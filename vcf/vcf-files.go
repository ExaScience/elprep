// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2017-2020 imec vzw.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version, and Additional Terms
// (see below).

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public
// License and Additional Terms along with this program. If not, see
// <https://github.com/ExaScience/elprep/blob/master/LICENSE.txt>.

package vcf

import (
	"bufio"
	"bytes"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"unicode/utf8"

	"github.com/exascience/pargo/pipeline"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/utils"
	"github.com/exascience/elprep/v5/utils/bgzf"
)

const (
	descriptionKey = "Description"
	idKey          = "ID"
	numberKey      = "Number"
	typeKey        = "Type"
)

var (
	specialToNormalString = strings.NewReplacer(
		"%3A", ":",
		"%3B", ";",
		"%3D", "=",
		"%25", "%",
		"%2C", ",",
		"%0D", "\r",
		"%0A", "\n",
		"%09", "\t",
	)
	normalToSpecialString = strings.NewReplacer(
		":", "%3A",
		";", "%3B",
		"=", "%3D",
		"%", "%25",
		",", "%2C",
		"\r", "%0D",
		"\n", "%0A",
		"\t", "%09",
	)
)

// ParseMetaField parses a VCF meta field
func (sc *StringScanner) ParseMetaField() (key, value string) {
	sc.SkipSpace()
	start := sc.index
	for ; sc.index < len(sc.data); sc.index++ {
		if c := sc.data[sc.index]; (c == ' ') || (c == '=') {
			break
		}
	}
	key = sc.data[start:sc.index]
	sc.SkipSpace()
	if sc.index >= len(sc.data) || sc.data[sc.index] != '=' {
		log.Panicf("invalid key=pair pair in a VCF meta-information line: %v", sc.data)
	}
	sc.index++
	sc.SkipSpace()
	start = sc.index
	if sc.data[sc.index] == '"' {
		start++
		sc.index++
		var buf strings.Builder
		for ; sc.index < len(sc.data); sc.index++ {
			switch sc.data[sc.index] {
			case '"':
				sc.index++
				return key, buf.String()
			case '\\':
				sc.index++
			}
			_ = buf.WriteByte(sc.data[sc.index])
		}
		sc.index = len(sc.data)
		log.Panicf("missing closing \" in a VCF meta-information line: %v", sc.data)
	}
	for ; sc.index < len(sc.data); sc.index++ {
		if c := sc.data[sc.index]; (c == ' ') || (c == ',') || (c == '>') {
			return key, sc.data[start:sc.index]
		}
	}
	log.Panicf("missing closing > in a VCF meta-information line: %v", sc.data)
	return "", ""
}

// ParseMetaInformation parses VCF meta information
func (sc *StringScanner) ParseMetaInformation() interface{} {
	if sc.data[sc.index] != '<' {
		start := sc.index
		sc.index = len(sc.data)
		return sc.data[start:]
	}
	sc.index++
	meta := NewMetaInformation()
	for {
		key, value := sc.ParseMetaField()
		switch key {
		case idKey:
			if meta.ID != nil {
				log.Panicf("multiple IDs in a VCF meta-information line: %v", sc.data)
			} else {
				meta.ID = utils.Intern(value)
			}
		case descriptionKey:
			if meta.Description != "" {
				log.Panicf("multiple Descriptions in a VCF meta-information line: %v", sc.data)
			} else {
				meta.Description = value
			}
		default:
			if !meta.Fields.SetUniqueEntry(key, value) {
				log.Panicf("duplicate field key %v in a VCF meta-information line: %v", key, sc.data)
			}
		}
		sc.SkipSpace()
		if c := sc.data[sc.index]; c == ',' {
			sc.index++
			continue
		} else if c == '>' {
			sc.index++
			break
		}
		log.Panicf("invalid syntax in a VCF meta-information line: %v", sc.data)
	}
	if meta.ID == nil {
		log.Panicf("missing ID in a VCF meta-information line: %v", sc.data)
	}
	return meta
}

// ParseFormatInformation parses VCF format information
func (sc *StringScanner) ParseFormatInformation() *FormatInformation {
	if sc.data[sc.index] != '<' {
		log.Panicf("Missing open angle bracket in a VCF INFO/FORMAT meta-information line: %v", sc.data)
	}
	sc.index++
	format := NewFormatInformation()
	for {
		key, value := sc.ParseMetaField()
		switch key {
		case idKey:
			if format.ID != nil {
				log.Panicf("multiple IDs in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			} else {
				format.ID = utils.Intern(value)
			}
		case descriptionKey:
			if format.Description != "" {
				log.Panicf("multiple Descriptions in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			} else {
				format.Description = value
			}
		case numberKey:
			if format.Number > InvalidNumber {
				log.Panicf("multiple Number entries in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			} else {
				switch value {
				case "a", "A":
					format.Number = NumberA
				case "r", "R":
					format.Number = NumberR
				case "g", "G":
					format.Number = NumberG
				case ".":
					format.Number = NumberDot
				default:
					format.Number = int32(internal.ParseInt(value, 10, 32))
				}
			}
		case typeKey:
			if format.Type != InvalidType {
				log.Panicf("Multiple types in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			} else {
				switch value {
				case "Integer":
					format.Type = Integer
				case "Float":
					format.Type = Float
				case "Flag":
					format.Type = Flag
				case "Character":
					format.Type = Character
				case "String":
					format.Type = String
				default:
					log.Panicf("Unknown type in a VCF INFO/FORMAT meta-information line: %v", sc.data)
				}
			}
		default:
			if !format.Fields.SetUniqueEntry(key, value) {
				log.Panicf("duplicate field key %v in a VCF meta-information line: %v", key, sc.data)
			}
		}
		sc.SkipSpace()
		if c := sc.data[sc.index]; c == ',' {
			sc.index++
			continue
		} else if c == '>' {
			sc.index++
			break
		}
		log.Panicf("invalid syntax in a VCF INFO/FORMAT meta-information line: %v", sc.data)
	}
	if format.ID == nil {
		log.Panicf("missing ID in a VCF INFO/FORMAT meta-information line: %v", sc.data)
	}
	if format.Number <= InvalidNumber {
		log.Panicf("missing number entry in a VCF INFO/FORMAT meta-information line: %v", sc.data)
	}
	if format.Type == InvalidType {
		log.Panicf("missing type in a VCF INFO/FORMAT meta-information line: %v", sc.data)
	}
	return format
}

func getLine(reader *bufio.Reader) string {
	line, err := reader.ReadString('\n')
	switch {
	case err == nil:
		if l := len(line); l > 1 && line[l-2] == '\r' {
			line = line[:l-2]
		} else {
			line = line[:l-1]
		}
	case err != io.EOF:
		log.Panic(err)
	}
	return line
}

// ParseHeader parses a VCF header
func ParseHeader(reader *bufio.Reader) (hdr *Header, lines int) {
	line := getLine(reader)
	lines++
	if line[:len(fileFormatVersionLinePrefix)] != fileFormatVersionLinePrefix {
		log.Panic("invalid first line in a VCF file")
	}
	hdr = NewHeader()
	hdr.FileFormat = line
	var sc StringScanner
	for {
		if data, err := reader.Peek(1); err != nil || data[0] != '#' {
			log.Panic("unexpected end of VCF header")
		}
		_, _ = reader.ReadByte()
		if data, err := reader.Peek(1); err != nil {
			log.Panic("unexpected end of VCF header")
		} else if data[0] != '#' {
			break
		}
		_, _ = reader.ReadByte()
		line = getLine(reader)
		lines++
		sc.Reset(line)
		if key, found := sc.readUntilByte('='); !found {
			log.Panic("invalid syntax in a VCF header")
		} else if key == "fileformat" {
			log.Panic("multiple file format meta-information lines in a VCF file")
		} else if key == "INFO" {
			hdr.Infos = append(hdr.Infos, sc.ParseFormatInformation())
		} else if key == "FORMAT" {
			hdr.Formats = append(hdr.Formats, sc.ParseFormatInformation())
		} else {
			hdr.Meta[key] = append(hdr.Meta[key], sc.ParseMetaInformation())
		}
	}
	line = getLine(reader)
	lines++
	sc.Reset(line)
	hdr.Columns = nil
	for sc.Len() > 0 {
		column, _ := sc.readUntilByte('\t')
		hdr.Columns = append(hdr.Columns, column)
	}
	return hdr, lines
}

func skipLine(reader *bufio.Reader) []byte {
	line, err := reader.ReadBytes('\n')
	if err != nil && err != io.EOF {
		log.Panic(err)
	}
	return line
}

// SkipHeader skips a VCF header. This is more efficient
// than calling ParseHeader and ignoring its result.
func SkipHeader(reader *bufio.Reader) (lines int) {
	line := skipLine(reader)
	lines++
	if string(line)[:len(fileFormatVersionLinePrefix)] != fileFormatVersionLinePrefix {
		log.Panic("invalid first line in a VCF file")
	}
	for {
		if data, err := reader.Peek(1); err != nil || data[0] != '#' {
			log.Panic("unexpected end of VCF header")
		}
		_, _ = reader.ReadByte()
		if data, err := reader.Peek(1); err != nil {
			log.Panic("unexpected end of VCF header")
		} else if data[0] != '#' {
			break
		}
		_, _ = reader.ReadByte()
		skipLine(reader)
		lines++
	}
	skipLine(reader)
	lines++
	return lines
}

// FieldParser is an abstraction for parsing VCF fields
type FieldParser func(*StringScanner) interface{}

func make1InfoParser(entryParser FieldParser) FieldParser {
	return func(sc *StringScanner) (result interface{}) {
		sc.SkipSpace()
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != '=') {
			log.Panicf("missing = in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			return nil
		}
		sc.SkipSpace()
		result = entryParser(sc)
		sc.SkipSpace()
		return
	}
}

func makeNInfoParser(entryParser FieldParser) FieldParser {
	return func(sc *StringScanner) interface{} {
		sc.SkipSpace()
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != '=') {
			log.Panicf("missing = in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			return nil
		}
		sc.SkipSpace()
		var result []interface{}
		for {
			result = append(result, entryParser(sc))
			sc.SkipSpace()
			if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ',') {
				break
			}
			sc.index++
		}
		return result
	}
}

var endOfInfoEntry = []byte{' ', ',', ';', '\t'}

// ParseGenericInfo parses a VCF info section without specific format information
func (sc *StringScanner) ParseGenericInfo() interface{} {
	sc.SkipSpace()
	if (sc.index < len(sc.data)) && (sc.data[sc.index] == '=') {
		var result []interface{}
		sc.index++
		sc.SkipSpace()
		for {
			result = append(result, sc.readUntilBytes(endOfInfoEntry))
			sc.SkipSpace()
			if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ',') {
				break
			}
			sc.index++
		}
		return result
	}
	return true
}

// ParseInfoInteger parses an integer in a VCF info section
func (sc *StringScanner) ParseInfoInteger() interface{} {
	return int(internal.ParseInt(sc.readUntilBytes(endOfInfoEntry), 10, 32))
}

// ParseInfoFloat parses a floating point number in a VCF info section
func (sc *StringScanner) ParseInfoFloat() interface{} {
	return internal.ParseFloat(sc.readUntilBytes(endOfInfoEntry), 64)
}

// ParseInfoFlag parses a boolean flag in a VCF info section (always returns true)
func (sc *StringScanner) ParseInfoFlag() interface{} {
	sc.SkipSpace()
	return true
}

// ParseInfoCharacter parses a rune in a VCF info section
func (sc *StringScanner) ParseInfoCharacter() interface{} {
	if sc.index >= len(sc.data) {
		log.Panic("missing Character entry in a VCF INFO meta-information line")
		return nil
	}
	if ch := sc.data[sc.index]; ch < utf8.RuneSelf {
		sc.index++
		return rune(ch)
	}
	rune, size := utf8.DecodeRune([]byte(sc.data[sc.index:]))
	if rune == utf8.RuneError {
		log.Panic("invalid rune encountered in a VCF INFO meta-information line")
	}
	sc.index += size
	return rune
}

// ParseInfoString parses a string in a VCF info section
func (sc *StringScanner) ParseInfoString() interface{} {
	return specialToNormalString.Replace(sc.readUntilBytes(endOfInfoEntry))
}

// CreateInfoParser creates a specific VCF info section parser for the given format information
func CreateInfoParser(format *FormatInformation) FieldParser {
	if format.Type == Flag {
		if format.Number != 0 {
			log.Panic("INFO Type Flag with Number != 0")
		}
		return (*StringScanner).ParseInfoFlag
	}
	if format.Number == 1 {
		switch format.Type {
		case Integer:
			return make1InfoParser((*StringScanner).ParseInfoInteger)
		case Float:
			return make1InfoParser((*StringScanner).ParseInfoFloat)
		case Character:
			return make1InfoParser((*StringScanner).ParseInfoCharacter)
		case String:
			return make1InfoParser((*StringScanner).ParseInfoString)
		default:
			log.Panic("invalid INFO Type")
			return nil
		}
	}
	switch format.Type {
	case Integer:
		return makeNInfoParser((*StringScanner).ParseInfoInteger)
	case Float:
		return makeNInfoParser((*StringScanner).ParseInfoFloat)
	case Character:
		return makeNInfoParser((*StringScanner).ParseInfoCharacter)
	case String:
		return makeNInfoParser((*StringScanner).ParseInfoString)
	default:
		log.Panic("invalid INFO Type")
		return nil
	}
}

var endOfFormatEntry = []byte{' ', ',', ':', '\t'}

// ParseFormatInteger parses an integer in a VCF format section
func (sc *StringScanner) ParseFormatInteger() interface{} {
	if sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		sc.index++
		return nil
	}
	return int(internal.ParseInt(sc.readUntilBytes(endOfFormatEntry), 10, 32))
}

// ParseFormatFloat parses a floating point number in a VCF format section
func (sc *StringScanner) ParseFormatFloat() interface{} {
	if sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		next := sc.index + 1
		if (next >= len(sc.data)) || (bytes.IndexByte(endOfFormatEntry, sc.data[next]) >= 0) {
			sc.index = next
			return nil
		}
	}
	return internal.ParseFloat(sc.readUntilBytes(endOfFormatEntry), 64)
}

// ParseFormatCharacter parses a rune in a VCF format section
func (sc *StringScanner) ParseFormatCharacter() interface{} {
	if sc.index >= len(sc.data) {
		return nil
	}
	if ch := sc.data[sc.index]; ch < utf8.RuneSelf {
		sc.index++
		if ch == '.' {
			return nil
		}
		return rune(ch)
	}
	rune, size := utf8.DecodeRune([]byte(sc.data[sc.index:]))
	if rune == utf8.RuneError {
		log.Panic("invalid rune encountered in a VCF FORMAT meta-information line")
	}
	sc.index += size
	return rune
}

// ParseFormatString parses a string in a VCF format section
func (sc *StringScanner) ParseFormatString() interface{} {
	if sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		next := sc.index + 1
		if (next >= len(sc.data)) || (bytes.IndexByte(endOfFormatEntry, sc.data[next]) >= 0) {
			sc.index = next
			return nil
		}
	}
	return specialToNormalString.Replace(sc.readUntilBytes(endOfFormatEntry))
}

func make1FormatParser(entryParser FieldParser) FieldParser {
	return func(sc *StringScanner) (result interface{}) {
		sc.SkipSpace()
		result = entryParser(sc)
		sc.SkipSpace()
		return
	}
}

func makeNFormatParser(entryParser FieldParser) FieldParser {
	return func(sc *StringScanner) interface{} {
		sc.SkipSpace()
		var result []interface{}
		for {
			result = append(result, entryParser(sc))
			sc.SkipSpace()
			if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ',') {
				break
			}
			sc.index++
		}
		return result
	}
}

// CreateFormatParser creates a specific VCF format section parser for the given format information
func CreateFormatParser(format *FormatInformation) FieldParser {
	if format.Number == 1 {
		switch format.Type {
		case Integer:
			return make1FormatParser((*StringScanner).ParseFormatInteger)
		case Float:
			return make1FormatParser((*StringScanner).ParseFormatFloat)
		case Character:
			return make1FormatParser((*StringScanner).ParseFormatCharacter)
		case String:
			return make1FormatParser((*StringScanner).ParseFormatString)
		default:
			log.Panic("invalid FORMAT Type")
			return nil
		}
	}
	switch format.Type {
	case Integer:
		return makeNFormatParser((*StringScanner).ParseFormatInteger)
	case Float:
		return makeNFormatParser((*StringScanner).ParseFormatFloat)
	case Character:
		return makeNFormatParser((*StringScanner).ParseFormatCharacter)
	case String:
		return makeNFormatParser((*StringScanner).ParseFormatString)
	default:
		log.Panic("invalid FORMAT Type")
		return nil
	}
}

func (sc *StringScanner) missingEntry() bool {
	if sc.index >= len(sc.data) {
		return true
	}
	if sc.data[sc.index] == '.' {
		next := sc.index + 1
		if (next >= len(sc.data)) || (sc.data[next] == '\t') {
			sc.index = next + 1
			return true
		}
	}
	return false
}

func (sc *StringScanner) scanTab() {
	if sc.index >= len(sc.data) || sc.data[sc.index] != '\t' {
		log.Panic("missing tabulator in VCF data line")
	}
	sc.index++
}

func (sc *StringScanner) maybeScanTab() bool {
	if sc.index >= len(sc.data) {
		return false
	}
	if sc.data[sc.index] != '\t' {
		log.Panic("missing tabulator in VCF data line")
	}
	sc.index++
	return true
}

func (sc *StringScanner) doString() string {
	if sc.missingEntry() {
		return "."
	}
	value, ok := sc.readUntilByte('\t')
	if !ok {
		log.Panic("missing tabulator in VCF data line")
		return ""
	}
	return value
}

func (sc *StringScanner) doInt32() int32 {
	if sc.missingEntry() {
		return -1
	}
	value, ok := sc.readUntilByte('\t')
	if !ok {
		log.Panic("missing tabulator in VCF data line")
		return -1
	}
	return int32(internal.ParseInt(value, 10, 32))
}

func (sc *StringScanner) doFloat() interface{} {
	if sc.missingEntry() {
		return nil
	}
	value, ok := sc.readUntilByte('\t')
	if !ok {
		log.Panic("missing tabulator in VCF data line")
		return nil
	}
	return internal.ParseFloat(value, 64)
}

func (sc *StringScanner) doStringList(separator []byte) (result []string) {
	if sc.missingEntry() {
		return nil
	}
	for {
		result = append(result, sc.readUntilBytes(separator))
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != separator[0]) {
			break
		}
		sc.index++
	}
	sc.scanTab()
	return result
}

var (
	filterSeparator = []byte{';', '\t'}
	passList        = []utils.Symbol{PASS}
)

func (sc *StringScanner) doFilter() []utils.Symbol {
	if sc.missingEntry() {
		return nil
	}
	str := sc.readUntilBytes(filterSeparator)
	if str == "PASS" {
		sc.scanTab()
		return passList
	}
	result := []utils.Symbol{utils.Intern(str)}
	for (sc.index < len(sc.data)) && (sc.data[sc.index] == ';') {
		sc.index++
		result = append(result, utils.Intern(sc.readUntilBytes(filterSeparator)))
	}
	sc.scanTab()
	return result
}

var formatSeparator = []byte{' ', ':', '\t'}

func (sc *StringScanner) doInfo(infoParsers utils.SmallMap) (result utils.SmallMap) {
	for {
		sc.SkipSpace()
		key := utils.Intern(sc.readUntilBytes(formatSeparator))
		var value interface{}
		if parser, ok := infoParsers.Get(key); ok {
			value = parser.(FieldParser)(sc)
		} else {
			value = sc.ParseGenericInfo()
		}
		result = append(result, utils.SmallMapEntry{Key: key, Value: value})
		sc.SkipSpace()
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ':') {
			return result
		}
		sc.index++
	}
}

func (sc *StringScanner) doSymbolList() (result []utils.Symbol) {
	for {
		sc.SkipSpace()
		str := sc.readUntilBytes(formatSeparator)
		result = append(result, utils.Intern(str))
		sc.SkipSpace()
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ':') {
			return result
		}
		sc.index++
	}
}

// VariantParser is an optimized parser for VCF variant lines.
//
// NSamples can be decreased as necessary to parse fewer samples, including down to zero.
type VariantParser struct {
	InfoParsers, FormatParsers utils.SmallMap
	NSamples                   int
}

// NewVariantParser creates a VariantParser for the given VCF header.
func (header *Header) NewVariantParser() *VariantParser {
	var vp VariantParser
	for _, format := range header.Infos {
		vp.InfoParsers = append(vp.InfoParsers, utils.SmallMapEntry{
			Key:   format.ID,
			Value: CreateInfoParser(format),
		})
	}
	for _, format := range header.Formats {
		vp.FormatParsers = append(vp.FormatParsers, utils.SmallMapEntry{
			Key:   format.ID,
			Value: CreateFormatParser(format),
		})
	}
	vp.NSamples = len(header.Columns) - len(DefaultHeaderColumns) - 1
	return &vp
}

var (
	idSeparator  = []byte{';', '\t'}
	altSeparator = []byte{',', '\t'}
)

// ParseVariant parses a VCF variant line
func (sc *StringScanner) ParseVariant(vp *VariantParser) Variant {
	var variant Variant
	variant.Chrom = sc.doString()
	variant.Pos = sc.doInt32()
	variant.ID = sc.doStringList(idSeparator)
	variant.Ref = sc.doString()
	variant.Alt = sc.doStringList(altSeparator)
	variant.Qual = sc.doFloat()
	variant.Filter = sc.doFilter()
	variant.Info = sc.doInfo(vp.InfoParsers)
	if vp.NSamples > 0 && sc.maybeScanTab() {
		variant.GenotypeFormat = sc.doSymbolList()
		if variant.GenotypeFormat[0] == GT {
			format := variant.GenotypeFormat[1:]
			parsers := make([]func(sc *StringScanner) interface{}, len(format))
			for p, f := range format {
				if parser, ok := vp.FormatParsers.Get(f); ok {
					parsers[p] = parser.(FieldParser)
				} else {
					parsers[p] = (*StringScanner).ParseFormatString
				}
			}
			for i := 0; i < vp.NSamples; i++ {
				if !sc.maybeScanTab() {
					break
				}
				var sample Genotype
				value := sc.ParseFormatString().(string)
				if len(value) > 0 {
					var sep byte
					for k := 0; k < len(value); k++ {
						if ch := value[k]; ch == '|' {
							sample.Phased = true
							sep = '|'
							break
						} else if ch == '/' {
							sample.Phased = false
							sep = '/'
							break
						}
					}
					sample.GT = make([]int32, 0, 2)
					for {
						end := strings.IndexByte(value, sep)
						if end < 0 {
							if value == "." {
								sample.GT = append(sample.GT, -1)
							} else {
								sample.GT = append(sample.GT, int32(internal.ParseInt(value, 10, 32)))
							}
							break
						}
						if sub := value[:end]; sub == "." {
							sample.GT = append(sample.GT, -1)
						} else {
							sample.GT = append(sample.GT, int32(internal.ParseInt(sub, 10, 32)))
						}
						value = value[end+1:]
					}
				}
				if (sc.index < len(sc.data)) && (sc.data[sc.index] == ':') {
					sc.index++
					sample.Data = make(utils.SmallMap, 0, len(parsers))
					for j := 0; j < len(parsers); j++ {
						key := format[j]
						value := parsers[j](sc)
						sample.Data = append(sample.Data, utils.SmallMapEntry{Key: key, Value: value})
						if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ':') {
							break
						}
						sc.index++
					}
				}
				variant.GenotypeData = append(variant.GenotypeData, sample)
			}
			return variant
		}
		parsers := make([]func(sc *StringScanner) interface{}, len(variant.GenotypeFormat))
		for p, f := range variant.GenotypeFormat {
			if parser, ok := vp.FormatParsers.Get(f); ok {
				parsers[p] = parser.(FieldParser)
			} else {
				parsers[p] = (*StringScanner).ParseFormatString
			}
		}
		for i := 0; i < vp.NSamples; i++ {
			if !sc.maybeScanTab() {
				break
			}
			var sample Genotype
			sample.Data = make(utils.SmallMap, 0, len(parsers))
			for j := 0; j < len(parsers); j++ {
				key := variant.GenotypeFormat[j]
				value := parsers[j](sc)
				sample.Data = append(sample.Data, utils.SmallMapEntry{Key: key, Value: value})
				if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ':') {
					break
				}
				sc.index++
			}
			variant.GenotypeData = append(variant.GenotypeData, sample)
		}
	}
	return variant
}

// FormatString outputs a string to a VCF file, adding necessary double quotes and escapes
func FormatString(out io.ByteWriter, str string) {
	internal.WriteByte(out, '"')
	for i := 0; i < len(str); i++ {
		b := str[i]
		if b == '"' || b == '\\' {
			internal.WriteByte(out, '\\')
		}
		internal.WriteByte(out, b)
	}
	internal.WriteByte(out, '"')
}

func needsQuotes(s string) bool {
	for i := 0; i < len(s); i++ {
		if ch := s[i]; ch == '"' || ch == ' ' {
			return true
		}
	}
	return false
}

// FormatMetaInformation outputs VCF meta information, which can be just a string or *MetaInformation
func FormatMetaInformation(out *bufio.Writer, meta interface{}) {
	switch m := meta.(type) {
	case string:
		internal.WriteString(out, m)
		internal.WriteByte(out, '\n')
	case *MetaInformation:
		internal.WriteString(out, "<ID=")
		internal.WriteString(out, *m.ID)
		for key, value := range m.Fields {
			internal.WriteByte(out, ',')
			internal.WriteString(out, key)
			internal.WriteByte(out, '=')
			if needsQuotes(value) {
				FormatString(out, value)
			} else {
				internal.WriteString(out, value)
			}
		}
		if m.Description != "" {
			internal.WriteString(out, ",Description=")
			FormatString(out, m.Description)
		}
		internal.WriteString(out, ">\n")
	default:
		log.Panic("invalid MetaInformation type")
	}
}

// FormatFormatInformation outputs VCF info or format information
func FormatFormatInformation(out *bufio.Writer, format *FormatInformation, infoNotFormat bool) {
	internal.WriteString(out, "<ID=")
	internal.WriteString(out, *format.ID)
	internal.WriteString(out, ",Number=")
	if format.Number >= 0 {
		internal.WriteString(out, strconv.FormatInt(int64(format.Number), 10))
	} else {
		switch format.Number {
		case NumberA:
			internal.WriteByte(out, 'A')
		case NumberR:
			internal.WriteByte(out, 'R')
		case NumberG:
			internal.WriteByte(out, 'G')
		case NumberDot:
			internal.WriteByte(out, '.')
		default:
			log.Panic("unknown Number kind in a VCF meta-information line")
		}
	}
	internal.WriteString(out, ",Type=")
	switch format.Type {
	case Integer:
		internal.WriteString(out, "Integer")
	case Float:
		internal.WriteString(out, "Float")
	case Flag:
		internal.WriteString(out, "Flag")
	case Character:
		internal.WriteString(out, "Character")
	case String:
		internal.WriteString(out, "String")
	default:
		log.Panic("invalid Type in a VCF meta-information line")
	}
	for key, value := range format.Fields {
		internal.WriteByte(out, ',')
		internal.WriteString(out, key)
		internal.WriteByte(out, '=')
		if (infoNotFormat && (key == "Source" || key == "Version")) || needsQuotes(value) {
			FormatString(out, value)
		} else {
			internal.WriteString(out, value)
		}
	}
	if format.Description != "" {
		internal.WriteString(out, ",Description=")
		FormatString(out, format.Description)
	}
	internal.WriteString(out, ">\n")
}

// Format outputs a VCF header
func (header *Header) Format(out *bufio.Writer) {
	internal.WriteString(out, header.FileFormat)
	internal.WriteByte(out, '\n')
	keys := []string{"FORMAT", "INFO"}
	for key := range header.Meta {
		keys = append(keys, key)
	}
	sort.Strings(keys)
	for _, key := range keys {
		switch key {
		case "FORMAT":
			for _, format := range header.Formats {
				internal.WriteString(out, "##FORMAT=")
				FormatFormatInformation(out, format, false)
			}
		case "INFO":
			for _, info := range header.Infos {
				internal.WriteString(out, "##INFO=")
				FormatFormatInformation(out, info, true)
			}
		default:
			for _, meta := range header.Meta[key] {
				internal.WriteString(out, "##")
				internal.WriteString(out, key)
				internal.WriteByte(out, '=')
				FormatMetaInformation(out, meta)
			}
		}
	}
	internal.WriteByte(out, '#')
	if len(header.Columns) > 0 {
		internal.WriteString(out, header.Columns[0])
		for _, col := range header.Columns[1:] {
			internal.WriteByte(out, '\t')
			internal.WriteString(out, col)
		}
	}
	internal.WriteByte(out, '\n')
}

func formatStringList(out []byte, list []string, separator byte) []byte {
	if len(list) == 0 {
		return append(out, '.', '\t')
	}
	out = append(out, list[0]...)
	for _, entry := range list[1:] {
		out = append(out, separator)
		out = append(out, entry...)
	}
	return append(out, '\t')
}

func formatSymbolList(out []byte, list []utils.Symbol, separator byte) []byte {
	if len(list) == 0 {
		return append(out, '.')
	}
	out = append(out, (*list[0])...)
	for _, sym := range list[1:] {
		out = append(out, separator)
		out = append(out, (*sym)...)
	}
	return out
}

func formatValue(out []byte, value interface{}) []byte {
	switch v := value.(type) {
	case int:
		return strconv.AppendInt(out, int64(v), 10)
	case float64:
		if v < 1 {
			if v < 0.01 {
				if math.Abs(v) < 1e-20 {
					return append(out, "0.00"...)
				}
				return strconv.AppendFloat(out, v, 'e', 3, 64)
			}
			return strconv.AppendFloat(out, v, 'f', 3, 64)
		}
		return strconv.AppendFloat(out, v, 'f', 2, 64)
	case rune:
		if v < utf8.RuneSelf {
			return append(out, byte(v))
		}
		pos := len(out)
		out = append(out, '1', '2', '3', '4', '5', '6')
		buf := out[pos:]
		return out[:pos+utf8.EncodeRune(buf, v)]
	case string:
		return append(out, normalToSpecialString.Replace(v)...)
	default:
		log.Panic("invalid value type")
		return nil
	}
}

func formatInfoEntry(out []byte, entry utils.SmallMapEntry) []byte {
	out = append(out, (*entry.Key)...)
	switch e := entry.Value.(type) {
	case bool:
		if !e {
			log.Panic("unexpected boolean value")
		}
		return out
	case []interface{}:
		out = append(out, '=')
		if len(e) == 0 {
			return out
		}
		out = formatValue(out, e[0])
		for _, v := range e[1:] {
			out = formatValue(append(out, ','), v)
		}
		return out
	default:
		return formatValue(append(out, '='), entry.Value)
	}
}

func formatInfo(out []byte, info utils.SmallMap) []byte {
	if len(info) == 0 {
		return append(out, '.')
	}
	out = formatInfoEntry(out, info[0])
	for _, entry := range info[1:] {
		out = formatInfoEntry(append(out, ';'), entry)
	}
	return out
}

func formatGenotypeDataEntry(out []byte, format utils.Symbol, data utils.SmallMap) ([]byte, bool) {
	switch value, _ := data.Get(format); val := value.(type) {
	case nil:
		return append(out, '.'), false
	case []interface{}:
		if len(val) == 0 {
			return out, true
		}
		if val[0] == nil {
			out = append(out, '.')
		} else {
			out = formatValue(out, val[0])
		}
		for _, v := range val[1:] {
			out = append(out, ',')
			if v == nil {
				out = append(out, '.')
			} else {
				out = formatValue(out, v)
			}
		}
		return out, true
	default:
		return formatValue(out, value), true
	}
}

func formatGenotypeData(out []byte, format []utils.Symbol, g Genotype) []byte {
	if len(format) == 0 {
		return out
	}
	var pos int
	var ok bool
	if format[0] == GT {
		var sep byte
		if g.Phased {
			sep = '|'
		} else {
			sep = '/'
		}
		if n := g.GT[0]; n < 0 {
			out = append(out, '.')
		} else {
			out = strconv.AppendInt(out, int64(n), 10)
		}
		for _, n := range g.GT[1:] {
			out = append(out, sep)
			if n < 0 {
				out = append(out, '.')
			} else {
				out = strconv.AppendInt(out, int64(n), 10)
			}
		}
		pos = len(out)
	} else {
		pos = len(out)
		out, ok = formatGenotypeDataEntry(out, format[0], g.Data)
		if ok {
			pos = len(out)
		}
	}
	for _, f := range format[1:] {
		out = append(out, ':')
		out, ok = formatGenotypeDataEntry(out, f, g.Data)
		if ok {
			pos = len(out)
		}
	}
	return out[:pos]
}

// Format outputs a VCF variant line
func (variant Variant) Format(out []byte) []byte {
	out = append(append(out, variant.Chrom...), '\t')
	if variant.Pos < 0 {
		out = append(out, '.', '\t')
	} else {
		out = append(strconv.AppendInt(out, int64(variant.Pos), 10), '\t')
	}
	out = formatStringList(out, variant.ID, ';')
	out = append(append(out, variant.Ref...), '\t')
	out = formatStringList(out, variant.Alt, ',')
	if value, ok := variant.Qual.(float64); ok {
		/*
			if math.Floor(value) != value {
				out = append(strconv.AppendFloat(out, value, 'f', 2, 64), '\t')
			} else {
				out = append(strconv.AppendFloat(out, value, 'f', -1, 64), '\t')
			}
		*/
		out = strconv.AppendFloat(out, value, 'f', 2, 64)
		if bytes.HasSuffix(out, []byte(".00")) {
			out = out[:len(out)-3]
		}
		out = append(out, '\t')
	} else {
		out = append(out, '.', '\t')
	}
	if len(variant.Filter) == 0 {
		out = append(out, '.', '\t')
	} else {
		out = append(formatSymbolList(out, variant.Filter, ';'), '\t')
	}
	out = formatInfo(out, variant.Info)
	if len(variant.GenotypeFormat) > 0 {
		out = append(out, '\t')
		out = formatSymbolList(out, variant.GenotypeFormat, ':')
		for _, data := range variant.GenotypeData {
			out = formatGenotypeData(append(out, '\t'), variant.GenotypeFormat, data)
		}
	}
	return append(out, '\n')
}

func FormatVariants(out *bufio.Writer, variants []Variant) {
	var p pipeline.Pipeline
	p.Source(variants)
	p.Add(
		pipeline.LimitedPar(0, pipeline.Receive(func(_ int, data interface{}) interface{} {
			variants := data.([]Variant)
			records := make([][]byte, 0, len(variants))
			var buf []byte
			for _, variant := range variants {
				buf = variant.Format(buf)
				records = append(records, append([]byte(nil), buf...))
				buf = buf[:0]
			}
			return records
		})),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			records := data.([][]byte)
			for _, record := range records {
				_, _ = out.Write(record)
			}
			return nil
		})),
	)
	internal.RunPipeline(&p)
}

// Format outputs a full VCF struct
func (vcf *Vcf) Format(out *bufio.Writer) {
	vcf.Header.Format(out)
	FormatVariants(out, vcf.Variants)
}

// The possible file extensions for VCF or gz-compressed VCF files
const (
	VcfExt = ".vcf"
	GzExt  = ".gz"
)

// InputFile represents a VCF or BCF file for input.
type InputFile struct {
	rc   io.ReadCloser
	bgzf *bgzf.Reader
	*bufio.Reader
}

// OutputFile represents a VCF or BCF file for output.
type OutputFile struct {
	wc   io.WriteCloser
	bgzf *bgzf.Writer
	*bufio.Writer
}

// Open a VCF file for input.
//
// Whether the format is gzipped or not is determined from the content
// of the input, not from any file extensions.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func Open(name string) *InputFile {
	var file io.ReadCloser
	if name == "/dev/stdin" {
		file = os.Stdin
	} else {
		file = internal.FileOpen(name)
	}
	buf := bufio.NewReader(file)
	ok, err := bgzf.IsGzip(buf)
	if err != nil {
		log.Panic(err)
	}
	if !ok {
		return &InputFile{
			rc:     file,
			Reader: buf,
		}
	}
	bgzf, err := bgzf.NewReader(buf)
	if err != nil {
		log.Panic(err)
	}
	return &InputFile{
		rc:     file,
		bgzf:   bgzf,
		Reader: bufio.NewReader(bgzf),
	}
}

// Open a VCF file for input, returning false if it doesn't exist.
//
// Whether the format is gzipped or not is determined from the content
// of the input, not from any file extensions.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func OpenIfExists(name string) (*InputFile, bool) {
	var file io.ReadCloser
	if name == "/dev/stdin" {
		file = os.Stdin
	} else if f, ok := internal.FileOpenIfExists(name); ok {
		file = f
	} else {
		return nil, false
	}
	buf := bufio.NewReader(file)
	ok, err := bgzf.IsGzip(buf)
	if err != nil {
		log.Panic(err)
	}
	if !ok {
		return &InputFile{
			rc:     file,
			Reader: buf,
		}, true
	}
	bgzf, err := bgzf.NewReader(buf)
	if err != nil {
		log.Panic(err)
	}
	return &InputFile{
		rc:     file,
		bgzf:   bgzf,
		Reader: bufio.NewReader(bgzf),
	}, true
}

// Create a VCF file for output.
//
// The format string can be "vcf" or "gz". If the format string
// is empty, the output format is determined by looking at the
// filename extension. If the filename extension is not .gz,
// then .vcf is always assumed.
//
// The format string will not become part of the resulting filename.
//
// Following zlib, levels range from 1 (BestSpeed) to 9 (BestCompression);
// higher levels typically run slower but compress more. Level 0
// (NoCompression) does not attempt any compression; it only adds the
// necessary DEFLATE framing.
// Level -1 (DefaultCompression) uses the default compression level.
// Level -2 (HuffmanOnly) will use Huffman compression only, giving
// a very fast compression for all types of input, but sacrificing considerable
// compression efficiency.
//
// If the name is "/dev/stdout", then the output is written to
// os.Stdout.
func Create(name string, format string, level int) *OutputFile {
	var file io.WriteCloser
	if name == "/dev/stdout" {
		file = os.Stdout
	} else {
		file = internal.FileCreate(name)
	}
	if format == "" {
		format = filepath.Ext(name)
	}
	switch strings.ToLower(format) {
	case "gz", ".gz", "vcf.gz", ".vcf.gz":
		bgzf := bgzf.NewWriter(file, level)
		return &OutputFile{
			wc:     file,
			bgzf:   bgzf,
			Writer: bufio.NewWriter(bgzf),
		}
	case "bcf", ".bcf":
		log.Panicf("BCF format not supported when creating %v", name)
		return nil
	default:
		return &OutputFile{
			wc:     file,
			Writer: bufio.NewWriter(file),
		}
	}
}

// Close the VCF input file.
func (input *InputFile) Close() {
	if input.bgzf != nil {
		internal.Close(input.bgzf)
	}
	if input.rc != os.Stdin {
		internal.Close(input.rc)
	}
}

// Close the VCF input file.
func (output *OutputFile) Close() {
	err := output.Flush()
	if err != nil {
		log.Panic(err)
	}
	if output.bgzf != nil {
		internal.Close(output.bgzf)
	}
	if output.wc != os.Stdout {
		internal.Close(output.wc)
	}
}

// ParseVariants parses VCF variant lines based on the given VCF header.
func (header *Header) ParseVariants(input *InputFile) []Variant {
	variantParser := header.NewVariantParser()
	var variants []Variant
	var p pipeline.Pipeline
	p.Source(pipeline.NewScanner(input.Reader))
	p.Add(
		pipeline.LimitedPar(0, pipeline.Receive(func(_ int, data interface{}) interface{} {
			lines := data.([]string)
			variants := make([]Variant, 0, len(lines))
			var sc StringScanner
			for _, line := range lines {
				sc.Reset(line)
				variants = append(variants, sc.ParseVariant(variantParser))
			}
			return variants
		})),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			variants = append(variants, data.([]Variant)...)
			return nil
		})),
	)
	internal.RunPipeline(&p)
	return variants
}

// Parse parseses a full VCF files.
func (input *InputFile) Parse() *Vcf {
	header, _ := ParseHeader(input.Reader)
	variants := header.ParseVariants(input)
	return &Vcf{header, variants}
}

// Format outputs a full VCF struct.
func (output *OutputFile) Format(vcf *Vcf) {
	vcf.Format(output.Writer)
}
