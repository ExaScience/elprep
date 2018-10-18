// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017, 2018 imec vzw.

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
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"unicode/utf8"

	"github.com/exascience/elprep/v4/utils"
)

const (
	descriptionKey = "Description"
	idKey          = "ID"
	numberKey      = "Number"
	typeKey        = "Type"
)

// ParseMetaField parses a VCF meta field
func (sc *StringScanner) ParseMetaField() (key, value string) {
	if sc.err != nil {
		return
	}
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
		if sc.err == nil {
			sc.err = fmt.Errorf("invalid key=pair pair in a VCF meta-information line: %v", sc.data)
		}
		return
	}
	sc.index++
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
		if sc.err == nil {
			sc.err = fmt.Errorf("missing closing \" in a VCF meta-information line: %v", sc.data)
		}
		return key, buf.String()
	}
	for ; sc.index < len(sc.data); sc.index++ {
		if c := sc.data[sc.index]; (c == ' ') || (c == ',') || (c == '>') {
			return key, sc.data[start:sc.index]
		}
	}
	if sc.err == nil {
		sc.err = fmt.Errorf("missing closing > in a VCF meta-information line: %v", sc.data)
	}
	return key, sc.data[start:]
}

// ParseMetaInformation parses VCF meta information
func (sc *StringScanner) ParseMetaInformation() interface{} {
	if sc.err != nil {
		return nil
	}
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
				if sc.err == nil {
					sc.err = fmt.Errorf("multiple IDs in a VCF meta-information line: %v", sc.data)
				}
			} else {
				meta.ID = utils.Intern(value)
			}
		case descriptionKey:
			if meta.Description != "" {
				if sc.err == nil {
					sc.err = fmt.Errorf("multiple Descriptions in a VCF meta-information line: %v", sc.data)
				}
			} else {
				meta.Description = value
			}
		default:
			if !meta.Fields.SetUniqueEntry(key, value) {
				if sc.err == nil {
					sc.err = fmt.Errorf("duplicate field key %v in a VCF meta-information line: %v", key, sc.data)
				}
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
		if sc.err == nil {
			sc.err = fmt.Errorf("invalid syntax in a VCF meta-information line: %v", sc.data)
		}
		break
	}
	if meta.ID == nil {
		if sc.err == nil {
			sc.err = fmt.Errorf("missing ID in a VCF meta-information line: %v", sc.data)
		}
	}
	return meta
}

// ParseFormatInformation parses VCF format information
func (sc *StringScanner) ParseFormatInformation() *FormatInformation {
	if sc.err != nil {
		return nil
	}
	if sc.data[sc.index] != '<' {
		sc.err = fmt.Errorf("Missing open angle bracket in a VCF INFO/FORMAT meta-information line: %v", sc.data)
		return nil
	}
	sc.index++
	format := NewFormatInformation()
	for {
		key, value := sc.ParseMetaField()
		switch key {
		case idKey:
			if format.ID != nil {
				if sc.err == nil {
					sc.err = fmt.Errorf("multiple IDs in a VCF INFO/FORMAT meta-information line: %v", sc.data)
				}
			} else {
				format.ID = utils.Intern(value)
			}
		case descriptionKey:
			if format.Description != "" {
				if sc.err == nil {
					sc.err = fmt.Errorf("multiple Descriptions in a VCF INFO/FORMAT meta-information line: %v", sc.data)
				}
			} else {
				format.Description = value
			}
		case numberKey:
			if format.Number > InvalidNumber {
				if sc.err == nil {
					sc.err = fmt.Errorf("multiple Number entries in a VCF INFO/FORMAT meta-information line: %v", sc.data)
				}
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
					n, err := strconv.ParseInt(value, 10, 32)
					if err != nil {
						if sc.err == nil {
							sc.err = err
						}
					} else {
						format.Number = int32(n)
					}
				}
			}
		case typeKey:
			if format.Type != InvalidType {
				if sc.err == nil {
					sc.err = fmt.Errorf("Multiple types in a VCF INFO/FORMAT meta-information line: %v", sc.data)
				}
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
					if sc.err == nil {
						sc.err = fmt.Errorf("Unknown type in a VCF INFO/FORMAT meta-information line: %v", sc.data)
					}
				}
			}
		default:
			if !format.Fields.SetUniqueEntry(key, value) {
				if sc.err == nil {
					sc.err = fmt.Errorf("duplicate field key %v in a VCF meta-information line: %v", key, sc.data)
				}
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
		if sc.err == nil {
			sc.err = fmt.Errorf("invalid syntax in a VCF INFO/FORMAT meta-information line: %v", sc.data)
		}
		break
	}
	if format.ID == nil {
		if sc.err == nil {
			sc.err = fmt.Errorf("missing ID in a VCF INFO/FORMAT meta-information line: %v", sc.data)
		}
	}
	if format.Number <= InvalidNumber {
		if sc.err == nil {
			sc.err = fmt.Errorf("missing number entry in a VCF INFO/FORMAT meta-information line: %v", sc.data)
		}
	}
	if format.Type == InvalidType {
		if sc.err == nil {
			sc.err = fmt.Errorf("missing type in a VCF INFO/FORMAT meta-information line: %v", sc.data)
		}
	}
	return format
}

func getLine(reader *bufio.Reader) (line string, err error) {
	line, err = reader.ReadString('\n')
	switch {
	case err == nil:
		line = line[:len(line)-1]
	case err == io.EOF:
		err = nil
	}
	return
}

// ParseHeader parses a VCF header
func ParseHeader(reader *bufio.Reader) (hdr *Header, lines int, err error) {
	line, err := getLine(reader)
	if err != nil {
		return nil, 0, err
	}
	lines++
	if line[:len(fileFormatVersionLinePrefix)] != fileFormatVersionLinePrefix {
		return nil, 0, errors.New("invalid first line in a VCF file")
	}
	hdr = NewHeader()
	hdr.FileFormat = line
	var sc StringScanner
	for {
		if data, e := reader.Peek(1); (e != nil) || (data[0] != '#') {
			return nil, 0, errors.New("unexpected end of VCF header")
		}
		_, _ = reader.ReadByte()
		if data, e := reader.Peek(1); e != nil {
			return nil, 0, errors.New("unexpected end of VCF header")
		} else if data[0] != '#' {
			break
		}
		_, _ = reader.ReadByte()
		line, err = getLine(reader)
		if err != nil {
			return nil, 0, err
		}
		lines++
		sc.Reset(line)
		if key, found := sc.readUntilByte('='); !found {
			return nil, 0, errors.New("invalid syntax in a VCF header")
		} else if key == "fileformat" {
			return nil, 0, errors.New("multiple file format meta-information lines in a VCF file")
		} else if key == "INFO" {
			hdr.Infos = append(hdr.Infos, sc.ParseFormatInformation())
		} else if key == "FORMAT" {
			hdr.Formats = append(hdr.Formats, sc.ParseFormatInformation())
		} else {
			hdr.Meta[key] = append(hdr.Meta[key], sc.ParseMetaInformation())
		}
		if sc.err != nil {
			return nil, 0, sc.err
		}
	}
	line, err = getLine(reader)
	if err != nil {
		return nil, 0, err
	}
	lines++
	sc.Reset(line)
	for sc.Len() > 0 {
		column, _ := sc.readUntilByte('\t')
		hdr.Columns = append(hdr.Columns, column)
	}
	if sc.err != nil {
		return nil, 0, sc.err
	}
	return hdr, lines, nil
}

// FieldParser is an abstraction for parsing VCF fields
type FieldParser func(*StringScanner) interface{}

func make1InfoParser(entryParser FieldParser) FieldParser {
	return func(sc *StringScanner) (result interface{}) {
		sc.SkipSpace()
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != '=') {
			if sc.err == nil {
				sc.err = fmt.Errorf("missing = in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			}
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
			if sc.err == nil {
				sc.err = fmt.Errorf("missing = in a VCF INFO/FORMAT meta-information line: %v", sc.data)
			}
			return nil
		}
		sc.SkipSpace()
		var result []interface{}
		for sc.err == nil {
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
	if sc.err != nil {
		return nil
	}
	sc.SkipSpace()
	if (sc.index < len(sc.data)) && (sc.data[sc.index] == '=') {
		var result []interface{}
		sc.index++
		sc.SkipSpace()
		for sc.err == nil {
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
	i, err := strconv.ParseInt(sc.readUntilBytes(endOfInfoEntry), 10, 32)
	if err != nil {
		if sc.err == nil {
			sc.err = err
		}
		return nil
	}
	return int(i)
}

// ParseInfoFloat parses a floating point number in a VCF info section
func (sc *StringScanner) ParseInfoFloat() interface{} {
	f, err := strconv.ParseFloat(sc.readUntilBytes(endOfInfoEntry), 64)
	if err != nil {
		if sc.err == nil {
			sc.err = err
		}
		return nil
	}
	return f
}

// ParseInfoFlag parses a boolean flag in a VCF info section (always returns true)
func (sc *StringScanner) ParseInfoFlag() interface{} {
	sc.SkipSpace()
	return true
}

// ParseInfoCharacter parses a rune in a VCF info section
func (sc *StringScanner) ParseInfoCharacter() interface{} {
	if sc.err != nil {
		return nil
	}
	if sc.index >= len(sc.data) {
		sc.err = errors.New("missing Character entry in a VCF INFO meta-information line")
		return nil
	}
	if ch := sc.data[sc.index]; ch < utf8.RuneSelf {
		sc.index++
		return rune(ch)
	}
	rune, size := utf8.DecodeRune([]byte(sc.data[sc.index:]))
	if rune == utf8.RuneError {
		if sc.err == nil {
			sc.err = errors.New("invalid rune encountered in a VCF INFO meta-information line")
		}
	}
	sc.index += size
	return rune
}

// ParseInfoString parses a string in a VCF info section
func (sc *StringScanner) ParseInfoString() interface{} {
	return sc.readUntilBytes(endOfInfoEntry)
}

// CreateInfoParser creates a specific VCF info section parser for the given format information
func CreateInfoParser(format *FormatInformation) (FieldParser, error) {
	if format.Type == Flag {
		if format.Number != 0 {
			return nil, errors.New("INFO Type Flag with Number != 0")
		}
		return (*StringScanner).ParseInfoFlag, nil
	}
	if format.Number == 1 {
		switch format.Type {
		case Integer:
			return make1InfoParser((*StringScanner).ParseInfoInteger), nil
		case Float:
			return make1InfoParser((*StringScanner).ParseInfoFloat), nil
		case Character:
			return make1InfoParser((*StringScanner).ParseInfoCharacter), nil
		case String:
			return make1InfoParser((*StringScanner).ParseInfoString), nil
		default:
			return nil, errors.New("invalid INFO Type")
		}
	}
	switch format.Type {
	case Integer:
		return makeNInfoParser((*StringScanner).ParseInfoInteger), nil
	case Float:
		return makeNInfoParser((*StringScanner).ParseInfoFloat), nil
	case Character:
		return makeNInfoParser((*StringScanner).ParseInfoCharacter), nil
	case String:
		return makeNInfoParser((*StringScanner).ParseInfoString), nil
	default:
		return nil, errors.New("invalid INFO Type")
	}
}

var endOfFormatEntry = []byte{' ', ',', ':', '\t'}

// ParseGenericFormat parses a VCF format section without specific format information
func (sc *StringScanner) ParseGenericFormat() interface{} {
	if sc.err != nil {
		return nil
	}
	sc.SkipSpace()
	if (sc.index < len(sc.data)) && (sc.data[sc.index] == '=') {
		var result []interface{}
		sc.index++
		sc.SkipSpace()
		for sc.err == nil {
			result = append(result, sc.readUntilBytes(endOfFormatEntry))
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

// ParseFormatInteger parses an integer in a VCF format section
func (sc *StringScanner) ParseFormatInteger() interface{} {
	if sc.err != nil || sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		sc.index++
		return nil
	}
	i, err := strconv.ParseInt(sc.readUntilBytes(endOfFormatEntry), 10, 32)
	if err != nil {
		if sc.err == nil {
			sc.err = err
		}
		return nil
	}
	return int(i)
}

func containsByte(b byte, bytes []byte) bool {
	for _, bb := range bytes {
		if b == bb {
			return true
		}
	}
	return false
}

// ParseFormatFloat parses a floating point number in a VCF format section
func (sc *StringScanner) ParseFormatFloat() interface{} {
	if sc.err != nil || sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		next := sc.index + 1
		if (next >= len(sc.data)) || containsByte(sc.data[next], endOfFormatEntry) {
			sc.index = next
			return nil
		}
	}
	f, err := strconv.ParseFloat(sc.readUntilBytes(endOfFormatEntry), 64)
	if err != nil {
		if sc.err == nil {
			sc.err = err
		}
		return nil
	}
	return f
}

// ParseFormatCharacter parses a rune in a VCF format section
func (sc *StringScanner) ParseFormatCharacter() interface{} {
	if sc.err != nil || sc.index >= len(sc.data) {
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
		if sc.err == nil {
			sc.err = errors.New("invalid rune encountered in a VCF FORMAT meta-information line")
		}
	}
	sc.index += size
	return rune
}

// ParseFormatString parses a string in a VCF format section
func (sc *StringScanner) ParseFormatString() interface{} {
	if sc.err != nil || sc.index >= len(sc.data) {
		return nil
	}
	if sc.data[sc.index] == '.' {
		next := sc.index + 1
		if (next >= len(sc.data)) || containsByte(sc.data[next], endOfFormatEntry) {
			sc.index = next
			return nil
		}
	}
	return sc.readUntilBytes(endOfFormatEntry)
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
		for sc.err == nil {
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
func CreateFormatParser(format *FormatInformation) (FieldParser, error) {
	if format.Number == 1 {
		switch format.Type {
		case Integer:
			return make1FormatParser((*StringScanner).ParseFormatInteger), nil
		case Float:
			return make1FormatParser((*StringScanner).ParseFormatFloat), nil
		case Character:
			return make1FormatParser((*StringScanner).ParseFormatCharacter), nil
		case String:
			return make1FormatParser((*StringScanner).ParseFormatString), nil
		default:
			return nil, errors.New("invalid FORMAT Type")
		}
	}
	switch format.Type {
	case Integer:
		return makeNFormatParser((*StringScanner).ParseFormatInteger), nil
	case Float:
		return makeNFormatParser((*StringScanner).ParseFormatFloat), nil
	case Character:
		return makeNFormatParser((*StringScanner).ParseFormatCharacter), nil
	case String:
		return makeNFormatParser((*StringScanner).ParseFormatString), nil
	default:
		return nil, errors.New("invalid FORMAT Type")
	}
}

func (sc *StringScanner) missingEntry() bool {
	if (sc.err != nil) || (sc.index >= len(sc.data)) {
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

func (sc *StringScanner) scanChar(ch byte) {
	if sc.err != nil {
		return
	}
	if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ch) {
		sc.err = errors.New("missing tabulator in VCF data line")
	}
	sc.index++
}

func (sc *StringScanner) doString() string {
	if sc.missingEntry() {
		return "."
	}
	value, ok := sc.readUntilByte('\t')
	if !ok {
		if sc.err == nil {
			sc.err = errors.New("missing tabulator in VCF data line")
		}
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
		if sc.err == nil {
			sc.err = errors.New("missing tabulator in VCF data line")
		}
		return -1
	}
	i, err := strconv.ParseInt(value, 10, 32)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return int32(i)
}

func (sc *StringScanner) doFloat() interface{} {
	if sc.missingEntry() {
		return nil
	}
	value, ok := sc.readUntilByte('\t')
	if !ok {
		if sc.err == nil {
			sc.err = errors.New("missing tabulator in VCF data line")
		}
		return nil
	}
	f, err := strconv.ParseFloat(value, 64)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return f
}

func (sc *StringScanner) doStringList(separator []byte) (result []string) {
	if sc.missingEntry() {
		return nil
	}
	for sc.err == nil {
		result = append(result, sc.readUntilBytes(separator))
		if (sc.index >= len(sc.data)) || (sc.data[sc.index] != separator[0]) {
			break
		}
		sc.index++
	}
	sc.scanChar('\t')
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
		sc.scanChar('\t')
		return passList
	}
	result := []utils.Symbol{utils.Intern(str)}
	for (sc.err == nil) && (sc.index < len(sc.data)) && (sc.data[sc.index] == ';') {
		sc.index++
		result = append(result, utils.Intern(sc.readUntilBytes(filterSeparator)))
	}
	sc.scanChar('\t')
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
		if sc.err != nil {
			return nil
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
		if sc.err != nil {
			return nil
		}
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
func (header *Header) NewVariantParser() (*VariantParser, error) {
	var vp VariantParser
	for _, format := range header.Infos {
		parser, err := CreateInfoParser(format)
		if err != nil {
			return nil, err
		}
		vp.InfoParsers = append(vp.InfoParsers, utils.SmallMapEntry{Key: format.ID, Value: parser})
	}
	for _, format := range header.Formats {
		parser, err := CreateFormatParser(format)
		if err != nil {
			return nil, err
		}
		vp.FormatParsers = append(vp.FormatParsers, utils.SmallMapEntry{Key: format.ID, Value: parser})
	}
	vp.NSamples = len(header.Columns) - len(DefaultHeaderColumns) - 1
	return &vp, nil
}

var (
	idSeparator  = []byte{';', '\t'}
	altSeparator = []byte{',', '\t'}
)

// ParseVariant parses a VCF variant line
func (sc *StringScanner) ParseVariant(vp *VariantParser) *Variant {
	var variant Variant
	variant.Chrom = sc.doString()
	variant.Pos = sc.doInt32()
	variant.ID = sc.doStringList(idSeparator)
	variant.Ref = sc.doString()
	variant.Alt = sc.doStringList(altSeparator)
	variant.Qual = sc.doFloat()
	variant.Filter = sc.doFilter()
	variant.Info = sc.doInfo(vp.InfoParsers)
	if vp.NSamples > 0 {
		sc.scanChar('\t')
		variant.GenotypeFormat = sc.doSymbolList()
		parsers := make([]func(sc *StringScanner) interface{}, len(variant.GenotypeFormat))
		for p, format := range variant.GenotypeFormat {
			if parser, ok := vp.FormatParsers.Get(format); ok {
				parsers[p] = parser.(FieldParser)
			} else {
				parsers[p] = (*StringScanner).ParseGenericFormat
			}
		}
		for i := 0; i < vp.NSamples; i++ {
			sample := make(utils.SmallMap, 0, len(parsers))
			sc.scanChar('\t')
			for j := 0; j < len(parsers); j++ {
				key := variant.GenotypeFormat[j]
				value := parsers[j](sc)
				if sc.err != nil {
					return nil
				}
				sample = append(sample, utils.SmallMapEntry{Key: key, Value: value})
				if (sc.index >= len(sc.data)) || (sc.data[sc.index] != ':') {
					break
				}
				sc.index++
			}
			variant.GenotypeData = append(variant.GenotypeData, sample)
		}
	}
	if sc.err != nil {
		return nil
	}
	return &variant
}

// FormatString outputs a string to a VCF file, adding necessary double quotes and escapes
func FormatString(out io.ByteWriter, str string) error {
	_ = out.WriteByte('"')
	for i := 0; i < len(str); i++ {
		b := str[i]
		if b == '"' || b == '\\' {
			_ = out.WriteByte('\\')
		}
		_ = out.WriteByte(b)
	}
	return out.WriteByte('"')
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
func FormatMetaInformation(out *bufio.Writer, meta interface{}) error {
	switch m := meta.(type) {
	case string:
		_, _ = out.WriteString(m)
		return out.WriteByte('\n')
	case *MetaInformation:
		_, _ = out.WriteString("<ID=")
		_, _ = out.WriteString(*m.ID)
		for key, value := range m.Fields {
			_ = out.WriteByte(',')
			_, _ = out.WriteString(key)
			_ = out.WriteByte('=')
			if needsQuotes(value) {
				_ = FormatString(out, value)
			} else {
				_, _ = out.WriteString(value)
			}
		}
		if m.Description != "" {
			_, _ = out.WriteString(",Description=")
			_ = FormatString(out, m.Description)
		}
		_, err := out.WriteString(">\n")
		return err
	default:
		return errors.New("invalid MetaInformation type")
	}
}

// FormatFormatInformation outputs VCF info or format information
func FormatFormatInformation(out *bufio.Writer, format *FormatInformation, infoNotFormat bool) error {
	_, _ = out.WriteString("<ID=")
	_, _ = out.WriteString(*format.ID)
	_, _ = out.WriteString(",Number=")
	if format.Number >= 0 {
		_, _ = out.WriteString(strconv.FormatInt(int64(format.Number), 10))
	} else {
		switch format.Number {
		case NumberA:
			_ = out.WriteByte('A')
		case NumberR:
			_ = out.WriteByte('R')
		case NumberG:
			_ = out.WriteByte('G')
		case NumberDot:
			_ = out.WriteByte('.')
		default:
			return errors.New("unknown Number kind in a VCF meta-information line")
		}
	}
	_, _ = out.WriteString(",Type=")
	switch format.Type {
	case Integer:
		_, _ = out.WriteString("Integer")
	case Float:
		_, _ = out.WriteString("Float")
	case Flag:
		_, _ = out.WriteString("Flag")
	case Character:
		_, _ = out.WriteString("Character")
	case String:
		_, _ = out.WriteString("String")
	default:
		return errors.New("invalid Type in a VCF meta-information line")
	}
	for key, value := range format.Fields {
		_ = out.WriteByte(',')
		_, _ = out.WriteString(key)
		_ = out.WriteByte('=')
		if (infoNotFormat && (key == "Source" || key == "Version")) || needsQuotes(value) {
			_ = FormatString(out, value)
		} else {
			_, _ = out.WriteString(value)
		}
	}
	if format.Description != "" {
		_, _ = out.WriteString(",Description=")
		_ = FormatString(out, format.Description)
	}
	_, err := out.WriteString(">\n")
	return err
}

// Format outputs a VCF header
func (header *Header) Format(out *bufio.Writer) (err error) {
	_, _ = out.WriteString(header.FileFormat)
	_ = out.WriteByte('\n')
	for _, info := range header.Infos {
		_, _ = out.WriteString("##INFO=")
		_ = FormatFormatInformation(out, info, true)
	}
	for _, format := range header.Formats {
		_, _ = out.WriteString("##FORMAT=")
		_ = FormatFormatInformation(out, format, false)
	}
	for key, metas := range header.Meta {
		for _, meta := range metas {
			_, _ = out.WriteString("##")
			_, _ = out.WriteString(key)
			_ = out.WriteByte('=')
			_ = FormatMetaInformation(out, meta)
		}
	}
	_ = out.WriteByte('#')
	if len(header.Columns) > 0 {
		_, _ = out.WriteString(header.Columns[0])
		for _, col := range header.Columns[1:] {
			_ = out.WriteByte('\t')
			_, _ = out.WriteString(col)
		}
	}
	return out.WriteByte('\n')
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

func formatValue(out []byte, value interface{}) ([]byte, error) {
	switch v := value.(type) {
	case int:
		return strconv.AppendInt(out, int64(v), 10), nil
	case float64:
		/*
			if math.Floor(v) != v {
				return strconv.AppendFloat(out, v, 'f', 2, 64)
			}
		*/
		return strconv.AppendFloat(out, v, 'f', -1, 64), nil
	case rune:
		if v < utf8.RuneSelf {
			return append(out, byte(v)), nil
		}
		pos := len(out)
		out = append(out, '1', '2', '3', '4', '5', '6')
		buf := out[pos:]
		return out[:pos+utf8.EncodeRune(buf, v)], nil
	case string:
		return append(out, v...), nil
	default:
		return nil, errors.New("invalid value type")
	}
}

func formatInfoEntry(out []byte, entry utils.SmallMapEntry) ([]byte, error) {
	out = append(out, (*entry.Key)...)
	switch e := entry.Value.(type) {
	case bool:
		if !e {
			return nil, errors.New("unexpected boolean value")
		}
		return out, nil
	case []interface{}:
		out = append(out, '=')
		if len(e) == 0 {
			return out, nil
		}
		var err error
		out, err = formatValue(out, e[0])
		if err != nil {
			return nil, err
		}
		for _, v := range e[1:] {
			out = append(out, ',')
			out, err = formatValue(out, v)
			if err != nil {
				return nil, err
			}
		}
		return out, nil
	default:
		out = append(out, '=')
		return formatValue(out, entry.Value)
	}
}

func formatInfo(out []byte, info utils.SmallMap) ([]byte, error) {
	if len(info) == 0 {
		return out, nil
	}
	var err error
	out, err = formatInfoEntry(out, info[0])
	if err != nil {
		return nil, err
	}
	for _, entry := range info[1:] {
		out = append(out, ';')
		out, err = formatInfoEntry(out, entry)
		if err != nil {
			return nil, err
		}
	}
	return out, nil
}

func formatGenotypeDataEntry(out []byte, format utils.Symbol, data utils.SmallMap) ([]byte, bool, error) {
	switch value, _ := data.Get(format); val := value.(type) {
	case nil:
		return append(out, '.'), false, nil
	case []interface{}:
		if len(val) == 0 {
			return out, true, nil
		}
		var err error
		if val[0] == nil {
			out = append(out, '.')
		} else {
			out, err = formatValue(out, val[0])
			if err != nil {
				return nil, false, err
			}
		}
		for _, v := range val[1:] {
			out = append(out, ',')
			if v == nil {
				out = append(out, '.')
			} else {
				out, err = formatValue(out, v)
				if err != nil {
					return nil, false, err
				}
			}
		}
		return out, true, nil
	default:
		var err error
		out, err = formatValue(out, value)
		if err != nil {
			return nil, false, err
		}
		return out, true, nil
	}
}

func formatGenotypeData(out []byte, format []utils.Symbol, data utils.SmallMap) ([]byte, error) {
	if len(format) == 0 {
		return out, nil
	}
	pos := len(out)
	out, ok, err := formatGenotypeDataEntry(out, format[0], data)
	if err != nil {
		return nil, err
	}
	if ok {
		pos = len(out)
	}
	for _, f := range format[1:] {
		out = append(out, ':')
		out, ok, err = formatGenotypeDataEntry(out, f, data)
		if err != nil {
			return nil, err
		}
		if ok {
			pos = len(out)
		}
	}
	if format[len(format)-1] == GT {
		return out, nil
	}
	return out[:pos], nil
}

// Format outputs a VCF variant line
func (variant *Variant) Format(out []byte) ([]byte, error) {
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
		out = append(strconv.AppendFloat(out, value, 'f', -1, 64), '\t')
	} else {
		out = append(out, '.', '\t')
	}
	if len(variant.Filter) == 0 {
		out = append(out, '.', '\t')
	} else {
		out = append(formatSymbolList(out, variant.Filter, ';'), '\t')
	}
	var err error
	out, err = formatInfo(out, variant.Info)
	if err != nil {
		return nil, err
	}
	if len(variant.GenotypeFormat) > 0 {
		out = append(out, '\t')
		out = formatSymbolList(out, variant.GenotypeFormat, ':')
		for _, data := range variant.GenotypeData {
			out = append(out, '\t')
			out, err = formatGenotypeData(out, variant.GenotypeFormat, data)
			if err != nil {
				return nil, err
			}
		}
	}
	return append(out, '\n'), nil
}

// Format outputs a full VCF struct
func (vcf *Vcf) Format(out *bufio.Writer) error {
	if err := vcf.Header.Format(out); err != nil {
		return err
	}
	var buf []byte
	var err error
	for _, variant := range vcf.Variants {
		if buf, err = variant.Format(buf); err != nil {
			return err
		}
		if _, err = out.Write(buf); err != nil {
			return err
		}
		buf = buf[:0]
	}
	return nil
}

// The possible file extensions for VCF or BCF files, or gz-compressed VCF files
const (
	VcfExt = ".vcf"
	BcfExt = ".bcf"
	GzExt  = ".gz"
)

// InputFile represents a VCF or BCF file for input.
type InputFile struct {
	rc io.ReadCloser
	*bufio.Reader
	*exec.Cmd
}

// OutputFile represents a VCF or BCF file for output.
type OutputFile struct {
	wc io.WriteCloser
	*bufio.Writer
	*exec.Cmd
}

// Reader is a bufio.Reader for a VCF or BCF InputFile.
type Reader bufio.Reader

// Writer is a bufio.Writer for a VCF or BCF OutputFile.
type Writer bufio.Writer

// VcfReader returns the reader for a VCF or BCF InputFile.
func (input *InputFile) VcfReader() *Reader {
	return (*Reader)(input.Reader)
}

// VcfWriter returns the Writer for a VCF or BCF OutputFile.
func (output *OutputFile) VcfWriter() *Writer {
	return (*Writer)(output.Writer)
}

// Open a VCF file for input.
//
// If the filename extension is .bcf or .gz, use bcftools view for
// input. Tell bcftools view to only return the header section for
// input when headerOnly is true.
//
// bcftools must be visible in the directories named by the PATH
// environment variable for .bcf or .gz input.
//
// If the filename extension is not .bcf or .gz, then .vcf is always
// assumed.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func Open(name string, headerOnly bool) (*InputFile, error) {
	switch filepath.Ext(name) {
	case BcfExt, GzExt:
		if _, err := os.Stat(name); os.IsNotExist(err) {
			return nil, err
		}
		args := []string{"view"}
		if headerOnly {
			args = append(args, "-h")
		}
		args = append(args, []string{"--threads", strconv.FormatInt(int64(runtime.GOMAXPROCS(0)), 10)}...)
		args = append(args, name)
		cmd := exec.Command("bcftools", args...)
		outPipe, err := cmd.StdoutPipe()
		if err != nil {
			return nil, err
		}
		err = cmd.Start()
		if err != nil {
			return nil, err
		}
		return &InputFile{outPipe, bufio.NewReader(outPipe), cmd}, nil
	default:
		if name == "/dev/stdin" {
			return &InputFile{os.Stdin, bufio.NewReader(os.Stdin), nil}, nil
		}
		file, err := os.Open(name)
		if err != nil {
			return nil, err
		}
		return &InputFile{file, bufio.NewReader(file), nil}, nil
	}
}

// Create a VCF file for output.
//
// If the filename extension is .bcf or .gz, use bcftools view for
// output.
//
// bcftools must be visible in the directories named by the PATH
// environment variable for .bcf or .gz output.
//
// If the filename extension is not .bcf or .gz, then .vcf is always
// assumed.
//
// If the name is "/dev/stdout", then the output is written to
// os.Stdout.
func Create(name string, compressed bool) (*OutputFile, error) {
	ext := filepath.Ext(name)
	if ext == BcfExt || ext == GzExt || compressed {
		args := []string{"view"}
		switch ext {
		case BcfExt:
			if compressed {
				args = append(args, "-Ob")
			} else {
				args = append(args, "-Ou")
			}
		case GzExt:
			args = append(args, "-Oz")
		default:
			if compressed {
				args = append(args, "-Oz")
			} else {
				args = append(args, "-Ov")
			}
		}
		args = append(args, []string{"--threads", strconv.FormatInt(int64(runtime.GOMAXPROCS(0)), 10)}...)
		args = append(args, []string{"-o", name, "-"}...)
		cmd := exec.Command("bcftools", args...)
		inPipe, err := cmd.StdinPipe()
		if err != nil {
			return nil, err
		}
		err = cmd.Start()
		if err != nil {
			return nil, err
		}
		return &OutputFile{inPipe, bufio.NewWriter(inPipe), cmd}, nil
	}
	if name == "/dev/stdout" {
		return &OutputFile{os.Stdout, bufio.NewWriter(os.Stdout), nil}, nil
	}
	file, err := os.Create(name)
	if err != nil {
		return nil, err
	}
	return &OutputFile{file, bufio.NewWriter(file), nil}, nil
}

// Close the VCF input file. If bcftools view is used for input, wait
// for its process to finish.
func (input *InputFile) Close() error {
	if input.rc != os.Stdin {
		err := input.rc.Close()
		if err != nil {
			return err
		}
	}
	if input.Cmd != nil {
		err := input.Wait()
		if err != nil {
			return err
		}
	}
	return nil
}

// Close the VCF input file. If bcftools view is used for input, wait
// for its process to finish.
func (output *OutputFile) Close() error {
	err := output.Flush()
	if err != nil {
		return err
	}
	if output.wc != os.Stdout {
		err := output.wc.Close()
		if err != nil {
			return err
		}
	}
	if output.Cmd != nil {
		err := output.Wait()
		if err != nil {
			return err
		}
	}
	return nil
}
