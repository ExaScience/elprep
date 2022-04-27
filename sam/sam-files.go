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

package sam

import (
	"bufio"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/exascience/elprep/v5/utils/nibbles"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/utils"
	"github.com/exascience/pargo/pipeline"
)

// parseSamHeaderField parses a field in a header line in a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
func (sc *stringScanner) parseSamHeaderField() (tag, value string) {
	tag, ok := sc.readUntil(':')
	if !ok || len(tag) != 2 {
		log.Panicf("invalid field tag %v", tag)
	}
	value, _ = sc.readUntil('\t')
	return tag, value
}

// parseSamHeaderLine parses a header line in a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
//
// The @ record type code must have already been
// scanned. parseSamHeaderLine cannot be used for @CO lines.
func (sc *stringScanner) parseSamHeaderLine() utils.StringMap {
	record := make(utils.StringMap)
	for sc.len() > 0 {
		tag, value := sc.parseSamHeaderField()
		if !record.SetUniqueEntry(tag, value) {
			log.Panicf("duplicate field tag %v in a SAM header line", tag)
		}
	}
	return record
}

// ParseSamHeader parses a complete header in a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
//
// Returns a freshly allocated header and a non-nil error value
// if an error occurred during parsing.
func ParseSamHeader(reader *bufio.Reader) (hdr *Header) {
	hdr = NewHeader()
	var sc stringScanner
	for first := true; ; first = false {
		switch data, err := reader.Peek(1); {
		case err == io.EOF:
			return hdr
		case err != nil:
			log.Panic(err)
		case data[0] != '@':
			return hdr
		}
		bytes, err := reader.ReadSlice('\n')
		length := len(bytes)
		switch {
		case err == nil:
			length--
		case err != io.EOF:
			log.Panic(err)
		}
		line := bytes[4:length]
		sc.reset(line)
		switch string(bytes[0:4]) {
		case "@HD\t":
			if !first {
				log.Panic("@HD line not in first line when parsing a SAM header")
			}
			hdr.HD = sc.parseSamHeaderLine()
		case "@SQ\t":
			hdr.SQ = append(hdr.SQ, sc.parseSamHeaderLine())
		case "@RG\t":
			hdr.RG = append(hdr.RG, sc.parseSamHeaderLine())
		case "@PG\t":
			hdr.PG = append(hdr.PG, sc.parseSamHeaderLine())
		case "@CO\t":
			hdr.CO = append(hdr.CO, string(line))
		default:
			switch code := string(bytes[0:3]); {
			case code == "@CO":
				hdr.CO = append(hdr.CO, string(bytes[3:]))
			case IsHeaderUserTag(code):
				if bytes[3] != '\t' {
					log.Panicf("header code %v not followed by a tab when parsing a SAM header", code)
				}
				hdr.AddUserRecord(code, sc.parseSamHeaderLine())
			default:
				log.Panicf("unknown SAM record type code %v", code)
			}
		}
	}
}

// SkipSamHeader skips the complete header in a SAM file. This is more
// efficient than calling ParseHeader and ignoring its result.
//
// Returns a non-nil error value if an error occurred.
func SkipSamHeader(reader *bufio.Reader) {
	for {
		data, err := reader.Peek(1)
		if err != nil {
			if err == io.EOF {
				return
			}
			log.Panic(err)
		}
		if data[0] != '@' {
			break
		}
		for {
			byte, err := reader.ReadByte()
			if err != nil {
				if err == io.EOF {
					return
				}
				log.Panic(err)
			}
			if byte == '\n' {
				break
			}
		}
	}
}

func splitHeaderField(field string) (tag, value string) {
	if field[2] != ':' {
		log.Panicf("incorrectly formatted SAM file field %v", field)
	}
	return field[:2], field[3:]
}

// ParseHeaderLineFromString parses a SAM header line from a string,
// except that entries are separated by white space instead of
// tabulators. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3.
//
// The @ record type code must have already been
// scanned. ParseHeaderLineFromString cannot be used for @CO lines.
func ParseHeaderLineFromString(line string) utils.StringMap {
	record := make(utils.StringMap)
	fields := strings.Fields(line)
	for _, field := range fields {
		tag, value := splitHeaderField(field)
		if !record.SetUniqueEntry(tag, value) {
			log.Panicf("duplicate field tag %v in a SAM header line", tag)
		}
	}
	return record
}

// samFieldParser is the signature for all parsers for optional fields in
// read alignment lines in SAM files.
type samFieldParser func(*stringScanner, utils.Symbol) (utils.Symbol, interface{})

// parseSamChar parses a single tab-delimited character and returns it as
// a byte. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 1.5.
func (sc *stringScanner) parseSamChar(tag utils.Symbol) (utils.Symbol, interface{}) {
	value, _ := sc.readByteUntil('\t')
	return tag, value
}

// parseSamInteger parses a single tab-delimited integer and returns it
// as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.5.
func (sc *stringScanner) parseSamInteger(tag utils.Symbol) (utils.Symbol, interface{}) {
	value, _ := sc.readUntil('\t')
	return tag, internal.ParseInt(value, 10, 64)
}

// parseSamFloat parses a single tab-delimited float and returns it as a
// float32. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.5.
func (sc *stringScanner) parseSamFloat(tag utils.Symbol) (utils.Symbol, interface{}) {
	value, _ := sc.readUntil('\t')
	return tag, float32(internal.ParseFloat(value, 32))
}

// parseSamString parses a single tab-delimited string and returns
// it. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 1.5.
func (sc *stringScanner) parseSamString(tag utils.Symbol) (utils.Symbol, interface{}) {
	value, _ := sc.readUntil('\t')
	switch tag {
	case CC, LB, PG, PU, RG:
		return tag, *utils.Intern(value)
	default:
		return tag, value
	}
}

// parseSamByteArray parses a byte array in the tab-delimited Hex format
// and returns it as a ByteArray. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
func (sc *stringScanner) parseSamByteArray(tag utils.Symbol) (utils.Symbol, interface{}) {
	value, _ := sc.readUntil('\t')
	result := ByteArray(make([]byte, 0, len(value)>>1))
	for i := 0; i < len(value); i += 2 {
		result = append(result, byte(internal.ParseUint(value[i:i+2], 16, 8)))
	}
	return tag, result
}

// parseSamNumericArray parses a typed, tab-delimited, and
// comma-separated integer or numeric array and returns it as a
// []int8, []uint8, []int16, []uint16, []int32, []uint32, or
// []float32. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.5.
func (sc *stringScanner) parseSamNumericArray(tag utils.Symbol) (utils.Symbol, interface{}) {
	ntype, ok := sc.readByteUntil(',')
	if !ok {
		log.Panic("missing entry in numeric array")
	}
	switch ntype {
	case 'c':
		var result []int8
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, int8(internal.ParseInt(entry, 10, 8)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'C':
		var result []uint8
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, uint8(internal.ParseUint(entry, 10, 8)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 's':
		var result []int16
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, int16(internal.ParseUint(entry, 10, 16)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'S':
		var result []uint16
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, uint16(internal.ParseUint(entry, 10, 16)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'i':
		var result []int32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, int32(internal.ParseInt(entry, 10, 32)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'I':
		var result []uint32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, uint32(internal.ParseUint(entry, 10, 32)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'f':
		var result []float32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			result = append(result, float32(internal.ParseFloat(entry, 32)))
			if sep != ',' {
				break
			}
		}
		return tag, result
	default:
		log.Panicf("invalid numeric array type %v", ntype)
		return nil, nil
	}
}

var optionalFieldParseTable = map[byte]samFieldParser{
	'A': (*stringScanner).parseSamChar,
	'i': (*stringScanner).parseSamInteger,
	'f': (*stringScanner).parseSamFloat,
	'Z': (*stringScanner).parseSamString,
	'H': (*stringScanner).parseSamByteArray,
	'B': (*stringScanner).parseSamNumericArray,
}

// parseSamOptionalField parses a single tab-delimited optional field in
// a SAM read alignment line and returns it as a tag/value pair. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
//
// The second return value is one of byte (representing an ASCII
// character), int32, float32, string, ByteArray, []int8, []uint8,
// []int16, []uint16, []int32, []uint32, or []float32.
func (sc *stringScanner) parseSamOptionalField() (tag utils.Symbol, value interface{}) {
	tagname, ok := sc.readUntil(':')
	if !ok || len(tagname) != 2 {
		log.Panicf("invalid field tag %v in SAM alignment line", tagname)
	}
	tag = utils.Intern(tagname)
	typebyte, ok := sc.readByteUntil(':')
	if !ok {
		log.Panicf("invalid field type %v in SAM alignment line", typebyte)
	}
	return optionalFieldParseTable[typebyte](sc, tag)
}

func (sc *stringScanner) doString() string {
	value, ok := sc.readUntil('\t')
	if !ok {
		log.Panic("missing tabulator in SAM alignment line")
	}
	return value
}

func (sc *stringScanner) doSeq() (seq Sequence) {
	var n nibbles.Nibbles
	for end := sc.index; end < len(sc.data); end++ {
		ch := sc.data[end]
		if ch == '\t' {
			sc.index = end + 1
			return Sequence(n)
		}
		nibble, ok := baseToNibble[ch]
		if !ok {
			nibble = 15
		}
		n = n.Append(nibble)
	}
	log.Panic("missing tabulator in SAM alignment line")
	return
}

func (sc *stringScanner) doInt32() int32 {
	return int32(internal.ParseInt(sc.doString(), 10, 32))
}

func (sc *stringScanner) doUint(bitSize int) uint64 {
	return internal.ParseUint(sc.doString(), 10, bitSize)
}

// parseSamAlignment parses a read alignment line in a SAM file and
// returns a freshly allocated alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and
// 1.5.
func (sc *stringScanner) parseSamAlignment() *Alignment {
	aln := new(Alignment)

	aln.QNAME = sc.doString()
	aln.FLAG = uint16(sc.doUint(16))
	aln.RNAME = *utils.Intern(sc.doString())
	aln.POS = sc.doInt32()
	aln.MAPQ = byte(sc.doUint(8))
	cigarString := sc.doString()
	aln.CIGAR = ScanCigarString(cigarString)
	aln.RNEXT = *utils.Intern(sc.doString())
	aln.PNEXT = sc.doInt32()
	aln.TLEN = sc.doInt32()
	aln.SEQ = sc.doSeq()
	aln.QUAL = sc.readBytes()
	for i := range aln.QUAL {
		aln.QUAL[i] -= 33
	}

	for sc.len() > 0 {
		aln.TAGS.Set(sc.parseSamOptionalField())
	}

	return aln
}

// formatSamString writes a SAM file TAG of type string.
func formatSamString(out []byte, tag, value string) []byte {
	out = append(out, '\t')
	out = append(out, tag...)
	out = append(out, ':')
	out = append(out, value...)
	return out
}

// formatSamHeaderLine writes a header line in a SAM file header
// section. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3.
func formatSamHeaderLine(out []byte, code string, record utils.StringMap) []byte {
	out = append(out, code...)
	// sort keys so we always get the same output
	keys := make([]string, 0, len(record))
	for k := range record {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, key := range keys {
		value := record[key]
		out = formatSamString(out, key, value)
	}
	out = append(out, '\n')
	return out
}

// formatSamComment writes a header comment line in a SAM file header
// section. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3.
func formatSamComment(out []byte, comment string) []byte {
	out = append(out, "@CO\t"...)
	out = append(out, comment...)
	out = append(out, '\n')
	return out
}

// FormatSam writes the header section of a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
func (hdr *Header) FormatSam(out []byte) []byte {
	if hdr.HD != nil {
		out = formatSamHeaderLine(out, "@HD", hdr.HD)
	}
	for _, record := range hdr.SQ {
		out = formatSamHeaderLine(out, "@SQ", record)
	}
	for _, record := range hdr.RG {
		out = formatSamHeaderLine(out, "@RG", record)
	}
	for _, record := range hdr.PG {
		out = formatSamHeaderLine(out, "@PG", record)
	}
	for _, comment := range hdr.CO {
		out = formatSamComment(out, comment)
	}
	for code, records := range hdr.UserRecords {
		for _, record := range records {
			out = formatSamHeaderLine(out, code, record)
		}
	}
	return out
}

// formatSamTag writes a SAM file TAG by appending its ASCII-string
// representation to out and returning the result, dispatching on the
// actual type of the given value. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
//
// The following types are accepted: byte (A), int64 (i), float32 (f),
// string (Z), ByteArray (H), []int8 (B:c), []uint8 (B:C), []int16
// (B:s), []uint16 (B:S), []int32 (B:i), []uint32 (B:I), and []float32
// (B:f).
func formatSamTag(out []byte, tag utils.Symbol, value interface{}) []byte {
	out = append(out, '\t')
	out = append(out, *tag...)

	switch val := value.(type) {
	case byte:
		out = append(append(out, ":A:"...), val)
	case int64:
		out = strconv.AppendInt(append(out, ":i:"...), val, 10)
	case float32:
		out = strconv.AppendFloat(append(out, ":f:"...), float64(val), 'g', -1, 32)
	case string:
		out = append(append(out, ":Z:"...), val...)
	case ByteArray:
		out = append(out, ":H:"...)
		for _, b := range val {
			if b < 16 {
				out = append(out, '0')
			}
			out = strconv.AppendUint(out, uint64(b), 16)
		}
	case []int8:
		out = append(out, ":B:c"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint8:
		out = append(out, ":B:C"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []int16:
		out = append(out, ":B:s"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint16:
		out = append(out, ":B:S"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []int32:
		out = append(out, ":B:i"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint32:
		out = append(out, ":B:I"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []float32:
		out = append(out, ":B:f"...)
		for _, v := range val {
			out = strconv.AppendFloat(append(out, ','), float64(v), 'g', -1, 32)
		}
	default:
		log.Panicf("unknown SAM alignment TAG type %v", value)
	}

	return out
}

func cigarToString(b []byte, cigar []CigarOperation) []byte {
	if len(cigar) == 0 {
		return append(b, '*')
	}
	for _, op := range cigar {
		b = strconv.AppendInt(b, int64(op.Length), 10)
		b = append(b, op.Operation)
	}
	return b
}

// FormatAlignment writes a SAM file read alignment line by appending its
// ASCII-string representation to out and returning the result. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and
// 1.5.
func FormatAlignment(aln *Alignment, out []byte) []byte {
	out = append(append(out, aln.QNAME...), '\t')
	out = append(strconv.AppendUint(out, uint64(aln.FLAG), 10), '\t')
	out = append(append(out, aln.RNAME...), '\t')
	out = append(strconv.AppendInt(out, int64(aln.POS), 10), '\t')
	out = append(strconv.AppendUint(out, uint64(aln.MAPQ), 10), '\t')
	out = append(cigarToString(out, aln.CIGAR), '\t')
	switch aln.RNEXT {
	case "=":
		out = append(out, "=\t"...)
	case "*":
		out = append(out, "*\t"...)
	default:
		if aln.RNEXT == aln.RNAME {
			out = append(out, "=\t"...)
		} else {
			out = append(append(out, aln.RNEXT...), '\t')
		}
	}
	out = append(strconv.AppendInt(out, int64(aln.PNEXT), 10), '\t')
	out = append(strconv.AppendInt(out, int64(aln.TLEN), 10), '\t')
	l := aln.SEQ.Len()
	for i := 0; i < l; i++ {
		out = append(out, aln.SEQ.Base(i))
	}
	out = append(out, '\t')
	for _, qual := range aln.QUAL {
		out = append(out, qual+33)
	}

	for _, entry := range aln.TAGS {
		out = formatSamTag(out, entry.Key, entry.Value)
	}

	return append(out, '\n')
}

// samReader is an AlignmentFileReader for a SAM InputFile.
type samReader struct {
	rc  io.Closer
	buf *bufio.Reader
	*pipeline.BytesScanner
}

// Close the SAM input file.
func (reader *samReader) Close() {
	if reader.rc != os.Stdin {
		internal.Close(reader.rc)
	}
}

const maxTokenSize = 16 * 64 * 1024

// ParseHeader implements the method of the AlignmentFileReader interface.
func (reader *samReader) ParseHeader() (hdr *Header) {
	hdr = ParseSamHeader(reader.buf)
	reader.BytesScanner = pipeline.NewBytesScanner(reader.buf)
	reader.BytesScanner.Buffer(nil, maxTokenSize)
	return
}

// SkipHeader implements the method of the AlignmentFileReader interface.
func (reader *samReader) SkipHeader() {
	SkipSamHeader(reader.buf)
	reader.BytesScanner = pipeline.NewBytesScanner(reader.buf)
	reader.BytesScanner.Buffer(nil, maxTokenSize)
}

// ParseAlignmentFromBytes parses a SAM alignment from its string representation
func ParseAlignmentFromBytes(s []byte) *Alignment {
	var sc stringScanner
	sc.reset(s)
	return sc.parseSamAlignment()
}

// ParseAlignment implements the method of the AlignmentFileReader interface.
func (*samReader) ParseAlignment(record []byte) *Alignment {
	return ParseAlignmentFromBytes(record)
}

// ParseAlignment parses a SAM alignment from its string representation
func ParseAlignment(s string) *Alignment {
	return ParseAlignmentFromBytes([]byte(s))
}

// samWriter is an AlignmentFileWriter for a SAM OutputFile.
type samWriter struct {
	wc io.WriteCloser
}

// Close the SAM output file.
func (writer *samWriter) Close() {
	if writer.wc != os.Stdout {
		internal.Close(writer.wc)
	}
}

// FormatHeader implements the method of the AlignmentFileWriter interface.
func (writer *samWriter) FormatHeader(hdr *Header) {
	internal.Write(writer.wc, hdr.FormatSam(nil))
}

// FormatAlignment implements the method of the AlignmentFileWriter interface.
func (*samWriter) FormatAlignment(aln *Alignment, out []byte) []byte {
	return FormatAlignment(aln, out)
}

// Write implements the method of the io.Writer interface.
func (writer *samWriter) Write(p []byte) int {
	n, err := writer.wc.Write(p)
	if err != nil {
		log.Panic(err)
	}
	return n
}
