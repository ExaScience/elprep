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
	"bytes"
	"context"
	"encoding/binary"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/utils"
	"github.com/exascience/elprep/v5/utils/bgzf"
	"github.com/exascience/elprep/v5/utils/nibbles"
)

// BAMReference is a an entry in a slice of BAM-encoded sequence dictionary entries.
// See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
type BAMReference struct {
	Name   string
	Length int32
}

func parseBamHeaderReferences(reader io.Reader, text []byte) (references []BAMReference) {
	var nRef int32
	internal.BinaryRead(reader, &nRef)
	for i := int32(0); i < nRef; i++ {
		var lName int32
		internal.BinaryRead(reader, &lName)
		for cap(text) < int(lName) {
			text = append(text[:cap(text)], 0)
		}
		text = text[:int(lName)]
		internal.ReadFull(reader, text)
		var lRef int32
		internal.BinaryRead(reader, &lRef)
		references = append(references, BAMReference{
			Name:   *utils.Intern(string(text[:len(text)-1])),
			Length: lRef,
		})
	}
	return
}

// bamMagic is the magic string for the BAM format. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
const bamMagic = "BAM\x01"

// ParseBamHeader parses a complete header in a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
//
// Returns a freshly allocated header, the BAM-encoded sequence dictionary,
// and a non-nil error value if an error occurred during parsing.
func ParseBamHeader(reader io.Reader) (*Header, []BAMReference) {
	text := make([]byte, 4)
	internal.ReadFull(reader, text)
	if string(text) != bamMagic {
		log.Panic("invalid BAM file header")
	}
	var lText int32
	internal.BinaryRead(reader, &lText)
	for cap(text) < int(lText) {
		text = append(text[:cap(text)], 0)
	}
	text = text[:int(lText)]
	internal.ReadFull(reader, text)
	for i, b := range text {
		if b == 0 {
			text = text[:i]
			break
		}
	}
	return ParseSamHeader(bufio.NewReader(bytes.NewReader(text))), parseBamHeaderReferences(reader, text)
}

// SkipBamHeader skips the complete header in a BAM file. This is more
// efficient than calling ParseBamHeader and ignoring its result.
//
// Returns the BAM-encoded sequence dictionary and a non-nil error value
// if an error occurred during parsing.
func SkipBamHeader(reader io.Reader) []BAMReference {
	text := make([]byte, 4)
	internal.ReadFull(reader, text)
	if string(text) != bamMagic {
		log.Panic("invalid BAM file header")
	}
	var lText int32
	internal.BinaryRead(reader, &lText)
	var err error
	if s, ok := reader.(io.Seeker); ok {
		_, err = s.Seek(int64(lText), io.SeekCurrent)
	} else {
		_, err = io.CopyN(ioutil.Discard, reader, int64(lText))
	}
	if err != nil {
		log.Panic(err)
	}
	return parseBamHeaderReferences(reader, text)
}

// bamFieldParser is the signature for all parsers for optional fields in
// read alignment records in BAM files.
type bamFieldParser func(record []byte, index int) (value interface{}, newIndex int)

// parseBamChar parses an A optional field in a BAM alignment record and returns
// it as a byte. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamChar(record []byte, index int) (value interface{}, newIndex int) {
	return record[index], index + 1
}

// parseBamI8 parses a c optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamI8(record []byte, index int) (value interface{}, newIndex int) {
	return int64(int8(record[index])), index + 1
}

// parseBamU8 parses a C optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamU8(record []byte, index int) (value interface{}, newIndex int) {
	return int64(record[index]), index + 1
}

// parseBamI16 parses an s optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamI16(record []byte, index int) (value interface{}, newIndex int) {
	return int64(int16(binary.LittleEndian.Uint16(record[index : index+2]))), index + 2
}

// parseBamU16 parses an S optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamU16(record []byte, index int) (value interface{}, newIndex int) {
	return int64(binary.LittleEndian.Uint16(record[index : index+2])), index + 2
}

// parseBamI32 parses an i optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamI32(record []byte, index int) (value interface{}, newIndex int) {
	return int64(int32(binary.LittleEndian.Uint32(record[index : index+4]))), index + 4
}

// parseBamU32 parses an I optional field in a BAM alignment record and returns
// it as an int64. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamU32(record []byte, index int) (value interface{}, newIndex int) {
	return int64(binary.LittleEndian.Uint32(record[index : index+4])), index + 4
}

// parseBamFloat parses an f optional field in a BAM alignment record and returns
// it as an float32. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamFloat(record []byte, index int) (value interface{}, newIndex int) {
	return math.Float32frombits(binary.LittleEndian.Uint32(record[index : index+4])), index + 4
}

// parseBamString parses a Z optional field in a BAM alignment record and returns
// it as a string. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamString(record []byte, index int) (value interface{}, newIndex int) {
	for end := index; end < len(record); end++ {
		if record[end] == 0 {
			value = string(record[index:end])
			newIndex = end + 1
			return
		}
	}
	log.Panic("missing NUL byte in an optional string field in a BAM alignment record")
	return nil, -1
}

// parseBamByteArray parses an H optional field in a BAM alignment record and returns
// it as a ByteArray. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamByteArray(record []byte, index int) (value interface{}, newIndex int) {
	for end := index; end < len(record); end++ {
		if record[end] == '0' {
			result := ByteArray(make([]byte, 0, (end-index)>>1))
			for i := index; i < end; i += 2 {
				result = append(result, byte(internal.ParseUint(string(record[i:i+2]), 16, 8)))
			}
			return result, end + 1
		}
	}
	log.Panic("missing NUL byte in an optional hex-formatted string field in a BAM alignment record")
	return nil, -1
}

// parseBamNumericArray parses a B optional field in a BAM alignment record and
// returns it as a []int8, []uint8, []int16, []uint16, []int32, []uint32, or
// []float32. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 4.2.4.
func parseBamNumericArray(record []byte, index int) (value interface{}, newIndex int) {
	ntype := record[index]
	index++
	count := int(int32(binary.LittleEndian.Uint32(record[index : index+4])))
	index += 4
	switch ntype {
	case 'c':
		result := make([]int8, count)
		for i := 0; i < count; i++ {
			result[i] = int8(record[index+i])
		}
		return result, index + count
	case 'C':
		result := make([]uint8, count)
		copy(result, record[index:index+count])
		return result, index + count
	case 's':
		result := make([]int16, count)
		for i, j := 0, 0; i < count; i, j = i+1, j+2 {
			result[i] = int16(binary.LittleEndian.Uint16(record[index+j : index+j+2]))
		}
		return result, index + (count << 1)
	case 'S':
		result := make([]uint16, count)
		for i, j := 0, 0; i < count; i, j = i+1, j+2 {
			result[i] = binary.LittleEndian.Uint16(record[index+j : index+j+2])
		}
		return result, index + (count << 1)
	case 'i':
		result := make([]int32, count)
		for i, j := 0, 0; i < count; i, j = i+1, j+4 {
			result[i] = int32(binary.LittleEndian.Uint32(record[index+j : index+j+4]))
		}
		return result, index + (count << 2)
	case 'I':
		result := make([]uint32, count)
		for i, j := 0, 0; i < count; i, j = i+1, j+4 {
			result[i] = binary.LittleEndian.Uint32(record[index+j : index+j+4])
		}
		return result, index + (count << 2)
	case 'f':
		result := make([]float32, count)
		for i, j := 0, 0; i < count; i, j = i+1, j+4 {
			result[i] = math.Float32frombits(binary.LittleEndian.Uint32(record[index+j : index+j+4]))
		}
		return result, index + (count << 2)
	default:
		log.Panic("invalid subtype in a numeric array in a BAM alignment record")
		return nil, -1
	}
}

var optionalBAMFieldParseTable = map[byte]bamFieldParser{
	'A': parseBamChar,
	'c': parseBamI8,
	'C': parseBamU8,
	's': parseBamI16,
	'S': parseBamU16,
	'i': parseBamI32,
	'I': parseBamU32,
	'f': parseBamFloat,
	'Z': parseBamString,
	'H': parseBamByteArray,
	'B': parseBamNumericArray,
}

var (
	star     = *utils.Intern("*")
	eq       = *utils.Intern("=")
	cg       = utils.Intern("CG")
	cigarOps = []byte("MIDNSHP=X")
	cigarMap = make(map[byte]byte)
)

func init() {
	for i, b := range cigarOps {
		cigarMap[b] = byte(i)
	}
}

const (
	refIDIndex     = 0
	posIndex       = 4
	lReadNameIndex = posIndex + 4
	mapqIndex      = lReadNameIndex + 1
	binIndex       = mapqIndex + 1
	nCigarOpIndex  = binIndex + 2
	flagIndex      = nCigarOpIndex + 2
	lSeqIndex      = flagIndex + 2
	nextRefIDIndex = lSeqIndex + 4
	nextPosIndex   = nextRefIDIndex + 4
	tlenIndex      = nextPosIndex + 4
	readNameIndex  = tlenIndex + 4
)

// parseBamAlignment parses a read alignment record in a BAM file and
// returns a freshly allocated alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 4.2.
func parseBamAlignment(record []byte, references []BAMReference) *Alignment {
	aln := new(Alignment)

	refID := int32(binary.LittleEndian.Uint32(record[refIDIndex : refIDIndex+4]))
	if refID < 0 {
		aln.RNAME = star
	} else {
		aln.RNAME = references[int(refID)].Name
	}

	aln.POS = int32(binary.LittleEndian.Uint32(record[posIndex:posIndex+4])) + 1

	lReadName := int(record[lReadNameIndex])

	aln.MAPQ = record[mapqIndex]

	nCigarOp := binary.LittleEndian.Uint16(record[nCigarOpIndex : nCigarOpIndex+2])

	aln.FLAG = binary.LittleEndian.Uint16(record[flagIndex : flagIndex+2])

	lSeq := int32(binary.LittleEndian.Uint32(record[lSeqIndex : lSeqIndex+4]))

	nextRefID := int32(binary.LittleEndian.Uint32(record[nextRefIDIndex : nextRefIDIndex+4]))
	if nextRefID < 0 {
		aln.RNEXT = star
	} else {
		aln.RNEXT = references[int(nextRefID)].Name
		if aln.RNEXT == aln.RNAME {
			aln.RNEXT = eq
		}
	}

	aln.PNEXT = int32(binary.LittleEndian.Uint32(record[nextPosIndex:nextPosIndex+4])) + 1

	aln.TLEN = int32(binary.LittleEndian.Uint32(record[tlenIndex : tlenIndex+4]))

	aln.QNAME = string(record[readNameIndex : readNameIndex+lReadName-1])

	index := readNameIndex + lReadName

	aln.CIGAR = make([]CigarOperation, nCigarOp)

	for i := uint16(0); i < nCigarOp; i, index = i+1, index+4 {
		cigar := binary.LittleEndian.Uint32(record[index : index+4])
		aln.CIGAR[i] = CigarOperation{
			Length:    int32(cigar >> 4),
			Operation: cigarOps[int(0xF&cigar)],
		}
	}

	nextIndex := index + ((int(lSeq) + 1) >> 1)
	aln.SEQ = Sequence(nibbles.ReflectMake(int(lSeq), 0, append([]byte(nil), record[index:nextIndex]...)))
	index = nextIndex

	aln.QUAL = append([]byte(nil), record[index:index+int(lSeq)]...)
	index += int(lSeq)

	for index < len(record) {
		tag := utils.Intern(string(record[index : index+2]))
		typebyte := record[index+2]
		index += 3
		value, newIndex := optionalBAMFieldParseTable[typebyte](record, index)
		if tag == cg && typebyte == 'B' {
			if cigars, ok := value.([]uint32); ok {
				if len(aln.CIGAR) > 0 {
					if op := aln.CIGAR[0]; op.Operation == 'S' && int(op.Length) == aln.SEQ.Len() {
						aln.CIGAR = aln.CIGAR[:0]
						for _, cigar := range cigars {
							aln.CIGAR = append(aln.CIGAR, CigarOperation{
								Length:    int32(cigar >> 4),
								Operation: cigarOps[int(0xF&cigar)],
							})
						}
					}
				}
				continue
			}
		}
		aln.TAGS = append(aln.TAGS, utils.SmallMapEntry{Key: tag, Value: value})
		index = newIndex
	}

	return aln
}

func enlarge(out []byte, by int) (int, []byte) {
	index := len(out)
	length := index + by
	for cap(out) < length {
		out = append(out[:cap(out)], 0)
	}
	out = out[:length]
	return index, out
}

// FormatBam writes the header section of a BAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
func (hdr *Header) FormatBam(out []byte) []byte {
	out = append(out, bamMagic...)
	lTextIndex := len(out)
	out = append(out, "0000"...)

	out = hdr.FormatSam(out)

	binary.LittleEndian.PutUint32(out[lTextIndex:lTextIndex+4], uint32(len(out)-lTextIndex-4))

	var index int
	index, out = enlarge(out, 4)
	binary.LittleEndian.PutUint32(out[index:], uint32(len(hdr.SQ)))

	for _, sq := range hdr.SQ {
		sn := sq["SN"]
		index, out = enlarge(out, 4+len(sn)+1+4)
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(sn)+1))
		index += 4
		copy(out[index:], sn)
		out[index+len(sn)] = 0
		index += len(sn) + 1
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(SQLN(sq)))
	}

	return out
}

var cigarConsumesReferenceBases = map[byte]int32{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}

func (aln *Alignment) bin() uint16 {
	beg := aln.POS - 1
	end := beg
	if !aln.IsUnmapped() {
		for _, op := range aln.CIGAR {
			end += cigarConsumesReferenceBases[op.Operation] * op.Length
		}
		end--
	}
	if beg>>14 == end>>14 {
		return uint16(((1<<15)-1)/7 + (beg >> 14))
	}
	if beg>>17 == end>>17 {
		return uint16(((1<<12)-1)/7 + (beg >> 17))
	}
	if beg>>20 == end>>20 {
		return uint16(((1<<9)-1)/7 + (beg >> 20))
	}
	if beg>>23 == end>>23 {
		return uint16(((1<<6)-1)/7 + (beg >> 23))
	}
	if beg>>26 == end>>26 {
		return uint16(((1<<3)-1)/7 + (beg >> 26))
	}
	return 0
}

const minus1 = 0xFFFFFFFF

// formatBamTag writes a BAM file TAG by appending its binary
// representation to out and returning the result, dispatching on the
// actual type of the given value. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.4.
//
// The following types are accepted: byte (A), int64 (c, C, s, S, i, I),
// float32 (f), string (Z), ByteArray (H), []int8 (B:c), []uint8 (B:C),
// []int16 (B:s), []uint16 (B:S), []int32 (B:i), []uint32 (B:I), and
// []float32 (B:f).
func formatBamTag(out []byte, tag utils.Symbol, value interface{}) []byte {
	var index int

	index, out = enlarge(out, 2)
	copy(out[index:], *tag)

	switch val := value.(type) {
	case byte:
		index, out = enlarge(out, 2)
		out[index] = 'A'
		out[index+1] = val
	case int64:
		if val < 0 {
			if val >= math.MinInt8 {
				index, out = enlarge(out, 2)
				out[index] = 'c'
				out[index+1] = byte(int8(val))
			} else if val >= math.MinInt16 {
				index, out = enlarge(out, 3)
				out[index] = 's'
				binary.LittleEndian.PutUint16(out[index+1:index+3], uint16(val))
			} else if val >= math.MinInt32 {
				index, out = enlarge(out, 5)
				out[index] = 'i'
				binary.LittleEndian.PutUint32(out[index+1:index+5], uint32(val))
			} else {
				log.Panicf("integer value too small in BAM alignment tag %v", value)
			}
		} else {
			if val <= math.MaxUint8 {
				index, out = enlarge(out, 2)
				out[index] = 'C'
				out[index+1] = uint8(val)
			} else if val <= math.MaxUint16 {
				index, out = enlarge(out, 3)
				out[index] = 'S'
				binary.LittleEndian.PutUint16(out[index+1:index+3], uint16(val))
			} else if val <= math.MaxUint32 {
				index, out = enlarge(out, 5)
				out[index] = 'I'
				binary.LittleEndian.PutUint32(out[index+1:index+5], uint32(val))
			} else {
				log.Panicf("integer value too large in BAM alignment tag %v", value)
			}
		}
	case float32:
		index, out = enlarge(out, 5)
		out[index] = 'f'
		binary.LittleEndian.PutUint32(out[index+1:index+5], math.Float32bits(val))
	case string:
		index, out = enlarge(out, 1+len(val)+1)
		out[index] = 'Z'
		index++
		copy(out[index:], val)
		out[index+len(val)] = 0
	case ByteArray:
		index, out = enlarge(out, 1+2*len(val)+1)
		out[index] = 'H'
		index++
		for _, b := range val {
			if hi := b >> 4; hi < 10 {
				out[index] = '0' + hi
			} else {
				out[index] = 'A' - 10 + hi
			}
			if lo := b & 0xF; lo < 10 {
				out[index+1] = '0' + lo
			} else {
				out[index+1] = 'A' - 10 + lo
			}
			index += 2
		}
		out[index] = 0
	case []int8:
		index, out = enlarge(out, 2+4+len(val))
		out[index] = 'B'
		out[index+1] = 'c'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		index += 4
		for _, v := range val {
			out[index] = byte(v)
			index++
		}
	case []uint8:
		index, out = enlarge(out, 2+4+len(val))
		out[index] = 'B'
		out[index+1] = 'C'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		index += 4
		copy(out[index:], val)
	case []int16:
		index, out = enlarge(out, 2+4+2*len(val))
		out[index] = 'B'
		out[index+1] = 's'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		index += 4
		for _, v := range val {
			binary.LittleEndian.PutUint16(out[index:index+2], uint16(v))
			index += 2
		}
	case []uint16:
		index, out = enlarge(out, 2+4+2*len(val))
		out[index] = 'B'
		out[index+1] = 'S'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		index += 4
		for _, v := range val {
			binary.LittleEndian.PutUint16(out[index:index+2], v)
			index += 2
		}
	case []int32:
		index, out = enlarge(out, 2+4+4*len(val))
		out[index] = 'B'
		out[index+1] = 'i'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		for _, v := range val {
			index += 4
			binary.LittleEndian.PutUint32(out[index:index+4], uint32(v))
		}
	case []uint32:
		index, out = enlarge(out, 2+4+4*len(val))
		out[index] = 'B'
		out[index+1] = 'I'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		for _, v := range val {
			index += 4
			binary.LittleEndian.PutUint32(out[index:index+4], v)
		}
	case []float32:
		index, out = enlarge(out, 2+4+4*len(val))
		out[index] = 'B'
		out[index+1] = 'f'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(val)))
		for _, v := range val {
			index += 4
			binary.LittleEndian.PutUint32(out[index:index+4], math.Float32bits(v))
		}
	default:
		log.Panicf("unknown BAM alignment TAG type %v", value)
	}

	return out
}

// formatBamAlignment writes a BAM file read alignment record by appending its
// binary representation to out and returning the result. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
func formatBamAlignment(aln *Alignment, out []byte, dictTable map[string]uint32) []byte {
	var index int

	index, out = enlarge(out, 4)
	blockSizeIndex := index

	index, out = enlarge(out, 4)
	refid, ok := dictTable[aln.RNAME]
	if ok {
		binary.LittleEndian.PutUint32(out[index:], refid)
	} else {
		binary.LittleEndian.PutUint32(out[index:], minus1)
	}

	index, out = enlarge(out, 4)
	binary.LittleEndian.PutUint32(out[index:], uint32(aln.POS-1))

	out = append(out, uint8(len(aln.QNAME)+1))

	out = append(out, aln.MAPQ)

	index, out = enlarge(out, 2)
	binary.LittleEndian.PutUint16(out[index:], aln.bin())

	index, out = enlarge(out, 2)
	if len(aln.CIGAR) <= math.MaxUint16 {
		binary.LittleEndian.PutUint16(out[index:], uint16(len(aln.CIGAR)))
	} else {
		binary.LittleEndian.PutUint16(out[index:], 2)
	}

	index, out = enlarge(out, 2)
	binary.LittleEndian.PutUint16(out[index:], aln.FLAG)

	seqLength := aln.SEQ.Len()
	index, out = enlarge(out, 4)
	binary.LittleEndian.PutUint32(out[index:], uint32(seqLength))

	index, out = enlarge(out, 4)
	if aln.RNEXT != "=" {
		refid, ok = dictTable[aln.RNEXT]
	}
	if ok {
		binary.LittleEndian.PutUint32(out[index:index+4], refid)
	} else {
		binary.LittleEndian.PutUint32(out[index:index+4], minus1)
	}

	index, out = enlarge(out, 4)
	binary.LittleEndian.PutUint32(out[index:], uint32(aln.PNEXT-1))

	index, out = enlarge(out, 4)
	binary.LittleEndian.PutUint32(out[index:], uint32(aln.TLEN))

	index, out = enlarge(out, len(aln.QNAME)+1)
	copy(out[index:], aln.QNAME)
	out[index+len(aln.QNAME)] = 0

	if len(aln.CIGAR) <= math.MaxUint16 {
		index, out = enlarge(out, len(aln.CIGAR)*4)
		for _, op := range aln.CIGAR {
			binary.LittleEndian.PutUint32(out[index:index+4], uint32((op.Length<<4)|int32(cigarMap[op.Operation])))
			index += 4
		}
	} else {
		index, out = enlarge(out, 2*4)
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(int32(aln.SEQ.Len()<<4)|int32(cigarMap['S'])))
		var m int32
		for _, op := range aln.CIGAR {
			m += cigarConsumesReferenceBases[op.Operation] * op.Length
		}
		binary.LittleEndian.PutUint32(out[index+4:index+8], uint32((m<<4)|int32(cigarMap['N'])))
	}

	index, out = enlarge(out, (seqLength+1)>>1)
	nib := nibbles.ReflectMake(seqLength, 0, out[index:])
	nib.Copy(nibbles.Nibbles(aln.SEQ))

	index, out = enlarge(out, len(aln.QUAL))
	copy(out[index:], aln.QUAL)

	for _, entry := range aln.TAGS {
		out = formatBamTag(out, entry.Key, entry.Value)
	}

	if len(aln.CIGAR) > math.MaxUint16 {
		index, out = enlarge(out, 2+2+4+4*len(aln.CIGAR))
		copy(out[index:], *cg)
		index += 2
		out[index] = 'B'
		out[index+1] = 'I'
		index += 2
		binary.LittleEndian.PutUint32(out[index:index+4], uint32(len(aln.CIGAR)))
		for _, op := range aln.CIGAR {
			index += 4
			binary.LittleEndian.PutUint32(out[index:index+4], uint32((op.Length<<4)|int32(cigarMap[op.Operation])))
		}
	}

	binary.LittleEndian.PutUint32(out[blockSizeIndex:blockSizeIndex+4], uint32(len(out)-blockSizeIndex-4))

	return out
}

// bamReader is an AlignmentFileReader for a BAM InputFile.
type bamReader struct {
	rc         io.Closer
	bgzf       *bgzf.Reader
	references []BAMReference
	buf        []byte
	data       interface{}
}

// Close the BAM input file.
func (reader *bamReader) Close() {
	internal.Close(reader.bgzf)
	if reader.rc != os.Stdin {
		internal.Close(reader.rc)
	}
}

// ParseHeader implements the method of the AlignmentFileReader interface.
func (reader *bamReader) ParseHeader() (hdr *Header) {
	hdr, reader.references = ParseBamHeader(reader.bgzf)
	reader.buf = make([]byte, 4)
	return
}

// SkipHeader implements the method of the AlignmentFileReader interface.
func (reader *bamReader) SkipHeader() {
	reader.references = SkipBamHeader(reader.bgzf)
	reader.buf = make([]byte, 4)
}

// Err implements the method of the pipeline.Source interface.
func (reader *bamReader) Err() error {
	return nil
}

// Prepare implements the method of the pipeline.Source interface.
func (*bamReader) Prepare(_ context.Context) (size int) {
	return -1
}

// Fetch implements the method of the pipeline.Source interface.
func (reader *bamReader) Fetch(size int) (fetched int) {
	var records [][]byte
	for fetched = 0; fetched < size; fetched++ {
		if _, err := io.ReadFull(reader.bgzf, reader.buf); err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		size := int(int32(binary.LittleEndian.Uint32(reader.buf)))
		for cap(reader.buf) < size {
			reader.buf = append(reader.buf[:cap(reader.buf)], 0)
		}
		reader.buf = reader.buf[:size]
		internal.ReadFull(reader.bgzf, reader.buf)
		records = append(records, append([]byte(nil), reader.buf...))
		reader.buf = reader.buf[:4]
	}
	reader.data = records
	return fetched
}

// Data implements the method of the pipeline.Source interface.
func (reader *bamReader) Data() interface{} {
	return reader.data
}

// ParseAlignment implements the method of the AlignmentFileReader interface.
func (reader *bamReader) ParseAlignment(record []byte) *Alignment {
	return parseBamAlignment(record, reader.references)
}

// bamWriter is an AlignmentFileWriter for a BAM OutputFile.
type bamWriter struct {
	dictTable map[string]uint32
	bgzf      *bgzf.Writer
	wc        io.Closer
}

func (writer *bamWriter) Close() {
	internal.Close(writer.bgzf)
	if writer.wc != os.Stdout {
		internal.Close(writer.wc)
	}
}

// FormatHeader implements the method of the AlignmentFileWriter interface.
func (writer *bamWriter) FormatHeader(hdr *Header) {
	dictTable := make(map[string]uint32)
	dictTable["*"] = minus1
	for index, entry := range hdr.SQ {
		dictTable[entry["SN"]] = uint32(index)
	}
	writer.dictTable = dictTable
	writer.Write(hdr.FormatBam(nil))
}

// FormatAlignment implements the method of the AlignmentFileWriter interface.
func (writer *bamWriter) FormatAlignment(aln *Alignment, out []byte) []byte {
	return formatBamAlignment(aln, out, writer.dictTable)
}

func (writer *bamWriter) Write(p []byte) int {
	n, err := writer.bgzf.Write(p)
	if err != nil {
		log.Panic(err)
	}
	return n
}
