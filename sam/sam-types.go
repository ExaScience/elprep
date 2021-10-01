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
	"log"
	"sort"
	"strconv"
	"sync"
	"unicode"
	"unsafe"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/utils/nibbles"

	psort "github.com/exascience/pargo/sort"

	"github.com/exascience/elprep/v5/utils"
)

// The SAM file format version and date strings supported by this
// library. This is entered by default in an @HD line in the header
// section of a SAM file, unless user code explicitly asks for a
// different version number. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
const (
	FileFormatVersion = "1.6"
	FileFormatDate    = "22 May 2018"
)

// IsHeaderUserTag determins whether this tag string represent a
// user-defined tag.
func IsHeaderUserTag(code string) bool {
	for i := 0; i < len(code); i++ {
		if c := code[i]; ('a' <= c) && (c <= 'z') {
			return true
		}
	}
	return false
}

// Header represents the information stored in the header section of a
// SAM file. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3.
//
// Each line (except for @CO) is represented as a map[string]string,
// mapping string tags to string values.
//
// The zero Header is valid and empty.
type Header struct {
	// The @HD line.
	HD utils.StringMap

	// The @SQ, @RG, and @PG lines, in the order they occur in the
	// header.
	SQ, RG, PG []utils.StringMap

	// The @CO lines in the order they occur in the header.
	CO []string

	// The lines with user-defined @ tags, for each tag in the order
	// they occur in the header.
	UserRecords map[string][]utils.StringMap
}

type (
	// SortingOrder represents the possible values for the SO tag stored
	// in the @HD line of a header.  See
	// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag
	// @HD.
	SortingOrder string

	// GroupingOrder represents the possible values for the GO tag stored
	// in the @HD line of a header.  See
	// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag
	// @HD.
	GroupingOrder string
)

// Sorting orders.
const (
	Keep       SortingOrder = "keep"
	Unknown    SortingOrder = "unknown"
	Unsorted   SortingOrder = "unsorted"
	Queryname  SortingOrder = "queryname"
	Coordinate SortingOrder = "coordinate"
)

// Grouping orders.
const (
	None      GroupingOrder = "none"
	Query     GroupingOrder = "query"
	Reference GroupingOrder = "reference"
)

// SQLN returns he LN field value, assuming that the given record
// represents an @SQ line in the the header section of a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
//
// If the LN field is present, error is nil unless the value cannot be
// successfully parsed into an int32. If the LN field is not present,
// SQLN returns the maximum possible value for LN and a non-nil error
// value.
func SQLN(record utils.StringMap) int32 {
	ln, found := record["LN"]
	if !found {
		log.Panic("LN entry in a SQ header line missing")
	}
	return int32(internal.ParseInt(ln, 10, 32))
}

// SetSQLN sets the LN field value, assumming that the given record
// represents an @SQ line in the header section of a SAM file. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
func SetSQLN(record utils.StringMap, value int32) {
	record["LN"] = strconv.FormatInt(int64(value), 10)
}

// NewHeader allocates and initializes an empty header.
func NewHeader() *Header { return &Header{} }

// EnsureHD ensures that an @HD line is present in the given
// header. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 1.3, Tag @HD.
//
// If an @HD line already exists, it is returned unchanged. Otherwise,
// the HD field is initialized with a default VN value.
func (hdr *Header) EnsureHD() utils.StringMap {
	if hdr.HD == nil {
		hdr.HD = utils.StringMap{"VN": FileFormatVersion}
	}
	return hdr.HD
}

// HDSO returns the sorting order (SO) stored in the @HD line of the
// given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD.
//
// If there is no @HD line, or the SO field is not set, returns
// "unknown".
func (hdr *Header) HDSO() SortingOrder {
	hd := hdr.EnsureHD()
	if sortingOrder, found := hd["SO"]; found {
		return SortingOrder(sortingOrder)
	}
	return Unknown
}

// SetHDSO sets the sorting order (SO) stored in the @HD line of the
// given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD.
//
// This also deletes the value for the GO field if it is set.
func (hdr *Header) SetHDSO(value SortingOrder) {
	hd := hdr.EnsureHD()
	delete(hd, "GO")
	hd["SO"] = string(value)
}

// HDGO returns the grouping order (GO) stored in the @HD line of the
// given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD.
//
// If there is no @HD line, or the GO field is not set, returns
// "none".
func (hdr *Header) HDGO() GroupingOrder {
	hd := hdr.EnsureHD()
	if groupingOrder, found := hd["GO"]; found {
		return GroupingOrder(groupingOrder)
	}
	return None
}

// SetHDGO sets the grouping order (GO) stored in the @HD line of the
// given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD.
//
// This also deletes the value for the SO field if it is set.
func (hdr *Header) SetHDGO(value GroupingOrder) {
	hd := hdr.EnsureHD()
	delete(hd, "SO")
	hd["GO"] = string(value)
}

// EnsureUserRecords ensures that a map for user-defined @ tags exists
// in the given header. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag
// @HD.
//
// If the map already exists, it is returned unchanged. Otherwise, the
// UserRecords field is initialized with an empty map.
func (hdr *Header) EnsureUserRecords() map[string][]utils.StringMap {
	if hdr.UserRecords == nil {
		hdr.UserRecords = make(map[string][]utils.StringMap)
	}
	return hdr.UserRecords
}

// AddUserRecord adds a header line for the given user-defined @ tag
// to the given header. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag
// @HD.
func (hdr *Header) AddUserRecord(code string, record utils.StringMap) {
	if records, found := hdr.UserRecords[code]; found {
		hdr.UserRecords[code] = append(records, record)
	} else {
		hdr.EnsureUserRecords()[code] = []utils.StringMap{record}
	}
}

var (
	nibbleToBase = []byte("=ACMGRSVTWYHKDBN")
	baseToNibble = make(map[byte]byte)
)

func init() {
	for i, b := range nibbleToBase {
		baseToNibble[b] = byte(i)
	}
}

func nibblesToByteSlice(seq nibbles.Nibbles) []byte {
	slice := seq.Expand()
	for i, b := range slice {
		slice[i] = nibbleToBase[b]
	}
	return slice
}

func nibblesToString(seq nibbles.Nibbles) string {
	slice := nibblesToByteSlice(seq)
	return *(*string)(unsafe.Pointer(&slice))
}

// Sequence encodes a SAM segment SEQuence as in the BAM format.
// See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 4.2.
type Sequence nibbles.Nibbles

// Len returns the length of a SAM segment SEQuence.
func (seq Sequence) Len() int {
	return nibbles.Nibbles(seq).Len()
}

// Slice slices a SAM segment SEQuence.
func (seq Sequence) Slice(low, high int) Sequence {
	return Sequence(nibbles.Nibbles(seq).Slice(low, high))
}

// Base returns the base in a SAM segment SEQuence at the given position.
func (seq Sequence) Base(i int) (base byte) {
	return nibbleToBase[int(nibbles.Nibbles(seq).Get(i))]
}

// SetBase sets the base in a SAM segment SEQuence at the given position.
func (seq Sequence) SetBase(i int, base byte) {
	nibble, ok := baseToNibble[base]
	if !ok {
		nibble = 15
	}
	nibbles.Nibbles(seq).Set(i, nibble)
}

// AsString returns the sequence representation as in the SAM format
func (seq Sequence) AsString() string {
	return nibblesToString(nibbles.Nibbles(seq))
}

// An Alignment represents a single read alignment with mandatory and
// optional fields that can be contained in a SAM file alignment
// line. See http://samtools.github.io/hts-specs/SAMv1.pdf - Sections
// 1.4 and 1.5. SEQ and QUAL are represented as in the BAM format, see
// Section 4.2.
type Alignment struct {
	// The Query template NAME.
	QNAME string

	// The Reference sequence NAME.
	RNAME string

	// The 1-based leftmost mapping POSition (as in the SAM format).
	POS int32

	// The bitwise FLAG.
	FLAG uint16

	// The MAPping Quality.
	MAPQ byte

	// The CIGAR string as a slice of CIGAR operations.
	CIGAR []CigarOperation

	// The Reference sequence name of the mate/NEXT read.
	RNEXT string

	// The 1-based leftmost mapping Position of the make/NEXT read (as in the SAM format).
	PNEXT int32

	// The observed Template LENgth.
	TLEN int32

	// The segment SEQuence (as in the BAM format).
	SEQ Sequence

	// The ASCII of Phred-scaled base QUALity.
	// A slice of the Phred-scaled base quality values (as in the BAM format,
	// without the increment of 33 to turn the values into printable ASCII characters).
	QUAL []byte

	// The optional fields in a read alignment.
	TAGS utils.SmallMap

	// Additional optional fields which are not stored in SAM files, but
	// reserved for temporary values in filters.
	Temps utils.SmallMap
}

// Symbols for some commonly used optional fields. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
var (
	CC = utils.Intern("CC")
	LB = utils.Intern("LB")
	PG = utils.Intern("PG")
	PU = utils.Intern("PU")
	RG = utils.Intern("RG")
)

// Symbols for some temporary fields.
var (
	LIBID     = utils.Intern("LIBID")
	REFID     = utils.Intern("REFID")
	NextREFID = utils.Intern("NextREFID")
)

// RG returns the (potentially empty) RG optional field.
func (aln *Alignment) RG() interface{} {
	rg, _ := aln.TAGS.Get(RG)
	return rg
}

// SetRG sets the RG optional field.
func (aln *Alignment) SetRG(rg interface{}) {
	aln.TAGS.Set(RG, rg)
}

// SetPG sets the PG optional field.
func (aln *Alignment) SetPG(pg interface{}) {
	aln.TAGS.Set(PG, pg)
}

// REFID returns the REFID temporary field.
//
// If REFID field is not set, this will panic with a log message. The
// AddREFID filter can be used to avoid this situation.  (The elPrep
// command line tool ensures that AddREFID is correctly used for its
// default pipelines.)
func (aln *Alignment) REFID() int32 {
	refid, ok := aln.Temps.Get(REFID)
	if !ok {
		log.Panic("REFID in SAM alignment ", aln.QNAME, " not set (use the AddREFID filter to fix this)")
	}
	return refid.(int32)
}

// SetREFID sets the REFID temporary field.
func (aln *Alignment) SetREFID(refid int32) {
	aln.Temps.Set(REFID, refid)
}

func (aln *Alignment) NextREFID() int32 {
	refid, ok := aln.Temps.Get(NextREFID)
	if !ok {
		log.Panic("NextREFID in SAM alignment ", aln.QNAME, " not set (use the AddREFID filter to fix this)")
	}
	return refid.(int32)
}

func (aln *Alignment) SetNextREFID(refid int32) {
	aln.Temps.Set(NextREFID, refid)
}

// LIBID returns the LIBID temporary field.
func (aln *Alignment) LIBID() interface{} {
	lb, _ := aln.Temps.Get(LIBID)
	return lb
}

// SetLIBID sets the LIBID temporary field.
func (aln *Alignment) SetLIBID(libid interface{}) {
	aln.Temps.Set(LIBID, libid)
}

func modFlag(flag uint16) uint16 {
	if flag&Multiple == 0 {
		flag &^= NextUnmapped
		flag &^= NextReversed
	}
	if flag&Unmapped != 0 {
		flag &^= Reversed
	}
	if flag&NextUnmapped != 0 {
		flag &^= NextReversed
	}
	return flag
}

// CoordinateLess compares two alignments according to their
// coordinate. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.3, Tag @HD, SO.
func CoordinateLess(aln1, aln2 *Alignment) bool {
	refid1 := aln1.REFID()
	refid2 := aln2.REFID()
	switch {
	case refid1 < refid2:
		return refid1 >= 0
	case refid2 < refid1:
		return refid2 < 0
	case aln1.POS < aln2.POS:
		return true
	case aln1.POS > aln2.POS:
		return false
	case aln1.IsReversed() != aln2.IsReversed():
		return !aln1.IsReversed()
	case aln1.QNAME != "" && aln2.QNAME != "":
		switch {
		case aln1.QNAME < aln2.QNAME:
			return true
		case aln1.QNAME > aln2.QNAME:
			return false
		}
	}
	flag1 := modFlag(aln1.FLAG)
	flag2 := modFlag(aln2.FLAG)
	switch {
	case flag1 < flag2:
		return true
	case flag1 > flag2:
		return false
	case aln1.MAPQ < aln2.MAPQ:
		return true
	case aln1.MAPQ > aln2.MAPQ:
		return false
	case aln1.IsMultiple() && aln2.IsMultiple():
		nextRefid1 := aln1.NextREFID()
		nextRefid2 := aln2.NextREFID()
		switch {
		case nextRefid1 < nextRefid2:
			return true // no special treatment of negative values!
		case nextRefid1 > nextRefid2:
			return false // no special treatment of negative values!
		case aln1.PNEXT < aln2.PNEXT:
			return true
		case aln1.PNEXT > aln2.PNEXT:
			return false
		}
	}
	return aln1.TLEN < aln2.TLEN
}

// QNAMELess compares two alignments according to their
// query template name. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section
// 1.3, Tag @HD, SO.
func QNAMELess(aln1, aln2 *Alignment) bool {
	return aln1.QNAME < aln2.QNAME
}

// Bit values for the FLAG field in the Alignment struct. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
const (
	// Template having multiple segments in sequencing.
	Multiple = 0x1

	// Each segment properly aligned according to the aligner.
	Proper = 0x2

	// Segment unmapped.
	Unmapped = 0x4

	// Next segment in the template unmapped.
	NextUnmapped = 0x8

	// SEQ being reversed complemented.
	Reversed = 0x10

	// SEQ of the next segment in the template being reverse
	// complemented.
	NextReversed = 0x20

	// The first segment in the template.
	First = 0x40

	// The last segment in the template.
	Last = 0x80

	// Secondary alignment.
	Secondary = 0x100

	// Not passing filters, such as platform/vendor quality controls.
	QCFailed = 0x200

	// PCR or optical duplicate.
	Duplicate = 0x400

	// Supplementary alignment.
	Supplementary = 0x800
)

// IsMultiple checks for template having multiple segments in
// sequencing. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
func (aln *Alignment) IsMultiple() bool { return (aln.FLAG & Multiple) != 0 }

// IsProper checks for each segment being properly aligned according
// to the aligner. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
func (aln *Alignment) IsProper() bool { return (aln.FLAG & Proper) != 0 }

// IsUnmapped checks for segment unmapped. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsUnmapped() bool { return (aln.FLAG & Unmapped) != 0 }

// IsNextUnmapped checks for next segment in the template
// unmapped. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.2.
func (aln *Alignment) IsNextUnmapped() bool { return (aln.FLAG & NextUnmapped) != 0 }

// IsReversed checks for SEQ being reversed complemented. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsReversed() bool { return (aln.FLAG & Reversed) != 0 }

// IsNextReversed check for SEQ of the next segment in the template
// being reverse complemented. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsNextReversed() bool { return (aln.FLAG & NextReversed) != 0 }

// IsFirst checks for being the first segment in the template. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsFirst() bool { return (aln.FLAG & First) != 0 }

// IsLast checks for being the last segment in the template. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsLast() bool { return (aln.FLAG & Last) != 0 }

// IsSecondary checks for secondary alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsSecondary() bool { return (aln.FLAG & Secondary) != 0 }

// IsQCFailed checks for not passing filters, such as platform/vendor
// quality controls. See http://samtools.github.io/hts-specs/SAMv1.pdf
// - Section 1.4.2.
func (aln *Alignment) IsQCFailed() bool { return (aln.FLAG & QCFailed) != 0 }

// IsDuplicate checks for PCR or optical duplicate. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsDuplicate() bool { return (aln.FLAG & Duplicate) != 0 }

// IsSupplementary checks for supplementary alignment. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
func (aln *Alignment) IsSupplementary() bool { return (aln.FLAG & Supplementary) != 0 }

// FlagEvery checks for every bit set in the given flag being also set
// in aln.FLAG.
func (aln *Alignment) FlagEvery(flag uint16) bool { return (aln.FLAG & flag) == flag }

// FlagSome checks for some bits set in the given flag being also set
// in aln.FLAG.
func (aln *Alignment) FlagSome(flag uint16) bool { return (aln.FLAG & flag) != 0 }

// FlagNotEvery checks for not every bit set in the given flag being
// also set in aln.FLAG.
func (aln *Alignment) FlagNotEvery(flag uint16) bool { return (aln.FLAG & flag) != flag }

// FlagNotAny checks for not any bit set in the given flag being also
// set in aln.FLAG.
func (aln *Alignment) FlagNotAny(flag uint16) bool { return (aln.FLAG & flag) == 0 }

// By is a type for comparison predicates on Alignment pointers.
type By func(aln1, aln2 *Alignment) bool

// AlignmentSorter is a helper for sorting Alignment slices that
// implements
// https://godoc.org/github.com/ExaScience/pargo/sort#StableSorter
type AlignmentSorter struct {
	alns []*Alignment
	by   By
}

// SequentialSort implements the method of the SequantialSorter interface.
func (s AlignmentSorter) SequentialSort(i, j int) {
	alns, by := s.alns[i:j], s.by
	sort.Slice(alns, func(i, j int) bool {
		return by(alns[i], alns[j])
	})
}

// NewTemp implements the method of the StableSorter interface
func (s AlignmentSorter) NewTemp() psort.StableSorter {
	return AlignmentSorter{make([]*Alignment, len(s.alns)), s.by}
}

// Len implements the method of the sort.Interface.
func (s AlignmentSorter) Len() int {
	return len(s.alns)
}

// Less implements the method of the sort.Interface.
func (s AlignmentSorter) Less(i, j int) bool {
	return s.by(s.alns[i], s.alns[j])
}

// Assign implements the method of the StableSorter interface.
func (s AlignmentSorter) Assign(p psort.StableSorter) func(i, j, len int) {
	dst, src := s.alns, p.(AlignmentSorter).alns
	return func(i, j, len int) {
		for k := 0; k < len; k++ {
			dst[i+k] = src[j+k]
		}
	}
}

// ParallelStableSort sorts a slice of alignments according to the
// given comparison predicate.
func (by By) ParallelStableSort(alns []*Alignment) {
	psort.StableSort(AlignmentSorter{alns, by})
}

// Sam represents a complete SAM data set that can be contained in a
// SAM or BAM file. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.
type Sam struct {
	Header     *Header
	Alignments []*Alignment
	nofBatches int
}

// NewSam allocates and initializes an empty SAM data set.
func NewSam() *Sam { return &Sam{Header: NewHeader()} }

// ByteArray is a representation for byte arrays as stored in optional
// fields of read alignments lines using type H. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
type ByteArray []byte

// cigarOperations contains all valid CIGAR operations.
const cigarOperations = "MmIiDdNnSsHhPpXx="

var cigarOperationsTable = make(map[byte]byte, len(cigarOperations))

func init() {
	for _, c := range cigarOperations {
		cigarOperationsTable[byte(c)] = byte(unicode.ToUpper(c))
	}
}

func isDigit(char byte) bool { return ('0' <= char) && (char <= '9') }

// CigarOperation represents a CIGAR operation. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.
type CigarOperation struct {
	Length    int32
	Operation byte // 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', or 'X'
}

func newCigarOperation(cigar string, i int) (op CigarOperation, j int) {
	for j = i; ; j++ {
		if char := cigar[j]; !isDigit(char) {
			length := internal.ParseInt(cigar[i:j], 10, 32)
			if operation := cigarOperationsTable[char]; operation != 0 {
				op = CigarOperation{int32(length), operation}
				j++
			} else {
				log.Panicf("invalid CIGAR operation %v", operation)
			}
			return
		}
	}
}

var (
	cigarSliceCache      = map[string][]CigarOperation{"*": []CigarOperation{}}
	cigarSliceCacheMutex = sync.RWMutex{}
)

func slowScanCigarString(cigar string) (slice []CigarOperation) {
	if len(cigar) == 0 {
		return nil
	}
	cigarOperation, i := newCigarOperation(cigar, 0)
	slice = []CigarOperation{cigarOperation}
	for i < len(cigar) {
		nextCigarOperation, j := newCigarOperation(cigar, i)
		if nextCigarOperation.Operation == cigarOperation.Operation {
			slice[len(slice)-1].Length += nextCigarOperation.Length
		} else {
			slice = append(slice, nextCigarOperation)
		}
		cigarOperation = nextCigarOperation
		i = j
	}
	cigarSliceCacheMutex.Lock()
	if value, found := cigarSliceCache[cigar]; found {
		slice = value
	} else {
		cigarSliceCache[cigar] = slice
	}
	cigarSliceCacheMutex.Unlock()
	return slice
}

// ScanCigarString converts a CIGAR string to a slice of
// CigarOperation. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Section 1.4.6.
//
// Uses an internal cache to reduce memory overhead. It is safe for
// multiple goroutines to call ScanCigarString concurrently.
func ScanCigarString(cigar string) []CigarOperation {
	cigarSliceCacheMutex.RLock()
	value, found := cigarSliceCache[cigar]
	cigarSliceCacheMutex.RUnlock()
	if found {
		return value
	}
	return slowScanCigarString(cigar)
}

var (
	// CigarOperatorConsumesReadBases maps operators that consume read bases to 1
	CigarOperatorConsumesReadBases = map[byte]int32{'M': 1, 'I': 1, 'S': 1, '=': 1, 'X': 1}
	// CigarOperatorConsumesReferenceBases maps operators that consume reference bases to 1
	CigarOperatorConsumesReferenceBases = map[byte]int32{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}
)

// ReadLengthFromCigar sums the lengths of all CIGAR operations that consume read bases.
func ReadLengthFromCigar(cigar []CigarOperation) int32 {
	var length int32
	for _, op := range cigar {
		length += CigarOperatorConsumesReadBases[op.Operation] * op.Length
	}
	return length
}

// ReferenceLengthFromCigar sums the lengths of all CIGAR operations that consume reference bases
func ReferenceLengthFromCigar(cigar []CigarOperation) int32 {
	var length int32
	for _, op := range cigar {
		length += CigarOperatorConsumesReferenceBases[op.Operation] * op.Length
	}
	return length
}

// End sums the lengths of all CIGAR operations that consume reference
// bases and adds them to the alignment's position - 1
func (aln *Alignment) End() int32 {
	var length int32
	for _, op := range aln.CIGAR {
		length += CigarOperatorConsumesReferenceBases[op.Operation] * op.Length
	}
	return aln.POS + length - 1
}
