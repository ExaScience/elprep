package sam

import (
	"errors"
	"fmt"
	"log"
	"sort"
	"strconv"
	"sync"
	"unicode"

	psort "github.com/exascience/pargo/sort"

	"github.com/exascience/elprep/utils"
)

/*
The SAM file format version and date strings supported by this
library. This is entered by default in an @HD line in the header
section of a SAM file, unless user code explicitly asks for a
different version number. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
const (
	FileFormatVersion = "1.5"
	FileFormatDate    = "1 Jun 2017"
)

/*
IsHeaderUserTag determins whether this tag string represent a
user-defined tag.
*/
func IsHeaderUserTag(code string) bool {
	for _, c := range code {
		if ('a' <= c) && (c <= 'z') {
			return true
		}
	}
	return false
}

/*
Header represents the information stored in the header section of a
SAM file. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section
1.3.

Each line (except for @CO) is represented as a map[string]string,
mapping string tags to string values.

The zero Header is valid and empty.
*/
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

/*
SQLN returns he LN field value, assuming that the given record
represents an @SQ line in the the header section of a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.

If the LN field is present, error is nil unless the value cannot be
successfully parsed into an int32. If the LN field is not present,
SQLN returns the maximum possible value for LN and a non-nil error
value.
*/
func SQLN(record utils.StringMap) (int32, error) {
	ln, found := record["LN"]
	if !found {
		return 0x7FFFFFFF, errors.New("LN entry in a SQ header line missing")
	}
	val, err := strconv.ParseInt(ln, 10, 32)
	if err != nil {
		return 0, err
	}
	return int32(val), nil
}

/*
SetSQLN sets the LN field value, assumming that the given record
represents an @SQ line in the header section of a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
func SetSQLN(record utils.StringMap, value int32) {
	record["LN"] = strconv.FormatInt(int64(value), 10)
}

/*
NewHeader allocates and initializes an empty header.
*/
func NewHeader() *Header { return &Header{} }

/*
EnsureHD ensures that an @HD line is present in the given header. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3, Tag @HD.

If an @HD line already exists, it is returned unchanged. Otherwise,
the HD field is initialized with a default VN value.
*/
func (hdr *Header) EnsureHD() utils.StringMap {
	if hdr.HD == nil {
		hdr.HD = utils.StringMap{"VN": FileFormatVersion}
	}
	return hdr.HD
}

/*
HDSO returns the sorting order (SO) stored in the @HD line of the
given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.

If there is no @HD line, or the SO field is not set, returns
"unknown".
*/
func (hdr *Header) HDSO() string {
	hd := hdr.EnsureHD()
	if sortingOrder, found := hd["SO"]; found {
		return sortingOrder
	}
	return "unknown"
}

/*
SetHDSO sets the sorting order (SO) stored in the @HD line of the
given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.

This also deletes the value for the GO field if it is set.
*/
func (hdr *Header) SetHDSO(value string) {
	hd := hdr.EnsureHD()
	delete(hd, "GO")
	hd["SO"] = value
}

/*
HDGO returns the grouping order (GO) stored in the @HD line of the
given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.

If there is no @HD line, or the GO field is not set, returns "none".
*/
func (hdr *Header) HDGO() string {
	hd := hdr.EnsureHD()
	if groupingOrder, found := hd["GO"]; found {
		return groupingOrder
	}
	return "none"
}

/*
SetHDGO sets the grouping order (GO) stored in the @HD line of the
given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.

This also deletes the value for the SO field if it is set.
*/
func (hdr *Header) SetHDGO(value string) {
	hd := hdr.EnsureHD()
	delete(hd, "SO")
	hd["GO"] = value
}

/*
EnsureUserRecords ensures that a map for user-defined @ tags exists in
the given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.

If the map already exists, it is returned unchanged. Otherwise, the
UserRecords field is initialized with an empty map.
*/
func (hdr *Header) EnsureUserRecords() map[string][]utils.StringMap {
	if hdr.UserRecords == nil {
		hdr.UserRecords = make(map[string][]utils.StringMap)
	}
	return hdr.UserRecords
}

/*
AddUserRecord adds a header line for the given user-defined @ tag to
the given header. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD.
*/
func (hdr *Header) AddUserRecord(code string, record utils.StringMap) {
	if records, found := hdr.UserRecords[code]; found {
		hdr.UserRecords[code] = append(records, record)
	} else {
		hdr.EnsureUserRecords()[code] = []utils.StringMap{record}
	}
}

/*
An Alignment represents a single read alignment with mandatory and
optional fields that can be contained in a SAM file alignment
line. See http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4
and 1.5.
*/
type Alignment struct {
	// The Query template NAME.
	QNAME string

	// The bitwise FLAG.
	FLAG uint16

	// The Reference sequence NAME.
	RNAME string

	// The 1-based leftmost mapping POSition.
	POS int32

	// The MAPping Quality.
	MAPQ byte

	// The CIGAR string.
	CIGAR string

	// The Reference sequence name of the mate/NEXT read.
	RNEXT string

	// The 1-based leftmost mapping Position of the make/NEXT read.
	PNEXT int32

	// The observed Template LENgth.
	TLEN int32

	// The segment SEQuence.
	SEQ string

	// The ASCII of Phred-scaled base QUALity+33.
	QUAL string

	// The optional fields in a read alignment.
	TAGS utils.SmallMap

	// Additional optional fields which are not stored in SAM files, but
	// resereved for temporary values in filters.
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
	LIBID = utils.Intern("LIBID")
	REFID = utils.Intern("REFID")
)

/*
RG returns the (potentially empty) RG optional field.
*/
func (aln *Alignment) RG() interface{} {
	rg, _ := aln.TAGS.Get(RG)
	return rg
}

/*
SetRG sets the RG optional field.
*/
func (aln *Alignment) SetRG(rg interface{}) {
	aln.TAGS.Set(RG, rg)
}

/*
REFID returns the REFID temporary field.

If REFID field is not set, this will panic with a log message. The
AddREFID filter can be used to avoid this situation.  (The elPrep
command line tool ensures that AddREFID is correctly used for its
default pipelines.)
*/
func (aln *Alignment) REFID() int32 {
	refid, ok := aln.Temps.Get(REFID)
	if !ok {
		log.Fatal("REFID in SAM alignment ", aln.QNAME, " not set (use the AddREFID filter to fix this)")
	}
	return refid.(int32)
}

/*
SetREFID sets the REFID temporary field.
*/
func (aln *Alignment) SetREFID(refid int32) {
	aln.Temps.Set(REFID, refid)
}

/*
LIBID returns the LIBID temporary field.
*/
func (aln *Alignment) LIBID() interface{} {
	lb, _ := aln.Temps.Get(LIBID)
	return lb
}

/*
SetLIBID sets the LIBID temporary field.
*/
func (aln *Alignment) SetLIBID(libid interface{}) {
	aln.Temps.Set(LIBID, libid)
}

/*
NewAlignment allocates and initializes an empty alignment.
*/
func NewAlignment() *Alignment {
	return &Alignment{
		TAGS:  make(utils.SmallMap, 0, 16),
		Temps: make(utils.SmallMap, 0, 4),
	}
}

/*
CoordinateLess compares two alignments according to their
coordinate. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.3, Tag @HD, SO.
*/
func CoordinateLess(aln1, aln2 *Alignment) bool {
	refid1 := aln1.REFID()
	refid2 := aln2.REFID()
	switch {
	case refid1 < refid2:
		return refid1 >= 0
	case refid2 < refid1:
		return refid2 < 0
	default:
		return aln1.POS < aln2.POS
	}
}

/*
Bit values for the FLAG field in the Alignment struct. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
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

/*
IsMultiple checks for template having multiple segments in
sequencing. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.4.2.
*/
func (aln *Alignment) IsMultiple() bool { return (aln.FLAG & Multiple) != 0 }

/*
IsProper checks for each segment being properly aligned according to
the aligner. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.4.2.
*/
func (aln *Alignment) IsProper() bool { return (aln.FLAG & Proper) != 0 }

/*
IsUnmapped checks for segment unmapped. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsUnmapped() bool { return (aln.FLAG & Unmapped) != 0 }

/*
IsNextUnmapped checks for next segment in the template unmapped. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsNextUnmapped() bool { return (aln.FLAG & NextUnmapped) != 0 }

/*
IsReversed checks for SEQ being reversed complemented. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsReversed() bool { return (aln.FLAG & Reversed) != 0 }

/*
IsNextReversed check for SEQ of the next segment in the template being
reverse complemented. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsNextReversed() bool { return (aln.FLAG & NextReversed) != 0 }

/*
IsFirst checks for being the first segment in the template. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsFirst() bool { return (aln.FLAG & First) != 0 }

/*
IsLast checks for being the last segment in the template. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsLast() bool { return (aln.FLAG & Last) != 0 }

/*
IsSecondary checks for secondary alignment. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsSecondary() bool { return (aln.FLAG & Secondary) != 0 }

/*
IsQCFailed checks for not passing filters, such as platform/vendor
quality controls. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.4.2.
*/
func (aln *Alignment) IsQCFailed() bool { return (aln.FLAG & QCFailed) != 0 }

/*
IsDuplicate checks for PCR or optical duplicate. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsDuplicate() bool { return (aln.FLAG & Duplicate) != 0 }

/*
IsSupplementary checks for supplementary alignment. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.2.
*/
func (aln *Alignment) IsSupplementary() bool { return (aln.FLAG & Supplementary) != 0 }

/*
FlagEvery checks for every bit set in the given flag being also set in
aln.FLAG.
*/
func (aln *Alignment) FlagEvery(flag uint16) bool { return (aln.FLAG & flag) == flag }

/*
FlagSome checks for some bits set in the given flag being also set in
aln.FLAG.
*/
func (aln *Alignment) FlagSome(flag uint16) bool { return (aln.FLAG & flag) != 0 }

/*
FlagNotEvery checks for not every bit set in the given flag being also
set in aln.FLAG.
*/
func (aln *Alignment) FlagNotEvery(flag uint16) bool { return (aln.FLAG & flag) != flag }

/*
FlagNotAny checks for not any bit set in the given flag being also set
in aln.FLAG.
*/
func (aln *Alignment) FlagNotAny(flag uint16) bool { return (aln.FLAG & flag) == 0 }

// By is a type for comparison predicates on Alignment pointers.
type By func(aln1, aln2 *Alignment) bool

/*
AlignmentSorter is a helper for sorting Alignment slices that
implements
https://godoc.org/github.com/ExaScience/pargo/sort#StableSorter
*/
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

/*
ParallelStableSort sorts a slice of alignments according to the given
comparison predicate.
*/
func (by By) ParallelStableSort(alns []*Alignment) {
	psort.StableSort(AlignmentSorter{alns, by})
}

/*
Sam represents a complete SAM data set that can be contained in a SAM
file. See http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.
*/
type Sam struct {
	Header     *Header
	Alignments []*Alignment
}

/*
NewSam allocates and initializes an empty SAM data set.
*/
func NewSam() *Sam { return &Sam{Header: NewHeader()} }

/*
ByteArray is a representation for byte arrays as stored in optional
fields of read alignments lines using type H. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
type ByteArray []byte

// CigarOperations contains all valid CIGAR operations.
const CigarOperations = "MmIiDdNnSsHhPpXx="

var cigarOperationsTable = make(map[byte]byte, len(CigarOperations))

func init() {
	for _, c := range CigarOperations {
		cigarOperationsTable[byte(c)] = byte(unicode.ToUpper(rune(c)))
	}
}

func isDigit(char byte) bool { return ('0' <= char) && (char <= '9') }

/*
CigarOperation represents a CIGAR operation. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.6.
*/
type CigarOperation struct {
	Length    int32
	Operation byte // 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', or 'X'
}

func newCigarOperation(cigar string, i int) (op CigarOperation, j int, err error) {
	for j = i; ; j++ {
		if char := cigar[j]; !isDigit(char) {
			length, nerr := strconv.ParseInt(cigar[i:j], 10, 32)
			if nerr != nil {
				err = nerr
				return
			}
			if operation := cigarOperationsTable[char]; operation != 0 {
				op = CigarOperation{int32(length), operation}
				j++
			} else {
				err = fmt.Errorf("Invalid CIGAR operation %v", operation)
			}
			return
		}
	}
}

var (
	cigarSliceCache      = map[string][]CigarOperation{"*": []CigarOperation{}}
	cigarSliceCacheMutex = sync.RWMutex{}
)

func slowScanCigarString(cigar string) (slice []CigarOperation, err error) {
	for i := 0; i < len(cigar); {
		cigarOperation, j, err := newCigarOperation(cigar, i)
		if err != nil {
			return nil, fmt.Errorf("%v, while scanning CIGAR string %v", err.Error(), cigar)
		}
		slice = append(slice, cigarOperation)
		i = j
	}
	cigarSliceCacheMutex.Lock()
	if value, found := cigarSliceCache[cigar]; found {
		slice = value
	} else {
		cigarSliceCache[cigar] = slice
	}
	cigarSliceCacheMutex.Unlock()
	return slice, nil
}

/*
ScanCigarString converts a CIGAR string to a slice of
CigarOperation. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.4.6.

Uses an internal cache to reduce memory overhead. It is safe for
multiple goroutines to call ScanCigarString concurrently.
*/
func ScanCigarString(cigar string) ([]CigarOperation, error) {
	cigarSliceCacheMutex.RLock()
	value, found := cigarSliceCache[cigar]
	cigarSliceCacheMutex.RUnlock()
	if found {
		return value, nil
	}
	return slowScanCigarString(cigar)
}
