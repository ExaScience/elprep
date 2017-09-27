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

const (
	FileFormatVersion = "1.5"
	FileFormatDate    = "1 Jun 2017"
)

func IsHeaderUserTag(code string) bool {
	for _, c := range code {
		if ('a' <= c) && (c <= 'z') {
			return true
		}
	}
	return false
}

type Header struct {
	HD          utils.StringMap
	SQ, RG, PG  []utils.StringMap
	CO          []string
	UserRecords map[string][]utils.StringMap
}

func SQ_LN(record utils.StringMap) (int32, error) {
	ln, found := record["LN"]
	if !found {
		return 0x7FFFFFFF, errors.New("LN entry in a SQ header line missing")
	}
	val, err := strconv.ParseInt(ln, 10, 32)
	return int32(val), err
}

func SetSQ_LN(record utils.StringMap, value int32) {
	record["LN"] = strconv.FormatInt(int64(value), 10)
}

func NewHeader() *Header { return &Header{} }

func (hdr *Header) EnsureHD() utils.StringMap {
	if hdr.HD == nil {
		hdr.HD = utils.StringMap{"VN": FileFormatVersion}
	}
	return hdr.HD
}

func (hdr *Header) HD_SO() string {
	hd := hdr.EnsureHD()
	if sortingOrder, found := hd["SO"]; found {
		return sortingOrder
	} else {
		return "unknown"
	}
}

func (hdr *Header) SetHD_SO(value string) {
	hd := hdr.EnsureHD()
	delete(hd, "GO")
	hd["SO"] = value
}

func (hdr *Header) HD_GO() string {
	hd := hdr.EnsureHD()
	if groupingOrder, found := hd["GO"]; found {
		return groupingOrder
	} else {
		return "none"
	}
}

func (hdr *Header) SetHD_GO(value string) {
	hd := hdr.EnsureHD()
	delete(hd, "SO")
	hd["GO"] = value
}

func (hdr *Header) EnsureUserRecords() map[string][]utils.StringMap {
	if hdr.UserRecords == nil {
		hdr.UserRecords = make(map[string][]utils.StringMap)
	}
	return hdr.UserRecords
}

func (hdr *Header) AddUserRecord(code string, record utils.StringMap) {
	if records, found := hdr.UserRecords[code]; found {
		hdr.UserRecords[code] = append(records, record)
	} else {
		hdr.EnsureUserRecords()[code] = []utils.StringMap{record}
	}
}

type Alignment struct {
	QNAME string
	FLAG  uint16
	RNAME string
	POS   int32
	MAPQ  byte
	CIGAR string
	RNEXT string
	PNEXT int32
	TLEN  int32
	SEQ   string
	QUAL  string
	TAGS  utils.SmallMap
	Temps utils.SmallMap
}

var (
	RG    = utils.Intern("RG")
	REFID = utils.Intern("REFID")
)

func (aln *Alignment) RG() interface{} {
	rg, ok := aln.TAGS.Get(RG)
	if !ok {
		log.Fatal("RG in SAM alignment ", aln.QNAME, " not set (use the AddOrReplaceReadGroup filter to fix this)")
	}
	return rg
}

func (aln *Alignment) SetRG(rg interface{}) {
	aln.TAGS.Set(RG, rg)
}

func (aln *Alignment) REFID() int32 {
	refid, ok := aln.Temps.Get(REFID)
	if !ok {
		log.Fatal("REFID in SAM alignment ", aln.QNAME, " not set (use the AddREFID filter to fix this)")
	}
	return refid.(int32)
}

func (aln *Alignment) SetREFID(refid int32) {
	aln.Temps.Set(REFID, refid)
}

func NewAlignment() *Alignment {
	return &Alignment{
		TAGS:  make(utils.SmallMap, 0, 16),
		Temps: make(utils.SmallMap, 0, 4),
	}
}

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

const (
	Multiple      = 0x1
	Proper        = 0x2
	Unmapped      = 0x4
	NextUnmapped  = 0x8
	Reversed      = 0x10
	NextReversed  = 0x20
	First         = 0x40
	Last          = 0x80
	Secondary     = 0x100
	QCFailed      = 0x200
	Duplicate     = 0x400
	Supplementary = 0x800
)

func (aln *Alignment) IsMultiple() bool      { return (aln.FLAG & Multiple) != 0 }
func (aln *Alignment) IsProper() bool        { return (aln.FLAG & Proper) != 0 }
func (aln *Alignment) IsUnmapped() bool      { return (aln.FLAG & Unmapped) != 0 }
func (aln *Alignment) IsNextUnmapped() bool  { return (aln.FLAG & NextUnmapped) != 0 }
func (aln *Alignment) IsReversed() bool      { return (aln.FLAG & Reversed) != 0 }
func (aln *Alignment) IsNextReversed() bool  { return (aln.FLAG & NextReversed) != 0 }
func (aln *Alignment) IsFirst() bool         { return (aln.FLAG & First) != 0 }
func (aln *Alignment) IsLast() bool          { return (aln.FLAG & Last) != 0 }
func (aln *Alignment) IsSecondary() bool     { return (aln.FLAG & Secondary) != 0 }
func (aln *Alignment) IsQCFailed() bool      { return (aln.FLAG & QCFailed) != 0 }
func (aln *Alignment) IsDuplicate() bool     { return (aln.FLAG & Duplicate) != 0 }
func (aln *Alignment) IsSupplementary() bool { return (aln.FLAG & Supplementary) != 0 }

func (aln *Alignment) FlagEvery(flag uint16) bool    { return (aln.FLAG & flag) == flag }
func (aln *Alignment) FlagSome(flag uint16) bool     { return (aln.FLAG & flag) != 0 }
func (aln *Alignment) FlagNotEvery(flag uint16) bool { return (aln.FLAG & flag) != flag }
func (aln *Alignment) FlagNotAny(flag uint16) bool   { return (aln.FLAG & flag) == 0 }

type (
	By func(aln1, aln2 *Alignment) bool

	AlignmentSorter struct {
		alns []*Alignment
		by   By
	}
)

func (s AlignmentSorter) SequentialSort(i, j int) {
	alns, by := s.alns[i:j], s.by
	sort.Slice(alns, func(i, j int) bool {
		return by(alns[i], alns[j])
	})
}

func (s AlignmentSorter) NewTemp() psort.StableSorter {
	return AlignmentSorter{make([]*Alignment, len(s.alns)), s.by}
}

func (s AlignmentSorter) Len() int {
	return len(s.alns)
}

func (s AlignmentSorter) Less(i, j int) bool {
	return s.by(s.alns[i], s.alns[j])
}

func (s AlignmentSorter) Assign(p psort.StableSorter) func(i, j, len int) {
	dst, src := s.alns, p.(AlignmentSorter).alns
	return func(i, j, len int) {
		for k := 0; k < len; k++ {
			dst[i+k] = src[j+k]
		}
	}
}

func (by By) ParallelStableSort(alns []*Alignment) {
	psort.StableSort(AlignmentSorter{alns, by})
}

type Sam struct {
	Header     *Header
	Alignments []*Alignment
}

func NewSam() *Sam { return &Sam{Header: NewHeader()} }

type ByteArray []byte

const CigarOperations = "MmIiDdNnSsHhPpXx="

var cigarOperationsTable = make(map[byte]byte, len(CigarOperations))

func init() {
	for _, c := range CigarOperations {
		cigarOperationsTable[byte(c)] = byte(unicode.ToUpper(rune(c)))
	}
}

func isDigit(char byte) bool { return ('0' <= char) && (char <= '9') }

type CigarOperation struct {
	Length    int32
	Operation byte
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

func ScanCigarString(cigar string) ([]CigarOperation, error) {
	cigarSliceCacheMutex.RLock()
	value, found := cigarSliceCache[cigar]
	cigarSliceCacheMutex.RUnlock()
	if found {
		return value, nil
	} else {
		return slowScanCigarString(cigar)
	}
}
