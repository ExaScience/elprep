package filters

import (
	"log"
	"runtime"
	"sync/atomic"
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/sam"
	"github.com/exascience/elprep/utils"
)

// Map Phred qualities to a reasonable range and an error flag
// indicating if it is outside a valid range.
var phredScoreTable [512]byte

func init() {
	for char := 0; char < 256; char++ {
		pos := char << 1
		if (char < 33) || (char > 126) {
			phredScoreTable[pos] = 0
			phredScoreTable[pos+1] = 1
		} else {
			qual := char - 33
			if qual >= 15 {
				phredScoreTable[pos] = byte(qual)
			} else {
				phredScoreTable[pos] = 0
			}
			phredScoreTable[pos+1] = 0
		}
	}
}

// computePhredScore sums the adapted Phred qualities of an alignment.
func computePhredScore(aln *sam.Alignment) (score int32) {
	var error int32
	for _, char := range aln.QUAL {
		pos := char << 1
		score += int32(phredScoreTable[pos])
		error |= int32(phredScoreTable[pos+1])
	}
	if error != 0 {
		log.Fatal("Invalid QUAL character in ", aln.QUAL)
	}
	return score
}

// Map CIGAR operations to flags indicating whether they are clipped
// and/or reference operations.
var (
	clippedTable   = map[byte]byte{'S': 1, 'H': 1}
	referenceTable = map[byte]byte{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}
)

// computeUnclippedPosition determines the unclipped position of an
// alignment, based on its FLAG, POS, and CIGAR string.
func computeUnclippedPosition(aln *sam.Alignment) (result int32) {
	cigar, err := sam.ScanCigarString(aln.CIGAR)
	if err != nil {
		log.Fatal(err, ", while scanning CIGAR string for ", aln.QNAME, " in computeUnclippedPosition")
	}

	result = aln.POS

	if len(cigar) == 0 {
		return result
	}

	if aln.IsReversed() {
		clipped := int32(1)
		result--
		for i := len(cigar) - 1; i >= 0; i-- {
			op := cigar[i]
			p := op.Operation
			c := int32(clippedTable[p])
			r := int32(referenceTable[p])
			clipped *= c
			result += (r | clipped) * op.Length
		}
	} else {
		for _, op := range cigar {
			p := op.Operation
			if clippedTable[p] == 0 {
				break
			}
			result -= op.Length
		}
	}

	return result
}

var (
	pos   = utils.Intern("pos")
	score = utils.Intern("score")
)

func adaptedPos(aln *sam.Alignment) int32 {
	p, ok := aln.Temps.Get(pos)
	if !ok {
		log.Fatal("Unclipped position not present in SAM alignment ", aln.QNAME)
	}
	return p.(int32)
}

func setAdaptedPos(aln *sam.Alignment, p int32) {
	aln.Temps.Set(pos, p)
}

func adaptedScore(aln *sam.Alignment) int32 {
	s, ok := aln.Temps.Get(score)
	if !ok {
		log.Fatal("Phred score not present in SAM alignment ", aln.QNAME)
	}
	return s.(int32)
}

func setAdaptedScore(aln *sam.Alignment, s int32) {
	aln.Temps.Set(score, s)
}

// Adapt the sam-alignment: Fill in library id; fill in unclipped
// position; fill in Phred score.
func adaptAlignment(aln *sam.Alignment, lbTable map[string]string) {
	rg := aln.RG()
	if rg != nil {
		lb, lbFound := lbTable[rg.(string)]
		if lbFound {
			aln.SetLIBID(lb)
		}
	}
	setAdaptedPos(aln, computeUnclippedPosition(aln))
	setAdaptedScore(aln, computePhredScore(aln))
}

// A handle wraps pointers in a box to enable using
// atomic.CompareAndSwapPointer safely.
type handle struct {
	object unsafe.Pointer
}

func newAlignmentHandle(aln *sam.Alignment) *handle {
	return &handle{unsafe.Pointer(aln)}
}

func (h *handle) alignment() *sam.Alignment {
	return (*sam.Alignment)(h.object)
}

func (h *handle) compareAndSwapAlignment(old, new *sam.Alignment) bool {
	return atomic.CompareAndSwapPointer(&h.object, unsafe.Pointer(old), unsafe.Pointer(new))
}

// Is this alignment definitely not part of a pair?
func isTrueFragment(aln *sam.Alignment) bool {
	return (aln.FLAG & (sam.Multiple | sam.NextUnmapped)) != sam.Multiple
}

// Is this alignment definitely part of pair?
func isTruePair(aln *sam.Alignment) bool {
	return (aln.FLAG & (sam.Multiple | sam.NextUnmapped)) == sam.Multiple
}

// The portion of an alignment that indicates its unclipped position
// and its direction.
type fragment struct {
	lb       interface{}
	refid    int32
	pos      int32
	reversed bool
}

func (f fragment) Hash() (hash uint64) {
	if f.lb != nil {
		hash = internal.StringHash(f.lb.(string))
	}
	return hash ^ uint64(f.refid) ^ uint64(f.pos) ^ internal.BoolHash(f.reversed)
}

// For each set of alignments with the same unclipped position and
// direction, all except the one with the highest score are marked as
// duplicates. If there are fragments in such a list that are actually
// part of pairs, all the true fragments are marked as duplicates and
// the pairs are left untouched.
//
// If multiple framents are tied for best score, and deterministic is
// true, all except the one with the lexicographically smallest QNAME
// are marked as duplicates.  If deterministic is false, the choice
// which of the tied fragments are marked as duplicates is random.
func classifyFragment(aln *sam.Alignment, fragments *sync.Map, deterministic bool) {
	entry, found := fragments.LoadOrStore(fragment{
		aln.LIBID(),
		aln.REFID(),
		adaptedPos(aln),
		aln.IsReversed(),
	}, newAlignmentHandle(aln))
	if !found {
		return
	}
	best := entry.(*handle)

	if isTrueFragment(aln) {
		alnScore := adaptedScore(aln)
		for {
			if bestAln := best.alignment(); isTruePair(bestAln) {
				aln.FLAG |= sam.Duplicate
				break
			} else if bestAlnScore := adaptedScore(bestAln); bestAlnScore > alnScore {
				aln.FLAG |= sam.Duplicate
				break
			} else if bestAlnScore == alnScore {
				if deterministic {
					if aln.QNAME > bestAln.QNAME {
						aln.FLAG |= sam.Duplicate
						break
					} else if best.compareAndSwapAlignment(bestAln, aln) {
						bestAln.FLAG |= sam.Duplicate
						break
					}
				} else {
					aln.FLAG |= sam.Duplicate
					break
				}
			} else if best.compareAndSwapAlignment(bestAln, aln) {
				bestAln.FLAG |= sam.Duplicate
				break
			}
		}
	} else {
		for {
			if bestAln := best.alignment(); isTruePair(bestAln) {
				break
			} else if best.compareAndSwapAlignment(bestAln, aln) {
				bestAln.FLAG |= sam.Duplicate
				break
			}
		}
	}
}

// The portion of an alignments that indicates the pair it belongs to.
type pairFragment struct {
	lb    interface{}
	qname string
}

func (f pairFragment) Hash() (hash uint64) {
	if f.lb != nil {
		hash = internal.StringHash(f.lb.(string))
	}
	return hash ^ internal.StringHash(f.qname)
}

// The portion of two alignments forming a pair that indicates their
// unclipped positions and their directions.
type pair struct {
	lb                   interface{}
	refid1, refid2       int32
	pos                  int64
	reversed1, reversed2 bool
}

func (p pair) Hash() (hash uint64) {
	if p.lb != nil {
		hash = internal.StringHash(p.lb.(string))
	}
	return hash ^ uint64(p.refid1) ^ uint64(p.refid2) ^ uint64(p.pos) ^ internal.BoolHash(p.reversed1) ^ internal.BoolHash(p.reversed2)
}

type samAlignmentPair struct {
	score      int32
	aln1, aln2 *sam.Alignment
}

func newPairHandle(score int32, aln1, aln2 *sam.Alignment) *handle {
	return &handle{unsafe.Pointer(&samAlignmentPair{score, aln1, aln2})}
}

func (h *handle) pair() *samAlignmentPair {
	return (*samAlignmentPair)(h.object)
}

func (h *handle) compareAndSwapPair(old, new *samAlignmentPair) bool {
	return atomic.CompareAndSwapPointer(&h.object, unsafe.Pointer(old), unsafe.Pointer(new))
}

// For each set of pairs with the same unclipped positions and
// directions, all except the one with the highest score are marked as
// duplicates.
//
// If multiple pairs are tied for best score, and deterministic is
// true, all except the one with the lexicographically smallest QNAME
// are marked as duplicates.  If deterministic is false, the choice
// which of the tied pairs are marked as duplicates is random.
func classifyPair(aln *sam.Alignment, fragments, pairs *sync.Map, deterministic bool) {
	if !isTruePair(aln) {
		return
	}

	aln1 := aln
	var aln2 *sam.Alignment
	if entry, deleted := fragments.DeleteOrStore(pairFragment{aln.LIBID(), aln.QNAME}, aln); deleted {
		aln2 = entry.(*sam.Alignment)
	} else {
		return
	}

	score := adaptedScore(aln1) + adaptedScore(aln2)
	aln1Pos := adaptedPos(aln1)
	aln2Pos := adaptedPos(aln2)
	if aln1Pos > aln2Pos {
		aln1, aln2 = aln2, aln1
		aln1Pos, aln2Pos = aln2Pos, aln1Pos
	}
	entry, found := pairs.LoadOrStore(pair{
		aln1.LIBID(),
		aln1.REFID(),
		aln2.REFID(),
		(int64(aln1Pos) << 32) + int64(aln2Pos),
		aln1.IsReversed(),
		aln2.IsReversed(),
	}, newPairHandle(score, aln1, aln2))
	if !found {
		return
	}
	best := entry.(*handle)

	var np *samAlignmentPair
	newPair := func() *samAlignmentPair {
		if np == nil {
			np = &samAlignmentPair{score, aln1, aln2}
		}
		return np
	}

	for {
		if bestPair := best.pair(); bestPair.score > score {
			aln1.FLAG |= sam.Duplicate
			aln2.FLAG |= sam.Duplicate
			break
		} else if bestPair.score == score {
			if deterministic {
				if aln1.QNAME > bestPair.aln1.QNAME {
					aln1.FLAG |= sam.Duplicate
					aln2.FLAG |= sam.Duplicate
					break
				} else if best.compareAndSwapPair(bestPair, newPair()) {
					bestPair.aln1.FLAG |= sam.Duplicate
					bestPair.aln2.FLAG |= sam.Duplicate
					break
				}
			} else {
				aln1.FLAG |= sam.Duplicate
				aln2.FLAG |= sam.Duplicate
				break
			}
		} else if best.compareAndSwapPair(bestPair, newPair()) {
			bestPair.aln1.FLAG |= sam.Duplicate
			bestPair.aln2.FLAG |= sam.Duplicate
			break
		}
	}
}

// MarkDuplicates returns a filter for marking duplicate
// alignments. Depends on the AddREFID filter being called before to
// fill in the refid.
//
// Duplicate marking is based on an adapted Phred score. In case of
// ties, if deterministic is true, the QNAME is used as a
// tie-breaker. Otherwise duplicate marking is random for alignments
// tied for best score.
func MarkDuplicates(deterministic bool) sam.Filter {
	return func(header *sam.Header) sam.AlignmentFilter {
		splits := 16 * runtime.GOMAXPROCS(0)
		fragments := sync.NewMap(splits)
		pairsFragments := sync.NewMap(splits)
		pairs := sync.NewMap(splits)
		// map read groups to library ids
		lbTable := make(map[string]string)
		for _, rgEntry := range header.RG {
			lb, found := rgEntry["LB"]
			if found {
				id, found := rgEntry["ID"]
				if !found {
					log.Fatal("Missing mandatory ID entry in an @RG line in a SAM file header.")
				}
				lbTable[id] = lb
			}
		}
		return func(aln *sam.Alignment) bool {
			if aln.FlagNotAny(sam.Unmapped | sam.Secondary | sam.Duplicate | sam.Supplementary) {
				adaptAlignment(aln, lbTable)
				classifyFragment(aln, fragments, deterministic)
				classifyPair(aln, pairsFragments, pairs, deterministic)
			}
			return true
		}
	}
}
