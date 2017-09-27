package sam

import (
	"log"
	"runtime"
	"sync/atomic"
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/utils"
)

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

func (aln *Alignment) ComputePhredScore() (score int32) {
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

var (
	clippedTable   = map[byte]byte{'S': 1, 'H': 1}
	referenceTable = map[byte]byte{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}
)

func (aln *Alignment) ComputeUnclippedPosition() (result int32) {
	cigar, err := ScanCigarString(aln.CIGAR)
	if err != nil {
		log.Fatal(err.Error(), ", while scanning CIGAR string for ", aln.QNAME, " in ComputeUnclippedPosition")
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

func adaptedPos(aln *Alignment) int32 {
	p, ok := aln.Temps.Get(pos)
	if !ok {
		log.Fatal("Unclipped position not present in SAM alignment ", aln.QNAME)
	}
	return p.(int32)
}

func setAdaptedPos(aln *Alignment, p int32) {
	aln.Temps.Set(pos, p)
}

func adaptedScore(aln *Alignment) int32 {
	s, ok := aln.Temps.Get(score)
	if !ok {
		log.Fatal("Phred score not present in SAM alignment ", aln.QNAME)
	}
	return s.(int32)
}

func setAdaptedScore(aln *Alignment, s int32) {
	aln.Temps.Set(score, s)
}

func adaptAlignment(aln *Alignment) {
	aln.SetRG(utils.Intern(aln.RG().(string)))
	setAdaptedPos(aln, aln.ComputeUnclippedPosition())
	setAdaptedScore(aln, aln.ComputePhredScore())
}

type handle struct {
	object unsafe.Pointer
}

func newAlignmentHandle(aln *Alignment) *handle {
	return &handle{unsafe.Pointer(aln)}
}

func (h *handle) alignment() *Alignment {
	return (*Alignment)(h.object)
}

func (h *handle) compareAndSwapAlignment(old, new *Alignment) bool {
	return atomic.CompareAndSwapPointer(&h.object, unsafe.Pointer(old), unsafe.Pointer(new))
}

func isTrueFragment(aln *Alignment) bool {
	return (aln.FLAG & (Multiple | NextUnmapped)) != Multiple
}

func isTruePair(aln *Alignment) bool {
	return (aln.FLAG & (Multiple | NextUnmapped)) == Multiple
}

type fragment struct {
	rg       utils.Symbol
	refid    int32
	pos      int32
	reversed bool
}

func (f fragment) Hash() uint64 {
	return utils.SymbolHash(f.rg) ^ uint64(f.refid) ^ uint64(f.pos) ^ internal.BoolHash(f.reversed)
}

func classifyFragment(aln *Alignment, fragments *sync.Map, deterministic bool) {
	entry, found := fragments.LoadOrStore(fragment{
		aln.RG().(utils.Symbol),
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
				aln.FLAG |= Duplicate
				break
			} else if bestAlnScore := adaptedScore(bestAln); bestAlnScore > alnScore {
				aln.FLAG |= Duplicate
				break
			} else if bestAlnScore == alnScore {
				if deterministic {
					if aln.QNAME > bestAln.QNAME {
						aln.FLAG |= Duplicate
						break
					} else if best.compareAndSwapAlignment(bestAln, aln) {
						bestAln.FLAG |= Duplicate
						break
					}
				} else {
					aln.FLAG |= Duplicate
					break
				}
			} else if best.compareAndSwapAlignment(bestAln, aln) {
				bestAln.FLAG |= Duplicate
				break
			}
		}
	} else {
		for {
			if bestAln := best.alignment(); isTruePair(bestAln) {
				break
			} else if best.compareAndSwapAlignment(bestAln, aln) {
				bestAln.FLAG |= Duplicate
				break
			}
		}
	}
}

type pairFragment struct {
	rg    utils.Symbol
	qname string
}

func (f pairFragment) Hash() uint64 {
	return utils.SymbolHash(f.rg) ^ internal.StringHash(f.qname)
}

type pair struct {
	rg                   utils.Symbol
	refid1, refid2       int32
	pos                  int64
	reversed1, reversed2 bool
}

func (p pair) Hash() uint64 {
	return utils.SymbolHash(p.rg) ^ uint64(p.refid1) ^ uint64(p.refid2) ^ uint64(p.pos) ^ internal.BoolHash(p.reversed1) ^ internal.BoolHash(p.reversed2)
}

type samAlignmentPair struct {
	score      int32
	aln1, aln2 *Alignment
}

func newPairHandle(score int32, aln1, aln2 *Alignment) *handle {
	return &handle{unsafe.Pointer(&samAlignmentPair{score, aln1, aln2})}
}

func (h *handle) pair() *samAlignmentPair {
	return (*samAlignmentPair)(h.object)
}

func (h *handle) compareAndSwapPair(old, new *samAlignmentPair) bool {
	return atomic.CompareAndSwapPointer(&h.object, unsafe.Pointer(old), unsafe.Pointer(new))
}

func classifyPair(aln *Alignment, fragments, pairs *sync.Map, deterministic bool) {
	if !isTruePair(aln) {
		return
	}

	aln1 := aln
	var aln2 *Alignment
	if entry, deleted := fragments.DeleteOrStore(pairFragment{aln.RG().(utils.Symbol), aln.QNAME}, aln); deleted {
		aln2 = entry.(*Alignment)
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
		aln1.RG().(utils.Symbol),
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
			aln1.FLAG |= Duplicate
			aln2.FLAG |= Duplicate
			break
		} else if bestPair.score == score {
			if deterministic {
				if aln1.QNAME > bestPair.aln1.QNAME {
					aln1.FLAG |= Duplicate
					aln2.FLAG |= Duplicate
					break
				} else if best.compareAndSwapPair(bestPair, newPair()) {
					bestPair.aln1.FLAG |= Duplicate
					bestPair.aln2.FLAG |= Duplicate
					break
				}
			} else {
				aln1.FLAG |= Duplicate
				aln2.FLAG |= Duplicate
				break
			}
		} else if best.compareAndSwapPair(bestPair, newPair()) {
			bestPair.aln1.FLAG |= Duplicate
			bestPair.aln2.FLAG |= Duplicate
			break
		}
	}
}

func MarkDuplicates(deterministic bool) Filter {
	return func(_ *Header) AlignmentFilter {
		splits := 16 * runtime.GOMAXPROCS(0)
		fragments := sync.NewMap(splits)
		pairsFragments := sync.NewMap(splits)
		pairs := sync.NewMap(splits)
		return func(aln *Alignment) bool {
			if aln.FlagNotAny(Unmapped | Secondary | Duplicate | Supplementary) {
				adaptAlignment(aln)
				classifyFragment(aln, fragments, deterministic)
				classifyPair(aln, pairsFragments, pairs, deterministic)
			}
			return true
		}
	}
}
