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

package filters

import (
	"log"
	"runtime"
	"sync/atomic"
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/v4/internal"
	"github.com/exascience/elprep/v4/sam"
	"github.com/exascience/elprep/v4/utils"
)

// Map Phred qualities to a reasonable range and an error flag
// indicating if it is outside a valid range.
var phredScoreTable [512]byte

func init() {
	for char := 0; char < 256; char++ {
		pos := char << 1
		if char > 126-33 {
			phredScoreTable[pos] = 0
			phredScoreTable[pos+1] = 1
		} else {
			qual := char
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

	cigar := aln.CIGAR

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

// Fill in the library id.
func addLIBID(aln *sam.Alignment, lbTable map[string]string) {
	rg := aln.RG()
	if rg != nil {
		lb, lbFound := lbTable[rg.(string)]
		if lbFound {
			aln.SetLIBID(lb)
		}
	}
}

// Adapt the sam-alignment: Fill in unclipped position and Phred score.
func adaptAlignment(aln *sam.Alignment) {
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
	return (*sam.Alignment)(atomic.LoadPointer(&h.object))
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
// If multiple framents are tied for best score, all except the one with
// the lexicographically smallest QNAME are marked as duplicates.
func classifyFragment(aln *sam.Alignment, fragments *sync.Map) {
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
				if aln.QNAME > bestAln.QNAME {
					aln.FLAG |= sam.Duplicate
					break
				} else if best.compareAndSwapAlignment(bestAln, aln) {
					bestAln.FLAG |= sam.Duplicate
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
	score             int32
	aln1, aln2        *sam.Alignment
	opticalDuplicates unsafe.Pointer
}

func newPairHandle(score int32, aln1, aln2 *sam.Alignment) *handle {
	return &handle{unsafe.Pointer(&samAlignmentPair{score: score, aln1: aln1, aln2: aln2})}
}

func (h *handle) pair() *samAlignmentPair {
	return (*samAlignmentPair)(atomic.LoadPointer(&h.object))
}

func (h *handle) compareAndSwapPair(old, new *samAlignmentPair) bool {
	return atomic.CompareAndSwapPointer(&h.object, unsafe.Pointer(old), unsafe.Pointer(new))
}

type alnCons struct {
	aln1, aln2 *sam.Alignment
	next       *alnCons
}

func (pair *samAlignmentPair) getOpticalDuplicates() *alnCons {
	return (*alnCons)(atomic.LoadPointer(&pair.opticalDuplicates))
}

func (pair *samAlignmentPair) addOpticalDuplicate(aln1, aln2 *sam.Alignment) {
	entry := &alnCons{aln1: aln1, aln2: aln2}
	for {
		head := atomic.LoadPointer(&pair.opticalDuplicates)
		entry.next = (*alnCons)(head)
		if atomic.CompareAndSwapPointer(&pair.opticalDuplicates, head, unsafe.Pointer(entry)) {
			break
		}
	}
}

// For each set of pairs with the same unclipped positions and
// directions, all except the one with the highest score are marked as
// duplicates.
//
// If multiple pairs are tied for best score, all except the one with
// the lexicographically smallest QNAME are marked as duplicates.
func classifyPair(aln *sam.Alignment, fragments, pairs *sync.Map) {
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
	aln1refid := aln1.REFID()
	aln2refid := aln2.REFID()
	aln1Pos := adaptedPos(aln1)
	aln2Pos := adaptedPos(aln2)
	if aln1refid > aln2refid ||
		(aln1refid == aln2refid && (aln1Pos > aln2Pos ||
			(aln1Pos == aln2Pos && aln1.IsReversed() && !aln2.IsReversed()))) {
		aln1, aln2 = aln2, aln1
		aln1refid, aln2refid = aln2refid, aln1refid
		aln1Pos, aln2Pos = aln2Pos, aln1Pos
	}
	entry, found := pairs.LoadOrStore(pair{
		aln1.LIBID(),
		aln1refid,
		aln2refid,
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
			np = &samAlignmentPair{score: score, aln1: aln1, aln2: aln2}
		}
		return np
	}

	for {
		if bestPair := best.pair(); bestPair.score > score {
			aln1.FLAG |= sam.Duplicate
			aln2.FLAG |= sam.Duplicate
			break
		} else if bestPair.score == score {
			if aln1.QNAME > bestPair.aln1.QNAME {
				aln1.FLAG |= sam.Duplicate
				aln2.FLAG |= sam.Duplicate
				break
			} else if best.compareAndSwapPair(bestPair, newPair()) {
				bestPair.aln1.FLAG |= sam.Duplicate
				bestPair.aln2.FLAG |= sam.Duplicate
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
// ties, the QNAME is used as a tie-breaker.
func MarkDuplicates(alsoOpticals bool) (sam.Filter, *sync.Map, *sync.Map) {
	splits := 16 * runtime.GOMAXPROCS(0)
	fragments := sync.NewMap(splits)
	pairsFragments := sync.NewMap(splits)
	pairs := sync.NewMap(splits)
	return func(header *sam.Header) sam.AlignmentFilter {
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
		if alsoOpticals {
			return func(aln *sam.Alignment) bool {
				addLIBID(aln, lbTable) // need this for all reads when marking optical duplicates
				if aln.FlagNotAny(sam.Unmapped | sam.Secondary | sam.Duplicate | sam.Supplementary) {
					adaptAlignment(aln)
					classifyFragment(aln, fragments)
					classifyPair(aln, pairsFragments, pairs)
				}
				return true
			}
		}
		return func(aln *sam.Alignment) bool {
			if aln.FlagNotAny(sam.Unmapped | sam.Secondary | sam.Duplicate | sam.Supplementary) {
				addLIBID(aln, lbTable) // don't need this for all reads when not marking optical duplicates
				adaptAlignment(aln)
				classifyFragment(aln, fragments)
				classifyPair(aln, pairsFragments, pairs)
			}
			return true
		}
	}, fragments, pairs
}
