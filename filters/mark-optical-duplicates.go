// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017-2019 imec vzw.

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
	"bytes"
	"encoding/gob"
	"fmt"
	"log"
	"math"
	"os"
	"path"
	"runtime"
	"sort"
	"strconv"
	"strings"
	osync "sync"
	"time"

	"github.com/exascience/elprep/v4/internal"
	"github.com/exascience/elprep/v4/sam"
	"github.com/exascience/pargo/parallel"
	"github.com/exascience/pargo/sync"
)

type (
	tileInfo struct {
		t, x, y int
	}

	tileInfoCache map[string]tileInfo
)

var once osync.Once

func computeTileInfo(aln *sam.Alignment) (info tileInfo, err error) {
	qnameInfo := strings.Split(aln.QNAME, ":")
	nColumns := len(qnameInfo)
	var t, x, y int64
	switch nColumns {
	case 7:
		t, err = strconv.ParseInt(qnameInfo[4], 10, 64)
		if err != nil {
			return
		}
		x, err = strconv.ParseInt(qnameInfo[5], 10, 64)
		if err != nil {
			return
		}
		y, err = strconv.ParseInt(qnameInfo[6], 10, 64)
		if err != nil {
			return
		}
	case 5:
		t, err = strconv.ParseInt(qnameInfo[2], 10, 64)
		if err != nil {
			return
		}
		x, err = strconv.ParseInt(qnameInfo[3], 10, 64)
		if err != nil {
			return
		}
		y, err = strconv.ParseInt(qnameInfo[4], 10, 64)
		if err != nil {
			return
		}
	default:
		unsupportedQnameFormatWarning := func() {
			log.Println("Warning: Unsupported qname format for extracting tile info. Cannot properly perform optical duplicate marking. ", len(qnameInfo), " columns.")
		}
		once.Do(unsupportedQnameFormatWarning)
		return tileInfo{-1, -1, -1}, nil
	}
	return tileInfo{int(t), int(x), int(y)}, nil
}

func (cache tileInfoCache) getTileInfo(aln *sam.Alignment) (tile tileInfo, err error) {
	if tile, ok := cache[aln.QNAME]; ok {
		return tile, nil
	}
	tile, err = computeTileInfo(aln)
	if err != nil {
		return
	}
	cache[aln.QNAME] = tile
	return tile, nil
}

func isOpticalDuplicateShort(tile1, tile2 tileInfo, deterministic bool, opticalPixelDistance int) bool {
	if deterministic {
		return absInt(int(int16(tile1.x))-int(int16(tile2.x))) <= opticalPixelDistance && absInt(int(int16(tile1.y))-int(int16(tile2.y))) <= opticalPixelDistance
	}
	return absInt(tile1.x-tile2.x) <= opticalPixelDistance && absInt(tile1.y-tile2.y) <= opticalPixelDistance
}

func isOpticalDuplicate(aln1 *sam.Alignment, tile1 tileInfo, aln2 *sam.Alignment, tile2 tileInfo, deterministic bool, opticalPixelDistance int) bool {
	if aln1.RG() != aln2.RG() {
		return false
	}
	if tile1.t == -1 || tile2.t == -1 { // no tile info available
		return false
	}
	if tile1.t != tile2.t {
		return false
	}
	return isOpticalDuplicateShort(tile1, tile2, deterministic, opticalPixelDistance)
}

// DuplicatesCtr implements a struct that stores metrics about reads such as the number of (optical) duplicates, unmapped reads, etc.
type DuplicatesCtr struct {
	UnpairedReadsExamined              int
	ReadPairsExamined                  int
	SecondaryOrSupplementaryReads      int
	UnmappedReads                      int
	UnpairedReadDuplicates             int
	ReadPairDuplicates                 int
	ReadPairOpticalDuplicates          int
	percentDuplication                 float64
	estimatedLibrarySize               int
	histogram                          []float64
	duplicatesCountHistogram           map[int]int
	nonOpticalDuplicatesCountHistogram map[int]int
	opticalDuplicatesCountHistogram    map[int]int
}

// DuplicatesCtrMap maps library names to duplicate counters.
type DuplicatesCtrMap struct {
	Map map[string]*DuplicatesCtr
	err error
}

// Err returns the error stored in this DuplicatesCtrMap.
func (ctrMap DuplicatesCtrMap) Err() error {
	return ctrMap.err
}

// DuplicatesCountHistograms keeps tracks of metrics for the number of pcr vs optical duplicates per list of duplicates
type DuplicatesCountsHistograms struct {
	duplicatesCountHistogram           map[int]int
	nonOpticalDuplicatesCountHistogram map[int]int
	opticalDuplicatesCountHistogram    map[int]int
}

// merge histograms2 into histograms1
func mergeDuplicatesCountsHistograms(histograms1, histograms2 DuplicatesCountsHistograms) (result DuplicatesCountsHistograms) {
	histogram1 := histograms1.duplicatesCountHistogram
	histogram2 := histograms2.duplicatesCountHistogram
	if len(histogram2) > len(histogram1) {
		histogram1, histogram2 = histogram2, histogram1
	}
	for i, v := range histogram2 {
		histogram1[i] += v
	}
	result.duplicatesCountHistogram = histogram1
	histogram1 = histograms1.nonOpticalDuplicatesCountHistogram
	histogram2 = histograms2.nonOpticalDuplicatesCountHistogram
	if len(histogram2) > len(histogram1) {
		histogram1, histogram2 = histogram2, histogram1
	}
	for i, v := range histogram2 {
		histogram1[i] += v
	}
	result.nonOpticalDuplicatesCountHistogram = histogram1
	histogram1 = histograms1.opticalDuplicatesCountHistogram
	histogram2 = histograms2.opticalDuplicatesCountHistogram
	if len(histogram2) > len(histogram1) {
		histogram1, histogram2 = histogram2, histogram1
	}
	for i, v := range histogram2 {
		histogram1[i] += v
	}
	result.opticalDuplicatesCountHistogram = histogram1
	return result
}

// indices is a slice of index1, index2, and index3
// index1 the length of the total number of duplicates for an origin
// index2 is the number of non optical duplicates in a list of duplicates for a given origin
// index3 is the number optical duplicates in a list of duplicates for a given origin
func incrementDuplicatesCountsHistograms(histograms DuplicatesCountsHistograms, indices []int) {
	for i, idx := range indices {
		histogram := histograms.duplicatesCountHistogram
		if i == 1 {
			if idx > 0 {
				histogram = histograms.nonOpticalDuplicatesCountHistogram
			} else {
				continue
			}
		}
		if i == 2 {
			if idx > 0 {
				histogram = histograms.opticalDuplicatesCountHistogram
			} else {
				continue
			}
		}
		histogram[idx] += 1
	}
}

func markOpticalDuplicatesFragment(aln *sam.Alignment, ctr *DuplicatesCtr) {
	if isTrueFragment(aln) {
		ctr.UnpairedReadDuplicates++
	}
}

func markOpticalDuplicatesPair(aln *sam.Alignment, pairFragments, pairs *sync.Map, ctr *DuplicatesCtr) error {
	if isTruePair(aln) {
		aln1 := aln
		var aln2 *sam.Alignment
		if entry, deleted := pairFragments.DeleteOrStore(pairFragment{aln.LIBID(), aln.QNAME}, aln); deleted {
			aln2 = entry.(*sam.Alignment)
		} else {
			return nil
		}
		ctr.ReadPairDuplicates++
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
		entry, found := pairs.Load(pair{
			aln1.LIBID(),
			aln1refid,
			aln2refid,
			(int64(aln1Pos) << 32) + int64(aln2Pos),
			aln1.IsReversed(),
			aln2.IsReversed(),
		})
		if !found {
			err := fmt.Errorf("origin for duplicate read pair %v:%v unknown", aln1.LIBID(), aln1.QNAME)
			return err
		}
		best := entry.(*handle)
		bestPair := best.pair()
		if bestPair.aln1 != aln1 {
			if aln1.IsFirst() {
				bestPair.addOpticalDuplicate(aln1) // map alns to origin duplicates
			} else {
				bestPair.addOpticalDuplicate(aln2)
			}
		}
		return nil
	}
	return nil
}

func fillGraphFromAGroup(tileCache tileInfoCache, duplicates []*sam.Alignment, group []int, opticalPixelDistance int, deterministic bool, opticalDistanceRelationGraph graph) error {
	for i, iIndex := range group {
		tileI, err := tileCache.getTileInfo(duplicates[iIndex])
		if err != nil {
			return err
		}
		for j := i + 1; j < len(group); j++ {
			jIndex := group[j]
			tileJ, err := tileCache.getTileInfo(duplicates[jIndex])
			if err != nil {
				return err
			}
			if isOpticalDuplicateShort(tileI, tileJ, deterministic, opticalPixelDistance) {
				opticalDistanceRelationGraph.addEdge(iIndex, jIndex)
			}
		}
	}
	return nil
}

type tileRGKey struct {
	rg interface{}
	t  int
}

func countOpticalDuplicatesWithGraph(duplicates []*sam.Alignment, deterministic bool, opticalPixelDistance int) (int, error) {
	tileCache := make(tileInfoCache)

	opticalDistanceRelationGraph := newGraph(len(duplicates))
	tileRGMap := make(map[tileRGKey][]int)
	for i, aln := range duplicates {
		tile, err := tileCache.getTileInfo(aln)
		if err != nil {
			return 0, err
		}
		if tile.t != -1 {
			key := tileRGKey{aln.RG(), tile.t}
			tileRGMap[key] = append(tileRGMap[key], i)
		}
	}

	for _, tileGroup := range tileRGMap {
		if len(tileGroup) > 1 {
			if err := fillGraphFromAGroup(tileCache, duplicates, tileGroup, opticalPixelDistance, deterministic, opticalDistanceRelationGraph); err != nil {
				return 0, err
			}
		}
	}

	var ctr int
	opticalDuplicateClusterMap := opticalDistanceRelationGraph.cluster()
	ctrPerCluster := make(map[int]int)
	for _, cluster := range opticalDuplicateClusterMap {
		ctrPerCluster[cluster] += 1
	}
	for _, clusterCtr := range ctrPerCluster {
		ctr += clusterCtr - 1
	}
	return ctr, nil
}

func countOpticalDuplicates(origin *samAlignmentPair, list *alnCons, deterministic bool, opticalPixelDistance int) (int, []int, error) {
	var forwardDuplicates, reverseDuplicates []*sam.Alignment
	var originAln *sam.Alignment
	if origin.aln1.IsFirst() {
		originAln = origin.aln1
	} else {
		originAln = origin.aln2
	}
	if originAln.IsReversed() {
		reverseDuplicates = append(reverseDuplicates, originAln)
	} else {
		forwardDuplicates = append(forwardDuplicates, originAln)
	}

	for entry := list; entry != nil; entry = entry.next {
		if entry.aln.IsReversed() {
			if len(reverseDuplicates) <= 300000 {
				reverseDuplicates = append(reverseDuplicates, entry.aln)
			}
		} else {
			if len(forwardDuplicates) <= 300000 {
				forwardDuplicates = append(forwardDuplicates, entry.aln)
			}
		}
	}

	var forwardCount, reverseCount int
	var forwardErr, reverseErr error
	parallel.Do(
		func() {
			forwardCount, forwardErr = countOpticalDuplicatesFromSlice(forwardDuplicates, deterministic, opticalPixelDistance)
		},
		func() {
			reverseCount, reverseErr = countOpticalDuplicatesFromSlice(reverseDuplicates, deterministic, opticalPixelDistance)
		},
	)
	if forwardErr != nil {
		return 0, nil, forwardErr
	}
	if reverseErr != nil {
		return 0, nil, reverseErr
	}
	// create histograms metrics counts
	opticalDuplicatesCount := forwardCount + reverseCount
	duplicatesCount := len(forwardDuplicates) + len(reverseDuplicates)

	index1 := duplicatesCount
	index2, index3 := 0, 0
	if duplicatesCount-opticalDuplicatesCount > 0 {
		index2 = duplicatesCount - opticalDuplicatesCount
	}
	if opticalDuplicatesCount > 0 {
		index3 = opticalDuplicatesCount + 1
	}

	histogramIndices := []int{index1, index2, index3}
	return forwardCount + reverseCount, histogramIndices, nil
}

func countOpticalDuplicatesFromSlice(duplicates []*sam.Alignment, deterministic bool, opticalPixelDistance int) (int, error) {
	if len(duplicates) > 300000 {
		return 0, nil
	}

	if len(duplicates) >= 4 {
		return countOpticalDuplicatesWithGraph(duplicates, deterministic, opticalPixelDistance)
	}

	if len(duplicates) < 2 {
		return 0, nil
	}

	tile0, err := computeTileInfo(duplicates[0])
	if err != nil {
		return 0, err
	}

	tile1, err := computeTileInfo(duplicates[1])
	if err != nil {
		return 0, err
	}

	var ctr int

	if isOpticalDuplicate(duplicates[0], tile0, duplicates[1], tile1, deterministic, opticalPixelDistance) {
		ctr++
	}

	if len(duplicates) < 3 {
		return ctr, nil
	}

	tile2, err := computeTileInfo(duplicates[2])
	if err != nil {
		return 0, err
	}

	if isOpticalDuplicate(duplicates[0], tile0, duplicates[2], tile2, deterministic, opticalPixelDistance) {
		ctr++
	}

	if ctr == 2 {
		return 2, nil
	}

	if isOpticalDuplicate(duplicates[1], tile1, duplicates[2], tile2, deterministic, opticalPixelDistance) {
		return ctr + 1, nil
	}

	return ctr, nil
}

type ctrsAndHistograms struct {
	ctrs       map[string]int
	histograms map[string]DuplicatesCountsHistograms
}

func countOpticalDuplicatesPairs(pairs *sync.Map, deterministic bool, opticalPixelDistance int) (ctrsAndHistograms, error) {
	result := pairs.ParallelReduce(
		func(alns map[interface{}]interface{}) interface{} {
			ctrs := make(map[string]int)
			histograms := make(map[string]DuplicatesCountsHistograms)
			for _, value := range alns {
				origin := value.(*handle).pair()
				ctr, histogramIndices, err := countOpticalDuplicates(origin, origin.getOpticalDuplicates(), deterministic, opticalPixelDistance)
				if err != nil {
					return err
				}
				libID := origin.aln1.LIBID()
				var libIDString string
				if libID != nil {
					libIDString = libID.(string)
					ctrs[libIDString] += ctr
				} else {
					libIDString = undefinedLibrary
					ctrs[undefinedLibrary] += ctr
				}
				// update histogram metrics
				histogramsForLibID, ok := histograms[libIDString]
				if !ok {
					histogramsForLibID = DuplicatesCountsHistograms{
						opticalDuplicatesCountHistogram:    make(map[int]int),
						nonOpticalDuplicatesCountHistogram: make(map[int]int),
						duplicatesCountHistogram:           make(map[int]int),
					}
					histograms[libIDString] = histogramsForLibID
				}
				incrementDuplicatesCountsHistograms(histogramsForLibID, histogramIndices)
			}
			return ctrsAndHistograms{ctrs, histograms}
		}, func(x, y interface{}) interface{} {
			var ctrsAndHistograms1, ctrsAndHistograms2 ctrsAndHistograms
			switch xt := x.(type) {
			case error:
				return xt
			case ctrsAndHistograms:
				ctrsAndHistograms1 = xt
			default:
				log.Fatal("invalid type during countOpticalDuplicatesPairs")
			}
			switch yt := y.(type) {
			case error:
				return yt
			case ctrsAndHistograms:
				ctrsAndHistograms2 = yt
			default:
				log.Fatal("invalid type during countOpticalDuplicatesPairs")
			}
			ctrs1 := ctrsAndHistograms1.ctrs
			ctrs2 := ctrsAndHistograms2.ctrs
			histograms1 := ctrsAndHistograms1.histograms
			histograms2 := ctrsAndHistograms2.histograms
			if len(ctrs1) < len(ctrs2) {
				ctrs1, ctrs2 = ctrs2, ctrs1
				histograms1, histograms2 = histograms2, histograms1
			}
			for library, ctr := range ctrs2 {
				ctrs1[library] += ctr
			}
			for library, histograms := range histograms2 {
				histograms1[library] = mergeDuplicatesCountsHistograms(histograms1[library], histograms)
			}
			ctrsAndHistograms1.ctrs = ctrs1
			ctrsAndHistograms1.histograms = histograms1
			return ctrsAndHistograms1
		})
	switch r := result.(type) {
	case error:
		return ctrsAndHistograms{}, r
	case ctrsAndHistograms:
		return r, nil
	default:
		log.Fatal("invalid type during countOpticalDuplicatesPairs")
		panic("Unreachable code.")
	}
}

const undefinedLibrary = "Unknown Library"

func initDuplicatesCtrMap(header *sam.Header) DuplicatesCtrMap {
	ctrs := DuplicatesCtrMap{Map: make(map[string]*DuplicatesCtr)}
	ctrs.Map[undefinedLibrary] = &DuplicatesCtr{}
	for _, entry := range header.RG {
		library, found := entry["LB"]
		if found {
			ctrs.Map[library] = &DuplicatesCtr{}
		}
	}
	return ctrs
}

func getDuplicatesCtr(aln *sam.Alignment, ctrs DuplicatesCtrMap) *DuplicatesCtr {
	libid := aln.LIBID()
	if libid == nil {
		return ctrs.Map[undefinedLibrary]
	}
	return ctrs.Map[libid.(string)]
}

func mergeDuplicatesCtrMaps(ctrs1, ctrs2 DuplicatesCtrMap) {
	for library, ctr2 := range ctrs2.Map {
		ctr1 := ctrs1.Map[library]
		if ctr1 == nil {
			ctr1 = new(DuplicatesCtr)
			ctrs1.Map[library] = ctr1
		}
		ctr1.UnpairedReadsExamined += ctr2.UnpairedReadsExamined
		ctr1.ReadPairsExamined += ctr2.ReadPairsExamined
		ctr1.SecondaryOrSupplementaryReads += ctr2.SecondaryOrSupplementaryReads
		ctr1.UnmappedReads += ctr2.UnmappedReads
		ctr1.UnpairedReadDuplicates += ctr2.UnpairedReadDuplicates
		ctr1.ReadPairDuplicates += ctr2.ReadPairDuplicates
		ctr1.ReadPairOpticalDuplicates += ctr2.ReadPairOpticalDuplicates
	}
	if ctrs1.err == nil {
		ctrs1.err = ctrs2.err
	}
}

// MarkOpticalDuplicates implements a function for calculating duplication metrics for a set of reads,
// optical pixel distance = 100
func MarkOpticalDuplicates(reads *sam.Sam, _, pairs *sync.Map, deterministic bool) DuplicatesCtrMap {
	return MarkOpticalDuplicatesWithPixelDistance(reads, pairs, deterministic, 100)
}

// MarkOpticalDuplicatesWithPixelDistance implements a function for calculating duplication metrics for a set of reads
func MarkOpticalDuplicatesWithPixelDistance(reads *sam.Sam, pairs *sync.Map, deterministic bool, opticalPixelDistance int) DuplicatesCtrMap {
	alns := reads.Alignments
	pairsFragments := sync.NewMap(16 * runtime.GOMAXPROCS(0))
	// Mark duplicates versus origin + collect for origins their duplicates.
	result := parallel.RangeReduce(0, len(alns), 0, func(low, high int) interface{} {
		ctrMap := initDuplicatesCtrMap(reads.Header)
		for _, aln := range alns[low:high] {
			ctr := getDuplicatesCtr(aln, ctrMap)
			if aln.IsUnmapped() {
				ctr.UnmappedReads++
				continue
			}
			if aln.FlagSome(sam.Secondary | sam.Supplementary) {
				ctr.SecondaryOrSupplementaryReads++
				continue
			}
			if isTrueFragment(aln) {
				ctr.UnpairedReadsExamined++
			}
			if isTruePair(aln) {
				ctr.ReadPairsExamined++
			}
			if aln.IsDuplicate() {
				markOpticalDuplicatesFragment(aln, ctr)
				if ctrMap.err = markOpticalDuplicatesPair(aln, pairsFragments, pairs, ctr); ctrMap.err != nil {
					return ctrMap
				}
			}
		}
		return ctrMap
	}, func(result1, result2 interface{}) interface{} {
		r1 := result1.(DuplicatesCtrMap)
		r2 := result2.(DuplicatesCtrMap)
		mergeDuplicatesCtrMaps(r1, r2)
		if r1.err == nil {
			r1.err = r2.err
		}
		return r1
	})
	ctrMap := result.(DuplicatesCtrMap)
	if ctrMap.err != nil {
		return ctrMap
	}
	for _, ctr := range ctrMap.Map {
		ctr.ReadPairsExamined = ctr.ReadPairsExamined / 2
	}
	// Now that for each "origin" we have the list of reads that are its duplicates, we check if among those duplicates are optical duplicates.
	//fnr := countOpticalDuplicatesFragments(fragments)
	pnr, err := countOpticalDuplicatesPairs(pairs, deterministic, opticalPixelDistance)
	if err != nil {
		ctrMap.err = err
		return ctrMap
	}
	// Combine ctrs
	for library, nr := range pnr.ctrs {
		ctr := ctrMap.Map[library]
		ctr.ReadPairOpticalDuplicates += nr
	}
	// Calculate derived metrics.
	// Fill in collected histograms
	for libID, ctr := range ctrMap.Map {
		calculateDerivedDuplicateMetrics(ctr)
		histograms := pnr.histograms[libID]
		ctr.opticalDuplicatesCountHistogram = histograms.opticalDuplicatesCountHistogram
		ctr.duplicatesCountHistogram = histograms.duplicatesCountHistogram
		ctr.nonOpticalDuplicatesCountHistogram = histograms.nonOpticalDuplicatesCountHistogram
	}
	return ctrMap
}

func calculateDerivedDuplicateMetrics(ctr *DuplicatesCtr) {
	ctr.estimatedLibrarySize = estimateLibrarySize(ctr.ReadPairsExamined-ctr.ReadPairOpticalDuplicates, ctr.ReadPairsExamined-ctr.ReadPairDuplicates)
	ctr.percentDuplication = float64(ctr.UnpairedReadDuplicates+ctr.ReadPairDuplicates*2) / float64(ctr.UnpairedReadsExamined+ctr.ReadPairsExamined*2)
	ctr.histogram = histogramRoi(ctr)
}

func f(x, c, n float64) float64 {
	return c/x - 1 + math.Exp(-n/x)
}

// Estimate the size of a library using the number of paired end molecules observed
// and the number of unique pairs observed
func estimateLibrarySize(nPairs, nUniquePairs int) int {
	n := float64(nPairs)
	c := float64(nUniquePairs)
	nReadPairDuplicates := nPairs - nUniquePairs
	if nPairs > 0 && nReadPairDuplicates > 0 {
		m := 1.0
		M := 100.0
		fd := f(M*c, c, n)
		for fd >= 0.0 {
			M *= 10.0
			fd = f(M*c, c, n)
		}
		for i := 0; i < 40; i++ {
			r := (m + M) / 2.0
			u := f(r*c, c, n)
			if u == 0.0 {
				break
			}
			if u > 0.0 {
				m = r
			}
			if u < 0.0 {
				M = r
			}
		}
		return int(c * ((m + M) / 2.0))
	}
	return 0
}

// Estimate the return on investment to be seen when a library was sequenced to x higher coverage than the observed coverage.
// Parameters: estimatedLibrarySize = nr of molecules in the library.
// x: the multiple of sequencing to be simulated.
// nPairs: the nr of pairs observed in actual sequencing.
// nUniquesPairs: the nr of uniques pairs, so without duplicate pairs.
// Returns a nr z <= x that if you sequenced pairs * x then you would get unique pairs = nUniquePairs * z.
func estimateRoi(estimatedLibrarySize, x, nPairs, nUniquePairs int) float64 {
	return float64(estimatedLibrarySize) * (1.0 - math.Exp(-float64(x*nPairs)/float64(estimatedLibrarySize))) / float64(nUniquePairs)
}

func histogramRoi(ctr *DuplicatesCtr) []float64 {
	histogram := make([]float64, 100)
	nUniquePairs := ctr.ReadPairsExamined - ctr.ReadPairDuplicates
	for x := 1; x <= 100; x++ {
		histogram[x-1] = estimateRoi(ctr.estimatedLibrarySize, x, ctr.ReadPairsExamined, nUniquePairs)
	}
	return histogram
}

func formatFloat(f float64) []byte {
	var s bytes.Buffer
	fmt.Fprintf(&s, "%.6f", f)
	b := s.Bytes()
	for i, c := range b {
		if c == '.' {
			for j := len(b) - 1; j > i; j-- {
				if b[j] != '0' {
					return b[:j+1]
				}
			}
			return b
		}
	}
	return b
}

// PrintDuplicatesMetrics writes the duplication metrics for a set of reads to a file.
func PrintDuplicatesMetrics(input, output, metrics string, removeDuplicates bool, ctrs DuplicatesCtrMap) (err error) {

	file, err := os.Create(metrics)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()
	// Header
	fmt.Fprintln(file, "## htsjdk.samtools.metrics.StringHeader")
	fmt.Fprintf(file, "# MarkDuplicates INPUT=[%v] OUTPUT=%v METRICS_FILE=%v    TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true REMOVE_DUPLICATES=%v ASSUME_SORTED=false DUPLICATE_SORTING_STRATEGY=SUM_OF_BASE_QUALITIES OPTICAL_DUPLICATE_PIXEL_DISTANCE=100\n", input, output, metrics, removeDuplicates)
	fmt.Fprintln(file, "## htsjdk.samtools.metrics.StringHeader")
	fmt.Fprintln(file, "# Started on:", time.Now().Format("Mon Jan 02 15:04:05 MST 2006"))
	fmt.Fprintln(file)
	fmt.Fprintln(file, "## METRICS CLASS\tpicard.sam.DuplicationMetrics")

	// Metrics
	fmt.Fprintln(file, "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE")
	for library, ctr := range ctrs.Map {
		if ctr.ReadPairsExamined > 0 {
			fmt.Fprintf(file, "%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%s\t%v\n", library, ctr.UnpairedReadsExamined, ctr.ReadPairsExamined, ctr.SecondaryOrSupplementaryReads, ctr.UnmappedReads, ctr.UnpairedReadDuplicates, ctr.ReadPairDuplicates, ctr.ReadPairOpticalDuplicates, formatFloat(ctr.percentDuplication), ctr.estimatedLibrarySize)
		}
	}
	fmt.Fprintln(file)

	// Histogram:
	// only print a histogram if metrics exist for exactly one library
	// otherwise print nothing
	var ctr *DuplicatesCtr
	for _, c := range ctrs.Map {
		if c.ReadPairsExamined > 0 {
			if ctr != nil {
				fmt.Fprintln(file)
				return nil
			}
			ctr = c
		}
	}
	if ctr == nil {
		fmt.Fprintln(file)
		return nil
	}
	fmt.Fprintln(file, "## HISTOGRAM\tjava.lang.Double")
	fmt.Fprintln(file, "BIN\tVALUE\tall_sets\toptical_sets\tnon_optical_sets")
	histogram := ctr.histogram
	for i := 0; i < len(histogram); i++ {
		fmt.Fprintf(file, "%v.0\t%s\t%v\t%v\t%v\n", i+1, formatFloat(histogram[i]),
			ctr.duplicatesCountHistogram[i+1],
			ctr.opticalDuplicatesCountHistogram[i+1],
			ctr.nonOpticalDuplicatesCountHistogram[i+1],
		)
	}

	// collect the keys > 100, so they can be sorted to be printed
	histogramUnprintedKeysMap := make(map[int]int)
	for k, _ := range ctr.nonOpticalDuplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	for k, _ := range ctr.opticalDuplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	for k, _ := range ctr.duplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	var histogramUnprintedKeys []int
	for k, _ := range histogramUnprintedKeysMap {
		histogramUnprintedKeys = append(histogramUnprintedKeys, k)
	}
	sort.Slice(histogramUnprintedKeys, func(i, j int) bool {
		return histogramUnprintedKeys[i] < histogramUnprintedKeys[j]
	})

	// plot the rest of the histogram counts
	for _, k := range histogramUnprintedKeys {
		fmt.Fprintf(file, "%v.0\t0\t%v\t%v\t%v\n", k,
			ctr.duplicatesCountHistogram[k],
			ctr.opticalDuplicatesCountHistogram[k],
			ctr.nonOpticalDuplicatesCountHistogram[k],
		)
	}
	fmt.Fprintln(file)

	return nil
}

// PrintDuplicatesMetricsToIntermediateFile writes the duplicate metrics to a gob file.
func PrintDuplicatesMetricsToIntermediateFile(name string, ctrs DuplicatesCtrMap) (err error) {
	file, err := os.Create(name)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()
	return gob.NewEncoder(file).Encode(ctrs)
}

// LoadAndCombineDuplicateMetrics loads partial duplication metrics from file and combines them
func LoadAndCombineDuplicateMetrics(metricsPath string) DuplicatesCtrMap {
	// create ctr
	ctrs := DuplicatesCtrMap{Map: make(map[string]*DuplicatesCtr)}
	// go through the files, loading intermediate metrics
	metricsPath, files, err := internal.Directory(metricsPath)
	if err != nil {
		ctrs.err = err
		return ctrs
	}
	for _, fileName := range files {
		partialResult := DuplicatesCtrMap{Map: make(map[string]*DuplicatesCtr)}
		file, err := os.Open(path.Join(metricsPath, fileName))
		if err != nil {
			ctrs.err = err
			return ctrs
		}
		if ctrs.err = gob.NewDecoder(file).Decode(&partialResult); ctrs.err != nil {
			_ = file.Close()
			return ctrs
		}
		if ctrs.err = file.Close(); ctrs.err != nil {
			return ctrs
		}
		// merge info
		mergeDuplicatesCtrMaps(ctrs, partialResult)
	}
	for _, ctr := range ctrs.Map {
		calculateDerivedDuplicateMetrics(ctr)
	}
	return ctrs
}
