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

package filters

import (
	"bytes"
	"encoding/gob"
	"fmt"
	"log"
	"math"
	"path"
	"runtime"
	"sort"
	"strings"
	osync "sync"
	"time"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/sam"
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

func computeTileInfo(aln *sam.Alignment) tileInfo {
	qnameInfo := strings.Split(aln.QNAME, ":")
	nColumns := len(qnameInfo)
	var t, x, y int64
	switch nColumns {
	case 7:
		t = internal.ParseInt(qnameInfo[4], 10, 64)
		x = internal.ParseInt(qnameInfo[5], 10, 64)
		y = internal.ParseInt(qnameInfo[6], 10, 64)
	case 5:
		t = internal.ParseInt(qnameInfo[2], 10, 64)
		x = internal.ParseInt(qnameInfo[3], 10, 64)
		y = internal.ParseInt(qnameInfo[4], 10, 64)
	default:
		unsupportedQnameFormatWarning := func() {
			log.Println("Warning: Unsupported qname format for extracting tile info. Cannot properly perform optical duplicate marking. ", len(qnameInfo), " columns.")
		}
		once.Do(unsupportedQnameFormatWarning)
		return tileInfo{-1, -1, -1}
	}
	return tileInfo{int(t), int(x), int(y)}
}

func (cache tileInfoCache) getTileInfo(aln *sam.Alignment) tileInfo {
	if tile, ok := cache[aln.QNAME]; ok {
		return tile
	}
	tile := computeTileInfo(aln)
	cache[aln.QNAME] = tile
	return tile
}

func isOpticalDuplicate(aln1 *sam.Alignment, tile1 tileInfo, aln2 *sam.Alignment, tile2 tileInfo, opticalPixelDistance int) bool {
	if aln1.RG() != aln2.RG() {
		return false
	}
	if tile1.t == -1 || tile2.t == -1 { // no tile info available
		return false
	}
	if tile1.t != tile2.t {
		return false
	}
	return isOpticalDuplicateShort(tile1, tile2, opticalPixelDistance)
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

// duplicatesCountHistograms keeps tracks of metrics for the number of pcr vs optical duplicates per list of duplicates
type duplicatesCountsHistograms struct {
	duplicatesCountHistogram           map[int]int
	nonOpticalDuplicatesCountHistogram map[int]int
	opticalDuplicatesCountHistogram    map[int]int
}

// merge histograms2 into histograms1
func mergeDuplicatesCountsHistograms(histograms1, histograms2 duplicatesCountsHistograms) (result duplicatesCountsHistograms) {
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
func incrementDuplicatesCountsHistograms(histograms duplicatesCountsHistograms, indices []int) {
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

func markOpticalDuplicatesPair(aln *sam.Alignment, pairFragments, pairs *sync.Map, ctr *DuplicatesCtr) {
	if isTruePair(aln) {
		aln1 := aln
		var aln2 *sam.Alignment
		if entry, deleted := pairFragments.DeleteOrStore(pairFragment{aln.LIBID(), aln.QNAME}, aln); deleted {
			aln2 = entry.(*sam.Alignment)
		} else {
			return
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
			log.Panicf("origin for duplicate read pair %v:%v unknown", aln1.LIBID(), aln1.QNAME)
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
	}
}

func fillGraphFromAGroup(tileCache tileInfoCache, duplicates []*sam.Alignment, group []int, opticalPixelDistance int, opticalDistanceRelationGraph graph) {
	for i, iIndex := range group {
		tileI := tileCache.getTileInfo(duplicates[iIndex])
		for j := i + 1; j < len(group); j++ {
			jIndex := group[j]
			tileJ := tileCache.getTileInfo(duplicates[jIndex])
			if isOpticalDuplicateShort(tileI, tileJ, opticalPixelDistance) {
				opticalDistanceRelationGraph.addEdge(iIndex, jIndex)
			}
		}
	}
}

type tileRGKey struct {
	rg interface{}
	t  int
}

func countOpticalDuplicatesWithGraph(duplicates []*sam.Alignment, opticalPixelDistance int) int {
	tileCache := make(tileInfoCache)

	opticalDistanceRelationGraph := newGraph(len(duplicates))
	tileRGMap := make(map[tileRGKey][]int)
	for i, aln := range duplicates {
		tile := tileCache.getTileInfo(aln)
		if tile.t != -1 {
			key := tileRGKey{aln.RG(), tile.t}
			tileRGMap[key] = append(tileRGMap[key], i)
		}
	}

	for _, tileGroup := range tileRGMap {
		if len(tileGroup) > 1 {
			fillGraphFromAGroup(tileCache, duplicates, tileGroup, opticalPixelDistance, opticalDistanceRelationGraph)
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
	return ctr
}

func countOpticalDuplicates(origin *samAlignmentPair, list *alnCons, opticalPixelDistance int) (int, []int) {
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
	parallel.Do(
		func() {
			forwardCount = countOpticalDuplicatesFromSlice(forwardDuplicates, opticalPixelDistance)
		},
		func() {
			reverseCount = countOpticalDuplicatesFromSlice(reverseDuplicates, opticalPixelDistance)
		},
	)
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
	return forwardCount + reverseCount, histogramIndices
}

func countOpticalDuplicatesFromSlice(duplicates []*sam.Alignment, opticalPixelDistance int) int {
	if len(duplicates) > 300000 {
		return 0
	}

	if len(duplicates) >= 4 {
		return countOpticalDuplicatesWithGraph(duplicates, opticalPixelDistance)
	}

	if len(duplicates) < 2 {
		return 0
	}

	tile0 := computeTileInfo(duplicates[0])
	tile1 := computeTileInfo(duplicates[1])

	var ctr int

	if isOpticalDuplicate(duplicates[0], tile0, duplicates[1], tile1, opticalPixelDistance) {
		ctr++
	}

	if len(duplicates) < 3 {
		return ctr
	}

	tile2 := computeTileInfo(duplicates[2])

	if isOpticalDuplicate(duplicates[0], tile0, duplicates[2], tile2, opticalPixelDistance) {
		ctr++
	}

	if ctr == 2 {
		return 2
	}

	if isOpticalDuplicate(duplicates[1], tile1, duplicates[2], tile2, opticalPixelDistance) {
		return ctr + 1
	}

	return ctr
}

type ctrsAndHistograms struct {
	ctrs       map[string]int
	histograms map[string]duplicatesCountsHistograms
}

func countOpticalDuplicatesPairs(pairs *sync.Map, opticalPixelDistance int) ctrsAndHistograms {
	result := pairs.ParallelReduce(
		func(alns map[interface{}]interface{}) interface{} {
			ctrs := make(map[string]int)
			histograms := make(map[string]duplicatesCountsHistograms)
			for _, value := range alns {
				origin := value.(*handle).pair()
				ctr, histogramIndices := countOpticalDuplicates(origin, origin.getOpticalDuplicates(), opticalPixelDistance)
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
					histogramsForLibID = duplicatesCountsHistograms{
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
			ctrsAndHistograms1 := x.(ctrsAndHistograms)
			ctrsAndHistograms2 := y.(ctrsAndHistograms)
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
	return result.(ctrsAndHistograms)
}

const undefinedLibrary = "Unknown Library"

func initDuplicatesCtrMap(header *sam.Header) map[string]*DuplicatesCtr {
	ctrs := make(map[string]*DuplicatesCtr)
	ctrs[undefinedLibrary] = &DuplicatesCtr{}
	for _, entry := range header.RG {
		library, found := entry["LB"]
		if found {
			ctrs[library] = &DuplicatesCtr{}
		}
	}
	return ctrs
}

func getDuplicatesCtr(aln *sam.Alignment, ctrs map[string]*DuplicatesCtr) *DuplicatesCtr {
	libid := aln.LIBID()
	if libid == nil {
		return ctrs[undefinedLibrary]
	}
	return ctrs[libid.(string)]
}

func mergeDuplicatesCtrMaps(ctrs1, ctrs2 map[string]*DuplicatesCtr) {
	for library, ctr2 := range ctrs2 {
		ctr1 := ctrs1[library]
		if ctr1 == nil {
			ctr1 = new(DuplicatesCtr)
			ctrs1[library] = ctr1
		}
		ctr1.UnpairedReadsExamined += ctr2.UnpairedReadsExamined
		ctr1.ReadPairsExamined += ctr2.ReadPairsExamined
		ctr1.SecondaryOrSupplementaryReads += ctr2.SecondaryOrSupplementaryReads
		ctr1.UnmappedReads += ctr2.UnmappedReads
		ctr1.UnpairedReadDuplicates += ctr2.UnpairedReadDuplicates
		ctr1.ReadPairDuplicates += ctr2.ReadPairDuplicates
		ctr1.ReadPairOpticalDuplicates += ctr2.ReadPairOpticalDuplicates
	}
}

// MarkOpticalDuplicates implements a function for calculating duplication metrics for a set of reads.
func MarkOpticalDuplicates(reads *sam.Sam, pairs *sync.Map, opticalPixelDistance int) map[string]*DuplicatesCtr {
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
				markOpticalDuplicatesPair(aln, pairsFragments, pairs, ctr)
			}
		}
		return ctrMap
	}, func(result1, result2 interface{}) interface{} {
		r1 := result1.(map[string]*DuplicatesCtr)
		r2 := result2.(map[string]*DuplicatesCtr)
		mergeDuplicatesCtrMaps(r1, r2)
		return r1
	})
	ctrMap := result.(map[string]*DuplicatesCtr)
	for _, ctr := range ctrMap {
		ctr.ReadPairsExamined = ctr.ReadPairsExamined / 2
	}
	// Now that for each "origin" we have the list of reads that are its duplicates, we check if among those duplicates are optical duplicates.
	//fnr := countOpticalDuplicatesFragments(fragments)
	pnr := countOpticalDuplicatesPairs(pairs, opticalPixelDistance)
	// Combine ctrs
	for library, nr := range pnr.ctrs {
		ctr := ctrMap[library]
		ctr.ReadPairOpticalDuplicates += nr
	}
	// Calculate derived metrics.
	// Fill in collected histograms
	for libID, ctr := range ctrMap {
		calculateDerivedDuplicateMetrics(ctr)
		histograms := pnr.histograms[libID]
		ctr.opticalDuplicatesCountHistogram = histograms.opticalDuplicatesCountHistogram
		ctr.duplicatesCountHistogram = histograms.duplicatesCountHistogram
		ctr.nonOpticalDuplicatesCountHistogram = histograms.nonOpticalDuplicatesCountHistogram
	}
	return ctrMap
}

func calculateDerivedDuplicateMetrics(ctr *DuplicatesCtr) {
	if ctr.ReadPairsExamined > 0 { // do not compute these metrics for single-end reads
		ctr.estimatedLibrarySize = estimateLibrarySize(ctr.ReadPairsExamined-ctr.ReadPairOpticalDuplicates, ctr.ReadPairsExamined-ctr.ReadPairDuplicates)
		ctr.histogram = histogramRoi(ctr)
	}
	ctr.percentDuplication = float64(ctr.UnpairedReadDuplicates+ctr.ReadPairDuplicates*2) / float64(ctr.UnpairedReadsExamined+ctr.ReadPairsExamined*2)
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
func PrintDuplicatesMetrics(metrics, commandLine string, ctrs map[string]*DuplicatesCtr) {
	file := internal.FileCreate(metrics)
	defer internal.Close(file)
	// Header
	fmt.Fprintln(file, "## htsjdk.samtools.metrics.StringHeader")
	fmt.Fprintln(file, "#", commandLine)
	fmt.Fprintln(file, "## htsjdk.samtools.metrics.StringHeader")
	fmt.Fprintln(file, "# Started on:", time.Now().Format("Mon Jan 02 15:04:05 MST 2006"))
	fmt.Fprintln(file)
	fmt.Fprintln(file, "## METRICS CLASS\tpicard.sam.DuplicationMetrics")

	// Metrics
	fmt.Fprintln(file, "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE")
	for library, ctr := range ctrs {
		if ctr.ReadPairsExamined > 0 {
			fmt.Fprintf(file, "%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%s\t%v\n", library, ctr.UnpairedReadsExamined, ctr.ReadPairsExamined, ctr.SecondaryOrSupplementaryReads, ctr.UnmappedReads, ctr.UnpairedReadDuplicates, ctr.ReadPairDuplicates, ctr.ReadPairOpticalDuplicates, formatFloat(ctr.percentDuplication), ctr.estimatedLibrarySize)
		} else {
			fmt.Fprintf(file, "%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%s\n", library, ctr.UnpairedReadsExamined, ctr.ReadPairsExamined, ctr.SecondaryOrSupplementaryReads, ctr.UnmappedReads, ctr.UnpairedReadDuplicates, ctr.ReadPairDuplicates, ctr.ReadPairOpticalDuplicates, formatFloat(ctr.percentDuplication))
		}
	}
	fmt.Fprintln(file)

	// Histogram:
	// only print a histogram if metrics exist for exactly one library
	// otherwise print nothing
	var ctr *DuplicatesCtr
	for _, c := range ctrs {
		if c.ReadPairsExamined > 0 {
			if ctr != nil {
				fmt.Fprintln(file)
				return
			}
			ctr = c
		}
	}
	if ctr == nil {
		fmt.Fprintln(file)
		return
	}
	fmt.Fprintln(file, "## HISTOGRAM\tjava.lang.Double")
	fmt.Fprintln(file, "BIN\tCoverageMult\tall_sets\toptical_sets\tnon_optical_sets")
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
	for k := range ctr.nonOpticalDuplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	for k := range ctr.opticalDuplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	for k := range ctr.duplicatesCountHistogram {
		if k > 100 {
			_, ok := histogramUnprintedKeysMap[k]
			if !ok {
				histogramUnprintedKeysMap[k] = 1
			}
		}
	}
	var histogramUnprintedKeys []int
	for k := range histogramUnprintedKeysMap {
		histogramUnprintedKeys = append(histogramUnprintedKeys, k)
	}
	sort.Ints(histogramUnprintedKeys)

	// plot the rest of the histogram counts
	for _, k := range histogramUnprintedKeys {
		fmt.Fprintf(file, "%v.0\t0\t%v\t%v\t%v\n", k,
			ctr.duplicatesCountHistogram[k],
			ctr.opticalDuplicatesCountHistogram[k],
			ctr.nonOpticalDuplicatesCountHistogram[k],
		)
	}
	fmt.Fprintln(file)
}

// PrintDuplicatesMetricsToIntermediateFile writes the duplicate metrics to a gob file.
func PrintDuplicatesMetricsToIntermediateFile(name string, ctrs map[string]*DuplicatesCtr) {
	file := internal.FileCreate(name)
	defer internal.Close(file)
	if err := gob.NewEncoder(file).Encode(ctrs); err != nil {
		log.Panic(err)
	}
}

// LoadAndCombineDuplicateMetrics loads partial duplication metrics from file and combines them
func LoadAndCombineDuplicateMetrics(metricsPath string) map[string]*DuplicatesCtr {
	// create ctr
	ctrs := make(map[string]*DuplicatesCtr)
	// go through the files, loading intermediate metrics
	metricsPath, files := internal.Directory(metricsPath)
	for _, fileName := range files {
		partialResult := make(map[string]*DuplicatesCtr)
		file := internal.FileOpen(path.Join(metricsPath, fileName))
		if err := gob.NewDecoder(file).Decode(&partialResult); err != nil {
			_ = file.Close()
			log.Panic(err)
		}
		internal.Close(file)
		// merge info
		mergeDuplicatesCtrMaps(ctrs, partialResult)
	}
	for _, ctr := range ctrs {
		calculateDerivedDuplicateMetrics(ctr)
	}
	return ctrs
}
