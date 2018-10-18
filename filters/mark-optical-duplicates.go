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
	"github.com/exascience/elprep/v4/utils"
	"github.com/exascience/pargo/parallel"
	"github.com/exascience/pargo/sync"
)

var optical = utils.Intern("optical")

func getOpticalDuplicate(aln *sam.Alignment) bool {
	_, ok := aln.Temps.Get(optical)
	return ok
}

func setOpticalDuplicate(aln *sam.Alignment) {
	aln.Temps.Set(optical, true)
}

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

const opticalPixelDistance = 100

func isOpticalDuplicate(aln1 *sam.Alignment, tile1 tileInfo, aln2 *sam.Alignment, tile2 tileInfo, deterministic bool) bool {
	if aln1.RG() != aln2.RG() {
		return false
	}
	if aln1.IsReversed() != aln2.IsReversed() {
		return false
	}
	if tile1.t == -1 || tile2.t == -1 { // no tile info available
		return false
	}
	if tile1.t != tile2.t {
		return false
	}
	if deterministic {
		return absInt(int(int16(tile1.x))-int(int16(tile2.x))) <= opticalPixelDistance && absInt(int(int16(tile1.y))-int(int16(tile2.y))) <= opticalPixelDistance
	}
	return absInt(tile1.x-tile2.x) <= opticalPixelDistance && absInt(tile1.y-tile2.y) <= opticalPixelDistance
}

// DuplicatesCtr implements a struct that stores metrics about reads such as the number of (optical) duplicates, unmapped reads, etc.
type DuplicatesCtr struct {
	UnpairedReadsExamined         int
	ReadPairsExamined             int
	SecondaryOrSupplementaryReads int
	UnmappedReads                 int
	UnpairedReadDuplicates        int
	ReadPairDuplicates            int
	ReadPairOpticalDuplicates     int
	percentDuplication            float64
	estimatedLibrarySize          int
	histogram                     []float64
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

func markOpticalDuplicatesFragment(aln *sam.Alignment, ctr *DuplicatesCtr) {
	if isTrueFragment(aln) {
		ctr.UnpairedReadDuplicates++
	}
}

func markOpticalDuplicatesPair(aln *sam.Alignment, pairFragments, pairs *sync.Map, ctr *DuplicatesCtr, deterministic bool) error {
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
		bestPair.addOpticalDuplicate(aln1, aln2) // map alns to origin duplicates
		// for aln2 same as other end, cause have same QNAME, so same tile info
		tile1, err := computeTileInfo(bestPair.aln1)
		if err != nil {
			return err
		}
		tile2, err := computeTileInfo(aln1)
		if err != nil {
			return err
		}
		if isOpticalDuplicate(bestPair.aln1, tile1, aln1, tile2, deterministic) {
			ctr.ReadPairOpticalDuplicates++
			setOpticalDuplicate(aln1)
			setOpticalDuplicate(aln2) // other end is duplicate too, cause have same QNAME, so same tile info
		}
		return nil
	}
	return nil
}

func pairOrientationLess(pair1, pair2 *alnCons) bool {
	if pair1.aln1.IsReversed() {
		if pair1.aln2.IsReversed() {
			return pair2.aln1.IsReversed() && !pair2.aln2.IsReversed()
		}
		return false
	}
	if pair1.aln2.IsReversed() {
		return pair2.aln1.IsReversed()
	}
	return pair2.aln1.IsReversed() || pair2.aln2.IsReversed()
}

func countOpticalDuplicates(list *alnCons, libraryTable map[string]int, deterministic bool) (int, error) {
	var pairs []*alnCons
	for entry := list; entry != nil; entry = entry.next {
		pairs = append(pairs, entry)
	}
	sort.Slice(pairs, func(i, j int) bool {
		pairI := pairs[i]
		pairJ := pairs[j]
		if libidI, libidJ := getLibIDIndex(pairI.aln1, libraryTable), getLibIDIndex(pairJ.aln1, libraryTable); libidI < libidJ {
			return true
		} else if libidI > libidJ {
			return false
		}
		if refidI, refidJ := pairI.aln1.REFID(), pairJ.aln1.REFID(); refidI < refidJ {
			return true
		} else if refidI > refidJ {
			return false
		}
		if pairIPos, pairJPos := adaptedPos(pairI.aln1), adaptedPos(pairJ.aln1); pairIPos < pairJPos {
			return true
		} else if pairIPos > pairJPos {
			return false
		}
		if pairOrientationLess(pairI, pairJ) {
			return true
		} else if pairOrientationLess(pairJ, pairI) {
			return false
		}
		if refidI, refidJ := pairI.aln2.REFID(), pairJ.aln2.REFID(); refidI < refidJ {
			return true
		} else if refidI > refidJ {
			return false
		}
		if pairIPos, pairJPos := adaptedPos(pairI.aln2), adaptedPos(pairJ.aln2); pairIPos < pairJPos {
			return true
		} else if pairIPos > pairJPos {
			return false
		}
		if pairI.aln1.POS < pairJ.aln1.POS {
			return true
		} else if pairI.aln1.POS > pairJ.aln1.POS {
			return false
		}
		if pairI.aln2.POS < pairJ.aln2.POS {
			return true
		} else if pairI.aln2.POS > pairJ.aln2.POS {
			return false
		}
		if !deterministic {
			return pairI.aln1.QNAME < pairJ.aln1.QNAME
		}
		if indexI, indexJ := pairI.aln1.FileIndex(), pairJ.aln1.FileIndex(); indexI < indexJ {
			return true
		} else if indexI > indexJ {
			return false
		}
		if indexI, indexJ := pairI.aln2.FileIndex(), pairJ.aln2.FileIndex(); indexI < indexJ {
			return true
		}
		return false
	})

	ctr := 0
	tileCache := make(tileInfoCache)
	for i := 0; i < len(pairs); i++ {
		opticalI := getOpticalDuplicate(pairs[i].aln1)
		tileI, err := tileCache.getTileInfo(pairs[i].aln1)
		if err != nil {
			return 0, err
		}
		for j := i + 1; j < len(pairs); j++ {
			if opticalJ := getOpticalDuplicate(pairs[j].aln1); !(opticalI && opticalJ) {
				tileJ, err := tileCache.getTileInfo(pairs[j].aln1)
				if err != nil {
					return 0, err
				}
				if isOpticalDuplicate(pairs[i].aln1, tileI, pairs[j].aln1, tileJ, deterministic) {
					ctr++
					if opticalJ {
						setOpticalDuplicate(pairs[i].aln1)
						setOpticalDuplicate(pairs[i].aln2)
					} else {
						setOpticalDuplicate(pairs[j].aln1)
						setOpticalDuplicate(pairs[j].aln2)
					}
				}
			}
		}
	}
	return ctr, nil
}

func countOpticalDuplicatesPairs(pairs *sync.Map, libraryTable map[string]int, deterministic bool) (map[string]int, error) {
	result := pairs.ParallelReduce(
		func(alns map[interface{}]interface{}) interface{} {
			ctrs := make(map[string]int)
			for _, value := range alns {
				origin := value.(*handle).pair()
				ctr, err := countOpticalDuplicates(origin.getOpticalDuplicates(), libraryTable, deterministic)
				if err != nil {
					return err
				}
				libID := origin.aln1.LIBID()
				if libID != nil {
					ctrs[libID.(string)] += ctr
				} else {
					ctrs[undefinedLibrary] += ctr
				}
			}
			return ctrs
		}, func(x, y interface{}) interface{} {
			var ctrs1, ctrs2 map[string]int
			switch xt := x.(type) {
			case error:
				return xt
			case map[string]int:
				ctrs1 = xt
			default:
				log.Fatal("invalid type during countOpticalDuplicatesPairs")
				panic("Unreachable code.")
			}
			switch yt := y.(type) {
			case error:
				return yt
			case map[string]int:
				ctrs2 = yt
			default:
				log.Fatal("invalid type during countOpticalDuplicatesPairs")
				panic("Unreachable code.")
			}
			if len(ctrs1) < len(ctrs2) {
				ctrs1, ctrs2 = ctrs2, ctrs1
			}
			for library, ctr := range ctrs2 {
				ctrs1[library] += ctr
			}
			return ctrs1
		})
	switch r := result.(type) {
	case error:
		return nil, r
	case map[string]int:
		return r, nil
	default:
		log.Fatal("invalid type during countOpticalDuplicatesPairs")
		panic("Unreachable code.")
	}
}

func initLibraryTable(header *sam.Header) map[string]int {
	table := make(map[string]int)
	ctr := 0
	for _, entry := range header.RG {
		library, found := entry["LB"]
		if found {
			table[library] = ctr
			ctr++
		}
	}
	return table
}

func getLibIDIndex(aln *sam.Alignment, libraryTable map[string]int) int {
	libid := aln.LIBID()
	if libid == nil {
		return -1
	}
	libIDIndex, found := libraryTable[libid.(string)]
	if !found {
		return -1
	}
	return libIDIndex
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

// MarkOpticalDuplicates implements a function for calculating duplication metrics for a set of reads
func MarkOpticalDuplicates(reads *sam.Sam, fragments, pairs *sync.Map, deterministic bool) DuplicatesCtrMap {
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
				if ctrMap.err = markOpticalDuplicatesPair(aln, pairsFragments, pairs, ctr, deterministic); ctrMap.err != nil {
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
	// We need to extract the order of the library ids from the header, since this is used for the second phase of optical duplicate marking.
	libraryTable := initLibraryTable(reads.Header)
	//fnr := countOpticalDuplicatesFragments(fragments)
	pnr, err := countOpticalDuplicatesPairs(pairs, libraryTable, deterministic)
	if err != nil {
		ctrMap.err = err
		return ctrMap
	}
	// Combine ctrs
	for library, nr := range pnr {
		ctr := ctrMap.Map[library]
		ctr.ReadPairOpticalDuplicates += nr
	}
	// Calculate derived metrics.
	for _, ctr := range ctrMap.Map {
		calculateDerivedDuplicateMetrics(ctr)
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
	fmt.Fprintln(file, "BIN\tVALUE")
	histogram := ctr.histogram
	for i := 0; i < len(histogram); i++ {
		fmt.Fprintf(file, "%v.0\t%s\n", i+1, formatFloat(histogram[i]))
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
	files, err := internal.Directory(metricsPath)
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
