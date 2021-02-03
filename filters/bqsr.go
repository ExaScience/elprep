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
	"log"
	"math"
	"sort"
	"sync"

	"github.com/exascience/elprep/v5/intervals"

	"github.com/exascience/pargo/parallel"

	"github.com/exascience/elprep/v5/fasta"
	"github.com/exascience/elprep/v5/sam"
)

func readGroupCovariate(hdr *sam.Header, aln *sam.Alignment) string {
	switch rg := aln.RG().(type) {
	default:
		log.Panic("Error: BQSR requires input with read groups. Alignment ", aln.QNAME, " has no read group. Please fix input, e.g. rerun with the --replace-read-group option.")
		panic("unreachable code")
	case string:
		for _, record := range hdr.RG {
			if id, ok := record["ID"]; ok && rg == id {
				if unit, ok := record["PU"]; ok {
					return unit
				}
				break
			}
		}
		return rg
	}
}

const lengthBits = 4

var simpleBaseToBaseIndexTable = map[byte]int32{'A': 0, 'a': 0, '*': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3}

func simpleBaseToBaseIndex(base byte) int32 {
	if index, ok := simpleBaseToBaseIndexTable[base]; ok {
		return index
	}
	return -1
}

func keyFromContext(dna []byte, start, end int) int32 {
	key := int32(end - start)
	var bitOffset uint = lengthBits
	for i := start; i < end; i++ {
		baseIndex := simpleBaseToBaseIndex(dna[i])
		if baseIndex == -1 {
			return -1
		}
		key |= baseIndex << bitOffset
		bitOffset += 2
	}
	return key
}

func createMask(contextSize int) (mask int32) {
	for i := 0; i < contextSize; i++ {
		mask = (mask << 2) | 3
	}
	return mask << lengthBits
}

const mismatchesContextSize = 2

func contextWith(bases []byte, contextSize int, keys []int32) []int32 {
	mask := createMask(contextSize)
	readLength := len(bases)
	keys = keys[:0]
	for i := 1; i < contextSize && i <= readLength; i++ {
		keys = append(keys, -1)
	}
	if readLength < contextSize {
		return keys
	}
	newBaseOffset := uint(2*(contextSize-1) + lengthBits)
	currentKey := keyFromContext(bases, 0, contextSize)
	keys = append(keys, currentKey)
	currentNPenalty := 0
	if currentKey == -1 {
		currentKey = 0
		currentNPenalty = contextSize - 1
		offset := newBaseOffset
		for bases[currentNPenalty] != 'N' {
			baseIndex := simpleBaseToBaseIndex(bases[currentNPenalty])
			currentKey |= baseIndex << offset
			offset -= 2
			currentNPenalty--
		}
	}
	contextSize32 := int32(contextSize)
	for currentIndex := contextSize; currentIndex < readLength; currentIndex++ {
		baseIndex := simpleBaseToBaseIndex(bases[currentIndex])
		if baseIndex == -1 {
			currentNPenalty = contextSize
			currentKey = 0
		} else {
			currentKey = (currentKey >> 2) & mask
			currentKey |= baseIndex << newBaseOffset
			currentKey |= contextSize32
		}
		if currentNPenalty == 0 {
			keys = append(keys, currentKey)
		} else {
			currentNPenalty--
			keys = append(keys, -1)
		}
	}
	return keys
}

func reverseInt32Slice(slice []int32) []int32 {
	for i, j := 0, len(slice)-1; i < j; i, j = i+1, j-1 {
		slice[i], slice[j] = slice[j], slice[i]
	}
	return slice
}

func computeBaseContextCovariate(aln *sam.Alignment, strandedClippedSeq []byte, keys []int32) []int32 {
	keys = contextWith(strandedClippedSeq, mismatchesContextSize, keys)
	if aln.IsReversed() {
		return reverseInt32Slice(keys)
	}
	return keys
}

func simpleBaseIndexToBase(baseIndex int32) (base byte) {
	switch baseIndex {
	case -1:
		base = 'N'
	case 0:
		base = 'A'
	case 1:
		base = 'C'
	case 2:
		base = 'G'
	case 3:
		base = 'T'
	default:
		log.Panic("invalid baseIndex")
	}
	return
}

func keyToString(key int32) string {
	length := int(key & 0xF)
	realKey := key >> 4
	result := make([]byte, 0, 4)
	for i := 0; i < length; i++ {
		result = append(result, simpleBaseIndexToBase(realKey&0x3))
		realKey >>= 2
	}
	if realKey != 0 {
		log.Panic("invalid key")
	}
	return string(result)
}

type (
	bqsrTableKey struct {
		Qual      uint8
		Covariate int32
		ReadGroup string
	}

	bqsrEntry struct {
		EmpiricalQuality uint8
		Observations     int
		Mismatches       int
	}

	bqsrTable map[bqsrTableKey]*bqsrEntry

	combinedBqsrEntry struct {
		reportedQuality float64
		bqsrEntry
	}
)

func (b bqsrTable) update(key bqsrTableKey, qual uint8, errorVal int) {
	if entry, ok := b[key]; ok {
		entry.Observations++
		entry.Mismatches += errorVal
		return
	}
	b[key] = &bqsrEntry{EmpiricalQuality: qual, Observations: 1, Mismatches: errorVal}
}

func (b bqsrTable) merge(c bqsrTable) bqsrTable {
	if len(b) < len(c) {
		b, c = c, b
	}
	for key, value := range c {
		if entry, ok := b[key]; ok {
			entry.Observations += value.Observations
			entry.Mismatches += value.Mismatches
		} else {
			b[key] = value
		}
	}
	return b
}

func recalibrateAln(hdr *sam.Header, aln *sam.Alignment) bool {
	_, found := aln.TAGS.Get(sr)
	if found {
		return false
	}
	if aln.MAPQ > 0 && aln.MAPQ < 255 &&
		aln.FlagNotAny(sam.Secondary|sam.Duplicate|sam.QCFailed) &&
		!isStrictUnmapped(aln) &&
		aln.POS > 0 &&
		aln.SEQ.Len() > 0 &&
		aln.SEQ.Len() == len(aln.QUAL) &&
		aln.RG() != nil &&
		alignmentAgreesWithHeader(hdr, aln) {
		cigarVec := aln.CIGAR
		return !cigarContainsN(cigarVec) &&
			sam.ReferenceLengthFromCigar(aln.CIGAR) >= 0 &&
			aln.SEQ.Len() == int(sam.ReadLengthFromCigar(cigarVec))
	}
	return false
}

// In the following table, 'N' and all IUPAC codes are 0.
var baseToIntMap = map[byte]int{
	'a': 1, 'A': 1, '*': 1,
	'c': 2, 'C': 2,
	'g': 3, 'G': 3,
	't': 4, 'T': 4,
}

func computeSnpEvents(aln *sam.Alignment, ref []byte, snps []int) []int {
	read := aln.SEQ
	readLength := read.Len()
	for cap(snps) < readLength {
		snps = append(snps[:cap(snps)], 0)
	}
	snps = snps[:readLength]
	for i := range snps {
		snps[i] = 0
	}
	i := 0
	j := int(aln.POS - 1)
	cigarVec := aln.CIGAR
	for _, cigarOp := range cigarVec {
		ln := int(cigarOp.Length)
		switch cigarOp.Operation {
		case 'M', '=', 'X':
			for k := 0; k < ln; k++ {
				if baseToIntMap[read.Base(i)] != baseToIntMap[ref[j]] {
					snps[i] = 1
				}
				i++
				j++
			}
		case 'D', 'N':
			j += ln
		case 'I', 'S':
			i += ln
		}
	}
	return snps
}

func readStartsWithInsertion(cigarVec []sam.CigarOperation) (int32, bool) {
	for _, element := range cigarVec {
		switch element.Operation {
		case 'I':
			return element.Length, true
		case 'H', 'S':
			continue
		default:
			return -1, false
		}
	}
	return -1, false
}

const lowQualityTail = 2

var baseComplementTable = map[byte]byte{'A': 'T', 'a': 'T', 'C': 'G', 'c': 'G', 'G': 'C', 'g': 'C', 'T': 'A', 't': 'A'}

func baseComplement(base byte) byte {
	if complement, ok := baseComplementTable[base]; ok {
		return complement
	}
	return base
}

func computeStrandedClippedSeq(aln *sam.Alignment, newSeq []byte) []byte {
	seq := aln.SEQ
	seqLength := seq.Len()
	leftPos := seqLength
	for i := 0; i < leftPos; i++ {
		if aln.QUAL[i] > lowQualityTail {
			leftPos = i
			break
		}
	}
	rightPos := leftPos - 1
	for i := seqLength - 1; i >= leftPos; i-- {
		if aln.QUAL[i] > lowQualityTail {
			rightPos = i
			break
		}
	}
	if leftPos > rightPos {
		return nil
	}
	for cap(newSeq) < seqLength {
		newSeq = append(newSeq[:cap(newSeq)], 0)
	}
	newSeq = newSeq[:seqLength]
	if aln.IsReversed() {
		j := -1
		for i := rightPos + 1; i < seqLength; i++ {
			j++
			newSeq[j] = 'N'
		}
		for i := rightPos; i >= leftPos; i-- {
			j++
			newSeq[j] = baseComplement(seq.Base(i))
		}
		for i := 0; i < leftPos; i++ {
			j++
			newSeq[j] = 'N'
		}
	} else {
		for i := 0; i < leftPos; i++ {
			newSeq[i] = 'N'
		}
		for i := leftPos; i <= rightPos; i++ {
			newSeq[i] = seq.Base(i)
		}
		for i := rightPos + 1; i < seqLength; i++ {
			newSeq[i] = 'N'
		}
	}
	return newSeq
}

func checkCycleCovariate(cycle, maxCycle int) int {
	if cycle > maxCycle || cycle < -maxCycle {
		log.Panic("cycle value exceeds maximum cycle value")
	}
	return cycle
}

const (
	reversedShift = 4
	lastShift     = 7
)

func prepareCycleCovariates(aln *sam.Alignment) (int, int) {
	reversed := int((aln.FLAG & sam.Reversed) >> reversedShift)
	last := int((aln.FLAG & sam.Last) >> lastShift)
	readOrderFactor := 1 - 2*last
	cycleFactor := readOrderFactor + reversed*(aln.SEQ.Len()-1)*readOrderFactor
	increment := (1 - 2*reversed) * readOrderFactor
	return cycleFactor, increment
}

func computeBaseCycleCovariate(cycleFactor, increment, index, maxCycle int) int {
	return checkCycleCovariate(cycleFactor+index*increment, maxCycle)
}

func calculateSkipSlice(aln *sam.Alignment, knownSites []intervals.Interval, skipSlice []bool) []bool {
	seqLength := aln.SEQ.Len()
	softStart := softStart(aln)
	softEnd := softEnd(aln)
	for cap(skipSlice) < seqLength {
		skipSlice = append(skipSlice[:cap(skipSlice)], false)
	}
	skipSlice = skipSlice[:seqLength]
	for i := range skipSlice {
		skipSlice[i] = false
	}
	for _, site := range intervals.Intersect(knownSites, int32(softStart), int32(softEnd)) {
		featureStartOnRead, ok := getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart, int(site.Start), left)
		if !ok || featureStartOnRead < 0 {
			featureStartOnRead = 0
		}
		featureEndOnRead, ok := getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart, int(site.End), left)
		if !ok || featureEndOnRead > seqLength-1 {
			featureEndOnRead = seqLength - 1
		}
		for i := featureStartOnRead; i <= featureEndOnRead; i++ {
			skipSlice[i] = true
		}
	}
	return skipSlice
}

// BaseRecalibrator implements the first step of base recalibration.
type BaseRecalibrator struct {
	forKnownIntervals sync.WaitGroup
	fastaMap          *fasta.MappedFasta
	knownIntervals    map[string][]intervals.Interval
}

// NewBaseRecalibrator returns a struct for the first step of base recalibration.
func NewBaseRecalibrator(knownSites []string, referenceFasta *fasta.MappedFasta) (recal *BaseRecalibrator) {
	recal = &BaseRecalibrator{
		fastaMap:       referenceFasta,
		knownIntervals: make(map[string][]intervals.Interval),
	}
	recal.forKnownIntervals.Add(1)
	go func() {
		defer recal.forKnownIntervals.Done()
		for _, knownSitesFile := range knownSites {
			for chrom, ivals := range intervals.FromElsitesFile(knownSitesFile) {
				recal.knownIntervals[chrom] = append(recal.knownIntervals[chrom], ivals...)
			}
		}
		for chrom, ivals := range recal.knownIntervals {
			intervals.ParallelSortByStart(ivals)
			recal.knownIntervals[chrom] = intervals.ParallelFlatten(ivals)
		}
	}()
	return recal
}

// BaseRecalibratorTables is the result of the base recalibration.
// All subsequent steps, including ApplyBQSR, are based on these tables.
type BaseRecalibratorTables struct {
	QualityScores, Cycles, Contexts bqsrTable
}

// NewBaseRecalibratorTables returns a struct for storing the result
// of the base recalibration.
func NewBaseRecalibratorTables() BaseRecalibratorTables {
	return BaseRecalibratorTables{
		QualityScores: make(map[bqsrTableKey]*bqsrEntry),
		Cycles:        make(map[bqsrTableKey]*bqsrEntry),
		Contexts:      make(map[bqsrTableKey]*bqsrEntry),
	}
}

func (recal *BaseRecalibrator) close() {
	recal.forKnownIntervals.Wait()
	recal.knownIntervals = nil
}

// Recalibrate implements the first step of base recalibration.
func (recal *BaseRecalibrator) Recalibrate(reads *sam.Sam, maxCycle int) *BaseRecalibratorTables {
	defer recal.close()
	hdr := reads.Header
	alns := reads.Alignments
	result := parallel.RangeReduce(0, len(alns), 0, func(low, high int) interface{} {
		result := NewBaseRecalibratorTables()
		var snps []int
		var strandedClippedSeq []byte
		var baseContextCovariate []int32
		var skipSlice []bool
		aln := new(sam.Alignment)
		for _, alignment := range alns[low:high] {
			*aln = *alignment
			if !recalibrateAln(hdr, aln) {
				continue
			}
			hardClipAdaptorSequence(aln)
			if aln.SEQ.Len() == 0 {
				continue
			}
			hardClipSoftClippedBases(aln)
			if aln.SEQ.Len() == 0 {
				continue
			}
			var waitForSkipSlice sync.WaitGroup
			waitForSkipSlice.Add(1)
			go func() {
				defer waitForSkipSlice.Done()
				recal.forKnownIntervals.Wait()
				skipSlice = calculateSkipSlice(aln, recal.knownIntervals[aln.RNAME], skipSlice)
			}()
			ref := recal.fastaMap.Seq(aln.RNAME)
			snps = computeSnpEvents(aln, ref, snps)
			readGroupCovariate := readGroupCovariate(hdr, aln)
			cycleFactor, cycleIncrement := prepareCycleCovariates(aln)
			strandedClippedSeq = computeStrandedClippedSeq(aln, strandedClippedSeq)
			baseContextCovariate = computeBaseContextCovariate(aln, strandedClippedSeq, baseContextCovariate)
			waitForSkipSlice.Wait()
			for i, l := 0, aln.SEQ.Len(); i < l; i++ {
				if skipSlice[i] {
					continue
				}
				if _, ok := simpleBaseToBaseIndexTable[aln.SEQ.Base(i)]; !ok {
					continue
				}
				qual := aln.QUAL[i]
				if qual < minInterestingQual {
					continue
				}
				// mismatch
				errVal := snps[i]
				key1 := bqsrTableKey{
					ReadGroup: readGroupCovariate,
					Qual:      qual,
				}
				result.QualityScores.update(key1, qual, errVal)
				baseCycleCovariate := computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, maxCycle)
				key2 := bqsrTableKey{
					Covariate: int32(baseCycleCovariate),
					ReadGroup: readGroupCovariate,
					Qual:      qual,
				}
				result.Cycles.update(key2, qual, errVal)
				if len(baseContextCovariate) > 0 && baseContextCovariate[i] >= 0 {
					key3 := bqsrTableKey{
						Covariate: baseContextCovariate[i],
						ReadGroup: readGroupCovariate,
						Qual:      qual,
					}
					result.Contexts.update(key3, qual, errVal)
				}
			}
		}
		return result
	}, func(result1, result2 interface{}) interface{} {
		r1 := result1.(BaseRecalibratorTables)
		r2 := result2.(BaseRecalibratorTables)
		r1.QualityScores = r1.QualityScores.merge(r2.QualityScores)
		r1.Cycles = r1.Cycles.merge(r2.Cycles)
		r1.Contexts = r1.Contexts.merge(r2.Contexts)
		return r1
	})
	tables := result.(BaseRecalibratorTables)
	return &tables
}

const (
	smoothingConstant        = 1
	maxQualityScore          = 93
	maxRecalibratedQualScore = 93
	maxReasonableQualScore   = 60
	maxNumberOfObservations  = math.MaxInt32 - 1
)

func qualityToErrorProbability(phred float64) float64 {
	return math.Pow(10, phred/-10)
}

func qualityToProbability(phred float64) float64 {
	return 1 - math.Pow(10, phred/-10)
}

var log10QualEmpiricalPriorCache = [...]float64{
	-0.045757490560675115,
	-0.9143464543671788,
	-3.5201133457866898,
	-7.863058164819208,
	-13.943180911464733,
	-21.760481585723266,
	-31.314960187594806,
	-42.606616717079355,
	-55.63545117417691,
	-70.40146355888747,
	-86.90465387121104,
	-105.14502211114761,
	-125.1225682786972,
	-146.83729237385978,
	-170.2891943966354,
	-195.47827434702398,
	-222.4045322250256,
	-251.06796803064023,
	-281.46858176386786,
	-313.60637342472336,
	-1.7976931348623157e308,
}

func log10QualEmpiricalPrior(empiricalQual, reportedQual float64) float64 {
	diff := minInt(absInt(int(empiricalQual-reportedQual)), len(log10QualEmpiricalPriorCache)-1)
	return log10QualEmpiricalPriorCache[diff]
}

func log10Gamma(n int) float64 {
	gammaN, _ := math.Lgamma(float64(n))
	return gammaN * math.Log10E
}

func log10BinomialCoefficient(n, k int) float64 {
	return log10Gamma(n+1) - log10Gamma(k+1) - log10Gamma(n-k+1)
}

func log10BinomialProbability(n, k int, log10p float64) float64 {
	if log10p == 0.0 {
		return -math.MaxFloat64
	}
	log10MinP := log10(1.0 - math.Pow(10, log10p))
	return log10BinomialCoefficient(n, k) + log10p*float64(k) + log10MinP*float64(n-k)
}

func log10QualEmpiricalLikelihood(empiricalQual float64, observations, mismatches int) float64 {
	if observations == 0 {
		return 0.0
	}
	qualToErrorProbLog10 := empiricalQual / -10.0
	return log10BinomialProbability(observations, mismatches, qualToErrorProbLog10)
}

func calculateBayesianEstimateOfEmpiricalQuality(observations, mismatches int, priorQuality float64) uint8 {
	if observations > maxNumberOfObservations {
		mismatches = int(math.Round(float64(mismatches) * (float64(maxNumberOfObservations) / float64(observations))))
		observations = maxNumberOfObservations
	}
	nbins := maxReasonableQualScore + 1
	max := -math.MaxFloat64
	var maxI uint8
	for i := 0; i < nbins; i++ {
		floati := float64(i)
		p1 := log10QualEmpiricalPrior(floati, priorQuality)
		p2 := log10QualEmpiricalLikelihood(floati, observations, mismatches)
		log10Posterior := p1 + p2
		if max < log10Posterior {
			max = log10Posterior
			maxI = uint8(i)
		}
	}
	return maxI
}

func (entry *bqsrEntry) calculateEmpiricalQuality(priorQual float64) uint8 {
	smoothMismatches := entry.Mismatches + smoothingConstant
	smoothObservations := entry.Observations + smoothingConstant + smoothingConstant
	empiricalQual := calculateBayesianEstimateOfEmpiricalQuality(smoothObservations, smoothMismatches, priorQual)
	return minUint8(empiricalQual, maxRecalibratedQualScore)
}

func calculateExpectedErrors(observations int, reportedQuality float64) float64 {
	return float64(observations) * qualityToErrorProbability(reportedQuality)
}

func initializeCombinedBQSRTable(bqsrTable bqsrTable) map[string]*combinedBqsrEntry {
	table := make(map[string]*combinedBqsrEntry)
	for key, entry := range bqsrTable {
		if tableEntry, ok := table[key.ReadGroup]; ok {
			sumErrors := calculateExpectedErrors(tableEntry.Observations, tableEntry.reportedQuality) + calculateExpectedErrors(entry.Observations, float64(key.Qual))
			tableEntry.Observations += entry.Observations
			tableEntry.Mismatches += entry.Mismatches
			tableEntry.reportedQuality = -10 * log10(sumErrors/float64(tableEntry.Observations))
		} else {
			table[key.ReadGroup] = &combinedBqsrEntry{
				reportedQuality: float64(key.Qual),
				bqsrEntry:       *entry,
			}
		}
	}
	for _, tableEntry := range table {
		tableEntry.EmpiricalQuality = tableEntry.calculateEmpiricalQuality(tableEntry.reportedQuality)
	}
	return table
}

// FinalizeBQSRTables finalizes the first step of base recalibration.
func (recal *BaseRecalibratorTables) FinalizeBQSRTables() {
	parallel.Do(
		func() {
			for key, entry := range recal.QualityScores {
				entry.EmpiricalQuality = entry.calculateEmpiricalQuality(float64(key.Qual))
			}
		},
		func() {
			for key, entry := range recal.Cycles {
				entry.EmpiricalQuality = entry.calculateEmpiricalQuality(float64(key.Qual))
			}
		},
		func() {
			for key, entry := range recal.Contexts {
				entry.EmpiricalQuality = entry.calculateEmpiricalQuality(float64(key.Qual))
			}
		})
}

const (
	globalQualityScorePrior = -1
	minInterestingQual      = 6
)

func errorProbabilityToQuality(prob float64) int {
	if prob == 0.0 {
		return maxQualityScore
	}
	return maxInt(minInt(int(math.Round(-10*log10(prob))), maxQualityScore), 1)
}

const maxStaticQuantizedQuality = 254

func initializeStaticQuantizedScores(quals []uint8) []uint8 {
	staticScores := make([]uint8, maxStaticQuantizedQuality)
	for i := 0; i < minInterestingQual; i++ {
		staticScores[i] = uint8(i)
	}
	if len(quals) == 1 {
		qual := quals[0]
		for i := minInterestingQual; i < maxStaticQuantizedQuality; i++ {
			staticScores[i] = qual
		}
		return staticScores
	}
	sort.Slice(quals, func(i, j int) bool { return quals[i] < quals[j] })
	firstQual := uint8(minInterestingQual)
	prevQual := firstQual
	prevProb := qualityToProbability(float64(prevQual))
	for _, nextQual := range quals {
		for i := prevQual; i < nextQual; i++ {
			nextProb := qualityToProbability(float64(nextQual))
			iProb := qualityToProbability(float64(i))
			if iProb-prevProb > nextProb-iProb {
				staticScores[i] = nextQual
			} else {
				staticScores[i] = prevQual
			}
			prevProb = nextProb
			prevQual = nextQual
		}
	}
	for i := prevQual; i < maxStaticQuantizedQuality; i++ {
		staticScores[i] = prevQual
	}
	return staticScores
}

type quantizationInterval struct {
	next      int
	errorRate float64
	nobs      int
	leafNobs  int
	nerrors   int
}

func leafInterval(index int, interval *quantizationInterval) bool {
	if interval.next < 0 {
		return index == maxQualityScore
	}
	return interval.next == index+1
}

func initializeQuantizationIntervals(quantizationMap []int) []quantizationInterval {
	intervals := make([]quantizationInterval, len(quantizationMap))
	for i, nobs := range quantizationMap {
		errorRate := qualityToErrorProbability(float64(i))
		nerrors := int(float64(nobs) * errorRate)
		j := i + 1
		if j == len(quantizationMap) {
			j = -1
		}
		intervals[i] = quantizationInterval{
			next:      j,
			errorRate: errorRate,
			nobs:      nobs,
			leafNobs:  nobs,
			nerrors:   nerrors,
		}
	}
	return intervals
}

func leafPenalty(k int, quantizationIntervals []quantizationInterval, globalErrorRate float64) float64 {
	interval := &quantizationIntervals[k]
	if k <= minInterestingQual {
		return 0.0
	}
	return math.Abs(log10(interval.errorRate)-log10(globalErrorRate)) * float64(interval.leafNobs)
}

func calculateErrorRate(nobs, nerrors int) float64 {
	if nobs == 0 {
		return 0.0
	}
	return float64(nerrors+1) / float64(nobs+1)
}

func computeMergePenalty(i, j int, quantizationIntervals []quantizationInterval) float64 {
	intervalI := &quantizationIntervals[i]
	intervalJ := &quantizationIntervals[j]
	mergedNobs := intervalI.nobs + intervalJ.nobs
	mergedNerrors := intervalI.nerrors + intervalJ.nerrors
	mergedErrorRate := calculateErrorRate(mergedNobs, mergedNerrors)
	if mergedErrorRate == 0 {
		return 0.0
	}
	var sumI, sumJ float64
	for k := i; k < j; k++ {
		sumI += leafPenalty(k, quantizationIntervals, mergedErrorRate)
	}
	var kend int
	if intervalJ.next >= 0 {
		kend = intervalJ.next
	} else {
		kend = len(quantizationIntervals)
	}
	for k := j; k < kend; k++ {
		sumJ += leafPenalty(k, quantizationIntervals, mergedErrorRate)
	}
	return sumI + sumJ
}

func mergeMinimalPenaltyQuantizationIntervals(quantizationIntervals []quantizationInterval) bool {
	i := 0
	intervalI := &quantizationIntervals[0]
	j := intervalI.next
	if j < 0 {
		return false
	}
	minI := i
	mergePenalty := computeMergePenalty(i, j, quantizationIntervals)
	for {
		i = j
		intervalI = &quantizationIntervals[i]
		j = intervalI.next
		if j < 0 {
			break
		}
		penalty := computeMergePenalty(i, j, quantizationIntervals)
		if penalty < mergePenalty {
			minI = i
			mergePenalty = penalty
		}
	}
	intervalMinI := &quantizationIntervals[minI]
	intervalMinJ := &quantizationIntervals[intervalMinI.next]
	mergedNobs := intervalMinI.nobs + intervalMinJ.nobs
	mergedNerrors := intervalMinI.nerrors + intervalMinJ.nerrors
	intervalMinI.next = intervalMinJ.next
	intervalMinI.nobs = mergedNobs
	intervalMinI.nerrors = mergedNerrors
	return true
}

func mergeQuantizationIntervals(quantizationIntervals []quantizationInterval, levels int) {
	n := len(quantizationIntervals)
	for n > levels {
		if mergeMinimalPenaltyQuantizationIntervals(quantizationIntervals) {
			n--
		} else {
			break
		}
	}
}

func initializeQuantizedQualityScores(bqsrTable bqsrTable, quantizationLevel int) ([]int, []uint8) {
	quantizationMap := make([]int, maxQualityScore+1)
	quantizedScores := make([]uint8, maxQualityScore+1)
	if quantizationLevel == 0 {
		for i := 0; i < maxQualityScore+1; i++ {
			quantizedScores[i] = uint8(i)
		}
		return quantizationMap, quantizedScores
	}
	for _, entry := range bqsrTable {
		quantizationMap[entry.EmpiricalQuality] += entry.Observations
	}
	quantizationIntervals := initializeQuantizationIntervals(quantizationMap)
	mergeQuantizationIntervals(quantizationIntervals, quantizationLevel)
	for i := 0; i >= 0; {
		intervalI := &quantizationIntervals[i]
		var quantizedScore uint8
		if leafInterval(i, intervalI) {
			quantizedScore = uint8(i)
		} else {
			errorRate := calculateErrorRate(intervalI.nobs, intervalI.nerrors)
			quantizedScore = uint8(errorProbabilityToQuality(errorRate))
		}
		j := intervalI.next
		var kend int
		if j >= 0 {
			kend = j
		} else {
			kend = len(quantizationIntervals)
		}
		for k := i; k < kend; k++ {
			quantizedScores[k] = quantizedScore
		}
		i = intervalI.next
	}
	return quantizationMap, quantizedScores
}

func estimateHierarchicalBayesianQuality(epsilon float64, empiricalReadGroupEntry *combinedBqsrEntry, empiricalBqEntry, empiricalContextEntry, empiricalCycleEntry *bqsrEntry) float64 {
	var deltaGlobalQuality, deltaReportedQuality float64
	if empiricalReadGroupEntry != nil {
		deltaGlobalQuality = float64(empiricalReadGroupEntry.calculateEmpiricalQuality(epsilon)) - epsilon
	}
	if empiricalBqEntry != nil {
		deltaReportedQuality = float64(empiricalBqEntry.calculateEmpiricalQuality(deltaGlobalQuality+epsilon)) - deltaGlobalQuality - epsilon
	}
	var deltaCovariatesQuality float64
	conditionalPrior := deltaReportedQuality + deltaGlobalQuality + epsilon
	if empiricalCycleEntry != nil {
		deltaCovariatesQuality = float64(empiricalCycleEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior
	}
	if empiricalContextEntry != nil {
		deltaCovariatesQuality += float64(empiricalContextEntry.calculateEmpiricalQuality(conditionalPrior)) - conditionalPrior
	}
	estimate := conditionalPrior + deltaCovariatesQuality
	return estimate
}

type (
	applyKey struct {
		readGroupCovariate               string
		qual                             uint8
		cycleCovariate, contextCovariate int32
	}

	applyBuffers struct {
		strandedClippedSeq        []byte
		baseContextCovariate      []int32
		recalibratedQualityScores map[applyKey]uint8
	}
)

// ApplyBQSR applies the base recalibration result to the QUAL strings of the given reads.
func (recal *BaseRecalibratorTables) ApplyBQSR(quantizeLevels int, sqqList []uint8, maxCycle int) sam.Filter {
	return func(hdr *sam.Header) sam.AlignmentFilter {
		combinedBQSRTable := initializeCombinedBQSRTable(recal.QualityScores)
		var staticQuantizedScores []uint8
		if len(sqqList) > 0 {
			staticQuantizedScores = initializeStaticQuantizedScores(sqqList)
		}
		_, quantizedQualityScores := initializeQuantizedQualityScores(recal.QualityScores, quantizeLevels)
		buffers := sync.Pool{New: func() interface{} {
			return &applyBuffers{recalibratedQualityScores: make(map[applyKey]uint8)}
		}}
		return func(aln *sam.Alignment) (result bool) {
			result = true
			readGroupCovariate := readGroupCovariate(hdr, aln)
			empiricalReadGroupEntry := combinedBQSRTable[readGroupCovariate]
			if empiricalReadGroupEntry == nil {
				return // no recalibration, bqsr table empty
			}
			buf := buffers.Get().(*applyBuffers)
			defer buffers.Put(buf)
			cycleFactor, cycleIncrement := prepareCycleCovariates(aln)
			buf.strandedClippedSeq = computeStrandedClippedSeq(aln, buf.strandedClippedSeq)
			buf.baseContextCovariate = computeBaseContextCovariate(aln, buf.strandedClippedSeq, buf.baseContextCovariate)
			var epsilon float64
			if globalQualityScorePrior > 0.0 {
				epsilon = globalQualityScorePrior
			} else {
				epsilon = empiricalReadGroupEntry.reportedQuality
			}
			// Adapt each base's quality score
			alnLength := aln.SEQ.Len()
			key := applyKey{readGroupCovariate: readGroupCovariate}
			for i := 0; i < alnLength; i++ {
				if qual := aln.QUAL[i]; qual >= minInterestingQual {
					key.qual = qual
					key.cycleCovariate = int32(computeBaseCycleCovariate(cycleFactor, cycleIncrement, i, maxCycle))
					key.contextCovariate = buf.baseContextCovariate[i]
					if recalibratedQualityScore, ok := buf.recalibratedQualityScores[key]; ok {
						aln.QUAL[i] = recalibratedQualityScore
						continue
					}
					bqKey := bqsrTableKey{
						ReadGroup: readGroupCovariate,
						Qual:      qual,
					}
					empiricalBqEntry := recal.QualityScores[bqKey]
					cycleKey := bqsrTableKey{
						Covariate: key.cycleCovariate,
						ReadGroup: readGroupCovariate,
						Qual:      qual,
					}
					empiricalCycleEntry := recal.Cycles[cycleKey]
					contextKey := bqsrTableKey{
						Covariate: key.contextCovariate,
						ReadGroup: readGroupCovariate,
						Qual:      qual,
					}
					empiricalContextEntry := recal.Contexts[contextKey]
					recalibratedQual := estimateHierarchicalBayesianQuality(epsilon, empiricalReadGroupEntry, empiricalBqEntry, empiricalContextEntry, empiricalCycleEntry)
					recalibratedQualityScore := quantizedQualityScores[maxInt(1, minInt(int(math.Round(recalibratedQual)), maxRecalibratedQualScore))]
					if len(staticQuantizedScores) != 0 {
						recalibratedQualityScore = staticQuantizedScores[recalibratedQualityScore]
					}
					buf.recalibratedQualityScores[key] = recalibratedQualityScore
					aln.QUAL[i] = recalibratedQualityScore
				}
			}
			return
		}
	}
}
