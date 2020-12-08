// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2020 imec vzw.

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
	"math"
	"sync"

	"github.com/exascience/elprep/v5/fasta"
	"github.com/exascience/elprep/v5/sam"
)

type smithWatermanOverhangStrategy int32

const (
	softclip smithWatermanOverhangStrategy = iota
	indel
	leadingIndel
	ignore
)

type int32Matrix struct {
	cols  int32
	array []int32
}

func (m *int32Matrix) ensureSize(rows, cols int32) {
	m.cols = cols
	totalSize := rows * cols
	if totalSize <= int32(cap(m.array)) {
		m.array = m.array[:totalSize]
		for i := int32(0); i < totalSize; i++ {
			m.array[i] = 0
		}
	} else {
		m.array = make([]int32, totalSize)
	}
}

func (m *int32Matrix) at(row, col int32) int32 {
	return m.array[row*m.cols+col]
}

func (m *int32Matrix) setAt(row, col, value int32) {
	m.array[row*m.cols+col] = value
}

func (m *int32Matrix) rowView(row int32) []int32 {
	offset := row * m.cols
	return m.array[offset : offset+m.cols]
}

type smithWatermanMatrices struct {
	sw, backtrack                          int32Matrix
	bestGapV, bestGapH, gapSizeV, gapSizeH []int32
}

var smithWatermanMatricesPool = sync.Pool{New: func() interface{} { return &smithWatermanMatrices{} }}

func getSmithWatermanMatrices() *smithWatermanMatrices {
	return smithWatermanMatricesPool.Get().(*smithWatermanMatrices)
}

func putSmithWatermanMatrices(sw *smithWatermanMatrices) {
	smithWatermanMatricesPool.Put(sw)
}

func ensureVector(v []int32, sz, initValue int32) (result []int32) {
	if sz <= int32(cap(v)) {
		result = v[:sz]
	} else {
		result = make([]int32, sz)
	}
	for i := int32(0); i < sz; i++ {
		result[i] = initValue
	}
	return
}

func lastIndex(ref, seq string) int32 {
	queryLength := int32(len(seq))
	for r := int32(len(ref)) - queryLength; r >= 0; r-- {
		q := int32(0)
		for q < queryLength && fasta.ToUpperAndN(ref[r+q]) == seq[q] {
			q++
		}
		if q == queryLength {
			return r
		}
	}
	return -1
}

func runSmithWaterman(reference, alternate string,
	matchValue, mismatchPenalty, gapOpenPenalty, gapExtendPenalty int32, strategy smithWatermanOverhangStrategy) ([]sam.CigarOperation, int32) {

	switch strategy {
	case softclip, ignore:
		if alignmentOffset := lastIndex(reference, alternate); alignmentOffset >= 0 {
			return []sam.CigarOperation{{int32(len(alternate)), 'M'}}, alignmentOffset
		}
	}

	sw := getSmithWatermanMatrices()
	defer putSmithWatermanMatrices(sw)

	refLength := int32(len(reference))
	altLength := int32(len(alternate))

	nrow := refLength + 1
	ncol := altLength + 1
	sw.sw.ensureSize(nrow, ncol)
	sw.backtrack.ensureSize(nrow, ncol)

	const (
		matrixMinCutoff = -1.0e8
		lowInitValue    = math.MinInt32 / 2
	)

	sw.bestGapV = ensureVector(sw.bestGapV, ncol+1, lowInitValue)
	sw.gapSizeV = ensureVector(sw.gapSizeV, ncol+1, 0)
	sw.bestGapH = ensureVector(sw.bestGapH, nrow+1, lowInitValue)
	sw.gapSizeH = ensureVector(sw.gapSizeH, nrow+1, 0)

	switch strategy {
	case indel, leadingIndel:
		topRow := sw.sw.rowView(0)
		topRow[1] = gapOpenPenalty
		currentValue := gapOpenPenalty
		for i := 2; i < len(topRow); i++ {
			currentValue += gapExtendPenalty
			topRow[i] = currentValue
		}
		sw.sw.setAt(1, 0, gapOpenPenalty)
		currentValue = gapOpenPenalty
		for i := int32(2); i < nrow; i++ {
			currentValue += gapExtendPenalty
			sw.sw.setAt(i, 0, currentValue)
		}
	}

	curRow := sw.sw.rowView(0)

	for i := int32(1); i < nrow; i++ {
		aBase := reference[i-1]
		lastRow := curRow
		curRow = sw.sw.rowView(i)
		curBacktrackRow := sw.backtrack.rowView(i)

		for j := int32(1); j < ncol; j++ {
			bBase := alternate[j-1]
			stepDiag := lastRow[j-1]
			if aBase == bBase {
				stepDiag += matchValue
			} else {
				stepDiag += mismatchPenalty
			}

			prevGap := lastRow[j] + gapOpenPenalty
			sw.bestGapV[j] += gapExtendPenalty
			if prevGap > sw.bestGapV[j] {
				sw.bestGapV[j] = prevGap
				sw.gapSizeV[j] = 1
			} else {
				sw.gapSizeV[j]++
			}

			stepDown := sw.bestGapV[j]
			kd := sw.gapSizeV[j]

			prevGap = curRow[j-1] + gapOpenPenalty
			sw.bestGapH[i] += gapExtendPenalty
			if prevGap > sw.bestGapH[i] {
				sw.bestGapH[i] = prevGap
				sw.gapSizeH[i] = 1
			} else {
				sw.gapSizeH[i]++
			}

			stepRight := sw.bestGapH[i]
			ki := sw.gapSizeH[i]

			if stepDiag >= stepDown && stepDiag >= stepRight {
				curRow[j] = maxInt32(matrixMinCutoff, stepDiag)
				curBacktrackRow[j] = 0
			} else if stepRight >= stepDown {
				curRow[j] = maxInt32(matrixMinCutoff, stepRight)
				curBacktrackRow[j] = -ki
			} else {
				curRow[j] = maxInt32(matrixMinCutoff, stepDown)
				curBacktrackRow[j] = kd
			}
		}
	}

	maxScore := math.MinInt32
	var segmentLength int32
	var p1 int32
	p2 := altLength

	if strategy == indel {
		p1 = refLength
	} else {
		for i := int32(1); i < nrow; i++ {
			curScore := int(sw.sw.at(i, altLength))
			if curScore >= maxScore {
				p1 = i
				maxScore = curScore
			}
		}
		if strategy != leadingIndel {
			bottomRow := sw.sw.rowView(refLength)
			for j := int32(1); j < ncol; j++ {
				if curScore := int(bottomRow[j]); curScore > maxScore || (curScore == maxScore && absInt32(refLength-j) < absInt32(p1-p2)) {
					p1 = refLength
					p2 = j
					maxScore = curScore
					segmentLength = altLength - j
				}
			}
		}
	}

	lce := make([]sam.CigarOperation, 0, 5)
	if segmentLength > 0 && strategy == softclip {
		lce = append(lce, sam.CigarOperation{int32(segmentLength), 'S'})
		segmentLength = 0
	}
	state := byte('M')
	for {
		stepLength := int32(1)
		btr := sw.backtrack.at(p1, p2)
		var newState byte
		if btr > 0 {
			newState = 'D'
			stepLength = btr
			p1 -= btr
		} else if btr < 0 {
			newState = 'I'
			stepLength = -btr
			p2 += btr
		} else {
			newState = 'M'
			p1--
			p2--
		}

		if newState == state {
			segmentLength += stepLength
		} else {
			lce = append(lce, sam.CigarOperation{int32(segmentLength), state})
			segmentLength = stepLength
			state = newState
		}

		if p1 <= 0 || p2 <= 0 {
			break
		}
	}

	var alignmentOffset int32
	switch strategy {
	case softclip:
		lce = append(lce, sam.CigarOperation{int32(segmentLength), state})
		if p2 > 0 {
			lce = append(lce, sam.CigarOperation{p2, 'S'})
		}
		alignmentOffset = p1
	case ignore:
		lce = append(lce, sam.CigarOperation{int32(segmentLength) + p2, state})
		alignmentOffset = p1 - p2
	default:
		lce = append(lce, sam.CigarOperation{int32(segmentLength), state})
		switch {
		case p1 > 0:
			lce = append(lce, sam.CigarOperation{p1, 'D'})
		case p2 > 0:
			lce = append(lce, sam.CigarOperation{p2, 'I'})
		}
		alignmentOffset = 0
	}

	for i, j := 0, len(lce)-1; i < j; i, j = i+1, j-1 {
		lce[i], lce[j] = lce[j], lce[i]
	}
	for i := 1; i < len(lce); {
		if lce[i-1].Length == 0 {
			lce = append(lce[:i-1], lce[i:]...)
		} else if lce[i-1].Operation == lce[i].Operation {
			lce[i-1].Length += lce[i].Length
			lce = append(lce[:i], lce[i+1:]...)
		} else {
			i++
		}
	}
	if l := len(lce) - 1; lce[l].Length == 0 {
		lce = lce[:l]
	}
	return lce, alignmentOffset
}

const swPad = "NNNNNNNNNN"

func isSWFailure(cigar []sam.CigarOperation, alignmentOffset int32) bool {
	if alignmentOffset > 0 {
		return true
	}
	for _, op := range cigar {
		if op.Operation == 'S' {
			return true
		}
	}
	return false
}

func trimCigarByBases(cigar []sam.CigarOperation, start, end int32) (newCigar []sam.CigarOperation) {
	pos := int32(0)
	for _, elt := range cigar {
		if elt.Operation == 'D' {
			if pos >= start {
				newCigar = append(newCigar, elt)
				continue
			}
		} else if pos > end {
			break
		}
		newCigar, pos = addCigarElement(newCigar, pos, start, end, elt)
	}
	for i := 1; i < len(newCigar); i++ {
		if newCigar[i-1].Operation == newCigar[i].Operation {
			newCigar[i-1].Length += newCigar[i].Length
			newCigar = append(newCigar[:i], newCigar[i+1:]...)
		} else {
			i++
		}
	}
	return newCigar
}

func leftAlignCigarSequentially(cigar []sam.CigarOperation, reference, alternate string) (newCigar []sam.CigarOperation) {
	var cigarToAlign []sam.CigarOperation
	var refIndex, readIndex int32
	for _, ce := range cigar {
		switch ce.Operation {
		case 'D', 'I':
			cigarToAlign = append(cigarToAlign, ce)
			newCigar = append(newCigar, leftAlignIndel(cigarToAlign, reference, alternate, refIndex, readIndex, false)...)
			refIndex += sam.ReferenceLengthFromCigar(cigarToAlign)
			readIndex += sam.ReadLengthFromCigar(cigarToAlign)
			cigarToAlign = cigarToAlign[:0]
		default:
			cigarToAlign = append(cigarToAlign, ce)
		}
	}
	newCigar = append(newCigar, cigarToAlign...)
	for len(newCigar) > 0 && newCigar[0].Length == 0 {
		newCigar = newCigar[1:]
	}
	for i := 1; i < len(newCigar); {
		if newCigar[i].Length == 0 {
			newCigar = append(newCigar[:i], newCigar[i+1:]...)
		} else if newCigar[i-1].Operation == newCigar[i].Operation {
			newCigar[i-1].Length += newCigar[i].Length
			newCigar = append(newCigar[:i], newCigar[i+1:]...)
		} else {
			i++
		}
	}
	return newCigar
}

func calculateCigar(reference, alternate, paddedRef string, strategy smithWatermanOverhangStrategy) []sam.CigarOperation {
	if len(reference) == len(alternate) {
		mismatches := 0
		for i := range reference {
			if reference[i] != alternate[i] {
				mismatches++
			}
		}
		if mismatches <= 2 {
			return []sam.CigarOperation{{int32(len(reference)), 'M'}}
		}
	}
	paddedAlt := swPad + alternate + swPad
	cigar, alignmentOffset := runSmithWaterman(paddedRef, paddedAlt, 200, -150, -260, -11, strategy)
	if isSWFailure(cigar, alignmentOffset) {
		return nil
	}
	baseStart := int32(len(swPad))
	baseEnd := int32(len(paddedAlt) - len(swPad) - 1)
	nonStandard := trimCigarByBases(cigar, baseStart, baseEnd)
	if refLength := int(sam.ReferenceLengthFromCigar(nonStandard)); refLength != len(reference) {
		nonStandard = append(nonStandard, sam.CigarOperation{int32(len(reference) - refLength), 'D'})
	}
	return leftAlignCigarSequentially(nonStandard, reference, alternate)
}
