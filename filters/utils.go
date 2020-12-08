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

	"github.com/exascience/elprep/v5/sam"
)

func minUint8(x, y uint8) uint8 {
	if x < y {
		return x
	}
	return y
}

func minInt32(x, y int32) int32 {
	if x < y {
		return x
	}
	return y
}

func maxInt32(x, y int32) int32 {
	if x > y {
		return x
	}
	return y
}

func minInt(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func maxInt(x, y int) int {
	if x > y {
		return x
	}
	return y
}

func absInt(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

func absInt32(x int32) int32 {
	if x < 0 {
		return -x
	}
	return x
}

var (
	operatorConsumesReadBases      = map[byte]bool{'M': true, 'I': true, 'S': true, '=': true, 'X': true}
	operatorConsumesReferenceBases = map[byte]bool{'M': true, 'D': true, 'N': true, '=': true, 'X': true}
)

func elementStradlessClippedRead(newCigar []sam.CigarOperation, operator byte, relativeClippingPosition, clippedBases int32) []sam.CigarOperation {
	if operatorConsumesReadBases[operator] {
		if operatorConsumesReferenceBases[operator] {
			if relativeClippingPosition > 0 {
				newCigar = append(newCigar, sam.CigarOperation{
					Length:    relativeClippingPosition,
					Operation: operator,
				})
			}
		} else {
			clippedBases += relativeClippingPosition
		}
	} else if relativeClippingPosition != 0 {
		log.Panic("Unexpected non-0 relative clipping position in CleanSam.")
	}
	return append(newCigar, sam.CigarOperation{
		Length:    clippedBases,
		Operation: 'S',
	})
}

func softClipEndOfRead(clipFrom int32, cigars []sam.CigarOperation) []sam.CigarOperation {
	var pos int32
	clipFrom--
	var newCigar []sam.CigarOperation
	for _, op := range cigars {
		endPos := pos
		endPos += sam.CigarOperatorConsumesReadBases[op.Operation] * op.Length
		if endPos < clipFrom {
			newCigar = append(newCigar, op)
		} else {
			clippedBases := sam.ReadLengthFromCigar(cigars) + clipFrom
			newCigar = elementStradlessClippedRead(newCigar, op.Operation, clipFrom-pos, clippedBases)
			break
		}
		pos += endPos
	}
	return newCigar
}

func cigarContainsN(cigarVec []sam.CigarOperation) bool {
	for _, cigarOp := range cigarVec {
		if cigarOp.Operation == 'N' {
			return true
		}
	}
	return false
}

func alignmentAgreesWithHeader(hdr *sam.Header, aln *sam.Alignment) bool {
	rname := aln.RNAME
	for _, sq := range hdr.SQ {
		if sq["SN"] == rname {
			return aln.POS <= sam.SQLN(sq)
		}
	}
	return false
}

func isStrictUnmapped(aln *sam.Alignment) bool {
	return aln.IsUnmapped() || aln.RNAME == "" || aln.RNAME == "*" || aln.POS == 0
}

func isStrictNextUnmapped(aln *sam.Alignment) bool {
	return aln.IsNextUnmapped() || aln.RNEXT == "" || aln.RNEXT == "*" || aln.PNEXT == 0
}

func hasWellDefinedFragmentSize(aln *sam.Alignment) (bool, int) {
	if aln.TLEN != 0 && aln.IsMultiple() && !isStrictUnmapped(aln) && !isStrictNextUnmapped(aln) && aln.IsReversed() != aln.IsNextReversed() {
		if aln.IsReversed() {
			alnEnd := aln.End()
			return alnEnd > aln.PNEXT, int(alnEnd)
		}
		return aln.POS <= aln.PNEXT+aln.TLEN, -1
	}
	return false, -1
}

func computeAdaptorBoundary(aln *sam.Alignment) (int, int, bool) {
	if wellDefinedSize, alnEnd := hasWellDefinedFragmentSize(aln); wellDefinedSize {
		var boundary int
		if aln.IsReversed() {
			boundary = int(aln.PNEXT) - 1
		} else {
			boundary = int(aln.POS) + absInt(int(aln.TLEN))
		}
		return boundary, alnEnd, true
	}
	return -1, -1, false
}

func isInsideRead(aln *sam.Alignment, alnEnd, refCoord int) bool {
	if refCoord >= int(aln.POS) {
		if alnEnd < 0 {
			alnEnd = int(aln.End())
		}
		return refCoord <= alnEnd
	}
	return false
}

func isInsideDeletion(cigar []sam.CigarOperation, offset int) bool {
	if offset < 0 {
		return false
	}

	var pos, prevPos int

	for _, ce := range cigar {
		switch ce.Operation {
		case 'I', 'S', 'D', 'M', '=', 'X':
			prevPos = pos
			pos += int(ce.Length)
		}

		if prevPos < offset && pos >= offset && ce.Operation == 'D' {
			return true
		}
	}

	return false
}

func nofAlignedBasesWithSoftClips(cigar []sam.CigarOperation) (n int) {
	for _, ce := range cigar {
		switch ce.Operation {
		case 'M', '=', 'X', 'S':
			n += int(ce.Length)
		}
	}
	return n
}

func hardClipAdaptorSequence(aln *sam.Alignment) {
	if adaptorBoundary, alnEnd, ok := computeAdaptorBoundary(aln); ok && isInsideRead(aln, alnEnd, adaptorBoundary) {
		if aln.IsReversed() {
			hardClipByReferenceCoordinatesLeftTail(aln, adaptorBoundary)
		} else {
			hardClipByReferenceCoordinatesRightTail(aln, adaptorBoundary)
		}
	}
}

func softStart(aln *sam.Alignment) int {
	start := aln.POS
	for _, op := range aln.CIGAR {
		if op.Operation == 'S' {
			start -= op.Length
		} else if op.Operation != 'H' {
			break
		}
	}
	return int(start)
}

func softEnd(aln *sam.Alignment) int {
	end := aln.End()
	softEnd := end
	for i := len(aln.CIGAR) - 1; i >= 0; i-- {
		op := aln.CIGAR[i]
		if op.Operation == 'S' {
			softEnd += op.Length
		} else if op.Operation != 'H' {
			return int(softEnd)
		}
	}
	return int(end)
}

func hardClipByReferenceCoordinatesLeftTail(aln *sam.Alignment, refStop int) {
	stop, ok := getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart(aln), refStop, left)
	if !ok {
		log.Panicf("reference coordinate matches a non-existing base in read %v", aln.QNAME)
	}
	hardClip(aln, 0, stop)
}

func hardClipByReferenceCoordinatesRightTail(aln *sam.Alignment, refStart int) {
	start, ok := getReadCoordinateForReferenceCoordinate(aln.CIGAR, softStart(aln), refStart, right)
	stop := aln.SEQ.Len() - 1
	if !ok {
		log.Panicf("reference coordinate matches a non-existing base in read %v", aln.QNAME)
	}
	hardClip(aln, start, stop)
}

func computeReadCoordinateForReferenceCoordinate(cigarVec []sam.CigarOperation, softStart, refIndex int) (int, bool) {
	goal := refIndex - softStart
	if goal < 0 {
		return -1, false
	}
	var readBases, refBases int
	fallsInsideDeletionOrSkippedRegion := false
	endsJustBeforeDeletionOrSkippedRegion := false
	fallsInsideOrJustBeforeDeletionOrSkippedRegion := false
	index := 0
	for refBases != goal && index < len(cigarVec) {
		element := cigarVec[index]
		index++
		elementLength := int(element.Length)
		var shift int
		if operatorConsumesReferenceBases[element.Operation] || element.Operation == 'S' {
			if refBases+elementLength < goal {
				shift = elementLength
			} else {
				shift = goal - refBases
			}
			refBases += shift
		}
		if refBases != goal {
			readBases += int(sam.CigarOperatorConsumesReadBases[element.Operation]) * elementLength
		} else {
			if shift >= elementLength && index == len(cigarVec) {
				return -1, false
			}
			var nextCigar sam.CigarOperation
			if shift < elementLength {
				fallsInsideDeletionOrSkippedRegion = element.Operation == 'D' || element.Operation == 'N'
			} else {
				nextCigar = cigarVec[index]
				index++
				if nextCigar.Operation == 'I' {
					readBases += int(nextCigar.Length)
					if index == len(cigarVec) {
						return -1, false
					}
					nextCigar = cigarVec[index]
					index++
				}
				endsJustBeforeDeletionOrSkippedRegion = nextCigar.Operation == 'D' || nextCigar.Operation == 'N'
			}
			fallsInsideOrJustBeforeDeletionOrSkippedRegion = endsJustBeforeDeletionOrSkippedRegion || fallsInsideDeletionOrSkippedRegion
			if !fallsInsideOrJustBeforeDeletionOrSkippedRegion {
				readBases += int(sam.CigarOperatorConsumesReadBases[element.Operation]) * shift
			} else if endsJustBeforeDeletionOrSkippedRegion {
				readBases += int(sam.CigarOperatorConsumesReadBases[element.Operation]) * (shift - 1)
			} else if fallsInsideDeletionOrSkippedRegion || (endsJustBeforeDeletionOrSkippedRegion && (nextCigar.Operation == 'D' || nextCigar.Operation == 'N')) {
				readBases--
			}
		}
	}
	if refBases != goal {
		return -1, false
	}
	return readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion
}

type clippingTail int

const (
	left clippingTail = iota
	right
)

func getReadCoordinateForReferenceCoordinate(cigarVec []sam.CigarOperation, softStart, refIndex int, tail clippingTail) (int, bool) {
	readBases, fallsInsideOrJustBeforeDeletionOrSkippedRegion := computeReadCoordinateForReferenceCoordinate(cigarVec, softStart, refIndex)
	if readBases == -1 {
		return -1, false
	}
	if tail == right && fallsInsideOrJustBeforeDeletionOrSkippedRegion {
		readBases++
	}
	if tail == left && readBases == 0 {
		if firstInsertionLength, firstInsertion := readStartsWithInsertion(cigarVec); firstInsertion {
			readBases = int(minInt32(firstInsertionLength, sam.ReadLengthFromCigar(cigarVec)-1))
		}
	}
	return readBases, true
}

func calculateHardSoftOffset(cigar []sam.CigarOperation) int32 {
	var size int32
	var i int
	for ; i < len(cigar); i++ {
		cigarOp := cigar[i]
		if cigarOp.Operation == 'H' {
			size += cigarOp.Length
		} else {
			break
		}
	}
	for ; i < len(cigar); i++ {
		cigarOp := cigar[i]
		if cigarOp.Operation == 'S' {
			size += cigarOp.Length
		} else {
			break
		}
	}
	return size
}

func calculateAlnStartShift(cigar, clippedCigar []sam.CigarOperation) int32 {
	return calculateHardSoftOffset(clippedCigar) - calculateHardSoftOffset(cigar)
}

func calculateHardClippingAlignmentShift(cigarOp sam.CigarOperation, cigarLength int) int {
	switch cigarOp.Operation {
	case 'I':
		return -cigarLength
	case 'D', 'N':
		return int(cigarOp.Length)
	default:
		return 0
	}
}

func hardClip(aln *sam.Alignment, start, stop int) {
	clippedCigar := hardClipCigar(aln, start, stop)
	readLength := aln.SEQ.Len()
	newLength := readLength - (stop - start + 1)
	copyStart := 0
	if start == 0 {
		copyStart = stop + 1
	}
	newBases := aln.SEQ.Slice(copyStart, copyStart+newLength)
	newQuals := aln.QUAL[copyStart : copyStart+newLength]
	cigarVec := aln.CIGAR
	aln.SEQ = newBases
	aln.QUAL = newQuals
	aln.CIGAR = clippedCigar
	if start == 0 && !isStrictUnmapped(aln) {
		aln.POS += calculateAlnStartShift(cigarVec, clippedCigar)
	}
}

func hardClipCigar(aln *sam.Alignment, start, stop int) []sam.CigarOperation {
	cigarVec := aln.CIGAR
	index := 0
	totalHardClipCount := stop - start + 1
	alignmentShift := 0
	var newCigar []sam.CigarOperation
	if start == 0 {
		var cigarOpIndex int
		var cigarOp sam.CigarOperation
		for cigarOpIndex, cigarOp = range cigarVec {
			if cigarOp.Operation != 'H' {
				break
			}
			totalHardClipCount += int(cigarOp.Length)
		}
		for ; index <= stop && cigarOpIndex < len(cigarVec); cigarOpIndex++ {
			cigarOp := cigarVec[cigarOpIndex]
			cigarOpLength := int(cigarOp.Length)
			shift := int(sam.CigarOperatorConsumesReadBases[cigarOp.Operation]) * cigarOpLength
			if index+shift == stop+1 {
				alignmentShift += calculateHardClippingAlignmentShift(cigarOp, cigarOpLength)
				newCigar = append(newCigar, sam.CigarOperation{Operation: 'H', Length: int32(totalHardClipCount + alignmentShift)})
			} else if index+shift > stop+1 {
				lengthAfterChopping := cigarOpLength - (stop - index + 1)
				alignmentShift += calculateHardClippingAlignmentShift(cigarOp, stop-index+1)
				newCigar = append(
					newCigar,
					sam.CigarOperation{Operation: 'H', Length: int32(totalHardClipCount + alignmentShift)},
					sam.CigarOperation{Operation: cigarOp.Operation, Length: int32(lengthAfterChopping)},
				)
			}
			index += shift
			alignmentShift += calculateHardClippingAlignmentShift(cigarOp, shift)
		}
		newCigar = append(newCigar, cigarVec[cigarOpIndex:]...)
	} else {
		var cigarOpIndex int
		for ; index < start && cigarOpIndex < len(cigarVec); cigarOpIndex++ {
			cigarOp := cigarVec[cigarOpIndex]
			cigarOpLength := int(cigarOp.Length)
			shift := int(sam.CigarOperatorConsumesReadBases[cigarOp.Operation]) * cigarOpLength
			if index+shift < start {
				newCigar = append(newCigar, cigarOp)
			} else {
				lengthAfterChopping := start - index
				alignmentShift += calculateHardClippingAlignmentShift(cigarOp, cigarOpLength-(start-index))
				if cigarOp.Operation == 'H' {
					totalHardClipCount += lengthAfterChopping
				} else {
					newCigar = append(newCigar, sam.CigarOperation{Operation: cigarOp.Operation, Length: int32(lengthAfterChopping)})
				}
			}
			index += shift
		}
		for ; cigarOpIndex < len(cigarVec); cigarOpIndex++ {
			cigarOp := cigarVec[cigarOpIndex]
			alignmentShift += calculateHardClippingAlignmentShift(cigarOp, int(cigarOp.Length))
			if cigarOp.Operation == 'H' {
				totalHardClipCount += int(cigarOp.Length)
			}
		}
		newCigar = append(newCigar, sam.CigarOperation{Operation: 'H', Length: int32(totalHardClipCount + alignmentShift)})
	}
	return cleanHardClippedCigar(newCigar)
}

func cleanHardClippedCigar(cigar []sam.CigarOperation) []sam.CigarOperation {
	totalHardClip := 0
	index := 0
forward:
	for ; index < len(cigar); index++ {
		switch element := cigar[index]; element.Operation {
		case 'H', 'D', 'N':
			totalHardClip += int(element.Length)
		default:
			break forward
		}
	}
	if index > 0 {
		cigar[0] = sam.CigarOperation{Operation: 'H', Length: int32(totalHardClip)}
		cigar = append(cigar[:1], cigar[index:]...)
	}
	totalHardClip = 0
	index = len(cigar) - 1
backward:
	for ; index >= 0; index-- {
		switch element := cigar[index]; element.Operation {
		case 'H', 'D', 'N':
			totalHardClip += int(element.Length)
		default:
			break backward
		}
	}
	if index < len(cigar)-1 {
		cigar = append(cigar[:index+1], sam.CigarOperation{Operation: 'H', Length: int32(totalHardClip)})
	}
	return cigar
}

func hardClipSoftClippedBases(aln *sam.Alignment) {
	cigarVec := aln.CIGAR
	readIndex := 0
	cutLeft := -1
	cutRight := -1
	rightTail := false
	for _, cigarOp := range cigarVec {
		key := cigarOp.Operation
		ln := int(cigarOp.Length)
		switch key {
		case 'S':
			if rightTail {
				cutRight = readIndex
			} else {
				cutLeft = readIndex + ln - 1
			}
		case 'H':
		default:
			rightTail = true
		}
		readIndex += int(sam.CigarOperatorConsumesReadBases[key]) * ln
	}
	if cutRight >= 0 {
		hardClip(aln, cutRight, aln.SEQ.Len()-1)
	}
	if cutLeft >= 0 {
		hardClip(aln, 0, cutLeft)
	}
}

func emptyRead(aln *sam.Alignment) {
	aln.FLAG |= sam.Unmapped
	aln.MAPQ = 0
	aln.CIGAR = nil
	aln.SEQ = sam.Sequence{}
	aln.QUAL = nil
	rg := aln.RG()
	aln.TAGS = nil
	if rg != nil {
		aln.SetRG(rg)
	}
}

func hardClipLowQualEnds(aln *sam.Alignment, lowQual byte) {
	length := aln.SEQ.Len()
	left, right := 0, length-1
	for right >= 0 && aln.QUAL[right] <= lowQual {
		right--
	}
	for left < length && aln.QUAL[left] <= lowQual {
		left++
	}
	if left > right {
		emptyRead(aln)
		return
	}
	if right < length-1 {
		hardClip(aln, right+1, length-1)
	}
	if left > 0 {
		hardClip(aln, 0, left-1)
	}
}

func revertSoftClippedBases(aln *sam.Alignment) {
	var unclipped []sam.CigarOperation
	var matches int32
	for _, op := range aln.CIGAR {
		switch op.Operation {
		case 'S', 'M':
			matches += op.Length
		default:
			if matches > 0 {
				unclipped = append(unclipped, sam.CigarOperation{
					Operation: 'M',
					Length:    matches,
				})
				matches = 0
			}
			unclipped = append(unclipped, op)
		}
	}
	if matches > 0 {
		unclipped = append(unclipped, sam.CigarOperation{
			Operation: 'M',
			Length:    matches,
		})
	}
	newStart := aln.POS + calculateAlnStartShift(aln.CIGAR, unclipped)
	aln.CIGAR = unclipped
	if newStart <= 0 {
		aln.POS = 1
		hardClip(aln, 0, -int(newStart))
		if !isStrictUnmapped(aln) {
			aln.POS = 1
		}
	} else {
		aln.POS = newStart
	}
}

func hardClipToRegion(aln *sam.Alignment, start, stop int) {
	if aln.SEQ.Len() == 0 || start-1 == stop+1 {
		emptyRead(aln)
		return
	}
	alnStart, alnStop := int(aln.POS), int(aln.End())
	if alnStart <= stop && alnStop >= start {
		if alnStop > stop {
			hardClipByReferenceCoordinatesRightTail(aln, stop+1)
			if alnStart < start && start-1 > int(aln.End()) {
				emptyRead(aln)
				return
			}
		}
		if alnStart < start {
			hardClipByReferenceCoordinatesLeftTail(aln, start-1)
		}
	} else {
		emptyRead(aln)
	}
}
