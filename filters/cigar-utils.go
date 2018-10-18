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

	"github.com/exascience/elprep/v4/sam"
)

var (
	cigarConsumesReadBases      = map[byte]int32{'M': 1, 'I': 1, 'S': 1, '=': 1, 'X': 1}
	cigarConsumesReferenceBases = map[byte]int32{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}
)

// Sums the lengths of all CIGAR operations that consume reference
// bases.
func end(aln *sam.Alignment, cigars []sam.CigarOperation) int32 {
	var length int32
	for _, op := range cigars {
		length += cigarConsumesReferenceBases[op.Operation] * op.Length
	}
	return aln.POS + length - 1
}

var (
	operatorConsumesReadBases      = map[byte]bool{'M': true, 'I': true, 'S': true, '=': true, 'X': true}
	operatorConsumesReferenceBases = map[byte]bool{'M': true, 'D': true, 'N': true, '=': true, 'X': true}
)

// Sums the lengths of all CIGAR operations that consume read bases.
func readLengthFromCigar(cigars []sam.CigarOperation) int32 {
	var length int32
	for _, op := range cigars {
		length += cigarConsumesReadBases[op.Operation] * op.Length
	}
	return length
}

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
		log.Fatal("Unexpected non-0 relative clipping position in CleanSam.")
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
		endPos += cigarConsumesReadBases[op.Operation] * op.Length
		if endPos < clipFrom {
			newCigar = append(newCigar, op)
		} else {
			clippedBases := readLengthFromCigar(cigars) + clipFrom
			newCigar = elementStradlessClippedRead(newCigar, op.Operation, clipFrom-pos, clippedBases)
			break
		}
		pos += endPos
	}
	return newCigar
}
