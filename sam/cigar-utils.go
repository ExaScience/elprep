package sam

import (
	"log"
	"strconv"

	"github.com/exascience/elprep/internal"
)

var cigarConsumesReferenceBases = map[byte]int32{'M': 1, 'D': 1, 'N': 1, '=': 1, 'X': 1}

// Sums the lengths of all CIGAR operations that consume reference
// bases.
func end(aln *Alignment, cigars []CigarOperation) int32 {
	var length int32
	for _, op := range cigars {
		length += cigarConsumesReferenceBases[op.Operation] * op.Length
	}
	return aln.POS + length - 1
}

func operatorConsumesReadBases(operator byte) bool {
	switch operator {
	case 'M', 'I', 'S', '=', 'X':
		return true
	default:
		return false
	}
}

func operatorConsumesReferenceBases(operator byte) bool {
	switch operator {
	case 'M', 'D', 'N', '=', 'X':
		return true
	default:
		return false
	}
}

// Sums the lengths of all CIGAR operations that consume read bases.
func readLengthFromCigar(cigars []CigarOperation) int32 {
	var length int32
	for _, op := range cigars {
		if operatorConsumesReadBases(op.Operation) {
			length += op.Length
		}
	}
	return length
}

func elementStradlessClippedRead(newCigar []byte, operator byte, relativeClippingPosition, clippedBases int32) []byte {
	if operatorConsumesReadBases(operator) {
		if operatorConsumesReferenceBases(operator) {
			if relativeClippingPosition > 0 {
				newCigar = strconv.AppendInt(append(newCigar, operator), int64(relativeClippingPosition), 10)
			}
		} else {
			clippedBases += relativeClippingPosition
		}
	} else if relativeClippingPosition != 0 {
		log.Fatal("Unexpected non-0 relative clipping position in CleanSam.")
	}
	return strconv.AppendInt(append(newCigar, 'S'), int64(clippedBases), 10)
}

func softClipEndOfRead(clipFrom int32, cigars []CigarOperation) string {
	var pos int32
	clipFrom--
	newCigar := internal.ReserveByteBuffer()
	defer internal.ReleaseByteBuffer(newCigar)
	for _, op := range cigars {
		endPos := pos
		if operatorConsumesReadBases(op.Operation) {
			endPos += op.Length
		}
		if endPos < clipFrom {
			*newCigar = strconv.AppendInt(append(*newCigar, op.Operation), int64(op.Length), 10)
		} else {
			clippedBases := readLengthFromCigar(cigars) + clipFrom
			*newCigar = elementStradlessClippedRead(*newCigar, op.Operation, clipFrom-pos, clippedBases)
			break
		}
		pos += endPos
	}
	return string(*newCigar)
}
