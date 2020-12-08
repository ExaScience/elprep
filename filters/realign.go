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
	"bytes"
	"math"

	"github.com/exascience/pargo/parallel"

	"github.com/exascience/elprep/v5/sam"
)

func (h *haplotype) priority() (result int) {
	if h.isRef {
		result = 1
	}
	if len(h.cigar) > 0 {
		result += 1 - len(h.cigar)
	}
	return result
}

type cigarTransform struct {
	op13                 byte
	advance12, advance23 int32
}

var cigarTransformers map[byte]map[byte]cigarTransform

func init() {
	cigarSets := map[byte][]byte{
		'M': []byte{'M', '=', 'X'},
		'I': []byte{'I', 'S'},
		'D': []byte{'D'},
	}

	cigarTransformers = make(map[byte]map[byte]cigarTransform)
	for _, op := range []byte{'M', '=', 'X', 'I', 'S', 'D'} {
		cigarTransformers[op] = make(map[byte]cigarTransform)
	}

	def := func(op12, op23, op13 byte, advance12, advance23 int32) {
		t := cigarTransform{op13, advance12, advance23}
		for _, op12 := range cigarSets[op12] {
			for _, op23 := range cigarSets[op23] {
				cigarTransformers[op12][op23] = t
			}
		}
	}

	def('M', 'M', 'M', 1, 1)
	def('M', 'I', 'I', 1, 1)
	def('M', 'D', 'D', 0, 1)

	def('D', 'M', 'D', 1, 1)
	def('D', 'D', 'D', 0, 1)
	def('D', 'I', 0, 1, 1)

	def('I', 'M', 'I', 1, 0)
	def('I', 'D', 'I', 1, 0)
	def('I', 'I', 'I', 1, 0)
}

func applyCigarToCigar(firstToSecond, secondToThird []sam.CigarOperation) []sam.CigarOperation {
	var result []sam.CigarOperation
	index := -1
	var addOp func(byte)
	addOp = func(op byte) {
		if op == 0 || op == 'D' {
			return
		}
		result = append(result, sam.CigarOperation{1, op})
		index = 0
		addOp = func(op byte) {
			if op == 0 {
				return
			}
			if result[index].Operation == op {
				result[index].Length++
			} else {
				result = append(result, sam.CigarOperation{1, op})
				index++
			}
		}
	}

	nElements12 := int32(len(firstToSecond))
	nElements23 := int32(len(secondToThird))

	var cigar12I, cigar23I, elt12I, elt23I int32

	for cigar12I < nElements12 && cigar23I < nElements23 {
		elt12 := firstToSecond[cigar12I]
		elt23 := secondToThird[cigar23I]

		transform := cigarTransformers[elt12.Operation][elt23.Operation]

		addOp(transform.op13)

		elt12I += transform.advance12
		elt23I += transform.advance23

		if elt12I == elt12.Length {
			cigar12I++
			elt12I = 0
		}
		if elt23I == elt23.Length {
			cigar23I++
			elt23I = 0
		}
	}

	return result
}

func createIndelString(alt []byte, cigar []sam.CigarOperation, indelIndex int, indel sam.CigarOperation, refSeq, readSeq string, refIndex, readIndex int32) []byte {
	var totalRefBases int32
	for _, ce := range cigar[:indelIndex] {
		switch ce.Operation {
		case 'M', '=', 'X':
			readIndex += ce.Length
			refIndex += ce.Length
			totalRefBases += ce.Length
		case 'S':
			readIndex += ce.Length
		case 'N':
			refIndex += ce.Length
			totalRefBases += ce.Length
		}
	}

	if refIndex > int32(len(refSeq)) {
		return nil
	}

	var indelLength int32

	if totalRefBases+indel.Length > int32(len(refSeq)) {
		indelLength = int32(len(refSeq)) - totalRefBases
	} else {
		indelLength = indel.Length
	}

	altLength := int32(len(refSeq))
	if indel.Operation == 'D' {
		altLength -= indelLength
	} else {
		altLength += indelLength
	}
	if refIndex > altLength {
		return nil
	}
	if altLength <= int32(cap(alt)) {
		alt = alt[:altLength]
	} else {
		alt = make([]byte, altLength)
	}
	copy(alt, refSeq[:refIndex])

	currentPos := refIndex

	if indel.Operation == 'D' {
		refIndex += indelLength
	} else {
		copy(alt[currentPos:currentPos+indelLength], readSeq[readIndex:readIndex+indelLength])
		currentPos += indelLength
	}

	if int32(len(refSeq))-refIndex > altLength-currentPos {
		return nil
	}
	copy(alt[currentPos:], refSeq[refIndex:])

	return alt
}

func moveCigarLeft(cigar []sam.CigarOperation, indelIndex int) []sam.CigarOperation {
	elements := make([]sam.CigarOperation, indelIndex-1, len(cigar))
	copy(elements, cigar)
	ce := cigar[indelIndex-1]
	elements = append(elements, sam.CigarOperation{maxInt32(ce.Length-1, 0), ce.Operation})
	elements = append(elements, cigar[indelIndex])
	if indelIndex+1 < len(cigar) {
		ce = cigar[indelIndex+1]
		elements = append(elements, sam.CigarOperation{ce.Length + 1, ce.Operation})
		elements = append(elements, cigar[indelIndex+2:]...)
	} else {
		elements = append(elements, sam.CigarOperation{1, 'M'})
	}
	return elements
}

func leftAlignIndel(cigar []sam.CigarOperation, refSeq, readSeq string, refIndex, readIndex int32, cleanupCigar bool) []sam.CigarOperation {
	indelIndex := -1
	var indel sam.CigarOperation
	for index, ce := range cigar {
		if ce.Operation == 'D' || ce.Operation == 'I' {
			if indelIndex != -1 {
				return cigar
			}
			indelIndex = index
			indel = ce
		}
	}
	if indelIndex <= 0 {
		return cigar
	}

	altString := createIndelString(nil, cigar, indelIndex, indel, refSeq, readSeq, refIndex, readIndex)
	if len(altString) == 0 {
		return cigar
	}

	newCigar := cigar
	var newAltString []byte
	for i := int32(0); i < indel.Length; i++ {
		newCigar = moveCigarLeft(newCigar, indelIndex)
		newAltString = createIndelString(newAltString, newCigar, indelIndex, indel, refSeq, readSeq, refIndex, readIndex)
		if bytes.Equal(altString, newAltString) {
			cigar = newCigar
			i = -1
			for _, ce := range newCigar {
				if ce.Length == 0 {
					if cleanupCigar {
						for i, ce := range cigar {
							if ce.Length != 0 && ce.Operation != 'D' {
								cigar = append(cigar[:0], cigar[i:]...)
								break
							}
						}
						for i := 1; i < len(cigar); {
							if cigar[i].Length == 0 {
								cigar = append(cigar[:i], cigar[i+1:]...)
							} else {
								i++
							}
						}
					}
					return cigar
				}
			}
		} else {
			for _, ce := range newCigar {
				if ce.Length == 0 {
					return cigar
				}
			}
		}
	}
	return cigar
}

func realignReadsToTheirBestHaplotype(likelihoods readLikelihoods, haplotypes []*haplotype) {
	var refHaplotype *haplotype
	for _, h := range haplotypes {
		if h.isRef {
			refHaplotype = h
			break
		}
	}
	parallel.Range(0, len(likelihoods.alns), 0, func(low, high int) {
		for r := low; r < high; r++ {
			aln := likelihoods.alns[r]
			var bestAllele, secondBestAllele *haplotype
			bestLikelihood, secondBestLikelihood := math.Inf(-1), math.Inf(-1)
			for _, h := range haplotypes {
				if likelihood := likelihoods.values[h][r]; likelihood > bestLikelihood {
					secondBestAllele, secondBestLikelihood = bestAllele, bestLikelihood
					bestAllele, bestLikelihood = h, likelihood
				} else if likelihood > secondBestLikelihood {
					secondBestAllele, secondBestLikelihood = h, likelihood
				}
			}
			if bestLikelihood-secondBestLikelihood < log10InformativeThreshold {
				bestPriority := bestAllele.priority()
				secondBestPriority := secondBestAllele.priority()
				for _, h := range haplotypes {
					if bestAllele == h {
						continue
					} else if likelihood := likelihoods.values[h][r]; bestLikelihood-likelihood > log10InformativeThreshold {
						// note: the test here is against bestLikelihood from the previous loop
						// it's not updated during the current loop!
						continue
					} else if priority := h.priority(); priority > bestPriority {
						secondBestAllele, secondBestPriority = bestAllele, bestPriority
						bestAllele, bestPriority = h, priority
					} else if priority > secondBestPriority {
						secondBestAllele, secondBestPriority = h, priority
					}
				}
			}
			seqString := aln.SEQ.AsString()
			cigar, alignmentOffset := runSmithWaterman(bestAllele.bases, seqString, 10, -15, -30, -5, softclip)
			if alignmentOffset < 0 {
				continue
			}

			var haplotypeCigar []sam.CigarOperation
			if last := len(bestAllele.cigar) - 1; bestAllele.cigar[last].Operation == 'M' {
				haplotypeCigar = make([]sam.CigarOperation, len(bestAllele.cigar))
				copy(haplotypeCigar, bestAllele.cigar)
				haplotypeCigar[last].Length += 1000
			} else {
				haplotypeCigar = make([]sam.CigarOperation, len(bestAllele.cigar)+1)
				copy(haplotypeCigar, bestAllele.cigar)
				haplotypeCigar[len(bestAllele.cigar)] = sam.CigarOperation{1000, 'M'}
			}

			var hapOffset, refOffset int32
		calcFirstMatchingBase:
			for _, ce := range haplotypeCigar {
				switch ce.Operation {
				case 'M', '=', 'X':
					if hapOffset >= alignmentOffset {
						break calcFirstMatchingBase
					}
					hapOffset += ce.Length
					refOffset += ce.Length
					if hapOffset > alignmentOffset {
						delta := hapOffset - alignmentOffset
						hapOffset -= delta
						refOffset -= delta
						break calcFirstMatchingBase
					}
				case 'I', 'S':
					hapOffset += ce.Length
				case 'D':
					refOffset += ce.Length
				}
			}
			readStartOnHaplotype := refOffset
			readStartOnReference := bestAllele.location + refOffset

			haplotypeCigarReadLength := sam.ReadLengthFromCigar(haplotypeCigar)

			haplotypeToRef := make([]sam.CigarOperation, 0, len(haplotypeCigar))
			var pos int32
			for _, ce := range haplotypeCigar {
				switch ce.Operation {
				case 'D':
					if pos >= alignmentOffset {
						haplotypeToRef = append(haplotypeToRef, ce)
					}
				default:
					if length := minInt32(pos+ce.Length, haplotypeCigarReadLength) - maxInt32(pos, alignmentOffset); length > 0 {
						haplotypeToRef = append(haplotypeToRef, sam.CigarOperation{length, ce.Operation})
					}
					pos += ce.Length
				}
			}
			for i := 1; i < len(haplotypeToRef); {
				if haplotypeToRef[i-1].Operation == haplotypeToRef[i].Operation {
					haplotypeToRef[i-1].Length += haplotypeToRef[i].Length
					haplotypeToRef = append(haplotypeToRef[:i], haplotypeToRef[i+1:]...)
				} else {
					i++
				}
			}

			readToRefCigarClean := applyCigarToCigar(cigar, haplotypeToRef)
			readToRefCigar := leftAlignIndel(readToRefCigarClean, refHaplotype.bases, seqString, readStartOnHaplotype, 0, true)
			leadingDeletions := sam.ReferenceLengthFromCigar(readToRefCigarClean) - sam.ReferenceLengthFromCigar(readToRefCigar)
			newAln := new(sam.Alignment)
			*newAln = *aln
			newAln.POS = int32(readStartOnReference) + leadingDeletions
			firstCE := aln.CIGAR[0]
			lastCE := aln.CIGAR[len(aln.CIGAR)-1]
			if firstCE.Operation == 'H' {
				newAln.CIGAR = append(aln.CIGAR[:1:1], readToRefCigar...)
			} else {
				newAln.CIGAR = readToRefCigar
			}
			if lastCE.Operation == 'H' {
				newAln.CIGAR = append(newAln.CIGAR, lastCE)
			}
			likelihoods.alns[r] = newAln
		}
	})
}
