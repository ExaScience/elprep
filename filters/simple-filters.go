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
	"math"
	"math/rand"
	"strconv"

	"github.com/exascience/elprep/v4/bed"
	"github.com/exascience/elprep/v4/intervals"
	"github.com/exascience/elprep/v4/sam"
	"github.com/exascience/elprep/v4/utils"
)

// ReplaceReferenceSequenceDictionary returns a filter for replacing
// the reference sequence dictionary in a Header.
func ReplaceReferenceSequenceDictionary(dict []utils.StringMap) sam.Filter {
	return func(header *sam.Header) sam.AlignmentFilter {
		if sortingOrder := sam.SortingOrder(header.HD["SO"]); sortingOrder == sam.Coordinate {
			previousPos := -1
			oldDict := header.SQ
			for _, entry := range dict {
				sn := entry["SN"]
				pos := utils.Find(oldDict, func(entry utils.StringMap) bool { return entry["SN"] == sn })
				if pos >= 0 {
					if pos > previousPos {
						previousPos = pos
					} else {
						header.SetHDSO(sam.Unknown)
						break
					}
				}
			}
		}
		dictTable := make(map[string]bool)
		dictTable["*"] = true
		for _, entry := range dict {
			dictTable[entry["SN"]] = true
		}
		header.SQ = dict
		return func(aln *sam.Alignment) bool { return dictTable[aln.RNAME] }
	}
}

// ReplaceReferenceSequenceDictionaryFromSamFile returns a filter for
// replacing the reference sequence dictionary in a Header with one
// parsed from the given SAM/DICT file.
func ReplaceReferenceSequenceDictionaryFromSamFile(samFile string) (f sam.Filter, err error) {
	input, err := sam.Open(samFile)
	if err != nil {
		return nil, err
	}
	defer func() {
		nerr := input.Close()
		if err == nil {
			err = nerr
		}
	}()
	header, err := input.ParseHeader()
	if err != nil {
		return nil, err
	}
	return ReplaceReferenceSequenceDictionary(header.SQ), nil
}

// RemoveUnmappedReads is a filter for removing unmapped sam-alignment
// instances, based on FLAG.
func RemoveUnmappedReads(_ *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool { return (aln.FLAG & sam.Unmapped) == 0 }
}

// RemoveUnmappedReadsStrict is a filter for removing unmapped
// sam-alignment instances, based on FLAG, or POS=0, or RNAME=*.
func RemoveUnmappedReadsStrict(_ *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool {
		return ((aln.FLAG & sam.Unmapped) == 0) && (aln.POS != 0) && (aln.RNAME != "*")
	}
}

var nonExactMappingOperator = map[byte]bool{'I': true, 'D': true, 'N': true, 'H': true, 'P': true, 'X': true, '=': true}

// RemoveNonExactMappingReads is a filter that removes all reads that
// are not exact matches with the reference (soft-clipping ok), based
// on CIGAR string (only M and S allowed).
func RemoveNonExactMappingReads(_ *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool {
		for _, op := range aln.CIGAR {
			if nonExactMappingOperator[op.Operation] {
				return false
			}
		}
		return true
	}
}

// Symbols for optional fields used for determining exact matches. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
var (
	X0 = utils.Intern("X0")
	X1 = utils.Intern("X1")
	XM = utils.Intern("XM")
	XO = utils.Intern("XO")
	XG = utils.Intern("XG")
)

// RemoveNonExactMappingReadsStrict is a filter that removes all reads
// that are not exact matches with the reference, based on the
// optional fields X0=1 (unique mapping), X1=0 (no suboptimal hit),
// XM=0 (no mismatch), XO=0 (no gap opening), XG=0 (no gap extension).
func RemoveNonExactMappingReadsStrict(header *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool {
		if x0, ok := aln.TAGS.Get(X0); !ok || x0.(int64) != 1 {
			return false
		}
		if x1, ok := aln.TAGS.Get(X1); !ok || x1.(int64) != 0 {
			return false
		}
		if xm, ok := aln.TAGS.Get(XM); !ok || xm.(int64) != 0 {
			return false
		}
		if xo, ok := aln.TAGS.Get(XO); !ok || xo.(int64) != 0 {
			return false
		}
		if xg, ok := aln.TAGS.Get(XG); !ok || xg.(int64) != 0 {
			return false
		}
		return true
	}
}

// RemoveDuplicateReads is a filter for removing duplicate
// sam-alignment instances, based on FLAG.
func RemoveDuplicateReads(_ *sam.Header) sam.AlignmentFilter {
	return func(aln *sam.Alignment) bool { return (aln.FLAG & sam.Duplicate) == 0 }
}

var sr = utils.Intern("sr")

// RemoveOptionalReads is a filter for removing alignments that
// represent optional information in elPrep.
func RemoveOptionalReads(header *sam.Header) sam.AlignmentFilter {
	if _, found := header.UserRecords["@sr"]; found {
		delete(header.UserRecords, "@sr")
		return func(aln *sam.Alignment) bool { _, found := aln.TAGS.Get(sr); return !found }
	}
	return nil
}

// AddOrReplaceReadGroup returns a filter for adding or replacing the
// read group both in the Header and in each Alignment.
func AddOrReplaceReadGroup(readGroup utils.StringMap) sam.Filter {
	return func(header *sam.Header) sam.AlignmentFilter {
		header.RG = []utils.StringMap{readGroup}
		id := readGroup["ID"]
		return func(aln *sam.Alignment) bool { aln.SetRG(id); return true }
	}
}

// AddPGLine returns a filter for adding a @PG tag to a Header, and
// ensuring that it is the first one in the chain.
func AddPGLine(newPG utils.StringMap) sam.Filter {
	return func(header *sam.Header) sam.AlignmentFilter {
		id := newPG["ID"]
		for utils.Find(header.PG, func(entry utils.StringMap) bool { return entry["ID"] == id }) >= 0 {
			id += " "
			id += strconv.FormatInt(rand.Int63n(0x10000), 16)
		}
		newPG["ID"] = id
		for _, PG := range header.PG {
			nextID := PG["ID"]
			if pos := utils.Find(header.PG, func(entry utils.StringMap) bool { return entry["PP"] == nextID }); pos < 0 {
				newPG["PP"] = nextID
				break
			}
		}
		header.PG = append(header.PG, newPG)
		return nil
	}
}

// RenameChromosomes is a filter for prepending "chr" to the reference
// sequence names in a Header, and in RNAME and RNEXT in each
// Alignment.
func RenameChromosomes(header *sam.Header) sam.AlignmentFilter {
	for _, entry := range header.SQ {
		if sn, found := entry["SN"]; found {
			entry["SN"] = "chr" + sn
		}
	}
	return func(aln *sam.Alignment) bool {
		if (aln.RNAME != "=") && (aln.RNAME != "*") {
			aln.RNAME = "chr" + aln.RNAME
		}
		if (aln.RNEXT != "=") && (aln.RNEXT != "*") {
			aln.RNEXT = "chr" + aln.RNEXT
		}
		return true
	}
}

// AddREFID is a filter for adding the refid (index in the reference
// sequence dictionary) to alignments as temporary values.
func AddREFID(header *sam.Header) sam.AlignmentFilter {
	dictTable := make(map[string]int32)
	dictTable["*"] = -1
	for index, entry := range header.SQ {
		dictTable[entry["SN"]] = int32(index)
	}
	return func(aln *sam.Alignment) bool {
		value, found := dictTable[aln.RNAME]
		if !found {
			value = -1
		}
		aln.SetREFID(value)
		return true
	}
}

// RemoveOptionalFields returns a filter for removing optional fields
// in an alignment.
func RemoveOptionalFields(tags []string) sam.Filter {
	if len(tags) == 0 {
		return nil
	}
	// Intern the tags once.
	var optionals []utils.Symbol
	for _, tag := range tags {
		optionals = append(optionals, utils.Intern(tag))
	}
	return func(header *sam.Header) sam.AlignmentFilter {
		return func(aln *sam.Alignment) bool {
			aln.TAGS, _ = aln.TAGS.DeleteIf(func(key utils.Symbol, val interface{}) bool {
				for _, tag := range optionals {
					if tag == key {
						return true
					}
				}
				return false
			})
			return true
		}
	}
}

// KeepOptionalFields returns a filter for removing all but a list of
// given optional fields in an alignment.
func KeepOptionalFields(tags []string) sam.Filter {
	if len(tags) == 0 {
		return func(header *sam.Header) sam.AlignmentFilter {
			return func(aln *sam.Alignment) bool {
				aln.TAGS = nil
				return true
			}
		}
	}
	// Intern the tags once.
	var optionals []utils.Symbol
	for _, tag := range tags {
		optionals = append(optionals, utils.Intern(tag))
	}
	return func(header *sam.Header) sam.AlignmentFilter {
		return func(aln *sam.Alignment) bool {
			aln.TAGS, _ = aln.TAGS.DeleteIf(func(key utils.Symbol, val interface{}) bool {
				for _, tag := range optionals {
					if tag == key {
						return false
					}
				}
				return true
			})
			return true
		}
	}
}

// CleanSam is a filter for soft-clipping an alignment at the end of a
// reference sequence, and setting MAPQ to 0 if unmapped.
func CleanSam(header *sam.Header) sam.AlignmentFilter {
	referenceSequenceTable := make(map[string]int32)
	for _, sn := range header.SQ {
		referenceSequenceTable[sn["SN"]], _ = sam.SQLN(sn)
	}
	return func(aln *sam.Alignment) bool {
		if aln.IsUnmapped() {
			aln.MAPQ = 0
		} else if length := referenceSequenceTable[aln.RNAME]; end(aln, aln.CIGAR) > length {
			clipFrom := length - aln.POS + 1
			aln.CIGAR = softClipEndOfRead(clipFrom, aln.CIGAR)
		}
		return true
	}
}

// RemoveNonOverlappingReads returns a filter for removing all reads
// that do not overlap with a set of regions specified by a bed file.
func RemoveNonOverlappingReads(bed *bed.Bed) sam.Filter {
	ivals := intervals.FromBed(bed)
	for chrom, ival := range ivals {
		intervals.ParallelSortByStart(ival)
		ivals[chrom] = intervals.ParallelFlatten(ival)
	}
	return func(header *sam.Header) sam.AlignmentFilter {
		return func(aln *sam.Alignment) bool {
			alnStart := aln.POS
			alnEnd := aln.POS
			if !aln.IsUnmapped() {
				if readLengthFromCigar(aln.CIGAR) > 0 {
					alnEnd = end(aln, aln.CIGAR)
				}
			}
			return intervals.Overlap(ivals[aln.RNAME], alnStart, alnEnd)
		}
	}
}

// RemoveMappingQualityLessThan is a filter for removing reads
// that do not match or exceed the given mapping quality.
func RemoveMappingQualityLessThan(mq int) sam.Filter {
	if mq == 0 {
		return nil // no need to add any filter because aln.MAPQ is always >= 0
	}
	if mq > math.MaxUint8 {
		return func(_ *sam.Header) sam.AlignmentFilter {
			return func(_ *sam.Alignment) bool {
				return false // no aln.MAPQ can be > math.MaxUint8
			}
		}
	}
	mapq := byte(mq)
	return func(_ *sam.Header) sam.AlignmentFilter {
		return func(aln *sam.Alignment) bool { return aln.MAPQ >= mapq }
	}
}
