package sam

import (
	"log"
	"math/rand"
	"strconv"
	"strings"

	"github.com/exascience/elprep/bed"
	"github.com/exascience/elprep/utils"
)

/*
A filter for replacing the reference sequence dictionary in a Header.
*/
func ReplaceReferenceSequenceDictionary(dict []utils.StringMap) Filter {
	return func(header *Header) AlignmentFilter {
		if sortingOrder := header.HD["SO"]; sortingOrder == "coordinate" {
			previousPos := -1
			oldDict := header.SQ
			for _, entry := range dict {
				sn := entry["SN"]
				pos := utils.Find(oldDict, func(entry utils.StringMap) bool { return entry["SN"] == sn })
				if pos >= 0 {
					if pos > previousPos {
						previousPos = pos
					} else {
						header.SetHD_SO("unknown")
						break
					}
				}
			}
		}
		dictTable := make(map[string]bool)
		for _, entry := range dict {
			dictTable[entry["SN"]] = true
		}
		header.SQ = dict
		return func(aln *Alignment) bool { return dictTable[aln.RNAME] }
	}
}

/*
A filter for replacing the reference sequence dictionary in a Header
with one parsed from the given SAM/DICT file.
*/
func ReplaceReferenceSequenceDictionaryFromSamFile(samFile string) (f Filter, err error) {
	input, err := Open(samFile, true)
	if err != nil {
		return nil, err
	}
	defer func() {
		nerr := input.Close()
		if err == nil {
			err = nerr
		}
	}()
	header, _, err := ParseHeader(input.Reader)
	if err != nil {
		return nil, err
	}
	return ReplaceReferenceSequenceDictionary(header.SQ), nil
}

/*
A filter for removing unmapped sam-alignment instances, based on FLAG.
*/
func FilterUnmappedReads(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool { return (aln.FLAG & Unmapped) == 0 }
}

/*
A filter for removing unmapped sam-alignment instances, based on FLAG,
or POS=0, or RNAME=*.
*/
func FilterUnmappedReadsStrict(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool {
		return ((aln.FLAG & Unmapped) == 0) && (aln.POS != 0) && (aln.RNAME != "*")
	}
}

/*
A filter that removes all reads that are not exact matches with the reference (soft-clipping ok),
based on CIGAR string (only M and S allowed).
*/
func FilterNonExactMappingReads(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool { return !strings.ContainsAny(aln.CIGAR, "IDNHPX=") }
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

/*
A filter that removes all reads that are not exact matches with the reference,
based on the optional fields X0=1 (unique mapping), X1=0 (no suboptimal hit),
XM=0 (no mismatch), XO=0 (no gap opening), XG=0 (no gap extension).
*/
func FilterNonExactMappingReadsStrict(header *Header) AlignmentFilter {

	return func(aln *Alignment) bool {
		if x0, ok := aln.TAGS.Get(X0); !ok || x0.(int32) != 1 {
			return false
		}
		if x1, ok := aln.TAGS.Get(X1); !ok || x1.(int32) != 0 {
			return false
		}
		if xm, ok := aln.TAGS.Get(XM); !ok || xm.(int32) != 0 {
			return false
		}
		if xo, ok := aln.TAGS.Get(XO); !ok || xo.(int32) != 0 {
			return false
		}
		if xg, ok := aln.TAGS.Get(XG); !ok || xg.(int32) != 0 {
			return false
		}
		return true
	}
}

/*
A filter for removing duplicate sam-alignment instances, based on
FLAG.
*/
func FilterDuplicateReads(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool { return (aln.FLAG & Duplicate) == 0 }
}

var sr = utils.Intern("sr")

/*
A filter for removing alignments that represent optional information
in elPrep.
*/
func FilterOptionalReads(header *Header) AlignmentFilter {
	if _, found := header.UserRecords["@sr"]; found {
		delete(header.UserRecords, "@sr")
		return func(aln *Alignment) bool { _, found := aln.TAGS.Get(sr); return !found }
	} else {
		return nil
	}
}

/*
A filter for adding or replacing the read group both in the Header and
each Alignment.
*/
func AddOrReplaceReadGroup(readGroup utils.StringMap) Filter {
	return func(header *Header) AlignmentFilter {
		header.RG = []utils.StringMap{readGroup}
		id := readGroup["ID"]
		return func(aln *Alignment) bool { aln.SetRG(id); return true }
	}
}

/*
A filter for adding a @PG tag to a Header, and ensuring that it is the
first one in the chain.
*/
func AddPGLine(newPG utils.StringMap) Filter {
	return func(header *Header) AlignmentFilter {
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

/*
A filter for prepending "chr" to the reference sequence names in a
Header, and in RNAME and RNEXT in each Alignment.
*/
func RenameChromosomes(header *Header) AlignmentFilter {
	for _, entry := range header.SQ {
		if sn, found := entry["SN"]; found {
			entry["SN"] = "chr" + sn
		}
	}
	return func(aln *Alignment) bool {
		if (aln.RNAME != "=") && (aln.RNAME != "*") {
			aln.RNAME = "chr" + aln.RNAME
		}
		if (aln.RNEXT != "=") && (aln.RNEXT != "*") {
			aln.RNEXT = "chr" + aln.RNEXT
		}
		return true
	}
}

/*
A filter for adding the refid (index in the reference sequence
dictionary) to alignments as temporary values.
*/
func AddREFID(header *Header) AlignmentFilter {
	dictTable := make(map[string]int32)
	for index, entry := range header.SQ {
		dictTable[entry["SN"]] = int32(index)
	}
	return func(aln *Alignment) bool {
		value, found := dictTable[aln.RNAME]
		if !found {
			value = -1
		}
		aln.SetREFID(value)
		return true
	}
}

/*
A filter for removing optional fields in an alignment.
*/
func RemoveOptionalFields(tags []string) Filter {
	if len(tags) == 0 {
		return nil
	} else {
		// Intern the tags once.
		var optionals []utils.Symbol
		for _, tag := range tags {
			optionals = append(optionals, utils.Intern(tag))
		}
		return func(header *Header) AlignmentFilter {
			return func(aln *Alignment) bool {
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
}

/*
A filter for removing all but a list of given optional fields in an alignment.
*/
func KeepOptionalFields(tags []string) Filter {
	if len(tags) == 0 {
		return func(header *Header) AlignmentFilter {
			return func(aln *Alignment) bool {
				aln.TAGS = nil
				return true
			}
		}
	} else {
		// Intern the tags once.
		var optionals []utils.Symbol
		for _, tag := range tags {
			optionals = append(optionals, utils.Intern(tag))
		}
		return func(header *Header) AlignmentFilter {
			return func(aln *Alignment) bool {
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
}

/*
A filter for soft-clipping an alignment at the end of a reference
sequence, and set MAPQ to 0 if unmapped.
*/
func CleanSam(header *Header) AlignmentFilter {
	referenceSequenceTable := make(map[string]int32)
	for _, sn := range header.SQ {
		referenceSequenceTable[sn["SN"]], _ = SQ_LN(sn)
	}
	return func(aln *Alignment) bool {
		if aln.IsUnmapped() {
			aln.MAPQ = 0
		} else if cigar, err := ScanCigarString(aln.CIGAR); err != nil {
			log.Fatal(err.Error(), ", while scanning a CIGAR string for ", aln.QNAME, " in CleanSam")
		} else if length := referenceSequenceTable[aln.RNAME]; end(aln, cigar) > length {
			clipFrom := length - aln.POS + 1
			aln.CIGAR = softClipEndOfRead(clipFrom, cigar)
		}
		return true
	}
}

/*
A filter for removing all reads that do not overlap with a set
of regions specified by a bed file.
*/
func FilterNonOverlappingReads(bed *bed.Bed) Filter {
	return func(header *Header) AlignmentFilter {
		return func(aln *Alignment) bool {
			alnStart := aln.POS
			cigar, err := ScanCigarString(aln.CIGAR)
			if err != nil {
				log.Fatal(err.Error(), ", while scanning a CIGAR string for", aln.QNAME, " in FilterNonOverlappingReads")
			}
			var alnEnd int32
			if aln.IsUnmapped() || readLengthFromCigar(cigar) <= 0 {
				alnEnd = aln.POS
			} else {
				alnEnd = end(aln, cigar)
			}
			regions := bed.RegionMap[utils.Intern(aln.RNAME)]
			left := 0
			right := len(regions) - 1
			for left <= right {
				mid := (left + right) / 2
				regionStart := regions[mid].Start
				regionEnd := regions[mid].End
				if regionStart > alnEnd-1 {
					right = mid - 1
				} else if regionEnd <= alnStart-1 {
					left = mid + 1
				} else {
					return true
				}
			}
			return false
		}
	}
}
