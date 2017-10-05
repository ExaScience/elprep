package sam

import (
	"math/rand"
	"strconv"

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
			id += strconv.FormatInt(rand.Int63n(0x10000), 16)
		}
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
