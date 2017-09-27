package sam

import (
	"math/rand"
	"strconv"

	"github.com/exascience/elprep/utils"
)

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

func FilterUnmappedReads(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool { return (aln.FLAG & Unmapped) == 0 }
}

func FilterUnmappedReadsStrict(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool {
		return ((aln.FLAG & Unmapped) == 0) && (aln.POS != 0) && (aln.RNAME != "*")
	}
}

func FilterDuplicateReads(_ *Header) AlignmentFilter {
	return func(aln *Alignment) bool { return (aln.FLAG & Duplicate) == 0 }
}

var sr = utils.Intern("sr")

func FilterOptionalReads(header *Header) AlignmentFilter {
	if _, found := header.UserRecords["@sr"]; found {
		delete(header.UserRecords, "@sr")
		return func(aln *Alignment) bool { _, found := aln.TAGS.Get(sr); return !found }
	} else {
		return nil
	}
}

func AddOrReplaceReadGroup(readGroup utils.StringMap) Filter {
	return func(header *Header) AlignmentFilter {
		header.RG = []utils.StringMap{readGroup}
		id := readGroup["ID"]
		return func(aln *Alignment) bool { aln.SetRG(id); return true }
	}
}

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
