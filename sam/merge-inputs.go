package sam

import (
	"log"
	"path"
	"strconv"
	"time"

	"github.com/exascience/elprep/v5/utils"
)

//mergeInputs handles merging of sam/bam files that are passed as inputs to the elPrep sfm command.
//It reads the headers from the files and attempts merging them into a new header object. If succesful, mergeInputs
//returns a list of filters that can be passed to the elPrep execution engine that take care of updating the header,
// correcting optional read tags that refer to the header, etc while executing an elPrep pipeline on the given inputs.
func MergeInputs(inputDir string, files []string) (*Header, []AlignmentFilter) {
	log.Println("Deteced multiple input files... attempting to merge headers")
	var newHeader *Header
	updateFilters := []AlignmentFilter{}
	for _, fileName := range files {
		sam := Open(path.Join(inputDir, fileName))
		header := sam.ParseHeader()
		sam.Close()
		if newHeader == nil {
			newHeader = header
			continue
		}
		filters := mergeHeaders(newHeader, header)
		for _, filter := range filters {
			updateFilters = append(updateFilters, filter)
		}
	}
	return newHeader, updateFilters
}

//mergeHeaders combines two headers into a new header, respecting the constraints set by the SAM specification.
//These constraints include:
//1. Order of the sequence dictionaries (@sq) must be kept. The sequence dictionaries must be combined into a new
//sequence dictionary in such a way that the order of the original sequence dictionaries is not modified. Any conflicts
//should result in an error. This is because the sequence dictionary determines the alignment sorting order and
//downstream tools may depend on that order.
//2. Read group identifiers must be unique (@rg). When merging read group lines, the IDs must be unique and collisions
//must be resolved by changing IDs to be unique. The optional @rg tags in the alignments must be updated accordingly.
//3. Program identifiers must be unique (@pg). If the IDs of @pg lines are modified to resolve collisions when merging, the
//PP tag must be updated accordingly, as well as the optional @pg tags in the alignments that refer to this ID.
//4. Comment lines must be merged (any order).
//5. When merging BAM files, the order in the header must be respected (GO entry). This may require either resorting the
//merged file or at least setting the order of the merged file to unknown.
func mergeHeaders(toHeader, fromHeader *Header) []AlignmentFilter {
	readUpdateFilters := []AlignmentFilter{}
	mergeMetaData(toHeader, fromHeader)
	mergeSequenceDictionariesOfHeaders(toHeader, fromHeader)
	filters := mergeReadGroups(toHeader, fromHeader)
	for _, filter := range filters {
		readUpdateFilters = append(readUpdateFilters, filter)
	}
	filters = mergeProgramLines(toHeader, fromHeader)
	for _, filter := range filters {
		readUpdateFilters = append(readUpdateFilters, filter)
	}
	mergeCommentLines(toHeader, fromHeader)
	mergeUserRecords(toHeader, fromHeader)
	toHeader.SetHDSO(Unknown)
	return readUpdateFilters
}

func mergeMetaData(toHeader, fromHeader *Header) {
	for k, v := range fromHeader.HD {
		toHeader.HD[k] = v
	}
}

func mergeUserRecords(toHeader, fromHeader *Header) {
	for k, v := range fromHeader.UserRecords {
		toHeader.UserRecords[k] = v
	}
}

func memberOfSequenceDictionary(dict []utils.StringMap, seq string) (int, bool) {
	for i, sq := range dict {
		if sq["SN"] == seq {
			return i, true
		}
	}
	return -1, false
}

//mergeSequenceDictionaries merges two sequence dictionaries while respecting the order expressed in both sequence
//dictionaries. The algorithm creates a new dictionary S3 and proceeds as follows:
//Iterate over all seq in dictionary S1:
//case 1: seq not in S2 and seq not in S3: S3 = S3 + seq
//case 2: seq in S2 and seq not in S3: merge (S3, S2[...index of seq in S2]), then S3 = S3 + seq
//case 3: seq in S2 and seq in S3: if index of seq in S3 == index of seq in S2 then ok else cannot merge
func mergeSequenceDictionaries(toDict, fromDict []utils.StringMap) ([]utils.StringMap, bool) {
	newDict := []utils.StringMap{}
	j := 0
	for iTo, sq := range toDict {
		sn := sq["SN"]
		if iFrom, ok := memberOfSequenceDictionary(fromDict, sn); ok {
			if iNew, ok := memberOfSequenceDictionary(newDict, sn); ok {
				if iNew != iTo {
					log.Panicln("Cannot merge sequence dictionaries.")
					return nil, false
				}
				j = iNew + 1
			} else {
				//merge everything between index j and index iFrom into the new sequence dictionary
				newDict, ok = mergeSequenceDictionaries(newDict, fromDict[j:iFrom])
				if !ok {
					log.Panicln("Cannot merge sequence dictionaries.")
					return nil, false
				}
				newDict = append(newDict, sq)
				j = iFrom + 1
			}
		} else {
			newDict = append(newDict, sq)
		}
	}
	for i := j; i < len(fromDict); i++ {
		newDict = append(newDict, fromDict[i])
	}
	return newDict, true
}

//mergeSequenceDictionariesOfHeaders merges the sequence dictionaries of two headers.
func mergeSequenceDictionariesOfHeaders(toHeader, fromHeader *Header) {
	var ok bool
	toHeader.SQ, ok = mergeSequenceDictionaries(toHeader.SQ, fromHeader.SQ)
	if !ok {
		log.Panicln("Cannot merge sequence dictionaries of given SAM headers.")
	}
}

//mergeReadGroups merges the read group lines from two headers. The read group IDs must be unique. If they are not,
//they are replaced to be made unique.
func mergeReadGroups(toHeader, fromHeader *Header) []AlignmentFilter {
	filters := []AlignmentFilter{}
	toHeaderIDS := map[string]bool{}
	for _, record := range toHeader.RG {
		id := record["ID"]
		toHeaderIDS[id] = true
	}
	for _, record := range fromHeader.RG {
		id := record["ID"]
		if _, ok := toHeaderIDS[id]; ok {
			//Read group ID not unique, make it unique
			newID := id + strconv.FormatInt(time.Now().Unix(), 10)
			record["ID"] = newID
			filters = append(filters, replaceReadGroupIfExists(id, newID))
		}
		toHeader.RG = append(toHeader.RG, record)
	}
	return filters
}

//replaceReadGroupIfExists replaces the RG optional tag of an alignment, if there was an RG tag filled in.
func replaceReadGroupIfExists(oldID, newID string) AlignmentFilter {
	return func(aln *Alignment) bool {
		if rID, ok := aln.TAGS.Get(RG); ok {
			if rID.(string) == oldID {
				aln.SetRG(newID)
			}
		}
		return true
	}
}

//mergeProgramLines merges the program header lines from two headers. The program IDs must be unique. If they are not,
//they are replaced to be made unique.
func mergeProgramLines(toHeader, fromHeader *Header) []AlignmentFilter {
	filters := []AlignmentFilter{}
	toHeaderPGIDs := map[string]bool{}
	replacedPGIDs := map[string]string{}
	for _, record := range toHeader.PG {
		id := record["ID"]
		toHeaderPGIDs[id] = true
	}
	for _, record := range fromHeader.PG {
		id := record["ID"]
		if toHeaderPGIDs[id] {
			//Program ID not unique, make it unique
			newID := id + strconv.FormatInt(time.Now().Unix(), 10)
			record["ID"] = newID
			filters = append(filters, replaceProgramIDIfExists(id, newID))
			replacedPGIDs[id] = newID
		}
		toHeader.PG = append(toHeader.PG, record)
	}
	for i := 1; i < len(fromHeader.RG); i++ {
		record := toHeader.RG[i]
		if pp, ok := record["PP"]; ok {
			if newID, ok := replacedPGIDs[pp]; ok {
				record["PP"] = newID
			}
		}
	}
	return filters
}

func replaceProgramIDIfExists(oldId, newID string) AlignmentFilter {
	return func(aln *Alignment) bool {
		if rid, ok := aln.TAGS.Get(PG); ok {
			if rid.(string) == oldId {
				aln.SetPG(newID)
			}
		}
		return true
	}
}

//mergeCommentLines combines the comment lines of two headers. There re no particular requirements in the SAM
//specification.
func mergeCommentLines(toHeader, fromHeader *Header) {
	//simply append the comments
	for _, comment := range fromHeader.CO {
		toHeader.CO = append(toHeader.CO, comment)
	}
}
