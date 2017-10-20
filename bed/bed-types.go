package bed

import (
	"fmt"
	"sort"
	"strconv"

	"github.com/exascience/elprep/utils"
)

/*
A struct for representing the contents of a BED file. See
https://genome.ucsc.edu/FAQ/FAQformat.html#format1
*/
type Bed struct {
	// Bed tracks defined in the file.
	Tracks []*BedTrack
	// Maps chromosome name onto bed regions.
	RegionMap map[utils.Symbol][]*BedRegion
}

/*
A struct for representing BED tracks. See
https://genome.ucsc.edu/FAQ/FAQformat.html#format1
*/
type BedTrack struct {
	// All track fields are optional.
	Fields map[string]string
	// The bed regions this track groups together.
	Regions []*BedRegion
}

/*
An interval as defined in a BED file. See
https://genome.ucsc.edu/FAQ/FAQformat.html#format1
*/
type BedRegion struct {
	Chrom          utils.Symbol
	Start          int32
	End            int32
	OptionalFields []interface{}
}

// Symbols for optional strand field of a BedRegion.
var (
	// Strand forward.
	SF = utils.Intern("+")
	// Strand reverse.
	SR = utils.Intern("-")
)

/*
Allocates and initializes a new BedRegion. Optional fields are given
in order. If a "later" field is entered, then the "earlier" field was
entered as well. See https://genome.ucsc.edu/FAQ/FAQformat.html#format1
*/
func NewBedRegion(chrom utils.Symbol, start int32, end int32, fields []string) (b *BedRegion, err error) {
	bedRegionFields, err := initializeBedRegionFields(fields)
	if err != nil {
		return nil, err
	}
	return &BedRegion{
		Chrom:          chrom,
		Start:          start,
		End:            end,
		OptionalFields: bedRegionFields,
	}, nil
}

// Valid bed region optional fields. See spec.
const (
	brName = iota
	brScore
	brStrand
	brThickStart
	brThickEnd
	brItemRgb
	brBlockCount
	brBlockSizes
	brBlockStarts
)

/*
Allocates a fresh SmallMap to initialize a BedRegion's optional fields.
*/
func initializeBedRegionFields(fields []string) ([]interface{}, error) {
	brFields := make([]interface{}, len(fields))
	for i, val := range fields {
		switch i {
		case brName:
			brFields[brName] = val
		case brScore:
			score, err := strconv.Atoi(val)
			if err != nil || score < 0 || score > 1000 {
				return nil, fmt.Errorf("Invalid Score field : %v", err)
			}
			brFields[brScore] = score
		case brStrand:
			if val != "+" && val != "-" {
				return nil, fmt.Errorf("Invalid Strand field: %v", val)
			}
			brFields[brStrand] = utils.Intern(val)
		case brThickStart:
			start, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid ThickStart field: %v", err)
			}
			brFields[brThickStart] = start
		case brThickEnd:
			end, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid ThickEnd field: %v", err)
			}
			brFields[brThickEnd] = end
		case brItemRgb:
			if val == "on" {
				brFields[brItemRgb] = true
			} else {
				brFields[brItemRgb] = false
			}
		case brBlockCount:
			count, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid BlockCount field: %v", err)
			}
			brFields[brBlockCount] = count
		case brBlockSizes:
			sizes, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid BlockSizes field: %v", err)
			}
			brFields[brBlockSizes] = sizes
		case brBlockStarts:
			start, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid BlockStarts field: %v", err)
			}
			brFields[brBlockStarts] = start
		default:
			return nil, fmt.Errorf("Invalid optional field: %v out of 0-8.", val)
		}
	}
	return brFields, nil
}

/*
Allocates and initializes a new BedTrack.
*/
func NewBedTrack(fields map[string]string) *BedTrack {
	return &BedTrack{
		Fields: fields,
	}
}

/*
Allocates and initializes an empty bed.
*/
func NewBed() *Bed {
	return &Bed{
		RegionMap: make(map[utils.Symbol][]*BedRegion),
	}
}

/*
Add region to the bed region map.
*/
func AddBedRegion(bed *Bed, region *BedRegion) {
	// append the region entry
	bed.RegionMap[region.Chrom] = append(bed.RegionMap[region.Chrom], region)
}

/*
A function for sorting the bed regions.
*/
func sortBedRegions(bed *Bed) {
	for _, regions := range bed.RegionMap {
		sort.SliceStable(regions, func(i, j int) bool {
			return regions[i].Start < regions[j].Start
		})
	}
}
