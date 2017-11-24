package bed

import (
	"fmt"
	"sort"
	"strconv"

	"github.com/exascience/elprep/utils"
)

// Bed is a struct for representing the contents of a BED file. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
type Bed struct {
	// Bed tracks defined in the file.
	Tracks []*Track
	// Maps chromosome name onto bed regions.
	RegionMap map[utils.Symbol][]*Region
}

// A Track is a struct for representing BED tracks. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
type Track struct {
	// All track fields are optional.
	Fields map[string]string
	// The bed regions this track groups together.
	Regions []*Region
}

// A Region is a struct for representing intervals as defined in a BED
// file. See https://genome.ucsc.edu/FAQ/FAQformat.html#format1
type Region struct {
	Chrom          utils.Symbol
	Start          int32
	End            int32
	OptionalFields []interface{}
}

// Symbols for optional strand field of a Region.
var (
	// Strand forward.
	SF = utils.Intern("+")
	// Strand reverse.
	SR = utils.Intern("-")
)

// NewRegion allocates and initializes a new Region. Optional fields
// are given in order. If a "later" field is entered, then the
// "earlier" field was entered as well. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
func NewRegion(chrom utils.Symbol, start int32, end int32, fields []string) (b *Region, err error) {
	regionFields, err := initializeRegionFields(fields)
	if err != nil {
		return nil, err
	}
	return &Region{
		Chrom:          chrom,
		Start:          start,
		End:            end,
		OptionalFields: regionFields,
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

// Allocates a fresh SmallMap to initialize a Region's optional
// fields.
func initializeRegionFields(fields []string) ([]interface{}, error) {
	brFields := make([]interface{}, len(fields))
	for i, val := range fields {
		switch i {
		case brName:
			brFields[brName] = val
		case brScore:
			score, err := strconv.Atoi(val)
			if err != nil || score < 0 || score > 1000 {
				return nil, fmt.Errorf("invalid Score field : %v", err.Error())
			}
			brFields[brScore] = score
		case brStrand:
			if val != "+" && val != "-" {
				return nil, fmt.Errorf("invalid Strand field: %v", val)
			}
			brFields[brStrand] = utils.Intern(val)
		case brThickStart:
			start, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("invalid ThickStart field: %v", err.Error())
			}
			brFields[brThickStart] = start
		case brThickEnd:
			end, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("invalid ThickEnd field: %v", err.Error())
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
				return nil, fmt.Errorf("invalid BlockCount field: %v", err.Error())
			}
			brFields[brBlockCount] = count
		case brBlockSizes:
			sizes, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("invalid BlockSizes field: %v", err.Error())
			}
			brFields[brBlockSizes] = sizes
		case brBlockStarts:
			start, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("invalid BlockStarts field: %v", err.Error())
			}
			brFields[brBlockStarts] = start
		default:
			return nil, fmt.Errorf("invalid optional field: %v out of 0-8", val)
		}
	}
	return brFields, nil
}

// NewTrack allocates and initializes a new Track.
func NewTrack(fields map[string]string) *Track {
	return &Track{
		Fields: fields,
	}
}

// NewBed allocates and initializes an empty bed.
func NewBed() *Bed {
	return &Bed{
		RegionMap: make(map[utils.Symbol][]*Region),
	}
}

// AddRegion adds a region to the bed region map.
func AddRegion(bed *Bed, region *Region) {
	// append the region entry
	bed.RegionMap[region.Chrom] = append(bed.RegionMap[region.Chrom], region)
}

// A function for sorting the bed regions.
func sortRegions(bed *Bed) {
	for _, regions := range bed.RegionMap {
		sort.SliceStable(regions, func(i, j int) bool {
			return regions[i].Start < regions[j].Start
		})
	}
}
