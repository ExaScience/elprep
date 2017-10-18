package bed

import (
	"fmt"
	"strconv"

	"github.com/exascience/elprep/utils"
)

type Bed struct {
	// Bed tracks defined in the file.
	Tracks []*BedTrack
	// Maps chromosome name onto bed regions.
	RegionMap map[utils.Symbol][]*BedRegion
}

type BedTrack struct {
	// All track fields are optional.
	Fields utils.SmallMap
	// The bed regions this track groups together.
	Regions []*BedRegion
}

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
Allocates and initializes a new BedRegion
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
Optional fields are given in order. If a "later" field is entered, then the
"earlier" field was entered as well. See spec.
*/
func initializeBedRegionFields(fields []string) ([]interface{}, error) {
	brFields := make([]interface{}, len(fields))
	for i, val := range fields {
		switch i {
		case brName:
			brFields[brName] = val
		case brScore:
			score, err := strconv.Atoi(val)
			if err != nil || score < 0 || score > 100 {
				return nil, fmt.Errorf("Invalid Score field : %v", err)
			}
			brFields[brScore] = score
		case brStrand:
			if val != "+" || val != "-" {
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
			return nil, fmt.Errorf("Invalid optional field: %v", val, " out of 0-8.")
		}
	}
	return brFields, nil
}

// Valid track line fields.
var (
	tName    = utils.Intern("name")
	tDescr   = utils.Intern("description")
	tPrior   = utils.Intern("priority")
	tColor   = utils.Intern("color")
	tUscore  = utils.Intern("useScore")
	tItemRgb = utils.Intern("intemRgb")
)

/*
Allocates a fresh SmallMap to initialize a BedTrack's optional fields.
Also checks that the given fields are valid.
*/
func initializeTrackFields(fields map[string]string) (utils.SmallMap, error) {
	trackFields := utils.SmallMap{}
	for key, val := range fields {
		switch key {
		case "name":
			trackFields.Set(tName, val)
		case "description":
			trackFields.Set(tDescr, val)
		case "priority":
			prior, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid Priority: %v", err)
			}
			trackFields.Set(tPrior, prior)
		case "color":
			trackFields.Set(tColor, val)
		case "usescore":
			score, err := strconv.Atoi(val)
			if err != nil {
				return nil, fmt.Errorf("Invalid UseScore: %v", err)
			}
			trackFields.Set(tUscore, score)
		case "itemRgb":
			rgb := false
			if val == "on" {
				rgb = true
			}
			trackFields.Set(tItemRgb, rgb)
		default:
			return nil, fmt.Errorf("Invalid track field: %v", val)
		}
	}
	return trackFields, nil
}

/*
Allocates and initializes a new BedTrack.
*/
func NewBedTrack(fields map[string]string) (*BedTrack, error) {
	trackFields, err := initializeTrackFields(fields)
	if err != nil {
		return nil, err
	}
	return &BedTrack{
		Fields: trackFields,
	}, nil
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
