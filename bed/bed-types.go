// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2017-2020 imec vzw.

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

package bed

import (
	"log"
	"sort"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/utils"
)

// Bed represents the contents of a BED file. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
// It maps chromosome names onto bed regions.
// utils.Symbol -> []*Region
type Bed = utils.SmallMap

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
func NewRegion(chrom utils.Symbol, start int32, end int32, fields []string) *Region {
	return &Region{
		Chrom:          chrom,
		Start:          start,
		End:            end,
		OptionalFields: initializeRegionFields(fields),
	}
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
func initializeRegionFields(fields []string) []interface{} {
	brFields := make([]interface{}, len(fields))
	for i, val := range fields {
		switch i {
		case brName:
			brFields[brName] = val
		case brScore:
			score := internal.ParseInt(val, 10, 64)
			if score < 0 || score > 1000 {
				log.Panicf("invalid Score field : %v", score)
			}
			brFields[brScore] = int(score)
		case brStrand:
			if val != "+" && val != "-" {
				log.Panicf("invalid Strand field: %v", val)
			}
			brFields[brStrand] = utils.Intern(val)
		case brThickStart:
			brFields[brThickStart] = int(internal.ParseInt(val, 10, 64))
		case brThickEnd:
			brFields[brThickEnd] = int(internal.ParseInt(val, 10, 64))
		case brItemRgb:
			if val == "on" {
				brFields[brItemRgb] = true
			} else {
				brFields[brItemRgb] = false
			}
		case brBlockCount:
			brFields[brBlockCount] = int(internal.ParseInt(val, 10, 64))
		case brBlockSizes:
			brFields[brBlockSizes] = int(internal.ParseInt(val, 10, 64))
		case brBlockStarts:
			brFields[brBlockStarts] = int(internal.ParseInt(val, 10, 64))
		default:
			log.Panicf("invalid optional field: %v out of 0-8", val)
		}
	}
	return brFields
}

// AddRegion adds a region to the bed region map.
func AddRegion(bed *Bed, region *Region) {
	// append the region entry
	if regions, found := bed.Get(region.Chrom); found {
		bed.Set(region.Chrom, append(regions.([]*Region), region))
	} else {
		bed.Set(region.Chrom, []*Region{region})
	}
}

// A function for sorting the bed regions.
func sortRegions(bed Bed) {
	for _, regions := range bed {
		regionsSlice := regions.Value.([]*Region)
		sort.SliceStable(regionsSlice, func(i, j int) bool {
			return regionsSlice[i].Start < regionsSlice[j].Start
		})
	}
}
