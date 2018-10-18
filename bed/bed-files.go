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

package bed

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/exascience/elprep/v4/utils"
)

// Helper function for parsing a track line field.
func splitTrackField(field string) (string, string) {
	split := strings.Split(field, "=")
	return split[0], split[1]
}

// ParseBed parses a BED file. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
func ParseBed(filename string) (b *Bed, err error) {

	bed := NewBed()

	// open file
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()

	scanner := bufio.NewScanner(file)

	var track *Track // for storing the current track

	for scanner.Scan() {
		line := scanner.Text()
		data := strings.Split(line, "\t")
		// check if the line is a new track
		if data[0] == "track" {
			// create new track, store the old one
			if track != nil {
				bed.Tracks = append(bed.Tracks, track)
			}
			// all track entries are optional
			// parse and collect those that are used
			fields := make(map[string]string)
			for _, field := range data[1:] {
				key, val := splitTrackField(field)
				fields[key] = val
			}
			track = NewTrack(fields)
		} else {
			// parse a region entry
			chrom := utils.Intern(data[0])
			var err error
			start, err := strconv.Atoi(data[1])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region start: %v ", err)
			}
			end, err := strconv.Atoi(data[2])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region end: %v ", err)
			}
			region, err := NewRegion(chrom, int32(start), int32(end), data[3:])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region: %v ", err)
			}
			AddRegion(bed, region)
			if track != nil {
				track.Regions = append(track.Regions, region)
			}
		}

	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error while reading bed file: %v ", err)
	}
	// Make sure bed regions are sorted.
	sortRegions(bed)
	return bed, nil
}
