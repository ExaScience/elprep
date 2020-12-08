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
	"bufio"
	"log"
	"strings"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/elprep/v5/utils"
)

// ParseBed parses a BED file. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
func ParseBed(filename string) Bed {
	var bed Bed

	// open file
	file := internal.FileOpen(filename)
	defer internal.Close(file)

	scanner := bufio.NewScanner(utils.HandleBGZF(bufio.NewReader(file)))

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") ||
			strings.HasPrefix(line, "browser") {
			continue
		}
		data := strings.Split(line, "\t")
		chrom := utils.Intern(data[0])
		start := internal.ParseInt(data[1], 10, 32)
		end := internal.ParseInt(data[2], 10, 32)
		region := NewRegion(chrom, int32(start), int32(end), data[3:])
		AddRegion(&bed, region)
	}
	if err := scanner.Err(); err != nil {
		log.Panic(err)
	}
	// Make sure bed regions are sorted.
	sortRegions(bed)
	return bed
}
