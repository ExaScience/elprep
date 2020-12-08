// +build !pedantic

// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2020 imec vzw.

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

package filters

import (
	"math"
	"strconv"
)

func log10(x float64) float64 {
	return math.Log10(x)
}

func isOpticalDuplicateShort(tile1, tile2 tileInfo, opticalPixelDistance int) bool {
	return absInt(tile1.x-tile2.x) <= opticalPixelDistance && absInt(tile1.y-tile2.y) <= opticalPixelDistance
}

func (hc *HaplotypeCaller) Close() {}

// formatf prints floating point numbers the way Go does
func formatf(value float64, precision int) string {
	return strconv.FormatFloat(value, 'f', precision, 64)
}
