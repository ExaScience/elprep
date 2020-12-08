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

package internal

import (
	"log"
	"strconv"
)

// ParseInt is strconv.ParseInt with panics in place of errors
func ParseInt(s string, base, bitSize int) int64 {
	result, err := strconv.ParseInt(s, base, bitSize)
	if err != nil {
		log.Panic(err)
	}
	return result
}

// ParseUint is strconv.ParseUint with panics in place of errors
func ParseUint(s string, base, bitSize int) uint64 {
	result, err := strconv.ParseUint(s, base, bitSize)
	if err != nil {
		log.Panic(err)
	}
	return result
}

// ParseFloat is strconv.ParseFloat with panics in place of errors
func ParseFloat(s string, bitSize int) float64 {
	result, err := strconv.ParseFloat(s, bitSize)
	if err != nil {
		log.Panic(err)
	}
	return result
}
