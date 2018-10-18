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

package utils

// A StringMap maps strings to strings.
type StringMap map[string]string

// Find returns the first index in a slice of StringMap where the
// predicate returns true, or -1 if predicate never returns true.
func Find(dict []StringMap, predicate func(record StringMap) bool) int {
	for index, record := range dict {
		if predicate(record) {
			return index
		}
	}
	return -1
}

// SetUniqueEntry checks if a mapping for the given key already exists
// in the StringMap. If this is the case, it returns false and the
// StringMap is not modified.  Otherwise, the given key/value pair is
// added to the StringMap.
func (record StringMap) SetUniqueEntry(key, value string) bool {
	if _, found := record[key]; found {
		return false
	}
	record[key] = value
	return true
}
