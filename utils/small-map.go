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

// SmallMapEntry is an entry in a SmallMap.
type SmallMapEntry struct {
	Key   Symbol
	Value interface{}
}

// A SmallMap maps keys to values, similar to Go's built-in maps. A
// SmallMap can be more efficient in terms of memory and runtime
// performance than a native map if it has only few entries. SmallMap
// keys are always symbols.
type SmallMap []SmallMapEntry

// Get returns the first entry in the SmallMap that has the same key
// as the given key.
//
// It returns the found value and true if the key was found, otherwise
// nil and false.
func (m SmallMap) Get(key Symbol) (interface{}, bool) {
	for _, entry := range m {
		if entry.Key == key {
			return entry.Value, true
		}
	}
	return nil, false
}

// Set associates the given value with the given key.
//
// It does so by either setting the value of the first entry that has
// the same key as the given key, or else by appending a new key/value
// pair to the end of the SmallMap if no entry already has that key.
func (m *SmallMap) Set(key Symbol, value interface{}) {
	for index := range *m {
		if (*m)[index].Key == key {
			(*m)[index].Value = value
			return
		}
	}
	*m = append(*m, SmallMapEntry{key, value})
}

// Delete returns a SmallMap from which the first entry has been
// removed that has the same key as the given key.
//
// It also returns true if an entry was removed, and false if no entry
// was removed because there was no entry for the given key.
func (m SmallMap) Delete(key Symbol) (SmallMap, bool) {
	for index, entry := range m {
		if entry.Key == key {
			return append(m[:index], m[index+1:]...), true
		}
	}
	return m, false
}

// DeleteIf returns a SmallMap from which all entries have been
// removed that satisfy the given test.
//
// It also returns true if any entry was removed, and false if no
// entry was removed because no entry matched the given test.
func (m SmallMap) DeleteIf(test func(key Symbol, val interface{}) bool) (SmallMap, bool) {
	i := 0
	for _, entry := range m {
		if !test(entry.Key, entry.Value) {
			m[i] = entry
			i++
		}
	}
	return m[:i], i < len(m)
}
