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

package vcf

// A StringScanner can be used scan/parse strings representing
// lines in VCF files.
//
// The zero StringScanner is valid and empty.
type StringScanner struct {
	index int
	data  string
}

// Reset resets the scanner, and initializes it with the given string.
func (sc *StringScanner) Reset(s string) {
	sc.index = 0
	sc.data = s
}

// Len returns the number of ASCII characters that still need to be
// scanned/parsed.
func (sc *StringScanner) Len() int {
	return len(sc.data) - sc.index
}

// SkipSpace skips ' ' runes
func (sc *StringScanner) SkipSpace() {
	for end := sc.index; end < len(sc.data); end++ {
		if sc.data[end] != ' ' {
			sc.index = end
			return
		}
	}
	sc.index = len(sc.data)
}

func (sc *StringScanner) readUntilByte(c byte) (s string, found bool) {
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if sc.data[end] == c {
			sc.index = end + 1
			return sc.data[start:end], true
		}
	}
	sc.index = len(sc.data)
	return sc.data[start:], false
}

func (sc *StringScanner) readUntilBytes(bytes []byte) string {
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		c := sc.data[end]
		for _, b := range bytes {
			if c == b {
				sc.index = end
				return sc.data[start:end]
			}
		}
	}
	sc.index = len(sc.data)
	return sc.data[start:]
}
