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

package sam

import (
	"fmt"
)

// A stringScanner can be used scan/parse ASCII strings representing
// lines in SAM files.
//
// The zero stringScanner is valid and empty.
type stringScanner struct {
	index int
	data  []byte
	err   error
}

// reset resets the scanner, and initializes it with the given string.
func (sc *stringScanner) reset(s []byte) {
	sc.index = 0
	sc.data = s
	sc.err = nil
}

// len returns the number of ASCII characters that still need to be
// scanned/parsed. Returns 0 if Err() would return a non-nil value.
func (sc *stringScanner) len() int {
	if sc.err != nil {
		return 0
	}
	return len(sc.data) - sc.index
}

func (sc *stringScanner) readByteUntil(c byte) (b byte, found bool) {
	if sc.err != nil {
		return 0, false
	}
	start := sc.index
	next := start + 1
	if next >= len(sc.data) {
		sc.index = len(sc.data)
		return sc.data[start], false
	} else if sc.data[next] != c {
		if sc.err == nil {
			sc.err = fmt.Errorf("unexpected character %v in stringScanner.ReadByteUntil", sc.data[next])
		}
		return 0, false
	} else {
		sc.index = next + 1
		return sc.data[start], true
	}
}

func (sc *stringScanner) readUntil(c byte) (s string, found bool) {
	if sc.err != nil {
		return "", false
	}
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if sc.data[end] == c {
			sc.index = end + 1
			return string(sc.data[start:end]), true
		}
	}
	sc.index = len(sc.data)
	return string(sc.data[start:]), false
}

func (sc *stringScanner) readBytes() []byte {
	if sc.err != nil {
		return nil
	}
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if sc.data[end] == '\t' {
			sc.index = end + 1
			return append([]byte(nil), sc.data[start:end]...)
		}
	}
	sc.index = len(sc.data)
	return append([]byte(nil), sc.data[start:]...)
}

func (sc *stringScanner) readUntil2(c1, c2 byte) (s string, b byte) {
	if sc.err != nil {
		return "", 0
	}
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if c := sc.data[end]; (c == c1) || (c == c2) {
			sc.index = end + 1
			return string(sc.data[start:end]), c
		}
	}
	sc.index = len(sc.data)
	return string(sc.data[start:]), 0
}
