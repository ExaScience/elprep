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

package utils

import (
	"bufio"
	"io"
	"log"

	"github.com/exascience/elprep/v5/utils/bgzf"
)

// HandleBGZF checks if the given reader produces a gzip file
// by looking at the initial byte. It then either returns
// a bgzf.Reader, or returns the given reader unchanged.
// HandleBGZF uses ReadByte und UnreadByte.
func HandleBGZF(buf *bufio.Reader) io.Reader {
	if ok, err := bgzf.IsGzip(buf); err != nil {
		log.Panic(err)
		return nil
	} else if ok {
		if r, err := bgzf.NewReader(buf); err != nil {
			log.Panic(err)
			return nil
		} else {
			return r
		}
	} else {
		return buf
	}
}
