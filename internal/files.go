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

package internal

import (
	"os"
	"path/filepath"
)

// Directory returns a slice of filenames. If the given filename
// refers to a directory, return a slice of names of files that are in
// this directory. If the given filename does not refer to a
// directory, return a slice with this filename as the only entry.
func Directory(file string) (files []string, err error) {
	info, err := os.Stat(file)
	if err != nil {
		return nil, err
	}
	if !info.IsDir() {
		return []string{filepath.Base(file)}, nil
	}
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer func() {
		nerr := f.Close()
		if err == nil {
			err = nerr
		}
	}()
	return f.Readdirnames(0)
}
