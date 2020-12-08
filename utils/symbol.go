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
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/v5/internal"
)

type symbolName string

// A Symbol is a unique pointer to a string.
type Symbol *string

// SymbolHash computes a hash value for the given Symbol.
func SymbolHash(s Symbol) uint64 {
	return uint64(uintptr(unsafe.Pointer(s)))
}

func (s symbolName) Hash() uint64 {
	return internal.StringHash(string(s))
}

var symbolTable = sync.NewMap(0)

// Intern returns a Symbol for the given string.
//
// It always returns the same pointer for strings that are equal, and
// different pointers for strings that are not equal. So for two
// strings s1 and s2, if s1 == s2, then Intern(s1) == Intern(s2), and
// if s1 != s2, then Intern(s1) != Intern(s2).
//
// Dereferencing the pointer always yields a string that is equal to
// the original string: *Intern(s) == s always holds.
//
// It is safe for multiple goroutines to call Intern concurrently.
func Intern(s string) Symbol {
	entry, _ := symbolTable.LoadOrStore(symbolName(s), Symbol(&s))
	return entry.(Symbol)
}
