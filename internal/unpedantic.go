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

package internal

import (
	"log"
	"math/rand"
)

const (
	// PedanticMode is a Boolean flag for conditional compilation
	PedanticMode = false

	// PedanticMessage can be added to the overall program message
	PedanticMessage = ""
)

type Rand = rand.Rand

// NewRand returns a Go-style random number generator.
func NewRand(seed int64) *Rand {
	return rand.New(rand.NewSource(seed))
}

// NewRandReflect works only in pedantic mode, so just panics unconditionally.
func NewRandReflect(seed int64) *Rand {
	log.Panic("NewRandReflect implemented only in pedantic mode.")
	return nil
}

// RandReflect works only in pedantic mode, so just panics unconditionally.
func RandReflect(r *Rand) int64 {
	log.Panic("RandReflect implemented only in pedantic mode.")
	return -1
}
