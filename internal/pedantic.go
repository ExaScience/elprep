// +build pedantic

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

const (
	// PedanticMode is a Boolean flag for conditional compilation
	PedanticMode = true

	// PedanticMessage can be added to the overall program message
	PedanticMessage = "pedantic mode "
)

// Rand produces random numbers,
// mimicking the behavior of the Java standard library.
type Rand struct {
	seed int64
}

const (
	multiplier = 0x5DEECE66D
	addend     = 0xB
	bits       = 31
	mask       = (1 << 48) - 1
)

// NewRand returns a Java-style random number generator.
func NewRand(seed int64) *Rand {
	return &Rand{seed: (seed ^ multiplier) & mask}
}

// NewRandReflect returns a Java-style random number generator
// without modifying the passed seed, but using it exactly as is.
func NewRandReflect(seed int64) *Rand {
	return &Rand{seed: seed}
}

// RandReflect returns the current seed value exactly as is.
func RandReflect(r *Rand) int64 {
	return r.seed
}

// Int31 produces the next int32.
func (r *Rand) Int31() int32 {
	r.seed = (r.seed*multiplier + addend) & mask
	b := uint(48 - bits)
	return int32((r.seed >> b) + (2 << ^b))
}

// Int31n produces the next int32 bounded by n.
func (r *Rand) Int31n(n int32) int32 {
	l := r.Int31()
	m := n - 1
	if (n & m) == 0 {
		l = int32((int(n) * int(l)) >> 31)
	} else {
		u := l
		for {
			l = u % n
			if u-l+m >= 0 {
				break
			}
			u = r.Int31()
		}
	}
	return l
}
