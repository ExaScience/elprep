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

/*
	elPrep's pedantic mode emulates some behaviours that are strictly
	not necessary for correct results, but guarantee exactly equal results at
    the binary level. Normally, it's not necessary to be that pedantic.
*/

package filters

// #cgo CFLAGS: -std=c99 -Wno-unknown-pragmas
// #cgo LDFLAGS: -lm
// #include <fenv.h>
// #include <math.h>
//
// #pragma STDC FENV_ACCESS_ON
// double x86_log10(double x) {
//	 register double result;
//	 fesetround(FE_TONEAREST);
//   __asm __volatile__ ("fldlg2; fxch; fyl2x" : "=t" (result) : "0" (x) : "st(1)");
//   return result;
// }
import "C"

import (
	"fmt"
	"strconv"

	"github.com/exascience/elprep/v5/internal"
)

func log10(x float64) float64 {
	return float64(C.x86_log10(C.double(x)))
}

func isOpticalDuplicateShort(tile1, tile2 tileInfo, opticalPixelDistance int) bool {
	return absInt(int(int16(tile1.x))-int(int16(tile2.x))) <= opticalPixelDistance && absInt(int(int16(tile1.y))-int(int16(tile2.y))) <= opticalPixelDistance
}

func (hc *HaplotypeCaller) Close() {
	if hc.randomSeedFile != "" {
		rsf := internal.FileCreate(hc.randomSeedFile)
		defer internal.Close(rsf)
		fmt.Fprintln(rsf, internal.RandReflect(hc.random))
	}
}

// formatf prints floating point numbers the way Java does
func formatf(value float64, precision int) string {
	formatted := strconv.AppendFloat(nil, value, 'f', -1, 64)
	var offset int
	if formatted[0] == '-' {
		offset = 1
	}
	for i := offset; i < len(formatted); i++ {
		if formatted[i] == '.' {
			if end := i + 1 + precision; end < len(formatted) {
				if formatted[end] >= '5' {
					overflow := true
					for j := end - 1; j >= offset; j-- {
						if c := formatted[j]; c == '9' {
							formatted[j] = '0'
						} else if c != '.' {
							formatted[j] = c + 1
							overflow = false
							break
						}
					}
					if overflow {
						formatted = formatted[:end+1]
						copy(formatted[offset+1:], formatted[offset:])
						formatted[offset] = '1'
					} else {
						formatted = formatted[:end]
					}
				} else {
					formatted = formatted[:end]
				}
			} else {
				for j := len(formatted); j < end; j++ {
					formatted = append(formatted, '0')
				}
			}
			return string(formatted)
		}
	}
	formatted = append(formatted, '.')
	for i := 0; i < precision; i++ {
		formatted = append(formatted, '0')
	}
	return string(formatted)
}
