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

package vcf

import (
	"fmt"
	"strconv"

	"github.com/exascience/elprep/v4/utils"
)

// The supported VCF file format version.
const (
	FileFormatVersion           = "VCFv4.3"
	FileFormatVersionLine       = "##fileformat=VCFv4.3"
	fileFormatVersionLinePrefix = "##fileformat=VCFv4."
)

// DefaultHeaderColumns for VCF files.
var DefaultHeaderColumns = []string{"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"}

// Type is an enumeration type for different VCF field types
type Type uint

// The different VCF field types
const (
	InvalidType Type = iota
	Integer          // represented as int (not int32, since that's the same as rune in Go)
	Float            // represented as float64 (parsing as float32 seems problematic in some cases in Go)
	Flag             // represented as bool with fixed value true
	Character        // represented as rune
	String           // represented as string
)

// Constants for format information Number entries.
const (
	NumberA int32 = -1 * (1 + iota)
	NumberR
	NumberG
	NumberDot
	InvalidNumber
)

// Commonly used VCF entries.
var (
	END  = utils.Intern("END")
	GT   = utils.Intern("GT")
	PASS = utils.Intern("PASS")
)

type (
	// MetaInformation in VCF files.
	MetaInformation struct {
		ID          utils.Symbol
		Description string // "" if not present
		Fields      utils.StringMap
	}

	// FormatInformation in VCF files.
	FormatInformation struct {
		ID          utils.Symbol
		Description string // "" if not present
		Number      int32  // > InvalidNumber
		Type        Type
		Fields      utils.StringMap
	}

	// Header section of a VCF files.
	Header struct {
		FileFormat string
		Infos      []*FormatInformation
		Formats    []*FormatInformation
		Meta       map[string][]interface{} // string or *MetaInformation
		Columns    []string
	}

	// Variant line in a VCF file.
	Variant struct {
		Chrom          string
		Pos            int32
		ID             []string // nil/empty if missing
		Ref            string
		Alt            []string       // nil/empty if missing
		Qual           interface{}    // float64, or nil if missing
		Filter         []utils.Symbol // nil/empty if missing
		Info           utils.SmallMap // values are int, float64, bool, rune, string, or []interface{}
		GenotypeFormat []utils.Symbol
		GenotypeData   []utils.SmallMap // values are nil (for missing entry), int, float64, rune, string, or []interface{}
	}

	// Vcf represents the full contents of a VCF file.
	Vcf struct {
		Header   *Header
		Variants []*Variant
	}
)

// NewMetaInformation creates an empty instance.
func NewMetaInformation() *MetaInformation {
	return &MetaInformation{Fields: make(utils.StringMap)}
}

// NewFormatInformation creates an empty instance.
func NewFormatInformation() *FormatInformation {
	return &FormatInformation{Number: InvalidNumber, Fields: make(utils.StringMap)}
}

// NewHeader creates an empty instance.
func NewHeader() *Header {
	return &Header{Meta: make(map[string][]interface{})}
}

// Start returns the start position of a VCF line in the reference.
func (v *Variant) Start() int32 {
	return v.Pos
}

// End returns the end position of a VCF line in the reference, determined either by the END field or len(v.Ref)
func (v *Variant) End() (int32, error) {
	if end, ok := v.Info.Get(END); ok {
		switch e := end.(type) {
		case int:
			return int32(e), nil
		case string:
			i, err := strconv.ParseInt(e, 10, 32)
			if err != nil {
				return v.Pos + 1, err
			}
			v.Info.Set(END, int(i))
			return int32(i), nil
		default:
			return v.Pos + 1, fmt.Errorf("invalid END value %v", end)
		}
	}
	return v.Pos - 1 + int32(len(v.Ref)), nil
}

// Pass determines whether the variant passed all filters.
func (v *Variant) Pass() bool {
	return len(v.Filter) == 1 && v.Filter[0] == PASS
}
