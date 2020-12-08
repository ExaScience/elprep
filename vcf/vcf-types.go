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

import (
	"log"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/utils"
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

	// Genotype is a structured representation of the GT entry in a VCF file.
	Genotype struct {
		Phased bool
		GT     []int32        // < 0 for unknown entries
		Data   utils.SmallMap // values are nil (for missing entry), int, float64, rune, string, or []interface{}
	}

	// Variant line in a VCF file.
	Variant struct {
		Source         string // this is not part of the VCF spec, but is needed in HaplotypeCaller
		Chrom          string
		Pos            int32    // < 0 if unknown
		ID             []string // nil/empty if missing
		Ref            string
		Alt            []string       // nil/empty if missing
		Qual           interface{}    // float64, or nil if missing
		Filter         []utils.Symbol // nil/empty if missing
		Info           utils.SmallMap // values are int, float64, bool, rune, string, or []interface{}
		GenotypeFormat []utils.Symbol
		GenotypeData   []Genotype
	}

	// Vcf represents the full contents of a VCF file.
	Vcf struct {
		Header   *Header
		Variants []Variant
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
	return &Header{
		FileFormat: FileFormatVersionLine,
		Meta:       make(map[string][]interface{}),
		Columns:    DefaultHeaderColumns,
	}
}

// Start returns the start position of a VCF line in the reference.
func (v Variant) Start() int32 {
	return v.Pos
}

// End returns the end position of a VCF line in the reference, determined either by the END field or len(v.Ref)
func (v *Variant) End() int32 {
	if end, ok := v.Info.Get(END); ok {
		switch e := end.(type) {
		case int:
			return int32(e)
		case string:
			i := internal.ParseInt(e, 10, 32)
			v.Info.Set(END, int(i))
			return int32(i)
		default:
			log.Panicf("invalid END value %v", end)
		}
	}
	return v.Pos - 1 + int32(len(v.Ref))
}

// SetEnd sets the end position of a VCF line in the reference by setting the END field.
// If the end position can be calculated from the start position and the length of Ref,
// delete the END field.
func (v *Variant) SetEnd(value int32) {
	if value == v.Pos-1+int32(len(v.Ref)) {
		v.Info.Delete(END)
	} else {
		v.Info.Set(END, int(value))
	}
}

// Pass determines whether the variant passed all filters.
func (v Variant) Pass() bool {
	return len(v.Filter) == 1 && v.Filter[0] == PASS
}
