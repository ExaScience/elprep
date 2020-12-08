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

package filters

import (
	"math"
	"sort"

	"github.com/exascience/elprep/v5/fasta"
	"github.com/exascience/elprep/v5/utils"
	"github.com/exascience/elprep/v5/vcf"
	"github.com/exascience/pargo/pipeline"
)

const maxGenotypeQual = 99

func (hc *HaplotypeCaller) findGQBand(gq int) (min, max int) {
	if gq > maxGenotypeQual {
		gq = maxGenotypeQual
	}
	index := sort.Search(len(hc.gqBands), func(index int) bool {
		return hc.gqBands[index] > gq
	})
	return hc.gqBands[index-1], hc.gqBands[index]
}

type (
	// a singleVariant is _either_ a fullVariant _or_ some
	// minimal information for a "reduced" variant
	// (never both)
	singleVariant struct {
		fullVariant *vcf.Variant
		location    int32
		dp          int
		ad          [2]int
		pls         [3]int
		gq          int
	}

	variantBlock struct {
		contig   string
		variants []singleVariant
	}

	variantCombiner struct {
		hc                    *HaplotypeCaller
		contig                string
		ref                   []byte
		first                 bool
		vPos, vEnd, nextStart int32
		vRef                  string
		dps                   []int
		minPLs                [3]int
		minGQ, maxGQ          int
	}
)

type variantSlice interface {
	addFullVariant(*vcf.Variant) variantSlice
	addFullVariants([]*vcf.Variant) variantSlice
	addReference(contig string, reference []byte, variant singleVariant) variantSlice
	makeVariantBlock(contig string) interface{}
}

type fullVariants []*vcf.Variant

func (variants fullVariants) addFullVariant(variant *vcf.Variant) variantSlice {
	return append(variants, variant)
}

func (variants fullVariants) addFullVariants(moreVariants []*vcf.Variant) variantSlice {
	return append(variants, moreVariants...)
}

func (variants fullVariants) addReference(contig string, reference []byte, variant singleVariant) variantSlice {
	refByte := fasta.ToUpperAndN(reference[variant.location-1])
	return append(variants, &vcf.Variant{
		Source:         "HC",
		Chrom:          contig,
		Pos:            variant.location,
		Ref:            singleAlleles[refByte],
		Alt:            vNonRefAlt,
		GenotypeFormat: noVariationFormatForNonGvcf,
		GenotypeData: []vcf.Genotype{{
			Phased: false,
			GT:     noVariationGT,
			Data: utils.SmallMap{
				{AD, []interface{}{variant.ad[0], variant.ad[1]}},
				{DP, variant.dp},
				{GQ, variant.gq},
				{PL, []interface{}{variant.pls[0], variant.pls[1], variant.pls[2]}},
			},
		}},
	})
}

func (variants fullVariants) makeVariantBlock(_ string) interface{} {
	return variants
}

type singleVariants []singleVariant

func (variants singleVariants) addFullVariant(variant *vcf.Variant) variantSlice {
	return append(variants, singleVariant{fullVariant: variant})
}

func (variants singleVariants) addFullVariants(moreVariants []*vcf.Variant) variantSlice {
	for _, variant := range moreVariants {
		variants = append(variants, singleVariant{fullVariant: variant})
	}
	return variants
}

func (variants singleVariants) addReference(_ string, _ []byte, variant singleVariant) variantSlice {
	return append(variants, variant)
}

func (variants singleVariants) makeVariantBlock(contig string) interface{} {
	return variantBlock{
		contig:   contig,
		variants: variants,
	}
}

func (hc *HaplotypeCaller) makeVariantSlice() variantSlice {
	if hc.confidenceMode == gvcf {
		return singleVariants(nil)
	}
	return fullVariants(nil)
}

func (hc *HaplotypeCaller) makeVariantCombiner() *variantCombiner {
	return &variantCombiner{hc: hc, first: true, nextStart: -1}
}

func (combiner *variantCombiner) setContig(variants []*vcf.Variant, contig string) []*vcf.Variant {
	if combiner.contig == contig {
		return variants
	}
	if !combiner.first {
		variants = append(variants, combiner.finalizeBlock())
	}
	combiner.contig = contig
	combiner.ref = combiner.hc.reference.Seq(contig)
	combiner.first = true
	combiner.nextStart = -1
	return variants
}

func (combiner *variantCombiner) initBlock(variant singleVariant) {
	combiner.vPos = variant.location
	combiner.vEnd = variant.location
	combiner.vRef = string(fasta.ToUpperAndN(combiner.ref[variant.location-1]))
	combiner.dps = append(combiner.dps[:0], variant.dp)
	combiner.minPLs = variant.pls
	combiner.minGQ, combiner.maxGQ = combiner.hc.findGQBand(variant.gq)
}

var vNonRefAlt = []string{nonRef}

func (combiner *variantCombiner) finalizeBlock() *vcf.Variant {
	sort.Ints(combiner.dps)
	var medianDP int
	if len(combiner.dps)%2 == 0 {
		half := len(combiner.dps) / 2
		medianDP = int(math.Round(float64(combiner.dps[half-1]+combiner.dps[half]) / 2))
	} else {
		half := (len(combiner.dps) + 1) / 2
		medianDP = combiner.dps[half-1]
	}
	return &vcf.Variant{
		Chrom:          combiner.contig,
		Pos:            combiner.vPos,
		Ref:            combiner.vRef,
		Alt:            vNonRefAlt,
		Info:           utils.SmallMap{{vcf.END, int(combiner.vEnd)}},
		GenotypeFormat: noVariationFormatForGvcf,
		GenotypeData: []vcf.Genotype{{
			Phased: false,
			GT:     noVariationGT,
			Data: utils.SmallMap{
				{DP, medianDP},
				{GQ, minInt(computeGQ(combiner.minPLs), maxGenotypeQual)},
				{MIN_DP, combiner.dps[0]},
				{PL, []interface{}{combiner.minPLs[0], combiner.minPLs[1], combiner.minPLs[2]}},
			},
		}},
	}
}

func (combiner *variantCombiner) mergeVariant(variants []*vcf.Variant, variant singleVariant) []*vcf.Variant {
	if variant.fullVariant != nil {
		if !combiner.first {
			variants = append(variants, combiner.finalizeBlock())
			combiner.first = true
		}
		combiner.nextStart = variant.fullVariant.End()
		return append(variants, variant.fullVariant)
	}
	if variant.location <= combiner.nextStart {
		return variants
	}
	if combiner.first {
		combiner.first = false
		combiner.initBlock(variant)
		return variants
	}
	if variant.location == combiner.vEnd+1 && variant.gq >= combiner.minGQ && variant.gq < combiner.maxGQ {
		combiner.vEnd = variant.location
		combiner.dps = append(combiner.dps, variant.dp)
		combiner.minPLs[0] = minInt(combiner.minPLs[0], variant.pls[0])
		combiner.minPLs[1] = minInt(combiner.minPLs[1], variant.pls[1])
		combiner.minPLs[2] = minInt(combiner.minPLs[2], variant.pls[2])
		return variants
	}
	variants = append(variants, combiner.finalizeBlock())
	combiner.initBlock(variant)
	return variants
}

func (combiner *variantCombiner) finalize() (*vcf.Variant, bool) {
	if combiner.first {
		return &vcf.Variant{}, false
	}
	return combiner.finalizeBlock(), true
}

func (hc *HaplotypeCaller) makeVariantCombinerFilter(finalVariant **vcf.Variant, finalVariantOK *bool) pipeline.Filter {
	if hc.confidenceMode == gvcf {
		combiner := hc.makeVariantCombiner()
		return pipeline.ReceiveAndFinalize(
			func(_ int, data interface{}) interface{} {
				block := data.(variantBlock)
				var variants fullVariants
				if block.contig == "" {
					return variants
				}
				variants = combiner.setContig(variants, block.contig)
				for _, singleVariant := range block.variants {
					variants = combiner.mergeVariant(variants, singleVariant)
				}
				return variants
			},
			func() { *finalVariant, *finalVariantOK = combiner.finalize() })
	}
	return pipeline.Identity
}
