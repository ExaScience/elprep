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
	"errors"
	"fmt"
	"runtime"
	"sort"

	"github.com/exascience/pargo/pipeline"
)

type (
	// An AlignmentFilter receives an Alignment which it can modify. It
	// returns true if the alignment should be kept, and false if the
	// alignment should be removed.
	AlignmentFilter func(*Alignment) bool

	// A Filter receives a Header and returns an AlignmentFilter or nil.
	Filter func(*Header) AlignmentFilter

	// A PipelineOutput can add nodes to the given pargo
	// pipeline. AddNodes also receives a header that should be added to
	// the output, and a sortingOrder. AddNodes should arrange for the
	// alignments that it receives to be sorted according to that
	// sortingOrder if possible, or report an error if it can't perform
	// such a sort. Any error should be reported to the pipeline by
	// calling p.Err(err) with a non-nil error value.
	PipelineOutput interface {
		AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder SortingOrder)
	}

	// A PipelineInput arranges for a pargo pipeline to be properly
	// initialized, arrange for the pipeline to run the given filters,
	// call output.AddNodes(...), and eventually run the pipeline. If
	// RunPipeline doesn't encounter an error of its own, it should
	// return the error of its pargo pipeline, if any.
	PipelineInput interface {
		RunPipeline(output PipelineOutput, filters []Filter, sortingOrder SortingOrder) error
	}
)

// AlignmentToBytes returns a pargo pipeline.Filter that formats
// slices of Alignment pointers into slices of bytes representing
// these alignments according to the SAM/BAM file format.
func AlignmentToBytes(writer *OutputFile) pipeline.Filter {
	return func(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
		receiver = func(_ int, data interface{}) interface{} {
			alns := data.([]*Alignment)
			records := make([][]byte, 0, len(alns))
			var buf []byte
			var err error
			for _, aln := range alns {
				buf, err = writer.FormatAlignment(aln, buf)
				if err != nil {
					p.SetErr(fmt.Errorf("%v in AlignmentToBytes", err))
				}
				records = append(records, append([]byte(nil), buf...))
				buf = buf[:0]
			}
			return records
		}
		return
	}
}

const (
	minBatchSize = 4096
	maxBatchSize = 262144
)

// BytesToAlignment returns a pargo pipeline.Filter that parses
// slices of bytes representing alignments according to the SAM/BAM file
// format into slices of pointers to freshly allocated Alignment
// values.
func BytesToAlignment(reader *InputFile) pipeline.Filter {
	return BytesToAlignmentFI(reader, false)
}

// BytesToAlignmentFI returns a pargo pipeline.Filter that parses
// slices of bytes representing alignments according to the SAM/BAM file
// format into slices of pointers to freshly allocated Alignment
// values, with an additional option to indicate whether a file index should
// be recorded with each alignment or not.
func BytesToAlignmentFI(reader *InputFile, setFileIndex bool) pipeline.Filter {
	if setFileIndex {
		return func(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
			receiver = func(serial int, data interface{}) interface{} {
				records := data.([][]byte)
				alns := make([]*Alignment, 0, len(records))
				offset := serial * maxBatchSize
				for index, record := range records {
					aln, err := reader.ParseAlignment(record)
					if err != nil {
						p.SetErr(fmt.Errorf("%v, while parsing SAM alignment %v", err, record))
						return alns
					}
					aln.setFileIndex(offset + index)
					alns = append(alns, aln)
				}
				return alns
			}
			return
		}
	}
	return func(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
		receiver = func(_ int, data interface{}) interface{} {
			records := data.([][]byte)
			alns := make([]*Alignment, 0, len(records))
			for _, record := range records {
				aln, err := reader.ParseAlignment(record)
				if err != nil {
					p.SetErr(fmt.Errorf("%v, while parsing SAM alignment %v", err, record))
					return alns
				}
				alns = append(alns, aln)
			}
			return alns
		}
		return
	}
}

// AddNodes implements the PipelineOutput interface for Sam values to
// represent complete SAM/BAM files in memory.
func (sam *Sam) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder SortingOrder) {
	sam.Header = header
	switch sortingOrder {
	case Keep, Unknown:
		p.Add(pipeline.StrictOrd(pipeline.Slice(&sam.Alignments)))
	case Coordinate:
		p.Add(pipeline.Seq(
			pipeline.Slice(&sam.Alignments),
			pipeline.Finalize(func() { By(CoordinateLess).ParallelStableSort(sam.Alignments) }),
		))
	case Queryname:
		p.Add(pipeline.Seq(
			pipeline.Slice(&sam.Alignments),
			pipeline.Finalize(func() { By(QNAMELess).ParallelStableSort(sam.Alignments) }),
		))
	case Unsorted:
		p.Add(pipeline.Seq(pipeline.Slice(&sam.Alignments)))
	default:
		p.SetErr(fmt.Errorf("unknown sorting order %v", sortingOrder))
	}
}

// AddNodes implements the PipelineOutput interface for SAM/BAM OutputFile values.
func (f *OutputFile) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder SortingOrder) {
	if err := f.FormatHeader(header); err != nil {
		p.SetErr(fmt.Errorf("%v, while writing a SAM header to output", err))
		return
	}
	var nodeCons func(...pipeline.Filter) pipeline.Node
	switch sortingOrder {
	case Keep, Unknown:
		nodeCons = pipeline.StrictOrd
	case Coordinate, Queryname:
		p.SetErr(errors.New("sorting on files not supported"))
		return
	case Unsorted:
		nodeCons = pipeline.Seq
	default:
		p.SetErr(fmt.Errorf("unknown sorting order %v", sortingOrder))
		return
	}
	p.Add(
		pipeline.LimitedPar(0, AlignmentToBytes(f)),
		nodeCons(pipeline.Receive(func(_ int, data interface{}) interface{} {
			var err error
			for _, aln := range data.([][]byte) {
				_, err = f.Write(aln)
			}
			if err != nil {
				p.SetErr(fmt.Errorf("%v, while writing SAM alignment strings to output", err))
			}
			return data
		})),
	)
}

// ComposeFilters takes a Header and a slice of Filter functions, and
// successively calls these functions to generate the corresponding
// AlignmentFilter predicates. It then returns a pargo
// pipeline.Receiver that applies these AlignmentFilter predicates on
// the slices of Alignment pointers it receives. ComposeFilters may
// return nil if all AlignmentFilters are nil.
func ComposeFilters(header *Header, hdrFilters []Filter) (receiver pipeline.Receiver) {
	var alnFilters []AlignmentFilter
	for _, f := range hdrFilters {
		if f != nil {
			if alnFilter := f(header); alnFilter != nil {
				alnFilters = append(alnFilters, alnFilter)
			}
		}
	}
	if len(alnFilters) > 0 {
		receiver = func(_ int, data interface{}) interface{} {
			alns := data.([]*Alignment)
			for i, aln := range alns {
				for _, alnFilter := range alnFilters {
					if !alnFilter(aln) {
						n := len(alns)
					jLoop:
						for j := i + 1; j < n; j++ {
							aln := alns[j]
							for _, alnFilter := range alnFilters {
								if !alnFilter(aln) {
									continue jLoop
								}
							}
							alns[i] = aln
							i++
						}
						return alns[0:i]
					}
				}
			}
			return alns
		}
	}
	return
}

// Determine effective sorting order: Some filters may destroy the
// sorting order recorded in the input. If this happens, and
// the requested sorting order is Keep, then we need to effectively
// sort the result according to the original sorting order.
// The reverse is also true: If the requested sorting order is
// Coordinate or Queryname, and the current sorting order already
// fulfills it, then we can just return Keep to avoid any additional
// sorting.
func effectiveSortingOrder(sortingOrder SortingOrder, header *Header, originalSortingOrder SortingOrder) SortingOrder {
	if sortingOrder == Keep {
		sortingOrder = originalSortingOrder
	}
	currentSortingOrder := header.HDSO()
	switch sortingOrder {
	case Coordinate, Queryname:
		if currentSortingOrder == sortingOrder {
			return Keep
		}
		header.SetHDSO(sortingOrder)
	case Unknown, Unsorted:
		if currentSortingOrder != sortingOrder {
			header.SetHDSO(sortingOrder)
		}
	}
	return sortingOrder
}

// NofBatches sets or gets the number of batches that are created from
// this Sam value for the next call of RunPipeline.
//
// NofBatches can be called safely by user programs before RunPipeline
// is called.
//
// If user programs do not call NofBatches, or call them with a value
// < 1, then the pipeline will choose a reasonable default value that
// takes runtime.GOMAXPROCS(0) into account.
func (sam *Sam) NofBatches(n int) {
	sam.nofBatches = n
}

// RunPipeline implements the PipelineInput interface for Sam values
// that represent complete SAM/BAM files in memory.
func (sam *Sam) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder SortingOrder) error {
	header := sam.Header
	alns := sam.Alignments
	sam.Header = NewHeader()
	sam.Alignments = nil
	originalSortingOrder := header.HDSO()
	alnFilter := ComposeFilters(header, hdrFilters)
	sortingOrder = effectiveSortingOrder(sortingOrder, header, originalSortingOrder)
	if out, ok := output.(*Sam); ok && (runtime.GOMAXPROCS(0) <= 3) {
		out.Header = header
		if alnFilter != nil {
			out.Alignments = alnFilter(0, alns).([]*Alignment)
		} else {
			out.Alignments = alns
		}
		switch sortingOrder {
		case Coordinate:
			sort.Slice(out.Alignments, func(i, j int) bool { return CoordinateLess(alns[i], alns[j]) })
		case Queryname:
			sort.Slice(out.Alignments, func(i, j int) bool { return QNAMELess(alns[i], alns[i]) })
		case Keep, Unknown, Unsorted:
			// nothing to do
		default:
			return fmt.Errorf("unknown sorting order %v", sortingOrder)
		}
		return nil
	}
	var p pipeline.Pipeline
	p.Source(alns)
	alns = nil
	if alnFilter != nil {
		p.Add(pipeline.LimitedPar(0, pipeline.Receive(alnFilter)))
	}
	output.AddNodes(&p, header, sortingOrder)
	p.NofBatches(sam.nofBatches)
	sam.nofBatches = 0
	p.Run()
	return p.Err()
}

// RunPipeline implements the PipelineInput interface for SAM/BAM InputFile values.
func (f *InputFile) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder SortingOrder) error {
	return f.RunPipelineFI(output, hdrFilters, sortingOrder, false)
}

// RunPipelineFI implements a variant of the PipelineInput interface for SAM/BAM InputFile
// values, with an additional option to indicate whether a file index should be recorded with
// each alignment or not.
func (f *InputFile) RunPipelineFI(output PipelineOutput, hdrFilters []Filter, sortingOrder SortingOrder, setFileIndex bool) error {
	header, err := f.ParseHeader()
	if err != nil {
		return err
	}
	originalSortingOrder := header.HDSO()
	alnFilter := ComposeFilters(header, hdrFilters)
	sortingOrder = effectiveSortingOrder(sortingOrder, header, originalSortingOrder)
	var p pipeline.Pipeline
	p.Source(f)
	if setFileIndex {
		// needed so that BytesToAlignment can compute file indexes for alignments
		p.SetVariableBatchSize(maxBatchSize, maxBatchSize)
	} else {
		p.SetVariableBatchSize(minBatchSize, maxBatchSize)
	}
	p.Add(pipeline.LimitedPar(0, BytesToAlignmentFI(f, setFileIndex)))
	if alnFilter != nil {
		p.Add(pipeline.LimitedPar(0, pipeline.Receive(alnFilter)))
	}
	output.AddNodes(&p, header, sortingOrder)
	p.Run()
	return p.Err()
}
