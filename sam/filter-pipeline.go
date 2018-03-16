package sam

import (
	"bufio"
	"errors"
	"fmt"
	"runtime"
	"sort"

	"github.com/exascience/pargo/pipeline"

	"github.com/exascience/elprep/internal"
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
		AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string)
	}

	// A PipelineInput arranges for a pargo pipeline to be properly
	// initialized, arrange for the pipeline to run the given filters,
	// call output.AddNodes(...), and eventually run the pipeline. If
	// RunPipeline doesn't encounter an error of its own, it should
	// return the error of its pargo pipeline, if any.
	PipelineInput interface {
		RunPipeline(output PipelineOutput, filters []Filter, sortingOrder string) error
	}
)

// AlignmentToString returns a pargo pipeline.Receiver that formats
// slices of Alignment pointers into slices of strings representing
// these alignments according to the SAM file format. See
// http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and
// 1.5.
func AlignmentToString(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
	receiver = func(_ int, data interface{}) interface{} {
		alns := data.([]*Alignment)
		strings := make([]string, 0, len(alns))
		buf := internal.ReserveByteBuffer()
		defer internal.ReleaseByteBuffer(buf)
		var err error
		for _, aln := range alns {
			*buf, err = aln.Format(*buf)
			if err != nil {
				p.SetErr(fmt.Errorf("%v in AlignmentToString", err))
			}
			strings = append(strings, string(*buf))
			*buf = (*buf)[:0]
		}
		return strings
	}
	return
}

// StringToAlignment returns a pargo pipeline.Receiver that parses
// slices of strings representing alignments according to the SAM file
// format into slices of pointers to freshly allocated Alignment
// values. See http://samtools.github.io/hts-specs/SAMv1.pdf -
// Sections 1.4 and 1.5.
func StringToAlignment(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
	receiver = func(_ int, data interface{}) interface{} {
		strings := data.([]string)
		alns := make([]*Alignment, 0, len(strings))
		var sc StringScanner
		for _, str := range strings {
			sc.Reset(str)
			aln := sc.ParseAlignment()
			if err := sc.Err(); err != nil {
				p.SetErr(fmt.Errorf("%v, while parsing SAM alignment %v", err, str))
				return alns
			}
			alns = append(alns, aln)
		}
		return alns
	}
	return
}

// AddNodes implements the PipelineOutput interface for Sam values to
// represent complete SAM files in memory.
func (sam *Sam) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string) {
	sam.Header = header
	switch sortingOrder {
	case Keep, Unknown:
		p.Add(pipeline.Ord(pipeline.Slice(&sam.Alignments)))
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

// AddNodes implements the PipelineOutput interface for Writer values
// to produce output in the SAM file format.
func (output *Writer) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string) {
	writer := (*bufio.Writer)(output)
	if err := header.Format(writer); err != nil {
		p.SetErr(fmt.Errorf("%v, while writing a SAM header to output", err))
		return
	}
	var kind pipeline.NodeKind
	switch sortingOrder {
	case Keep, Unknown:
		kind = pipeline.Ordered
	case Coordinate, Queryname:
		p.SetErr(errors.New("sorting on files not supported"))
		return
	case Unsorted:
		kind = pipeline.Sequential
	default:
		p.SetErr(fmt.Errorf("unknown sorting order %v", sortingOrder))
		return
	}
	p.Add(
		pipeline.Par(AlignmentToString),
		pipeline.NewNode(kind, pipeline.Receive(func(_ int, data interface{}) interface{} {
			var err error
			for _, aln := range data.([]string) {
				_, err = writer.WriteString(aln)
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
func effectiveSortingOrder(sortingOrder string, header *Header, originalSortingOrder string) string {
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

// RunPipeline implements the PipelineInput interface for Sam values
// that represent complete SAM files in memory.
func (sam *Sam) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder string) error {
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
	if alnFilter != nil {
		p.Add(pipeline.Par(pipeline.Receive(alnFilter)))
	}
	output.AddNodes(&p, header, sortingOrder)
	p.Run()
	return p.Err()
}

// RunPipeline implements the PipelineInput interface for Reader
// values that produce input in the SAM file format.
func (input *Reader) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder string) error {
	reader := (*bufio.Reader)(input)
	header, _, err := ParseHeader(reader)
	if err != nil {
		return err
	}
	originalSortingOrder := header.HDSO()
	alnFilter := ComposeFilters(header, hdrFilters)
	sortingOrder = effectiveSortingOrder(sortingOrder, header, originalSortingOrder)
	if out, ok := output.(*Writer); ok && (alnFilter == nil) {
		switch sortingOrder {
		case Coordinate, Queryname:
			return errors.New("sorting on files not supported")
		case Keep, Unknown, Unsorted:
			// nothing to do
		default:
			return fmt.Errorf("unknown sorting order %v", sortingOrder)
		}
		writer := (*bufio.Writer)(out)
		if err := header.Format(writer); err != nil {
			return fmt.Errorf("%v, while writing a SAM header to output", err)
		}
		_, err := reader.WriteTo(writer)
		return fmt.Errorf("%v, while writing SAM alignments to output", err)
	}
	var p pipeline.Pipeline
	p.Source(pipeline.NewScanner(reader))
	p.Add(pipeline.Par(StringToAlignment))
	if alnFilter != nil {
		p.Add(pipeline.Par(pipeline.Receive(alnFilter)))
	}
	output.AddNodes(&p, header, sortingOrder)
	p.Run()
	return p.Err()
}
