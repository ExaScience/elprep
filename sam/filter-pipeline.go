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
	AlignmentFilter func(*Alignment) bool
	Filter          func(*Header) AlignmentFilter

	PipelineOutput interface {
		AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string)
	}

	PipelineInput interface {
		RunPipeline(output PipelineOutput, filters []Filter, sortingOrder string) error
	}
)

func AlignmentToString(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
	receiver = func(_ int, data interface{}) interface{} {
		alns := data.([]*Alignment)
		strings := make([]string, 0, len(alns))
		buf := internal.ReserveByteBuffer()
		defer internal.ReleaseByteBuffer(buf)
		var err error
		for _, aln := range alns {
			buf, err = aln.Format(buf)
			if err != nil {
				p.Err(fmt.Errorf("%v in AlignmentToString", err.Error()))
			}
			strings = append(strings, string(buf))
			buf = buf[:0]
		}
		return strings
	}
	return
}

func StringToAlignment(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
	receiver = func(_ int, data interface{}) interface{} {
		strings := data.([]string)
		alns := make([]*Alignment, 0, len(strings))
		var sc StringScanner
		for _, str := range strings {
			sc.Reset(str)
			aln := sc.ParseAlignment()
			if err := sc.Err(); err != nil {
				p.Err(fmt.Errorf("%v, while parsing SAM alignment %v", err.Error(), str))
				return alns
			}
			alns = append(alns, aln)
		}
		return alns
	}
	return
}

func (output *Sam) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string) {
	output.Header = header
	switch sortingOrder {
	case "keep", "unknown":
		p.Add(pipeline.Ord(pipeline.Slice(&output.Alignments)))
	case "coordinate":
		p.Add(pipeline.Seq(
			pipeline.Slice(&output.Alignments),
			pipeline.Finalize(func() { By(CoordinateLess).ParallelStableSort(output.Alignments) }),
		))
	case "queryname":
		p.Add(pipeline.Seq(
			pipeline.Slice(&output.Alignments),
			pipeline.Finalize(func() {
				By(func(aln1, aln2 *Alignment) bool { return aln1.QNAME < aln2.QNAME }).ParallelStableSort(output.Alignments)
			}),
		))
	case "unsorted":
		p.Add(pipeline.Seq(pipeline.Slice(&output.Alignments)))
	default:
		p.Err(fmt.Errorf("Unknown sorting order %v", sortingOrder))
	}
}

func (output *Writer) AddNodes(p *pipeline.Pipeline, header *Header, sortingOrder string) {
	writer := (*bufio.Writer)(output)
	header.Format(writer)
	var kind pipeline.NodeKind
	switch sortingOrder {
	case "keep", "unknown":
		kind = pipeline.Ordered
	case "coordinate", "queryname":
		p.Err(errors.New("Sorting on files not supported"))
		return
	case "unsorted":
		kind = pipeline.Sequential
	default:
		p.Err(fmt.Errorf("Unknown sorting order %v", sortingOrder))
		return
	}
	p.Add(
		pipeline.Par(AlignmentToString),
		pipeline.NewNode(kind, pipeline.Receive(func(_ int, data interface{}) interface{} {
			for _, aln := range data.([]string) {
				writer.WriteString(aln)
			}
			return data
		})),
	)
}

func ComposeFilters(header *Header, hdrFilters []Filter) (receiver pipeline.Receiver) {
	var alnFilters []AlignmentFilter
	for _, f := range hdrFilters {
		if alnFilter := f(header); alnFilter != nil {
			alnFilters = append(alnFilters, alnFilter)
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

func effectiveSortingOrder(sortingOrder string, header *Header, originalSortingOrder string) string {
	if sortingOrder == "keep" {
		currentSortingOrder := header.HD_SO()
		if currentSortingOrder == originalSortingOrder {
			return "keep"
		}
		header.SetHD_SO(originalSortingOrder)
		return originalSortingOrder
	}
	header.SetHD_SO(sortingOrder)
	return sortingOrder
}

func (input *Sam) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder string) error {
	header := input.Header
	alns := input.Alignments
	*input = Sam{}
	originalSortingOrder := header.HD_SO()
	alnFilter := ComposeFilters(header, hdrFilters)
	sortingOrder = effectiveSortingOrder(sortingOrder, header, originalSortingOrder)
	if out, ok := output.(*Sam); ok && (runtime.GOMAXPROCS(0) <= 3) {
		out.Header = header
		out.Alignments = alns
		if alnFilter != nil {
			alnFilter(0, &out.Alignments)
		}
		switch sortingOrder {
		case "coordinate":
			sort.Slice(out.Alignments, func(i, j int) bool { return CoordinateLess(alns[i], alns[j]) })
		case "queryname":
			sort.Slice(out.Alignments, func(i, j int) bool { return alns[i].QNAME < alns[j].QNAME })
		case "keep", "unknown", "unsorted":
			// nothing to do
		default:
			return fmt.Errorf("Unknown sorting order %v", sortingOrder)
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
	return p.Err(nil)
}

func (input *Reader) RunPipeline(output PipelineOutput, hdrFilters []Filter, sortingOrder string) error {
	reader := (*bufio.Reader)(input)
	header, _, err := ParseHeader(reader)
	if err != nil {
		return err
	}
	originalSortingOrder := header.HD_SO()
	alnFilter := ComposeFilters(header, hdrFilters)
	sortingOrder = effectiveSortingOrder(sortingOrder, header, originalSortingOrder)
	if out, ok := output.(*Writer); ok && (alnFilter == nil) {
		switch sortingOrder {
		case "coordinate", "queryname":
			return errors.New("Sorting on files not supported")
		case "keep", "unknown", "unsorted":
			// nothing to do
		default:
			return fmt.Errorf("Unknown sorting order %v", sortingOrder)
		}
		writer := (*bufio.Writer)(out)
		header.Format(writer)
		_, err := reader.WriteTo(writer)
		return err
	}
	var p pipeline.Pipeline
	p.Source(pipeline.NewScanner(reader))
	p.Add(pipeline.Par(StringToAlignment))
	if alnFilter != nil {
		p.Add(pipeline.Par(pipeline.Receive(alnFilter)))
	}
	output.AddNodes(&p, header, sortingOrder)
	p.Run()
	return p.Err(nil)
}
