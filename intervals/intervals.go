// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017-2019 imec vzw.

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

package intervals

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"

	"github.com/exascience/elprep/v4/bed"
	"github.com/exascience/elprep/v4/vcf"
	"github.com/exascience/pargo/parallel"
	"github.com/exascience/pargo/pipeline"
	psort "github.com/exascience/pargo/sort"
)

// Interval is a generic struct with a start and an end position.
type Interval struct {
	Start, End int32
}

// SortByStart sorts a slice of Interval by Start position.
func SortByStart(intervals []Interval) {
	sort.SliceStable(intervals, func(i, j int) bool {
		return intervals[i].Start < intervals[j].Start
	})
}

type stableIntervalSorter []Interval

func (s stableIntervalSorter) SequentialSort(i, j int) {
	SortByStart(s[i:j])
}

func (s stableIntervalSorter) NewTemp() psort.StableSorter {
	return stableIntervalSorter(make([]Interval, len(s)))
}

func (s stableIntervalSorter) Len() int {
	return len(s)
}

func (s stableIntervalSorter) Less(i, j int) bool {
	return s[i].Start < s[j].Start
}

func (s stableIntervalSorter) Assign(source psort.StableSorter) func(i, j, len int) {
	dst, src := s, source.(stableIntervalSorter)
	return func(i, j, len int) {
		copy(dst[i:i+len], src[j:j+len])
	}
}

// ParallelSortByStart sorts a slice of Interval by Start position using
// a parallel stable sort.
func ParallelSortByStart(intervals []Interval) {
	psort.StableSort(stableIntervalSorter(intervals))
}

// Extend makes interval1 larger if it overlaps with interval2,
// by storing max(interval1.End, interval2.End) in interval1.End;
// otherwise, interval1 remains unchanged.
// Returns true if the two intervals overlap, false otherwise.
// interval2.Start >= interval1.Start must be true before
// calling Extend.
func (interval1 *Interval) Extend(interval2 Interval) bool {
	if interval2.Start > interval1.End {
		return false
	}
	if interval2.End > interval1.End {
		interval1.End = interval2.End
	}
	return true
}

// Flatten merges overlapping intervals into larger intervals.
// intervals must be sorted by Start before calling Flatten.
// The resulting slice is sorted by Start, and no two
// intervals in the result overlap with each other.
// The result shares memory with the intervals argument.
func Flatten(intervals []Interval) []Interval {
	for i, n := 0, len(intervals)-1; i < n; i++ {
		if intervals[i].Extend(intervals[i+1]) {
			n++
			for j := i + 1; j < n; j++ {
				if !intervals[i].Extend(intervals[j]) {
					i++
					intervals[i] = intervals[j]
				}
			}
			return intervals[:i+1]
		}
	}
	return intervals
}

const parallelFlattenGrainSize = 0x1000

// ParallelFlatten merges overlapping intervals into larger intervals,
// using a parallel algorithm.
// intervals must be sorted by Start before calling Flatten.
// The resulting slice is sorted by Start, and no two
// intervals in the result overlap with each other.
// The result shares memory with the intervals argument.
func ParallelFlatten(intervals []Interval) []Interval {
	if len(intervals) < parallelFlattenGrainSize {
		return Flatten(intervals)
	}
	half := len(intervals) >> 1
	left, right := intervals[:half], intervals[half:]
	parallel.Do(
		func() { left = ParallelFlatten(left) },
		func() { right = ParallelFlatten(right) },
	)
	for left[len(left)-1].Extend(right[0]) {
		right = right[1:]
	}
	return append(left, right...)
}

// Overlap determines whether the given start/end range overlaps
// with any of the given intervals.
// intervals must be Flattened and sorted by Start.
func Overlap(intervals []Interval, start, end int32) bool {
	for left, right := 0, len(intervals)-1; left <= right; {
		mid := (left + right) / 2
		intervalStart := intervals[mid].Start
		intervalEnd := intervals[mid].End
		if intervalStart > end-1 {
			right = mid - 1
		} else if intervalEnd <= start-1 {
			left = mid + 1
		} else {
			return true
		}
	}
	return false
}

// Intersect returns a slice of all intervals that overlap with the
// given start/end range.
// intervals must be Flattened and sorted by Start.
// The result shares memory with the intervals argument.
func Intersect(intervals []Interval, start, end int32) []Interval {
	n := len(intervals)
	return intervals[sort.Search(n, func(i int) bool {
		return intervals[i].End >= start
	}):sort.Search(n, func(i int) bool {
		return intervals[i].Start > end
	})]
}

// ElsitesHeader is the header line that every .elsites file starts with.
const ElsitesHeader = "# elsites format version 1.0\n"

// ToElsitesFile stores intervals in an elPrep-defined .elsites file.
func ToElsitesFile(intervals map[string][]Interval, filename string) (err error) {
	pathname, err := filepath.Abs(filename)
	if err != nil {
		return err
	}
	output, err := os.Create(pathname)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := output.Close(); nerr != nil {
			if err == nil {
				err = nerr
			}
		}
	}()
	if _, err = output.WriteString(ElsitesHeader); err != nil {
		return err
	}
	for chrom, ivals := range intervals {
		var buf []byte
		for _, ival := range ivals {
			buf = append(buf, chrom...)
			buf = append(buf, '\t')
			buf = strconv.AppendInt(buf, int64(ival.Start), 10)
			buf = append(buf, '\t')
			buf = strconv.AppendInt(buf, int64(ival.End), 10)
			buf = append(buf, '\n')
		}
		if _, err := output.Write(buf); err != nil {
			return err
		}
	}
	return nil
}

// FromElsitesFile loads intervals from an elPrep-defined .elsites file.
func FromElsitesFile(filename string) (intervals map[string][]Interval, err error) {
	pathname, err := filepath.Abs(filename)
	if err != nil {
		return nil, err
	}
	in, err := os.Open(pathname)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := in.Close(); nerr != nil {
			if err == nil {
				err = nerr
			}
		}
	}()
	input := bufio.NewReader(in)
	header, err := input.ReadString('\n')
	if err != nil {
		return nil, err
	}
	if header != ElsitesHeader {
		return nil, fmt.Errorf("%v is not a .elsites file - invalid header", filename)
	}
	var p pipeline.Pipeline
	p.Source(pipeline.NewScanner(input))
	p.Add(pipeline.LimitedPar(0, pipeline.Receive(func(_ int, data interface{}) interface{} {
		strs := data.([]string)
		intervals := make(map[string][]Interval)
		for _, str := range strs {
			i := 0
			for ; i < len(str); i++ {
				if str[i] == '\t' {
					break
				}
			}
			if i == 0 || i == len(str) {
				p.SetErr(fmt.Errorf("invalid sites line %v", str))
				return intervals
			}
			chrom := str[:i]
			j := i + 1
			for i = j; i < len(str); i++ {
				if str[i] == '\t' {
					break
				}
			}
			if i == 0 || i == len(str) {
				p.SetErr(fmt.Errorf("invalid sites line %v", str))
				return intervals
			}
			var interval Interval
			value, err := strconv.ParseInt(str[j:i], 10, 32)
			if err != nil {
				p.SetErr(err)
				return intervals
			}
			interval.Start = int32(value)
			value, err = strconv.ParseInt(str[i+1:], 10, 32)
			if err != nil {
				p.SetErr(err)
				return intervals
			}
			interval.End = int32(value)
			intervals[chrom] = append(intervals[chrom], interval)
		}
		return intervals
	})))
	intervals = make(map[string][]Interval)
	p.Add(pipeline.Ord(pipeline.Receive(func(_ int, data interface{}) interface{} {
		for chrom, ivals := range data.(map[string][]Interval) {
			intervals[chrom] = append(intervals[chrom], ivals...)
		}
		return data
	})))
	p.Run()
	if err = p.Err(); err != nil {
		return nil, err
	}
	return
}

// FromBed returns a slice of intervals that correspond to the BED file entries.
func FromBed(bed *bed.Bed) (intervals map[string][]Interval) {
	intervals = make(map[string][]Interval)
	for _, regions := range bed.RegionMap {
		for _, region := range regions {
			intervals[*region.Chrom] = append(intervals[*region.Chrom], Interval{Start: region.Start, End: region.End})
		}
	}
	return
}

// FromBedFile returns a slice of intervals that correspond to the BED file entries.
func FromBedFile(filename string) (map[string][]Interval, error) {
	pathname, err := filepath.Abs(filename)
	if err != nil {
		return nil, err
	}
	bed, err := bed.ParseBed(pathname)
	if err != nil {
		return nil, err
	}
	return FromBed(bed), nil
}

// FromVcf returns a slice of intervals that correspond to the Vcf file entries.
func FromVcf(vcf *vcf.Vcf) (intervals map[string][]Interval, err error) {
	intervals = make(map[string][]Interval)
	for _, variant := range vcf.Variants {
		start := variant.Start()
		end, err := variant.End()
		if err != nil {
			return nil, err
		}
		intervals[variant.Chrom] = append(intervals[variant.Chrom], Interval{Start: start, End: end})
	}
	return
}

// FromVcfFile returns a slice of intervals that correspond to the Vcf file entries.
func FromVcfFile(filename string) (intervals map[string][]Interval, err error) {
	pathname, err := filepath.Abs(filename)
	if err != nil {
		return nil, err
	}
	input, err := vcf.Open(pathname, false)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := input.Close(); nerr != nil {
			if err == nil {
				intervals = nil
				err = nerr
			}
		}
	}()
	reader := (*bufio.Reader)(input.VcfReader())
	header, _, err := vcf.ParseHeader(reader)
	if err != nil {
		return nil, err
	}
	variantParser, err := header.NewVariantParser()
	if err != nil {
		return nil, err
	}
	variantParser.NSamples = 0 // no need to parse the samples just to retrieve the region information
	var p pipeline.Pipeline
	p.Source(pipeline.NewScanner(reader))
	p.Add(pipeline.LimitedPar(0, func(p *pipeline.Pipeline, _ pipeline.NodeKind, _ *int) (receiver pipeline.Receiver, _ pipeline.Finalizer) {
		receiver = func(_ int, data interface{}) interface{} {
			strings := data.([]string)
			intervals := make(map[string][]Interval)
			var sc vcf.StringScanner
			for _, str := range strings {
				sc.Reset(str)
				variant := sc.ParseVariant(variantParser)
				if err := sc.Err(); err != nil {
					p.SetErr(fmt.Errorf("%v, while parsing VCF variant %v", err, str))
					return intervals
				}
				start := variant.Start()
				end, err := variant.End()
				if err != nil {
					p.SetErr(fmt.Errorf("%v, while parsing VCF variant %v", err, str))
					return intervals
				}
				intervals[variant.Chrom] = append(intervals[variant.Chrom], Interval{Start: start, End: end})
			}
			return intervals
		}
		return
	}))
	intervals = make(map[string][]Interval)
	p.Add(pipeline.Ord(pipeline.Receive(func(_ int, data interface{}) interface{} {
		for chrom, ivals := range data.(map[string][]Interval) {
			intervals[chrom] = append(intervals[chrom], ivals...)
		}
		return data
	})))
	p.Run()
	if err = p.Err(); err != nil {
		return nil, err
	}
	return
}
