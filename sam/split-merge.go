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

package sam

import (
	"context"
	"fmt"
	"log"
	"math"
	"path/filepath"
	"sort"

	"github.com/exascience/pargo/parallel"

	"github.com/exascience/pargo/pipeline"

	"github.com/exascience/elprep/v5/internal"
	"github.com/exascience/elprep/v5/utils"
)

func minInt(i, j int) int {
	if i < j {
		return i
	}
	return j
}

func maxInt(i, j int) int {
	if i > j {
		return i
	}
	return j
}

const mergeGrainSize = 0x3000

// ParallelStableMerge merges two slices of reads in the given order
// by using an efficient parallel merge sort.
func (by By) ParallelStableMerge(alns1, alns2 []*Alignment) (result []*Alignment) {
	result = make([]*Alignment, len(alns1)+len(alns2))
	t := func(index int) *Alignment {
		if index < len(alns1) {
			return alns1[index]
		}
		return alns2[index-len(alns1)]
	}
	assign := func(i, j, sz int) {
		s1 := minInt(j, len(alns1))
		e1 := minInt(j+sz, len(alns1))
		copy(result[i:], alns1[s1:e1])
		s2 := maxInt(j-len(alns1), 0)
		e2 := maxInt(j-len(alns1)+sz, 0)
		copy(result[i+e1-s1:], alns2[s2:e2])
	}
	sMerge := func(p1, r1, p2, r2, p3 int) {
		for {
			if p2 > r2 {
				assign(p3, p1, r1+1-p1)
				return
			}

			q1 := p1
			tp2 := t(p2)
			for (p1 <= r1) && !by(tp2, t(p1)) {
				p1++
			}
			n1 := p1 - q1
			assign(p3, q1, n1)
			p3 += n1

			if p1 > r1 {
				assign(p3, p2, r2+1-p2)
				return
			}

			q2 := p2
			tp1 := t(p1)
			for (p2 <= r2) && by(t(p2), tp1) {
				p2++
			}
			n2 := p2 - q2
			assign(p3, q2, n2)
			p3 += n2
		}
	}
	binarySearchEq := func(x, p, r int) int {
		low, high := p, r+1
		if low > high {
			return low
		}
		tx := t(x)
		for low < high {
			mid := (low + high) / 2
			if !by(t(mid), tx) {
				high = mid
			} else {
				low = mid + 1
			}
		}
		return high
	}
	binarySearchNeq := func(x, p, r int) int {
		low, high := p, r+1
		if low > high {
			return low
		}
		tx := t(x)
		for low < high {
			mid := (low + high) / 2
			if by(tx, t(mid)) {
				high = mid
			} else {
				low = mid + 1
			}
		}
		return high
	}
	var pMerge func(p1, r1, p2, r2, p3 int)
	pMerge = func(p1, r1, p2, r2, p3 int) {
		n1 := r1 - p1 + 1
		n2 := r2 - p2 + 1
		if (n1 + n2) < mergeGrainSize {
			sMerge(p1, r1, p2, r2, p3)
			return
		}
		if n1 > n2 {
			if n1 == 0 {
				return
			}
			q1 := (p1 + r1) / 2
			q2 := binarySearchEq(q1, p2, r2)
			q3 := p3 + (q1 - p1) + (q2 - p2)
			result[q3] = t(q1)
			parallel.Do(
				func() { pMerge(p1, q1-1, p2, q2-1, p3) },
				func() { pMerge(q1+1, r1, q2, r2, q3+1) },
			)
		} else {
			if n2 == 0 {
				return
			}
			q2 := (p2 + r2) / 2
			q1 := binarySearchNeq(q2, p1, r1)
			q3 := p3 + (q1 - p1) + (q2 - p2)
			result[q3] = t(q2)
			parallel.Do(
				func() { pMerge(p1, q1-1, p2, q2-1, p3) },
				func() { pMerge(q1, r1, q2+1, r2, q3+1) },
			)
		}
	}
	pMerge(0, len(alns1)-1, len(alns1), len(alns1)+len(alns2)-1, 0)
	return
}

var sr = utils.Intern("sr")

func formatGroup(groupIndex int) string {
	return fmt.Sprintf("group%05d", groupIndex)
}

func computeContigGroups(SQ []utils.StringMap, contigGroupSize int) (groups []string, contigToGroup map[string]string, groupToContigs map[string][]string) {
	if contigGroupSize <= 0 {
		for _, sn := range SQ {
			if ln := int(SQLN(sn)); ln > contigGroupSize {
				contigGroupSize = ln
			}
		}
		if contigGroupSize <= 0 {
			log.Panic("no valid contig group size")
		}
	}
	groups = []string{"unmapped"}
	contigToGroup = make(map[string]string)
	contigToGroup["*"] = "unmapped"
	groupToContigs = make(map[string][]string)
	groupToContigs["unmapped"] = []string{"*"}
	currentContigGroupIndex := 1
	currentContigGroupSize := 0
	currentContigGroup := formatGroup(currentContigGroupIndex)
	for _, sn := range SQ {
		ln := int(SQLN(sn))
		if currentContigGroupSize > 0 && currentContigGroupSize+ln > contigGroupSize {
			currentContigGroupIndex++
			currentContigGroupSize = 0
			currentContigGroup = formatGroup(currentContigGroupIndex)
		}
		contig := sn["SN"]
		contigToGroup[contig] = currentContigGroup
		groupToContigs[currentContigGroup] = append(groupToContigs[currentContigGroup], contig)
		if groups[len(groups)-1] != currentContigGroup {
			groups = append(groups, currentContigGroup)
		}
		currentContigGroupSize += ln
	}
	return groups, contigToGroup, groupToContigs
}

func (hdr *Header) Contigs() (contigs []string, ok bool) {
	if records := hdr.UserRecords["@cs"]; len(records) > 0 {
		for _, record := range records {
			contigs = append(contigs, record["cs"])
		}
		return contigs, true
	}
	return nil, false
}

// SplitFilePerChromosome splits a SAM file into: a file containing
// all unmapped reads, a file containing all pairs where reads map to
// different chromosomes, and a file per chromosome containing all
// pairs where the reads map to that chromosome. There are no
// requirements on the input file for splitting.
func SplitFilePerChromosome(input, outputPath, outputPrefix, outputExtension string, contigGroupSize int) {
	inputPath, files := internal.Directory(input)
	var header *Header
	var filters []AlignmentFilter
	firstFile := filepath.Join(inputPath, files[0])
	firstIn := Open(firstFile)
	if len(files) > 1 {
		header, filters = MergeInputs(inputPath, files)
		firstIn.SkipHeader() // skip header because it won't be used as it is replaced by the merged header
	} else {
		header = firstIn.ParseHeader()
	}
	groups, contigToGroup, groupToContigs := computeContigGroups(header.SQ, contigGroupSize)
	splitsPath := filepath.Join(outputPath, "splits")
	internal.MkdirAll(splitsPath, 0700)
	header.EnsureUserRecords()["@cs"] = []utils.StringMap{}
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split."})
	groupFiles := make(map[string]*OutputFile)
	for _, group := range groups {
		groupname := filepath.Join(splitsPath, outputPrefix+"-"+group+"."+outputExtension)
		out := Create(groupname, outputExtension)
		defer out.Close()
		var csRecords []utils.StringMap
		for _, contig := range groupToContigs[group] {
			csRecords = append(csRecords, utils.StringMap{"cs": contig})
		}
		header.UserRecords["@cs"] = csRecords
		out.FormatHeader(header)
		groupFiles[group] = out
	}
	delete(header.UserRecords, "@cs")
	spreadname := filepath.Join(outputPath, outputPrefix+"-spread."+outputExtension)
	spreadReads := Create(spreadname, outputExtension)
	defer spreadReads.Close()
	spreadReads.FormatHeader(header)

	var buf []byte

	processFile := func(in *InputFile) {
		var p pipeline.Pipeline
		p.Source(in)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(in)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([]*Alignment) {
					//execute the filters that may update optional rg and pg tags in reads
					for _, filter := range filters {
						filter(aln)
					}
					group := contigToGroup[aln.RNAME]
					out := groupFiles[group]
					if out == nil {
						log.Panic("Attempting to output a read mapped to a region not present in the header: QNAME: ", aln.QNAME, "RNAME: ", aln.RNAME)
					}
					buf = out.FormatAlignment(aln, buf[:0])
					if aln.RNEXT == "=" || aln.RNAME == "*" || contigToGroup[aln.RNEXT] == group {
						out.Write(buf) // untagged
					} else {
						spreadReads.Write(buf) // untagged
						aln.TAGS.Set(sr, int64(1))
						buf = out.FormatAlignment(aln, buf[:0])
						out.Write(buf) // tagged
					}
				}
				return nil
			})),
		)
		internal.RunPipeline(&p)
	}
	processFile(firstIn)
	firstIn.Close()
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		func() {
			in := Open(inFile)
			defer in.Close()
			in.SkipHeader()
			processFile(in)
		}()
	}
}

type (
	groupOut struct {
		channel chan []*Alignment
		block   []*Alignment
		err     chan error
	}

	groupOutSlice struct {
		dictTable     map[string]int32
		currentIndex  int
		currentContig string
		slice         []*groupOut
		data          interface{}
		err           error
	}
)

func (groups *groupOutSlice) Err() error {
	return groups.err
}

func (groups *groupOutSlice) Prepare(ctx context.Context) int {
	return -1
}

func (groups *groupOutSlice) Fetch(size int) (fetched int) {
	for groups.currentIndex < len(groups.slice) {
		group := groups.slice[groups.currentIndex]
		if group.block == nil {
			if block, ok := <-group.channel; ok {
				group.block = block
			} else {
				if err := <-group.err; err != nil {
					if groups.err == nil {
						groups.err = err
					}
					groups.data = nil
					return 0
				}
				groups.slice = append(groups.slice[:groups.currentIndex], groups.slice[groups.currentIndex+1:]...)
				break
			}
		}
		if group.block[0].RNAME == groups.currentContig {
			fetched = len(group.block)
			groups.data = group.block
			group.block = nil
			return
		}
		break
	}
	minRefid := int32(math.MaxInt32)
	minIndex := len(groups.slice)
	for i := 0; i < len(groups.slice); {
		group := groups.slice[i]
		if group.block == nil {
			if block, ok := <-group.channel; ok {
				group.block = block
			} else {
				if err := <-group.err; err != nil {
					if groups.err == nil {
						groups.err = err
					}
					groups.data = nil
					return 0
				}
				groups.slice = append(groups.slice[:i], groups.slice[i+1:]...)
				continue
			}
		}
		if refid, ok := groups.dictTable[group.block[0].RNAME]; !ok {
			log.Panicf("Invalid RNAME %v.", group.block[0].RNAME)
		} else if refid < minRefid {
			minRefid = refid
			minIndex = i
		}
		i++
	}
	if len(groups.slice) == 0 {
		groups.data = nil
		return 0
	}
	group := groups.slice[minIndex]
	groups.currentIndex = minIndex
	groups.currentContig = group.block[0].RNAME
	fetched = len(group.block)
	groups.data = group.block
	group.block = nil
	return
}

func (groups *groupOutSlice) Data() interface{} {
	return groups.data
}

// MergeSortedFilesSplitPerChromosome merges files that were split
// with SplitFilePerChromosome and sorted in coordinate order.
func MergeSortedFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension, outputFormat string, header *Header) {
	dictTable := make(map[string]int32)
	dictTable["*"] = -1
	for index, entry := range header.SQ {
		dictTable[entry["SN"]] = int32(index)
	}

	coordinateLess := func(aln1, aln2 *Alignment) bool {
		refid1, ok := dictTable[aln1.RNAME]
		if !ok {
			refid1 = -1
		}
		refid2, ok := dictTable[aln2.RNAME]
		if !ok {
			refid2 = -1
		}
		switch {
		case refid1 < refid2:
			return refid1 >= 0
		case refid2 < refid1:
			return refid2 < 0
		default:
			return aln1.POS < aln2.POS
		}
	}

	spreadReadsName := filepath.Join(inputPath, inputPrefix+"-spread."+inputExtension)
	spreadReads := Open(spreadReadsName)
	defer spreadReads.Close()
	spreadReads.SkipHeader()
	spreadReads.Prepare(context.Background())

	var spreadRead *Alignment

	fetchSpreadRead := func() {
		spreadRead = nil
		if spreadReads.Fetch(1) == 0 {
			if err := spreadReads.Err(); err != nil {
				log.Panic(err)
			}
			return
		}
		spreadRead = spreadReads.ParseAlignment(spreadReads.Data().([][]byte)[0])
	}

	fetchSpreadRead()

	out := Create(output, outputFormat)
	defer out.Close()

	delete(header.UserRecords, "@cs")
	out.FormatHeader(header)

	groups := &groupOutSlice{dictTable: dictTable}

	for currentContigGroupIndex := 1; ; currentContigGroupIndex++ {
		path := filepath.Join(inputPath, inputPrefix+"-"+formatGroup(currentContigGroupIndex)+"."+inputExtension)
		file, ok := OpenIfExists(path)
		if !ok {
			break
		}

		group := &groupOut{
			channel: make(chan []*Alignment, 1),
			err:     make(chan error),
		}

		groups.slice = append(groups.slice, group)

		go func() {
			defer func() {
				file.Close()
				close(group.err)
				close(group.channel)
			}()

			file.SkipHeader()

			var p pipeline.Pipeline

			p.Source(file)
			p.SetVariableBatchSize(512, 4096)
			p.Add(
				pipeline.LimitedPar(0, BytesToAlignment(file)),
				pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
					alns := data.([]*Alignment)
					for len(alns) > 0 {
						currentContig := alns[0].RNAME
						i := sort.Search(len(alns), func(i int) bool {
							return alns[i].RNAME != currentContig
						})
						// specifying the capacity in the next line is important
						// because insertion of spread reads may otherwise
						// spill reads from this block into the next one
						group.channel <- alns[:i:i]
						alns = alns[i:]
					}
					return nil
				})),
			)
			p.Run()
			if err := p.Err(); err != nil {
				group.err <- err
			}
		}()

	}

	var p pipeline.Pipeline

	p.Source(groups)
	p.Add(
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			alns := data.([]*Alignment)
			if spreadRead == nil {
				return alns
			}
			for i := 0; i < len(alns); i++ {
				if coordinateLess(spreadRead, alns[i]) {
					alns = append(alns[:i+1], alns[i:]...)
					alns[i] = spreadRead
					fetchSpreadRead()
					if spreadRead == nil {
						return alns
					}
				}
			}
			return alns
		})),
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			for _, aln := range data.([][]byte) {
				out.Write(aln)
			}
			return nil
		})),
	)
	internal.RunPipeline(&p)

	var buf []byte

	for spreadRead != nil {
		buf = out.FormatAlignment(spreadRead, buf[:0])
		out.Write(buf)
		fetchSpreadRead()
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile := Open(unmappedName)
	defer unmappedFile.Close()
	unmappedFile.SkipHeader()

	p = pipeline.Pipeline{}

	p.Source(unmappedFile)
	p.SetVariableBatchSize(512, 4096)
	p.Add(
		pipeline.LimitedPar(0, BytesToAlignment(unmappedFile)),
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			for _, aln := range data.([][]byte) {
				out.Write(aln)
			}
			return nil
		})),
	)
	internal.RunPipeline(&p)
}

// MergeUnsortedFilesSplitPerChromosome merges files that were split
// with SplitFilePerChromosome and are unsorted.
func MergeUnsortedFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension, outputFormat string, header *Header) {
	out := Create(output, outputFormat)
	defer out.Close()

	delete(header.UserRecords, "@cs")
	out.FormatHeader(header)

	processFile := func(file *InputFile) {
		defer file.Close()
		file.SkipHeader()

		var p pipeline.Pipeline

		p.Source(file)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(file)),
			pipeline.LimitedPar(0, AlignmentToBytes(out)),
			pipeline.Seq(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([][]byte) {
					out.Write(aln)
				}
				return nil
			})),
		)
		internal.RunPipeline(&p)
	}

	processFile(Open(filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)))
	processFile(Open(filepath.Join(inputPath, inputPrefix+"-spread."+inputExtension)))

	for currentContigGroupIndex := 1; ; currentContigGroupIndex++ {
		if file, ok := OpenIfExists(filepath.Join(inputPath, inputPrefix+"-"+formatGroup(currentContigGroupIndex)+"."+inputExtension)); ok {
			processFile(file)
		} else {
			break
		}
	}
}

func MergeSortedFilesSplitPerChromosomeWithoutSpreadFile(inputPath, output, inputPrefix, inputExtension, outputFormat string, header *Header) {
	out := Create(output, outputFormat)
	defer out.Close()

	delete(header.UserRecords, "@cs")
	out.FormatHeader(header)

	processFile := func(file *InputFile) {
		defer file.Close()
		file.SkipHeader()

		var p pipeline.Pipeline

		p.Source(file)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(file)),
			pipeline.LimitedPar(0, AlignmentToBytes(out)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([][]byte) {
					out.Write(aln)
				}
				return nil
			})),
		)
		internal.RunPipeline(&p)
	}

	processFile(Open(filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)))

	for currentContigGroupIndex := 1; ; currentContigGroupIndex++ {
		if file, ok := OpenIfExists(filepath.Join(inputPath, inputPrefix+"-"+formatGroup(currentContigGroupIndex)+"."+inputExtension)); ok {
			processFile(file)
		} else {
			break
		}
	}
}

// SplitSingleEndFilePerChromosome splits a SAM file containing
// single-end reads into a file for the unmapped reads, and a file per
// chromosome, containing all reads that map to that chromosome. There
// are no requirements on the input file for splitting.
func SplitSingleEndFilePerChromosome(input, outputPath, outputPrefix, outputExtension string, contigGroupSize int) {
	inputPath, files := internal.Directory(input)
	var header *Header
	var filters []AlignmentFilter
	firstFile := filepath.Join(inputPath, files[0])
	firstIn := Open(firstFile)
	if len(files) > 1 {
		header, filters = MergeInputs(inputPath, files)
		firstIn.SkipHeader()
	} else {
		header = firstIn.ParseHeader()
	}
	groups, contigToGroup, _ := computeContigGroups(header.SQ, contigGroupSize)
	groupFiles := make(map[string]*OutputFile)
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split --single-end."})
	for _, group := range groups {
		groupname := filepath.Join(outputPath, outputPrefix+"-"+group+"."+outputExtension)
		out := Create(groupname, outputExtension)
		defer out.Close()
		out.FormatHeader(header)
		groupFiles[group] = out
	}

	var buf []byte

	processFile := func(in *InputFile) {
		var p pipeline.Pipeline
		p.Source(in)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(in)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([]*Alignment) {
					//execute the filters that may update optional rg and pg tags in reads
					for _, filter := range filters {
						filter(aln)
					}
					out := groupFiles[contigToGroup[aln.RNAME]]
					if out == nil {
						log.Panic("Attempting to output a read mapped to a region not present in the header: QNAME: ", aln.QNAME, " RNAME: ", aln.RNAME)
					}
					buf = out.FormatAlignment(aln, buf[:0])
					out.Write(buf) // untagged
				}
				return nil
			})),
		)
		internal.RunPipeline(&p)
	}
	processFile(firstIn)
	firstIn.Close()
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		func() {
			in := Open(inFile)
			defer in.Close()
			in.SkipHeader()
			processFile(in)
		}()
	}
}

// MergeSingleEndFilesSplitPerChromosome merges files containing
// single-end reads that were split with
// SplitSingleEndFilePerChromosome.
func MergeSingleEndFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension, outputFormat string, header *Header) {
	dictTable := make(map[string]int32)
	dictTable["*"] = -1
	for index, entry := range header.SQ {
		dictTable[entry["SN"]] = int32(index)
	}

	out := Create(output, outputFormat)
	defer out.Close()

	delete(header.UserRecords, "@cs")
	out.FormatHeader(header)

	groups := &groupOutSlice{dictTable: dictTable}

	for currentContigGroupIndex := 1; ; currentContigGroupIndex++ {
		path := filepath.Join(inputPath, inputPrefix+"-"+formatGroup(currentContigGroupIndex)+"."+inputExtension)
		file, ok := OpenIfExists(path)
		if !ok {
			break
		}

		group := &groupOut{
			channel: make(chan []*Alignment, 1),
			err:     make(chan error),
		}

		groups.slice = append(groups.slice, group)

		go func() {
			defer func() {
				file.Close()
				close(group.err)
				close(group.channel)
			}()

			file.SkipHeader()

			var p pipeline.Pipeline

			p.Source(file)
			p.SetVariableBatchSize(512, 4096)
			p.Add(
				pipeline.LimitedPar(0, BytesToAlignment(file)),
				pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
					alns := data.([]*Alignment)
					for len(alns) > 0 {
						currentContig := alns[0].RNAME
						i := sort.Search(len(alns), func(i int) bool {
							return alns[i].RNAME != currentContig
						})
						group.channel <- alns[:i]
						alns = alns[i:]
					}
					return nil
				})),
			)
			p.Run()
			if err := p.Err(); err != nil {
				group.err <- err
			}
		}()
	}

	var p pipeline.Pipeline

	p.Source(groups)
	p.Add(
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			for _, aln := range data.([][]byte) {
				out.Write(aln)
			}
			return nil
		})),
	)
	internal.RunPipeline(&p)

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile := Open(unmappedName)
	defer unmappedFile.Close()
	unmappedFile.SkipHeader()

	p = pipeline.Pipeline{}

	p.Source(unmappedFile)
	p.SetVariableBatchSize(512, 4096)
	p.Add(
		pipeline.LimitedPar(0, BytesToAlignment(unmappedFile)),
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			for _, aln := range data.([][]byte) {
				out.Write(aln)
			}
			return nil
		})),
	)
	internal.RunPipeline(&p)
}
