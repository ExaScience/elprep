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
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/exascience/pargo/pipeline"

	"github.com/exascience/elprep/v4/internal"
	"github.com/exascience/elprep/v4/utils"
)

var sr = utils.Intern("sr")

func computeContigGroups(SQ []utils.StringMap, contigGroupSize int) (groups []string, contigMap map[string]string, err error) {
	if contigGroupSize <= 0 {
		for _, sn := range SQ {
			if ln32, err := SQLN(sn); err == nil {
				if ln := int(ln32); ln > contigGroupSize {
					contigGroupSize = ln
				}
			}
		}
		if contigGroupSize <= 0 {
			return nil, nil, errors.New("no valid contig group size")
		}
	}
	contigMap = make(map[string]string)
	contigMap["*"] = "unmapped"
	groups = []string{"unmapped"}
	currentContigGroupIndex := 1
	currentContigGroupSize := 0
	currentContigGroup := fmt.Sprintf("group%v", currentContigGroupIndex)
	for _, sn := range SQ {
		ln32, _ := SQLN(sn)
		ln := int(ln32)
		if currentContigGroupSize > 0 && currentContigGroupSize+ln > contigGroupSize {
			currentContigGroupIndex++
			currentContigGroupSize = 0
			currentContigGroup = fmt.Sprintf("group%v", currentContigGroupIndex)
		}
		contigMap[sn["SN"]] = currentContigGroup
		if groups[len(groups)-1] != currentContigGroup {
			groups = append(groups, currentContigGroup)
		}
		currentContigGroupSize += ln
	}
	return groups, contigMap, nil
}

// SplitFilePerChromosome splits a SAM file into: a file containing
// all unmapped reads, a file containing all pairs where reads map to
// different chromosomes, and a file per chromosome containing all
// pairs where the reads map to that chromosome. There are no
// requirements on the input file for splitting.
func SplitFilePerChromosome(input, outputPath, outputPrefix, outputExtension string, contigGroupSize int) (funcErr error) {
	files, err := internal.Directory(input)
	if err != nil {
		return fmt.Errorf("%v, while attempting to fetch file(s) %v in SplitFilePerChromosome", err, input)
	}
	inputPath := filepath.Dir(input)
	firstFile := filepath.Join(inputPath, files[0])
	firstIn, err := Open(firstFile)
	if err != nil {
		return err
	}
	header, err := firstIn.ParseHeader()
	if err != nil {
		return fmt.Errorf("%v, while parsing header of %v in SplitFilePerChromosome", err, firstFile)
	}
	groups, contigMap, err := computeContigGroups(header.SQ, contigGroupSize)
	if err != nil {
		return fmt.Errorf("%v, while splitting file %v", err, input)
	}
	splitsPath := filepath.Join(outputPath, "splits")
	if err = os.MkdirAll(splitsPath, 0700); err != nil {
		return err
	}
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split."})
	groupFiles := make(map[string]*OutputFile)
	for _, group := range groups {
		groupname := filepath.Join(splitsPath, outputPrefix+"-"+group+"."+outputExtension)
		out, err := Create(groupname)
		if err != nil {
			return err
		}
		defer func() {
			if err := out.Close(); funcErr == nil {
				funcErr = err
			}
		}()
		if err = out.FormatHeader(header); err != nil {
			return fmt.Errorf("%v, while writing header to %v", err, groupname)
		}
		groupFiles[group] = out
	}
	spreadname := filepath.Join(outputPath, outputPrefix+"-spread."+outputExtension)
	spreadReads, err := Create(spreadname)
	if err != nil {
		return err
	}
	defer func() {
		if err := spreadReads.Close(); funcErr == nil {
			funcErr = err
		}
	}()
	if err = spreadReads.FormatHeader(header); err != nil {
		return fmt.Errorf("%v, while writing header to %v", err, spreadname)
	}

	var buf []byte

	processFile := func(in *InputFile) error {
		var p pipeline.Pipeline
		p.Source(in)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(in)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([]*Alignment) {
					group := contigMap[aln.RNAME]
					out := groupFiles[group]
					untagged, err := out.FormatAlignment(aln, buf[:0])
					if err != nil {
						p.SetErr(err)
						return nil
					}
					if aln.RNEXT == "=" || aln.RNAME == "*" || contigMap[aln.RNEXT] == group {
						if _, err := out.Write(untagged); err != nil {
							p.SetErr(err)
							return nil
						}
					} else {
						if _, err := spreadReads.Write(untagged); err != nil {
							p.SetErr(err)
							return nil
						}
						aln.TAGS.Set(sr, int64(1))
						tagged, err := out.FormatAlignment(aln, buf[:0])
						if err != nil {
							p.SetErr(err)
							return nil
						}
						if _, err := out.Write(tagged); err != nil {
							p.SetErr(err)
							return nil
						}
					}
				}
				return nil
			})),
		)
		p.Run()
		return p.Err()
	}
	if err = processFile(firstIn); err != nil {
		return fmt.Errorf("%v, while processing file %v in SplitFilePerChromosome", err, firstFile)
	}
	if err = firstIn.Close(); err != nil {
		return err
	}
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		err = func() (funcErr error) {
			in, err := Open(inFile)
			if err != nil {
				return err
			}
			defer func() {
				if err := in.Close(); funcErr == nil {
					funcErr = err
				}
			}()
			if err = in.SkipHeader(); err != nil {
				return fmt.Errorf("%v, while skipping header of file %v in SplitFilePerChromosome", err, inFile)
			} else if err = processFile(in); err != nil {
				return fmt.Errorf("%v, while processing file %v in SplitFilePerChromosome", err, inFile)
			}
			return nil
		}()
		if err != nil {
			return err
		}
	}
	return nil
}

// MergeSortedFilesSplitPerChromosome merges files that were split
// with SplitFilePerChromosome and sorted in coordinate order.
func MergeSortedFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension string, header *Header, contigGroupSize int) (funcErr error) {

	// Extract the header to identify the files names.  Assume that all
	// files are sorted per cooordinate order, i.e. first sorted on
	// refid entry according to sequence dictionary, then sorted on
	// position entry. There is a file per chromosome in the sequence
	// dictionary. These contain all reads that map to that chromosome.
	// On top of that, there is a file that contains the unmapped (or *)
	// reads and a file that contains the reads that map to different
	// chromosomes. Merge these files in the order of the sequence
	// dictionary. Put the unmapped reads as the last entries. When
	// merging a particular chromosome file into the merged file, make
	// sure that reads that map to different chromosomes are merged in
	// correctly. So while merging a particular chromosome file, pop
	// and compare against reads in the file for reads that map to
	// different chromosomes until the next chromosome is encountered on
	// the refid position. When a file is empty, close it and remove it
	// from the list of files to merge. Loop for identifying and
	// opening the files to merge.

	_, contigMap, err := computeContigGroups(header.SQ, contigGroupSize)
	if err != nil {
		return fmt.Errorf("%v, while merging files from %v", err, inputPath)
	}

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
	spreadReads, err := Open(spreadReadsName)
	if err != nil {
		return err
	}
	defer func() {
		if err := spreadReads.Close(); funcErr == nil {
			funcErr = err
		}
	}()
	err = spreadReads.SkipHeader()
	if err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err, spreadReadsName)
	}

	spreadReads.Prepare(context.Background())

	var spreadRead *Alignment

	fetchSpreadRead := func() (err error) {
		spreadRead = nil
		if spreadReads.Fetch(1) == 0 {
			return spreadReads.Err()
		}
		if spreadRead, err = spreadReads.ParseAlignment(spreadReads.Data().([][]byte)[0]); err != nil {
			spreadRead = nil
		}
		return err
	}

	if err = fetchSpreadRead(); err != nil {
		return fmt.Errorf("%v, while fetching a read from %v", err, spreadReadsName)
	}

	out, err := Create(output)
	if err != nil {
		return err
	}
	defer func() {
		if err := out.Close(); funcErr == nil {
			funcErr = err
		}
	}()

	if err = out.FormatHeader(header); err != nil {
		return fmt.Errorf("%v, while writing header to %v", err, output)
	}

	var buf []byte

	processGroup := func(group, fullInputPath string) (funcErr error) {
		file, err := Open(fullInputPath)
		if err != nil {
			if os.IsNotExist(err) {
				return nil
			}
			return err
		}
		defer func() {
			if err := file.Close(); funcErr == nil {
				funcErr = err
			}
		}()

		if err = file.SkipHeader(); err != nil {
			return fmt.Errorf("%v, while skipping SAM file header in file %v", err, fullInputPath)
		}

		var p pipeline.Pipeline

		p.Source(file)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(file)),
			pipeline.StrictOrd(pipeline.Receive(func(serial int, data interface{}) interface{} {
				if spreadRead == nil || contigMap[spreadRead.RNAME] != group {
					return data
				}
				alns := data.([]*Alignment)
				for i := 0; i < len(alns); i++ {
					if coordinateLess(spreadRead, alns[i]) {
						alns = append(alns[:i+1], alns[i:]...)
						alns[i] = spreadRead
						if err := fetchSpreadRead(); err != nil {
							p.SetErr(err)
							return nil
						}
						if spreadRead == nil || contigMap[spreadRead.RNAME] != group {
							return alns
						}
					}
				}
				return alns
			})),
			pipeline.LimitedPar(0, AlignmentToBytes(out)),
			pipeline.StrictOrd(pipeline.Receive(func(serial int, data interface{}) interface{} {
				var err error
				for _, aln := range data.([][]byte) {
					_, err = out.Write(aln)
				}
				if err != nil {
					p.SetErr(err)
				}
				return nil
			})),
		)
		p.Run()
		if err = p.Err(); err != nil {
			return err
		}
		for spreadRead != nil && contigMap[spreadRead.RNAME] == group {
			if buf, err = out.FormatAlignment(spreadRead, buf[:0]); err != nil {
				return err
			}
			if _, err = out.Write(buf); err != nil {
				return err
			}
			if err = fetchSpreadRead(); err != nil {
				return err
			}
		}
		return nil
	}

	groupsEncountered := make(map[string]bool)

	for _, sn := range header.SQ {
		group := contigMap[sn["SN"]]
		if groupsEncountered[group] {
			continue
		}
		fullInputPath := filepath.Join(inputPath, inputPrefix+"-"+group+"."+inputExtension)
		if err = processGroup(group, fullInputPath); err != nil {
			return fmt.Errorf("%v, while processing %v", err, fullInputPath)
		}
		groupsEncountered[group] = true
	}

	for spreadRead != nil {
		if buf, err = out.FormatAlignment(spreadRead, buf[:0]); err != nil {
			break
		}
		if _, err = out.Write(buf); err != nil {
			break
		}
		if err = fetchSpreadRead(); err != nil {
			break
		}
	}
	if err != nil {
		return fmt.Errorf("%v, while writing alignments to %v", err, output)
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile, err := Open(unmappedName)
	if err != nil {
		return err
	}
	defer func() {
		if err := unmappedFile.Close(); funcErr == nil {
			funcErr = err
		}
	}()
	if err = unmappedFile.SkipHeader(); err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err, unmappedName)
	}

	var p pipeline.Pipeline

	p.Source(unmappedFile)
	p.SetVariableBatchSize(512, 4096)
	p.Add(
		pipeline.LimitedPar(0, BytesToAlignment(unmappedFile)),
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			var err error
			for _, aln := range data.([][]byte) {
				_, err = out.Write(aln)
			}
			if err != nil {
				p.SetErr(err)
			}
			return nil
		})),
	)
	p.Run()
	if err = p.Err(); err != nil {
		return fmt.Errorf("%v, while processing %v", err, unmappedName)
	}
	return nil
}

// MergeUnsortedFilesSplitPerChromosome merges files that were split
// with SplitFilePerChromosome and are unsorted.
func MergeUnsortedFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension string, header *Header, contigGroupSize int) (funcErr error) {

	_, contigMap, err := computeContigGroups(header.SQ, contigGroupSize)
	if err != nil {
		return fmt.Errorf("%v, while merging files from %v", err, inputPath)
	}

	out, err := Create(output)
	if err != nil {
		return err
	}
	defer func() {
		if err := out.Close(); funcErr == nil {
			funcErr = err
		}
	}()

	if err = out.FormatHeader(header); err != nil {
		return fmt.Errorf("%v, while writing header to %v", err, output)
	}

	processFile := func(filename string, missingOK bool) (funcErr error) {
		file, err := Open(filename)
		if err != nil {
			if missingOK && os.IsNotExist(err) {
				return nil
			}
			return err
		}
		defer func() {
			if err := file.Close(); funcErr == nil {
				funcErr = err
			}
		}()
		if err = file.SkipHeader(); err != nil {
			return err
		}

		var p pipeline.Pipeline

		p.Source(file)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(file)),
			pipeline.LimitedPar(0, AlignmentToBytes(out)),
			pipeline.Seq(pipeline.Receive(func(_ int, data interface{}) interface{} {
				var err error
				for _, aln := range data.([][]byte) {
					_, err = out.Write(aln)
				}
				if err != nil {
					p.SetErr(err)
				}
				return nil
			})),
		)
		p.Run()
		return p.Err()
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	if err = processFile(unmappedName, false); err != nil {
		return fmt.Errorf("%v, while processing %v", err, unmappedName)
	}
	spreadName := filepath.Join(inputPath, inputPrefix+"-spread."+inputExtension)
	if err = processFile(spreadName, false); err != nil {
		return fmt.Errorf("%v, while processing %v", err, spreadName)
	}

	groupsEncountered := make(map[string]bool)

	for _, sn := range header.SQ {
		group := contigMap[sn["SN"]]
		if groupsEncountered[group] {
			continue
		}
		groupName := filepath.Join(inputPath, inputPrefix+"-"+group+"."+inputExtension)
		if err = processFile(groupName, true); err != nil {
			return fmt.Errorf("%v, while processing %v", err, groupName)
		}
		groupsEncountered[group] = true
	}

	return nil
}

// SplitSingleEndFilePerChromosome splits a SAM file containing
// single-end reads into a file for the unmapped reads, and a file per
// chromosome, containing all reads that map to that chromosome. There
// are no requirements on the input file for splitting.
func SplitSingleEndFilePerChromosome(input, outputPath, outputPrefix, outputExtension string, contigGroupSize int) (funcErr error) {

	files, err := internal.Directory(input)
	if err != nil {
		return fmt.Errorf("%v, while attempting to fetch file(s) %v in SplitSingleEndFilePerChromosome", err, input)
	}
	inputPath := filepath.Dir(input)
	firstFile := filepath.Join(inputPath, files[0])
	firstIn, err := Open(firstFile)
	if err != nil {
		return err
	}
	header, err := firstIn.ParseHeader()
	if err != nil {
		return fmt.Errorf("%v, while parsing header of %v in SplitSingleEndFilePerChromosome", err, firstFile)
	}
	groups, contigMap, err := computeContigGroups(header.SQ, contigGroupSize)
	if err != nil {
		return fmt.Errorf("%v, while splitting file %v", err, input)
	}
	groupFiles := make(map[string]*OutputFile)
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split --single-end."})
	for _, group := range groups {
		groupname := filepath.Join(outputPath, outputPrefix+"-"+group+"."+outputExtension)
		out, err := Create(groupname)
		if err != nil {
			return err
		}
		defer func() {
			if err := out.Close(); funcErr == nil {
				funcErr = err
			}
		}()
		if err = out.FormatHeader(header); err != nil {
			return fmt.Errorf("%v, while writing header to %v", err, groupname)
		}
		groupFiles[group] = out
	}

	var buf []byte

	processFile := func(in *InputFile) (funcErr error) {
		var p pipeline.Pipeline
		p.Source(in)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(in)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				for _, aln := range data.([]*Alignment) {
					out := groupFiles[contigMap[aln.RNAME]]
					untagged, err := out.FormatAlignment(aln, buf[:0])
					if err != nil {
						p.SetErr(err)
						return nil
					}
					if _, err = out.Write(untagged); err != nil {
						p.SetErr(err)
						return nil
					}
				}
				return nil
			})),
		)
		p.Run()
		return p.Err()
	}
	if err = processFile(firstIn); err != nil {
		return fmt.Errorf("%v, while processing file %v in SplitSingleEndFilePerChromosome", err, firstFile)
	}
	if err = firstIn.Close(); err != nil {
		return err
	}
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		err = func() (funcErr error) {
			in, err := Open(inFile)
			if err != nil {
				return err
			}
			defer func() {
				if err := in.Close(); funcErr == nil {
					funcErr = err
				}
			}()
			if err := in.SkipHeader(); err != nil {
				return fmt.Errorf("%v, while skipping header of file %v in SplitSingleEndFilePerChromosome", err, inFile)
			} else if err = processFile(in); err != nil {
				return fmt.Errorf("%v, while processing file %v in SplitSingleEndFilePerChromosome", err, inFile)
			}
			return nil
		}()
		if err != nil {
			return err
		}
	}
	return nil
}

// MergeSingleEndFilesSplitPerChromosome merges files containing
// single-end reads that were split with
// SplitSingleEndFilePerChromosome.
func MergeSingleEndFilesSplitPerChromosome(inputPath, output, inputPrefix, inputExtension string, header *Header, contigGroupSize int) (funcErr error) {

	_, contigMap, err := computeContigGroups(header.SQ, contigGroupSize)
	if err != nil {
		return fmt.Errorf("%v, while merging files from %v", err, inputPath)
	}

	out, err := Create(output)
	if err != nil {
		return err
	}
	defer func() {
		if err := out.Close(); funcErr == nil {
			funcErr = err
		}
	}()

	if err = out.FormatHeader(header); err != nil {
		return fmt.Errorf("%v, while writing header to %v", err, output)
	}

	processGroup := func(group, fullInputPath string) (funcErr error) {
		file, err := Open(fullInputPath)
		if err != nil {
			return err
		}
		defer func() {
			if err := file.Close(); funcErr == nil {
				funcErr = err
			}
		}()

		if err = file.SkipHeader(); err != nil {
			return fmt.Errorf("%v, while skipping header of file %v", err, fullInputPath)
		}

		var p pipeline.Pipeline
		p.Source(file)
		p.SetVariableBatchSize(512, 4096)
		p.Add(
			pipeline.LimitedPar(0, BytesToAlignment(file)),
			pipeline.LimitedPar(0, AlignmentToBytes(out)),
			pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
				var err error
				for _, aln := range data.([][]byte) {
					_, err = out.Write(aln)
				}
				if err != nil {
					p.SetErr(err)
				}
				return nil
			})),
		)
		p.Run()
		return p.Err()
	}

	groupsEncountered := make(map[string]bool)

	for _, sn := range header.SQ {
		group := contigMap[sn["SN"]]
		if groupsEncountered[group] {
			continue
		}
		fullInputPath := filepath.Join(inputPath, inputPrefix+"-"+group+"."+inputExtension)
		if err = processGroup(group, fullInputPath); err != nil {
			return fmt.Errorf("%v, while processing %v", err, fullInputPath)
		}
		groupsEncountered[group] = true
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile, err := Open(unmappedName)
	if err != nil {
		return err
	}
	defer func() {
		if err := unmappedFile.Close(); funcErr == nil {
			funcErr = err
		}
	}()
	if err = unmappedFile.SkipHeader(); err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err, unmappedName)
	}

	var p pipeline.Pipeline

	p.Source(unmappedFile)
	p.SetVariableBatchSize(512, 4096)
	p.Add(
		pipeline.LimitedPar(0, BytesToAlignment(unmappedFile)),
		pipeline.LimitedPar(0, AlignmentToBytes(out)),
		pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
			var err error
			for _, aln := range data.([][]byte) {
				_, err = out.Write(aln)
			}
			if err != nil {
				p.SetErr(err)
			}
			return nil
		})),
	)
	p.Run()
	if err = p.Err(); err != nil {
		return fmt.Errorf("%v, while processing %v", err, unmappedName)
	}
	return nil
}
