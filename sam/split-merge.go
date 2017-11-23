package sam

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/utils"
)

const (
	// _qname = 1
	// _flag  = 2
	_rname = 3
	_pos   = 4
	// _mapq  = 5
	// _cigar = 6
	_rnext = 7
	// _pnext = 8
	// _tlen  = 9
	// _seq   = 10
	// _qual  = 11
)

type lineScanner struct {
	_err          error
	reader        *bufio.Reader
	nfield, index int
	bytes         []byte
	line          int
	filename      string
}

func newLineScanner(reader *bufio.Reader, filename string, line int) *lineScanner {
	return &lineScanner{
		reader:   reader,
		filename: filename,
		line:     line,
	}
}

func (s *lineScanner) scan() bool {
	if s._err != nil {
		return false
	}
	s.nfield, s.index = 0, 0
	s.bytes, s._err = s.reader.ReadSlice('\n')
	switch {
	case s._err == nil:
		s.bytes = s.bytes[:len(s.bytes)-1]
		fallthrough
	case s._err == io.EOF:
		if len(s.bytes) == 0 {
			return false
		}
		s.line++
		return true
	default:
		s._err = fmt.Errorf("%v, while scanning line %v in file %v", s._err, s.line, s.filename)
		return false
	}
}

func (s *lineScanner) field(n int) string {
	if s.nfield < n {
		start := s.index
		for end := start; end < len(s.bytes); end++ {
			if s.bytes[end] == '\t' {
				s.nfield++
				if s.nfield == n {
					s.index = end + 1
					return string(s.bytes[start:end])
				}
				start = end + 1
			}
		}
	}
	s._err = fmt.Errorf("Invalid index %v, while scanning line %v in file %v", n, s.line, s.filename)
	return ""
}

func (s *lineScanner) parseInt(n int) int64 {
	field := s.field(n)
	value, err := strconv.ParseInt(field, 10, 32)
	if (err != nil) && ((s._err == nil) || (s._err == io.EOF)) {
		s._err = fmt.Errorf("%v, while parsing integer %v in line %v of file %v", err.Error(), field, s.line, s.filename)
	}
	return value
}

func (s *lineScanner) err() error {
	if s._err == io.EOF {
		return nil
	}
	return s._err
}

/*
SplitFilePerChromosome splits a SAM file into: a file containing all
unmapped reads, a file containing all pairs where reads map to
different chromosomes, and a file per chromosome containing all pairs
where the reads map to that chromosome. There are no requirements on
the input file for splitting.
*/
func SplitFilePerChromosome(input, outputPath, outputPrefix, outputExtension, fai, fasta string) (err error) {
	files, err := internal.Directory(input)
	if err != nil {
		return fmt.Errorf("%v, while attempting to fetch file(s) %v in SplitFilePerChromosome", err.Error(), input)
	}
	inputPath := filepath.Dir(input)
	firstFile := filepath.Join(inputPath, files[0])
	firstIn, err := Open(firstFile, false)
	if err != nil {
		return err
	}
	header, lines, err := ParseHeader(firstIn.Reader)
	if err != nil {
		return fmt.Errorf("%v, while parsing header of %v in SplitFilePerChromosome", err.Error(), firstFile)
	}
	splitsPath := filepath.Join(outputPath, "splits")
	chromsEncountered := make(map[string]*OutputFile)
	if err = os.MkdirAll(splitsPath, 0700); err != nil {
		return err
	}
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split."})
	out, err := Create(filepath.Join(splitsPath, outputPrefix+"-unmapped."+outputExtension), fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := out.Close(); err == nil {
			err = nerr
		}
	}()
	header.Format(out.Writer)
	chromsEncountered["*"] = out
	for _, sn := range header.SQ {
		chrom := sn["SN"]
		out, err := Create(filepath.Join(splitsPath, outputPrefix+"-"+chrom+"."+outputExtension), fai, fasta)
		if err != nil {
			return err
		}
		defer func() {
			if nerr := out.Close(); err == nil {
				err = nerr
			}
		}()
		header.Format(out.Writer)
		chromsEncountered[chrom] = out
	}
	spreadReads, err := Create(filepath.Join(outputPath, outputPrefix+"-spread."+outputExtension), fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := spreadReads.Close(); err == nil {
			err = nerr
		}
	}()
	header.Format(spreadReads.Writer)
	processFile := func(in *bufio.Reader, filename string, lines int) error {
		s := newLineScanner(in, filename, lines)
		for s.scan() {
			rname := s.field(_rname)
			rnext := s.field(_rnext)
			out := chromsEncountered[rname]
			if (rnext == "=") || (rname == rnext) || (rname == "*") {
				out.Write(s.bytes)
				out.WriteByte('\n')
			} else {
				spreadReads.Write(s.bytes)
				spreadReads.WriteByte('\n')
				out.Write(s.bytes)
				out.WriteString("\tsr:i:1\n")
			}
		}
		return s.err()
	}
	if err = processFile(firstIn.Reader, firstFile, lines); err != nil {
		return fmt.Errorf("%v, while processing file %v in SplitFilePerChromosome", err.Error(), err)
	}
	if err = firstIn.Close(); err != nil {
		return err
	}
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		err = func() (err error) {
			in, err := Open(inFile, false)
			if err != nil {
				return err
			}
			defer func() {
				if nerr := in.Close(); err == nil {
					err = nerr
				}
			}()
			if lines, err := SkipHeader(in.Reader); err != nil {
				return fmt.Errorf("%v, while skipping header of file %v in SplitFilePerChromosome", err.Error(), inFile)
			} else if err = processFile(in.Reader, inFile, lines); err != nil {
				return fmt.Errorf("%v, while processing file %v in SplitFilePerChromosome", err.Error(), inFile)
			}
			return nil
		}()
		if err != nil {
			return err
		}
	}
	return nil
}

/*
MergeSortedFilesSplitPerChromosome merges files that were split with
SplitFilePerChromosome and sorted in coordinate order.
*/
func MergeSortedFilesSplitPerChromosome(inputPath, output, fai, fasta, inputPrefix, inputExtension string, header *Header) (err error) {

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

	spreadReadsName := filepath.Join(inputPath, inputPrefix+"-spread."+inputExtension)
	spreadReads, err := Open(spreadReadsName, false)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := spreadReads.Close(); err == nil {
			err = nerr
		}
	}()
	spreadReadsLines, err := SkipHeader(spreadReads.Reader)
	if err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err.Error(), spreadReadsName)
	}

	spreadReadsScanner := newLineScanner(spreadReads.Reader, spreadReadsName, spreadReadsLines)
	spreadReadsReady := spreadReadsScanner.scan()

	var spreadReadRname string

	if spreadReadsReady {
		spreadReadRname = spreadReadsScanner.field(_rname)
	} else if err = spreadReadsScanner.err(); err != nil {
		return err
	}

	out, err := Create(output, fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := out.Close(); err == nil {
			err = nerr
		}
	}()

	header.Format(out.Writer)

	finishSpreadReads := func(chrom string) error {
		out.Write(spreadReadsScanner.bytes)
		out.WriteByte('\n')
		for {
			spreadReadsReady = spreadReadsScanner.scan()
			if !spreadReadsReady {
				break
			}
			spreadReadRname = spreadReadsScanner.field(_rname)
			if spreadReadRname != chrom {
				return nil
			}
			out.Write(spreadReadsScanner.bytes)
			out.WriteByte('\n')
		}
		return spreadReadsScanner.err()
	}

	processSN := func(chrom, fullInputPath string) (err error) {
		file, err := Open(fullInputPath, false)
		if err != nil {
			if os.IsNotExist(err) {
				return nil
			}
			return err
		}
		defer func() {
			if nerr := file.Close(); err == nil {
				err = nerr
			}
		}()

		chromosomeLines, err := SkipHeader(file.Reader)
		if err != nil {
			return fmt.Errorf("%v, while skipping SAM file header in file %v", err.Error(), fullInputPath)
		}

		if !spreadReadsReady || (spreadReadRname != chrom) {
			_, err = file.WriteTo(out)
			return err
		}

		chromosomeScanner := newLineScanner(file.Reader, fullInputPath, chromosomeLines)
		chromosomeReady := chromosomeScanner.scan()

		if !chromosomeReady {
			if err = chromosomeScanner.err(); err != nil {
				return err
			}
			return finishSpreadReads(chrom)
		}

		finishChromosomeReads := func() error {
			out.Write(chromosomeScanner.bytes)
			out.WriteByte('\n')
			_, err := file.WriteTo(out)
			return err
		}

		chromosomeReadRname := chromosomeScanner.field(_rname)
		if chromosomeReadRname != chrom {
			return fmt.Errorf("Invalid RNAME %v, expected %v, in file %v", chromosomeReadRname, chrom, fullInputPath)
		}

		spreadReadPos := spreadReadsScanner.parseInt(_pos)
		if err = spreadReadsScanner.err(); err != nil {
			return err
		}

		chromosomeReadPos := chromosomeScanner.parseInt(_pos)
		if err = chromosomeScanner.err(); err != nil {
			return err
		}

		for {
			if spreadReadPos < chromosomeReadPos {
				out.Write(spreadReadsScanner.bytes)
				out.WriteByte('\n')
				spreadReadsReady = spreadReadsScanner.scan()
				if spreadReadsReady {
					spreadReadRname = spreadReadsScanner.field(_rname)
					if spreadReadRname == chrom {
						spreadReadPos = spreadReadsScanner.parseInt(_pos)
						if err = spreadReadsScanner.err(); err != nil {
							return err
						}
					} else {
						return finishChromosomeReads()
					}
				} else if err = spreadReadsScanner.err(); err != nil {
					return err
				} else {
					return finishChromosomeReads()
				}
			} else {
				out.Write(chromosomeScanner.bytes)
				out.WriteByte('\n')
				chromosomeReady = chromosomeScanner.scan()
				if chromosomeReady {
					chromosomeReadRname = chromosomeScanner.field(_rname)
					chromosomeReadPos = chromosomeScanner.parseInt(_pos)
				}
				if err = chromosomeScanner.err(); err != nil {
					return err
				}
				if !chromosomeReady {
					return finishSpreadReads(chrom)
				}
				if chromosomeReadRname != chrom {
					return fmt.Errorf("Invalid RNAME %v, expected %v, in file %v", chromosomeReadRname, chrom, fullInputPath)
				}
			}
		}
	}

	for _, sn := range header.SQ {
		chrom := sn["SN"]
		fullInputPath := filepath.Join(inputPath, inputPrefix+"-"+chrom+"."+inputExtension)
		if err = processSN(chrom, fullInputPath); err != nil {
			return err
		}
	}

	if spreadReadsReady {
		out.Write(spreadReadsScanner.bytes)
		out.WriteByte('\n')
		if _, err = spreadReads.WriteTo(out); err != nil {
			return err
		}
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile, err := Open(unmappedName, false)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := unmappedFile.Close(); err == nil {
			err = nerr
		}
	}()
	if _, err = SkipHeader(unmappedFile.Reader); err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err, unmappedName)
	}
	_, err = unmappedFile.WriteTo(out)
	return err
}

/*
MergeUnsortedFilesSplitPerChromosome merges files that were split with
SplitFilePerChromosome and are unsorted.
*/
func MergeUnsortedFilesSplitPerChromosome(inputPath, output, fai, fasta, inputPrefix, inputExtension string, header *Header) (err error) {
	out, err := Create(output, fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := out.Close(); err == nil {
			err = nerr
		}
	}()

	header.Format(out.Writer)

	processFile := func(filename string, missingOK bool) (err error) {
		file, err := Open(filename, false)
		if err != nil {
			if missingOK && os.IsNotExist(err) {
				return nil
			}
			return err
		}
		defer func() {
			if nerr := file.Close(); err == nil {
				err = nerr
			}
		}()
		if _, err = SkipHeader(file.Reader); err != nil {
			return err
		}
		_, err = file.WriteTo(out)
		return err
	}

	if err = processFile(filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension), false); err != nil {
		return err
	}
	if err = processFile(filepath.Join(inputPath, inputPrefix+"-spread."+inputExtension), false); err != nil {
		return err
	}
	for _, sn := range header.SQ {
		chrom := sn["SN"]
		if err = processFile(filepath.Join(inputPath, inputPrefix+"-"+chrom+"."+inputExtension), true); err != nil {
			return err
		}
	}
	return nil
}

/*
SplitSingleEndFilePerChromosome splits a SAM file containing
single-end reads into a file for the unmapped reads, and a file per
chromosome, containing all reads that map to that chromosome. There
are no requirements on the input file for splitting.
*/
func SplitSingleEndFilePerChromosome(input, outputPath, outputPrefix, outputExtension, fai, fasta string) (err error) {
	files, err := internal.Directory(input)
	if err != nil {
		return fmt.Errorf("%v, while attempting to fetch file(s) %v in SplitSingleEndFilePerChromosome", err.Error(), input)
	}
	inputPath := filepath.Dir(input)
	firstFile := filepath.Join(inputPath, files[0])
	firstIn, err := Open(firstFile, false)
	if err != nil {
		return err
	}
	header, lines, err := ParseHeader(firstIn.Reader)
	if err != nil {
		return fmt.Errorf("%v, while parsing header of %v in SplitSingleEndFilePerChromosome", err.Error(), firstFile)
	}
	chromsEncountered := make(map[string]*OutputFile)
	header.AddUserRecord("@sr", utils.StringMap{"co": "This file was created using elprep split --single-end."})
	out, err := Create(filepath.Join(outputPath, outputPrefix+"-unmapped."+outputExtension), fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := out.Close(); err == nil {
			err = nerr
		}
	}()
	header.Format(out.Writer)
	chromsEncountered["*"] = out
	for _, sn := range header.SQ {
		chrom := sn["SN"]
		out, err := Create(filepath.Join(outputPath, outputPrefix+"-"+chrom+"."+outputExtension), fai, fasta)
		if err != nil {
			return err
		}
		defer func() {
			if nerr := out.Close(); err == nil {
				err = nerr
			}
		}()
		header.Format(out.Writer)
		chromsEncountered[chrom] = out
	}
	processFile := func(in *bufio.Reader, filename string, lines int) error {
		s := newLineScanner(in, filename, lines)
		for s.scan() {
			rname := s.field(_rname)
			out := chromsEncountered[rname]
			out.Write(s.bytes)
			out.WriteByte('\n')
		}
		return s.err()
	}
	if err = processFile(firstIn.Reader, firstFile, lines); err != nil {
		return fmt.Errorf("%v, while processing file %v in SplitSingleEndFilePerChromosome", err.Error(), err)
	}
	if err = firstIn.Close(); err != nil {
		return err
	}
	for _, name := range files[1:] {
		inFile := filepath.Join(inputPath, name)
		err = func() (err error) {
			in, err := Open(inFile, false)
			if err != nil {
				return err
			}
			defer func() {
				if nerr := in.Close(); err == nil {
					err = nerr
				}
			}()
			if lines, err := SkipHeader(in.Reader); err != nil {
				return fmt.Errorf("%v, while skipping header of file %v in SplitSingleEndFilePerChromosome", err.Error(), inFile)
			} else if err = processFile(in.Reader, inFile, lines); err != nil {
				return fmt.Errorf("%v, while processing file %v in SplitSingleEndFilePerChromosome", err.Error(), inFile)
			}
			return nil
		}()
		if err != nil {
			return err
		}
	}
	return nil
}

/*
MergeSingleEndFilesSplitPerChromosome merges files containing
single-end reads that were split with SplitSingleEndFilePerChromosome.
*/
func MergeSingleEndFilesSplitPerChromosome(inputPath, output, fai, fasta, inputPrefix, inputExtension string, header *Header) (err error) {

	out, err := Create(output, fai, fasta)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := out.Close(); err == nil {
			err = nerr
		}
	}()

	header.Format(out.Writer)

	for _, sn := range header.SQ {
		chrom := sn["SN"]
		chromName := filepath.Join(inputPath, inputPrefix+"-"+chrom+"."+inputExtension)
		chromFile, err := Open(chromName, false)
		if err != nil {
			return err
		}
		defer func() {
			if nerr := chromFile.Close(); err == nil {
				err = nerr
			}
		}()
		if _, err = SkipHeader(chromFile.Reader); err != nil {
			return fmt.Errorf("%v, while skipping header of file %v", err, chromName)
		}
		_, err = chromFile.WriteTo(out)
		if err != nil {
			return err
		}
	}

	unmappedName := filepath.Join(inputPath, inputPrefix+"-unmapped."+inputExtension)
	unmappedFile, err := Open(unmappedName, false)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := unmappedFile.Close(); err == nil {
			err = nerr
		}
	}()
	if _, err = SkipHeader(unmappedFile.Reader); err != nil {
		return fmt.Errorf("%v, while skipping header of file %v", err, unmappedName)
	}
	_, err = unmappedFile.WriteTo(out)
	return err
}
