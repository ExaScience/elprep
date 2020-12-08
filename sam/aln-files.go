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
	"bufio"
	"context"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/exascience/elprep/v5/utils/bgzf"

	"github.com/exascience/elprep/v5/internal"

	"github.com/exascience/pargo/pipeline"
)

type (
	// alignmentReader is a common interface for reading both SAM and BAM files.
	alignmentReader interface {
		ParseHeader() *Header
		SkipHeader()
		ParseAlignment([]byte) *Alignment
		Close()
		pipeline.Source
	}

	// InputFile represents a SAM or BAM file for input.
	InputFile struct {
		reader alignmentReader
	}
)

// Close closes the SAM/BAM input file.
func (f *InputFile) Close() {
	f.reader.Close()
}

// ParseHeader fetches a header from a SAM or BAM file.
func (f *InputFile) ParseHeader() *Header {
	return f.reader.ParseHeader()
}

// SkipHeader skips the header section of a SAM or BAM file.
// This is more efficient than calling ParseHeader and ignoring its result.
func (f *InputFile) SkipHeader() {
	f.reader.SkipHeader()
}

// ParseAlignment parses a block of bytes into an alignment.
// For example in a SAM file, each block of bytes must be
// one line from the alignment section.
func (f *InputFile) ParseAlignment(block []byte) *Alignment {
	return f.reader.ParseAlignment(block)
}

// Err implements the method of the pipeline.Source interface.
func (f *InputFile) Err() error {
	return f.reader.Err()
}

// Prepare implements the method of the pipeline.Source interface.
func (f *InputFile) Prepare(ctx context.Context) int {
	return f.reader.Prepare(ctx)
}

// Fetch implements the method of the pipeline.Source interface.
func (f *InputFile) Fetch(size int) int {
	return f.reader.Fetch(size)
}

// Data implements the method of the pipeline.Source interface.
func (f *InputFile) Data() interface{} {
	return f.reader.Data()
}

type (
	// alignmentWriter is a common interface for writing both SAM and BAM files.
	alignmentWriter interface {
		FormatHeader(hdr *Header)
		FormatAlignment(aln *Alignment, out []byte) []byte
		Close()
		Write(p []byte) int
	}

	// OutputFile represents a SAM or BAM file for output.
	OutputFile struct {
		writer alignmentWriter
	}
)

// Close closes a SAM or BAM output file.
func (f *OutputFile) Close() {
	f.writer.Close()
}

// FormatHeader writes the header to a SAM or BAM file.
func (f *OutputFile) FormatHeader(hdr *Header) {
	f.writer.FormatHeader(hdr)
}

// FormatAlignment formats an alignment into a block of bytes for a SAM or BAM file.
func (f *OutputFile) FormatAlignment(aln *Alignment, out []byte) []byte {
	return f.writer.FormatAlignment(aln, out)
}

// Write can be used to write the blocks of bytes from FormatAlignment
// to the underlying SAM or BAM file.
func (f *OutputFile) Write(p []byte) int {
	return f.writer.Write(p)
}

// SAM file extensions.
const (
	SamExt = ".sam"
	BamExt = ".bam"
)

// Open a SAM or BAM file for input.
//
// Whether the format is SAM or BAM is determined from the content
// of the input, not from any file extensions.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func Open(name string) *InputFile {
	var file io.ReadCloser
	if name == "/dev/stdin" {
		file = os.Stdin
	} else {
		file = internal.FileOpen(name)
	}
	buf := bufio.NewReader(file)
	ok, err := bgzf.IsGzip(buf)
	if err != nil {
		log.Panic(err)
	}
	if !ok {
		return &InputFile{
			reader: &samReader{
				rc:  file,
				buf: buf,
			},
		}
	}
	bgzf, err := bgzf.NewReader(buf)
	if err != nil {
		log.Panic(err)
	}
	return &InputFile{
		reader: &bamReader{
			rc:   file,
			bgzf: bgzf,
		},
	}
}

// OpenIfExists a SAM or BAM file for input, returning false if it doesn't exist.
//
// Whether the format is SAM or BAM is determined from the content
// of the input, not from any file extensions.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func OpenIfExists(name string) (*InputFile, bool) {
	var file io.ReadCloser
	if name == "/dev/stdin" {
		file = os.Stdin
	} else if f, ok := internal.FileOpenIfExists(name); ok {
		file = f
	} else {
		return nil, false
	}
	buf := bufio.NewReader(file)
	ok, err := bgzf.IsGzip(buf)
	if err != nil {
		log.Panic(err)
	}
	if !ok {
		return &InputFile{
			reader: &samReader{
				rc:  file,
				buf: buf,
			},
		}, true
	}
	bgzf, err := bgzf.NewReader(buf)
	if err != nil {
		log.Panic(err)
	}
	return &InputFile{
		reader: &bamReader{
			rc:   file,
			bgzf: bgzf,
		},
	}, true
}

// Create a SAM or BAM file for output.
//
// The format string can be "sam" or "bam". If the format string
// is empty, the output format is determined by looking at the
// filename extension. If the filename extension is not .bam,
// then .sam is always assumed.
//
// The format string will not become part of the resulting filename.
//
// If the name is "/dev/stdout", then the output is written to
// os.Stdout.
func Create(name string, format string) *OutputFile {
	var file io.WriteCloser
	if name == "/dev/stdout" {
		file = os.Stdout
	} else {
		file = internal.FileCreate(name)
	}
	if format == "" {
		format = filepath.Ext(name)
	}
	switch strings.ToLower(format) {
	case "bam", ".bam":
		return &OutputFile{
			writer: &bamWriter{
				wc:   file,
				bgzf: bgzf.NewWriter(file, -1),
			},
		}
	case "cram", ".cram":
		log.Panicf("CRAM format not supported when creating %v", name)
		return nil
	default:
		return &OutputFile{writer: &samWriter{wc: file}}
	}
}
