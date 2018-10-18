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
	"bufio"
	"context"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/exascience/pargo/pipeline"
)

type (
	// alignmentReader is a common interface for reading both SAM and BAM files.
	alignmentReader interface {
		ParseHeader() (*Header, error)
		SkipHeader() error
		ParseAlignment([]byte) (*Alignment, error)
		pipeline.Source
		io.Closer
	}

	// InputFile represents a SAM or BAM file for input.
	InputFile struct {
		reader alignmentReader
	}
)

// Close closes the SAM/BAM input file.
func (f *InputFile) Close() error {
	return f.reader.Close()
}

// ParseHeader fetches a header from a SAM or BAM file.
func (f *InputFile) ParseHeader() (*Header, error) {
	return f.reader.ParseHeader()
}

// SkipHeader skips the header section of a SAM or BAM file.
// This is more efficient than calling ParseHeader and ignoring its result.
func (f *InputFile) SkipHeader() error {
	return f.reader.SkipHeader()
}

// ParseAlignment parses a block of bytes into an alignment.
// For example in a SAM file, each block of bytes must be
// one line from the alignment section.
func (f *InputFile) ParseAlignment(block []byte) (*Alignment, error) {
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
		FormatHeader(hdr *Header) error
		FormatAlignment(aln *Alignment, out []byte) ([]byte, error)
		io.WriteCloser
	}

	// OutputFile represents a SAM or BAM file for output.
	OutputFile struct {
		writer alignmentWriter
	}
)

// Close closes a SAM or BAM output file.
func (f *OutputFile) Close() error {
	return f.writer.Close()
}

// FormatHeader writes the header to a SAM or BAM file.
func (f *OutputFile) FormatHeader(hdr *Header) error {
	return f.writer.FormatHeader(hdr)
}

// FormatAlignment formats an alignment into a block of bytes for a SAM or BAM file.
func (f *OutputFile) FormatAlignment(aln *Alignment, out []byte) ([]byte, error) {
	return f.writer.FormatAlignment(aln, out)
}

// Write can be used to write the blocks of bytes from FormatAlignment
// to the underlying SAM or BAM file.
func (f *OutputFile) Write(p []byte) (int, error) {
	return f.writer.Write(p)
}

// SAM file extensions.
const (
	SamExt  = ".sam"
	BamExt  = ".bam"
	cramExt = ".cram"
)

// Open a SAM or BAM file for input.
//
// If the filename extension is not .bam, then .sam is always
// assumed.
//
// If the name is "/dev/stdin", then the input is read from os.Stdin
func Open(name string) (*InputFile, error) {
	switch filepath.Ext(name) {
	case BamExt:
		file, err := os.Open(name)
		if err != nil {
			return nil, err
		}
		bgzf, err := NewBGZFReader(bufio.NewReader(file))
		if err != nil {
			return nil, err
		}
		return &InputFile{
			reader: &bamReader{
				rc:   file,
				bgzf: bgzf,
			},
		}, nil
	case cramExt:
		return nil, fmt.Errorf("CRAM format not supported when opening %v", name)
	default:
		if name == "/dev/stdin" {
			return &InputFile{
				reader: &samReader{
					rc:  os.Stdin,
					buf: bufio.NewReader(os.Stdin),
				},
			}, nil
		}
		file, err := os.Open(name)
		if err != nil {
			return nil, err
		}
		return &InputFile{
			reader: &samReader{
				rc:  file,
				buf: bufio.NewReader(file),
			},
		}, nil
	}
}

// Create a SAM or BAM file for output.
//
// If the filename extension is not .bam, then .sam is always
// assumed.
//
// If the name is "/dev/stdout", then the output is written to
// os.Stdout.
func Create(name string) (*OutputFile, error) {
	switch filepath.Ext(name) {
	case BamExt:
		file, err := os.Create(name)
		if err != nil {
			return nil, err
		}
		return &OutputFile{
			writer: &bamWriter{
				wc:   file,
				bgzf: NewBGZFWriter(file),
			},
		}, nil
	case cramExt:
		return nil, fmt.Errorf("CRAM format not supported when opening %v", name)
	default:
		if name == "/dev/stdout" {
			return &OutputFile{writer: &samWriter{wc: os.Stdout}}, nil
		}
		file, err := os.Create(name)
		if err != nil {
			return nil, err
		}
		return &OutputFile{writer: &samWriter{wc: file}}, nil
	}
}
