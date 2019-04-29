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

package fasta

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"strconv"
	"sync"
	"unicode"

	"golang.org/x/sys/unix"
)

// FaiReference represents an entry in an FAI file.
type FaiReference struct {
	Length    int32
	Offset    int64
	LineBases int32
	LineWidth int32
}

func atoi(b []byte) (int32, error) {
	i, err := strconv.ParseInt(string(b), 10, 32)
	return int32(i), err
}

// ParseFai parses an FAI file.
func ParseFai(filename string) (fai map[string]FaiReference, err error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := f.Close(); err == nil {
			err = nerr
		}
	}()

	fai = make(map[string]FaiReference)

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		b := bytes.Split(scanner.Bytes(), []byte("\t"))
		if len(b) != 5 {
			return nil, fmt.Errorf("badly formatted fai file %v - invalid number of entries", filename)
		}

		var r FaiReference
		r.Length, err = atoi(b[1])
		if err == nil {
			r.Offset, err = strconv.ParseInt(string(b[2]), 10, 64)
		}
		if err == nil {
			r.LineBases, err = atoi(b[3])
		}
		if err == nil {
			r.LineWidth, err = atoi(b[4])
		}
		if err != nil {
			return nil, fmt.Errorf("badly formatted fai file %v - invalid integer values", filename)
		}

		fai[string(b[0])] = r
	}

	if err = scanner.Err(); err != nil {
		return nil, fmt.Errorf("%v, while parsing fai file %v", err, filename)
	}

	return fai, nil
}

func contigFromHeader(b []byte) string {
	i := 1
	for ; i < len(b); i++ {
		if c := b[i]; c >= '!' && c <= '~' {
			break
		}
	}
	j := i + 1
	for ; j < len(b); j++ {
		if c := b[j]; c < '!' || c > '~' {
			break
		}
	}
	return string(b[i:j])
}

func initSeq(contig string, fai map[string]FaiReference) []byte {
	if fai != nil {
		if ref, ok := fai[contig]; ok {
			return make([]byte, 0, ref.Length)
		}
	}
	return nil
}

var iupacTable = map[byte]byte{
	'A': 'A', 'a': 'a',
	'C': 'C', 'c': 'c',
	'G': 'G', 'g': 'g',
	'T': 'T', 't': 't',
	'N': 'N', 'n': 'N',
	'R': 'N', 'r': 'N',
	'Y': 'N', 'y': 'N',
	'M': 'N', 'm': 'N',
	'K': 'N', 'k': 'N',
	'W': 'N', 'w': 'N',
	'S': 'N', 's': 'N',
	'B': 'N', 'b': 'N',
	'D': 'N', 'd': 'N',
	'H': 'N', 'h': 'N',
	'V': 'N', 'v': 'N',
}

// ParseFasta sequentially parses a FASTA file.
//
// If fai is given, the sequences can be pre-allocated
// to reduce pressure on the garbage collector.
func ParseFasta(filename string, fai map[string]FaiReference, toUpper, toN bool) (fasta map[string][]byte, err error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := f.Close(); err == nil {
			err = nerr
		}
	}()

	scanner := bufio.NewScanner(f)

	if !scanner.Scan() {
		return nil, fmt.Errorf("empty fasta file %v", filename)
	}
	b := scanner.Bytes()
	for len(b) == 0 {
		if !scanner.Scan() {
			return nil, fmt.Errorf("empty fasta file %v", filename)
		}
		b = scanner.Bytes()
	}
	if b[0] != '>' {
		return nil, fmt.Errorf("invalid fasta file %v - missing first header", filename)
	}

	contig := contigFromHeader(b)
	seq := initSeq(contig, fai)
	fasta = make(map[string][]byte)

scanLoop:
	for scanner.Scan() {
		b := scanner.Bytes()
		if len(b) == 0 {
			if !scanner.Scan() {
				break scanLoop
			}
			b = scanner.Bytes()
			for len(b) == 0 {
				if !scanner.Scan() {
					break scanLoop
				}
				b = scanner.Bytes()
			}
			if b[0] != '>' {
				return nil, fmt.Errorf("invalid fasta file %v - empty line", filename)
			}
		}
		if b[0] == '>' {
			fasta[contig] = seq
			contig = contigFromHeader(b)
			seq = initSeq(contig, fai)
		} else {
			if toUpper {
				for i, c := range b {
					b[i] = byte(unicode.ToUpper(rune(c)))
				}
			}
			if toN {
				for i, c := range b {
					if n, ok := iupacTable[c]; ok {
						b[i] = n
					}
				}
			}
			seq = append(seq, b...)
		}
	}

	fasta[contig] = seq

	if err = scanner.Err(); err != nil {
		return nil, fmt.Errorf("%v, while parsing fasta file %v", err, filename)
	}

	return fasta, nil
}

type offsetTableEntry struct {
	contig string
	offset int
}

// ElfastaMagic is the magic byte sequence that every .elfasta file starts with.
var ElfastaMagic = []byte{0x31, 0xFA, 0x57, 0xA1} // 31FA57A1 => ELFASTA1

// ToElfasta stores fasta data into an elprep-defined mmappable .elfasta file.
func ToElfasta(fasta map[string][]byte, filename string) (err error) {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()
	offset, err := file.Write(ElfastaMagic)
	if err != nil {
		return err
	}
	var offsetTable []offsetTableEntry
	for contig := range fasta {
		n, err := file.WriteString(contig)
		if err != nil {
			return err
		}
		t, err := file.WriteString("\t")
		if err != nil {
			return err
		}
		offset += n + t
		offsetTable = append(offsetTable, offsetTableEntry{contig: contig, offset: offset})
		offset += 2 * binary.MaxVarintLen64
		if _, err := file.Seek(int64(offset), 0); err != nil {
			return err
		}
	}
	n, err := file.WriteString("\n")
	if err != nil {
		return err
	}
	offset += n
	offsetMap := make(map[string]int)
	for contig, ref := range fasta {
		offsetMap[contig] = offset
		n, err := file.Write(ref)
		if err != nil {
			return err
		}
		offset += n
	}
	data, err := unix.Mmap(int(file.Fd()), 0, offset, unix.PROT_WRITE, unix.MAP_SHARED)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := unix.Munmap(data); err == nil {
			err = nerr
		}
	}()
	for _, entry := range offsetTable {
		binary.PutVarint(data[entry.offset:entry.offset+binary.MaxVarintLen64], int64(offsetMap[entry.contig]))
		binary.PutVarint(data[entry.offset+binary.MaxVarintLen64:entry.offset+2*binary.MaxVarintLen64], int64(len(fasta[entry.contig])))
	}
	return nil
}

// MappedFasta represents the contents of a .elfasta file.
type MappedFasta struct {
	fasta map[string][]byte
	data  []byte
	file  *os.File
}

// OpenElfasta opens a .elfasta file.
func OpenElfasta(filename string) (*MappedFasta, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	stat, err := file.Stat()
	if err != nil {
		_ = file.Close()
		return nil, err
	}
	data, err := unix.Mmap(int(file.Fd()), 0, int(stat.Size()), unix.PROT_READ, unix.MAP_SHARED)
	if err != nil {
		_ = file.Close()
		return nil, err
	}
	for i, b := range ElfastaMagic {
		if data[i] != b {
			_ = file.Close()
			return nil, fmt.Errorf("%v is not a .elfasta file - invalid magic byte sequence", filename)
		}
	}
	fasta := make(map[string][]byte)
	index := len(ElfastaMagic)
	for data[index] != '\n' {
		start := index
		for ; data[index] != '\t'; index++ {
		}
		contig := string(data[start:index])
		index++
		offset, n := binary.Varint(data[index : index+binary.MaxVarintLen64])
		if n <= 0 {
			_ = unix.Munmap(data)
			_ = file.Close()
			return nil, fmt.Errorf("bad number of bytes while parsing offset in elfasta file %v", filename)
		}
		size, n := binary.Varint(data[index+binary.MaxVarintLen64 : index+2*binary.MaxVarintLen64])
		if n <= 0 {
			_ = unix.Munmap(data)
			_ = file.Close()
			return nil, fmt.Errorf("bad number of bytes while parsing size in elfasta file %v", filename)
		}
		fasta[contig] = data[int(offset):int(offset+size)]
		index += 2 * binary.MaxVarintLen64
	}
	return &MappedFasta{fasta: fasta, data: data, file: file}, nil
}

// Close closes a .elfasta file.
func (fasta *MappedFasta) Close() (err error) {
	err = unix.Munmap(fasta.data)
	fasta.data = nil
	if nerr := fasta.file.Close(); err == nil {
		err = nerr
	}
	fasta.file = nil
	fasta.fasta = nil
	return err
}

// Seq fetches a sequence for the given contig
// from the .elfasta file.
func (fasta *MappedFasta) Seq(contig string) []byte {
	return fasta.fasta[contig]
}

type (
	// ConcurrentFastaEntry represents an entry
	// in a FASTA file.
	ConcurrentFastaEntry struct {
		wait sync.WaitGroup
		seq  []byte
	}

	// ConcurrentFasta represents a FASTA file
	// that can be parsed concurrently.
	ConcurrentFasta map[string]*ConcurrentFastaEntry
)

// Seq fetches a sequence for the given contig
// from the FASTA file representation. This method
// blocks if the contig has not been parsed yet.
func (fasta ConcurrentFasta) Seq(contig string) []byte {
	entry := fasta[contig]
	entry.wait.Wait()
	return entry.seq
}

// OpenConcurrentFasta prepares the FASTA representation for
// concurrent parsing. The fai parameter must be passed.
func OpenConcurrentFasta(filename, fai string, toUpper, toN bool) (fasta ConcurrentFasta, err error) {
	faiMap, err := ParseFai(fai)
	if err != nil {
		return nil, err
	}
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	fasta = make(map[string]*ConcurrentFastaEntry)
	for contig := range faiMap {
		entry := &ConcurrentFastaEntry{}
		entry.wait.Add(1)
		fasta[contig] = entry
	}

	scanner := bufio.NewScanner(file)

	if !scanner.Scan() {
		return nil, fmt.Errorf("empty fasta file %v", filename)
	}
	b := scanner.Bytes()
	for len(b) == 0 {
		if !scanner.Scan() {
			return nil, fmt.Errorf("empty fasta file %v", filename)
		}
		b = scanner.Bytes()
	}
	if b[0] != '>' {
		return nil, fmt.Errorf("invalid fasta file %v - missing first header", filename)
	}

	go func() {
		defer func() {
			if err := file.Close(); err != nil {
				log.Fatal(err)
			}
		}()

		contig := contigFromHeader(b)
		ref := faiMap[contig]
		seq := make([]byte, 0, ref.Length)

	scanLoop:
		for scanner.Scan() {
			b := scanner.Bytes()
			if len(b) == 0 {
				if !scanner.Scan() {
					break scanLoop
				}
				b = scanner.Bytes()
				for len(b) == 0 {
					if !scanner.Scan() {
						break scanLoop
					}
					b = scanner.Bytes()
				}
				if b[0] != '>' {
					log.Fatal(fmt.Errorf("invalid fasta file %v - empty line", filename))
				}
			}
			if b[0] == '>' {
				entry := fasta[contig]
				entry.seq = seq
				entry.wait.Done()
				contig = contigFromHeader(b)
				ref = faiMap[contig]
				seq = make([]byte, 0, ref.Length)
			} else {
				if toUpper {
					for i, c := range b {
						b[i] = byte(unicode.ToUpper(rune(c)))
					}
				}
				if toN {
					for i, c := range b {
						if n, ok := iupacTable[c]; ok {
							b[i] = n
						}
					}
				}
				seq = append(seq, b...)
			}
		}

		entry := fasta[contig]
		entry.seq = seq
		entry.wait.Done()

		if err = scanner.Err(); err != nil {
			log.Fatal(err)
		}
	}()
	return
}
