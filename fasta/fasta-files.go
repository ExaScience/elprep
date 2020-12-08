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

package fasta

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"log"
	"os"
	"sync"
	"unicode"

	"github.com/exascience/elprep/v5/utils"

	"github.com/exascience/elprep/v5/internal"

	"golang.org/x/sys/unix"
)

// FaiReference represents an entry in an FAI file.
type FaiReference struct {
	Length    int32
	Offset    int64
	LineBases int32
	LineWidth int32
}

// ParseFai parses an FAI file.
func ParseFai(filename string) (fai map[string]FaiReference) {
	f := internal.FileOpen(filename)
	defer internal.Close(f)

	fai = make(map[string]FaiReference)

	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		b := bytes.Split(scanner.Bytes(), []byte("\t"))
		if len(b) != 5 {
			log.Panicf("badly formatted fai file %v - invalid number of entries", filename)
		}

		fai[string(b[0])] = FaiReference{
			Length:    int32(internal.ParseInt(string(b[1]), 10, 32)),
			Offset:    internal.ParseInt(string(b[2]), 10, 64),
			LineBases: int32(internal.ParseInt(string(b[3]), 10, 32)),
			LineWidth: int32(internal.ParseInt(string(b[4]), 10, 32)),
		}
	}

	if err := scanner.Err(); err != nil {
		log.Panic(err)
	}

	return fai
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

// ToN can be used to normalize ambiguity codes in FASTA references.
func ToN(base byte) byte {
	if n, ok := iupacTable[base]; ok {
		return n
	}
	return base
}

var iupacUpperTable = map[byte]byte{
	'A': 'A', 'a': 'A',
	'C': 'C', 'c': 'C',
	'G': 'G', 'g': 'G',
	'T': 'T', 't': 'T',
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

// ToUpperAndN can be used to normalize ambiguity codes in FASTA references,
// and convert all codes to upper case.
func ToUpperAndN(base byte) byte {
	if n, ok := iupacUpperTable[base]; ok {
		return n
	}
	return base
}

// ParseFasta sequentially parses a FASTA file.
//
// If fai is given, the sequences can be pre-allocated
// to reduce pressure on the garbage collector.
// If toUpper is true, the contents are converetd to upper case.
// If toN is true, ambiguity codes are normalized.
func ParseFasta(filename string, fai map[string]FaiReference, toUpper, toN bool) (fasta map[string][]byte) {
	f := internal.FileOpen(filename)
	defer internal.Close(f)

	scanner := bufio.NewScanner(utils.HandleBGZF(bufio.NewReader(f)))

	if !scanner.Scan() {
		log.Panicf("empty fasta file %v", filename)
	}
	b := scanner.Bytes()
	for len(b) == 0 {
		if !scanner.Scan() {
			log.Panicf("empty fasta file %v", filename)
		}
		b = scanner.Bytes()
	}
	if b[0] != '>' {
		log.Panicf("invalid fasta file %v - missing first header", filename)
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
				log.Panicf("invalid fasta file %v - empty line", filename)
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

	if err := scanner.Err(); err != nil {
		log.Panic(err)
	}

	return fasta
}

type offsetTableEntry struct {
	contig string
	offset int
}

// ElfastaMagic is the magic byte sequence that every .elfasta file starts with.
var ElfastaMagic = []byte{0x31, 0xFA, 0x57, 0xA1} // 31FA57A1 => ELFASTA1

// ToElfasta stores fasta data into an elprep-defined mmappable .elfasta file.
func ToElfasta(fasta map[string][]byte, filename string) {
	file := internal.FileCreate(filename)
	defer internal.Close(file)
	offset := internal.Write(file, ElfastaMagic)
	var offsetTable []offsetTableEntry
	for contig := range fasta {
		n := internal.WriteString(file, contig)
		t := internal.WriteString(file, "\t")
		offset += n + t
		offsetTable = append(offsetTable, offsetTableEntry{contig: contig, offset: offset})
		offset += 2 * binary.MaxVarintLen64
		if _, err := file.Seek(int64(offset), 0); err != nil {
			log.Panic(err)
		}
	}
	n := internal.WriteString(file, "\n")
	offset += n
	offsetMap := make(map[string]int)
	for contig, ref := range fasta {
		offsetMap[contig] = offset
		offset += internal.Write(file, ref)
	}
	data, err := unix.Mmap(int(file.Fd()), 0, offset, unix.PROT_WRITE, unix.MAP_SHARED)
	if err != nil {
		log.Panic(err)
	}
	defer func() {
		if err := unix.Munmap(data); err != nil {
			log.Panic(err)
		}
	}()
	for _, entry := range offsetTable {
		binary.PutVarint(data[entry.offset:entry.offset+binary.MaxVarintLen64], int64(offsetMap[entry.contig]))
		binary.PutVarint(data[entry.offset+binary.MaxVarintLen64:entry.offset+2*binary.MaxVarintLen64], int64(len(fasta[entry.contig])))
	}
}

// MappedFasta represents the contents of an .elfasta file.
type MappedFasta struct {
	wait  sync.WaitGroup
	fasta map[string][]byte
	data  []byte
	file  *os.File
}

// OpenElfasta opens a .elfasta file.
func OpenElfasta(filename string) (result *MappedFasta) {
	result = new(MappedFasta)
	result.wait.Add(1)
	go func() {
		defer result.wait.Done()
		file := internal.FileOpen(filename)
		stat, err := file.Stat()
		if err != nil {
			_ = file.Close()
			log.Panic(err)
		}
		data, err := unix.Mmap(int(file.Fd()), 0, int(stat.Size()), unix.PROT_READ, unix.MAP_SHARED)
		if err != nil {
			_ = file.Close()
			log.Panic(err)
		}
		for i, b := range ElfastaMagic {
			if data[i] != b {
				_ = file.Close()
				log.Panicf("%v is not a .elfasta file - invalid magic byte sequence", filename)
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
				log.Panicf("bad number of bytes while parsing offset in elfasta file %v", filename)
			}
			size, n := binary.Varint(data[index+binary.MaxVarintLen64 : index+2*binary.MaxVarintLen64])
			if n <= 0 {
				_ = unix.Munmap(data)
				_ = file.Close()
				log.Panicf("bad number of bytes while parsing size in elfasta file %v", filename)
			}
			fasta[contig] = data[int(offset):int(offset+size)]
			index += 2 * binary.MaxVarintLen64
		}
		result.fasta = fasta
		result.data = data
		result.file = file
	}()
	return result
}

// Close closes the .elfasta file.
func (fasta *MappedFasta) Close() {
	fasta.wait.Wait()
	err := unix.Munmap(fasta.data)
	fasta.data = nil
	if nerr := fasta.file.Close(); err == nil {
		err = nerr
	}
	fasta.file = nil
	fasta.fasta = nil
	if err != nil {
		log.Panic(err)
	}
}

// Seq fetches a sequence for the given contig
// from the .elfasta file.
func (fasta *MappedFasta) Seq(contig string) []byte {
	fasta.wait.Wait()
	return fasta.fasta[contig]
}
