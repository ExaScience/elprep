// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2020 imec vzw.

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

package internal

import (
	"encoding/binary"
	"io"
	"log"
	"os"
	"path/filepath"
)

// FilepathAbs is filepath.Abs with panics in place of errors
func FilepathAbs(path string) string {
	result, err := filepath.Abs(path)
	if err != nil {
		log.Panic(err)
	}
	return result
}

// MkdirAll is os.MkdirAll with panics in place of errors
func MkdirAll(path string, perm os.FileMode) {
	if err := os.MkdirAll(path, perm); err != nil {
		log.Panic(err)
	}
}

// FileCreate is os.Create with panics in place of errors
func FileCreate(name string) *os.File {
	f, err := os.Create(name)
	if err != nil {
		log.Panic(err)
	}
	return f
}

// FileCreateOpenIfExists Creates a new file or opens if it already exists for appending
func FileCreateOpenIfExists(name string) *os.File {
	f, err := os.OpenFile(name, os.O_APPEND|os.O_CREATE|os.O_RDWR, 0666)
	if err != nil {
		log.Panic(err)
	}
	return f
}

// FileOpen is os.Open with panics in place of errors
func FileOpen(name string) *os.File {
	f, err := os.Open(name)
	if err != nil {
		log.Panic(err)
	}
	return f
}

// FileOpenIfExists is os.Oppen with panics in place of errors.
// Returns an additional boolean flag to indicate whether the
// actually exists or not.
func FileOpenIfExists(name string) (*os.File, bool) {
	f, err := os.Open(name)
	switch {
	case err == nil:
		return f, true
	case os.IsNotExist(err):
		return nil, false
	default:
		log.Panic(err)
		return nil, false
	}
}

// Close is f.Close with panics in place of errors
func Close(f io.Closer) {
	if err := f.Close(); err != nil {
		log.Panic(err)
	}
}

// Write is f.Write with panics in place of errors
func Write(f io.Writer, b []byte) int {
	n, err := f.Write(b)
	if err != nil {
		log.Panic(err)
	}
	return n
}

// WriteByte is f.WriteByte with panics in place of errors
func WriteByte(f io.ByteWriter, b byte) {
	if err := f.WriteByte(b); err != nil {
		log.Panic(err)
	}
}

// WriteString is f.WriteString with panics in place of errors
func WriteString(f io.StringWriter, s string) int {
	n, err := f.WriteString(s)
	if err != nil {
		log.Panic(err)
	}
	return n
}

// BinaryRead is binary.Read with panics in place of errors
func BinaryRead(r io.Reader, data interface{}) {
	if err := binary.Read(r, binary.LittleEndian, data); err != nil {
		log.Panic(err)
	}
}

// ReadFull is io.ReadFull with panics in place of errors
func ReadFull(r io.Reader, buf []byte) int {
	n, err := io.ReadFull(r, buf)
	if err != nil {
		log.Panic(err)
	}
	return n
}

// Directory returns a slice of filenames. If the given filename
// refers to a directory, return a slice of names of files that are in
// this directory. If the given filename does not refer to a
// directory, return a slice with this filename as the only entry.
// Directory returns a slice of filenames. If the given filename
// refers to a directory, return a slice of names of files that are in
// this directory. If the given filename does not refer to a
// directory, return a slice with this filename as the only entry.
func Directory(file string) (directory string, files []string) {
	info, err := os.Stat(file)
	if err != nil {
		log.Panic(err)
	}
	if !info.IsDir() {
		return filepath.Dir(file), []string{filepath.Base(file)}
	}
	f := FileOpen(file)
	defer Close(f)
	files, err = f.Readdirnames(0)
	if err != nil {
		log.Panic(err)
	}
	return filepath.Clean(file), files
}

// Copy is io.Copy with panics in place of errors
func Copy(dst io.Writer, src io.Reader) int64 {
	written, err := io.Copy(dst, src)
	if err != nil {
		log.Panic(err)
	}
	return written
}
