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
	"bytes"
	"compress/flate"
	"compress/gzip"
	"context"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"sync"

	"github.com/exascience/pargo/pipeline"
)

// maxBgzfBlockSize defines the maximum block size for BGZF files.
const maxBgzfBlockSize = 65536

var bgzfEOF []byte

func init() {
	bgzfEOF = []byte{
		0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
		0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
		0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
		0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
		0x00, 0x00, 0x00, 0x00,
	}
}

type (
	// bgzfBlock is one block of compressed data in a BGZF file.
	bgzfBlock struct {
		Data  []byte
		Crc32 uint32
		Size  uint32
	}

	// BGZFReader reads in parallel from a BGZF file.
	BGZFReader struct {
		err     error
		r       io.Reader
		gz      *gzip.Reader
		p       pipeline.Pipeline
		w       sync.WaitGroup
		channel chan *bgzfBlock
		ctx     context.Context
		cancel  func()
		data    interface{}
		index   int
		block   *bgzfBlock
	}

	internalBGZFReader BGZFReader
)

var blockPool = sync.Pool{New: func() interface{} {
	return &bgzfBlock{Data: make([]byte, 0, maxBgzfBlockSize)}
}}

func (bgzf *internalBGZFReader) readBgzfBlock() (block *bgzfBlock, err error) {
	var slen int
	for i := 0; i < len(bgzf.gz.Extra); i += 4 + slen {
		if bgzf.gz.Extra[i] == 66 && bgzf.gz.Extra[i+1] == 67 {
			if slen = int(binary.LittleEndian.Uint16(bgzf.gz.Extra[i+2 : i+4])); slen == 2 {
				bsize := int(binary.LittleEndian.Uint16(bgzf.gz.Extra[i+4 : i+6]))
				block = blockPool.Get().(*bgzfBlock)
				block.Data = block.Data[:bsize-len(bgzf.gz.Extra)-19]
				if _, err = io.ReadFull(bgzf.r, block.Data); err != nil {
					return
				}
				var tail [8]byte
				if _, err = io.ReadFull(bgzf.r, tail[:]); err != nil {
					return
				}
				block.Crc32 = binary.LittleEndian.Uint32(tail[0:4])
				block.Size = binary.LittleEndian.Uint32(tail[4:8])
				err = bgzf.gz.Reset(bgzf.r)
				if err == io.EOF {
					if len(block.Data) != 2 || block.Data[0] != 3 || block.Data[1] != 0 || block.Crc32 != 0 || block.Size != 0 {
						err = errors.New("invalid BGZF file: does not end in proper EOF marker")
					}
				}
				return
			}
		}
	}
	err = errors.New("missing BC extra subfield in BGZF header")
	return
}

// Err implements the corresponding method of pipeline.Source
func (bgzf *internalBGZFReader) Err() error {
	if bgzf.err != io.EOF {
		return bgzf.err
	}
	return nil
}

// Prepare implements the corresponding method of pipeline.Source
func (bgzf *internalBGZFReader) Prepare(_ context.Context) (size int) {
	return -1
}

// Fetch implements the corresponding method of pipeline.Source
func (bgzf *internalBGZFReader) Fetch(size int) (fetched int) {
	if bgzf.err != nil {
		return 0
	}
	block, err := bgzf.readBgzfBlock()
	if err != nil {
		bgzf.err = err
		bgzf.data = nil
		return 0
	}
	bgzf.data = block
	return 1
}

// Data implements the corresponding method of pipeline.Source
func (bgzf *internalBGZFReader) Data() interface{} {
	return bgzf.data
}

var flateReaderPool sync.Pool

// NewBGZFReader returns a BGZFReader for the given flate.Reader
func NewBGZFReader(r flate.Reader) (*BGZFReader, error) {
	gz, err := gzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	ctx, cancel := context.WithCancel(context.Background())
	bgzf := &BGZFReader{
		r:       r,
		gz:      gz,
		channel: make(chan *bgzfBlock, 1),
		ctx:     ctx,
		cancel:  cancel,
	}
	bgzf.p.Source((*internalBGZFReader)(bgzf))
	bgzf.p.Add(pipeline.LimitedPar(0, pipeline.Receive(func(_ int, data interface{}) interface{} {
		block := data.(*bgzfBlock)
		blockReader := bytes.NewReader(block.Data)
		var flateReader io.ReadCloser
		if pooled := flateReaderPool.Get(); pooled == nil {
			flateReader = flate.NewReader(blockReader)
		} else {
			flateReader = pooled.(io.ReadCloser)
			if err := flateReader.(flate.Resetter).Reset(blockReader, nil); err != nil {
				flateReader = flate.NewReader(blockReader)
			}
		}
		uncompressed := blockPool.Get().(*bgzfBlock)
		uncompressed.Data = uncompressed.Data[:int(block.Size)]
		if _, err := io.ReadFull(flateReader, uncompressed.Data); err == io.EOF {
			bgzf.p.SetErr(io.ErrUnexpectedEOF)
		} else if err != nil {
			bgzf.p.SetErr(err)
		} else if crc32.ChecksumIEEE(uncompressed.Data) != block.Crc32 {
			bgzf.p.SetErr(errors.New("invalid CRC-32 value for a data block in a BGZF file"))
		}
		if err := flateReader.Close(); err != nil {
			bgzf.p.SetErr(err)
		}
		flateReaderPool.Put(flateReader)
		blockPool.Put(block)
		return uncompressed
	})), pipeline.StrictOrd(pipeline.ReceiveAndFinalize(func(_ int, data interface{}) interface{} {
		select {
		case <-bgzf.ctx.Done():
		case bgzf.channel <- data.(*bgzfBlock):
		}
		return nil
	}, func() {
		close(bgzf.channel)
	})))
	bgzf.w.Add(1)
	go func() {
		defer bgzf.w.Done()
		bgzf.p.Run()
	}()
	return bgzf, nil
}

// Close implements the corresponding method of io.Closer
func (bgzf *BGZFReader) Close() error {
	bgzf.cancel()
	bgzf.w.Wait()
	if err := bgzf.gz.Close(); err != nil {
		return err
	}
	return bgzf.p.Err()
}

func (bgzf *BGZFReader) fetchBlock() (err error) {
	select {
	case <-bgzf.ctx.Done():
		if bgzf.err != nil {
			return bgzf.err
		}
		return bgzf.ctx.Err()
	case b, ok := <-bgzf.channel:
		if !ok {
			return bgzf.err
		}
		bgzf.index = 0
		bgzf.block = b
		return nil
	}
}

// Read implements the corresponding method of io.Reader
func (bgzf *BGZFReader) Read(p []byte) (n int, err error) {
	if bgzf.block == nil {
		if err = bgzf.fetchBlock(); err != nil {
			return
		}
	} else if bgzf.index == len(bgzf.block.Data) {
		blockPool.Put(bgzf.block)
		bgzf.block = nil
		if err = bgzf.fetchBlock(); err != nil {
			return
		}
	}
	n = copy(p, bgzf.block.Data[bgzf.index:])
	bgzf.index += n
	return
}

type (
	bytesBlock struct {
		bytes []byte
	}

	// BGZFWriter writes in parallel to a BGZF file.
	BGZFWriter struct {
		w       io.Writer
		p       pipeline.Pipeline
		wait    sync.WaitGroup
		block   *bytesBlock
		channel chan *bytesBlock
		data    interface{}
	}

	internalBGZFWriter BGZFWriter
)

func (*internalBGZFWriter) Err() error {
	return nil
}

func (writer *internalBGZFWriter) Prepare(_ context.Context) (size int) {
	return -1
}

func (writer *internalBGZFWriter) Fetch(size int) (fetched int) {
	if block, ok := <-writer.channel; ok {
		writer.data = block
		return 1
	}
	writer.data = nil
	return 0
}

func (writer *internalBGZFWriter) Data() interface{} {
	return writer.data
}

var (
	bytesPool = sync.Pool{New: func() interface{} {
		return &bytesBlock{bytes: make([]byte, 0, maxBgzfBlockSize)}
	}}

	flateWriterPool sync.Pool
)

// NewBGZFWriter returns a BGZFWriter for the given io.Writer.
func NewBGZFWriter(w io.Writer) *BGZFWriter {
	bgzf := &BGZFWriter{
		w:       w,
		block:   bytesPool.Get().(*bytesBlock),
		channel: make(chan *bytesBlock, 1),
	}
	bgzf.p.Source((*internalBGZFWriter)(bgzf))
	bgzf.p.Add(pipeline.LimitedPar(0, pipeline.Receive(func(n int, data interface{}) interface{} {
		block := data.(*bytesBlock)
		gzBytes := bytesPool.Get().(*bytesBlock)
		gzBuf := bytes.NewBuffer(gzBytes.bytes)

		gzBuf.Write([]byte{
			0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00,
			0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
			0x42, 0x43, 0x02, 0x00, 0x00, 0x00,
		})

		var flateWriter *flate.Writer
		if pooled := flateWriterPool.Get(); pooled != nil {
			flateWriter = pooled.(*flate.Writer)
			flateWriter.Reset(gzBuf)
		} else {
			var err error
			flateWriter, err = flate.NewWriter(gzBuf, -1)
			if err != nil {
				bgzf.p.SetErr(err)
			}
		}
		if _, err := flateWriter.Write(block.bytes); err != nil {
			bgzf.p.SetErr(err)
		} else if err := flateWriter.Close(); err != nil {
			bgzf.p.SetErr(err)
		}
		gzBytes.bytes = gzBuf.Bytes()
		index := len(gzBytes.bytes)
		gzBytes.bytes = gzBytes.bytes[:index+8]
		binary.LittleEndian.PutUint32(gzBytes.bytes[index:index+4], crc32.ChecksumIEEE(block.bytes))
		binary.LittleEndian.PutUint32(gzBytes.bytes[index+4:index+8], uint32(len(block.bytes)))
		binary.LittleEndian.PutUint16(gzBytes.bytes[16:18], uint16(len(gzBytes.bytes)-1))
		block.bytes = block.bytes[:0]
		bytesPool.Put(block)
		flateWriterPool.Put(flateWriter)
		return gzBytes
	})), pipeline.StrictOrd(pipeline.Receive(func(_ int, data interface{}) interface{} {
		gzBytes := data.(*bytesBlock)
		if _, err := w.Write(gzBytes.bytes); err != nil {
			bgzf.p.SetErr(err)
		}
		gzBytes.bytes = gzBytes.bytes[:0]
		bytesPool.Put(gzBytes)
		return nil
	})))
	bgzf.wait.Add(1)
	go func() {
		defer bgzf.wait.Done()
		bgzf.p.Run()
	}()
	return bgzf
}

func (bgzf *BGZFWriter) sendBlock() (err error) {
	defer func() {
		if x := recover(); x != nil {
			err = errors.New(fmt.Sprint(x))
		}
	}()
	bgzf.channel <- bgzf.block
	return nil
}

// Close closes this BGZFWriter.
func (bgzf *BGZFWriter) Close() error {
	if bgzf.block != nil && len(bgzf.block.bytes) > 0 {
		if err := bgzf.sendBlock(); err != nil {
			return err
		}
	}
	close(bgzf.channel)
	bgzf.wait.Wait()
	if err := bgzf.p.Err(); err != nil {
		return err
	}
	_, err := bgzf.w.Write(bgzfEOF)
	return err
}

// Write implements the corresponding method of io.Writer.
func (bgzf *BGZFWriter) Write(p []byte) (n int, err error) {
	n = len(p)
	for {
		blockIndex := len(bgzf.block.bytes)
		newBlockLength := blockIndex + len(p)
		if newBlockLength >= maxBgzfBlockSize {
			bgzf.block.bytes = bgzf.block.bytes[:maxBgzfBlockSize]
			k := copy(bgzf.block.bytes[blockIndex:], p)
			p = p[k:]
			if err := bgzf.sendBlock(); err != nil {
				return n - len(p), err
			}
			bgzf.block = bytesPool.Get().(*bytesBlock)
		} else {
			bgzf.block.bytes = bgzf.block.bytes[:newBlockLength]
			copy(bgzf.block.bytes[blockIndex:], p)
			return
		}
	}
}
