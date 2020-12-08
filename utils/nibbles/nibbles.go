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

package nibbles

import (
	"log"
	"strconv"
)

// Nibbles is a slice-like data structure for storing
// sequences of 4-bit values.
type Nibbles struct {
	info  int
	bytes []byte
}

// Len returns the number of 4-bit values stored in these nibbles.
func (n Nibbles) Len() int {
	return n.info >> 1
}

// Cap returns the capacity of these nibbles.
func (n Nibbles) Cap() int {
	return cap(n.bytes) << 1
}

func (n Nibbles) offset() int {
	return n.info & 1
}

// Make creates nibbles of the given length.
func Make(n int) Nibbles {
	return Nibbles{
		info:  n << 1,
		bytes: make([]byte, (n+1)>>1),
	}
}

// Make2 creates nibbles of the given length and capacity.
func Make2(n, m int) Nibbles {
	return Nibbles{
		info:  n << 1,
		bytes: make([]byte, (n+1)>>1, (m+1)>>1),
	}
}

// ReflectMake creates nibbles of the given length, offset, and raw byte slice.
func ReflectMake(len, offset int, bytes []byte) Nibbles {
	return Nibbles{
		info:  (len << 1) | (offset & 1),
		bytes: bytes,
	}
}

// ReflectValue returns the underlying representation of the nibbles.
func (n Nibbles) ReflectValue() (len, offset int, bytes []byte) {
	return n.Len(), n.offset(), n.bytes
}

// Expand returns a byte slice with the same contents, but where each entry is stored in a byte
func (n Nibbles) Expand() []byte {
	length := n.Len()
	offset := n.offset()
	result := make([]byte, length)
	for k := 0; k < length; k++ {
		index := k + offset
		i := index >> 1
		bit := index & 1
		result[k] = 0xF & (n.bytes[i] >> uint((1^bit)<<2))
	}
	return result
}

// Get returns the nibble at the given index.
func (n Nibbles) Get(index int) byte {
	if index >= n.Len() {
		log.Panic("index out of range")
	}
	index += n.offset()
	i := index >> 1
	bit := index & 1
	return 0xF & (n.bytes[i] >> uint((1^bit)<<2))
}

// Set sets the nibble at the given index.
func (n Nibbles) Set(index int, value byte) {
	if index >= n.Len() {
		log.Panic("index out of range")
	}
	index += n.offset()
	i := index >> 1
	bit := index & 1
	n.bytes[i] = ((0xF << uint(bit<<2)) & n.bytes[i]) | ((0xF & value) << uint((1^bit)<<2))
}

// String returns a string representation of the given nibbles.
func (n Nibbles) String() string {
	if len := n.Len(); len > 0 {
		b := []byte("[")
		b = strconv.AppendInt(b, int64(n.Get(0)), 10)
		for i := 1; i < len; i++ {
			b = append(b, ' ')
			b = strconv.AppendInt(b, int64(n.Get(i)), 10)
		}
		return string(append(b, ']'))
	}
	return "[]"
}

// Slice returns a subsequence of the given nibbles.
func (n Nibbles) Slice(low, high int) Nibbles {
	offset := n.offset()
	return Nibbles{
		info:  ((high - low) << 1) | (offset ^ (low & 1)),
		bytes: n.bytes[(low+offset)>>1 : (high+offset+1)>>1],
	}
}

// Slice3 returns a subsequence of the given nibbles,
// also adapting the capacity.
func (n Nibbles) Slice3(low, high, max int) Nibbles {
	offset := n.offset()
	return Nibbles{
		info:  ((high - low) << 1) | (offset ^ (low & 1)),
		bytes: n.bytes[(low+offset)>>1 : (high+offset+1)>>1 : (max+offset+1)>>1],
	}
}

// Append appends the given nibble.
func (n Nibbles) Append(value byte) Nibbles {
	len := n.Len()
	offset := n.offset()
	index := len + offset
	if index&1 == 1 {
		n := Nibbles{
			info:  ((len + 1) << 1) | offset,
			bytes: n.bytes,
		}
		i := index >> 1
		n.bytes[i] = ((0xF << 4) & n.bytes[i]) | (0xF & value)
		return n
	}
	return Nibbles{
		info:  ((len + 1) << 1) | offset,
		bytes: append(n.bytes, (0xF&value)<<4),
	}
}

// AppendSlice appends the given nibbles.
func (n Nibbles) AppendSlice(m Nibbles) Nibbles {
	mLen := m.Len()
	if mLen == 0 {
		return n
	}
	if mLen == 1 {
		return n.Append(m.Get(0))
	}
	mOffset := m.offset()
	len := n.Len()
	offset := n.offset()
	index := len + offset
	if index&1 == 1 {
		if mOffset == 1 {
			i := index >> 1
			n.bytes[i] = ((0xF << 4) & n.bytes[i]) | (0xF & m.bytes[0])
			return Nibbles{
				info:  ((len + mLen) << 1) | offset,
				bytes: append(n.bytes, m.bytes[1:]...),
			}
		}
	} else {
		if mOffset == 0 {
			return Nibbles{
				info:  ((len + mLen) << 1) | offset,
				bytes: append(n.bytes, m.bytes...),
			}
		}
	}
	for i := 0; i < mLen; i++ {
		n = n.Append(m.Get(i))
	}
	return n
}

// Copy copies the given nibbles.
func (n Nibbles) Copy(m Nibbles) int {
	copyLen := n.Len()
	if mLen := m.Len(); mLen < copyLen {
		copyLen = mLen
	}
	if copyLen == 0 {
		return 0
	}
	offset := n.offset()
	mOffset := m.offset()
	if offset == 0 {
		if mOffset == 0 {
			index := copyLen >> 1
			copy(n.bytes[:index], m.bytes[:index])
			if (copyLen & 1) == 1 {
				n.bytes[index] = (0xF & n.bytes[index]) | ((0xF << 4) & m.bytes[index])
			}
			return copyLen
		}
	} else {
		if mOffset == 1 {
			n.bytes[0] = ((0xF << 4) & n.bytes[0]) | (0xF & m.bytes[0])
			cLen := copyLen + 1
			index := cLen >> 1
			copy(n.bytes[1:index], m.bytes[1:index])
			if (cLen & 1) == 1 {
				n.bytes[index] = (0xF & n.bytes[index]) | ((0xF << 4) & m.bytes[index])
			}
			return copyLen
		}
	}
	for i := 0; i < copyLen; i++ {
		n.Set(i, m.Get(i))
	}
	return copyLen
}
