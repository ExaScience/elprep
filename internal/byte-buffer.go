package internal

import "sync"

var bufPool = sync.Pool{New: func() interface{} {
	b := make([]byte, 0)
	return &b
}}

// ReserveByteBuffer uses a sync.Pool to either reuse or make a slice
// of bytes of length 0, but of capacity potentially larger than 0.
//
// Use ReleaseByteBuffer to return slices of bytes to the internal
// pool.
func ReserveByteBuffer() *[]byte {
	b := bufPool.Get().(*[]byte)
	*b = (*b)[:0]
	return b
}

// ReleaseByteBuffer returns the given slice of bytes to the internal
// sync.Pool from which ReserveByteBuffer can fetch it again.
func ReleaseByteBuffer(buf *[]byte) {
	bufPool.Put(buf)
}
