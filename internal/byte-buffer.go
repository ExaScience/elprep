package internal

import "sync"

var bufPool = sync.Pool{New: func() interface{} {
	return []byte(nil)
}}

func ReserveByteBuffer() []byte    {
	return bufPool.Get().([]byte)[:0]
}

func ReleaseByteBuffer(buf []byte) {
	bufPool.Put(buf)
}
