package sam

import (
	"fmt"
)

type StringScanner struct {
	index int
	data  string
	err   error
}

func (sc *StringScanner) Err() error {
	return sc.err
}

func (sc *StringScanner) Reset(s string) {
	sc.index = 0
	sc.data = s
	sc.err = nil
}

func (sc *StringScanner) Len() int {
	if sc.err != nil {
		return 0
	}
	return len(sc.data) - sc.index
}

func (sc *StringScanner) readByteUntil(c byte) (b byte, found bool) {
	if sc.err != nil {
		return 0, false
	}
	start := sc.index
	next := start + 1
	if next >= len(sc.data) {
		sc.index = len(sc.data)
		return sc.data[start], false
	} else if sc.data[next] != c {
		if sc.err == nil {
			sc.err = fmt.Errorf("Unexpected character %v in StringScanner.ReadByteUntil", sc.data[next])
		}
		return 0, false
	} else {
		sc.index = next + 1
		return sc.data[start], true
	}
}

func (sc *StringScanner) readUntil(c byte) (s string, found bool) {
	if sc.err != nil {
		return "", false
	}
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if sc.data[end] == c {
			sc.index = end + 1
			return sc.data[start:end], true
		}
	}
	sc.index = len(sc.data)
	return sc.data[start:], false
}

func (sc *StringScanner) readUntil2(c1, c2 byte) (s string, b byte) {
	if sc.err != nil {
		return "", 0
	}
	start := sc.index
	for end := sc.index; end < len(sc.data); end++ {
		if c := sc.data[end]; (c == c1) || (c == c2) {
			sc.index = end + 1
			return sc.data[start:end], c
		}
	}
	sc.index = len(sc.data)
	return sc.data[start:], 0
}
