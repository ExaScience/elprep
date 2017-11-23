package utils

import (
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/internal"
)

type symbolName string

// A Symbol is a unique pointer to a string.
type Symbol *string

/*
SymbolHash computes a hash value for the given Symbol.
*/
func SymbolHash(s Symbol) uint64 {
	return uint64(uintptr(unsafe.Pointer(s)))
}

func (s symbolName) Hash() uint64 {
	return internal.StringHash(string(s))
}

var symbolTable = sync.NewMap(0)

/*
Intern returns a Symbol for the given string.

It always returns the same pointer for strings that are equal, and
different pointers for strings that are not equal. So for two strings
s1 and s2, if s1 == s2, then Intern(s1) == Intern(s2), and if s1 !=
s2, then Intern(s1) != Intern(s2).

Dereferencing the pointer always yields a string that is equal to the
original string: *Intern(s) == s always holds.

It is safe for multiple goroutines to call Intern concurrently.
*/
func Intern(s string) Symbol {
	entry, _ := symbolTable.LoadOrStore(symbolName(s), Symbol(&s))
	return entry.(Symbol)
}
