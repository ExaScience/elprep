package utils

import (
	"unsafe"

	"github.com/exascience/pargo/sync"

	"github.com/exascience/elprep/internal"
)

type (
	symbolName string
	Symbol     *string
)

func SymbolHash(s Symbol) uint64 {
	return uint64(uintptr(unsafe.Pointer(s)))
}

func (s symbolName) Hash() uint64 {
	return internal.StringHash(string(s))
}

var symbolTable = sync.NewMap(0)

func Intern(s string) Symbol {
	entry, _ := symbolTable.LoadOrStore(symbolName(s), Symbol(&s))
	return entry.(Symbol)
}
