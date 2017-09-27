package internal

import (
	"os"
	"path/filepath"
)

func Directory(file string) (files []string, err error) {
	info, err := os.Stat(file)
	if err != nil {
		return nil, err
	}
	if !info.IsDir() {
		return []string{filepath.Base(file)}, nil
	}
	f, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer func() {
		nerr := f.Close()
		if err == nil {
			err = nerr
		}
	}()
	return f.Readdirnames(0)
}

func FullPathname(filename string) (string, error) {
	if filepath.IsAbs(filename) {
		return filename, nil
	}
	wd, err := os.Getwd()
	return filepath.Join(wd, filename), err
}
