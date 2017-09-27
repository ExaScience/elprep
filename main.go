package main

import (
	"fmt"
	"log"
	"os"

	"github.com/exascience/elprep/cmd"
)

func printHelp() {
	fmt.Fprintln(os.Stderr, "Available commands: filter, split, merge")
	fmt.Fprint(os.Stderr, "\n", cmd.FilterHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.SplitHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.MergeHelp)
}

func main() {
	fmt.Fprintln(os.Stderr, cmd.ProgramMessage)
	if len(os.Args) < 2 {
		log.Println("Incorrect number of parameters.")
		printHelp()
		os.Exit(1)
	} else {
		var err error
		switch os.Args[1] {
		case "filter":
			err = cmd.Filter()
		case "split":
			err = cmd.Split()
		case "merge":
			err = cmd.Merge()
		case "help", "-help", "--help", "-h", "--h":
			printHelp()
		default:
			err = cmd.DeprecatedFilter()
		}
		if err != nil {
			log.Fatal(err.Error())
		}
	}
}
