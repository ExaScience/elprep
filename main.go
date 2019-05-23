// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2017-2019 imec vzw.

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

// elPrep is a high-performance tool for preparing .sam/.bam
// files for variant calling in sequencing pipelines.
//
// Please see https://github.com/exascience/elprep for a documentation
// of the tool, and below (and/or
// https://godoc.org/github.com/ExaScience/elprep) for the API
// documentation.
package main

import (
	"fmt"
	"log"
	"os"

	"github.com/exascience/elprep/v4/cmd"
)

func printHelp() {
	fmt.Fprintln(os.Stderr, "Available commands: filter, sfm, vcf-to-elsites, bed-to-elsites, fasta-to-elfasta")
	fmt.Fprint(os.Stderr, "\n", cmd.CombinedSfmFilterHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.VcfToElsitesHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.BedToElsitesHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.FastaToElfastaHelp)
}

func prinExtendedHelp() {
	fmt.Fprintln(os.Stderr, "Available commands: filter, split, merge, merge-optical-duplicates-metrics, sfm, vcf-to-elsites, bed-to-elsites, fasta-to-elfasta")
	fmt.Fprint(os.Stderr, "\n", cmd.FilterExtendedHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.SplitHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.MergeHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.SfmHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.MergeOpticalDuplicatesMetricsHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.VcfToElsitesHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.BedToElsitesHelp)
	fmt.Fprint(os.Stderr, "\n", cmd.FastaToElfastaHelp)
}

func main() {
	fmt.Fprintln(os.Stderr, cmd.ProgramMessage)
	if len(os.Args) < 2 {
		log.Println("Incorrect number of parameters.")
		fmt.Fprintln(os.Stderr, cmd.HelpMessage)
		printHelp()
		os.Exit(1)
	}

	var err error
	switch os.Args[1] {
	case "filter":
		err = cmd.Filter()
	case "split":
		err = cmd.Split()
	case "merge":
		err = cmd.Merge()
	case "merge-optical-duplicates-metrics":
		err = cmd.MergeOpticalDuplicatesMetrics()
	case "vcf-to-elsites":
		err = cmd.VcfToElsites()
	case "bed-to-elsites":
		err = cmd.BedToElsites()
	case "fasta-to-elfasta":
		err = cmd.FastaToElfasta()
	case "sfm":
		err = cmd.Sfm()
	case "help", "-help", "--help", "-h", "--h":
		printHelp()
	case "help-extended", "-help-extended", "--help-extended", "-he", "--he":
		prinExtendedHelp()
	default:
		err = cmd.DeprecatedFilter()
	}
	if err != nil {
		log.Fatal(err)
	}

}
