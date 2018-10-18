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

package cmd

import (
	"flag"
	"fmt"
	"os"

	"github.com/exascience/elprep/v4/fasta"
	"github.com/exascience/elprep/v4/intervals"
)

// VcfToElsitesHelp is the help string for this command.
const VcfToElsitesHelp = "vcf-to-elsites parameters:\n" +
	"elprep vcf-to-elsites vcf-file elsites-file\n" +
	"[--log-path path]\n"

// VcfToElsites implements the elprep vcf-to-elsites command.
func VcfToElsites() error {

	var logPath string

	var flags flag.FlagSet

	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	var input, output string

	if len(os.Args) < 4 {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprint(os.Stderr, VcfToElsitesHelp)
		os.Exit(1)
	}

	input = getFilename(os.Args[2], VcfToElsitesHelp)
	output = getFilename(os.Args[3], VcfToElsitesHelp)

	setLogOutput(logPath)

	inter, err := intervals.FromVcfFile(input)

	if err != nil {
		return err
	}

	for chrom, ivals := range inter {
		intervals.ParallelSortByStart(ivals)
		inter[chrom] = intervals.ParallelFlatten(ivals)
	}

	return intervals.ToElsitesFile(inter, output)
}

// BedToElsitesHelp is the help string for this command.
const BedToElsitesHelp = "\nbed-to-elsites parameters:\n" +
	"elprep bed-to-elsites bed-file elsites-file\n" +
	"[--log-path path]\n"

// BedToElsites implements the elprep bed-to-elsites command.
func BedToElsites() error {

	var logPath string

	var flags flag.FlagSet

	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	parseFlags(flags, 4, BedToElsitesHelp)

	input := getFilename(os.Args[2], BedToElsitesHelp)
	output := getFilename(os.Args[3], BedToElsitesHelp)

	setLogOutput(logPath)

	inter, err := intervals.FromBedFile(input)

	if err != nil {
		return err
	}

	for chrom, ivals := range inter {
		intervals.ParallelSortByStart(ivals)
		inter[chrom] = intervals.ParallelFlatten(ivals)
	}

	return intervals.ToElsitesFile(inter, output)
}

// FastaToElfastaHelp is the help string for this command.
const FastaToElfastaHelp = "fasta-to-elfasta parameters:\n" +
	"elprep fasta-to-elfasta fasta-file elfasta-file\n" +
	"[--log-path path]\n"

// FastaToElfasta implements the elprep fasta-to-elfasta command.
func FastaToElfasta() error {

	var logPath string

	var flags flag.FlagSet

	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")

	var input, output string

	if len(os.Args) < 4 {
		fmt.Fprintln(os.Stderr, "Incorrect number of parameters.")
		fmt.Fprintln(os.Stderr, FastaToElfastaHelp)
		os.Exit(1)
	}

	input = getFilename(os.Args[2], FastaToElfastaHelp)
	output = getFilename(os.Args[3], FastaToElfastaHelp)

	setLogOutput(logPath)

	fst, err := fasta.ParseFasta(input, nil, false, false)

	if err != nil {
		return err
	}

	return fasta.ToElfasta(fst, output)
}
