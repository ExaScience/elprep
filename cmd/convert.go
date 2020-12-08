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

package cmd

import (
	"flag"
	"os"

	"github.com/exascience/elprep/v5/fasta"
	"github.com/exascience/elprep/v5/intervals"
)

// VcfToElsitesHelp is the help string for this command.
const VcfToElsitesHelp = "vcf-to-elsites parameters:\n" +
	"elprep vcf-to-elsites vcf-file elsites-file\n" +
	"[--log-path path]\n"

// VcfToElsites implements the elprep vcf-to-elsites command.
func VcfToElsites() {
	var logPath string

	var flags flag.FlagSet
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	parseFlags(flags, 4, VcfToElsitesHelp)

	input := getFilename(os.Args[2], VcfToElsitesHelp)
	output := getFilename(os.Args[3], VcfToElsitesHelp)

	setLogOutput(logPath)

	inter := intervals.FromVcfFile(input)
	for chrom, ivals := range inter {
		intervals.ParallelSortByStart(ivals)
		inter[chrom] = intervals.ParallelFlatten(ivals)
	}
	intervals.ToElsitesFile(inter, output)
}

// BedToElsitesHelp is the help string for this command.
const BedToElsitesHelp = "\nbed-to-elsites parameters:\n" +
	"elprep bed-to-elsites bed-file elsites-file\n" +
	"[--log-path path]\n"

// BedToElsites implements the elprep bed-to-elsites command.
func BedToElsites() {
	var logPath string

	var flags flag.FlagSet
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	parseFlags(flags, 4, BedToElsitesHelp)

	input := getFilename(os.Args[2], BedToElsitesHelp)
	output := getFilename(os.Args[3], BedToElsitesHelp)

	setLogOutput(logPath)

	inter := intervals.FromBedFile(input)
	for chrom, ivals := range inter {
		intervals.ParallelSortByStart(ivals)
		inter[chrom] = intervals.ParallelFlatten(ivals)
	}
	intervals.ToElsitesFile(inter, output)
}

// FastaToElfastaHelp is the help string for this command.
const FastaToElfastaHelp = "fasta-to-elfasta parameters:\n" +
	"elprep fasta-to-elfasta fasta-file elfasta-file\n" +
	"[--log-path path]\n"

// FastaToElfasta implements the elprep fasta-to-elfasta command.
func FastaToElfasta() {
	var logPath string

	var flags flag.FlagSet
	flags.StringVar(&logPath, "log-path", "", "write log files to the specified directory")
	parseFlags(flags, 4, FastaToElfastaHelp)

	input := getFilename(os.Args[2], FastaToElfastaHelp)
	output := getFilename(os.Args[3], FastaToElfastaHelp)

	setLogOutput(logPath)

	fasta.ToElfasta(fasta.ParseFasta(input, nil, false, false), output)
}
