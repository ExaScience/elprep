// elPrep: a high-performance tool for analyzing SAM/BAM files.
// Copyright (c) 2020 imec vzw.

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

package internal

import (
	"log"
	"os/exec"

	"github.com/exascience/pargo/pipeline"
)

var BQSRTablenamePrefix string

// RunPipeline is p.Run() with panics in place of errors
func RunPipeline(p *pipeline.Pipeline) {
	p.Run()
	if err := p.Err(); err != nil {
		log.Panic(err)
	}
}

// RunCmd is cmd.Run() with panics in place of errors
func RunCmd(cmd *exec.Cmd) {
	if err := cmd.Run(); err != nil {
		log.Panic(err)
	}
}
