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

package filters

import (
	"log"
	"sync"
)

/*
 A "deletion" is a piece of information we need in calculateGenotype.
 A challenge is that an assembly region may be influenced by deletions from
 another assembly region left to it, which strictly speaking creates a serial
 dependency between assembly regions. These serial dependencies are something
 we want to avoid. The "side channel" construction below helps to minimize
 the dependencies, so that even in the general case, some of the computations
 per assembly region can overlap.
*/

type (
	deletion struct {
		start, end int32
	}

	sideChannel struct {
		// input produces deletions from other assembly regions to the left
		// output receive deletions for other assembly regions to the right
		input, output chan interface{}
	}

	// Once we know we will need deletion information in an assembly region,
	// we can use a deletionsHandler to ensure even more computational overlap.
	deletionsHandler struct {
		wg    sync.WaitGroup
		ch    *sideChannel
		slice []deletion
	}
)

// Used when an assembly region is the first (left-most) in a contig.
// Closing the input ensures that the input produces nil.
func (ch *sideChannel) makeInitial() {
	ch.input = make(chan interface{})
	ch.output = make(chan interface{}, 1)
	close(ch.input)
}

// Used when an assembly region is not the first in a contig.
// Previous is from the assembly region directly to the left of
// the current one.
func (ch *sideChannel) linkFrom(previous sideChannel) {
	ch.input = previous.output
	ch.output = make(chan interface{}, 1)
}

// Receive deletions from assembly regions to the left.
// An assembly region may send its own input channel to
// its output channel in case it doesn't need any deletion
// information for itself, which will forward deletion
// information further from the left to the next
// assembly region to the right.
func (ch *sideChannel) receiveDeletions() []deletion {
	for {
		if item := <-ch.input; item == nil {
			return nil
		} else {
			switch it := item.(type) {
			case []deletion:
				return it
			case chan interface{}:
				ch.input = it
			default:
				log.Panicf("Invalid value %v received from side channel.", item)
			}
		}
	}
}

// When an assembly region is done processing and creating new deletions,
// it can send the new deletions to the next assembly region.
func (ch *sideChannel) sendDeletions(deletions []deletion) {
	ch.output <- deletions
	close(ch.output)
}

// When an assembly region determines that it doesn't need deletion information
// and won't create any of its own, it can call this method. The deletions
// from the left are then just forwarded to the right.
func (ch *sideChannel) noDeletions() {
	ch.output <- ch.input
	close(ch.output)
}

// Create a deletions handler for the current assembly region.
// Waiting for the deletions is handled in a separate goroutine to
// avoid unnecessary blocking.
// Deletions can be inspected and modified by calling handler.wg.Wait()
// and then just accessing handler.slice.
func (ch *sideChannel) handleDeletions() *deletionsHandler {
	w := &deletionsHandler{ch: ch}
	w.wg.Add(1)
	go func() {
		defer w.wg.Done()
		w.slice = ch.receiveDeletions()
	}()
	return w
}

// Closing a deletionsHandler ensures that the current deletions slice
// is properly sent to the next assembly region.
func (h *deletionsHandler) close() {
	h.wg.Wait()
	h.ch.sendDeletions(h.slice)
	h.ch = nil
	h.slice = nil
}
