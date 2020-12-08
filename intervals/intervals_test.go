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

package intervals

import (
	"math/rand"
	"testing"
)

func intervalsEqual(intervals1, intervals2 []Interval) bool {
	if len(intervals1) != len(intervals2) {
		return false
	}
	for i, interval1 := range intervals1 {
		if interval1 != intervals2[i] {
			return false
		}
	}
	return true
}

func makeLargeIntervalsSlice() (result []Interval) {
	result = make([]Interval, 0x30000)
	result[0].Start = 0
	result[0].End = 3
	for i := 1; i < len(result); i++ {
		if rand.Intn(100) < 20 {
			result[i].Start = result[i-1].End - 1
		} else {
			result[i].Start = result[i-1].End + 1
		}
		result[i].End = result[i].Start + 3
	}
	return result
}

func TestFlatten(t *testing.T) {
	if Flatten(nil) != nil {
		t.Error("empty Flatten failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 3}, {3, 4}}), []Interval{{2, 4}}) {
		t.Error("Flatten 1 failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 3}, {4, 5}}), []Interval{{2, 3}, {4, 5}}) {
		t.Error("Flatten 2 failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 4}, {3, 5}, {4, 6}}), []Interval{{2, 6}}) {
		t.Error("Flatten 3 failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 4}, {3, 5}, {4, 6}, {7, 9}}), []Interval{{2, 6}, {7, 9}}) {
		t.Error("Flatten 4 failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 3}, {3, 4}, {5, 6}, {6, 7}}), []Interval{{2, 4}, {5, 7}}) {
		t.Error("Flatten 5 failed")
	}
	if !intervalsEqual(Flatten([]Interval{{2, 3}, {2, 5}, {2, 4}, {2, 3}, {2, 6}, {2, 7}}), []Interval{{2, 7}}) {
		t.Error("Flatten 6 failed")
	}
	intervals := Flatten(makeLargeIntervalsSlice())
	if intervals[0].Start > intervals[0].End {
		t.Error("Flatten 7a failed")
	}
	for i := 1; i < len(intervals); i++ {
		interval := intervals[i]
		if interval.Start > interval.End || interval.Start <= intervals[i-1].End {
			t.Error("Flatten 7b failed")
		}
	}
}
func TestParallelFlatten(t *testing.T) {
	if ParallelFlatten(nil) != nil {
		t.Error("empty ParallelFlatten failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 3}, {3, 4}}), []Interval{{2, 4}}) {
		t.Error("ParallelFlatten 1 failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 3}, {4, 5}}), []Interval{{2, 3}, {4, 5}}) {
		t.Error("ParallelFlatten 2 failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 4}, {3, 5}, {4, 6}}), []Interval{{2, 6}}) {
		t.Error("ParallelFlatten 3 failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 4}, {3, 5}, {4, 6}, {7, 9}}), []Interval{{2, 6}, {7, 9}}) {
		t.Error("ParallelFlatten 4 failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 3}, {3, 4}, {5, 6}, {6, 7}}), []Interval{{2, 4}, {5, 7}}) {
		t.Error("ParallelFlatten 5 failed")
	}
	if !intervalsEqual(ParallelFlatten([]Interval{{2, 3}, {2, 5}, {2, 4}, {2, 3}, {2, 6}, {2, 7}}), []Interval{{2, 7}}) {
		t.Error("ParallelFlatten 6 failed")
	}
	intervals := ParallelFlatten(makeLargeIntervalsSlice())
	if intervals[0].Start > intervals[0].End {
		t.Error("ParallelFlatten 7a failed")
	}
	for i := 1; i < len(intervals); i++ {
		interval := intervals[i]
		if interval.Start > interval.End || interval.Start <= intervals[i-1].End {
			t.Error("ParallelFlatten 7b failed")
		}
	}
}

func BenchmarkFlatten(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		intervals := makeLargeIntervalsSlice()
		b.StartTimer()
		intervals = Flatten(intervals)
	}
}
func BenchmarkParallelFlatten(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		intervals := makeLargeIntervalsSlice()
		b.StartTimer()
		intervals = ParallelFlatten(intervals)
	}
}

func TestOverlap(t *testing.T) {
	if Overlap(nil, 2, 3) {
		t.Error("empty Overlap failed")
	}
	if Overlap([]Interval{{1, 3}, {7, 8}}, 4, 6) {
		t.Error("Overlap 1 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 1, 3) {
		t.Error("Overlap 2 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 2, 3) {
		t.Error("Overlap 3 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 2, 5) {
		t.Error("Overlap 4 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 2, 6) {
		t.Error("Overlap 5 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 3, 7) {
		t.Error("Overlap 6 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 5, 7) {
		t.Error("Overlap 7 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 6, 8) {
		t.Error("Overlap 8 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 6, 9) {
		t.Error("Overlap 9 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 5, 9) {
		t.Error("Overlap 10 failed")
	}
	if !Overlap([]Interval{{2, 4}, {6, 8}}, 1, 10) {
		t.Error("Overlap 11 failed")
	}
}

func TestIntersect(t *testing.T) {
	if !intervalsEqual(Intersect(nil, 2, 3), nil) {
		t.Error("empty Intersect failed")
	}
	if !intervalsEqual(Intersect([]Interval{{1, 3}, {7, 8}}, 4, 6), nil) {
		t.Error("Intersect 1 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 1, 3), []Interval{{2, 4}}) {
		t.Error("Intersect 2 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 2, 3), []Interval{{2, 4}}) {
		t.Error("Intersect 3 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 2, 5), []Interval{{2, 4}}) {
		t.Error("Intersect 4 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 2, 6), []Interval{{2, 4}, {6, 8}}) {
		t.Error("Intersect 5 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 3, 7), []Interval{{2, 4}, {6, 8}}) {
		t.Error("Intersect 6 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 5, 7), []Interval{{6, 8}}) {
		t.Error("Intersect 7 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 6, 8), []Interval{{6, 8}}) {
		t.Error("Intersect 8 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 6, 9), []Interval{{6, 8}}) {
		t.Error("Intersect 9 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 5, 9), []Interval{{6, 8}}) {
		t.Error("Intersect 10 failed")
	}
	if !intervalsEqual(Intersect([]Interval{{2, 4}, {6, 8}}, 1, 10), []Interval{{2, 4}, {6, 8}}) {
		t.Error("Intersect 11 failed")
	}
}
