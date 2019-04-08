// elPrep: a high-performance tool for preparing SAM/BAM files.
// Copyright (c) 2019 imec vzw.

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

// A graph data structure for determining optical distance relations
// in the graph-based algorithm for marking optical duplicates.

type graph [][]int // adjacency lists for the edges

func newGraph(size int) graph {
	return make([][]int, size)
}

func (g graph) addNeighbor(from, to int) {
	n := g[from]
	for _, t := range n {
		if t == to {
			return
		}
	}
	g[from] = append(n, to)
}

func (g graph) addEdge(left, right int) {
	if left == right {
		return
	}
	g.addNeighbor(left, right)
	g.addNeighbor(right, left)
}

func findRepNode(grouping []int, nodeId int) int {
	representativeUmi := nodeId
	for representativeUmi != grouping[representativeUmi] {
		representativeUmi = grouping[representativeUmi]
	}
	for nodeId != representativeUmi {
		newUmiId := grouping[nodeId]
		grouping[nodeId] = representativeUmi
		nodeId = newUmiId
	}
	return representativeUmi
}

func joinNodes(grouping []int, nodeId1, nodeId2 int) {
	repNode1 := findRepNode(grouping, nodeId1)
	repNode2 := findRepNode(grouping, nodeId2)
	if repNode1 == repNode2 {
		return
	}
	grouping[repNode1] = repNode2
}

func (g graph) cluster() map[int]int {
	cluster := make([]int, len(g))
	for i := range cluster {
		cluster[i] = i
	}
	for i := range g {
		for _, j := range g[i] {
			joinNodes(cluster, j, i)
		}
	}
	result := make(map[int]int, len(g))
	for index := range g {
		result[index] = findRepNode(cluster, index)
	}
	return result
}
