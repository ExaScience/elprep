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
	"math"
	"sort"
	"strings"

	"github.com/exascience/elprep/v5/fasta"

	"github.com/exascience/elprep/v5/sam"
)

func makeReferenceHaplotype(bases string, location int32) *haplotype {
	return &haplotype{
		bases:    bases,
		location: location,
		cigar:    []sam.CigarOperation{{int32(len(bases)), 'M'}},
		isRef:    true,
		score:    math.NaN(),
	}
}

type kmer struct {
	bases       string
	start, stop int32
	isRef       bool
}

func (hc *HaplotypeCaller) baseUseableForAssembly(base, qual byte) bool {
	return base != 'N' && qual >= hc.minBaseQual
}

func (hc *HaplotypeCaller) addSequencesForKmers(sequences []kmer, aln *sam.Alignment, kmerSize int32) []kmer {
	start := int32(-1)
	alnSequence := aln.SEQ.AsString()
	end := int32(len(alnSequence))
	for stop := int32(0); stop < end; stop++ {
		if !hc.baseUseableForAssembly(alnSequence[stop], aln.QUAL[stop]) {
			if start != -1 && stop-start >= kmerSize {
				sequences = append(sequences, kmer{
					bases: alnSequence,
					start: start,
					stop:  stop,
				})
			}
			start = -1
		} else if start == -1 {
			start = stop
		}
	}
	if start != -1 && end-start >= kmerSize {
		sequences = append(sequences, kmer{
			bases: alnSequence,
			start: start,
			stop:  end,
		})
	}
	return sequences
}

func makeSequenceForKmersFromReference(sequence string, kmerSize int32) kmer {
	return kmer{
		bases: sequence,
		start: 0,
		stop:  int32(len(sequence)),
		isRef: true,
	}
}

type (
	vertexInfo struct {
		id    int32
		bases string
	}

	edgeInfo struct {
		id, multiplicity int32
		isRef            bool
	}

	kmerGraph struct {
		verticesId     int32
		kmerSize       int32
		vertices       map[int32]*vertexInfo
		uniqueKmers    map[string]*vertexInfo
		nonUniqueKmers map[string]bool
		outgoingEdges  map[int32][]*edgeInfo
		incomingEdges  map[int32][]*edgeInfo
	}
)

func (g *kmerGraph) vertexOutDegree(vertex *vertexInfo) int {
	return len(g.outgoingEdges[vertex.id])
}

func (g *kmerGraph) vertexInDegree(vertex *vertexInfo) int {
	return len(g.incomingEdges[vertex.id])
}

func newKmerGraph(kmerSize int32) *kmerGraph {
	return &kmerGraph{
		kmerSize:       kmerSize,
		vertices:       make(map[int32]*vertexInfo),
		uniqueKmers:    make(map[string]*vertexInfo),
		nonUniqueKmers: make(map[string]bool),
		outgoingEdges:  make(map[int32][]*edgeInfo),
		incomingEdges:  make(map[int32][]*edgeInfo),
	}
}

func (g *kmerGraph) newVertexId() int32 {
	g.verticesId++
	return g.verticesId
}

func (g *kmerGraph) addVertex(vertex *vertexInfo) {
	vertex.id = g.newVertexId()
	g.vertices[vertex.id] = vertex
}

func (g *kmerGraph) updateVertexId(vertex *vertexInfo) {
	if vertex.id == g.verticesId {
		return
	}
	oldId := vertex.id
	g.verticesId++
	newId := g.verticesId
	for _, incoming := range g.incomingEdges[oldId] {
		for _, outgoing := range g.outgoingEdges[incoming.id] {
			if outgoing.id == oldId {
				outgoing.id = newId
			}
		}
	}
	for _, outgoing := range g.outgoingEdges[oldId] {
		for _, incoming := range g.incomingEdges[outgoing.id] {
			if incoming.id == oldId {
				incoming.id = newId
			}
		}
	}
	incomingEdges := g.incomingEdges[oldId]
	delete(g.incomingEdges, oldId)
	g.incomingEdges[newId] = incomingEdges
	outgoingEdges := g.outgoingEdges[oldId]
	delete(g.outgoingEdges, oldId)
	g.outgoingEdges[newId] = outgoingEdges
	vertex.id = newId
	delete(g.vertices, oldId)
	g.vertices[newId] = vertex
}

type kmerGraphCopy struct {
	nofVertices   int
	vertexBases   map[string]bool
	outgoingEdges [][2]int32
	incomingEdges [][2]int32
}

func int32PairLess(pairs [][2]int32, i, j int) bool {
	x := pairs[i]
	y := pairs[j]
	if x[0] < y[0] {
		return true
	}
	if x[0] > y[0] {
		return false
	}
	return x[1] < y[1]
}

func partialCopyEdges(edges map[int32][]*edgeInfo) (result [][2]int32) {
	for id, edgeSlice := range edges {
		for _, edge := range edgeSlice {
			result = append(result, [2]int32{id, edge.id})
		}
	}
	sort.Slice(result, func(i, j int) bool {
		return int32PairLess(result, i, j)
	})
	return
}

func (g *kmerGraph) partialCopy() (result kmerGraphCopy) {
	result.nofVertices = len(g.vertices)
	result.vertexBases = make(map[string]bool)
	for _, vertex := range g.vertices {
		result.vertexBases[vertex.bases] = true
	}
	result.outgoingEdges = partialCopyEdges(g.outgoingEdges)
	result.incomingEdges = partialCopyEdges(g.incomingEdges)
	return
}

func (prev kmerGraphCopy) partialEqual(g kmerGraphCopy) bool {
	if prev.nofVertices != g.nofVertices {
		return false
	}
	if len(prev.outgoingEdges) != len(g.outgoingEdges) {
		return false
	}
	if len(prev.incomingEdges) != len(g.incomingEdges) {
		return false
	}
	for prevBases := range prev.vertexBases {
		if !g.vertexBases[prevBases] {
			return false
		}
	}
	for i, prevEdge := range prev.outgoingEdges {
		if prevEdge != g.outgoingEdges[i] {
			return false
		}
	}
	for i, prevEdge := range prev.incomingEdges {
		if prevEdge != g.incomingEdges[i] {
			return false
		}
	}
	return true
}

func (g *kmerGraph) setOutgoingEdges(vertex *vertexInfo, edges []*edgeInfo) {
	if len(edges) == 0 {
		delete(g.outgoingEdges, vertex.id)
	} else {
		g.outgoingEdges[vertex.id] = edges
	}
}

func (g *kmerGraph) setIncomingEdges(vertex *vertexInfo, edges []*edgeInfo) {
	if len(edges) == 0 {
		delete(g.incomingEdges, vertex.id)
	} else {
		g.incomingEdges[vertex.id] = edges
	}
}

func (g *kmerGraph) getOutgoingEdge(source, target *vertexInfo) (*edgeInfo, bool) {
	for _, edge := range g.outgoingEdges[source.id] {
		if edge.id == target.id {
			return edge, true
		}
	}
	return nil, false
}

func (g *kmerGraph) addEdge(vertex1, vertex2 *vertexInfo, multiplicity int32, isRef bool) (incoming, outgoing *edgeInfo) {
	if _, ok := g.getOutgoingEdge(vertex1, vertex2); ok {
		return
	}
	incoming = &edgeInfo{vertex1.id, multiplicity, isRef}
	g.incomingEdges[vertex2.id] = append(g.incomingEdges[vertex2.id], incoming)
	outgoing = &edgeInfo{vertex2.id, multiplicity, isRef}
	g.outgoingEdges[vertex1.id] = append(g.outgoingEdges[vertex1.id], outgoing)
	return
}

func (g *kmerGraph) getHeaviestOutgoingEdge(vertex *vertexInfo) *edgeInfo {
	edges := g.outgoingEdges[vertex.id]
	maxEdge := edges[0]
	for _, edge := range edges[1:] {
		if edge.multiplicity > maxEdge.multiplicity {
			maxEdge = edge
		}
	}
	return maxEdge
}

func (g *kmerGraph) isSingletonVertex(vertex *vertexInfo) bool {
	return g.vertexInDegree(vertex) == 0 && g.vertexOutDegree(vertex) == 0
}

func (g *kmerGraph) removeSingletonVertex(vertex *vertexInfo) {
	if vertex.id == -1 {
		return
	}
	delete(g.vertices, vertex.id)
	delete(g.uniqueKmers, vertex.bases)
	vertex.id = -1
}

func (g *kmerGraph) removeEdgeRaw(source, target *vertexInfo) {
	var newOutgoingEdges []*edgeInfo
	for _, edge := range g.outgoingEdges[source.id] {
		if edge.id != target.id {
			newOutgoingEdges = append(newOutgoingEdges, edge)
		}
	}
	g.setOutgoingEdges(source, newOutgoingEdges)
	var newIncomingEdges []*edgeInfo
	for _, edge := range g.incomingEdges[target.id] {
		if edge.id != source.id {
			newIncomingEdges = append(newIncomingEdges, edge)
		}
	}
	g.setIncomingEdges(target, newIncomingEdges)
}

func (g *kmerGraph) removeEdge(source, target *vertexInfo) {
	if target == nil {
		return
	}
	g.removeEdgeRaw(source, target)
	if g.isSingletonVertex(target) {
		g.removeSingletonVertex(target)
	}
	if g.isSingletonVertex(source) {
		if len(g.vertices) != 1 {
			g.removeSingletonVertex(source)
		}
	}
}

func (g *kmerGraph) removeAllOutgoingEdges(source *vertexInfo) {
	for _, edge := range g.outgoingEdges[source.id] {
		target := g.vertices[edge.id]
		var newIncomingEdges []*edgeInfo
		for _, edge := range g.incomingEdges[target.id] {
			if edge.id != source.id {
				newIncomingEdges = append(newIncomingEdges, edge)
			}
		}
		g.setIncomingEdges(target, newIncomingEdges)
		if g.isSingletonVertex(target) {
			g.removeSingletonVertex(target)
		}
	}
	delete(g.outgoingEdges, source.id)
	if g.isSingletonVertex(source) {
		if len(g.vertices) != 1 {
			g.removeSingletonVertex(source)
		}
	}
}

func (g *kmerGraph) removeAllIncomingEdges(target *vertexInfo) {
	for _, edge := range g.incomingEdges[target.id] {
		source := g.vertices[edge.id]
		var newOutgoingEdges []*edgeInfo
		for _, edge := range g.outgoingEdges[source.id] {
			if edge.id != target.id {
				newOutgoingEdges = append(newOutgoingEdges, edge)
			}
		}
		g.setOutgoingEdges(source, newOutgoingEdges)
		if g.isSingletonVertex(source) {
			g.removeSingletonVertex(source)
		}
	}
	delete(g.incomingEdges, target.id)
	if g.isSingletonVertex(target) {
		if len(g.vertices) != 1 {
			g.removeSingletonVertex(target)
		}
	}
}

func (g *kmerGraph) removeVertex(vertex *vertexInfo) {
	g.removeAllOutgoingEdges(vertex)
	g.removeAllIncomingEdges(vertex)
	g.removeSingletonVertex(vertex)
}

func (vertex *vertexInfo) getSuffix() byte {
	return vertex.bases[len(vertex.bases)-1]
}

func (g *kmerGraph) getVertices(predicate func(*vertexInfo) bool) (result []*vertexInfo) {
	for _, vertex := range g.vertices {
		if predicate(vertex) {
			result = append(result, vertex)
		}
	}
	sort.Slice(result, func(i, j int) bool {
		return result[i].id < result[j].id
	})
	return
}

func (g *kmerGraph) getAllVertices() []*vertexInfo {
	result := make([]*vertexInfo, 0, len(g.vertices))
	for _, vertex := range g.vertices {
		result = append(result, vertex)
	}
	sort.Slice(result, func(i, j int) bool {
		return result[i].id < result[j].id
	})
	return result
}

func (g *kmerGraph) getVertex(predicate func(*vertexInfo) bool) *vertexInfo {
	vertices := g.getAllVertices()
	for _, vertex := range vertices {
		if predicate(vertex) {
			return vertex
		}
	}
	return nil
}

func (g *kmerGraph) getNonReferenceDestinations() []*vertexInfo {
	return g.getVertices(func(vertex *vertexInfo) bool {
		return g.vertexOutDegree(vertex) == 0 && !g.vertexIsReferenceSink(vertex)
	})
}

func (g *kmerGraph) getNonReferenceStarts() []*vertexInfo {
	return g.getVertices(func(vertex *vertexInfo) bool {
		return g.vertexInDegree(vertex) == 0 && !g.vertexIsReferenceSource(vertex)
	})
}

func (g *kmerGraph) getReferenceSourceVertex() *vertexInfo {
	return g.getVertex(g.vertexIsReferenceSource)
}

func (g *kmerGraph) getReferenceSinkVertex() *vertexInfo {
	return g.getVertex(g.vertexIsReferenceSink)
}

const (
	processing = 1
	done       = 2
)

type cycleDetector struct {
	g        *kmerGraph
	vertices []*vertexInfo
	seen     map[*vertexInfo]int
	stack    []*vertexInfo
	path     []*vertexInfo
}

func newCycleDetector(g *kmerGraph) *cycleDetector {
	return &cycleDetector{
		g:        g,
		vertices: g.getAllVertices(),
		seen:     make(map[*vertexInfo]int),
	}
}

func (cd *cycleDetector) popVertex() (vertex *vertexInfo) {
	vertex = cd.vertices[0]
	cd.vertices = cd.vertices[1:]
	return
}

func (cd *cycleDetector) peekStack() (vertex *vertexInfo) {
	return cd.stack[len(cd.stack)-1]
}

func (cd *cycleDetector) popStack() (vertex *vertexInfo) {
	last := len(cd.stack) - 1
	vertex = cd.stack[last]
	cd.stack = cd.stack[:last]
	return
}

func (cd *cycleDetector) doublePopStack() {
	last := len(cd.stack) - 2
	cd.stack = cd.stack[:last]
}

func (cd *cycleDetector) encounterVertex(vertex *vertexInfo) {
	cd.seen[vertex] = processing
	cd.stack = append(cd.stack, vertex)
}

func (cd *cycleDetector) componentExhausted() bool {
	for {
		if len(cd.stack) == 0 {
			return true
		}
		if cd.peekStack() != nil {
			return false
		}
		cd.doublePopStack()
	}
}

func (cd *cycleDetector) hasNext() bool {
	if cd.componentExhausted() {
		for len(cd.vertices) > 0 {
			if v := cd.popVertex(); cd.seen[v] == 0 {
				cd.encounterVertex(v)
				return true
			}
		}
		return false
	}
	return true
}

func (cd *cycleDetector) nextVertex() *vertexInfo {
	var v *vertexInfo
	for {
		v = cd.popStack()
		if v == nil {
			cd.popStack()
		} else {
			break
		}
	}
	cd.stack = append(cd.stack, v, nil)
	cd.seen[v] = done
outer:
	for i := len(cd.path) - 1; i >= 0; i-- {
		for _, edge := range cd.g.outgoingEdges[cd.path[i].id] {
			if edge.id == v.id {
				break outer
			}
		}
		cd.path = cd.path[:i]
	}
	cd.path = append(cd.path, v)
	return v
}

func (cd *cycleDetector) encounterVertexAgain(vertex *vertexInfo) (hasCycle bool) {
	for _, v := range cd.path {
		if v == vertex {
			return true
		}
	}
	if cd.seen[vertex] == processing {
		for i := len(cd.stack) - 1; i >= 0; i-- {
			if cd.stack[i] == vertex {
				cd.stack = append(cd.stack[:i], cd.stack[i+1:]...)
				break
			}
		}
		cd.stack = append(cd.stack, vertex)
	}
	return false
}

func (cd *cycleDetector) addNextVertices(vertex *vertexInfo) (hasCycle bool) {
	for _, edge := range cd.g.outgoingEdges[vertex.id] {
		if target := cd.g.vertices[edge.id]; cd.seen[target] > 0 {
			if cd.encounterVertexAgain(target) {
				return true
			}
		} else {
			cd.encounterVertex(target)
		}
	}
	return false
}

func (cd *cycleDetector) processNext() (hasCycle bool) {
	return cd.addNextVertices(cd.nextVertex())
}

func (g *kmerGraph) hasCycle() bool {
	cd := newCycleDetector(g)
	cd.encounterVertex(cd.popVertex())
	for cd.hasNext() {
		if cd.processNext() {
			return true
		}
	}
	return false
}

type chainInfo struct {
	source *vertexInfo
	edges  []edgeInfo
}

func (g *kmerGraph) findChain(start *vertexInfo, edge *edgeInfo) []edgeInfo {
	end := g.vertices[edge.id]
	chain := []edgeInfo{*edge}
	for {
		if g.vertexInDegree(end) > 1 || start == end {
			return chain
		}
		endEdges := g.outgoingEdges[end.id]
		if len(endEdges) != 1 {
			return chain
		}
		nextEdge := endEdges[0]
		chain = append(chain, *nextEdge)
		end = g.vertices[nextEdge.id]
	}
}

func (g *kmerGraph) findAllChains() (chains []chainInfo) {
	sources := g.getVertices(g.isSourceVertex)
	seenVertices := make(map[int32]bool)
	for _, vertex := range sources {
		seenVertices[vertex.id] = true
	}
	for len(sources) > 0 {
		last := len(sources) - 1
		source := sources[last]
		sources = sources[:last]
		for _, edge := range g.outgoingEdges[source.id] {
			chain := g.findChain(source, edge)
			chains = append(chains, chainInfo{source, chain})
			lastEdgeId := chain[len(chain)-1].id
			if !seenVertices[lastEdgeId] {
				sources = append(sources, g.vertices[lastEdgeId])
				seenVertices[lastEdgeId] = true
			}
		}
	}
	return
}

const minPruningFactor = 2

func (g *kmerGraph) pruneChain(chain chainInfo) bool {
	for _, edge := range chain.edges {
		if edge.multiplicity >= minPruningFactor || edge.isRef {
			return false
		}
	}
	return true
}

func (g *kmerGraph) removeChain(chain chainInfo) {
	source := chain.source
	for _, edge := range chain.edges {
		target := g.vertices[edge.id]
		g.removeEdgeRaw(source, target)
		source = target
	}
	for _, vertex := range g.getAllVertices() {
		if len(g.vertices) == 1 {
			break
		}
		if g.isSingletonVertex(vertex) {
			g.removeSingletonVertex(vertex)
		}
	}
}

func (g *kmerGraph) pruneChainsWithLowWeight() {
	for _, chain := range g.findAllChains() {
		if g.pruneChain(chain) {
			g.removeChain(chain)
		}
	}
}

func (vertex *vertexInfo) isDeleted() bool {
	return vertex.id == -1
}

func (g *kmerGraph) mergeOutgoingEdges(vertex1, vertex2 *vertexInfo) {
	edges2 := g.outgoingEdges[vertex2.id]
	for _, edge := range edges2 {
		incomingEdges := g.incomingEdges[edge.id]
		for i, inEdge := range incomingEdges {
			if inEdge.id == vertex2.id {
				inEdge.id = vertex1.id
				g.incomingEdges[edge.id] = append(append(incomingEdges[:i], incomingEdges[i+1:]...), inEdge)
				break
			}
		}
	}
	for _, edge := range g.incomingEdges[vertex1.id] {
		outgoingEdges := g.outgoingEdges[edge.id]
		for i, outEdge := range outgoingEdges {
			if outEdge.id == vertex1.id {
				g.outgoingEdges[edge.id] = append(append(outgoingEdges[:i], outgoingEdges[i+1:]...), outEdge)
				break
			}
		}
	}
	g.setOutgoingEdges(vertex1, edges2)
	delete(g.outgoingEdges, vertex2.id)
	delete(g.incomingEdges, vertex2.id)
	g.removeSingletonVertex(vertex2)
}

func (g *kmerGraph) mergeVertices(vertex1, vertex2 *vertexInfo) {
	vertex1.bases = vertex1.bases + vertex2.bases
	g.mergeOutgoingEdges(vertex1, vertex2)
	g.updateVertexId(vertex1)
}

func (g *kmerGraph) mergeLinearChains() (modified bool) {
	for _, vertex := range g.getAllVertices() {
		if vertex.isDeleted() || !g.vertexIsLinearChainStart(vertex) {
			continue
		}
		prevVertex := vertex
		prevVertexIsReference := g.vertexIsReferenceNode(vertex)
		for edges := g.outgoingEdges[vertex.id]; len(edges) == 1; edges = g.outgoingEdges[vertex.id] {
			nextVertex := g.vertices[edges[0].id]
			if prevVertex == nextVertex ||
				g.vertexInDegree(nextVertex) != 1 ||
				prevVertexIsReference != g.vertexIsReferenceNode(nextVertex) {
				break
			}
			modified = true
			g.mergeVertices(vertex, nextVertex)
		}
	}
	return
}

func nonUniqueKmersExist(bases string, kmerSize int32) bool {
	kmersSeenSoFar := make(map[string]bool)
	for i, end := int32(0), int32(len(bases))-kmerSize; i <= end; i++ {
		kmer := bases[i : i+kmerSize]
		if kmersSeenSoFar[kmer] {
			return true
		}
		kmersSeenSoFar[kmer] = true
	}
	return false
}

func (g *kmerGraph) determineNonUniqueKmersForSequence(sequence kmer, kmerSize int32) {
	kmersSeenSoFar := make(map[string]bool)
	for i, end := int32(0), sequence.stop-kmerSize; i <= end; i++ { // from 0!!!
		kmer := sequence.bases[i : i+kmerSize]
		if kmersSeenSoFar[kmer] {
			g.nonUniqueKmers[kmer] = true
		} else {
			kmersSeenSoFar[kmer] = true
		}
	}
}

func (g *kmerGraph) initializeNonUniqueKmers(seqs []kmer, kmerSize int32) {
	for _, seq := range seqs {
		g.determineNonUniqueKmersForSequence(seq, kmerSize)
	}
}

func (g *kmerGraph) isStartingKmer(sequence string) bool {
	return !g.nonUniqueKmers[sequence]
}

func (g *kmerGraph) findStartOfKmers(sequence kmer) int32 {
	if sequence.isRef {
		return 0
	}
	for i, end := sequence.start, sequence.stop-g.kmerSize; i < end; i++ {
		if g.isStartingKmer(sequence.bases[i : i+g.kmerSize]) {
			return i
		}
	}
	return -1
}

func (g *kmerGraph) incrementOutgoingEdgeMultiplicity(from, to *vertexInfo) {
	for _, edge := range g.outgoingEdges[from.id] {
		if edge.id == to.id {
			edge.multiplicity += 1
			return
		}
	}
}

func (g *kmerGraph) incrementIncomingEdgeMultiplicity(from, to *vertexInfo) {
	for _, edge := range g.incomingEdges[to.id] {
		if edge.id == from.id {
			edge.multiplicity += 1
			return
		}
	}
}

func (g *kmerGraph) increaseCountsMatchedKmers(sequence kmer, original *vertexInfo) {
	var recur func(*vertexInfo, int32)
	recur = func(vertex *vertexInfo, offset int32) {
		if offset == -1 {
			return
		}
		if g.vertexInDegree(vertex) == 1 {
			for _, edge := range g.incomingEdges[vertex.id] {
				previous := g.vertices[edge.id]
				if previous.getSuffix() == original.bases[offset] {
					edge.multiplicity += 1
					g.incrementOutgoingEdgeMultiplicity(previous, vertex)
					recur(previous, offset-1)
				}
			}
		}
	}
	recur(original, g.kmerSize-2)
}

func (g *kmerGraph) newKmerVertex(kmerSequence string) (vertex *vertexInfo) {
	vertex = &vertexInfo{bases: kmerSequence}
	g.addVertex(vertex)
	if !g.nonUniqueKmers[kmerSequence] && g.uniqueKmers[kmerSequence] == nil {
		g.uniqueKmers[kmerSequence] = vertex
	}
	return
}

func (g *kmerGraph) getKmerVertex(sequence kmer, start int32) (vertex *vertexInfo) {
	kmerSequence := sequence.bases[start : start+g.kmerSize]
	if vertex = g.uniqueKmers[kmerSequence]; vertex == nil {
		vertex = g.newKmerVertex(kmerSequence)
	}
	return
}

func (g *kmerGraph) getKmerVertexButNotRefSource(refSource string, sequence kmer, start int32) (vertex *vertexInfo) {
	kmerSequence := sequence.bases[start : start+g.kmerSize]
	if kmerSequence == refSource {
		vertex = g.newKmerVertex(kmerSequence)
	} else if vertex = g.uniqueKmers[kmerSequence]; vertex == nil {
		vertex = g.newKmerVertex(kmerSequence)
	}
	return
}

func (g *kmerGraph) vertexIsReferenceSource(vertex *vertexInfo) bool {
	if len(g.vertices) == 1 {
		return true
	}
	for _, incoming := range g.incomingEdges[vertex.id] {
		if incoming.isRef {
			return false
		}
	}
	for _, outgoing := range g.outgoingEdges[vertex.id] {
		if outgoing.isRef {
			return true
		}
	}
	return false
}

func (g *kmerGraph) vertexIsReferenceSink(vertex *vertexInfo) bool {
	if len(g.vertices) == 1 {
		return true
	}
	for _, outgoing := range g.outgoingEdges[vertex.id] {
		if outgoing.isRef {
			return false
		}
	}
	for _, incoming := range g.incomingEdges[vertex.id] {
		if incoming.isRef {
			return true
		}
	}
	return false
}

func (g *kmerGraph) vertexIsLinearChainStart(vertex *vertexInfo) bool {
	if g.vertexOutDegree(vertex) != 1 {
		return false
	}
	if g.vertexInDegree(vertex) != 1 {
		return true
	}
	incoming := g.incomingEdges[vertex.id]
	return len(incoming) == 1 && g.vertexOutDegree(g.vertices[incoming[0].id]) > 1
}

func (g *kmerGraph) vertexIsReferenceNode(vertex *vertexInfo) bool {
	if len(g.vertices) == 1 {
		return true
	}
	for _, incoming := range g.incomingEdges[vertex.id] {
		if incoming.isRef {
			return true
		}
	}
	for _, outgoing := range g.outgoingEdges[vertex.id] {
		if outgoing.isRef {
			return true
		}
	}
	return false
}

func reversePath(path []*vertexInfo) []*vertexInfo {
	for i, j := 0, len(path)-1; i < j; i, j = i+1, j-1 {
		path[i], path[j] = path[j], path[i]
	}
	return path
}

func (g *kmerGraph) findPathUpwardsToLowestCommonAncestor(vertex *vertexInfo) (path []*vertexInfo) {
	currentVertex := vertex
	for g.vertexInDegree(currentVertex) == 1 && g.vertexOutDegree(currentVertex) < 2 {
		edge := g.incomingEdges[currentVertex.id][0]
		targetVertex := g.vertices[edge.id]
		if edge.multiplicity < minPruningFactor {
			path = path[:0]
		} else {
			path = append(path, currentVertex)
		}
		currentVertex = targetVertex
	}
	if g.vertexOutDegree(currentVertex) > 1 {
		return reversePath(append(path, currentVertex))
	}
	return nil
}

func (g *kmerGraph) findReferencePath(path []*vertexInfo) []*vertexInfo {
	vertex := path[0]
	var maxEdge *edgeInfo
	{
		edges := g.incomingEdges[path[1].id]
		maxEdge = edges[0]
		for _, edge := range edges[1:] {
			if edge.multiplicity > maxEdge.multiplicity {
				maxEdge = edge
			}
		}
		if maxEdge.id == vertex.id {
			maxEdge, _ = g.getOutgoingEdge(vertex, path[1])
		} else {
			maxEdge = nil
		}
	}

	path = nil
loop:
	for {
		path = append(path, vertex)
		edges := g.outgoingEdges[vertex.id]
		if len(edges) == 0 {
			return path
		}
		for _, edge := range edges {
			if edge.isRef {
				vertex = g.vertices[edge.id]
				continue loop
			}
		}
		if maxEdge == nil {
			if len(edges) == 1 {
				vertex = g.vertices[edges[0].id]
				continue loop
			}
			return path
		}
		var nextVertex *vertexInfo
		for _, edge := range edges {
			if edge != maxEdge {
				if nextVertex == nil {
					nextVertex = g.vertices[edge.id]
				} else {
					return path
				}
			}
		}
		if nextVertex != nil {
			vertex = nextVertex
			continue loop
		}
		return path
	}
}

func getReferenceBasesForPath(path []*vertexInfo) string {
	var result strings.Builder
	for _, vertex := range path {
		result.WriteByte(vertex.getSuffix())
	}
	return result.String()
}

func (g *kmerGraph) isSourceVertex(vertex *vertexInfo) bool {
	return g.vertexInDegree(vertex) == 0
}

func (g *kmerGraph) getReferenceBasesForPathWithExpandedSources(path []*vertexInfo) string {
	var result strings.Builder
	for _, vertex := range path {
		if g.isSourceVertex(vertex) {
			for i := len(vertex.bases) - 1; i >= 0; i-- {
				result.WriteByte(vertex.bases[i])
			}
		} else {
			result.WriteByte(vertex.getSuffix())
		}
	}
	return result.String()
}

func longestSuffixMatch(sequence, kmer string, start int32) int32 {
	kmerLength := int32(len(kmer))
	for length := int32(1); length <= kmerLength; length++ {
		if seqi := start - length + 1; seqi < 0 || sequence[seqi] != kmer[kmerLength-length] {
			return length - 1
		}
	}
	return kmerLength
}

const maxCigarComplexity = 3

func (g *kmerGraph) mergeDanglingTail(altPath, refPath []*vertexInfo, altBases, refBases string, cigar []sam.CigarOperation) bool {
	lastRefIndex := sam.ReferenceLengthFromCigar(cigar) - 1
	matchingSuffix := minInt32(longestSuffixMatch(refBases, altBases, lastRefIndex), cigar[len(cigar)-1].Length)
	if matchingSuffix == 0 {
		return false
	}
	altIndexToMerge := maxInt32(sam.ReadLengthFromCigar(cigar)-matchingSuffix-1, 0)
	refIndexToMerge := lastRefIndex - matchingSuffix + 1
	if cigar[0].Operation == 'D' && cigar[0].Length+matchingSuffix == lastRefIndex+1 {
		refIndexToMerge++
	}
	if refIndexToMerge == 0 {
		return false
	}
	g.addEdge(altPath[altIndexToMerge], refPath[refIndexToMerge], 1, false)
	return true
}

func (g *kmerGraph) recoverDanglingTail(vertex *vertexInfo) bool {
	altPath := g.findPathUpwardsToLowestCommonAncestor(vertex)
	if len(altPath) < 5 || g.vertexIsReferenceSource(altPath[0]) {
		return false
	}
	refPath := g.findReferencePath(altPath)
	altBases := getReferenceBasesForPath(altPath)
	refBases := getReferenceBasesForPath(refPath)
	cigar, _ := runSmithWaterman(refBases, altBases, 25, -50, -110, -6, leadingIndel)
	if len(cigar) > 0 && cigar[len(cigar)-1].Operation == 'D' {
		cigar = cigar[:len(cigar)-1]
	}
	if length := len(cigar); length == 0 || length > maxCigarComplexity || cigar[length-1].Operation != 'M' {
		return false
	}
	return g.mergeDanglingTail(altPath, refPath, altBases, refBases, cigar)
}

func (g *kmerGraph) recoverDanglingTails() {
	for _, vertex := range g.getNonReferenceDestinations() {
		g.recoverDanglingTail(vertex)
	}
}

func (g *kmerGraph) findPathDownwardsToHighestCommonDescendant(vertex *vertexInfo) (path []*vertexInfo) {
	currentVertex := vertex
	for !g.vertexIsReferenceNode(currentVertex) && g.vertexOutDegree(currentVertex) == 1 {
		edge := g.outgoingEdges[currentVertex.id][0]
		targetVertex := g.vertices[edge.id]
		if edge.multiplicity < minPruningFactor {
			path = path[:0]
		} else {
			path = append(path, currentVertex)
		}
		currentVertex = targetVertex
	}
	if g.vertexIsReferenceNode(currentVertex) {
		return reversePath(append(path, currentVertex))
	}
	return nil
}

func (g *kmerGraph) findReferencePathUp(path []*vertexInfo) (newPath []*vertexInfo) {
	vertex := path[0]
loop:
	for {
		newPath = append(newPath, vertex)
		for _, edge := range g.incomingEdges[vertex.id] {
			if vertex = g.vertices[edge.id]; g.vertexIsReferenceNode(vertex) {
				continue loop
			}
		}
		return
	}
}

func bestPrefixMatch(sequence1, sequence2 string, maxIndex, kmerSize int32) int32 {
	maxMismatches := maxInt32(1, maxIndex/kmerSize)
	mismatches := int32(0)
	lastGoodIndex := int32(-1)
	for index := int32(0); index < maxIndex; index++ {
		if sequence1[index] != sequence2[index] {
			if mismatches++; mismatches > maxMismatches {
				return -1
			}
			lastGoodIndex = index
		}
	}
	return lastGoodIndex
}

func (g *kmerGraph) extendPathAgainstReference(altPath, refPath []*vertexInfo, nofNodesToExtend, kmerSize int32) ([]*vertexInfo, bool) {
	indexLastDanglingNode := len(altPath) - 1
	indexRefNode := indexLastDanglingNode + int(nofNodesToExtend)
	if indexRefNode >= len(refPath) {
		return altPath, false
	}
	danglingSource := altPath[indexLastDanglingNode]
	altPath = append(altPath[:indexLastDanglingNode], altPath[indexLastDanglingNode+1:]...)
	refSourceSequence := refPath[indexRefNode].bases
	sequenceToExtend := refSourceSequence[:nofNodesToExtend] + danglingSource.bases
	sourceEdge := g.getHeaviestOutgoingEdge(danglingSource)
	sourceTarget := g.vertices[sourceEdge.id]
	g.removeEdge(danglingSource, sourceTarget)
	for i := nofNodesToExtend; i >= 1; i-- {
		newVertex := &vertexInfo{bases: sequenceToExtend[i:minInt32(i+kmerSize, int32(len(sequenceToExtend)))]}
		g.addVertex(newVertex)
		g.addEdge(newVertex, sourceTarget, sourceEdge.multiplicity, false)
		altPath = append(altPath, newVertex)
		sourceTarget = newVertex
	}
	return altPath, true
}

func (g *kmerGraph) mergeDanglingHead(altPath, refPath []*vertexInfo, altBases, refBases string, cigar []sam.CigarOperation, kmerSize int32) bool {
	indexToMerge := bestPrefixMatch(refBases, altBases, cigar[0].Length, kmerSize)
	if indexToMerge <= 0 || indexToMerge >= int32(len(refPath)-1) {
		return false
	}
	if indexToMerge >= int32(len(altPath)) {
		newAltPath, extended := g.extendPathAgainstReference(altPath, refPath, indexToMerge-int32(len(altPath))+2, kmerSize)
		if !extended {
			return false
		}
		altPath = newAltPath
	}
	g.addEdge(refPath[indexToMerge+1], altPath[indexToMerge], 1, false)
	return true
}

func (g *kmerGraph) recoverDanglingHead(vertex *vertexInfo, kmerSize int32) bool {
	altPath := g.findPathDownwardsToHighestCommonDescendant(vertex)
	if len(altPath) < 5 || g.vertexIsReferenceSink(altPath[0]) {
		return false
	}
	refPath := g.findReferencePathUp(altPath)
	altBases := g.getReferenceBasesForPathWithExpandedSources(altPath)
	refBases := g.getReferenceBasesForPathWithExpandedSources(refPath)
	cigar, _ := runSmithWaterman(refBases, altBases, 25, -50, -110, -6, leadingIndel)
	if len(cigar) > 0 && cigar[len(cigar)-1].Operation == 'D' {
		cigar = cigar[:len(cigar)-1]
	}
	if len(cigar) == 0 || len(cigar) > maxCigarComplexity || cigar[0].Operation != 'M' {
		return false
	}
	return g.mergeDanglingHead(altPath, refPath, altBases, refBases, cigar, kmerSize)
}

func (g *kmerGraph) recoverDanglingHeads(kmerSize int32) {
	for _, vertex := range g.getNonReferenceStarts() {
		g.recoverDanglingHead(vertex, kmerSize)
	}
}

func (g *kmerGraph) removePathsNotConnectedToReference() {
	referenceSourceVertex := g.getReferenceSourceVertex()
	visitedVerticesFromReferenceSource := make(map[int32]bool)
	visitedVerticesFromReferenceSink := make(map[int32]bool)
	verticesToVisit := []*vertexInfo{referenceSourceVertex}
	for len(verticesToVisit) > 0 {
		last := len(verticesToVisit) - 1
		vertex := verticesToVisit[last]
		verticesToVisit = verticesToVisit[:last]
		if visitedVerticesFromReferenceSource[vertex.id] {
			continue
		}
		visitedVerticesFromReferenceSource[vertex.id] = true
		for _, edge := range g.outgoingEdges[vertex.id] {
			verticesToVisit = append(verticesToVisit, g.vertices[edge.id])
		}
	}
	referenceSinkVertex := g.getReferenceSinkVertex()
	verticesToVisit = append(verticesToVisit, referenceSinkVertex)
	for len(verticesToVisit) > 0 {
		last := len(verticesToVisit) - 1
		vertex := verticesToVisit[last]
		verticesToVisit = verticesToVisit[:last]
		if visitedVerticesFromReferenceSink[vertex.id] {
			continue
		}
		visitedVerticesFromReferenceSink[vertex.id] = true
		for _, edge := range g.incomingEdges[vertex.id] {
			verticesToVisit = append(verticesToVisit, g.vertices[edge.id])
		}
	}
	for _, vertex := range g.getAllVertices() {
		if !(visitedVerticesFromReferenceSource[vertex.id] && visitedVerticesFromReferenceSink[vertex.id]) {
			g.removeVertex(vertex)
		}
	}
}

func (g *kmerGraph) extendChainByOne(refSource string, vertex *vertexInfo, kmersSeq kmer, i int32) *vertexInfo {
	nextPos := i + g.kmerSize - 1
	lastBase := kmersSeq.bases[nextPos]
	for _, edge := range g.outgoingEdges[vertex.id] {
		connectingVertex := g.vertices[edge.id]
		if lastBase == connectingVertex.getSuffix() {
			edge.multiplicity += 1
			g.incrementIncomingEdgeMultiplicity(vertex, connectingVertex)
			return connectingVertex
		}
	}
	newVertex := g.getKmerVertexButNotRefSource(refSource, kmersSeq, i)
	g.addEdge(vertex, newVertex, 1, kmersSeq.isRef)
	return newVertex
}

func (g *kmerGraph) convertToSequenceGraph() {
	for _, vertex := range g.getAllVertices() {
		if !g.isSourceVertex(vertex) {
			vertex.bases = string(vertex.getSuffix())
		}
	}
}

func (g *kmerGraph) removeNonReferenceComponents() {
	referenceSourceVertex := g.getReferenceSourceVertex()
	visitedVertices := make(map[int32]bool)
	verticesToVisit := []*vertexInfo{referenceSourceVertex}
	for len(verticesToVisit) > 0 {
		last := len(verticesToVisit) - 1
		vertex := verticesToVisit[last]
		verticesToVisit = verticesToVisit[:last]
		if visitedVertices[vertex.id] {
			continue
		}
		visitedVertices[vertex.id] = true
		for _, edge := range g.incomingEdges[vertex.id] {
			verticesToVisit = append(verticesToVisit, g.vertices[edge.id])
		}
		for _, edge := range g.outgoingEdges[vertex.id] {
			verticesToVisit = append(verticesToVisit, g.vertices[edge.id])
		}
	}
	for _, vertex := range g.getAllVertices() {
		if !visitedVertices[vertex.id] {
			g.removeVertex(vertex)
		}
	}
}

func commonMaximumPrefixLength(vertices []*vertexInfo, min int) (length int) {
	for i := 0; i < min; i++ {
		char := vertices[0].bases[i]
		for _, seq := range vertices[1:] {
			if seq.bases[i] != char {
				return
			}
		}
		length++
	}
	return
}

func commonMaximumSuffixLength(vertices []*vertexInfo, min int) (length int) {
	for i := 1; i <= min; i++ {
		bases0 := vertices[0].bases
		char := bases0[len(bases0)-i]
		for _, seq := range vertices[1:] {
			if seq.bases[len(seq.bases)-i] != char {
				return
			}
		}
		length++
	}
	return
}

func commonPrefixAndSuffixOfVertices(vertices []*vertexInfo) (string, string) {
	min := math.MaxInt64
	for _, vertex := range vertices {
		if len(vertex.bases) < min {
			min = len(vertex.bases)
		}
	}
	prefixLength := commonMaximumPrefixLength(vertices, min)
	suffixLength := commonMaximumSuffixLength(vertices, min-prefixLength)
	seq := vertices[0].bases
	return seq[:prefixLength], seq[len(seq)-suffixLength:]
}

func sequenceWithoutPrefixAndSuffix(seq string, prefixLength, suffixLength int) string {
	if len(seq)-prefixLength-suffixLength <= 0 {
		return ""
	}
	return seq[prefixLength : len(seq)-suffixLength]
}

func (g *kmerGraph) mergeDiamondSequences(top, bottom *vertexInfo, middles []*vertexInfo) bool {
	prefix, suffix := commonPrefixAndSuffixOfVertices(middles)
	if len(prefix) == 0 && len(suffix) == 0 {
		return false
	}
	for _, middle := range middles {
		g.updateVertexId(middle)
	}
	prefixVertex := top
	if len(prefix) > 0 {
		prefixVertex = &vertexInfo{bases: prefix}
		g.addVertex(prefixVertex)
		var anyRef bool
		for _, outgoingEdge := range g.outgoingEdges[top.id] {
			if outgoingEdge.isRef {
				anyRef = true
			}
			for _, incomingEdge := range g.incomingEdges[outgoingEdge.id] {
				if incomingEdge.id == top.id {
					incomingEdge.id = prefixVertex.id
				}
			}
		}
		g.outgoingEdges[prefixVertex.id] = g.outgoingEdges[top.id]
		delete(g.outgoingEdges, top.id)
		g.addEdge(top, prefixVertex, 1, anyRef)
	}
	suffixVertex := bottom
	if len(suffix) > 0 {
		suffixVertex = &vertexInfo{bases: suffix}
		g.addVertex(suffixVertex)
		var anyRef bool
		for _, incomingEdge := range g.incomingEdges[bottom.id] {
			if incomingEdge.isRef {
				anyRef = true
			}
			for _, outgoingEdge := range g.outgoingEdges[incomingEdge.id] {
				if outgoingEdge.id == bottom.id {
					outgoingEdge.id = suffixVertex.id
				}
			}
		}
		g.incomingEdges[suffixVertex.id] = g.incomingEdges[bottom.id]
		delete(g.incomingEdges, bottom.id)
		g.addEdge(suffixVertex, bottom, 1, anyRef)
	}

	var newIncomingEdges, newOutgoingEdges []*edgeInfo
	var directIncoming, directOutgoing *edgeInfo

	for _, edge := range g.outgoingEdges[prefixVertex.id] {
		middle := g.vertices[edge.id]
		if remainingSeq := sequenceWithoutPrefixAndSuffix(middle.bases, len(prefix), len(suffix)); remainingSeq != "" {
			middle.bases = remainingSeq
			newOutgoingEdges = append(newOutgoingEdges, edge)
		} else {
			incoming := g.incomingEdges[middle.id][0]
			outgoing := g.outgoingEdges[middle.id][0]
			inAndOutMultiplicity := incoming.multiplicity + outgoing.multiplicity
			inAndOutIsRef := incoming.isRef || outgoing.isRef
			if directOutgoing == nil {
				directIncoming, directOutgoing = g.addEdge(prefixVertex, suffixVertex, inAndOutMultiplicity, inAndOutIsRef)
				newIncomingEdges = append(newIncomingEdges, directIncoming)
				newOutgoingEdges = append(newOutgoingEdges, directOutgoing)
			} else {
				directIncoming.multiplicity += inAndOutMultiplicity
				directOutgoing.multiplicity += inAndOutMultiplicity
				if inAndOutIsRef {
					directIncoming.isRef = true
					directOutgoing.isRef = true
				}
			}
			delete(g.incomingEdges, middle.id)
			delete(g.outgoingEdges, middle.id)
			delete(g.vertices, middle.id)
			middle.id = -1
		}
	}
	for _, edge := range g.incomingEdges[suffixVertex.id] {
		if middle := g.vertices[edge.id]; middle != nil && middle.id != prefixVertex.id {
			newIncomingEdges = append(newIncomingEdges, edge)
		}
	}
	g.outgoingEdges[prefixVertex.id] = newOutgoingEdges
	g.incomingEdges[suffixVertex.id] = newIncomingEdges
	return true
}

func (g *kmerGraph) mergeDiamond(vertex *vertexInfo) bool {
	var middles []*vertexInfo
	for _, edge := range g.outgoingEdges[vertex.id] {
		middles = append(middles, g.vertices[edge.id])
	}
	if len(middles) <= 1 {
		return false
	}
	var bottom *vertexInfo
	for _, middle := range middles {
		if g.vertexOutDegree(middle) < 1 || g.vertexInDegree(middle) != 1 {
			return false
		}
		for _, edge := range g.outgoingEdges[middle.id] {
			target := g.vertices[edge.id]
			if bottom == nil {
				bottom = target
			} else if bottom != target {
				return false
			}
		}
	}
	if len(g.incomingEdges[bottom.id]) != len(middles) {
		return false
	}
	return g.mergeDiamondSequences(vertex, bottom, middles)
}

func (g *kmerGraph) mergeDiamonds() (merged bool) {
	for foundNodes := true; foundNodes; {
		for _, vertex := range g.getAllVertices() {
			if vertex.id != -1 {
				if foundNodes = g.mergeDiamond(vertex); foundNodes {
					merged = true
					break
				}
			}
		}
	}
	return
}

func (g *kmerGraph) mergeTailSequences(top *vertexInfo, tails []*vertexInfo) bool {
	prefix, suffix := commonPrefixAndSuffixOfVertices(tails)
	if len(suffix) < 10 {
		return false
	}
	prefixVertex := top
	if len(prefix) > 0 {
		prefixVertex = &vertexInfo{bases: prefix}
		g.addVertex(prefixVertex)
		var anyRef bool
		for _, outgoingEdge := range g.outgoingEdges[top.id] {
			if outgoingEdge.isRef {
				anyRef = true
			}
			for _, incomingEdge := range g.incomingEdges[outgoingEdge.id] {
				if incomingEdge.id == top.id {
					incomingEdge.id = prefixVertex.id
				}
			}
		}
		g.outgoingEdges[prefixVertex.id] = g.outgoingEdges[top.id]
		delete(g.outgoingEdges, top.id)
		g.addEdge(top, prefixVertex, 1, anyRef)
	}
	suffixVertex := &vertexInfo{bases: suffix}
	g.addVertex(suffixVertex)
	multiplicity := int32(0)
	var anyRef bool
	for _, tail := range tails {
		if remainingSeq := sequenceWithoutPrefixAndSuffix(tail.bases, len(prefix), len(suffix)); remainingSeq != "" {
			tail.bases = remainingSeq
		} else {
			incoming := g.incomingEdges[tail.id][0]
			if incoming.isRef {
				anyRef = true
			}
			multiplicity += incoming.multiplicity
			g.removeEdge(prefixVertex, tail)
		}
	}
	if multiplicity > 0 {
		g.addEdge(prefixVertex, suffixVertex, multiplicity, anyRef)
	}
	return true
}

func (g *kmerGraph) mergeTail(vertex *vertexInfo) bool {
	var tails []*vertexInfo
	for _, edge := range g.outgoingEdges[vertex.id] {
		tails = append(tails, g.vertices[edge.id])
	}
	if len(tails) <= 1 {
		return false
	}
	for _, tail := range tails {
		if g.vertexOutDegree(tail) != 0 || g.vertexInDegree(tail) > 1 {
			return false
		}
	}
	return g.mergeTailSequences(vertex, tails)
}

func (g *kmerGraph) mergeTails() (merged bool) {
	for foundNodes := true; foundNodes; {
		for _, vertex := range g.getAllVertices() {
			if vertex.id != -1 {
				if foundNodes = g.mergeTail(vertex); foundNodes {
					merged = true
					break
				}
			}
		}
	}
	return
}

func (g *kmerGraph) safeToSplit(bottom *vertexInfo, tops []*vertexInfo) bool {
	bottomConnections := g.outgoingEdges[bottom.id]
	for _, top := range tops {
		if top.id == bottom.id {
			return false
		}
		if middleEdges := g.outgoingEdges[top.id]; len(middleEdges) != 1 || middleEdges[0].id != bottom.id {
			return false
		}
		for _, edge := range bottomConnections {
			if edge.id == top.id {
				return false
			}
		}
	}
	return true
}

func commonSuffixOfVertices(vertices []*vertexInfo) string {
	min := math.MaxInt64
	for _, vertex := range vertices {
		if len(vertex.bases) < min {
			min = len(vertex.bases)
		}
	}
	seq := vertices[0].bases
	return seq[len(seq)-commonMaximumSuffixLength(vertices, min):]
}

func eliminatesReferenceSource(suffix string, referenceSource *vertexInfo) bool {
	return referenceSource != nil && len(referenceSource.bases) == len(suffix)
}

func allVerticesAreTheCommonSuffix(suffix string, tops []*vertexInfo) bool {
	for _, top := range tops {
		if len(top.bases) != len(suffix) {
			return false
		}
	}
	return true
}

func (g *kmerGraph) commonSuffix(vertex *vertexInfo, tops []*vertexInfo) (string, bool) {
	if len(tops) < 2 || !g.safeToSplit(vertex, tops) {
		return "", false
	}
	suffix := commonSuffixOfVertices(tops)
	if suffix == "" {
		return "", false
	}
	var referenceSource *vertexInfo
	for _, top := range tops {
		if g.vertexIsReferenceSource(top) {
			referenceSource = top
			break
		}
	}
	if eliminatesReferenceSource(suffix, referenceSource) || allVerticesAreTheCommonSuffix(suffix, tops) {
		return "", false
	}
	return suffix, true
}

func sequenceWithoutSuffix(seq string, suffixLength int) string {
	length := len(seq) - suffixLength
	if length < 0 {
		return ""
	}
	return seq[:length]
}

func (g *kmerGraph) splitCommonSuffixesOfVertex(vertex *vertexInfo) bool {
	incomingEdges := g.incomingEdges[vertex.id]
	var tops []*vertexInfo
	for _, edge := range incomingEdges {
		tops = append(tops, g.vertices[edge.id])
	}
	suffix, ok := g.commonSuffix(vertex, tops)
	if !ok {
		return false
	}
	topEdges := make([]*edgeInfo, len(incomingEdges))
	copy(topEdges, incomingEdges)
	for _, topEdge := range topEdges {
		top := g.vertices[topEdge.id]
		out := g.outgoingEdges[top.id][0]
		topMultiplicity := topEdge.multiplicity
		suffixVertex := &vertexInfo{bases: suffix}
		g.addVertex(suffixVertex)
		remainingSeq := sequenceWithoutSuffix(top.bases, len(suffix))
		var targetVertex *vertexInfo
		if remainingSeq == "" {
			targetVertex = suffixVertex
		} else {
			targetVertex = &vertexInfo{bases: remainingSeq}
			g.addVertex(targetVertex)
			g.addEdge(targetVertex, suffixVertex, 1, out.isRef)
		}
		g.addEdge(suffixVertex, vertex, topMultiplicity, out.isRef)
		topConnectingEdges := g.incomingEdges[top.id]
		for _, edge := range topConnectingEdges {
			g.addEdge(g.vertices[edge.id], targetVertex, edge.multiplicity, edge.isRef)
		}
		g.removeVertex(top)
	}
	return true
}

func (g *kmerGraph) splitCommonSuffixes() (split bool) {
	alreadySplit := make(map[*vertexInfo]bool)
	for foundNodes := true; foundNodes; {
		for _, vertex := range g.getAllVertices() {
			if vertex.id != -1 {
				if !alreadySplit[vertex] {
					alreadySplit[vertex] = true
					if foundNodes = g.splitCommonSuffixesOfVertex(vertex); foundNodes {
						split = true
						break
					}
				}
			}
		}
	}
	return
}

func (g *kmerGraph) mergeCommonSequences(bottom *vertexInfo) bool {
	var tops []*vertexInfo
	for _, edge := range g.incomingEdges[bottom.id] {
		tops = append(tops, g.vertices[edge.id])
	}
	if len(tops) == 0 {
		return false
	}
	seq := tops[0].bases
	for _, top := range tops {
		if top.bases != seq ||
			g.vertexOutDegree(top) != 1 ||
			g.vertexInDegree(top) == 0 ||
			g.outgoingEdges[top.id][0].id != bottom.id {
			return false
		}
	}
	bottom.bases = seq + bottom.bases
	g.updateVertexId(bottom)
	for _, top := range tops {
		for _, edge := range g.incomingEdges[top.id] {
			g.addEdge(g.vertices[edge.id], bottom, edge.multiplicity, edge.isRef)
		}
		g.removeVertex(top)
	}
	for _, edge := range g.outgoingEdges[bottom.id] {
		incomingEdges := g.incomingEdges[edge.id]
		for i, inEdge := range incomingEdges {
			if inEdge.id == bottom.id {
				g.incomingEdges[edge.id] = append(append(incomingEdges[:i], incomingEdges[i+1:]...), inEdge)
				break
			}
		}
	}
	return true
}

func (g *kmerGraph) mergeCommonIncomingSequences() (merged bool) {
	for foundNodes := true; foundNodes; {
		for _, vertex := range g.getAllVertices() {
			if vertex.id != -1 {
				if foundNodes = g.mergeCommonSequences(vertex); foundNodes {
					merged = true
					break
				}
			}
		}
	}
	return
}

func (g *kmerGraph) simplifyOnce() (modified bool) {
	if g.mergeDiamonds() {
		modified = true
	}
	if g.mergeTails() {
		modified = true
	}
	if g.splitCommonSuffixes() {
		modified = true
	}
	if g.mergeCommonIncomingSequences() {
		modified = true
	}
	if g.mergeLinearChains() {
		modified = true
	}
	return
}

func (g *kmerGraph) simplify() {
	g.mergeLinearChains()
	for i := 0; i <= 6; i++ {
		if !g.simplifyOnce() {
			return
		}
	}
	previousGraph := g.partialCopy()
	for i := 7; i <= 100; i++ {
		if !g.simplifyOnce() {
			return
		}
		currentGraph := g.partialCopy()
		if previousGraph.partialEqual(currentGraph) {
			return
		}
		previousGraph = currentGraph
	}
}

func (g *kmerGraph) cleanSequenceGraph() {
	g.mergeLinearChains()
	g.removeNonReferenceComponents()
	g.simplify()
	g.removeNonReferenceComponents()
	g.simplify()
	if len(g.vertices) == 1 {
		for _, vertex := range g.vertices {
			dummy := new(vertexInfo)
			g.addVertex(dummy)
			g.addEdge(vertex, dummy, 0, true)
			break
		}
	}
}

type (
	haplotypePath struct {
		vertices []*vertexInfo
		score    float64
		isRef    bool
	}

	priorityQueue []*haplotypePath
)

func (pq priorityQueue) siftUp(k int, x *haplotypePath) {
	for k > 0 {
		parent := (k - 1) >> 1
		e := pq[parent]
		if x.score <= e.score {
			break
		}
		pq[k] = e
		k = parent
	}
	pq[k] = x
}

func (pq *priorityQueue) enqueue(path *haplotypePath) {
	if len(*pq) == 0 {
		*pq = append(*pq, path)
		return
	}
	*pq = append(*pq, nil)
	pq.siftUp(len(*pq)-1, path)
}

func (pq priorityQueue) siftDown(k int, x *haplotypePath) {
	half := len(pq) >> 1
	for k < half {
		child := (k << 1) + 1
		c := pq[child]
		right := child + 1
		if right < len(pq) && c.score < pq[right].score {
			child = right
			c = pq[child]
		}
		if x.score >= c.score {
			break
		}
		pq[k] = c
		k = child
	}
	pq[k] = x
}

func (pq *priorityQueue) dequeue() *haplotypePath {
	s := len(*pq) - 1
	result := (*pq)[0]
	x := (*pq)[s]
	*pq = (*pq)[:s]
	if s != 0 {
		pq.siftDown(0, x)
	}
	return result
}

func (path *haplotypePath) vertexToExtend() *vertexInfo {
	return path.vertices[len(path.vertices)-1]
}

func (path *haplotypePath) extend(toVertex *vertexInfo, edge *edgeInfo, totalMultiplicityLog10 float64) *haplotypePath {
	extendedScore := path.score + log10(float64(edge.multiplicity)) - totalMultiplicityLog10
	return &haplotypePath{
		vertices: append(path.vertices[:len(path.vertices):len(path.vertices)], toVertex),
		score:    extendedScore,
		isRef:    false,
	}
}

func (path *haplotypePath) convertToHaplotype() *haplotype {
	var bases strings.Builder
	for _, vertex := range path.vertices {
		bases.WriteString(vertex.bases)
	}
	return &haplotype{
		bases: bases.String(),
		score: path.score,
		isRef: path.isRef,
	}
}

type assemblyResultSet []*haplotype

func (set assemblyResultSet) contains(h *haplotype) bool {
	for _, r := range set {
		if h.bases == r.bases {
			return true
		}
	}
	return false
}

func pathTooDivergent(cigar []sam.CigarOperation) bool {
	for _, op := range cigar {
		if op.Operation == 'N' {
			return true
		}
	}
	return false
}

const (
	maxHaplotypes               = 128
	minHaplotypeReferenceLength = 30
)

func (g *kmerGraph) addBestHaplotypes(returnHaplotypes assemblyResultSet, referenceHaplotype *haplotype, paddedReferenceBases string, regionStart int32) assemblyResultSet {
	source := g.getReferenceSourceVertex()
	sink := g.getReferenceSinkVertex()
	var pq priorityQueue
	vertexCounts := make(map[*vertexInfo]int32)
	pq.enqueue(&haplotypePath{vertices: []*vertexInfo{source}})
	var results []*haplotypePath
	for len(pq) > 0 && len(results) < maxHaplotypes {
		pathToExtend := pq.dequeue()
		vertexToExtend := pathToExtend.vertexToExtend()
		if vertexToExtend == sink {
			results = append(results, pathToExtend)
			continue
		}
		count := vertexCounts[vertexToExtend]
		vertexCounts[vertexToExtend] = count + 1
		if count >= maxHaplotypes {
			continue
		}
		var totalMultiplicity int32
		edges := g.outgoingEdges[vertexToExtend.id]
		for _, edge := range edges {
			totalMultiplicity += edge.multiplicity
		}
		totalMultiplicityLog10 := log10(float64(totalMultiplicity))
		for _, edge := range edges {
			pq.enqueue(pathToExtend.extend(g.vertices[edge.id], edge, totalMultiplicityLog10))
		}
	}
	for _, path := range results {
		haplotype := path.convertToHaplotype()
		if returnHaplotypes.contains(haplotype) {
			continue
		}
		cigar := calculateCigar(referenceHaplotype.bases, haplotype.bases, paddedReferenceBases, softclip)
		if len(cigar) == 0 || pathTooDivergent(cigar) || sam.ReferenceLengthFromCigar(cigar) < minHaplotypeReferenceLength {
			continue
		}
		haplotype.cigar = cigar
		haplotype.location = regionStart
		returnHaplotypes = append(returnHaplotypes, haplotype)
	}
	return returnHaplotypes
}

func (region *assemblyRegion) referenceBases() string {
	referenceStart := region.paddedStart() - 1
	referenceEnd := region.paddedEnd()
	referenceLength := referenceEnd - referenceStart
	referenceBytes := make([]byte, referenceLength)
	copy(referenceBytes, region.reference[referenceStart:referenceEnd])
	for i := range referenceBytes {
		referenceBytes[i] = fasta.ToUpperAndN(referenceBytes[i])
	}
	return string(referenceBytes)
}

func (g *kmerGraph) isLowComplexity() bool {
	return len(g.nonUniqueKmers)*4 > len(g.uniqueKmers)
}

func (hc *HaplotypeCaller) assembleReads(region *assemblyRegion) (result assemblyResultSet) {
	referenceBases := region.referenceBases()
	paddedReferenceBases := swPad + referenceBases + swPad
	paddedStart := region.paddedStart()
	referenceHaplotype := makeReferenceHaplotype(referenceBases, paddedStart)
	result = append(result, referenceHaplotype)

	processKmerSize := func(kmerSize int32, lastAttempt bool) bool {
		if !lastAttempt && nonUniqueKmersExist(referenceBases, kmerSize) {
			return false
		}
		graph := newKmerGraph(kmerSize)
		seqs := []kmer{makeSequenceForKmersFromReference(referenceBases, kmerSize)}
		refSource := referenceBases[:kmerSize]
		for _, aln := range region.alns {
			seqs = hc.addSequencesForKmers(seqs, aln, kmerSize)
		}
		graph.initializeNonUniqueKmers(seqs, kmerSize)
		for _, kmersSeq := range seqs {
			uniqueStart := graph.findStartOfKmers(kmersSeq)
			if uniqueStart == -1 {
				continue
			}
			vertex := graph.getKmerVertex(kmersSeq, uniqueStart)
			graph.increaseCountsMatchedKmers(kmersSeq, vertex)
			for i, end := uniqueStart+1, kmersSeq.stop-kmerSize; i <= end; i++ {
				vertex = graph.extendChainByOne(refSource, vertex, kmersSeq, i)
			}
		}
		if len(graph.vertices) == 0 {
			return false
		}

		graph.pruneChainsWithLowWeight()

		if graph.hasCycle() {
			return false
		}
		if !lastAttempt && graph.isLowComplexity() {
			return false
		}
		graph.recoverDanglingTails()
		graph.recoverDanglingHeads(kmerSize)
		graph.removePathsNotConnectedToReference()
		graph.convertToSequenceGraph()
		graph.cleanSequenceGraph()
		result = graph.addBestHaplotypes(result, referenceHaplotype, paddedReferenceBases, paddedStart)
		return true
	}

	var graphSeen bool

	for _, kmerSize := range []int32{10, 25} {
		if len(referenceBases) < int(kmerSize) {
			return result
		}
		if processKmerSize(kmerSize, false) {
			graphSeen = true
		}
	}

	if graphSeen {
		return result
	}

	kmerSize := int32(35)

	for attempt := 1; attempt < 6; attempt++ {
		if len(referenceBases) < int(kmerSize) {
			return result
		}
		if processKmerSize(kmerSize, false) {
			return result
		}
		kmerSize += 10
	}

	if len(referenceBases) < int(kmerSize) {
		return result
	}
	processKmerSize(kmerSize, true)
	return result
}
