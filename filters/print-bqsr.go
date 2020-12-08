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

package filters

import (
	"encoding/gob"
	"fmt"
	"io"
	"log"
	"path"
	"sort"
	"strconv"

	"github.com/exascience/elprep/v5/internal"
)

const quantizationLevel = 16

const (
	countString              = "Count"
	covariateNameString      = "CovariateName"
	covariateValueString     = "CovariateValue"
	empiricalQualityString   = "EmpiricalQuality"
	errorsString             = "Errors"
	estimatedQReportedString = "EstimatedQReported"
	eventTypeString          = "EventType"
	observationsString       = "Observations"
	qualityScoreString       = "QualityScore"
	quantizedScoreString     = "QuantizedScore"
	readGroupString          = "ReadGroup"
)

func (recal *BaseRecalibratorTables) printQuantizationTable(file io.Writer) {
	observations, scores := initializeQuantizedQualityScores(recal.QualityScores, quantizationLevel)

	fmt.Fprintf(file, "#:%sTable:3:%d:%%d:%%d:%%d:;\n", internal.BQSRTablenamePrefix, len(observations))
	fmt.Fprintf(file, "#:%sTable:Quantized:Quality quantization map\n", internal.BQSRTablenamePrefix)

	maxLenQualityScore := len(qualityScoreString)
	maxLenCount := len(countString)
	maxLenQuantizedScore := len(quantizedScoreString)

	for i, observation := range observations {
		maxLenQualityScore = maxInt(maxLenQualityScore, len(strconv.FormatInt(int64(i), 10)))
		maxLenCount = maxInt(maxLenCount, len(strconv.FormatInt(int64(observation), 10)))
		maxLenQuantizedScore = maxInt(maxLenQuantizedScore, len(strconv.FormatInt(int64(scores[i]), 10)))
	}

	fmt.Fprintf(file, "%-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenCount, countString)
	fmt.Fprintf(file, "  %-[1]*[2]s\n", maxLenQuantizedScore, quantizedScoreString)

	for i, observation := range observations {
		fmt.Fprintf(file, "%[1]*[2]d", maxLenQualityScore, i)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenCount, observation)
		fmt.Fprintf(file, "  %[1]*[2]d\n", maxLenQuantizedScore, scores[i])
	}

	fmt.Fprintln(file)
}

func (recal *BaseRecalibratorTables) printCombinedBQSRTable(file io.Writer) {
	table := initializeCombinedBQSRTable(recal.QualityScores)

	fmt.Fprintf(file, "#:%sTable:6:%d:%%s:%%s:%%.4f:%%.4f:%%d:%%.2f:;\n", internal.BQSRTablenamePrefix, len(table))
	fmt.Fprintf(file, "#:%sTable:RecalTable0:\n", internal.BQSRTablenamePrefix)

	maxLenReadGroup := len(readGroupString)
	maxLenEventType := len(eventTypeString)
	maxLenEmpiricalQuality := len(empiricalQualityString)
	maxLenEstimatedQReported := len(estimatedQReportedString)
	maxLenObservations := len(observationsString)
	maxLenErrors := len(errorsString)

	var readGroups []string

	for rg, entry := range table {
		readGroups = append(readGroups, rg)
		maxLenReadGroup = maxInt(maxLenReadGroup, len(rg))
		maxLenEmpiricalQuality = maxInt(maxLenEmpiricalQuality, len(strconv.FormatInt(int64(entry.EmpiricalQuality), 10))+5)
		maxLenEstimatedQReported = maxInt(maxLenEstimatedQReported, len(strconv.FormatFloat(entry.reportedQuality, 'f', 4, 64)))
		maxLenObservations = maxInt(maxLenObservations, len(strconv.FormatInt(int64(entry.Observations), 10)))
		maxLenErrors = maxInt(maxLenErrors, len(strconv.FormatInt(int64(entry.Mismatches), 10))+3)
	}

	fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEstimatedQReported, estimatedQReportedString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	fmt.Fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

	sort.Strings(readGroups)

	for _, rg := range readGroups {
		entry := table[rg]
		fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, rg)
		fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		fmt.Fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		fmt.Fprintf(file, "  %[1]*.4[2]f", maxLenEstimatedQReported, entry.reportedQuality)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		fmt.Fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	fmt.Fprintln(file)
}

func (recal *BaseRecalibratorTables) printBQSRTable(file io.Writer) {
	fmt.Fprintf(file, "#:%sTable:6:%d:%%s:%%d:%%s:%%.4f:%%d:%%.2f:;\n", internal.BQSRTablenamePrefix, len(recal.QualityScores))
	fmt.Fprintf(file, "#:%sTable:RecalTable1:\n", internal.BQSRTablenamePrefix)

	maxLenReadGroup := len(readGroupString)
	maxLenQualityScore := len(qualityScoreString)
	maxLenEventType := len(eventTypeString)
	maxLenEmpiricalQuality := len(empiricalQualityString)
	maxLenObservations := len(observationsString)
	maxLenErrors := len(errorsString)

	var keys []bqsrTableKey

	for key, entry := range recal.QualityScores {
		keys = append(keys, key)
		maxLenReadGroup = maxInt(maxLenReadGroup, len(key.ReadGroup))
		maxLenQualityScore = maxInt(maxLenQualityScore, len(strconv.FormatInt(int64(key.Qual), 10)))
		maxLenEmpiricalQuality = maxInt(maxLenEmpiricalQuality, len(strconv.FormatInt(int64(entry.EmpiricalQuality), 10))+5)
		maxLenObservations = maxInt(maxLenObservations, len(strconv.FormatInt(int64(entry.Observations), 10)))
		maxLenErrors = maxInt(maxLenErrors, len(strconv.FormatInt(int64(entry.Mismatches), 10))+3)
	}

	fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	fmt.Fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

	sort.Slice(keys, func(i, j int) bool {
		e1, e2 := keys[i], keys[j]
		if e1.ReadGroup < e2.ReadGroup {
			return true
		}
		if e1.ReadGroup == e2.ReadGroup {
			return e1.Qual < e2.Qual
		}
		return false
	})

	for _, key := range keys {
		entry := recal.QualityScores[key]
		fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, key.ReadGroup)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenQualityScore, key.Qual)
		fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		fmt.Fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		fmt.Fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	fmt.Fprintln(file)
}

type covariateBqsrTableKey struct {
	cycleNotContext bool
	text            string
	bqsrTableKey
}

func (recal *BaseRecalibratorTables) printOtherCovariateTable(file io.Writer) {
	fmt.Fprintf(file, "#:%sTable:8:%d:%%s:%%d:%%s:%%s:%%s:%%.4f:%%d:%%.2f:;\n", internal.BQSRTablenamePrefix, len(recal.Cycles)+len(recal.Contexts))
	fmt.Fprintf(file, "#:%sTable:RecalTable2:\n", internal.BQSRTablenamePrefix)

	maxLenReadGroup := len(readGroupString)
	maxLenQualityScore := len(qualityScoreString)
	maxLenCovariateValue := len(covariateValueString)
	maxLenCovariateName := maxInt(maxInt(len(covariateNameString), len("Cycle")), len("Context"))
	maxLenEventType := len(eventTypeString)
	maxLenEmpiricalQuality := len(empiricalQualityString)
	maxLenObservations := len(observationsString)
	maxLenErrors := len(errorsString)

	var keys []covariateBqsrTableKey

	for key, entry := range recal.Cycles {
		text := strconv.FormatInt(int64(key.Covariate), 10)
		keys = append(keys, covariateBqsrTableKey{cycleNotContext: true, text: text, bqsrTableKey: key})
		maxLenReadGroup = maxInt(maxLenReadGroup, len(key.ReadGroup))
		maxLenQualityScore = maxInt(maxLenQualityScore, len(strconv.FormatInt(int64(key.Qual), 10)))
		maxLenCovariateValue = maxInt(maxLenCovariateValue, len(text))
		maxLenEmpiricalQuality = maxInt(maxLenEmpiricalQuality, len(strconv.FormatInt(int64(entry.EmpiricalQuality), 10))+5)
		maxLenObservations = maxInt(maxLenObservations, len(strconv.FormatInt(int64(entry.Observations), 10)))
		maxLenErrors = maxInt(maxLenErrors, len(strconv.FormatInt(int64(entry.Mismatches), 10))+3)
	}

	for key, entry := range recal.Contexts {
		text := keyToString(key.Covariate)
		keys = append(keys, covariateBqsrTableKey{cycleNotContext: false, text: text, bqsrTableKey: key})
		maxLenReadGroup = maxInt(maxLenReadGroup, len(key.ReadGroup))
		maxLenQualityScore = maxInt(maxLenQualityScore, len(strconv.FormatInt(int64(key.Qual), 10)))
		maxLenCovariateValue = maxInt(maxLenCovariateValue, len(text))
		maxLenEmpiricalQuality = maxInt(maxLenEmpiricalQuality, len(strconv.FormatInt(int64(entry.EmpiricalQuality), 10))+5)
		maxLenObservations = maxInt(maxLenObservations, len(strconv.FormatInt(int64(entry.Observations), 10)))
		maxLenErrors = maxInt(maxLenErrors, len(strconv.FormatInt(int64(entry.Mismatches), 10))+3)
	}

	fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenCovariateValue, covariateValueString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenCovariateName, covariateNameString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	fmt.Fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	fmt.Fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

	sort.Slice(keys, func(i, j int) bool {
		e1, e2 := keys[i], keys[j]
		if e1.ReadGroup < e2.ReadGroup {
			return true
		}
		if e1.ReadGroup == e2.ReadGroup {
			if e1.Qual < e2.Qual {
				return true
			}
			if e1.Qual == e2.Qual {
				return e1.text < e2.text
			}
		}
		return false
	})

	for _, key := range keys {
		var entry *bqsrEntry
		var name string
		if key.cycleNotContext {
			entry = recal.Cycles[key.bqsrTableKey]
			name = "Cycle"
		} else {
			entry = recal.Contexts[key.bqsrTableKey]
			name = "Context"
		}
		fmt.Fprintf(file, "%-[1]*[2]s", maxLenReadGroup, key.ReadGroup)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenQualityScore, key.Qual)
		fmt.Fprintf(file, "  %-[1]*[2]s", maxLenCovariateValue, key.text)
		fmt.Fprintf(file, "  %-[1]*[2]s", maxLenCovariateName, name)
		fmt.Fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		fmt.Fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		fmt.Fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		fmt.Fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	fmt.Fprintln(file)
}

// PrintBQSRTables creates a recalibration report file.
func (recal *BaseRecalibratorTables) PrintBQSRTables(name string) {
	file := internal.FileCreate(name)
	defer internal.Close(file)
	fmt.Fprintf(file, "#:%sReport.v1.1:5\n", internal.BQSRTablenamePrefix)
	fmt.Fprintf(file, "#:%sTable:2:17:%%s:%%s:;\n", internal.BQSRTablenamePrefix)
	fmt.Fprintf(file, "#:%sTable:Arguments:Recalibration argument collection values used in this run\n", internal.BQSRTablenamePrefix)
	fmt.Fprintln(file, "Argument                    Value                                                                   ")
	fmt.Fprintln(file, "binary_tag_name             null                                                                    ")
	fmt.Fprintln(file, "covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate")
	fmt.Fprintln(file, "default_platform            null                                                                    ")
	fmt.Fprintln(file, "deletions_default_quality   45                                                                      ")
	fmt.Fprintln(file, "force_platform              null                                                                    ")
	fmt.Fprintln(file, "indels_context_size         3                                                                       ")
	fmt.Fprintln(file, "insertions_default_quality  45                                                                      ")
	fmt.Fprintln(file, "low_quality_tail            2                                                                       ")
	fmt.Fprintln(file, "maximum_cycle_value         500                                                                     ")
	fmt.Fprintln(file, "mismatches_context_size     2                                                                       ")
	fmt.Fprintln(file, "mismatches_default_quality  -1                                                                      ")
	fmt.Fprintln(file, "no_standard_covs            false                                                                   ")
	fmt.Fprintln(file, "quantizing_levels           16                                                                      ")
	fmt.Fprintln(file, "recalibration_report        null                                                                    ")
	fmt.Fprintln(file, "run_without_dbsnp           false                                                                   ")
	fmt.Fprintln(file, "solid_nocall_strategy       THROW_EXCEPTION                                                         ")
	fmt.Fprintln(file, "solid_recal_mode            SET_Q_ZERO                                                              ")
	fmt.Fprintln(file)
	recal.printQuantizationTable(file)
	recal.printCombinedBQSRTable(file)
	recal.printBQSRTable(file)
	recal.printOtherCovariateTable(file)
}

// PrintBQSRTablesToIntermediateFile prints the recalibration tables to a gob file.
func (recal *BaseRecalibratorTables) PrintBQSRTablesToIntermediateFile(name string) {
	file := internal.FileCreate(name)
	defer internal.Close(file)
	if err := gob.NewEncoder(file).Encode(recal); err != nil {
		log.Panic(err)
	}
}

// LoadAndCombineBQSRTables loads and merges multiple recalibration tables from file into a single, new recalibration table.
func LoadAndCombineBQSRTables(bqsrPath string) *BaseRecalibratorTables {
	// create bqsr tables
	result := NewBaseRecalibratorTables()
	// go through the files, loading intermediate tables
	bqsrPath, files := internal.Directory(bqsrPath)
	for _, fileName := range files {
		partialResult := BaseRecalibratorTables{}
		file := internal.FileOpen(path.Join(bqsrPath, fileName))
		if err := gob.NewDecoder(file).Decode(&partialResult); err != nil {
			_ = file.Close()
			log.Panic(err)
		}
		internal.Close(file)
		// BqsrTable
		result.QualityScores = result.QualityScores.merge(partialResult.QualityScores)
		result.Contexts = result.Contexts.merge(partialResult.Contexts)
		result.Cycles = result.Cycles.merge(partialResult.Cycles)
	}
	return &result
}
