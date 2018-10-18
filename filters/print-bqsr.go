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

package filters

import (
	"encoding/gob"
	"fmt"
	"io"
	"os"
	"path"
	"sort"
	"strconv"

	"github.com/exascience/elprep/v4/internal"
)

const quantizationLevel = 16

func (recal *BaseRecalibratorTables) fprintln(w io.Writer, a ...interface{}) {
	if recal.err != nil {
		return
	}
	fmt.Fprintln(w, a...)
}

func (recal *BaseRecalibratorTables) fprintf(w io.Writer, format string, a ...interface{}) {
	if recal.err != nil {
		return
	}
	fmt.Fprintf(w, format, a...)
}

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

	recal.fprintln(file, "#:GATKTable:3:94:\045d:\045d:\045d:;")
	recal.fprintln(file, "#:GATKTable:Quantized:Quality quantization map")

	maxLenQualityScore := len(qualityScoreString)
	maxLenCount := len(countString)
	maxLenQuantizedScore := len(quantizedScoreString)

	for i, observation := range observations {
		maxLenQualityScore = maxInt(maxLenQualityScore, len(strconv.FormatInt(int64(i), 10)))
		maxLenCount = maxInt(maxLenCount, len(strconv.FormatInt(int64(observation), 10)))
		maxLenQuantizedScore = maxInt(maxLenQuantizedScore, len(strconv.FormatInt(int64(scores[i]), 10)))
	}

	recal.fprintf(file, "%-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenCount, countString)
	recal.fprintf(file, "  %-[1]*[2]s\n", maxLenQuantizedScore, quantizedScoreString)

	for i, observation := range observations {
		recal.fprintf(file, "%[1]*[2]d", maxLenQualityScore, i)
		recal.fprintf(file, "  %[1]*[2]d", maxLenCount, observation)
		recal.fprintf(file, "  %[1]*[2]d\n", maxLenQuantizedScore, scores[i])
	}

	recal.fprintln(file)
}

func (recal *BaseRecalibratorTables) printCombinedBQSRTable(file io.Writer) {
	table := initializeCombinedBQSRTable(recal.QualityScores)

	recal.fprintln(file, "#:GATKTable:6:1:\045s:\045s:\045.4f:\045.4f:\045d:\045.2f:;")
	recal.fprintln(file, "#:GATKTable:RecalTable0:")

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

	recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEstimatedQReported, estimatedQReportedString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	recal.fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

	sort.Strings(readGroups)

	for _, rg := range readGroups {
		entry := table[rg]
		recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, rg)
		recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		recal.fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		recal.fprintf(file, "  %[1]*.4[2]f", maxLenEstimatedQReported, entry.reportedQuality)
		recal.fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		recal.fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	recal.fprintln(file)
}

func (recal *BaseRecalibratorTables) printBQSRTable(file io.Writer) {
	recal.fprintln(file, "#:GATKTable:6:36:\045s:\045d:\045s:\045.4f:\045d:\045.2f:;")
	recal.fprintln(file, "#:GATKTable:RecalTable1:")

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

	recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	recal.fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

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
		recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, key.ReadGroup)
		recal.fprintf(file, "  %[1]*[2]d", maxLenQualityScore, key.Qual)
		recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		recal.fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		recal.fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		recal.fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	recal.fprintln(file)
}

type covariateBqsrTableKey struct {
	cycleNotContext bool
	text            string
	bqsrTableKey
}

func (recal *BaseRecalibratorTables) printOtherCovariateTable(file io.Writer) {
	recal.fprintln(file, "#:GATKTable:8:7756:\045s:\045d:\045s:\045s:\045s:\045.4f:\045d:\045.2f:;")
	recal.fprintln(file, "#:GATKTable:RecalTable2:")

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

	recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, readGroupString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenQualityScore, qualityScoreString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenCovariateValue, covariateValueString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenCovariateName, covariateNameString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, eventTypeString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenEmpiricalQuality, empiricalQualityString)
	recal.fprintf(file, "  %-[1]*[2]s", maxLenObservations, observationsString)
	recal.fprintf(file, "  %-[1]*[2]s\n", maxLenErrors, errorsString)

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
		recal.fprintf(file, "%-[1]*[2]s", maxLenReadGroup, key.ReadGroup)
		recal.fprintf(file, "  %[1]*[2]d", maxLenQualityScore, key.Qual)
		recal.fprintf(file, "  %-[1]*[2]s", maxLenCovariateValue, key.text)
		recal.fprintf(file, "  %-[1]*[2]s", maxLenCovariateName, name)
		recal.fprintf(file, "  %-[1]*[2]s", maxLenEventType, "M")
		recal.fprintf(file, "  %[1]*[2]d.0000", maxLenEmpiricalQuality-5, entry.EmpiricalQuality)
		recal.fprintf(file, "  %[1]*[2]d", maxLenObservations, entry.Observations)
		recal.fprintf(file, "  %[1]*[2]d.00\n", maxLenErrors-3, entry.Mismatches)
	}

	recal.fprintln(file)
}

// PrintBQSRTables creates a recalibration report file.
func (recal *BaseRecalibratorTables) PrintBQSRTables(name string) error {
	file, err := os.Create(name)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()
	recal.fprintln(file, "#:GATKReport.v1.1:5")
	recal.fprintln(file, "#:GATKTable:2:17:\045s:\045s:;")
	recal.fprintln(file, "#:GATKTable:Arguments:Recalibration argument collection values used in this run")
	recal.fprintln(file, "Argument                    Value                                                                   ")
	recal.fprintln(file, "binary_tag_name             null                                                                    ")
	recal.fprintln(file, "covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate")
	recal.fprintln(file, "default_platform            null                                                                    ")
	recal.fprintln(file, "deletions_default_quality   45                                                                      ")
	recal.fprintln(file, "force_platform              null                                                                    ")
	recal.fprintln(file, "indels_context_size         3                                                                       ")
	recal.fprintln(file, "insertions_default_quality  45                                                                      ")
	recal.fprintln(file, "low_quality_tail            2                                                                       ")
	recal.fprintln(file, "maximum_cycle_value         500                                                                     ")
	recal.fprintln(file, "mismatches_context_size     2                                                                       ")
	recal.fprintln(file, "mismatches_default_quality  -1                                                                      ")
	recal.fprintln(file, "no_standard_covs            false                                                                   ")
	recal.fprintln(file, "quantizing_levels           16                                                                      ")
	recal.fprintln(file, "recalibration_report        null                                                                    ")
	recal.fprintln(file, "run_without_dbsnp           false                                                                   ")
	recal.fprintln(file, "solid_nocall_strategy       THROW_EXCEPTION                                                         ")
	recal.fprintln(file, "solid_recal_mode            SET_Q_ZERO                                                              ")
	recal.fprintln(file)
	recal.printQuantizationTable(file)
	recal.printCombinedBQSRTable(file)
	recal.printBQSRTable(file)
	recal.printOtherCovariateTable(file)
	return recal.err
}

// PrintBQSRTablesToIntermediateFile prints the recalibration tables to a gob file.
func (recal *BaseRecalibratorTables) PrintBQSRTablesToIntermediateFile(name string) error {
	file, err := os.Create(name)
	if err != nil {
		return err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()
	return gob.NewEncoder(file).Encode(recal)
}

// LoadAndCombineBQSRTables loads and merges multiple recalibration tables from file into a single, new recalibration table.
func LoadAndCombineBQSRTables(bqsrPath string) (BaseRecalibratorTables, error) {
	// create bqsr tables
	result := NewBaseRecalibratorTables()
	// go through the files, loading intermediate tables
	files, err := internal.Directory(bqsrPath)
	if err != nil {
		return BaseRecalibratorTables{}, err
	}
	for _, fileName := range files {
		partialResult := BaseRecalibratorTables{}
		file, err := os.Open(path.Join(bqsrPath, fileName))
		if err != nil {
			return BaseRecalibratorTables{}, err
		}
		if err = gob.NewDecoder(file).Decode(&partialResult); err != nil {
			_ = file.Close()
			return BaseRecalibratorTables{}, err
		}
		if err = file.Close(); err != nil {
			return BaseRecalibratorTables{}, err
		}
		// BqsrTable
		result.QualityScores = result.QualityScores.merge(partialResult.QualityScores)
		result.Contexts = result.Contexts.merge(partialResult.Contexts)
		result.Cycles = result.Cycles.merge(partialResult.Cycles)
	}
	return result, nil
}
