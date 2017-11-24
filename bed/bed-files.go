package bed

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/exascience/elprep/utils"
)

// Helper function for parsing a track line field.
func splitTrackField(field string) (string, string) {
	split := strings.Split(field, "=")
	return split[0], split[1]
}

// ParseBed parses a BED file. See
// https://genome.ucsc.edu/FAQ/FAQformat.html#format1
func ParseBed(filename string) (b *Bed, err error) {

	bed := NewBed()

	// open file
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer func() {
		if nerr := file.Close(); err == nil {
			err = nerr
		}
	}()

	scanner := bufio.NewScanner(file)

	var track *Track // for storing the current track

	for scanner.Scan() {
		line := scanner.Text()
		data := strings.Split(line, "\t")
		// check if the line is a new track
		if data[0] == "track" {
			// create new track, store the old one
			if track != nil {
				bed.Tracks = append(bed.Tracks, track)
			}
			// all track entries are optional
			// parse and collect those that are used
			fields := make(map[string]string)
			for _, field := range data[1:] {
				key, val := splitTrackField(field)
				fields[key] = val
			}
			track = NewTrack(fields)
		} else {
			// parse a region entry
			chrom := utils.Intern(data[0])
			var err error
			start, err := strconv.Atoi(data[1])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region start: %v ", err.Error())
			}
			end, err := strconv.Atoi(data[2])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region end: %v ", err.Error())
			}
			region, err := NewRegion(chrom, int32(start), int32(end), data[3:])
			if err != nil {
				return nil, fmt.Errorf("invalid bed region: %v ", err.Error())
			}
			AddRegion(bed, region)
			if track != nil {
				track.Regions = append(track.Regions, region)
			}
		}

	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error while reading bed file: %v ", err.Error())
	}
	// Make sure bed regions are sorted.
	sortRegions(bed)
	return bed, nil
}
