package utils

// A StringMap maps strings to strings.
type StringMap map[string]string

/*
Find returns the first index in a slice of StringMap where the
predicate returns true, or -1 if predicate never returns true.
*/
func Find(dict []StringMap, predicate func(record StringMap) bool) int {
	for index, record := range dict {
		if predicate(record) {
			return index
		}
	}
	return -1
}

/*
SetUniqueEntry checks if a mapping for the given key already exists in
the StringMap. If this is the case, it returns false and the StringMap
is not modified.  Otherwise, the given key/value pair is added to the
StringMap.
*/
func (record StringMap) SetUniqueEntry(key, value string) bool {
	if _, found := record[key]; found {
		return false
	}
	record[key] = value
	return true
}
