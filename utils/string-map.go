package utils

type StringMap map[string]string

func Find(dict []StringMap, predicate func(record StringMap) bool) int {
	for index, record := range dict {
		if predicate(record) {
			return index
		}
	}
	return -1
}

func (record StringMap) SetUniqueEntry(key, value string) bool {
	if _, found := record[key]; found {
		return false
	} else {
		record[key] = value
		return true
	}
}
