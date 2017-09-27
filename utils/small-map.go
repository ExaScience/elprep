package utils

type (
	SmallMapEntry struct {
		Key   Symbol
		Value interface{}
	}

	SmallMap []SmallMapEntry
)

func (m SmallMap) Get(key Symbol) (interface{}, bool) {
	for _, entry := range m {
		if entry.Key == key {
			return entry.Value, true
		}
	}
	return nil, false
}

func (m *SmallMap) Set(key Symbol, value interface{}) {
	for index := range *m {
		if (*m)[index].Key == key {
			(*m)[index].Value = value
			return
		}
	}
	*m = append(*m, SmallMapEntry{key, value})
}

func (m SmallMap) Delete(key Symbol) (SmallMap, bool) {
	for index, entry := range m {
		if entry.Key == key {
			return append(m[:index], m[index+1:]...), true
		}
	}
	return m, false
}
