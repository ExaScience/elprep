package internal

func BoolHash(b bool) uint64 {
	if b {
		return (1 << 35) - 1
	} else {
		return ((1 << 29) - 1) << 35
	}
}

func StringHash(s string) (hash uint64) {
	// DJBX33A
	hash = 5381
	for _, b := range s {
		hash = ((hash << 5) + hash) + uint64(b)
	}
	return
}
