package internal

// BoolHash returns a hash value for the given boolean value.
func BoolHash(b bool) uint64 {
	if b {
		return (1 << 35) - 1
	}
	return ((1 << 29) - 1) << 35
}

// StringHash returns a hash value for the given string value.
func StringHash(s string) (hash uint64) {
	// DJBX33A
	hash = 5381
	for _, b := range s {
		hash = ((hash << 5) + hash) + uint64(b)
	}
	return
}
