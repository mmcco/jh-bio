package bioutils

// K is the kmer size (in basepairs).
const K = 31
// M is the minimizer size, which must be <= k.
const M = 15

// KMask and MMask contain k and m consecutive right aligned 1 bits
// respectively (e.g. "0000011111111" for k=8).
const KMask = (1 << (2 * K)) - 1
const MMask = (1 << (2 * M)) - 1
