/*
   General sequence-manipulation functions.
*/

package repeatgenome

import (
	"unsafe"
)

/*
   Converts a TextSeq to the more memory-efficient Seq type.
   Upper- and lower-case base bytes are currently supported, but stable code should immediately convert to lower-case.
   The logic works and is sane, but could be altered in the future for brevity and efficiency.
*/
func GetSeq(textSeq []byte) Seq {
	numBytes := ceilDiv(len(textSeq), 4)
	seq := Seq{
		Bytes: make([]byte, numBytes, numBytes),
		Len:   uint64(len(textSeq)),
	}

	// determines how much to shift the current byte of interest
	// starts at 6, wraps around to 254, which mods to 0
	// therefore loops 6 -> 4 -> 2 -> 0 -> ...
	// see the last line of the for-loop
	var shift uint8 = 6

	for i := 0; i < len(textSeq); i++ {

		switch textSeq[i] {
		case 'a':
			// already zero
			break
		case 'c':
			seq.Bytes[i/4] |= 1 << shift
			break
		case 'g':
			seq.Bytes[i/4] |= 2 << shift
			break
		case 't':
			seq.Bytes[i/4] |= 3 << shift
			break
		case 'A':
			// already zero
			break
		case 'C':
			seq.Bytes[i/4] |= 1 << shift
			break
		case 'G':
			seq.Bytes[i/4] |= 2 << shift
			break
		case 'T':
			seq.Bytes[i/4] |= 3 << shift
			break
		default:
			panic("TextSeq.GetSeq(): byte other than 'a', 'c', 'g', or 't' encountered")
		}

		shift = (shift - 2) % 8
	}

	return seq
}

/*
   Return the subsequence of the supplied Seq from a (inclusive) to b (exclusive), like a slice.
*/
func (seq Seq) Subseq(a, b uint64) Seq {
	numBytes := ceilDiv_U64(b-a, 4)
	subseq := Seq{
		Bytes: make([]byte, numBytes, numBytes),
		Len:   b - a,
	}

	// How much to left-shift the byte currently being altered.
	// For more details, see the identical usage in TextSeq.Seq()
	var shift uint8 = 6

	var i uint64
	for i = 0; i < subseq.Len; i++ {
		subseq.Bytes[i/4] |= seq.GetBase(i+a) << shift
	}
	shift = (shift - 2) % 8
	return subseq
}

/*
   Return the i-th byte of the Seq (zero-indexed).
*/
func (seq Seq) GetBase(i uint64) uint8 {
	thisByte := seq.Bytes[i/4]
	return thisByte >> (6 - 2*(i%4))
}

/*
   A more declarative and modifiable accessor function.
   While it would almost certainly be inlined, this is such a performance-critical operation that this function isn't currently used.
*/
func (kmer Kmer) Int() uint64 {
	return *(*uint64)(unsafe.Pointer(&kmer))
}

/*
   A more declarative and modifiable accessor function.
   While it would almost certainly be inlined, this is such a performance-critical operation that this function isn't currently used.
*/
func (kmer Kmer) ClassID() ClassID {
	return *(*ClassID)(unsafe.Pointer(&kmer[8]))
}

func (kmer *Kmer) SetInt(kmerInt uint64) {
	*(*uint64)(unsafe.Pointer(kmer)) = kmerInt
}

func (kmer *Kmer) SetClassID(classID ClassID) {
	*(*ClassID)(unsafe.Pointer(&kmer[8])) = classID
}
