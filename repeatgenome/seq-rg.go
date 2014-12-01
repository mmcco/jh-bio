/*
   General sequence-manipulation functions.
*/

package repeatgenome

import (
    "fmt"
    "unsafe"
)

func (seq TextSeq) kmerInt() KmerInt {
    if len(seq) < 1 || len(seq) > 31 {
        panic("TextSeq.kmerInt() can only int-ize sequences where 0 < length < 32")
    }
    var kmerInt KmerInt = 0
    for _, c := range seq {
        kmerInt = kmerInt << 2
        switch c {
        case 'a':
            break
        case 'c':
            kmerInt |= 1
            break
        case 'g':
            kmerInt |= 2
            break
        case 't':
            kmerInt |= 3
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to TextSeq.kmerInt(): %c", c))
        }
    }
    return kmerInt
}

func (seq TextSeq) minInt() MinInt {
    if len(seq) < 1 || len(seq) > 15 {
        panic("TextSeq.minInt() can only int-ize sequences where 0 < length < 32")
    }
    var minInt MinInt = 0
    for _, c := range seq {
        minInt = minInt << 2
        switch c {
        case 'a':
            break
        case 'c':
            minInt |= 1
            break
        case 'g':
            minInt |= 2
            break
        case 't':
            minInt |= 3
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to TextSeq.minInt(): %c", c))
        }
    }
    return minInt
}

/*
   Returns the reverse complement of a KmerInt.
*/
func (kmerInt KmerInt) revComp() KmerInt {
    var revComp KmerInt
    var i uint8
    for i = 0; i < k; i++ {
        revComp <<= 2
        revComp |= 3 &^ kmerInt // BEWARE: &^ is not commutative
        kmerInt >>= 2
    }
    return revComp
}

/*
   Returns the reverse complement of a MinInt.
*/
func (minInt MinInt) revComp() MinInt {
    var revComp MinInt
    var i uint8
    for i = 0; i < m; i++ {
        revComp <<= 2
        revComp |= 3 &^ minInt // BEWARE: &^ is not commutative
        minInt >>= 2
    }
    return revComp
}

/*
   Return the minimizer of a KmerInt using canonical representation (testing both the positive-strand and its reverse complement)
*/
func (kmerInt KmerInt) Minimize() MinInt {
    rcKmerInt := kmerInt.revComp()
    // Despite being minimizers (and therefore expected to be
    // MinInts), we make possMin and currMin KmerInts for ease of
    // manipulation and then convert upon returning.
    // Initialize the current best minimizer candidate as MAX_INT
    currMin := ^KmerInt(0)
    var possMin KmerInt

    var i uint8 // the positive-strand index of the minimizers tested on the iteration
    for i = 0; i <= k-m; i++ {
        possMin = mMask & kmerInt
        if possMin < currMin {
            currMin = possMin
        }

        possMin = mMask & rcKmerInt
        if possMin < currMin {
            currMin = possMin
        }

        kmerInt >>= 2
        rcKmerInt >>= 2
    }

    return MinInt(currMin)
}

/*
   Converts a TextSeq to the more memory-efficient Seq type.
   Upper- and lower-case base bytes are currently supported, but stable code should immediately convert to lower-case.
   The logic works and is sane, but could be altered in the future for brevity and efficiency.
*/
func (textSeq TextSeq) Seq() Seq {
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
   Returns the index of the KmerInt's minimizer.
   The positive-strand index of the minimizer is minKey % k.
   If the MinKey is >= k, the minimizer is the reverse complement of the indexed minimizer.
*/
func (kmerInt KmerInt) MinKey() MinKey {
    rcKmerInt := kmerInt.revComp()
    numPossMins := k - m + 1
    // Despite being minimizers (and therefore expected to be
    // MinInts), we make possMin and currMin KmerInts for ease of
    // manipulation and then convert upon returning.
    currMin := ^KmerInt(0) // initialize the current best minimizer candidate as MAX_INT
    var currKey MinKey = 0

    var i uint8 // the positive-strand start index of the minimizers tested in this iteration
    for i = 0; i < numPossMins; i++ {
        possMin := minKmerInt(mMask&kmerInt, mMask&rcKmerInt)

        possMin = mMask & kmerInt
        if possMin < currMin {
            currMin = possMin
            currKey = MinKey(i)
        }

        possMin = mMask & rcKmerInt
        if possMin < currMin {
            currMin = possMin
            currKey = MinKey(k + i)
        }

        kmerInt >>= 2
        rcKmerInt >>= 2
    }

    return currKey
}

/*
   The canonical representation of a sequence is the lexicographically smaller of its positive-strand and its reverse complement.
*/
func (kmerInt KmerInt) canonicalRepr() KmerInt {
    revComp := kmerInt.revComp()
    if revComp < kmerInt {
        return revComp
    } else {
        return kmerInt
    }
}

/*
   The canonical representation of a sequence is the lexicographically smaller of its positive-strand and its reverse complement.
*/
func (minInt MinInt) canonicalRepr() MinInt {
    revComp := minInt.revComp()
    if revComp < minInt {
        return revComp
    } else {
        return minInt
    }
}

/*
   A more declarative and modifiable accessor function.
   While it would almost certainly be inlined, this is such a performance-critical operation that this function isn't currently used.
*/
func (kmer Kmer) Int() KmerInt {
    return *(*KmerInt)(unsafe.Pointer(&kmer))
}

/*
   A more declarative and modifiable accessor function.
   While it would almost certainly be inlined, this is such a performance-critical operation that this function isn't currently used.
*/
func (kmer Kmer) ClassID() ClassID {
    return *(*ClassID)(unsafe.Pointer(&kmer[8]))
}

func (kmer *Kmer) SetInt(kmerInt KmerInt) {
    *(*KmerInt)(unsafe.Pointer(kmer)) = kmerInt
}

func (kmer *Kmer) SetClassID(classID ClassID) {
    *(*ClassID)(unsafe.Pointer(&kmer[8])) = classID
}
