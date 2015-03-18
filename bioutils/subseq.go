package bioutils

import (
    "fmt"
)

/*
   The canonical representation of a sequence is the lexicographically smaller of its positive-strand and its reverse complement.
*/
func CanonicalRepr32(seq uint32) uint32 {
    revComp := RevComp32(seq)
    if revComp < seq {
        return revComp
    } else {
        return seq
    }
}

/*
   The canonical representation of a sequence is the lexicographically smaller of its positive-strand and its reverse complement.
*/
func CanonicalRepr64(seq uint64) uint64 {
    revComp := RevComp64(seq)
    if revComp < seq {
        return revComp
    } else {
        return seq
    }
}


/*
   Returns the index of the uint64's minimizer.
   The positive-strand index of the minimizer is minKey % K.
   If the MinKey is >= K, the minimizer is the reverse complement of the indexed minimizer.
*/
func MinKey(seq uint64) uint8 {
    rcKmerInt := RevComp64(seq)
    numPossMins := uint8(K - M + 1)
    // Despite being minimizers (and therefore expected to be
    // uint32s), we make possMin and currMin uint64s for ease of
    // manipulation and then convert upon returning.
    currMin := ^uint64(0) // initialize the current best minimizer candidate as MAX_INT
    var currKey uint8 = 0

    var i uint8 // the positive-strand start index of the minimizers tested in this iteration
    for i = 0; i < numPossMins; i++ {
        possMin := minU64(MMask&seq, MMask&rcKmerInt)

        possMin = MMask & seq
        if possMin < currMin {
            currMin = possMin
            currKey = i
        }

        possMin = MMask & rcKmerInt
        if possMin < currMin {
            currMin = possMin
            currKey = K + i
        }

        seq >>= 2
        rcKmerInt >>= 2
    }

    return currKey
}


func BytesToU64(seq []byte) uint64 {
    if len(seq) < 1 || len(seq) > 31 {
        panic("TextSeq.kmerInt() can only int-ize sequences where 0 < length < 32")
    }
    var kmerInt uint64 = 0
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

func BytesToU32(seq []byte) uint32 {
    if len(seq) < 1 || len(seq) > 15 {
        panic("TextSeq.minInt() can only int-ize sequences where 0 < length < 16")
    }
    var minInt uint32 = 0
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
   Return the minimizer of a uint64 using canonical representation (testing both the positive-strand and its reverse complement)
*/
func Minimize(seq uint64) uint32 {
    rcKmerInt := RevComp64(seq)
    // Despite being minimizers (and therefore expected to be
    // uint32s), we make possMin and currMin KmerInts for ease of
    // manipulation and then convert upon returning.
    // Initialize the current best minimizer candidate as MAX_INT
    currMin := ^uint64(0)
    var possMin uint64

    var i uint8 // the positive-strand index of the minimizers tested on the iteration
    for i = 0; i <= K-M; i++ {
        possMin = MMask & seq
        if possMin < currMin {
            currMin = possMin
        }

        possMin = MMask & rcKmerInt
        if possMin < currMin {
            currMin = possMin
        }

        seq >>= 2
        rcKmerInt >>= 2
    }

    return uint32(currMin)
}

/*
   Returns the reverse complement of a sequence stored in
   a uint32.
*/
func RevComp32(seq uint32) uint32 {
    var revComp uint32
    var i uint8
    for i = 0; i < M; i++ {
        revComp <<= 2
        revComp |= 3 &^ seq // BEWARE: &^ is not commutative
        seq >>= 2
    }
    return revComp
}

/*
   Returns the reverse complement of a sequence stored in
   a uint64.
*/
func RevComp64(seq uint64) uint64 {
    var revComp uint64
    var i uint8
    for i = 0; i < K; i++ {
        revComp <<= 2
        revComp |= 3 &^ seq // BEWARE: &^ is not commutative
        seq >>= 2
    }
    return revComp
}

func U64ToBytes(kmerInt uint64) []byte {
    b := make([]byte, K)
    var i uint8
    for i = 0; i < K; i++ {
        // this tricky bit arithmetic shifts the two bits of interests
        // to the two rightmost positions, then selects them with the
        // and statement
        switch (kmerInt >> (2 * (K - i - 1))) & 3 {
        case 0:
            b = append(b, 'a')
            break
        case 1:
            b = append(b, 'c')
            break
        case 2:
            b = append(b, 'g')
            break
        case 3:
            b = append(b, 't')
            break
        default:
            panic("error in printSeqInt() base selection")
        }
    }

    return b
}

func U32ToBytes(minInt uint32) []byte {
    b := make([]byte, K)
    var i uint8
    for i = 0; i < M; i++ {
        // this tricky bit arithmetic shifts the two bits of interests
        // to the two rightmost positions, then selects them with the
        // and statement
        switch (minInt >> (2 * (M - i - 1))) & 3 {
        case 0:
            b = append(b, 'a')
            break
        case 1:
            b = append(b, 'c')
            break
        case 2:
            b = append(b, 'g')
            break
        case 3:
            b = append(b, 't')
            break
        default:
            panic("error in printSeqInt() base selection")
        }
    }

    return b
}
