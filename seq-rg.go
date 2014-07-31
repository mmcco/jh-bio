package repeatgenome

// General sequence-manipulation functions.

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

func (kmerInt KmerInt) revComp() KmerInt {
    var revComp KmerInt
    var i uint8
    for i = 0; i < k; i++ {
        revComp <<= 2
        revComp |= 3 &^ kmerInt    // BEWARE: &^ is not commutative
        kmerInt >>= 2
    }
    return revComp
}

func (minInt MinInt) revComp() MinInt {
    var revComp MinInt
    var i uint8
    for i = 0; i < m; i++ {
        revComp <<= 2
        revComp |= 3 &^ minInt    // BEWARE: &^ is not commutative
        minInt >>= 2
    }
    return revComp
}

func (kmerInt KmerInt) Minimize() MinInt {
    rcKmerInt := kmerInt.revComp()
    // despite being minimizers (and therefore expected to be MinInts), we make possMin and currMin KmerInts for ease of manipulation and then convert upon returning
    // initialize the current best minimizer candidate as MAX_INT
    currMin := ^KmerInt(0)
    // i is the index of the offset
    var possMin KmerInt
    var i uint8
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

func (textSeq TextSeq) Seq() Seq {
    numBases := len(textSeq)
    numBytes := ceilDiv(numBases, 4)
    seq := Seq{
        Bytes: make([]byte, numBytes, numBytes),
        Len: uint64(numBases),
    }

    // determines how much to shift the byte of interest
    var shift uint8 = 6
    for i := 0; i < len(textSeq); i++ {
        byteInd := i / 4

        switch textSeq[i] {
        case 'a':
            // already zero
            break
        case 'c':
            seq.Bytes[byteInd] |= 1 << shift
            break
        case 'g':
            seq.Bytes[byteInd] |= 2 << shift
            break
        case 't':
            seq.Bytes[byteInd] |= 3 << shift
            break
        case 'A':
            // already zero
            break
        case 'C':
            seq.Bytes[byteInd] |= 1 << shift
            break
        case 'G':
            seq.Bytes[byteInd] |= 2 << shift
            break
        case 'T':
            seq.Bytes[byteInd] |= 3 << shift
            break
        default:
            panic("TextSeq.GetSeq(): byte other than 'a', 'c', 'g', or 't' encountered")
        }

        // starts at 6, wraps around to 254, which mods to 0
        // therefore loops 6 -> 4 -> 2 -> 0 -> ...
        shift = (shift - 2) % 8
    }

    return seq
}

func (seq Seq) Subseq(a, b uint64) Seq {
    numBytes := ceilDiv_U64(b-a, 4)
    subseq := Seq{
        Bytes: make([]byte, numBytes, numBytes),
        Len: b - a,
    }
    var i uint64
    for i = 0; i < subseq.Len; i++ {
        subseq.Bytes[i / 4] |= seq.GetBase(i + a) << (6 - 2*(i%4))
    }
    return subseq
}

func (seq Seq) GetBase(i uint64) uint8 {
    thisByte := seq.Bytes[i / 4]
    return thisByte >> (6 - 2*(i%4))
}

func (kmerInt KmerInt) MinKey() MinKey {
    rcKmerInt := kmerInt.revComp()
    numPossMins := k - m + 1
    // despite being minimizers (and therefore expected to be MinInts), we make possMin and currMin KmerInts for ease of manipulation and then convert upon returning
    // initialize the current best minimizer candidate as MAX_INT
    currMin := ^KmerInt(0)
    var currKey MinKey = 0
    // i is the index of the offset
    var i uint8
    for i = 0; i < numPossMins; i++ {
        possMin := minKmerInt(mMask & kmerInt, mMask & rcKmerInt)

        possMin = mMask & kmerInt
        if possMin < currMin {
            currMin = possMin
            currKey = MinKey(k - m - i)
        }

        possMin = mMask & rcKmerInt
        if possMin < currMin {
            currMin = possMin
            currKey = MinKey(i)
        }

        kmerInt >>= 2
        rcKmerInt >>= 2
    }

    return currKey
}

func (kmerInt KmerInt) canonicalRepr() KmerInt {
    revComp := kmerInt.revComp()
    if revComp < kmerInt {
        return revComp
    } else {
        return kmerInt
    }
}

func (minInt MinInt) canonicalRepr() MinInt {
    revComp := minInt.revComp()
    if revComp < minInt {
        return revComp
    } else {
        return minInt
    }
}

func (kmer Kmer) Int() KmerInt {
    return *(*KmerInt)(unsafe.Pointer(&kmer))
}
