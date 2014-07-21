package repeatgenome

import (
    "fmt"
)
// General sequence-manipulation functions.

func seqToInt(seq TextSeq) uint64 {
    if len(seq) < 1 || len(seq) > 31 {
        panic("seqToInt() can only int-ize sequences where 0 < length < 32")
    }
    var seqInt uint64 = 0
    for _, c := range seq {
        seqInt = seqInt << 2
        switch c {
        case 'a':
            break
        case 'c':
            seqInt |= 1
            break
        case 'g':
            seqInt |= 2
            break
        case 't':
            seqInt |= 3
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to seqToInt(): %c", c))
        }
    }
    return seqInt
}

func revCompToInt(seq TextSeq) uint64 {
    if len(seq) < 1 || len(seq) > 31 {
        panic("revCompToInt() can only int-ize sequences where 0 < length < 32")
    }
    var seqInt uint64 = 0
    for i := range seq {
        seqInt = seqInt << 2
        switch seq[len(seq)-(i+1)] {
        case 'a':
            seqInt |= 3
            break
        case 'c':
            seqInt |= 2
            break
        case 'g':
            seqInt |= 1
            break
        case 't':
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to revCompToInt(): %c", seq[i]))
        }
    }
    return seqInt
}

func intRevComp(posInt uint64, seqLen uint8) uint64 {
    var revComp uint64 = 0
    var i uint8
    for i = 0; i < seqLen; i++ {
        // should execute before every loop but the first
        if i != 0 {
            revComp <<= 2
            posInt >>= 2
        }
        switch posInt & 3 {
        case 0:
            revComp |= 3
            break
        case 1:
            revComp |= 2
            break
        case 2:
            revComp |= 1
            break
        case 3:
            break
        default:
            panic("bit manipulation logic error in intRevComp()")
        }
    }
    return revComp
}

func getMinimizer(kmer uint64, k, m uint8) (uint8, uint64) {
    if m > k || m < 1 {
        panic("getMinimizer(): m must be <= k and > 0")
    }

    // stores the index of the leftmost base included in the minimizer
    var currOffset uint8 = 0
    numExtraBases := 32 - k
    revCompKmer := intRevComp(kmer, k)
    currMin := kmer >> uint64(64-2*(32-k+m))
    possMin := currMin
    numHangingBases := k - m
    var i uint8
    // i is the index of the offset
    for i = 0; i <= numHangingBases; i++ {
        // overflow off the first excluded base
        possMin = kmer << (2 * (numExtraBases + i))
        // return to proper alignment
        possMin >>= 64 - 2*m

        if possMin < currMin {
            currMin = possMin
            currOffset = i
        }

        possMin = revCompKmer << (2 * (numExtraBases + i))
        possMin >>= 64 - 2*m

        if possMin < currMin {
            currMin = possMin
            currOffset = numHangingBases - i
        }
    }

    return currOffset, currMin
}

func revComp(seq TextSeq) TextSeq {
    var revCompSeq = make(TextSeq, 0, len(seq))
    for i := range seq {
        switch seq[len(seq)-i-1] {
        case 'a':
            revCompSeq = append(revCompSeq, 't')
            break
        case 't':
            revCompSeq = append(revCompSeq, 'a')
            break
        case 'c':
            revCompSeq = append(revCompSeq, 'g')
            break
        case 'g':
            revCompSeq = append(revCompSeq, 'c')
            break
        default:
            panic("byte other than 'a', 'c', 'g', or 't' supplied to revComp")
        }
    }
    return revCompSeq
}

// The logic for determining the minimizer
// Currently, it uses simple lexicographic ordering
func seqLessThan(a, b TextSeq) bool {
    // min function manually inlined for speed - dubious
    var size int
    if len(a) < len(b) {
        size = len(a)
    } else {
        size = len(b)
    }
    for i := 0; i < size; i++ {
        if a[i] < b[i] {
            return true
        }
        if a[i] > b[i] {
            return false
        }
    }
    return false
}

func GetSeq(seqStr TextSeq) Seq {
    // ceiling division of len(seqStr) by 4
    var numBytes uint64 = 1 + ((uint64(len(seqStr)) - 1) / 4)
    seq := Seq{make([]byte, numBytes, numBytes), uint64(len(seqStr))}
    basesInFirstByte := 1 + ((seq.Len - 1) % numBytes)

    seqStrInd := 0
    var i uint64
    for i = 0; i < basesInFirstByte; i++ {
        seq.Bases[0] <<= 2

        switch seqStr[seqStrInd] {
        case 'a':
            break
        case 'c':
            seq.Bases[0] |= 1
            break
        case 'g':
            seq.Bases[0] |= 2
            break
        case 't':
            seq.Bases[0] |= 3
            break
        default:
            panic("byte other than 'a', 'c', 'g', or 't' supplied to revComp")
        }

        seqStrInd++
    }

    for i := 1; i < len(seq.Bases); i++ {
        for j := 0; j < 4; j++ {
            seq.Bases[i] <<= 2

            switch seqStr[seqStrInd] {
            case 'a':
                break
            case 'c':
                seq.Bases[i] |= 1
                break
            case 'g':
                seq.Bases[i] |= 2
                break
            case 't':
                seq.Bases[i] |= 3
                break
            default:
                panic("byte other than 'a', 'c', 'g', or 't' supplied to revComp")
            }

            seqStrInd++
        }
    }

    return seq
}
