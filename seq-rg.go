package repeatgenome

import (
    "fmt"
)
// General sequence-manipulation functions.

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

func (seq TextSeq) revCompKmerInt() KmerInt {
    if len(seq) < 1 || len(seq) > 31 {
        panic("TextSeq.revCompKmerInt() can only int-ize sequences where 0 < length < 32")
    }
    var kmerInt KmerInt = 0
    for i := range seq {
        kmerInt = kmerInt << 2
        switch seq[len(seq)-(i+1)] {
        case 'a':
            kmerInt |= 3
            break
        case 'c':
            kmerInt |= 2
            break
        case 'g':
            kmerInt |= 1
            break
        case 't':
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to TextSeq.revCompKmerInt(): %c", seq[i]))
        }
    }
    return kmerInt
}

func (seq TextSeq) revCompMinInt() MinInt {
    if len(seq) < 1 || len(seq) > 15 {
        panic("TextSeq.revCompMinInt() can only int-ize sequences where 0 < length < 15")
    }
    var minInt MinInt = 0
    for i := range seq {
        minInt = minInt << 2
        switch seq[len(seq)-(i+1)] {
        case 'a':
            minInt |= 3
            break
        case 'c':
            minInt |= 2
            break
        case 'g':
            minInt |= 1
            break
        case 't':
            break
        default:
            panic(fmt.Errorf("byte other than 'a', 'c', 'g', or 't' passed to TextSeq.revCompMinInt(): %c", seq[i]))
        }
    }
    return minInt
}

func (kmerInt KmerInt) revComp(seqLen uint8) KmerInt {
    if seqLen > 31 {
        panic("seqLen provided to KmerInt.revComp too large")
    }

    var revComp KmerInt = 0
    var i uint8
    for i = 0; i < seqLen; i++ {
        // should execute before every loop but the first
        if i != 0 {
            revComp <<= 2
            kmerInt >>= 2
        }
        switch kmerInt & 3 {
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
            panic("bit manipulation logic error in KmerInt.revComp()")
        }
    }
    return revComp
}

func (minInt MinInt) revComp(seqLen uint8) MinInt {
    if seqLen > 15 {
        panic("seqLen provided to MinInt.revComp too large")
    }

    var revComp MinInt = 0
    var i uint8
    for i = 0; i < seqLen; i++ {
        // should execute before every loop but the first
        if i != 0 {
            revComp <<= 2
            minInt >>= 2
        }
        switch minInt & 3 {
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
            panic("bit manipulation logic error in KmerInt.revComp()")
        }
    }
    return revComp
}

func (kmerInt KmerInt) Minimize(k, m uint8) MinInt {
    if m > k || m < 1 {
        panic("KmerInt.Minimize(): m must be <= k and > 0")
    }

    // stores the index of the leftmost base included in the minimizer
    numExtraBases := 32 - k
    numHangingBases := k - m
    // despite being minimizers (and therefore expected to be MinInts), we make possMin and currMin KmerInts for ease of manipulation and then convert upon returning
    // initialize the current best minimizer candidate as MAX_INT
    currMin := ^KmerInt(0)
    // i is the index of the offset
    var i uint8
    for i = 0; i <= numHangingBases; i++ {
        // overflow off the first excluded base
        possMin := kmerInt << (2 * (numExtraBases + i))
        // return to proper alignment
        possMin >>= 64 - 2*m
        possMin = minKmerInt(possMin, possMin.revComp(m))

        if possMin < currMin {
            currMin = possMin
        }
    }

    return MinInt(currMin)
}

func (seq TextSeq) revComp() TextSeq {
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
// The function's name allows it to be used in a sort.Interface implementation, if we ever need one.
func (a TextSeq) Less(b TextSeq) bool {
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

func (seqStr TextSeq) Seq() Seq {
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
            panic("byte other than 'a', 'c', 'g', or 't' supplied to TextSeq.Seq()")
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

func (seq Seq) subseq(start, end uint64) Seq {
    subSeqLen := end - start
    if subSeqLen > seq.Len {
        panic("Seq.subseq(): subsequence larger than parent sequence requested")
    }
    // the actual number of bytes needed in Seq.Seq, calculated with ceiling division
    numBytes := ((2*subSeqLen - 1) / 8) + 1
    subSeq := Seq{make([]byte, numBytes, numBytes), subSeqLen}

    //seqExtraBits := seq.numExtraBits()
    subSeqExtraBits := subSeq.numExtraBits()
    // the byte offset (floored, as there is likely a bit offset) of the subseq into seq
    //byteOffset := (seq.Len - subSeqLen) / 4
    //var bitOffset uint8 = seqExtraBits - subSeqExtraBits

    var i uint64
    for i = 0; i < subSeqLen; i++ {
        byteIndex := (i*2 + uint64(subSeqExtraBits)) / 8
        subSeq.Bases[byteIndex] <<= 2
        // this method call is probably a performance drag
        subSeq.Bases[byteIndex] |= seq.getBase(i + start)
    }
    
    return subSeq
}

func (seq Seq) numExtraBits() uint8 {
    return uint8(8 * uint64(len(seq.Bases)) - 2 * seq.Len)
}

func (seq Seq) getBase(i uint64) uint8 {
    // the bit offset of this base in seq.Seq
    bitOffset := i*2 + uint64(seq.numExtraBits())
    // the byte index of seq.Seq in which our base exists
    byteIndex := bitOffset / 8
    var base uint8 = seq.Bases[byteIndex]
    // the number of bits that must be shifted to move the base of interest into the least significant bits
    bitsToShift := 6 - (bitOffset - 8*byteIndex)

    return (base >> bitsToShift) | 3
}
