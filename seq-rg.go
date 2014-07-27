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

func (kmerInt KmerInt) revComp() KmerInt {
    var revComp KmerInt = 0
    var i uint8
    for i = 0; i < k; i++ {
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

func (minInt MinInt) revComp() MinInt {
    var revComp MinInt = 0
    var i uint8
    for i = 0; i < m; i++ {
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

func (kmerInt KmerInt) Minimize() MinInt {
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
        possMin = minKmerInt(possMin, possMin.revComp())

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

func (textSeq TextSeq) Seq() Seq {
    numBases := uint64(len(textSeq))
    numBytes := ceilDiv_U64(numBases, 4)
    seq := Seq{
        Bytes: make([]byte, numBytes, numBytes),
        Len: numBases,
    }
    var i uint64
    for i = 0; i < uint64(len(textSeq)); i++ {
        byteInd := i / 4
        shift := 6 - 2*(i%4)
        switch textSeq[i] {
        case 'a':
            // already zero
            break
        case 'c':
            seq.Bytes[byteInd] |= (uint8(1) << shift)
            break
        case 'g':
            seq.Bytes[byteInd] |= (uint8(2) << shift)
            break
        case 't':
            seq.Bytes[byteInd] |= (uint8(3) << shift)
            break
        case 'A':
            // already zero
            break
        case 'C':
            seq.Bytes[byteInd] |= (uint8(1) << shift)
            break
        case 'G':
            seq.Bytes[byteInd] |= (uint8(2) << shift)
            break
        case 'T':
            seq.Bytes[byteInd] |= (uint8(3) << shift)
            break
        default:
            panic("TextSeq.GetSeq(): byte other than 'a', 'c', 'g', or 't' encountered")
        }
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

func Debug() {
    ts := TextSeq("tgaatatacgtagctctagctagcgcttatatg")
    fmt.Println("ts:", ts)
    s := ts.Seq()
    fmt.Print("s: ")
    s.Print()
    fmt.Println()
    fmt.Println("s[4]: ", s.GetBase(4))
    fmt.Println("s[-1]: ", s.GetBase(s.Len - 1))
    fmt.Print("s[:7]: ")
    s.Subseq(0, 7).Print()
    fmt.Println()
    fmt.Print("s[7:] ")
    s.Subseq(7, s.Len).Print()
    fmt.Println()
}

/*
func (kmerInt KmerInt) MinKey() MinKey {
    var currMin MinInt = ^MinInt(0)
    var currKey MinKey = ^MinKey(0)

    var mask KmerInt = ^KmerInt(0)
    mask >>= 32 - m

    var i uint8
    for i = 0; i <= k - m; i++ {
        // we start with the last minimizer and work back - this makes masking easier and faster
        possMin := MinInt((kmerInt | mask) >> i)
        rcPossMin := possMin.revComp()

        if possMin < currMin {
            currMin = possMin
            currKey = MinKey(k - m - i)
        }
        if rcPossMin < currMin {
            currMin = rcPossMin
            currKey = MinKey(k - m - i + 32)
        }

        mask <<= 2
    }

    return currKey
}
*/

/*
func (minPair MinPair) getMin() MinInt {
    minKey := *(*MinKey)(unsafe.Pointer(&minPair[8]))
    isRevComp := false
    if minKey > 31 {
        isRevComp = true
        minKey -= 32
    }
    minInt := MinInt(*(*KmerInt)(unsafe.Pointer(&minPair)))
    minInt <<= (64 - 2*k + uint8(minKey))
    minInt >>= (64 - 2*m)
    
    if isRevComp {
        return minInt.revComp()
    } else {
        return minInt
    }
}
*/
