package repeatgenome

import (
    "sort"
//    "fmt"
    "sync"
)

type matchSpan struct {
    start uint64
    end uint64
    className string
}

type matchSpans []matchSpan

func (matchSpans matchSpans) Len() int {
    return len(matchSpans)
}

func (matchSpans matchSpans) Swap(i, j int) {
    matchSpans[i], matchSpans[j] = matchSpans[j], matchSpans[i]
}

func (matchSpans matchSpans) Less(i, j int) bool {
    return matchSpans[i].start < matchSpans[j].start
}


func (rg *RepeatGenome) GetMatchSpans() map[string]matchSpans {
    // maps a seq's name to its matchSpans
    spanMap := make(map[string]matchSpans)

    for _, match := range rg.Matches {
        /*
        if _, exists := spanMap[match.SeqName]; !exists {
            spanMap[match.SeqName] = make
        */
        seqName := match.SeqName
        matchSpan := matchSpan{match.SeqStart, match.SeqEnd, match.RepeatName}
        spanMap[seqName] = append(spanMap[seqName], matchSpan)
    }

    for _, matchSpans := range spanMap {
        sort.Sort(matchSpans)
    }

    return spanMap
}


type KRespPair struct {
    KmerInt KmerInt
    ClassName string
}


type MRespPair struct {
    MinInt MinInt
    ClassName string
}


func (rg *RepeatGenome) SplitChromsK() (chan KRespPair, chan KRespPair) {
    spanMap := rg.GetMatchSpans()

    repChan, nonrepChan := make(chan KRespPair, 5000), make(chan KRespPair, 5000)
    var wg = new(sync.WaitGroup)

    for chromName, matchSpans := range spanMap {
        // used to determine where in seq each iteration starts
        seq := rg.chroms[chromName][chromName]
        var start uint64 = 0
        for _, matchSpan := range matchSpans {
            mStart, mEnd, className := matchSpan.start, matchSpan.end, matchSpan.className
            // necessary because of overlapping matches
            if start < mStart {
                wg.Add(1)
                go sendKmers(seq[start:mStart], className, k, nonrepChan, wg)
            }
            wg.Add(1)
            go sendKmers(seq[mStart:mEnd], className, k, repChan, wg)
            start = mEnd
        }
        wg.Add(1)
        go sendKmers(seq[start:], "NON_REPEAT", k, nonrepChan, wg)
    }

    go func() {
        wg.Wait()
        close(repChan)
        close(nonrepChan)
    }()

    return repChan, nonrepChan
}


func sendKmers(seq TextSeq, className string, k uint8, c chan KRespPair, wg *sync.WaitGroup) {
    defer wg.Done()
    var numKmers = len(seq) - int(k) + 1
KmerLoop:
    for i := 0; i < numKmers; i++ {
        for j := int(k) + i - 1; j >= i; j-- {
            if seq[j] == byte('n') {
                i += j - i
                continue KmerLoop
            }
        }
        kmerInt := seq[i:i+int(k)].kmerInt()
        c <- KRespPair{kmerInt, className}
        c <- KRespPair{kmerInt.revComp(), className}
    }
}


/*
func (rg *RepeatGenome) getMatchKmers(matchChan chan *bioutils.Match, respChan chan ResponsePair, wg *sync.WaitGroup) {
    k_ := uint64(k)
    for match := range matchChan {
        start, end := match.SeqStart, match.SeqEnd
        matchSeq := rg.chroms[match.SeqName][match.SeqName][start:end]

        if len(matchSeq) < int(k) {
            continue
        }
        // includes kmers containing n's, which are ignored
        numRawKmers := end - start - k_ + 1

        var i uint64
    KmerLoop:
        for i = 0; i < numRawKmers; i++ {

            var j int64
            for j = int64(k_) - 1; j >= 0; j-- {
                if matchSeq[i+uint64(j)] == 'n' {
                    i += uint64(j)
                    continue KmerLoop
                }
            }

            kmerInt := TextSeq(matchSeq[i : i+k_]).kmerInt().canonicalRepr()
            var kmer Kmer
            (&kmer).SetInt(kmerInt)
            (&kmer).SetClassID(rg.matchNodes[match].ID)
            respChan <- ResponsePair{kmer, kmerInt.Minimize()}
        }
    }
    wg.Done()
}


/*
    Returns a channel to which all repeat minimizers are pushed.
/*
func (rg *RepeatGenome) minChan() chan ResponsePair {
    matchChan := make(chan *bioutils.Match, 500)
    go func() {
        for i := range rg.Matches {
            matchChan <- &rg.Matches[i]
        }
        close(matchChan)
    }()

    var wg = new(sync.WaitGroup)
    wg.Add(numCPU)
    respChan := make(chan ResponsePair) // do not buffer without additional WaitGroup

    go func() {
        wg.Wait()
        close(respChan)
    }()

    for i := 0; i < numCPU; i++ {
        go rg.getMatchKmers(matchChan, respChan, wg)
    }

    return respChan
}
*/


func (rg *RepeatGenome) SplitChromsM() (chan MRespPair, chan MRespPair) {
    spanMap := rg.GetMatchSpans()

    repChan, nonrepChan := make(chan MRespPair, 5000), make(chan MRespPair, 5000)
    var wg = new(sync.WaitGroup)

    for chromName, matchSpans := range spanMap {
        // used to determine where in seq each iteration starts
        seq := rg.chroms[chromName][chromName]
        var start uint64 = 0
        for _, matchSpan := range matchSpans {
            mStart, mEnd, className := matchSpan.start, matchSpan.end, matchSpan.className
            // necessary because of overlapping matches
            if start < mStart {
                wg.Add(1)
                go sendMins(seq[start:mStart], className, m, nonrepChan, wg)
            }
            wg.Add(1)
            go sendMins(seq[mStart:mEnd], className, m, repChan, wg)
            start = mEnd
        }
        wg.Add(1)
        go sendMins(seq[start:], "NON_REPEAT", m, nonrepChan, wg)
    }

    go func() {
        wg.Wait()
        close(repChan)
        close(nonrepChan)
    }()

    return repChan, nonrepChan
}


func sendMins(seq TextSeq, className string, m uint8, c chan MRespPair, wg *sync.WaitGroup) {
    defer wg.Done()
    var numMins = len(seq) - int(m) + 1
MinLoop:
    for i := 0; i < numMins; i++ {
        for j := int(m) + i - 1; j >= i; j-- {
            if seq[j] == byte('n') {
                i += j - i
                continue MinLoop
            }
        }
        minInt := seq[i:i+int(m)].minInt()
        c <- MRespPair{minInt, className}
        c <- MRespPair{minInt.revComp(), className}
    }
}
