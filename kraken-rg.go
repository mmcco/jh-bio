package repeatgenome

import (
    "fmt"
    "runtime"
    "sort"
    "sync"
    "unsafe"
)

func (classTree *ClassTree) getLCA(cnA, cnB *ClassNode) *ClassNode {
    for cnBWalker := cnB; cnBWalker != classTree.Root; cnBWalker = cnBWalker.Parent {
        for cnAWalker := cnA; cnAWalker != classTree.Root; cnAWalker = cnAWalker.Parent {
            if cnBWalker == cnAWalker {
                return cnBWalker
            }
        }
    }
    return classTree.Root
}

func (rg *RepeatGenome) getKmer(kmerInt KmerInt) *Kmer {
    minimizer := kmerInt.Minimize()
    
    // simple binary search within the range of RepeatGenome.Kmers that has this kmer's minimizer
    i, exists := rg.MinOffsets[minimizer]
    if !exists {
        return nil
    }
    j := i + uint64(rg.MinCounts[minimizer])
    for i < j {
        x := (j+i)/2
        thisKmerInt := *(*KmerInt)(unsafe.Pointer(&rg.Kmers[x]))
        if thisKmerInt == kmerInt {
            return &rg.Kmers[x]
        } else if thisKmerInt < kmerInt {
            i = x + 1
        } else {
            j = x
        }
    }

    return nil
}

// Could probably use a map to prevent duplicates in the future, although this wouldn't make a difference if the slice size were still manually chosen
func (rg *RepeatGenome) getRawMins(matchChan chan *Match, respChan chan MinInt, wg *sync.WaitGroup) {
    k_ := uint64(k)
    for match := range matchChan {
        start, end := match.SeqStart, match.SeqEnd
        matchSeq := rg.chroms[match.SeqName][match.SeqName][start : end]

        if len(matchSeq) < int(k) {
            continue
        }
        // includes kmers containing n's, which are ignored
        numKmers := end - start - k_ + 1

        var i uint64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            var j int64
            for j = int64(k_) - 1; j >= 0; j-- {
                if matchSeq[i+uint64(j)] == 'n' {
                    i += uint64(j)
                    continue KmerLoop
                }
            }

            respChan <- TextSeq(matchSeq[i : i+k_]).kmerInt().canonicalRepr().Minimize()
        }
    }
    wg.Done()
}

type ResponsePair struct {
    Kmer Kmer
    MinInt MinInt
}

// Could probably use a map to prevent duplicates in the future, although this wouldn't make a difference if the slice size were still manually chosen
func (rg *RepeatGenome) getMatchKmers(matchChan chan *Match, respChan chan ResponsePair, wg *sync.WaitGroup) {
    k_ := uint64(k)
    for match := range matchChan {
        start, end := match.SeqStart, match.SeqEnd
        matchSeq := rg.chroms[match.SeqName][match.SeqName][start : end]

        if len(matchSeq) < int(k) {
            continue
        }
        // includes kmers containing n's, which are ignored
        numKmers := end - start - k_ + 1

        var i uint64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            var j int64
            for j = int64(k_) - 1; j >= 0; j-- {
                if matchSeq[i+uint64(j)] == 'n' {
                    i += uint64(j)
                    continue KmerLoop
                }
            }

            kmerInt := TextSeq(matchSeq[i : i+k_]).kmerInt().canonicalRepr()
            var kmer Kmer
            *(*KmerInt)(unsafe.Pointer(&kmer)) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = match.ClassNode.ID
            respChan <- ResponsePair{kmer, kmerInt.Minimize()}
        }
    }
    wg.Done()
}

// we first count the number of non-unique kmers and the number of non-unique minimizers associated with each kmer
func (rg *RepeatGenome) krakenFirstPass() (numRawKmers uint64, rawMinCounts map[MinInt]uint32) {
    matchChan := make(chan *Match)
    go func() {
        for i := range rg.Matches {
            matchChan <- &rg.Matches[i]
        }
        close(matchChan)
    }()

    wg := new(sync.WaitGroup)
    wg.Add(numCPU)
    respChan := make(chan MinInt)   // cannot buffer without additional WaitGroup

    go func() {
        wg.Wait()
        close(respChan)
    }()

    for i := 0; i < numCPU; i++ {
        go rg.getRawMins(matchChan, respChan, wg)
    }

    rawMinCounts = make(map[MinInt]uint32)
    numRawKmers = 0
    for minInt := range respChan {
        rawMinCounts[minInt]++
        numRawKmers++
    }

    wg.Wait()

    return numRawKmers, rawMinCounts
}

func (rg *RepeatGenome) sortRawKmers(rawKmers Kmers, rawMinCounts map[MinInt]uint32) {
    wg := new(sync.WaitGroup)
    wg.Add(numCPU)
    kmersChan := make(chan Kmers)
    for i := 0; i < numCPU; i++ {
        go func() {
            for kmers := range kmersChan {
                sort.Sort(kmers)
                /*
                // a simple insertion sort, which is faster than sort.Sort's QuickSort for such small lists
                for i := 1; i < len(kmers); i++ {
                    for j := i; j > 0 && kmers.Less(j, j-1); j-- {
                        kmers.Swap(j, j-1)
                    }
                }
                */
            }
            wg.Done()
        }()
    }

    var start, end uint64 = 0, 0
    for _, minInt := range rg.SortedMins {
        cnt, exists := rawMinCounts[minInt]
        if !exists {
            minInt.print(); fmt.Println("does not exist in rawMinCounts")
        }
        end += uint64(cnt)
        kmersChan <- rawKmers[start : end]
        start = end
    }

    close(kmersChan)
    wg.Wait()
}

func (rg *RepeatGenome) getRawKmers(numRawKmers uint64, rawMinCounts map[MinInt]uint32) Kmers {
    // we need a map telling us the next place to insert a kmer associated with any given minimizer
    // we simply increment a minimizer's offset when we insert a kmer associated with it
    locMap := make(map[MinInt]uint64, len(rawMinCounts))
    if len(rg.SortedMins) > 0 {
        locMap[rg.SortedMins[0]] = 0

        for i := 1; i < len(rg.SortedMins); i++ {
            thisMin := rg.SortedMins[i]
            lastMin := rg.SortedMins[i-1]
            locMap[thisMin] = locMap[lastMin] + uint64(rawMinCounts[lastMin])
        }
    }

    // we now populate the raw kmers list, filing each according to its minimizer
    rawKmers := make(Kmers, numRawKmers, numRawKmers)
    matchChan := make(chan *Match, 500)
    go func() {
        for i := range rg.Matches {
            matchChan <- &rg.Matches[i]
        }
        close(matchChan)
    }()

    var wg = new(sync.WaitGroup)
    wg.Add(numCPU)
    respChan := make(chan ResponsePair)   // do not buffer without additional WaitGroup

    go func() {
        wg.Wait()
        close(respChan)
    }()

    for i := 0; i < numCPU; i++ {
        go rg.getMatchKmers(matchChan, respChan, wg)
    }

    for respPair := range respChan {
        rawKmers[locMap[respPair.MinInt]] = respPair.Kmer
        locMap[respPair.MinInt]++
    }


    return rawKmers
}

type ReducePair struct {
    LcaPtr *uint16
    Set Kmers
}

func (rg *RepeatGenome) uniqKmers(rawKmers Kmers) {
    // we first count the uniques so that we don't waste any capacity in the rg.Kmers slice
    var numUniqs uint64 = 0
    if len(rawKmers) > 0 {
        numUniqs++     // account for the first one, which is skipped
        for i := 1; i < len(rawKmers); i++ {
            lastKmerInt := *(*KmerInt)(unsafe.Pointer(&rawKmers[i-1]))
            kmerInt := *(*KmerInt)(unsafe.Pointer(&rawKmers[i]))
            if kmerInt != lastKmerInt {
                numUniqs++
            }
        }
    }

    rg.Kmers = make(Kmers, 0, numUniqs)

    pairChan := make(chan ReducePair)
    var wg = new(sync.WaitGroup)
    wg.Add(numCPU)

    go func() {
        wg.Wait()
        close(pairChan)
    }()

    for i := 0; i < numCPU; i++ {
        // this function finds the LCA of a set of kmers and updates the supplied LCA_ID pointer accordingly
        go func() {
            for pair := range pairChan {
                // grab the LCA of the first kmer in this set, and use it as our starting point
                lcaID := *(*uint16)(unsafe.Pointer(&pair.Set[0][8]))
                currLCA := rg.ClassTree.NodesByID[lcaID]

                // loop through the rest, updating currLCA
                for i := 1; i < len(pair.Set); i++ {
                    lcaID = *(*uint16)(unsafe.Pointer(&pair.Set[i][8]))
                    currLCA = rg.ClassTree.getLCA(currLCA, rg.ClassTree.NodesByID[lcaID])
                }

                *pair.LcaPtr = currLCA.ID
            }
            wg.Done()
        }()
    }

    numKmers := uint64(len(rawKmers))
    var start, end uint64 = 0, 0
    for end < numKmers {
        start = end
        end++
        for end < numKmers && *(*uint64)(unsafe.Pointer(&rawKmers[end-1])) == *(*uint64)(unsafe.Pointer(&rawKmers[end])) {
            end++
        }
        rg.Kmers = append(rg.Kmers, rawKmers[start])
        lcaPtr := (*uint16)(unsafe.Pointer(&rg.Kmers[len(rg.Kmers)-1][8]))
        pairChan <- ReducePair{lcaPtr, rawKmers[start : end]}
    }

}

func (rg *RepeatGenome) genKrakenLib() {

    fmt.Println("beginning first pass")
    numRawKmers, rawMinCounts := rg.krakenFirstPass()

    fmt.Println("expecting", comma(uint64(len(rawMinCounts))), "unique minimizers")
    fmt.Println("expecting", comma(numRawKmers), "non-unique kmers")
    fmt.Println()

    // populate a list of sorted minimizers
    rg.SortedMins = make(MinInts, 0, len(rawMinCounts))
    for minInt, _ := range rawMinCounts {
        rg.SortedMins = append(rg.SortedMins, minInt)
    }
    fmt.Println("RepeatGenome.SortedMins generated")
    fmt.Println()
    fmt.Println("Beginning sorting RepeatGenome.SortedMins")
    sort.Sort(rg.SortedMins)
    fmt.Println("RepeatGenome.SortedMins sorted")
    fmt.Println("len(rg.SortedMins):", len(rg.SortedMins))
    fmt.Println()
    
    fmt.Println("generating rawKmers and locMap")
    rawKmers := rg.getRawKmers(numRawKmers, rawMinCounts)
    runtime.GC()    // manual memory clear
    fmt.Println("raw kmers slice generated - len =", comma(uint64(len(rawKmers))))

    if debug {
        rawKmers.checkIntegrity()
    }

    fmt.Println("len(rawMinCounts):", comma(uint64(len(rawMinCounts))))

    fmt.Println("sorting rawKmers")
    rg.sortRawKmers(rawKmers, rawMinCounts)
    fmt.Println("rawKmers sorted")

    fmt.Println("generating RepeatGenome.Kmers")
    rg.uniqKmers(rawKmers)
    fmt.Println("done generating RepeatGenome.Kmers")

    fmt.Println(comma(uint64(len(rg.Kmers))), "unique kmers processed")

    rawKmers = nil
    runtime.GC()

    fmt.Println("generating RepeatGenome.MinCounts")
    rg.MinCounts = make(map[MinInt]uint32)
    for _, kmer := range rg.Kmers {
        kmerInt := *(*KmerInt)(unsafe.Pointer(&kmer))
        rg.MinCounts[kmerInt.Minimize()]++
    }
    fmt.Println("RepeatGenome.MinCounts generated")

    fmt.Println("generating RepeatGenome.MinOffsets")
    rg.populateMinOffsets()
    fmt.Println("rg.MinOffsets generated - done!")
}
