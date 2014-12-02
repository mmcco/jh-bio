/*
   Code involved with Kraken library generation.
*/

package repeatgenome

import (
    "fmt"
    "github.com/plsql/jh-bio/bioutils"
    "runtime"
    "sort"
    "sync"
    "unsafe"
)

/*
   Returns a pointer the the LCA of the two supplied ClassNode in the ClassTree.
   Currently uses a simple pointer-walking algorithm.
   This is a performance-critical function, and a location for potential optimization.
*/
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

/*
   Finds the ClassNode associated with a Kmer in the Kraken library.
   Part of the core Kraken read-classification functionality.
*/
func (rg *RepeatGenome) getKmerLCA(kmerInt KmerInt) *ClassNode {
    minimizer := kmerInt.Minimize()

    // simple binary search within the range of RepeatGenome.Kmers that has this kmer's minimizer
    i := rg.MinOffsets[minimizer]
    if i < 0 {
        return nil
    }
    j := i + int64(rg.MinCounts[minimizer])
    for i < j {
        x := (j + i) / 2
        thisKmerInt := rg.Kmers[x].Int()
        if thisKmerInt == kmerInt {
            return rg.ClassTree.NodesByID[rg.Kmers[x].ClassID()] // get kmer's LCA ID and return the corresponding ClassNode
        } else if thisKmerInt < kmerInt {
            i = x + 1
        } else {
            j = x
        }
    }

    return nil
}

/*
   Pushes all minimizers of all kmers of each Match in the supplied slice to the supplied chan.
*/
func (rg *RepeatGenome) getRawMins(matchChan chan *bioutils.Match, respChan chan MinInt, wg *sync.WaitGroup) {
    k_ := uint64(k)
    for match := range matchChan {
        start, end := match.SeqStart, match.SeqEnd
        matchSeq := rg.chroms[match.SeqName][match.SeqName][start:end]

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

            respChan <- TextSeq(matchSeq[i : i+k_]).kmerInt().Minimize()
        }
    }
    wg.Done()
}

/*
   The type returned by RepeatGenome.getMatchKmers(), which process
   raw kmers. The LCA contained in the Kmer value is not the Kmer's
   final LCA, but simply the ClassNode ID of the match this instance
   of the Kmer came from.
*/
type ResponsePair struct {
    Kmer   Kmer
    MinInt MinInt
}

/*
   Pushes all Kmers, along with their minimizers and the ClassNode ID
   of the match instance the kmer came from. As mentioned above, note
   that at this point the Kmer's LCA field contains the ClassNode ID
   of the match processed, not the Kmer's final LCA ID.
*/
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
   This is the first step of Kraken library generation. While it
   minimizes all repeat sequences genome, it only returns various
   counts, which are used for memory-allocation efficiency. It counts
   the number of raw (non-unique and possibly n-containing) kmers and
   the number of unique minimizers. It also returns a slice containing
   the number of times each minimizer appears in the genome.
*/
func (rg *RepeatGenome) krakenFirstPass() (numRawKmers uint64, numRawMins uint64, rawMinCounts []uint32) {
    matchChan := make(chan *bioutils.Match)
    go func() {
        for i := range rg.Matches {
            matchChan <- &rg.Matches[i]
        }
        close(matchChan)
    }()

    wg := new(sync.WaitGroup)
    wg.Add(numCPU)
    respChan := make(chan MinInt) // cannot buffer without additional WaitGroup

    go func() {
        wg.Wait()
        close(respChan)
    }()

    for i := 0; i < numCPU; i++ {
        go rg.getRawMins(matchChan, respChan, wg)
    }

    rawMinCounts = make([]uint32, minSliceSize, minSliceSize)
    for minInt := range respChan {
        if rawMinCounts[minInt] == 0 {
            numRawMins++
        }
        rawMinCounts[minInt]++
        numRawKmers++
    }

    wg.Wait()

    return numRawKmers, numRawMins, rawMinCounts
}

/*
   Sorts a slice of Kmers that are already sorted by minimizer. The
   resulting slice is therefore sorted primarily by minimizer and
   secondarily by kmer. This is done concurrently, using goroutines to
   sort each minimizer's kmers as a unique slice.
*/
func (rg *RepeatGenome) sortRawKmers(rawKmers Kmers, rawMinCounts []uint32) {
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
        end += uint64(rawMinCounts[minInt])
        kmersChan <- rawKmers[start:end]
        start = end
    }

    close(kmersChan)
    wg.Wait()
}

/*
   Returns a slice of all kmers in the genome's repeat sequences,
   unsorted and non-unique.
*/
func (rg *RepeatGenome) getRawKmers(numRawKmers uint64, rawMinCounts []uint32) Kmers {
    // we need a map telling us the next place to insert a kmer associated with any given minimizer
    // we simply increment a minimizer's offset when we insert a kmer associated with it
    locMap := make([]int64, minSliceSize, minSliceSize)

    for i := 1; i < len(rg.SortedMins); i++ {
        thisMin := rg.SortedMins[i]
        lastMin := rg.SortedMins[i-1]
        locMap[thisMin] = locMap[lastMin] + int64(rawMinCounts[lastMin])
    }
    // could set all unoccupied locMap indices to -1 here for debugging

    // we now populate the raw kmers list, filing each according to its minimizer
    rawKmers := make(Kmers, numRawKmers)
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

    for respPair := range respChan {
        rawKmers[locMap[respPair.MinInt]] = respPair.Kmer
        locMap[respPair.MinInt]++
    }

    return rawKmers
}

type ReducePair struct {
    LcaPtr *ClassID
    Set    Kmers
}

func (rg *RepeatGenome) uniqKmers(rawKmers Kmers) {
    // we first count the uniques so that we don't waste any capacity
    // in the rg.Kmers slice
    var numUniqs uint64 = 0
    if len(rawKmers) > 0 {
        numUniqs++ // account for the first one, which is skipped
        for i := 1; i < len(rawKmers); i++ {
            lastKmerInt := rawKmers[i-1].Int()
            kmerInt := rawKmers[i].Int()
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
        // this function finds the LCA of a set of kmers and updates
        // the supplied ClassID pointer accordingly
        go func() {
            for pair := range pairChan {
                // grab the LCA of the first kmer in this set, and use
                // it as our starting point
                classID := pair.Set[0].ClassID()
                currLCA := rg.ClassTree.NodesByID[classID]

                // loop through the rest, updating currLCA
                for i := 1; i < len(pair.Set); i++ {
                    classID = pair.Set[i].ClassID()
                    currLCA = rg.ClassTree.getLCA(currLCA, rg.ClassTree.NodesByID[classID])
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
        for end < numKmers && rawKmers[end-1].Int() == rawKmers[end].Int() {
            end++
        }
        rg.Kmers = append(rg.Kmers, rawKmers[start])
        classPtr := (*ClassID)(unsafe.Pointer(&rg.Kmers[len(rg.Kmers)-1][8]))
        pairChan <- ReducePair{classPtr, rawKmers[start:end]}
    }

}

func (rg *RepeatGenome) genKrakenLib() {

    fmt.Println()
    fmt.Println("beginning first pass")
    numRawKmers, numRawMins, rawMinCounts := rg.krakenFirstPass()
    fmt.Println("first pass complete")

    //fmt.Println("len(rawMinCounts):", comma(uint64(len(rawMinCounts))))

    fmt.Println("expecting", comma(int(numRawMins)), "unique minimizers")
    fmt.Println("expecting", comma(int(numRawKmers)), "non-unique kmers")

    // populate a list of sorted minimizers
    fmt.Println()
    fmt.Println("generating RepeatGenome.SortedMins")
    rg.SortedMins = make(MinInts, 0, numRawMins)
    for minInt, cnt := range rawMinCounts {
        if cnt > 0 {
            rg.SortedMins = append(rg.SortedMins, MinInt(minInt))
        }
    }
    fmt.Println("RepeatGenome.SortedMins generated")
    fmt.Println()
    fmt.Println("Beginning sorting RepeatGenome.SortedMins")
    sort.Sort(rg.SortedMins)
    fmt.Println("RepeatGenome.SortedMins sorted")
    fmt.Println("len(rg.SortedMins):", comma(len(rg.SortedMins)))

    fmt.Println()
    fmt.Println("generating rawKmers and locMap")
    rawKmers := rg.getRawKmers(numRawKmers, rawMinCounts)
    runtime.GC() // manual memory clear
    fmt.Println("raw kmers slice generated")
    fmt.Println("len(rawKmers) =", comma(len(rawKmers)))

    if debug {
        rawKmers.checkIntegrity()
    }

    fmt.Println("sorting rawKmers")
    rg.sortRawKmers(rawKmers, rawMinCounts)
    fmt.Println("rawKmers sorted")

    fmt.Println("generating RepeatGenome.Kmers")
    rg.uniqKmers(rawKmers)
    fmt.Println("done generating RepeatGenome.Kmers")

    fmt.Println(comma(len(rg.Kmers)), "unique kmers processed")

    rawKmers = nil
    runtime.GC()

    fmt.Println("generating RepeatGenome.MinCounts")
    rg.MinCounts = make([]uint32, minSliceSize, minSliceSize)
    for _, kmer := range rg.Kmers {
        rg.MinCounts[kmer.Int().Minimize()]++
    }
    fmt.Println("RepeatGenome.MinCounts generated")

    fmt.Println("generating RepeatGenome.MinOffsets")
    rg.populateMinOffsets()
    fmt.Println("rg.MinOffsets generated - done!")
}
