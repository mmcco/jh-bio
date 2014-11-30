/*
   General debugging functions.
   Most effort is put into checking the integrity and correctness of data.
*/

package repeatgenome

import (
    "fmt"
    "os"
    "reflect"
    "unsafe"
)

func (rg *RepeatGenome) printSample(numMins, numKmers int) {
    for i := 0; i < len(rg.SortedMins); i += len(rg.SortedMins) / numMins {
        thisMin := rg.SortedMins[i]
        thisMin.print()
        fmt.Println()
        theseKmers := rg.getMinsKmers(thisMin)
        if len(theseKmers) > numKmers {
            theseKmers = theseKmers[:numKmers]
        }
        for _, kmer := range theseKmers {
            fmt.Print("\t")
            kmer.Int().print()
            fmt.Printf("\t%s\n", rg.ClassTree.NodesByID[kmer.ClassID()].Name)
        }
    }
}

func (rg *RepeatGenome) RunDebugTests() {
    fmt.Println()
    rg.checkIntegrity()

    fmt.Println()
    for k, v := range rg.chroms {
        for k_, v_ := range v {
            fmt.Printf("chrom: %s\tseq: %s\t%s...%s\n", k, k_, v_[:20], v_[len(v_)-20:])
        }
    }
    fmt.Println()

    fmt.Println("number of chromosomes parsed:", len(rg.chroms))
    fmt.Println()

    fmt.Println("total number of bases in genome:", rg.Size())

    rg.ClassTree.PrintBranches()
    fmt.Println()
    fmt.Println("number of ClassNodes:", len(rg.ClassTree.ClassNodes))
    fmt.Println()

    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['DNA/TcMar-Mariner'], rg.ClassTree.ClassNodes['DNA/TcMar-Tc1']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["DNA/TcMar-Mariner"], rg.ClassTree.ClassNodes["DNA/TcMar-Tc1"]).Name)
    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['ARTEFACT'], rg.ClassTree.ClassNodes['DNA/TcMar-Tc1']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["ARTEFACT"], rg.ClassTree.ClassNodes["DNA/TcMar-Tc1"]).Name)
    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['LINE/LOA'], rg.ClassTree.ClassNodes['root']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["LINE/LOA"], rg.ClassTree.ClassNodes["root"]).Name)
    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['Simple_repeat/(T)n'], rg.ClassTree.ClassNodes['Simple_repeat/(T)n']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["Simple_repeat/(T)n"], rg.ClassTree.ClassNodes["Simple_repeat/(T)n"]).Name)
    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['LTR/Gypsy/MICROPIA_I-int'], rg.ClassTree.ClassNodes['LTR/Gypsy']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["LTR/Gypsy/MICROPIA_I-int"], rg.ClassTree.ClassNodes["LTR/Gypsy"]).Name)
    fmt.Println("rg.ClassTree.getLCA(rg.ClassTree.ClassNodes['LTR/Gypsy'], rg.ClassTree.ClassNodes['LTR/Gypsy/MICROPIA_I-int']).Name:", rg.ClassTree.getLCA(rg.ClassTree.ClassNodes["LTR/Gypsy"], rg.ClassTree.ClassNodes["LTR/Gypsy/MICROPIA_I-int"]).Name)
    fmt.Println()

    fmt.Println("min(5, 7):", min(5, 7))
    fmt.Println("max64(int64(5), int64(7)):", max64(int64(5), int64(7)))
    fmt.Println()

    testSeq := TextSeq("atgtttgtgtttttcataaagacgaaagatg")
    thisKmer := testSeq.kmerInt()
    thisMin := thisKmer.Minimize()
    fmt.Println("minimizer of 'tgctcctgtcatgcatacgcaggtcatgcat': ")
    thisMin.print()
    fmt.Println()

    fmt.Print("revComp of ")
    thisKmer.print()
    fmt.Println(":")
    thisKmer.revComp().print()
    fmt.Println()
    fmt.Printf("Kmer struct size: %d\n", unsafe.Sizeof(Kmer{}))
    fmt.Println()

    fmt.Printf("GetReadClassChan() using %d CPUs\n", numCPU)
    fmt.Println()

    rg.printSample(5, 5)
}

func DebugSeq() {
    ts := TextSeq("tgaatatacgtagctctagctagcgcttatatg")
    fmt.Println("ts:", ts)
    s := ts.Seq()
    fmt.Print("s: ")
    s.Print()
    fmt.Println()
    fmt.Println("s[4]: ", s.GetBase(4))
    fmt.Println("s[-1]: ", s.GetBase(s.Len-1))
    fmt.Print("s[:7]: ")
    s.Subseq(0, 7).Print()
    fmt.Println()
    fmt.Print("s[7:] ")
    s.Subseq(7, s.Len).Print()
    fmt.Println()
}

func (rg *RepeatGenome) checkIntegrity() {

    for i, match := range rg.Matches {
        if uint64(i) != match.ID {
            fmt.Println("ERROR: match at index", i, "has an ID that does not match its index")
            os.Exit(1)
        }
    }

    for i := 1; i < len(rg.SortedMins); i++ {
        if rg.SortedMins[i] <= rg.SortedMins[i-1] {
            fmt.Println("ERROR: rg.SortedMins not sorted")
            os.Exit(1)
        }
        if rg.SortedMins[i] == rg.SortedMins[i-1] {
            fmt.Println("ERROR: rg.SortedMins contains multiple copies of a minimizer")
            os.Exit(1)
        }
    }
    fmt.Println("rg.SortedMins is indeed sorted and unique")

    mcLen := 0
    for _, cnt := range rg.MinCounts {
        if cnt > 0 {
            mcLen++
        }
    }

    ofstLen := 0
    for _, offset := range rg.MinOffsets {
        if offset >= 0 {
            ofstLen++
        }
    }

    if len(rg.SortedMins) != mcLen || mcLen != ofstLen {
        fmt.Println("ERROR: rg.SortedMins, rg.MinCounts, and rg.MinOffsets of inconsistent size")
        os.Exit(1)
    } else {
        fmt.Println("rg.SortedMins, rg.MinCounts, and rg.MinOffsets of compatible size")
    }

    for _, minInt := range rg.SortedMins {
        cnt := rg.MinCounts[minInt]
        if cnt < 1 {
            fmt.Println("minimizer in rg.SortedMins absent from rg.MinCounts")
            os.Exit(1)
        }
        offset := rg.MinOffsets[minInt]
        if offset < 0 {
            fmt.Println("minimizer in rg.SortedMins absent from rg.MinOffsets")
            os.Exit(1)
        }
    }
    fmt.Println("all minimizers in rg.SortedMins present in rg.MinCounts and rg.MinOffsets")

    for _, kmer := range rg.Kmers {
        lca := rg.getKmerLCA(kmer.Int())
        if lca == nil {
            fmt.Print("ERROR: kmer ")
            kmer.Int().print()
            fmt.Println("not found by rg.getKmer()")
            os.Exit(1)
        }
    }
}

/*
   Checks whether every index in rawKmers was set when it was populated.
   rawKmers is initialized as a list of zero-kmers.
   If there was a mismatch in the supplied minimizer counts and the actual minimizer counts, at least one kmer will almost certainly be left unassigned.
   Because no reference repeats map to root (whose ClassID is 0), no rawKmer index should be a zero kmer.
*/
func (rawKmers Kmers) checkIntegrity() {
    var zeroKmer Kmer // works because no raw kmer should have an LCA of zero, which represents root
    for i, kmer := range rawKmers {
        if reflect.DeepEqual(kmer, zeroKmer) {
            fmt.Println("rawKmer at index", i, "not set")
            os.Exit(1)
        }
    }
}
