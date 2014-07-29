package repeatgenome

import (
    "fmt"
    "os"
    "unsafe"
)

func (rg *RepeatGenome) printSample(numMins, numKmers int) {
    for i := 0; i < len(rg.SortedMins); i += len(rg.SortedMins) / numMins {
        thisMin := rg.SortedMins[i]
        thisMin.print(); fmt.Println()
        theseKmers := rg.getMinsKmers(thisMin)
        if len(theseKmers) > numKmers {
            theseKmers = theseKmers[:numKmers]
        }
        for _, kmer := range theseKmers {
            fmt.Print("\t")
            (*(*KmerInt)(unsafe.Pointer(&kmer))).print()
            lca_ID := *(*uint16)(unsafe.Pointer(&kmer[8]))
            fmt.Printf("\t%s\n", rg.ClassTree.NodesByID[lca_ID].Name)
        }
    }
}

func (rg *RepeatGenome) RunDebugTests() {
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

    fmt.Print("revComp of "); thisKmer.print(); fmt.Println(":")
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
    fmt.Println("s[-1]: ", s.GetBase(s.Len - 1))
    fmt.Print("s[:7]: ")
    s.Subseq(0, 7).Print()
    fmt.Println()
    fmt.Print("s[7:] ")
    s.Subseq(7, s.Len).Print()
    fmt.Println()
}

func (rg *RepeatGenome) testDataStructIntegrity() {

    for i := 1; i < len(rg.SortedMins); i++ {
        if rg.SortedMins[i] > rg.SortedMins[i-1] {
            fmt.Println("ERROR: rg.SortedMins not sorted")
            os.Exit(1)
        }
    }
    fmt.Println("rg.SortedMins is indeed sorted")

    if len(rg.SortedMins) != len(rg.MinCounts) || len(rg.MinCounts) != len(rg.MinOffsets) {
        fmt.Println("ERROR: rg.SortedMins, rg.MinCounts, and rg.MinOffsets of inconsistent size")
    } else {
        fmt.Println("rg.SortedMins, rg.MinCounts, and rg.MinOffsets of compatible size")
    }
}
