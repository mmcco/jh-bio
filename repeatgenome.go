package repeatgenome

/*
   WARNING!!! This program is currently under development and may be buggy or broken.

   A barebones (at the moment) Go script for parsing and minimizing RepeatMasker output files alongside FASTA reference genomes.

   This script expects there to be a subdirectory of the current directory named after the reference genome used (e.g. "dm3") that contains the following files:
       * a RepeatMasker library containing:
           - the match library (e.g. "dm3.fa.out")
           - the alignment information (e.g. "dm3.fa.align")
       * one or more reference genome files in FASTA format with the suffix ".fa"

   Premature commenting is the root of all evil, and I have sinned. Please read comments skeptically - they will eventually be audited.

   Should probably make a file solely for type defs.

   Error handling should be updated with a custom ParseError type and removal fo checkError()

   For portability's sake, the flags should be used as args to Generate() rather than globals.

   Minimizers could be stored as uint32's.

   Should switch from strings to byte slices.

   The concurrent read-kmer generator could be reintroduced using a select statement.

   Should probably restrict activity of chans with directionals

   Kmer counting should be re-added eventually - it's currently excluded for performance reasons because we aren't using it.

   We should test a version that doesn't cache minimizers, as that seems to be a needless bottleneck. It could also be conditional on the number of CPUs available.

   Minimizers are currently not written to file in any order. This is for memory efficiency, and can be changed in the future.

   I should probably change some variable names, like repeatGenome, to less verbose variants, and use more aliases.

   Ranges should be changed to use actual values instead of indexes.

   Slice sizes should be specified in the make() call when the size is known.

   seqToInt and revCompToInt need minor cleanup and a potential name-change.

   All sequences containing Ns are currently ignored.

   We should consider taking end minimizers once the code base is more mature.

   We should also review how to deal with m <= len(match) < k.

   For caching efficiency, we should change the minimizer data structure to a map-indexed 1D slice of Kmers (not *Kmers). (This technique originated in Kraken.)

   Int sizes should be reviewed for memory efficiency.

   The sole command line argument is the name of the reference genome (e.g. "dm3").
*/

import (
    "bytes"
    "fmt"
    "io/ioutil"
    "log"
    //"mapset"
    "os"
    "runtime"
    "runtime/pprof"
    "sort"
    "strconv"
    "strings"
    "sync"
    "unsafe"
)

type Flags struct {
    Debug       bool
    CPUProfile  bool
    MemProfile  bool
    Minimize    bool
    WriteKraken bool
    WriteJSON   bool
}

// Match.SW_Score - Smith-Waterman score, describing the likeness to the repeat reference sequence
// Match.PercDiv
// Match.PercDel
// Match.PercIns
// Match.SeqName -  the reference genome file this match came from (typically the chromosome)
// Match.SeqStart -  the starting index (inclusive) in the reference genome
// Match.SeqEnd -  the ending index (exclusive) in the reference genome
// Match.SeqRemains
// Match.IsRevComp -  the match may be for the complement of the reference sequence
// Match.RepeatClass - the repeat's full ancestry, including its repeat class and repeat name (which are listed separately in the RepeatMasker output file)
// Match.RepeatStart-  the starting index in the repeat consensus sequence
// Match.RepeatEnd -  the ending sequence (exclusive) in the repeat consensus sequence
// Match.RepeatRemains
// Match.InsertionID - a numerical ID for the repeat type (starts at 1)

//     below are not in parsed data file
// Match.RepeatName - simply repeatClass concatenated - used for printing and map indexing
// Match.ClassNode
// Match.Repeat
type Match struct {
    SW_Score    int32
    PercDiv     float64
    PercDel     float64
    PercIns     float64
    SeqName     string
    SeqStart    uint64
    SeqEnd      uint64
    SeqRemains  uint64
    IsRevComp   bool
    RepeatClass []string
    // in weird cases, RepeatStart can be negative, so they must be signed
    RepeatStart   int64
    RepeatEnd     int64
    RepeatRemains int64
    InsertionID   uint64

    // these are generated, not parsed
    RepeatName string
    ClassNode  *ClassNode
    Repeat     *Repeat
}

type RepeatGenome struct {
    Name       string
    Flags Flags
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    chroms        map[string](map[string]TextSeq)
    K             uint8
    M             uint8
    Kmers         Kmers
    // stores the offset of each minimizer's first kmer in RepeatGenome.Kmers
    OffsetsToMin  map[uint64]uint64
    // stores the number of kmers that each minimizer is associated with
    MinCounts     map[uint64]uint64
    SortedMins    Uint64Slice
    Matches       Matches
    ClassTree     ClassTree
    Repeats       Repeats
    RepeatMap     map[string]*Repeat
    memProfFile   *os.File
}

type ClassTree struct {
    // maps all class names to a pointer to their node struct
    // we must use pointers because of this foible in golang: https://code.google.com/p/go/issues/detail?id=3117
    // if we didn't use pointers, we'd have to completely reassign the struct when adding parents, children, etc.
    ClassNodes map[string](*ClassNode)
    NodesByID  []*ClassNode
    // a pointer to the the class tree's root, used for recursive descent etc.
    // we explicitly create the root (because RepeatMatcher doesn't)
    Root *ClassNode
}

// Used to differentiate sequence representations with one base per byte (type TextSeq) from those with four bases per byte (type Seq).
type TextSeq []byte

// can store a kmer where k <= 32
// the value of k is not stored in the struct, but rather in the RepeatGenome, for memory efficiency
// first eight bits are the int representation of the sequence
// the last two are the LCA ID
type Kmer [10]byte

// as with the Kmer type, each base is represented by two bits
// any excess bits are the first bits of the first byte (seq is right-justified)
// remember that len(Seq.Bases) is not the actual number of bases, but rather the number of bytes necessary to represent them
type Seq struct {
    Bases []byte
    Len   uint64
}

type Seqs []Seq

type repeatLoc struct {
    SeqName  string
    StartInd uint64
    EndInd   uint64
}

type Repeat struct {
    // assigned in simple incremented order starting from 1
    // they are therefore not compatible across genomes
    // we give root ID = 0
    ID uint64
    // a list containing the repeat's ancestry path, from top down
    // root is implicit, and is therefore excluded from the list
    ClassList []string
    ClassNode *ClassNode
    Name  string
    Instances []*Match
    Locations []repeatLoc
}

type ClassNode struct {
    Name     string
    ID       uint16
    Class    []string
    Parent   *ClassNode
    Children []*ClassNode
    IsRoot   bool
    Repeat   *Repeat
}

// type synonyms, necessary to implement interfaces (e.g. sort) and methods
type Kmers []Kmer
type PKmers []*Kmer
type MinMap map[uint64]PKmers
type Repeats []Repeat
type Matches []Match

//type Chroms map[string](map[string]TextSeq)

type MinPair struct {
    KmerInt   Kmer
    Minimizer uint64
}

type ThreadResponse struct {
    KmerInt   uint64
    Minimizer uint64
    Relative *ClassNode
}

type MinCache struct {
    sync.Mutex
    Cache map[uint64]uint64
}

func parseMatches(genomeName string) (error, Matches) {
    // "my_genome_name"  ->  "my_genome_name/my_genome_name.fa.out"
    filepath := strings.Join([]string{genomeName, "/", genomeName, ".fa.out"}, "")
    err, matchLines := fileLines(filepath)
    if err != nil {
        return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
    }
    // drop header
    matchLines = matchLines[3:]

    var matches Matches
    var sw_Score int64

    for _, matchLine := range matchLines {
        rawVals := strings.Fields(string(matchLine))
        if len(rawVals) != 15 {
            return ParseError{"repeatgenome.parseMatches()",
                           filepath,
                           fmt.Errorf("supplied match line is not 15 fields long (has %d fields and length %d):\n", len(rawVals), len(matchLine))},
                   nil
        }
        var match Match
        match.IsRevComp = rawVals[8] == "C"

        // remove enclosing parentheses
        // !!! in the future, checkes to ensure that the parentheses exist should be added
        // !!! it would also be sensible to check that rawVals[8] is either "C" or "+"
        rawVals[7] = rawVals[7][1 : len(rawVals[7])-1]
        if match.IsRevComp {
            rawVals[11] = rawVals[11][1 : len(rawVals[11])-1]
        } else {
            rawVals[13] = rawVals[13][1 : len(rawVals[13])-1]
        }

        // everything in this block is just vanilla trimming, converting, and error checking
        sw_Score, err = strconv.ParseInt(rawVals[0], 10, 32)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.SW_Score = int32(sw_Score)
        match.PercDiv, err = strconv.ParseFloat(rawVals[1], 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.PercDel, err = strconv.ParseFloat(rawVals[2], 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.PercIns, err = strconv.ParseFloat(rawVals[3], 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.SeqName = strings.TrimSpace(rawVals[4])
        match.SeqStart, err = strconv.ParseUint(rawVals[5], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.SeqEnd, err = strconv.ParseUint(rawVals[6], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.SeqRemains, err = strconv.ParseUint(rawVals[7], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        // match.IsComplement, rawVals[8], moved above
        match.RepeatClass = append(strings.Split(strings.TrimSpace(rawVals[10]), "/"), strings.TrimSpace(rawVals[9]))
        match.RepeatStart, err = strconv.ParseInt(rawVals[11], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.RepeatEnd, err = strconv.ParseInt(rawVals[12], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.RepeatRemains, err = strconv.ParseInt(rawVals[13], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }
        match.InsertionID, err = strconv.ParseUint(rawVals[14], 10, 64)
        if err != nil {
            return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
        }

        // necessary swaps to convert reverse complement repeat indexes to positive-strand indexes
        if match.IsRevComp {
            match.RepeatStart = match.RepeatRemains
            match.RepeatEnd = match.RepeatStart
            match.RepeatRemains = match.RepeatRemains + (match.RepeatEnd - match.RepeatStart)
        }

        // decrement match.SeqStart and match.RepeatStart so that they work from a start index of 0 rather than 1
        // that way, we can use them without modification in slices
        match.SeqStart--
        match.RepeatStart--

        // "Other" and "Unknown" classes are heirarchically meaningless and really just mean "root", so we remove them
        if match.RepeatClass[0] == "Other" || match.RepeatClass[0] == "Unknown" {
            match.RepeatClass = match.RepeatClass[1:]
        }

        match.RepeatName = strings.Join(match.RepeatClass, "/")

        matches = append(matches, match)
    }
    return nil, matches
}

func parseGenome(genomeName string) (error, map[string](map[string]TextSeq)) {
    chromFileInfos, err := ioutil.ReadDir(genomeName)
    if err != nil {
        return IOError{"repeatgenome.parseGenome()", err}, nil
    }
    warned := false
    chroms := make(map[string](map[string]TextSeq))
    // used below to store the two keys for RepeatGenome.chroms
    for i := range chromFileInfos {
        // "my_genome_name", "my_chrom_name"  ->  "my_genome_name/my_chrom_name"
        chromFilename := chromFileInfos[i].Name()
        chromFilepath := strings.Join([]string{genomeName, chromFilename}, "/")
        // process the ref genome files (*.fa), not the repeat ref files (*.fa.out and *.fa.align) or anything else
        if strings.HasSuffix(chromFilepath, ".fa") {
            err, seqLines := fileLines(chromFilepath)
            if err != nil {
                // error should already be IOError
                return err, nil
            }

            // maps each sequence name in this chrom to a slice of its sequence's lines
            // the list is concatenated at the end for efficiency's sake
            seqMap := make(map[string][][]byte)
            numLines := uint64(len(seqLines))
            var seqName string = ""        // forward initialization necessary
            var i uint64
            for i = 0; i < numLines; i++ {
                seqLine := bytes.TrimSpace(seqLines[i])
                if seqLine[0] == byte('>') {
                    seqName = string(bytes.TrimSpace(seqLine[1:]))
                    if !warned && seqName != chromFilename[:len(chromFilename)-3] {
                        fmt.Println("WARNING: reference genome is two-dimensional, containing sequences not named after their chromosome.")
                        fmt.Println("Because RepeatMasker supplied only one-dimensional indexing, this may cause unexpected behavior or program failure.")
                        fmt.Println("seqName:", seqName, "\tlen(seqName):", len(seqName))
                        fmt.Println("chrom name:", chromFilename[:len(chromFilename)-3], "\tlen(chrom name):", len(chromFilename)-3)
                        warned = true
                    }
                } else {
                    if seqName == "" {
                        return ParseError{"repeatgenome.parseGenome", chromFilepath, fmt.Errorf("Empty or missing sequence name")}, nil
                    }
                    seqMap[seqName] = append(seqMap[seqName], seqLine)
                }
            }
            // finally, we insert this map into the full map
            chromName := chromFilepath[len(genomeName)+1 : len(chromFilepath)-3]
            chroms[chromName] = make(map[string]TextSeq)
            for seqName, seqLines := range seqMap {
                chroms[chromName][seqName] = TextSeq(bytes.ToLower(bytes.Join(seqLines, []byte{})))
            }
        }
    }
    return nil, chroms
}

func Generate(genomeName string, k, m uint8, rgFlags Flags) (error, *RepeatGenome) {
    var err error
    // we popoulate the RepeatGenome mostly with helper functions
    // we should consider whether it makes more sense for them to alter the object directly, than to return their results
    rg := new(RepeatGenome)
    rg.Name = genomeName
    rg.Flags = rgFlags

    if rg.Flags.MemProfile {
        os.Mkdir("profiles", os.ModeDir)
        rg.memProfFile, err = os.Create("profiles/" + rg.Name + ".memprof")
        if err != nil {
            return IOError{"RepeatGenome.getKrakenSlice()", err}, nil
        }
        pprof.WriteHeapProfile(rg.memProfFile)
        defer rg.memProfFile.Close()
    }

    err, rg.chroms = parseGenome(genomeName)
    if err != nil {
        return err, nil
    }
    err, rg.Matches = parseMatches(genomeName)
    if err != nil {
        return err, nil
    }
    rg.getRepeats()
    rg.getClassTree()
    rg.K = k
    rg.M = m

    if rg.Flags.Minimize {
        // calling the parallel minimizer and writing the result
        rg.getKrakenSlice()
    }

    if rg.Flags.WriteJSON {
        rg.WriteClassJSON(false, false)
    }

    if rg.Flags.Debug {
        rg.RunDebugTests()
    }

    return nil, rg
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

    fmt.Println()
    fmt.Println("min(5, 7):", min(5, 7))
    fmt.Println("max64(int64(5), int64(7)):", max64(int64(5), int64(7)))
    fmt.Println()

    testSeq := TextSeq("atgtttgtgtttttcataaagacgaaagatg")
    offset, thisMin := getMinimizer(seqToInt(testSeq), uint8(len(testSeq)), 15)
    fmt.Println("getMinimizer('tgctcctgtcatgcatacgcaggtcatgcat', 15) offset :", offset)
    printSeqInt(thisMin, 15)
    fmt.Println()

    fmt.Printf("Kmer struct size: %d\n", unsafe.Sizeof(Kmer{}))
}

func (rg *RepeatGenome) getRepeats() {
    // we now populate a list of unique repeat types
    // repeats are stored in the below slice, indexed by their ID
    // we first determine the necessary size of the slice - we can't use append because matches are not sorted by repeatID

    rg.RepeatMap = make(map[string]*Repeat)

    // DON'T use the second field of the range - this causes the Match struct to be copied
    // creating an alias struct (match := rg.Matches[i]) of type Match rather than *Match causes the repeat.Instance item to point to a copy, not the original Match struct
    for i := range rg.Matches {
        match := &rg.Matches[i]
        // don't bother overwriting
        if repeat, exists := rg.RepeatMap[match.RepeatName]; exists {
            repeat.Instances = append(repeat.Instances, match)
            match.Repeat = repeat
        } else {
            var repeat Repeat
            repeat.ID = uint64(len(rg.Repeats))
            repeat.ClassList = match.RepeatClass
            repeat.Name = match.RepeatName
            repeat.Instances = []*Match{match}
            repeat.Locations = append(repeat.Locations, repeatLoc{match.SeqName, match.SeqStart, match.SeqEnd})

            rg.Repeats = append(rg.Repeats, repeat)
            rg.RepeatMap[repeat.Name] = &repeat

            match.Repeat = &repeat
        }
    }
}

func (rg *RepeatGenome) getClassTree() {
    // mapping to pointers allows us to make references (i.e. pointers) to values
    tree := &rg.ClassTree
    tree.ClassNodes = make(map[string](*ClassNode))
    tree.NodesByID = make([]*ClassNode, 0)
    // would be prettier if expanded
    tree.Root = &ClassNode{"root", 0, []string{"root"}, nil, nil, true, nil}
    tree.ClassNodes["root"] = tree.Root
    tree.NodesByID = append(tree.NodesByID, tree.Root)

    for _, repeat := range rg.Repeats {
        // process every heirarchy level (e.g. for "DNA/LINE/TiGGER", process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
        for j := 1; j <= len(repeat.ClassList); j++ {
            thisClass := repeat.ClassList[:j]
            thisClassName := strings.Join(thisClass, "/")
            _, keyExists := tree.ClassNodes[thisClassName]
            if !keyExists {
                if len(tree.NodesByID) > 65534 {
                    panic("RepeatGenome.getClassTree(): more than 65,536 class nodes - ID is overflowed")
                }
                classNode := new(ClassNode)
                classNode.Name = thisClassName
                classNode.ID = uint16(len(tree.NodesByID))
                classNode.Class = thisClass
                classNode.IsRoot = false
                if repeat, exists := rg.RepeatMap[thisClassName]; exists {
                    classNode.Repeat = repeat
                }

                tree.ClassNodes[thisClassName] = classNode
                tree.NodesByID = append(tree.NodesByID, classNode)
                // first case handles primary classes, as root is implicit and not listed in thisClass
                if j == 1 {
                    classNode.Parent = tree.Root
                } else {
                    classNode.Parent = tree.ClassNodes[strings.Join(thisClass[:len(thisClass)-1], "/")]
                }
                if classNode.Parent.Children == nil {
                    classNode.Parent.Children = make([]*ClassNode, 0)
                }
                classNode.Parent.Children = append(classNode.Parent.Children, tree.ClassNodes[thisClassName])
            }
        }
    }

    // MUST NOT USE RANGE - the struct will be copied!
    for i := 0; i < len(rg.Repeats); i++ {
        repeat := &rg.Repeats[i]
        repeat.ClassNode = tree.ClassNodes[repeat.Name]

        if repeat.ClassNode == nil {
            fmt.Println(repeat.Name)
            log.Fatal("getClassTree(): nil Repeat.ClassNode")
        }
    }

    // MUST NOT USE RANGE - the struct will be copied!
    for i := 0; i < len(rg.Matches); i++ {
        match := &rg.Matches[i]
        match.ClassNode = tree.ClassNodes[match.RepeatName]

        if match.ClassNode == nil {
            fmt.Println(match.RepeatName)
            log.Fatal("getClassTree(): nil Match.ClassNode")
        }
    }
}

// returns ancestry list, from self to root
func (cn *ClassNode) getAncestry() []*ClassNode {
    // condition to prevent nil dereference of Root.Parent
    if cn.IsRoot {
        return []*ClassNode{}
    }
    ancestry := []*ClassNode{cn.Parent}
    for !ancestry[len(ancestry)-1].IsRoot {
        ancestry = append(ancestry, ancestry[len(ancestry)-1].Parent)
    }
    return ancestry
}

func (classTree *ClassTree) getLCA(cnA, cnB *ClassNode) *ClassNode {
    ancestryA := cnA.getAncestry()
    cnBWalker := cnB
    for cnBWalker != classTree.Root {
        for i := 0; i < len(ancestryA); i++ {
            if cnBWalker == ancestryA[i] {
                return cnBWalker
            }
        }
        cnBWalker = cnBWalker.Parent
    }
    // necessary for compilation - Root should be in the ancestry paths
    return classTree.Root
}

// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func (rg *RepeatGenome) minimizeThread(matchStart, matchEnd uint64, c chan ThreadResponse) {
    var match *Match
    k := rg.K
    k_ := uint64(k)
    m := rg.M
    var matchLen, kmerInt, thisMin uint64

    for i := matchStart; i < matchEnd; i++ {
        match = &rg.Matches[i]
        seq := rg.chroms[match.SeqName][match.SeqName]
        matchLen = match.SeqEnd - match.SeqStart
        // for now, we will ignore matches too short to be traditionally minimized
        if matchLen < k_ {
            continue
        }

    KmerLoop:
        for j := match.SeqStart; j <= match.SeqEnd-k_; j++ {
            // we begin by skipping any kmers containing n's
            // we start checking from the end for maximum skipping efficiency
            for x := k_ - 1; x >= 0; x-- {
                if seq[j+x] == byte('n') {
                    j += x
                    continue KmerLoop
                }
                // prevents negative overflow - a bit of a hack, but seqs can be big, so we need uint64's capacity
                if x == 0 { break }
            }

            kmerInt = seqToInt(seq[j : j+k_])
            // make the sequence strand-agnostic
            kmerInt = minU64(kmerInt, intRevComp(kmerInt, k))

            _, thisMin = getMinimizer(kmerInt, k, m)
            if match.ClassNode == nil {
                fmt.Println("current match's repeat:", match.RepeatName)
                panic("minimizeThread(): match has nil ClassNode")
            }
            c <- ThreadResponse{kmerInt, thisMin, match.ClassNode}
        }
    }

    close(c)
}

func (rg *RepeatGenome) CoordMinimizerThreads() {
}

func (rg *RepeatGenome) getKrakenSlice() error {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    numCPU := runtime.NumCPU()
    if rg.Flags.Debug {
        fmt.Printf("getKrakenSlice() using %d CPUs\n", numCPU)
    }
    runtime.GOMAXPROCS(numCPU)
    var mStart, mEnd uint64

    if rg.Flags.Debug {
        numKmers := rg.numKmers()
        fmt.Printf("expecting >= %d million kmers\n", numKmers/1000000)
    }

    var threadChans [](chan ThreadResponse)
    for i := 0; i < numCPU; i++ {
        mStart = uint64(i * len(rg.Matches) / numCPU)
        mEnd = uint64((i + 1) * len(rg.Matches) / numCPU)
        c := make(chan ThreadResponse, 1000)
        threadChans = append(threadChans, c)
        go rg.minimizeThread(mStart, mEnd, c)
    }

    // below is the atomic section of minimizing, which is parallel
    // this seems to be the rate-limiting section, as threads use only ~6-7 CPU-equivalents
    // it should therefore be optimized before other sections
    var kmer *Kmer
    var kmersProcessed, kmerInt, minimizer uint64
    var lca_ID uint16
    var lca, relative *ClassNode
    var exists bool
    minCache := make(map[uint64]uint64)
    kmerMap := make(map[uint64]*Kmer)

    for response := range mergeThreadResp(threadChans) {
        if kmersProcessed%5000000 == 0 {
            fmt.Println(comma(kmersProcessed/1000000), "million kmers processed...")
        }
        kmersProcessed++

        kmerInt, minimizer, relative = response.KmerInt, response.Minimizer, response.Relative
        minCache[kmerInt] = minimizer

        kmer, exists = kmerMap[kmerInt]
        if exists {
            lca_ID = *(*uint16)(unsafe.Pointer(&kmer[8]))
            lca = rg.ClassTree.getLCA(rg.ClassTree.NodesByID[lca_ID], relative)
            *(*uint16)(unsafe.Pointer(&kmer[8])) = lca.ID
        } else {
            kmer = new(Kmer)
            *(*uint64)(unsafe.Pointer(&kmer[0])) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = relative.ID
            kmerMap[kmerInt] = kmer
        }
    }
    fmt.Println("...all kmers processed")
    fmt.Println()

    if rg.Flags.MemProfile {
        pprof.WriteHeapProfile(rg.memProfFile)
    }

    return rg.populateKraken(minCache, kmerMap)
}

func (rg *RepeatGenome) populateKraken(minCache map[uint64]uint64, kmerMap map[uint64]*Kmer) error {
    if len(kmerMap) != len(minCache) {
        panic("lengths of kmerMap and minCache are incompatible")
    }
    numUniqKmers := uint64(len(kmerMap))

    fmt.Println(comma(numUniqKmers), "unique kmers generated")

    minMap := make(MinMap)

    for kmerInt, kmer := range kmerMap {
        minimizer := minCache[kmerInt]

        if kmers, exists := minMap[minimizer]; exists {
            minMap[minimizer] = append(kmers, kmer)
        } else {
            minMap[minimizer] = PKmers{kmer}
        }
        delete(kmerMap, kmerInt)
    }

    var minMapCount uint64 = 0
    for i := range minMap {
        minMapCount += uint64(len(minMap[i]))
    }
    if minMapCount != numUniqKmers {
        panic("error populated minMap - length not equal to that of kmerMap")
    }

    numUniqMins := uint64(len(minMap))
    fmt.Println(comma(numUniqMins), "unique minimizers used")

    // manual deletion to save memory
    kmerMap = nil
    runtime.GC()

    rg.MinCounts = make(map[uint64]uint64, numUniqMins)
    rg.SortedMins = make(Uint64Slice, 0, numUniqMins)
    for thisMin, kmers := range minMap {
        if len(kmers) == 0 {
            panic("empty kmer list stored in minMap")
        }
        rg.MinCounts[thisMin] += uint64(len(kmers))
        rg.SortedMins = append(rg.SortedMins, thisMin)
        sort.Sort(kmers)
    }
    sort.Sort(rg.SortedMins)

    if uint64(len(rg.SortedMins)) != numUniqMins {
        panic(fmt.Errorf("error populating RepeatGenome.SortedMins - %d minimizers inserted rather than expected %d", len(rg.SortedMins), numUniqMins))
    }

    if rg.Flags.WriteKraken {
        err := rg.WriteMins(minMap)
        if err != nil {
            return err
        }
    }

    var currOffset uint64 = 0
    rg.OffsetsToMin = make(map[uint64]uint64)
    rg.Kmers = make(Kmers, 0, numUniqKmers)
    for _, thisMin := range rg.SortedMins {
        rg.OffsetsToMin[thisMin] = currOffset
        currOffset += uint64(len(minMap[thisMin]))
        for _, kmer := range minMap[thisMin] {
            rg.Kmers = append(rg.Kmers, *kmer)
        }
        delete(minMap, thisMin)
    }

    if uint64(len(rg.Kmers)) != numUniqKmers {
        panic(fmt.Errorf("error populating RepeatGenome.Kmers - %d kmers inserted rather than expected %d", len(rg.Kmers), numUniqKmers))
    }

    if rg.Flags.MemProfile {
        pprof.WriteHeapProfile(rg.memProfFile)
    }

    return nil
}

func (rg *RepeatGenome) numKmers() uint64 {
    var k = int(rg.K)
    var numKmers uint64 = 0

    splitOnN := func(c rune) bool { return c == 'n' }

    for i := range rg.Matches {
        match := &rg.Matches[i]
        seq := rg.chroms[match.SeqName][match.SeqName][match.SeqStart:match.SeqEnd]
        seqs := bytes.FieldsFunc([]byte(seq), splitOnN)
        for j := range seqs {
            if len(seqs[j]) >= k {
                numKmers += uint64(len(seqs[j]) - k + 1)
            }
        }
    }
    return numKmers
}

func (rg *RepeatGenome) getKmer(kmerInt uint64) *Kmer {
    _, minimizer := getMinimizer(kmerInt, rg.K, rg.M)
    startInd, minExists := rg.OffsetsToMin[minimizer]
    if !minExists {
        return nil
    }
    var endInd uint64
    minCount, minExists := rg.MinCounts[minimizer]
    if minExists {
        endInd = startInd + minCount
    } else {
        panic("minimizer exists in RepeatGenome.OffsetsToMin but not RepeatGenome.MinCounts")
    }
    
    if endInd > uint64(len(rg.Kmers)) {
        panic(fmt.Errorf("getKmer(): out-of-bounds RepeatGenome.Kmers access (len(rg.Kmers) = %d, endInd = %d)", len(rg.Kmers), endInd))
    }

    if !sort.IsSorted(rg.Kmers[startInd:endInd]) {
        panic("minimizer's kmers not sorted")
    }

    xs := make([]uint64, 0, 0)
    i, j := startInd, endInd
    for i < j {
        x := (j+i)/2
        xs = append(xs, x)
        thisKmerInt := *(*uint64)(unsafe.Pointer(&rg.Kmers[x][0]))
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

type SeqAndClass struct {
    Seq   Seq
    Class *ClassNode
}

/*
// not that the caller is responsible for closing the channel
func (rg *RepeatGenome) kmerSeqFeed(seq TextSeq) chan uint64 {
    c := make(chan uint64)

    go func() {
        numKmers := uint8(len(seq)) - (rg.K - 1)
        var i uint8
    KmerLoop:
        for i = 0; i < numKmers; i++ {
            kmerSeq := seq[i : i+rg.K]

            // an ugly but necessary n-skipper
            for j := rg.K - 1; j >= 0; j-- {
                if kmerSeq[j] == byte('n') {
                    i += j
                    continue KmerLoop
                }
                // necessary for j to be unsigned
                if j == 0 { break }
            }

            kmerInt := seqToInt(string(kmerSeq))
            c <- minU64(kmerInt, intRevComp(kmerInt, rg.K))
        }
    }()

    return c
}
*/

type ReadResponse struct {
    Seq       []byte
    ClassNode *ClassNode
}

func (rg *RepeatGenome) ClassifyReads(reads []TextSeq, responseChan chan ReadResponse) {
    var kmerSet map[uint64]bool
    var byteBuf TextSeq
    if rg.Flags.Debug {
        byteBuf = make(TextSeq, rg.K, rg.K)
        kmerSet = make(map[uint64]bool, len(rg.Kmers))
        for _, kmer := range rg.Kmers {
            kmerSeq := *(*uint64)(unsafe.Pointer(&kmer[0]))
            kmerSet[kmerSeq] = true
        }
    }


ReadLoop:
    for _, readSeq := range reads {
        k_ := int64(rg.K)
        numKmers := int64(len(readSeq)) - k_ + 1

        var i int64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            for j := k_ + i - 1; j >= i ; j-- {
                if readSeq[j] == byte('n') {
                    i += j - i
                    continue KmerLoop
                }
            }

            kmerBytes := readSeq[i : i+k_]
            kmerInt := seqToInt((kmerBytes))
            kmerInt = minU64(kmerInt, intRevComp(kmerInt, rg.K))
            kmer := rg.getKmer(kmerInt)

            if rg.Flags.Debug && kmer == nil && kmerSet[kmerInt] {
                fillKmerBuf(byteBuf, kmerInt)
                panic(fmt.Errorf("RepeatGenome.getKmer() returned nil for %s, but kmer exists\n", byteBuf))
            }

            if kmer != nil {
                fillKmerBuf(byteBuf, kmerInt)
                lcaID := *(*uint16)(unsafe.Pointer(&kmer[8]))
                responseChan <- ReadResponse{readSeq, rg.ClassTree.NodesByID[lcaID]}
                // only use the first matched kmer
                continue ReadLoop
            }
        }
        responseChan <- ReadResponse{readSeq, nil}
    }
    close(responseChan)
}

func (rg *RepeatGenome) GetReadClassChan(reads []TextSeq) chan ReadResponse {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    numCPU := uint64(runtime.NumCPU())
    if rg.Flags.Debug {
        fmt.Printf("GetReadClassChan() using %d CPUs\n", numCPU)
    }
    runtime.GOMAXPROCS(int(numCPU))

    responseChans := make([]chan ReadResponse, 0, numCPU)

    numReads := uint64(len(reads))
    var i uint64
    for i = 0; i < numCPU; i++ {
        responseChans = append(responseChans, make(chan ReadResponse, 50))
        startInd := i * (numReads / numCPU)
        endInd := ((i + 1) * numReads) / numCPU
        go rg.ClassifyReads(reads[startInd : endInd], responseChans[i])
    }

    var wg sync.WaitGroup
    //fmt.Println("GetReadClassChan() using %d chans for %d reads\n", len(responseChans), len(reads))
    //os.Exit(0)
    wg.Add(len(responseChans))
    master := make(chan ReadResponse)

    for _, respChan := range responseChans {
        go func(respChan chan ReadResponse) {
            for resp := range respChan {
                master <- resp
            }
            wg.Done()
        }(respChan)
    }

    go func() {
        wg.Wait()
        close(master)
    }()

    return master
}

func (rg *RepeatGenome) ProcessReads() (error, chan ReadResponse) {
    workingDirName, err := os.Getwd()
    if err != nil {
        return err, nil
    }
    readsDirName := workingDirName + "/" + rg.Name + "-reads"
    currDir, err := os.Open(readsDirName)
    if err != nil {
        return err, nil
    }
    fileinfos, err := currDir.Readdir(-1)
    if err != nil {
        return err, nil
    }
    processedFiles := []os.FileInfo{}
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 5 && fileinfo.Name()[len(fileinfo.Name())-5 : ] == ".proc" {
            processedFiles = append(processedFiles, fileinfo)
        }
    }
    var reads []TextSeq
    for _, fileinfo := range processedFiles {
        _, theseReadsBytes := fileLines(readsDirName + "/" + fileinfo.Name())
        for _, bytesLine := range theseReadsBytes {
            reads = append(reads, TextSeq(bytesLine))
        }
    }

    return nil, rg.GetReadClassChan(reads)
}
