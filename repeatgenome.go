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

   Functions for pointer accesses?

   We should probably just store minimizers fully - they're only 32 bits

   KmerInt.Minimize() logic could be changed now that minimizers are 32 bits

   Should a Seq's first field be a *byte to discard the extra two fields?

   Should probably make a file solely for type defs.

   Reads are currently kept in TextSeq form until the bitter end because, with Go's referenced based slices, there's no compelling reason not to, and because they're easier (and probably faster) to manipulate than Seqs. This may change at some point, though.

   If a minimizer is associated with a single repeat type, can we use that heuristically?

   Error handling should be updated with a custom ParseError type - panics should be removed, excepting performance-cricial sequence manipulation functions

   Should consider splitting at hyphenated class names like TcMar-Tc1

   The concurrent read-kmer generator could be reintroduced using a select statement.

   Should probably restrict activity of chans with directionals

   It would make sense to discard kmers associated with ClassNodes greater than a certain size.

   Kmer counting should be re-added eventually - it's currently excluded for performance reasons because we aren't using it.

   We should test a version that doesn't cache minimizers, as that seems to be a needless bottleneck. It could also be conditional on the number of CPUs available.

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

// These variables are used ubiquitously, especially in performance-critical functions, so we grudgingly make them globals
var k, m uint8

type Flags struct {
    Debug       bool
    CPUProfile  bool
    MemProfile  bool
    Minimize    bool
    WriteKraken bool
    WriteJSON   bool
}

// Match.SW_Score - Smith-Waterman score, describing the likeness to the repeat reference sequence
// Match.PercDiv - "% substitutions in matching region compared to the consensus" - RepeatMasker docs
// Match.PercDel - "% of bases opposite a gap in the query sequence (deleted bp)" - RepeatMasker docs
// Match.PercIns - "% of bases opposite a gap in the repeat consensus (inserted bp)" - RepeatMasker docs
// Match.SeqName -  the reference genome file this match came from (typically the chromosome)
// Match.SeqStart -  the starting index (inclusive) in the reference genome
// Match.SeqEnd -  the ending index (exclusive) in the reference genome
// Match.SeqRemains - the number of bases past the end of the match in the relevant reference seqence
// Match.IsRevComp -  the match may be for the complement of the reference sequence
// Match.RepeatClass - the repeat's full ancestry, including its repeat class and repeat name (which are listed separately in the RepeatMasker output file)
// Match.RepeatStart-  the starting index in the repeat consensus sequence
// Match.RepeatEnd -  the ending sequence (exclusive) in the consensus repeat  sequence
// Match.RepeatRemains - the number of bases past the end of the match in the consensus repeat sequence
// Match.InsertionID - a numerical ID for the repeat type (starts at 1)
//
//     below are not in parsed data file
// Match.RepeatName - simply repeatClass concatenated - used for printing and map indexing
// Match.ClassNode - pointer to corresponding ClassNode in RepeatGenome.ClassTree
// Match.Repeat - pointer to corresponding Repeat struct in RepeatGenome.Repeats
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
    ID         uint64
}

type RepeatGenome struct {
    Name       string
    Flags Flags
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    chroms         map[string](map[string]TextSeq)
    Kmers          Kmers
    FullKmers      FullKmers
    // stores the offset of each minimizer's first kmer in RepeatGenome.Kmers - indexed by the minimizer's index in SortedMins
    OffsetsToMin   []uint64
    // stores the number of kmers that each minimizer is associated with
    SortedMins MinInts
    Matches        Matches
    ClassTree      ClassTree
    Repeats        Repeats
    RepeatMap      map[string]*Repeat
    memProfFile    *os.File
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

type MuxKmers struct {
    sync.Mutex
    Kmers Kmers
}

// Used to differentiate sequence representations with one base per byte (type TextSeq) from those with four bases per byte (type Seq).
type TextSeq []byte

// Used to clarify context of integers, and to differentiate full Kmers (which include an LCA ID) from integer-represented kmer sequences.
type KmerInt uint64
type KmerInts []KmerInt

type MinInt uint32
type MinInts []MinInt

/*
// Indexes the base of the associated kmer that is the starting index of its minimizer
// If < 32, the minimizer is the positive strand representation
// Otherwise, the minimizer is the reverse complement of kmer[minkey-32 : minkey+m-32]
type MinKey uint8

// Stores a KmerInt and a MinKey in that order
type MinPair [9]byte
type MinPairs []MinPair
*/

type MinPair [12]byte
type MinPairs []MinPair

// can store a kmer where k <= 32
// the value of k is not stored in the struct, but rather in the RepeatGenome, for memory efficiency
// first eight bits are the int representation of the sequence
// the last two are the LCA ID
type Kmer [10]byte

// Like a Kmer, but with a MinInt between the KmerInt and the LCA ID.
type FullKmer [14]byte
type FullKmers []FullKmer

// as with the Kmer type, each base is represented by two bits
// any excess bits are the first bits of the first byte (seq is right-justified)
// remember that len(Seq.Bytes) is not the actual number of bases, but rather the number of bytes necessary to represent them
type Seq struct {
    Bytes []byte
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
type MinMap map[MinInt]Kmers
type Repeats []Repeat
type Matches []Match

//type Chroms map[string](map[string]TextSeq)

/*
type MinPair struct {
    KmerInt   Kmer
    Minimizer MinInt
}
*/

type ThreadResponse struct {
    KmerInt   KmerInt
    MinInt    MinInt
    Relative *ClassNode
}

type MinCache map[KmerInt]MinInt

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

        match.ID = uint64(len(matches))

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

func Generate(genomeName string, k_arg, m_arg uint8, rgFlags Flags) (error, *RepeatGenome) {
    var err error
    if m > k {
        return fmt.Errorf("m must be <= k"), nil
    }
    k = k_arg
    m = m_arg
    fmt.Println("k =", k)
    fmt.Println("m =", m)
    fmt.Println()
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

    /*
    if rg.Flags.Minimize {
        // calling the parallel minimizer and writing the result
        rg.getKrakenSlice()
    }
    */
    if rg.Flags.Minimize {
        rg.getMinPairSlices()
        os.Exit(0)
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
    thisMin := testSeq.kmerInt().Minimize()
    fmt.Println("getMinimizer('tgctcctgtcatgcatacgcaggtcatgcat'): ")
    thisMin.print()
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

/*
// returns ancestry list, beginning at (and including) self, excluding root
func (cn *ClassNode) getAncestry() []*ClassNode {
    ancestry := []*ClassNode{}
    walker := cn
    for !walker.IsRoot {
        ancestry = append(ancestry, walker)
        walker = walker.Parent
    }
    return ancestry
}
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
// some of the logic in here is deeply nested or non-obvious for efficiency's sake
// specifically, we made sure not to make any heap allocations, which means reverse complements can never be explicitly evaluated
func (rg *RepeatGenome) minimizeThread(matchStart, matchEnd uint64, c chan ThreadResponse) {
    k_ := uint64(k)

    for i := matchStart; i < matchEnd; i++ {
        match := &rg.Matches[i]
        seq := rg.chroms[match.SeqName][match.SeqName]
        matchLen := match.SeqEnd - match.SeqStart
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

            kmerInt := TextSeq(seq[j : j+k_]).kmerInt()
            // make the sequence strand-agnostic
            kmerInt = minKmerInt(kmerInt, kmerInt.revComp())

            thisMin := kmerInt.Minimize()

            c <- ThreadResponse{kmerInt, thisMin, match.ClassNode}
        }
    }

    close(c)
}

type UpdateInfo struct {
    muxKmers *MuxKmers
    kmerInt KmerInt
    relative *ClassNode
}

func (rg *RepeatGenome) krakenUpdateThread(wg *sync.WaitGroup, updateChan chan UpdateInfo) {
    for updateInfo := range updateChan {
        muxKmers, kmerInt, relative := updateInfo.muxKmers, updateInfo.kmerInt, updateInfo.relative
        muxKmers.Lock()
        kmerExists := false

        for i, kmer := range muxKmers.Kmers {
            // the case that we've already processed this exact kmer - we just update the LCA
            if kmerInt == *(*KmerInt)(unsafe.Pointer(&kmer[0])) {
                kmerExists = true
                prev_LCA_ID := *(*uint16)(unsafe.Pointer(&kmer[8]))
                lca := rg.ClassTree.getLCA(rg.ClassTree.NodesByID[prev_LCA_ID], relative)

                *(*uint16)(unsafe.Pointer(&kmer[8])) = lca.ID
                muxKmers.Kmers[i] = kmer

                break
            }
        }

        if !kmerExists {
            var kmer Kmer
            *(*KmerInt)(unsafe.Pointer(&kmer[0])) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = relative.ID
            muxKmers.Kmers = append(muxKmers.Kmers, kmer)
        }
        muxKmers.Unlock()
        wg.Done()
    }
}

func (rg *RepeatGenome) getKrakenSlice() error {
    // a rudimentary way of deciding how many threads to allow, should eventually be improved
    numCPU := runtime.NumCPU()
    if rg.Flags.Debug {
        fmt.Printf("getKrakenSlice() using %d CPUs\n", ceilDiv(numCPU, 2))
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

    minToKmers := make(map[MinInt]*MuxKmers)
    numUpdateThreads := numCPU / 2
    updateChan := make(chan UpdateInfo, 1000)
    var wg sync.WaitGroup
    for i := 0; i < numUpdateThreads; i++ {
        go rg.krakenUpdateThread(&wg, updateChan)
    }
    
    var kmersProcessed uint64 = 0

    // below is the atomic section of minimizing
    // this seems to be the rate-limiting section, as 24+ goroutines use only ~9-10 CPU-equivalents
    // it should therefore be optimized first
    for resp := range mergeThreadResp(threadChans) {
        if kmersProcessed%5000000 == 0 {
            fmt.Println(comma(kmersProcessed/1000000), "million kmers processed...")
        }
        kmersProcessed++

        kmerInt, minInt, relative := resp.KmerInt, resp.MinInt, resp.Relative
        
        if muxKmers, minExists := minToKmers[minInt]; minExists {
            wg.Add(1)
            updateChan <- UpdateInfo{muxKmers, kmerInt, relative}
        // ...otherwise we initialize it in the kmerMap
        } else {
            var kmer Kmer
            *(*KmerInt)(unsafe.Pointer(&kmer[0])) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = relative.ID
            // we don't need to lock because the update is atomic with the addition
            minToKmers[minInt] = &MuxKmers{ sync.Mutex{}, Kmers{kmer} }
        }
    }

    wg.Wait()
    close(updateChan)

    fmt.Println("...all kmers processed")
    fmt.Println()

    if rg.Flags.MemProfile {
        pprof.WriteHeapProfile(rg.memProfFile)
    }

    return rg.populateKraken(minToKmers)
}

func (rg *RepeatGenome) populateKraken(minToKmers map[MinInt]*MuxKmers) error {
    var numUniqKmers uint64 = 0

    for minInt, muxKmers := range minToKmers {
        kmers := muxKmers.Kmers
        numUniqKmers += uint64(len(kmers))
        sort.Sort(kmers)
        rg.SortedMins = append(rg.SortedMins, minInt)
    }
    sort.Sort(rg.SortedMins)

    fmt.Println(comma(numUniqKmers), "unique kmers generated")

    numUniqMins := uint64(len(rg.SortedMins))
    fmt.Println(comma(numUniqMins), "unique minimizers used")

    var currOffset uint64 = 0
    rg.OffsetsToMin = make([]uint64, 0, numUniqMins)
    rg.Kmers = make(Kmers, 0, numUniqKmers)
    for _, thisMin := range rg.SortedMins {
        rg.OffsetsToMin = append(rg.OffsetsToMin, currOffset)
        currOffset += uint64(len(minToKmers[thisMin].Kmers))
        for _, kmer := range minToKmers[thisMin].Kmers {
            rg.Kmers = append(rg.Kmers, kmer)
        }
        // delete(minMap, thisMin)     // useless because it's going to be almost immediately nilled anyway
    }

    if uint64(len(rg.Kmers)) != numUniqKmers {
        panic(fmt.Errorf("error populating RepeatGenome.Kmers - %d kmers inserted rather than expected %d", len(rg.Kmers), numUniqKmers))
    }

    if rg.Flags.WriteKraken {
        err := rg.WriteMins()
        if err != nil {
            return err
        }
    }

    if rg.Flags.MemProfile {
        pprof.WriteHeapProfile(rg.memProfFile)
    }

    return nil
}
*/

func (rg *RepeatGenome) getMinIndex(minInt MinInt) (bool, uint64) {
    var i uint64 = 0
    j := uint64(len(rg.SortedMins))

    for i < j {
        x := (i+j)/2

        if minInt == rg.SortedMins[x] {
            return true, x
        } else if minInt < rg.SortedMins[x] {
            j = x
        } else {
            i = x + 1
        }
    }

    return false, 0
}

func (rg *RepeatGenome) getKmer(kmerInt KmerInt) *Kmer {
    minimizer := kmerInt.Minimize()
    minExists, minIndex := rg.getMinIndex(minimizer)
    if !minExists {
        return nil
    }
    startInd := rg.OffsetsToMin[minIndex]
    var endInd uint64
    if minIndex == uint64(len(rg.SortedMins)) - 1 {
        endInd = uint64(len(rg.Kmers))
    } else {
        endInd = rg.OffsetsToMin[minIndex+1]
    }
    
    if endInd > uint64(len(rg.Kmers)) {
        panic(fmt.Errorf("getKmer(): out-of-bounds RepeatGenome.Kmers access (len(rg.Kmers) = %d, endInd = %d)", len(rg.Kmers), endInd))
    }

    /*
    if !sort.IsSorted(rg.Kmers[startInd:endInd]) {
        panic("minimizer's kmers not sorted")
    }
    */

    // simple binary search within the range of RepeatGenome.Kmers that has this kmer's minimizer
    i, j := startInd, endInd
    for i < j {
        x := (j+i)/2
        thisKmerInt := *(*KmerInt)(unsafe.Pointer(&rg.Kmers[x][0]))
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
func (rg *RepeatGenome) kmerSeqFeed(seq TextSeq) chan uint64 {
    c := make(chan uint64)

    go func() {
        numKmers := uint8(len(seq)) - (k - 1)
        var i uint8
    KmerLoop:
        for i = 0; i < numKmers; i++ {
            kmerSeq := seq[i : i+k]

            // an ugly but necessary n-skipper
            for j := k - 1; j >= 0; j-- {
                if kmerSeq[j] == byte('n') {
                    i += j
                    continue KmerLoop
                }
                // necessary for j to be unsigned
                if j == 0 { break }
            }

            kmerInt := seqToInt(string(kmerSeq))
            c <- minU64(kmerInt, intRevComp(kmerInt, k))
        }
    }()

    return c
}
*/

type ReadResponse struct {
    Seq       TextSeq
    ClassNode *ClassNode
}

// This function assumes that the Seqs in readSeqs do not contain 'n's.
// The output reads of sequencing simulators will generally contain 'n's if the input reference genome does.
// They must therefore be filtered upstream.
func (rg *RepeatGenome) ClassifyReads(readTextSeqs []TextSeq, responseChan chan ReadResponse) {
    var kmerSet map[KmerInt]bool
    var byteBuf TextSeq
    if rg.Flags.Debug {
        byteBuf = make(TextSeq, k, k)
        kmerSet = make(map[KmerInt]bool, len(rg.Kmers))
        for _, kmer := range rg.Kmers {
            kmerSeq := *(*KmerInt)(unsafe.Pointer(&kmer[0]))
            kmerSet[kmerSeq] = true
        }
    }


ReadLoop:
    for _, read := range readTextSeqs {
        // we use sign int64s in this triple-nested loop because the innermost one counts down and would otherwise overflow
        // this isn't a significant limitation because reads are never big enough to overflow one
        k_ := int64(k)
        numKmers := int64(len(read)) - k_ + 1

        var i int64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            // if this is ever revived, all the loop ints must be made signed again
            for j := k_ + i - 1; j >= i ; j-- {
                if read[j] == byte('n') {
                    i += j - i
                    continue KmerLoop
                }
            }

            kmerBytes := read[i : i+k_]
            kmerInt := kmerBytes.kmerInt()
            kmerInt = minKmerInt(kmerInt, kmerInt.revComp())
            kmer := rg.getKmer(kmerInt)

            if rg.Flags.Debug && kmer == nil && kmerSet[kmerInt] {
                fillKmerBuf(byteBuf, kmerInt)
                panic(fmt.Errorf("RepeatGenome.getKmer() returned nil for %s, but kmer exists\n", byteBuf))
            }

            if kmer != nil {
                fillKmerBuf(byteBuf, kmerInt)
                lcaID := *(*uint16)(unsafe.Pointer(&kmer[8]))
                responseChan <- ReadResponse{read, rg.ClassTree.NodesByID[lcaID]}
                // only use the first matched kmer
                continue ReadLoop
            }
        }
        responseChan <- ReadResponse{read, nil}
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

type MinPairResponse struct {
    ClassNodeID uint16
    MinPairs MinPairs
}

// WARNING: Will return nil rather than empty slice in response.MinPairs
// Could probably use a map to prevent duplicates in the future, although this wouldn't make a difference if the slice size were still manually chosen
func (rg *RepeatGenome) getMinPairs(match *Match) MinPairResponse {
    k_ := uint64(k)
    start, end := match.SeqStart, match.SeqEnd
    matchSeq := rg.chroms[match.SeqName][match.SeqName][start : end]
    if len(matchSeq) < int(k) {
        return MinPairResponse{match.ClassNode.ID, nil}
    }
    numKmers := end - start - k_ + 1
    var minPairs MinPairs

    var i uint64
KmerLoop:
    for i = 0; i < numKmers; i++ {

        // loop logic uses manual condition to prevent unsigned int overflow
        for j := k_ - 1; ; j-- {
            if matchSeq[i+j] == 'n' {
                i += uint64(j)
                continue KmerLoop
            }
            if j == 0 { break }
        }

        kmerInt := TextSeq(matchSeq[i : i+k_]).kmerInt()
        kmerInt = minKmerInt(kmerInt, kmerInt.revComp())
        var minPair MinPair
        *(*KmerInt)(unsafe.Pointer(&minPair)) = kmerInt
        *(*MinInt)(unsafe.Pointer(&minPair[8])) = kmerInt.Minimize()
        minPairs = append(minPairs, minPair)
    }

    return MinPairResponse{match.ClassNode.ID, minPairs}
}

func (rg *RepeatGenome) getMinPairSlices() {
    numCPU := runtime.NumCPU()
    runtime.GOMAXPROCS(numCPU)

    matchChan := make(chan *Match)
    go func() {
        for i := range rg.Matches {
            matchChan <- &rg.Matches[i]
        }
        close(matchChan)
    }()

    var wg sync.WaitGroup
    wg.Add(numCPU)
    respChan := make(chan MinPairResponse)
    for i := 0; i < numCPU; i++ {
        go func() {
            for match := range matchChan {
                minPairs := rg.getMinPairs(match)
                respChan <- minPairs
            }
            wg.Done()
        }()
    }

    var minPairResps []MinPairResponse
    go func() {
        for resp := range respChan {
            if resp.MinPairs != nil {
                minPairResps = append(minPairResps, resp)
            }
        }
    }()

    wg.Wait()
    close(respChan)

    numMinPairs := 0
    for _, resp := range minPairResps {
        numMinPairs += len(resp.MinPairs)
    }
    fmt.Println("total number of (non-unique) kmers processed:", comma(uint64(numMinPairs)))

    fullKmers := make(FullKmers, 0, numMinPairs)
    for _, resp := range minPairResps {
        for _, minPair := range resp.MinPairs {
            var fullKmer FullKmer
            *(*[12]byte)(unsafe.Pointer(&fullKmer)) = minPair
            *(*uint16)(unsafe.Pointer(&fullKmer[12])) = resp.ClassNodeID
            fullKmers = append(fullKmers, fullKmer)
        }
        resp.MinPairs = nil    // aid garbage collection
    }

    /*
    minPairResps = nil
    runtime.GC()
    */

    sort.Sort(fullKmers)
    fmt.Println("raw kmer list sorted")
    /*
    // make rg.MinPairs unique
    if len(fullKmers) > 0 {
        rg.FullKmers = append(rg.FullKmers, fullKmers[0])
    }
    for i := 1; i < len(fullKmers); i++ {
        if *(*uint64)(unsafe.Pointer(&fullKmers[i])) != *(*uint64)(unsafe.Pointer(&fullKmers[i-1])) {
            rg.FullKmers = append(rg.FullKmers, fullKmers[i])
        } else {
            lastLCA := rg.ClassTree.NodesByID[*(*uint16)(unsafe.Pointer(&rg.FullKmers[len(rg.FullKmers)-1][9]))]
            newNode := rg.ClassTree.NodesByID[*(*uint16)(unsafe.Pointer(&fullKmers[i][9]))]
            *(*uint16)(unsafe.Pointer(&rg.FullKmers[len(rg.FullKmers)-1])) = rg.ClassTree.getLCA(lastLCA, newNode).ID
        }
    }
    */
    rg.parallelFullKmers(fullKmers)

    /*
    fullKmers = nil
    runtime.GC()
    */

    fmt.Println(comma(uint64(len(rg.FullKmers))), "unique kmers processed")
}

type ReducePair struct {
    Loc *FullKmer
    Set FullKmers
}

func (rg *RepeatGenome) parallelFullKmers(fullKmers FullKmers) {
    numCPU := runtime.NumCPU()
    runtime.GOMAXPROCS(numCPU)

    pairChan := make(chan ReducePair, 500)
    var wg sync.WaitGroup
    wg.Add(numCPU)

    for i := 0; i < numCPU; i++ {
        go func() {
            for pair := range pairChan {
                lcaID := *(*uint16)(unsafe.Pointer(&pair.Set[0][12]))
                currLCA := rg.ClassTree.NodesByID[lcaID]

                for i := 1; i < len(pair.Set); i++ {
                    lcaID = *(*uint16)(unsafe.Pointer(&pair.Set[i][12]))
                    currLCA = rg.ClassTree.getLCA(currLCA, rg.ClassTree.NodesByID[lcaID])
                }

                *(*uint16)(unsafe.Pointer(&pair.Loc[12])) = currLCA.ID
            }
            wg.Done()
        }()
    }

    go func() {
        numKmers := uint64(len(fullKmers))
        var start, end uint64 = 0, 0
        for end < numKmers {
            start = end
            end++
            for end < numKmers && *(*uint64)(unsafe.Pointer(&fullKmers[end-1])) == *(*uint64)(unsafe.Pointer(&fullKmers[end])) {
                end++
            }
            rg.FullKmers = append(rg.FullKmers, FullKmer{})
            pairChan <- ReducePair{&rg.FullKmers[len(rg.FullKmers)-1], fullKmers[start : end]}
        }

        close(pairChan)
    }()

    wg.Wait()
}
