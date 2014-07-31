package repeatgenome

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
// masks contain k and m consecutive right aligned 1 bits respectively (e.g. "0000011111111")
var kMask, mMask KmerInt
var debug bool

// a rudimentary way of deciding how many threads to allow, should eventually be improved
var numCPU = runtime.NumCPU()

type Config struct {
    Name        string
    K           uint8
    M           uint8
    Debug       bool
    CPUProfile  bool
    MemProfile  bool
    WriteLib    bool
    ForceGen    bool
}

/*
   Match.SW_Score - Smith-Waterman score, describing the likeness to the repeat reference sequence
   Match.PercDiv - "% substitutions in matching region compared to the consensus" - RepeatMasker docs
   Match.PercDel - "% of bases opposite a gap in the query sequence (deleted bp)" - RepeatMasker docs
   Match.PercIns - "% of bases opposite a gap in the repeat consensus (inserted bp)" - RepeatMasker docs
   Match.SeqName -  the reference genome file this match came from (typically the chromosome)
   Match.SeqStart -  the starting index (inclusive) in the reference genome
   Match.SeqEnd -  the ending index (exclusive) in the reference genome
   Match.SeqRemains - the number of bases past the end of the match in the relevant reference seqence
   Match.IsRevComp -  the match may be for the complement of the reference sequence
   Match.RepeatClass - the repeat's full ancestry, including its repeat class and repeat name (which are listed separately in the RepeatMasker output file)
   Match.RepeatStart-  the starting index in the repeat consensus sequence
   Match.RepeatEnd -  the ending sequence (exclusive) in the consensus repeat  sequence
   Match.RepeatRemains - the number of bases past the end of the match in the consensus repeat sequence
   Match.InsertionID - a numerical ID for the repeat type (starts at 1)
  
       below are not in parsed data file
   Match.RepeatName - simply repeatClass concatenated - used for printing and map indexing
   Match.ClassNode - pointer to corresponding ClassNode in RepeatGenome.ClassTree
   Match.Repeat - pointer to corresponding Repeat struct in RepeatGenome.Repeats
*/
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
    // maps a chromosome name to a map of its sequences
    // as discussed above, though, matches only contain 1D sequence indexes
    chroms         Chroms
    Kmers          Kmers
    // stores the offset of each minimizer's first kmer in RepeatGenome.Kmers - indexed by the minimizer's index in SortedMins
    MinOffsets     map[MinInt]uint64
    MinCounts      map[MinInt]uint32
    SortedMins     MinInts
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

// Indexes the base of the associated kmer that is the starting index of its minimizer
// If < 32, the minimizer is the positive strand representation
// Otherwise, the minimizer is the reverse complement of kmer[minkey-32 : minkey+m-32]
type MinKey uint8

type MinPair [12]byte
type MinPairs []MinPair

// can store a kmer where k <= 32
// the value of k is not stored in the struct, but rather in the RepeatGenome, for memory efficiency
// first eight bits are the int representation of the sequence
// the last two are the LCA ID
type Kmer [10]byte

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

type Chroms map[string](map[string]TextSeq)

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

func New(config Config) (error, *RepeatGenome) {
    debug = config.Debug

    if config.K < config.M {
        return fmt.Errorf("m must be <= k"), nil
    }
    k = config.K
    m = config.M
    kMask = (1 << (2*k)) - 1
    mMask = (1 << (2*m)) - 1
    fmt.Println("k =", k)
    fmt.Println("m =", m)
    fmt.Println()

    runtime.GOMAXPROCS(numCPU)

    // we popoulate the RepeatGenome mostly with helper functions
    // we should consider whether it makes more sense for them to alter the object directly, than to return their results
    rg := new(RepeatGenome)
    rg.Name = config.Name

    var err error

    if config.MemProfile {
        os.Mkdir("profiles", os.ModeDir)
        rg.memProfFile, err = os.Create("profiles/" + rg.Name + ".memprof")
        if err != nil {
            return IOError{"RepeatGenome.getKrakenSlice()", err}, nil
        }
        pprof.WriteHeapProfile(rg.memProfFile)
        defer rg.memProfFile.Close()
    }

    err, rg.chroms = parseGenome(rg.Name)
    if err != nil {
        return err, nil
    }
    err, rg.Matches = parseMatches(rg.Name)
    if err != nil {
        return err, nil
    }
    rg.getRepeats()
    rg.getClassTree()

    if config.ForceGen {
        rg.genKrakenLib()
        err = rg.WriteKraken()
        if err != nil {
            return err, nil
        }
    } else {
        minsFile, err := os.OpenFile(rg.Name + ".mins", os.O_RDONLY, 0400)
        // implies the file already exists
        if err == nil {
            fmt.Println("Kraken library file exists - using contents\n")
            err = rg.ReadKraken(minsFile)
            if err != nil {
                return IOError{"Kmers.ReadKraken()", err}, nil
            }
        } else if os.IsNotExist(err) {    // the case that there isn't a written file yet
            fmt.Println("Kraken library file doesn't exist - generating library\n")
            rg.genKrakenLib()
            err = rg.WriteKraken()
            if err != nil {
                return err, nil
            }
        } else {    // otherwise we're dealing with a generic error of some sort
            return err, nil
        }
    }


    if debug {
        rg.RunDebugTests()
    }

    return nil, rg
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

type ReadResponse struct {
    Seq       TextSeq
    ClassNode *ClassNode
}

// This function is a mess and needs cleaning up
func (rg *RepeatGenome) ClassifyReads(readTextSeqs []TextSeq, responseChan chan ReadResponse) {
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
            kmerInt := kmerBytes.kmerInt().canonicalRepr()
            kmer := rg.getKmer(kmerInt)

            if kmer != nil {
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
    responseChans := make([]chan ReadResponse, 0, numCPU)

    numReads := uint64(len(reads))
    numCPU_ := uint64(numCPU)
    var i uint64
    for i = 0; i < uint64(numCPU_); i++ {
        responseChans = append(responseChans, make(chan ReadResponse, 50))
        startInd := i * (numReads / numCPU_)
        endInd := ((i + 1) * numReads) / numCPU_
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
