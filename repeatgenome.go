package repeatgenome

import (
    "bytes"
    "fmt"
    "io/ioutil"
    "log"
    "os"
    "runtime"
    "runtime/pprof"
    "strconv"
    "strings"
)

// These variables are used ubiquitously, especially in performance-critical functions, so we grudgingly make them globals.
// K is the kmer size (in basepairs).
// M is the minimizer size, which must be <= k.
var k, m uint8

// kMask and mMask contain k and m consecutive right aligned 1 bits respectively (e.g. "0000011111111").
var kMask, mMask KmerInt

// minSliceSize is the size of a slice that contains indexes for all possible mMers.
// It is calculated as: 1 << m
var minSliceSize uint64
var debug bool

// This is used as a rudimentary way of determining how many goroutines to spawn in concurrent sections.
var numCPU = runtime.NumCPU()

// Used to write heap memory profile.
var memProfFile *os.File

// A value of type Config is passed to the New() function, which constructs and returns a new RepeatGenome.
type Config struct {
    Name       string
    K          uint8
    M          uint8
    Debug      bool
    CPUProfile bool
    MemProfile bool
    WriteLib   bool
    ForceGen   bool
    WriteStats bool
}

/*
   RepeatMasker is a program that takes as input a set of reference repeat sequences and a reference genome.
   It outputs "matches", specific instances of supplied reference repeats in the supplied reference genome.
   These are stored in the file <genome-name>.fa.out, and are parsed line-by-line into values of this type.

   Match.SW_Score - Smith-Waterman score, describing the likeness of this match to the repeat reference sequence.
   Match.PercDiv - "% substitutions in matching region compared to the consensus" - RepeatMasker docs
   Match.PercDel - "% of bases opposite a gap in the query sequence (deleted bp)" - RepeatMasker docs
   Match.PercIns - "% of bases opposite a gap in the repeat consensus (inserted bp)" - RepeatMasker docs
   Match.SeqName - The name (without ".fa") of the reference genome FASTA file this match came from.
                   It is typically the chromosome name, such as "chr2L".
                   This is inherently unsound, as RepeatMasker gives only a 1-dimensional qualification, but FASTA-formatted reference genomes are 2-dimensional, using both filename and sequence name.
                   Reference genome FASTA files generally contain only a single sequence, with the same name as the file.
                   If this is not the case when parsing a FASTA reference genome, we print an ominous warning to stdout and use the sequence name.
   Match.SeqStart - The match's start index (inclusive and zero-indexed) in the reference genome. Note that RepeatMasker's output is one-indexed.
   Match.SeqEnd - The end index (exclusive and zero-indexed) in the reference genome.
   Match.SeqRemains - The number of bases past the end of the match in the relevant reference sequence.
   Match.IsRevComp - Whether the match was for the reverse complement of the reference repeat sequence. In this case, we manually adjust some location fields, as RepeatMasker's output gives indexes for the reverse complement of the reference sequence. This allows us to treat all matches with the same logic.
   Match.RepeatClass - The repeat's full ancestry in a slice of strings.
                       This includes its repeat class and repeat name, which are listed separately in the RepeatMasker output file.
                       Root is implicit and excluded.
   Match.RepeatStart - The start index (inclusive and zero-indexed) of this match in the repeat consensus sequence.
                       A signed integer is used because it can be negative in weird cases.
   Match.RepeatEnd - The end sequence (exclusive and zero-indexed) of this match in the consensus repeat sequence.
                     A signed integer is used, in agreement with Match.RepeatStart.
   Match.RepeatRemains - The number of bases at the end of the consensus repeat sequence that this match excludes.
   Match.InsertionID - A numerical ID that is the same only for matches of the same long terminal repeat (LTR) instance.
                       The sequences classified as <LTR name>_I or <LTR name>_int are the internal sequences of LTRs. These are less well defined than the core LTR sequence.
                       These IDs begin at 1.

   The below fields are not parsed, but rather calculated.
   Match.RepeatName - Simply Match.RepeatClass's items concatenated.
                      It is used for quick printing, and is not necessarily going to remain in the long-term.
   Match.ClassNode - A pointer to the match's corresponding ClassNode in RepeatGenome.ClassTree.
   Match.Repeat - A pointer to the match's corresponding Repeat in RepeatGenome.Repeats.
   Match.ID - A unique ID, used as a quick way of referencing and indexing a repeat.
              Pointers can generally be used, so this isn't necessarily going to remain in the long-term.
*/
type Match struct {
    SW_Score      int32
    PercDiv       float64
    PercDel       float64
    PercIns       float64
    SeqName       string
    SeqStart      uint64
    SeqEnd        uint64
    SeqRemains    uint64
    IsRevComp     bool
    RepeatClass   []string
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

type Matches []Match

/*
   RepeatGenome.Name - The name of the reference genome, such as "dm3" or "hg38".
                       This is used to name created directories, and to find directories and files that may be read from, such as a stored Kraken library and reference sequences.
   RepeatGenome.chroms - A 2-dimensional map mapping a chromosome name to a map of its sequence names to their sequences (in text form).
                         Actual 2-dimensional mapping is currently impossible because of RepeatMasker's 1-dimensional output.
   RepeatGenome.Kmers - A slice of all Kmers, sorted primarily by minimizer and secondarily by lexicographical value.
   RepeatGenome.MinOffsets - Maps a minimizer to its offset in the Kmers slice, or -1 if no kmers of this minimizer were stored.
   RepeatGenome.MinCounts - Maps a minimizer to the number of stored kmers associated with it.
   RepeatGenome.SortedMins - A sorted slice of all minimizers of stored kmers.
   RepeatGenome.Matches - All matches, indexed by their assigned IDs.
   RepeatGenome.ClassTree - Contains all information used for LCA determination and read classification.
                            It may eventually be collapsed into RepeatGenome, as accessing it is rather verbose.
   RepeatGenome.Repeats - A slice of all repeats, indexed by their assigned IDs.
   RepeatGenome.RepeatMap - Maps a fully qualified repeat name, excluding root, to its struct.
*/
type RepeatGenome struct {
    Name       string
    chroms     Chroms
    Kmers      Kmers
    MinOffsets []int64
    MinCounts  []uint32
    SortedMins MinInts
    Matches    Matches
    ClassTree  ClassTree
    Repeats    Repeats
    RepeatMap  map[string]*Repeat
}

/*
   Repeat.ID - A unique ID that we assign (not included in RepeatMasker output).
               Because these are assigned in the order in which they are encountered in <genome name>.fa.out, they are not compatible across even different versions of the same reference genome. This may change.
   Repeat.Name - The repeat's fully qualified name, excluding root.
   Repeat.ClassList - A slice of this Repeat's class ancestry from the top of the tree down, excluding root.
   Repeat.ClassNode - A pointer to the ClassNode which corresponds to this repeat.
   Repeat.Instances - A slice of pointers to all matches that are instances of this repeat.
*/
type Repeat struct {
    ID        uint64
    Name      string
    ClassList []string
    ClassNode *ClassNode
    Instances []*Match
}

type Repeats []Repeat

/*
   ClassNode.Name - This ClassNode's fully qualified name, excluding root.
   ClassNode.ID - A unique ID starting at 0 that we assign (not included in RepeatMasker output).
                  Root has ID 0.
   ClassNode.Class - This ClassNode's name cut on "/".
                     This likely isn't necessary, and may be removed in the future.
   ClassNode.Parent - A pointer to this ClassNode's parent in the ancestry tree.
                      It should be nil for root and only for root.
   ClassNode.Children - A slice containing pointers to all of this ClassNode's children in the tree.
   ClassNode.Repeat - A pointer to this ClassNode's corresonding Repeat, if it has one.
                      This field is of dubious value.
*/
type ClassNode struct {
    Name     string
    ID       uint16
    Class    []string
    Parent   *ClassNode
    Children []*ClassNode
    Repeat   *Repeat
}

type ClassNodes []*ClassNode

/*
   ClassTree.ClassNodes - Maps a fully qualified class name (excluding root) to that class's ClassNode struct, if it exists.
                          This is slower than ClassTree.NodesByID, and should only be used when necessary.
   ClassTree.NodesByID - A slice of pointers to all ClassNode structs, indexed by ID.
                         This should be the default means of accessing a ClassNode.
   ClassTree.Root - A pointer to the ClassTree's root, which has name "root" and ID 0.
                    We explicitly create this - it isn't present in the RepeatMasker output.
*/
type ClassTree struct {
    ClassNodes map[string](*ClassNode)
    NodesByID  []*ClassNode
    Root       *ClassNode
}

/* A two-bits-per-base sequence of up to 31 bases, with low-order bits occupied first.
   00 = 'a'
   01 = 'c'
   10 = 'g'
   11 = 't'
*/
type KmerInt uint64
type KmerInts []KmerInt

// A two-bits-per-base sequence of up to 15 bases, with low-bits occupied first.
type MinInt uint32
type MinInts []MinInt

/*
   Indexes the base of a kmer that is the starting index of its minimizer.
   If less than 32, the minimizer is the positive strand representation
   Otherwise, the minimizer is the reverse complement of kmer[minkey%32 : minkey + (m%32)]
*/
type MinKey uint8

/*
   This is what is stored by the main Kraken data structure: RepeatGenome.Kmers
   The first eight bits are the integer representation of the kmer's sequence (type KmerInt).
   The last two are the LCA ID (type uint16).
*/
type Kmer [10]byte
type Kmers []Kmer
type PKmers []*Kmer

/*
   Each base is represented by two bits.
   High-order bits are occupied first.
   Remember that Seq.Len is the number of bases contained, while len(Seq.Bytes) is the number of bytes necessary to represent them.
*/
type Seq struct {
    Bytes []byte
    Len   uint64
}

type Seqs []Seq

// A 2-dimensional map used to represent a newly-parsed FASTA-formatted reference genome.
type Chroms map[string](map[string]TextSeq)

// parseMatches() is kept as a solo function rather than a method of type RepeatGenome.
// This allows the function to be as portable as possible, which is important because some may come here wanting RepeatMasker-output-parsing code and nothing else.
func parseMatches(genomeName string) (error, Matches) {
    // "my_genome_name"  ->  "my_genome_name/my_genome_name.fa.out"
    filepath := strings.Join([]string{genomeName, "/", genomeName, ".fa.out"}, "")

    err, matchLines := fileLines(filepath)
    if err != nil {
        return ParseError{"repeatgenome.parseMatches()", filepath, err}, nil
    }

    matchLines = matchLines[3:] // drop the header lines

    var matches Matches

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
        sw_Score, err := strconv.ParseInt(rawVals[0], 10, 32)
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
            var seqName string = "" // forward initialization necessary
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
    minSliceSize = 1 << (2 * m)
    kMask = (1 << (2 * k)) - 1
    mMask = (1 << (2 * m)) - 1
    fmt.Println("k =", k)
    fmt.Println("m =", m)

    runtime.GOMAXPROCS(numCPU)

    // we popoulate the RepeatGenome mostly with helper functions
    // we should consider whether it makes more sense for them to alter the object directly, than to return their results
    rg := new(RepeatGenome)
    rg.Name = config.Name

    var err error

    if config.MemProfile {
        os.Mkdir("profiles", os.ModeDir)
        memProfFile, err = os.Create("profiles/" + rg.Name + ".memprof")
        if err != nil {
            return IOError{"RepeatGenome.getKrakenSlice()", err}, nil
        }
        pprof.WriteHeapProfile(memProfFile)
        defer memProfFile.Close()
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
        minsFile, err := os.OpenFile(rg.Name+".mins", os.O_RDONLY, 0400)

        if err == nil {    // implies the file already exists
            fmt.Println("\nKraken library file exists - using contents")
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

    if config.WriteStats {
        err := rg.WriteStatData()
        if err != nil {
            return err, nil
        }

        err = rg.WriteClassJSON(false, false)
        if err != nil {
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
            repeat := Repeat{
                ID:        uint64(len(rg.Repeats)),
                ClassList: match.RepeatClass,
                Name:      match.RepeatName,
                Instances: []*Match{match},
            }

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
    tree.Root = &ClassNode{
        Name:     "root",
        ID:       0,
        Class:    []string{"root"},
        Parent:   nil,
        Children: nil,
        Repeat:   nil,
    }
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
