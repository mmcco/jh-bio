/*
   The core of the package, from which other functions are dispatched. Includes
   RepeatMasker-related data structure definitions and parsers, as well as the
   RepeatGenome definition and initializer.
*/

package repeatgenome

import (
	"bytes"
	"fmt"
	"github.com/mmcco/jh-bio/bioutils"
	"io/ioutil"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strings"
)

// minSliceSize is the size of a slice that contains an index for each possible
// mMer.
const minSliceSize uint64 = 1 << (2 * bioutils.M)

var debug bool

// This is used as a rudimentary way of determining how many goroutines to
// spawn in concurrent sections.
var numCPU = runtime.NumCPU()

// Used to write heap memory profile.
var memProfFile *os.File

// A value of type Config is passed to the New() function, which
// constructs and returns a new RepeatGenome.
type Config struct {
	Name       string
	Debug      bool
	CPUProfile bool
	MemProfile bool
	WriteLib   bool
	ForceGen   bool
	WriteStats bool
}

/*
   RepeatGenome.Name - The name of the reference genome, such as "dm3"
       or "hg38". This is used to name created directories, and to find
       directories and files that may be read from, such as a stored
       Kraken library and reference sequences.
   RepeatGenome.chroms - A 2-dimensional map mapping a chromosome name
       to a map of its sequence names to their sequences (in text form).
       Actual 2-dimensional mapping is currently impossible because of
       RepeatMasker's 1-dimensional output.
   RepeatGenome.Kmers - A slice of all Kmers, sorted primarily by
       minimizer and secondarily by lexicographical value.
   RepeatGenome.MinOffsets - Maps a minimizer to its offset in the
       Kmers slice, or -1 if no kmers of this minimizer were stored.
   RepeatGenome.MinCounts - Maps a minimizer to the number of stored
       kmers associated with it.
   RepeatGenome.SortedMins - A sorted slice of all minimizers of
       stored kmers.
   RepeatGenome.Matches - All matches, indexed by their assigned IDs.
   RepeatGenome.ClassTree - Contains all information used for LCA
       determination and read classification. It may eventually be
       collapsed into RepeatGenome, as accessing it is rather verbose.
   RepeatGenome.Repeats - A slice of all repeats, indexed by their
       assigned IDs.
   RepeatGenome.RepeatMap - Maps a fully qualified repeat name,
       excluding root, to its struct.
*/
type RepeatGenome struct {
	Name       string
	chroms     Chroms
	Kmers      Kmers
	MinOffsets []int64
	MinCounts  []uint32
	SortedMins MinInts
	Matches    bioutils.Matches
	ClassTree  ClassTree
	Repeats    Repeats
	RepeatMap  map[string]*Repeat
	// Maps a Match to its associated ClassNode.
	// Currently used in read classification logic.
	// Not exposed in API.
	matchNodes map[*bioutils.Match]*ClassNode
}

/*
   Repeat.ID - A unique ID that we assign (not included in
       RepeatMasker output). Because these are assigned in the order in
       which they are encountered in <genome name>.fa.out, they are not
       compatible across even different versions of the same reference
       genome. This may change.
   Repeat.Name - The repeat's fully qualified name, excluding root.
   Repeat.ClassList - A slice of this Repeat's class ancestry from the
       top of the tree down, excluding root.
   Repeat.ClassNode - A pointer to the ClassNode which corresponds to
       this repeat.
   Repeat.Instances - A slice of pointers to all matches that are
       instances of this repeat.
*/
type Repeat struct {
	ID        uint64
	Name      string
	ClassList []string
	ClassNode *ClassNode
	Instances []*bioutils.Match
}

type Repeats []*Repeat

/*
   ClassNode.Name - This ClassNode's fully qualified name, excluding
       root.
   ClassNode.ID - A unique ID starting at 0 that we assign (not
       included in RepeatMasker output). Root has ID 0.
   ClassNode.Class - This ClassNode's name cut on "/". This likely
       isn't necessary, and may be removed in the future.
   ClassNode.Parent - A pointer to this ClassNode's parent in the
       ancestry tree. It should be nil for root and only for root.
   ClassNode.Children - A slice containing pointers to all of this
       ClassNode's children in the tree.
   ClassNode.Repeat - A pointer to this ClassNode's corresonding
       Repeat, if it has one. This field is of dubious value.
*/
type ClassNode struct {
	Name     string
	ID       ClassID
	Class    []string
	Parent   *ClassNode
	Children []*ClassNode
	Repeat   *Repeat
}

type ClassNodes []*ClassNode

/*
   ClassTree.ClassNodes - Maps a fully qualified class name (excluding
       root) to that class's ClassNode struct, if it exists. This is
       slower than ClassTree.NodesByID, and should only be used when
       necessary.
   ClassTree.NodesByID - A slice of pointers to all ClassNode structs,
       indexed by ID. This should be the default means of accessing a
       ClassNode.
   ClassTree.Root - A pointer to the ClassTree's root, which has name
       "root" and ID 0. We explicitly create this - it isn't present in
       the RepeatMasker output.
*/
type ClassTree struct {
	ClassNodes map[string](*ClassNode)
	NodesByID  []*ClassNode
	Root       *ClassNode
}

/* A two-bits-per-base sequence of up to 31 bases, with low-order bits
    occupied first.
   00 = 'a'
   01 = 'c'
   10 = 'g'
   11 = 't'

   The definitions of KmerInt was previously here, but I reverted to uint64 for
   simplicity.
*/
type KmerInts []uint64

/* A two-bits-per-base sequence of up to 15 bases, with low-bits
   occupied first.

   The definitions of MinInt was previously here, but I reverted to uint32 for
   simplicity.
*/
type MinInts []uint32

/*
   Indexes the base of a kmer that is the starting index of its
   minimizer. If less than 32, the minimizer is the positive strand
   representation Otherwise, the minimizer is the reverse complement
   of kmer[minkey%32 : minkey + (m%32)]

   The definition of MinKey previously existed here, but I
   reverted to uint8 for simplicity.
*/

/*
   A type synonym representing a ClassNode by ID. Used to
   space-efficiently store a read's classification.
*/
type ClassID uint16

/*
   This is what is stored by the main Kraken data structure:
   RepeatGenome.Kmers The first eight bits are the integer
   representation of the kmer's sequence (type KmerInt). The last two
   are the class ID (type ClassID).
*/
type Kmer [10]byte
type Kmers []Kmer
type PKmers []*Kmer

/*
   Each base is represented by two bits. High-order bits are occupied
   first. Remember that Seq.Len is the number of bases contained,
   while len(Seq.Bytes) is the number of bytes necessary to represent
   them.
*/
type Seq struct {
	Bytes []byte
	Len   uint64
}

type Seqs []Seq

// A 2-dimensional map used to represent a newly-parsed
// FASTA-formatted reference genome.
type Chroms map[string](map[string][]byte)

func parseGenome(genomeName string) (error, Chroms) {
	chromFileInfos, err := ioutil.ReadDir(genomeName + "-fasta")
	if err != nil {
		return fmt.Errorf("repeatgenome.parseGenome():" + err.Error()), nil
	}
	warned := false
	chroms := make(map[string](map[string][]byte))
	// used below to store the two keys for RepeatGenome.chroms
	for i := range chromFileInfos {
		chromFilename := chromFileInfos[i].Name()
		chromName := chromFilename[:len(chromFilename)-3]
		// "my_genome_name", "my_chrom_name"  ->  "my_genome_name/my_chrom_name"
		chromFilepath := strings.Join([]string{genomeName + "-fasta", chromFilename}, "/")
		// process the ref genome files (*.fa), not the repeat ref
		// files (*.fa.out and *.fa.align) or anything else
		infile, err := os.Open(chromFilepath)
		if err != nil {
			return err, nil
		}
		err, seqMap := bioutils.ReadFASTA(infile)
		if err != nil {
			return err, nil
		}

		for seqName := range seqMap {
			if warned {
				break
			}
			if seqName != chromName {
				fmt.Println("WARNING: reference genome is two-dimensional, containing sequences not named after their chromosome.")
				fmt.Println("Because RepeatMasker supplied only one-dimensional indexing, this may cause unexpected behavior or program failure.")
				fmt.Printf("seqName: %s\tlen(seqName): %d\n", seqName, len(seqName))
				fmt.Printf("chrom name: %s\tlen(chrom name): %d\n", chromName, len(chromName))
				warned = true
			}
		}

		chroms[chromName] = make(map[string][]byte)
		for seqName, seqBytes := range seqMap {
			chroms[chromName][seqName] = bytes.ToLower(seqBytes)
		}
	}
	return nil, chroms
}

func New(config Config) (error, *RepeatGenome) {
	debug = config.Debug

	runtime.GOMAXPROCS(numCPU)

	// We popoulate the RepeatGenome mostly with helper functions.
	// We should consider whether it makes more sense for them to
	// alter the object directly, than to return their results.
	rg := new(RepeatGenome)
	rg.Name = config.Name

	var err error

	if config.MemProfile {
		os.Mkdir("profiles", os.ModePerm)
		memProfFile, err = os.Create("profiles/" + rg.Name + ".memprof")
		if err != nil {
			return fmt.Errorf("repeatgenome.New():" + err.Error()), nil
		}
		pprof.WriteHeapProfile(memProfFile)
		defer memProfFile.Close()
	}

	err, rg.chroms = parseGenome(rg.Name)
	if err != nil {
		return err, nil
	}
	repeatMaskerFile, err := os.Open(rg.Name + "/" + rg.Name + ".fa.out")
	if err != nil {
		return err, nil
	}
	err, rg.Matches = bioutils.ParseMatches(repeatMaskerFile)
	if err != nil {
		return err, nil
	}
	rg.getRepeats()
	rg.getClassTree()

	if config.ForceGen {
		rg.genKrakenLibVerbose()
		err = rg.WriteKraken()
		if err != nil {
			return err, nil
		}
	} else {
		krakenFile, err := os.Open(rg.Name + "-lib/" + rg.Name + ".kraken")

		if err == nil { // implies the file already exists
			fmt.Println()
			fmt.Println("Kraken library file exists - using contents")
			err = rg.ReadKraken(krakenFile)
			if err != nil {
				return fmt.Errorf("repeatgenome.New():" + err.Error()), nil
			}
		} else if os.IsNotExist(err) { // the case that there isn't a written file yet
			fmt.Println("Kraken library file doesn't exist - generating library")
			fmt.Println()
			rg.genKrakenLibVerbose()
			err = rg.WriteKraken()
			if err != nil {
				return err, nil
			}
		} else { // otherwise we're dealing with a generic error of some sort
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
	// We first determine the necessary size of the slice - we can't
	// use append because matches are not sorted by repeatID.

	rg.RepeatMap = make(map[string]*Repeat)

	// DON'T use the second field of the range - this causes the Match
	// struct to be copied.
	// Creating an alias struct (match := rg.Matches[i]) of type Match
	// rather than *Match causes the repeat.Instance item to point to
	// a copy, not the original Match struct.
	for i := range rg.Matches {
		match := &rg.Matches[i]
		// don't bother overwriting
		if repeat, exists := rg.RepeatMap[match.RepeatName]; exists {
			repeat.Instances = append(repeat.Instances, match)
		} else {
			repeat := Repeat{
				ID:        uint64(len(rg.Repeats)),
				ClassList: match.RepeatClass,
				Name:      match.RepeatName,
				Instances: []*bioutils.Match{match},
			}

			rg.Repeats = append(rg.Repeats, &repeat)
			rg.RepeatMap[repeat.Name] = &repeat
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
		// process every heirarchy level (e.g. for "DNA/LINE/TiGGER",
		// process "DNA", then "DNA/LINE", then "DNA/LINE/TiGGER")
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
				classNode.ID = ClassID(len(tree.NodesByID))
				classNode.Class = thisClass
				if repeat, exists := rg.RepeatMap[thisClassName]; exists {
					classNode.Repeat = repeat
				}

				tree.ClassNodes[thisClassName] = classNode
				tree.NodesByID = append(tree.NodesByID, classNode)
				// first case handles primary classes, as root is
				// implicit and not listed in thisClass
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
		repeat := rg.Repeats[i]
		repeat.ClassNode = tree.ClassNodes[repeat.Name]

		// cond. mostly for debugging - remove
		if repeat.ClassNode == nil {
			fmt.Println(repeat.Name)
			log.Fatal("getClassTree(): nil Repeat.ClassNode")
		}
	}

	rg.matchNodes = make(map[*bioutils.Match]*ClassNode)
	// MUST NOT USE RANGE - the struct will be copied!
	for i := 0; i < len(rg.Matches); i++ {
		match := &rg.Matches[i]
		rg.matchNodes[match] = tree.ClassNodes[match.RepeatName]

		// cond. mostly for debugging - remove
		if rg.matchNodes[match] == nil {
			fmt.Println(match.RepeatName)
			log.Fatal("getClassTree(): nil bioutils.Match.ClassNode")
		}
	}
}
