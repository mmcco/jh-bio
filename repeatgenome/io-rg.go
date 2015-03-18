/*
   A set of functions and datatypes used for reading input data, Kraken library writing, and JSON data writing.
*/

package repeatgenome

import (
    "bufio"
    "bytes"
    "encoding/binary"
    "encoding/json"
    "fmt"
    "github.com/plsql/jh-bio/bioutils"
    "io"
    "os"
    "strconv"
    "unsafe"
)

// A magic value, used to give an explicit, declarative error message
// if the user tries to parse a non-Kraken file.
const magicVal int = 13377331

// Used only for recursively writing the JSON representation of the
// ClassTree.
type JSONNode struct {
    Name     string      `json:"name"`
    Size     uint64      `json:"size"`
    Children []*JSONNode `json:"children"`
}

/*
   Writes a JSON representation of the class tree.
   Used by the Javascript visualization, among other things.
   Currently, each node is associated with a value "size", the number of kmers associated with it.
   useCumSize determines whether the kmer count is cumulative, counting all kmers in its subtree.
*/
func (rg *RepeatGenome) WriteClassJSON(useCumSize, printLeaves bool) error {
    tree := &rg.ClassTree

    statDirName := rg.Name + "-stats"
    err := os.Mkdir(statDirName, os.ModePerm)
    if err != nil && !os.IsExist(err) {
        return fmt.Errorf("RepeatGenome.WriteClassJSON():" + err.Error())
    }
    filename := statDirName + "/" + rg.Name + ".classtree.json"
    outfile, err := os.Create(filename)
    defer outfile.Close()
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteClassJSON():" + err.Error())
    }

    classToCount := make(map[ClassID]uint64)
    for i := range rg.Kmers {
        classID := rg.Kmers[i].ClassID()
        classToCount[classID]++
    }

    root := JSONNode{tree.Root.Name, classToCount[0], nil}
    tree.jsonRecPopulate(&root, classToCount)

    if useCumSize {
        root.jsonRecSize()
    }

    if !printLeaves {
        root.deleteLeaves()
    }

    jsonBytes, err := json.MarshalIndent(root, "", "\t")
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteClassJSON():" + err.Error())
    }
    fmt.Fprint(outfile, string(jsonBytes))

    return nil
}

/*
   Recursively populates a JSONNode tree with its kmer counts.
*/
func (classTree *ClassTree) jsonRecPopulate(jsonNode *JSONNode, classToCount map[ClassID]uint64) {
    classNode := classTree.ClassNodes[jsonNode.Name]
    for i := range classNode.Children {
        child := new(JSONNode)
        child.Name = classNode.Children[i].Name
        child.Size = classToCount[classTree.ClassNodes[child.Name].ID]
        jsonNode.Children = append(jsonNode.Children, child)
        classTree.jsonRecPopulate(child, classToCount)
    }
}

/*
   Populates each node in the JSONNode tree with its cumulative kmer count.
*/
func (jsonNode *JSONNode) jsonRecSize() uint64 {
    for i := range jsonNode.Children {
        jsonNode.Size += jsonNode.Children[i].jsonRecSize()
    }
    return jsonNode.Size
}

/*
   Deletes the leaves from a JSONNode tree.
   Used to ignore leaves for visualizations.
*/
func (jsonNode *JSONNode) deleteLeaves() {
    branchChildren := []*JSONNode{}
    for i := range jsonNode.Children {
        child := jsonNode.Children[i]
        if child.Children != nil && len(child.Children) > 0 {
            branchChildren = append(branchChildren, child)
            child.deleteLeaves()
        }
    }
    jsonNode.Children = branchChildren
}

func (refGenome *RepeatGenome) PrintChromInfo() {
    fmt.Println()
    for k, v := range refGenome.chroms {
        for k_, v_ := range v {
            fmt.Printf("refGenome.chroms[%s][%s] = %s . . . %s\n", k, k_, v_[:10], v_[len(v_)-10:])
            fmt.Printf("len(refGenome.chroms[%s][%s]) = %d\n", k, k_, len(v_))
        }
        fmt.Println()
    }
}

// assumes that all bytes in the slice to be filled are initialized
// (a.k.a initialize buffer with make([]byte, k, k))
func fillKmerBuf(slice []byte, seqInt uint64) {
    if len(slice) > 31 {
        panic("slice of length greater than 31 passed to fillKmerBuf()")
    }
    for i := range slice {
        switch (seqInt >> uint((2 * (len(slice) - i - 1)))) & 3 {
        case 0:
            slice[i] = byte('a')
            break
        case 1:
            slice[i] = byte('c')
            break
        case 2:
            slice[i] = byte('g')
            break
        case 3:
            slice[i] = byte('t')
            break
        default:
            panic("error in fillKmerBuf() base selection")
        }
    }
}

func fillMinBuf(slice []byte, seqInt uint32) {
    if len(slice) > 15 {
        panic("slice of length greater than 15 passed to fillMinBuf()")
    }
    for i := range slice {
        switch (seqInt >> uint((2 * (len(slice) - i - 1)))) & 3 {
        case 0:
            slice[i] = byte('a')
            break
        case 1:
            slice[i] = byte('c')
            break
        case 2:
            slice[i] = byte('g')
            break
        case 3:
            slice[i] = byte('t')
            break
        default:
            panic("error in fillMinBuf() base selection")
        }
    }
}

func (rg *RepeatGenome) WriteStatData() error {
    statDirName := rg.Name + "-stats/"
    err := os.Mkdir(statDirName, os.ModePerm)
    if err != nil && !os.IsExist(err) {
        return fmt.Errorf("RepeatGenome.WriteStatData():" + err.Error())
    }

    err = rg.Repeats.Write(statDirName + "repeats.txt")
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteStatData():" + err.Error())
    }

    err = rg.Matches.Write(statDirName + "matches.txt")
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteStatData():" + err.Error())
    }

    err = ClassNodes(rg.ClassTree.NodesByID).Write(statDirName + "class-nodes.txt")
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteStatData():" + err.Error())
    }

    return nil
}

func (repeats Repeats) Write(filename string) error {
    outfile, err := os.Create(filename)
    if err != nil {
        return fmt.Errorf("Repeats.Write():" + err.Error())
    }
    defer outfile.Close()

    fmt.Fprintln(outfile, "ID\tName\tNumInstances\tDepth")
    fmt.Fprintln(outfile)
    for _, repeat := range repeats {
        var depth int
        if repeat.Name == "root" {
            depth = 0
        } else {
            depth = len(repeat.ClassList)
        }
        fmt.Fprintf(outfile, "%d\t%s\t%d\t%d", repeat.ID, repeat.Name, len(repeat.Instances), depth)
        fmt.Fprintln(outfile)
    }

    return nil
}

func (classNodes ClassNodes) Write(filename string) error {
    outfile, err := os.Create(filename)
    if err != nil {
        return fmt.Errorf("ClassNodes.Write():" + err.Error())
    }
    defer outfile.Close()

    fmt.Fprintln(outfile, "ID\tName\tDepth\tNumChildren\tIsRepeat")
    fmt.Fprintln(outfile)
    for _, cn := range classNodes {
        depth := 0
        walker := cn
        for walker.Parent != nil {
            depth++
            walker = walker.Parent
        }
        fmt.Fprintf(outfile, "%d\t%s\t%d\t%d\t%v", cn.ID, cn.Name, depth, len(cn.Children), cn.Repeat != nil)
        fmt.Fprintln(outfile)
    }

    return nil
}

func (rg *RepeatGenome) WriteKraken() error {
    buf := make([]byte, 10) // varints can occupy at most 10 bytes
    libDirName := rg.Name + "-lib"
    err := os.Mkdir(libDirName, os.ModePerm)
    if err != nil && !os.IsExist(err) {
        return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
    }
    outfile, err := os.Create(libDirName + "/" + rg.Name + ".kraken")
    if err != nil {
        return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
    }
    defer outfile.Close()

    // write magic bytes, k, m, numMins, and numKmers
    write_vals := []interface{}{
        magicVal,
        int(bioutils.K),
        int(bioutils.M),
        len(rg.SortedMins),
        len(rg.Kmers),
    }

    for _, val := range write_vals {
        numBytes := binary.PutUvarint(buf, uint64(val.(int)))
        bytesWritten, err := outfile.Write(buf[:numBytes])
        if bytesWritten != numBytes {
            return fmt.Errorf("RepeatGenome.WriteKraken(): did not write expected number of bytes")
        }
        if err != nil {
            return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
        }
    }

    // write mins in order
    for _, minInt := range rg.SortedMins {
        numBytes := binary.PutUvarint(buf, uint64(minInt))
        bytesWritten, err := outfile.Write(buf[:numBytes])
        if bytesWritten != numBytes {
            return fmt.Errorf("RepeatGenome.WriteKraken(): did not write expected number of bytes")
        }
        if err != nil {
            return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
        }
    }

    // write their counts in order
    for _, minInt := range rg.SortedMins {
        numBytes := binary.PutUvarint(buf, uint64(rg.MinCounts[minInt]))
        bytesWritten, err := outfile.Write(buf[:numBytes])
        if bytesWritten != numBytes {
            return fmt.Errorf("RepeatGenome.WriteKraken(): did not write expected number of bytes")
        }
        if err != nil {
            return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
        }
    }

    // write kmers in order
    kmerSize := len(Kmer{})
    var kmerBuf []byte
    for _, kmer := range rg.Kmers {
        kmerBuf = kmer[:]
        numBytes, err := outfile.Write(kmerBuf)
        if numBytes != kmerSize {
            return fmt.Errorf("RepeatGenome.WriteKraken(): wrote too few kmer bytes")
        }
        if err != nil {
            return fmt.Errorf("RepeatGenome.WriteKraken():" + err.Error())
        }
    }

    return nil
}

// has a lot of error handling, but pretty simple logic
func (rg *RepeatGenome) ReadKraken(infile *os.File) error {
    bufioReader := bufio.NewReader(infile)
    defer infile.Close()

    // check magic value
    thisMagicVal, err := binary.ReadUvarint(bufioReader)
    if err != nil {
        return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
    }
    if int(thisMagicVal) != magicVal {
        return fmt.Errorf("RepeatGenome.ReadKraken(): incorrect magic value - are you sure" + infile.Name() + "is the right file?")
    }

    // check K
    thisK, err := binary.ReadUvarint(bufioReader)
    if err != nil {
        return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
    }
    if uint8(thisK) != bioutils.K {
        return fmt.Errorf("RepeatGenome.ReadKraken(): incompatible k value - are you sure" + infile.Name() + "is the right file?")
    }

    // check M
    thisM, err := binary.ReadUvarint(bufioReader)
    if err != nil {
        return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
    }
    if uint8(thisM) != bioutils.M {
        return fmt.Errorf("RepeatGenome.ReadKraken(): incompatible m value - are you sure" + infile.Name() + "is the right file?")
    }

    // get number of minimizers
    numMins, err := binary.ReadUvarint(bufioReader)
    if err != nil {
        return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
    }

    // get number of kmers
    numKmers, err := binary.ReadUvarint(bufioReader)
    if err != nil {
        return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
    }

    // populate rg.SortedMins
    if len(rg.SortedMins) > 0 {
        fmt.Println("!!! WARNING !!! RepeatGenome.ReadKraken() overwriting RepeatGenome.SortedMins")
    }
    rg.SortedMins = make(MinInts, 0, numMins)

    var i uint64
    for i = 0; i < numMins; i++ {
        minInt, err := binary.ReadUvarint(bufioReader)
        if err != nil {
            return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
        }
        rg.SortedMins = append(rg.SortedMins, uint32(minInt))
    }

    // populate rg.MinCounts
    if rg.MinCounts != nil && len(rg.MinCounts) > 0 {
        fmt.Println("!!! WARNING !!! RepeatGenome.ReadKraken() overwriting RepeatGenome.MinCounts")
    }
    rg.MinCounts = make([]uint32, minSliceSize, minSliceSize)

    for _, minInt := range rg.SortedMins {
        minCount, err := binary.ReadUvarint(bufioReader)
        if err != nil {
            return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
        }
        rg.MinCounts[minInt] = uint32(minCount)
    }

    // populate rg.Kmers
    if len(rg.Kmers) > 0 {
        fmt.Println("!!! WARNING !!! RepeatGenome.ReadKraken() overwriting RepeatGenome.Kmers")
    }
    rg.Kmers = make(Kmers, 0, numKmers)
    kmerBuf := make([]byte, len(Kmer{}))

    for i = 0; i < numKmers; i++ {
        numBytes, err := io.ReadFull(bufioReader, kmerBuf)
        if err != nil {
            return fmt.Errorf("RepeatGenome.ReadKraken():" + err.Error())
        }
        if numBytes != len(kmerBuf) {
            return fmt.Errorf("RepeatGenome.ReadKraken(): only %d kmer bytes read instead of expected %d - on kmer index %d - file %s incorrectly written or corrupted",
                numBytes, len(kmerBuf), i, infile.Name())
        }

        rg.Kmers = append(rg.Kmers, *(*Kmer)(unsafe.Pointer(&kmerBuf[0])))
    }

    rg.populateMinOffsets()

    return nil
}

func (repeat *Repeat) Print() {
    for j := range repeat.ClassList {
        for j_ := 0; j_ < j; j_++ {
            fmt.Printf("\t")
        }
        fmt.Printf("%s\n", repeat.ClassList[j])
    }
    fmt.Println()
}

func (classTree *ClassTree) PrintTree() {
    classTree.Root.printTreeRec(0, true)
}

// Doesn't print leaves.
// Prevents the terminal from being flooded with Unknowns, Others, and
// Simple Repeats.
func (classTree *ClassTree) PrintBranches() {
    classTree.Root.printTreeRec(0, false)
}

func (classNode *ClassNode) printTreeRec(indent int, printLeaves bool) {
    for i := 0; i < indent; i++ {
        fmt.Printf("\t")
    }
    fmt.Println(classNode.Class[len(classNode.Class)-1])
    for i := range classNode.Children {
        if printLeaves || len(classNode.Children[i].Children) > 0 {
            classNode.Children[i].printTreeRec(indent+1, printLeaves)
        }
    }
}

func readSimSeqReads(filepath string) (error, Seqs) {
    err, lines := fileLines(filepath)
    if err != nil {
        return err, nil
    }

    simReads := make(Seqs, 0, len(lines))
    for _, line := range lines {
        simReads = append(simReads, GetSeq(line))
    }
    return nil, simReads
}

type ReadSAM struct {
    TextSeq []byte
    SeqName  string
    StartInd uint64
}

type ReadSAMResponse struct {
    ReadSAM   ReadSAM
    ClassNode *ClassNode
}

func parseReadSAMs(filepath string) (error, []ReadSAM) {
    err, lines := fileLines(filepath)
    if err != nil {
        return err, nil
    }

    if len(lines) < 3 {
        err := fmt.Errorf("repeatgenome.parseReadSAMs(): only %d lines in SAM file %s - need at least three", len(lines), filepath)
        return err, nil
    }

    // drop header
    lines = lines[3:]
    readSAMs := make([]ReadSAM, len(lines), len(lines))

    for i, line := range lines {
        fields := bytes.Fields(line)
        if len(fields) != 12 {
            err := fmt.Errorf("parseReadSAMs(): SAM line in %s has %d fields - 12 expected", filepath, len(fields))
            return err, nil
        }

        readSAMs[i].TextSeq = bytes.ToLower(fields[9])
        readSAMs[i].SeqName = string(fields[2])
        readSAMs[i].StartInd, err = strconv.ParseUint(string(fields[3]), 10, 64)
        if err != nil {
            return err, nil
        }
        // start index is 1 in the SAM standard
        readSAMs[i].StartInd--
    }

    return nil, readSAMs
}

// Passes all file names in the dir to parseReadSAMs and returns the
// concatenated results.
func GetReadSAMs(readsDirPath string) (error, []ReadSAM) {
    currDir, err := os.Open(readsDirPath)
    if err != nil {
        return err, nil
    }

    fileinfos, err := currDir.Readdir(-1)
    if err != nil {
        return err, nil
    }

    var samFiles []os.FileInfo
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 10 && fileinfo.Name()[len(fileinfo.Name())-10:] == ".fasta.sam" {
            samFiles = append(samFiles, fileinfo)
        }
    }
    readSAMs := []ReadSAM{}
    for _, samFile := range samFiles {
        err, theseReadSAMs := parseReadSAMs(readsDirPath + "/" + samFile.Name())
        if err != nil {
            return err, nil
        }
        for _, readSAM := range theseReadSAMs {
            readSAMs = append(readSAMs, readSAM)
        }
    }

    return nil, readSAMs
}

func (seq Seq) Print() {
    var i uint64
    for i = 0; i < seq.Len; i++ {
        thisByte := seq.Bytes[i/4]
        switch (thisByte >> (6 - 2*(i%4))) & 3 {
        case 0:
            fmt.Print("a")
            break
        case 1:
            fmt.Print("c")
            break
        case 2:
            fmt.Print("g")
            break
        case 3:
            fmt.Print("t")
            break
        default:
            panic("Seq.Print(): logic error in switch statement")
        }
    }
}
