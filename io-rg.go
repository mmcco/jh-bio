package repeatgenome

import (
    "bufio"
    "bytes"
    "encoding/json"
    "fmt"
    "io"
    "os"
    "strconv"
    "strings"
    "unsafe"
)

// used only for recursive JSON writing
type JSONNode struct {
    Name     string      `json:"name"`
    Size     uint64      `json:"size"`
    Children []*JSONNode `json:"children"`
}

func (repeatGenome *RepeatGenome) WriteClassJSON(useCumSize, printLeaves bool) error {
    tree := &repeatGenome.ClassTree

    filename := repeatGenome.Name + ".classtree.json"
    outfile, err := os.Create(filename)
    defer outfile.Close()
    if err != nil {
        return IOError{"RepeatGenome.WriteClassJSON()", err}
    }

    classToCount := make(map[uint16]uint64)
    for i := range repeatGenome.Kmers {
        lcaID := *(*uint16)(unsafe.Pointer(&repeatGenome.Kmers[i][8]))
        classToCount[lcaID]++
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
        return IOError{"RepeatGenome.WriteClassJSON()", err}
    }
    fmt.Fprint(outfile, string(jsonBytes))

    return nil
}

func (classTree *ClassTree) jsonRecPopulate(jsonNode *JSONNode, classToCount map[uint16]uint64) {
    classNode := classTree.ClassNodes[jsonNode.Name]
    for i := range classNode.Children {
        child := new(JSONNode)
        child.Name = classNode.Children[i].Name
        child.Size = classToCount[classTree.ClassNodes[child.Name].ID]
        jsonNode.Children = append(jsonNode.Children, child)
        classTree.jsonRecPopulate(child, classToCount)
    }
}

func (jsonNode *JSONNode) jsonRecSize() uint64 {
    for i := range jsonNode.Children {
        jsonNode.Size += jsonNode.Children[i].jsonRecSize()
    }
    return jsonNode.Size
}

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

func (repeatGenome *RepeatGenome) WriteMins(minMap MinMap) error {
    k := repeatGenome.K
    m := repeatGenome.M
    kmerBuf := make(TextSeq, k, k)
    minBuf := make(TextSeq, m, m)
    filename := strings.Join([]string{repeatGenome.Name, ".mins"}, "")
    outfile, err := os.Create(filename)
    if err != nil {
        return err
    }
    defer outfile.Close()
    writer := bufio.NewWriter(outfile)
    defer writer.Flush()

    var kmers PKmers
    var thisMin, kmerSeqInt uint64
    var lcaID uint16

    for thisMin, kmers = range minMap {
        fillKmerBuf(minBuf, thisMin)
        _, err = fmt.Fprintf(writer, ">%s\n", minBuf)
        if err != nil {
            return err
        }

        for i := range kmers {
            kmerSeqInt = *(*uint64)(unsafe.Pointer(&kmers[i][0]))
            lcaID = *(*uint16)(unsafe.Pointer(&kmers[i][8]))
            fillKmerBuf(kmerBuf, kmerSeqInt)
            _, err = fmt.Fprintf(writer, "\t%s %s\n", kmerBuf, repeatGenome.ClassTree.NodesByID[lcaID].Name)
            if err != nil {
                return err
            }
        }
    }
    return nil
}

func printSeqInt(seqInt uint64, seqLen uint8) {
    var i uint8
    for i = 0; i < seqLen; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (seqLen - i - 1))) & 3 {
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
            panic("error in printSeqInt() base selection")
        }
    }
}

func writeSeqInt(writer io.ByteWriter, seqInt uint64, seqLen uint8) error {
    var i uint8
    var err error
    for i = 0; i < seqLen; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (seqLen - i - 1))) & 3 {
        case 0:
            err = writer.WriteByte('a')
            break
        case 1:
            err = writer.WriteByte('c')
            break
        case 2:
            err = writer.WriteByte('g')
            break
        case 3:
            err = writer.WriteByte('t')
            break
        default:
            err = fmt.Errorf("error in printSeqInt() base selection")
        }
        if err != nil {
            return IOError{"repeatgenome.writeSeqInt()", err}
        }
    }
    return nil
}

// assumes that all bytes in the slice to be filled are initialized
// (a.k.a initialize buffer with make(TextSeq, k, k))
func fillKmerBuf(slice TextSeq, seqInt uint64) {
    if len(slice) > 32 {
        panic("slice of length greater than 32 passed to fillKmerBuf()")
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
            panic("error in printSeqInt() base selection")
        }
    }
}

func (repeats *Repeats) Write(filename string) error {
    outfile, err := os.Create(filename)
    if err != nil {
        return IOError{"Repeats.Write()", err}
    }
    defer outfile.Close()

    for i := range *repeats {
        if int((*repeats)[i].ID) == i {
            fmt.Fprintf(outfile, "%d %s\n", (*repeats)[i].ID, (*repeats)[i].Name)
        }
    }
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

// doesn't print leaves
// prevents the terminal from being flooded with Unknowns, Others, and Simple Repeats
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

func (seq *Seq) Print() {
    var numBases uint64 = uint64(len(seq.Bases))
    basesInFirstByte := 1 + ((seq.Len - 1) % numBases)
    var byteCopy byte = seq.Bases[0]
    // we start by truncating the extra bits from the first byte and printing the rest
    byteCopy <<= (4 - basesInFirstByte) * 2
    var i uint64

    for i = 0; i < basesInFirstByte; i++ {
        // \xc0 is 11000000 in binary - it selects the first base
        switch byteCopy & byte('\xc0') {
        case 0:
            fmt.Print("a")
            break
        case byte('\x40'):
            fmt.Print("c")
            break
        case byte('\x80'):
            fmt.Print("g")
            break
        case byte('\xc0'):
            fmt.Print("t")
            break
        default:
            panic("Seq.Print(): bad value in switch")
        }
        byteCopy <<= 2
    }

    // we then print the bases in the rest of the bytes
    for i := 1; i < len(seq.Bases); i++ {
        byteCopy = seq.Bases[i]

        for j := 0; j < 4; j++ {
            switch byteCopy & byte('\xc0') {
            case 0:
                fmt.Print("a")
                break
            case byte('\x40'):
                fmt.Print("c")
                break
            case byte('\x80'):
                fmt.Print("g")
                break
            case byte('\xc0'):
                fmt.Print("t")
                break
            default:
                panic("Seq.Print(): bad value in switch")
            }

            byteCopy <<= 2
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
    Seq      []byte
    SeqName  string
    StartInd uint64
}

type ReadSAMResponse struct {
    ReadSAM ReadSAM
    ClassNode *ClassNode
}

func ParseReadSAMs(filepath string) (error, []ReadSAM) {
    err, lines := fileLines(filepath)
    if err != nil {
        return err, nil
    }

    if len(lines) < 3 {
        err := ParseError{"ParseReadSAMs()", filepath, fmt.Errorf("only %d lines in SAM file - need at least three", len(lines))}
        return err, nil
    }

    // drop header
    lines = lines[3:]
    readSAMs := make([]ReadSAM, len(lines), len(lines))

    for i, line := range lines {
        fields := bytes.Fields(line)
        if len(fields) != 12 {
            err := ParseError{"ParseReadSAMs()", filepath, fmt.Errorf("SAM line has %d fields - 12 expected", len(fields))}
            return err, nil
        }
        
        readSAMs[i].Seq = bytes.ToLower(fields[9])
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
