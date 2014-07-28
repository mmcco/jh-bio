package repeatgenome

import (
    "bufio"
    "bytes"
    "encoding/binary"
    "encoding/json"
    "fmt"
    "io"
    "os"
    "reflect"
    "sort"
    "strconv"
    "strings"
    "unsafe"
)

var magicBytes [3]byte = [3]byte{'\x49', '\x47', '\x4d'}

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

func (rg *RepeatGenome) WriteMins() error {
    filename := strings.Join([]string{rg.Name, ".mins"}, "")
    outfile, err := os.Create(filename)
    if err != nil {
        return err
    }
    defer outfile.Close()
    writer := bufio.NewWriter(outfile)
    defer writer.Flush()

    kmerBuf := make(TextSeq, k, k)
    minBuf := make(TextSeq, m, m)

    for i, thisMin := range rg.SortedMins {
        fillMinBuf(minBuf, thisMin)
        _, err = fmt.Fprintf(writer, ">%s\n", minBuf)
        if err != nil {
            return err
        }
        for _, kmer := range rg.getMinsKmers(uint64(i)) {
            kmerSeqInt := *(*KmerInt)(unsafe.Pointer(&kmer[0]))
            lcaID := *(*uint16)(unsafe.Pointer(&kmer[8]))
            fillKmerBuf(kmerBuf, kmerSeqInt)
            _, err = fmt.Fprintf(writer, "\t%s %s\n", kmerBuf, rg.ClassTree.NodesByID[lcaID].Name)
            if err != nil {
                return err
            }
        }
    }
    return nil
}

func (kmerInt KmerInt) print() {
    var i uint8
    for i = 0; i < k; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (kmerInt >> (2 * (k - i - 1))) & 3 {
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

func (minInt MinInt) print() {
    var i uint8
    for i = 0; i < m; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (minInt >> (2 * (m - i - 1))) & 3 {
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

func writeKmerInt(writer io.ByteWriter, seqInt KmerInt) error {
    var i uint8
    var err error
    for i = 0; i < k; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (k - i - 1))) & 3 {
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

func writeMinInt(writer io.ByteWriter, seqInt MinInt) error {
    var i uint8
    var err error
    for i = 0; i < m; i++ {
        // this tricky bit arithmetic shifts the two bits of interests to the two rightmost positions, then selects them with the and statement
        switch (seqInt >> (2 * (m - i - 1))) & 3 {
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
func fillKmerBuf(slice TextSeq, seqInt KmerInt) {
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

func fillMinBuf(slice TextSeq, seqInt MinInt) {
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

func (rg *RepeatGenome) WriteFullKmers(filepath string) error {
    outfile, err := os.Create(filepath)
    if err != nil {
        return IOError{"FullKmers.WriteFullKmers()", err}
    }
    defer outfile.Close()

    // write magic bytes, k, m, numMins, and numKmers
    for _, val := range [...]interface{}{magicBytes, k, m, uint64(len(rg.SortedMins)), uint64(len(rg.FullKmers))} {
        err := binary.Write(outfile, binary.LittleEndian, val)
        if err != nil {
            return ParseError{"RepeatGenome.WriteFullKmers()", filepath, err}
        }
    }

    fmt.Println("writing", comma(uint64(len(rg.FullKmers))), "kmers")

    for _, fullKmer := range rg.FullKmers {
        err := binary.Write(outfile, binary.LittleEndian, fullKmer)
        if err != nil {
            return ParseError{"RepeatGenome.WriteFullKmers()", filepath, err}
        }
    }

    return nil
}

func (rg *RepeatGenome) ReadFullKmers(filepath string) error {
    outfile, err := os.OpenFile(filepath, os.O_RDONLY, 0400)
    if err != nil {
        return IOError{"FullKmers.ReadFullKmers()", err}
    }
    defer outfile.Close()

    headerBuf := make([]byte, 21, 21)
    bytesRead, err := outfile.Read(headerBuf)
    if err == io.EOF || bytesRead < len(headerBuf) {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("file too short")}
    }
    if err != nil {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, err}
    }
    if !reflect.DeepEqual(magicBytes[:], headerBuf[:3]) {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("magic bytes incorrect - are you sure this is the right file?")}
    }

    k_, bytesRead := binary.Uvarint(headerBuf[3:4])
    if bytesRead != 1 {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("error parsing k")}
    }
    if uint8(k_) != k {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("k value inconsistent with the supplied k value")}
    }

    m_, bytesRead := binary.Uvarint(headerBuf[4:5])
    if bytesRead != 1 {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("error parsing m")}
    }
    if uint8(m_) != m {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("m value inconsistent with the supplied m value")}
    }

    numMins, bytesRead := binary.Uvarint(headerBuf[5:13])
    if bytesRead != 8 {
        return ParseError{"FullMins.ReadFullMins()", filepath, fmt.Errorf("error parsing numMins")}
    }
    fmt.Println("expecting to read", comma(numMins), "mins")

    numKmers, bytesRead := binary.Uvarint(headerBuf[13:21])
    if bytesRead != 8 {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("error parsing numKmers")}
    }
    fmt.Println("expecting to read", comma(numKmers), "kmers")

    if len(rg.FullKmers) > 0 {
        fmt.Println("WARNING: RepeatGenome.ReadFullKmers() overwriting RepeatGenome.FullKmers"); fmt.Println()
    }
    rg.FullKmers = make(FullKmers, numKmers, numKmers)

    kmerBuf := make([]byte, 14)
    var fullKmer FullKmer
    var i uint64
    for i = 0; i < numKmers; i++ {
        bytesRead, err := outfile.Read(kmerBuf)
        if err == io.EOF || bytesRead < len(kmerBuf) {
            return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("fewer kmers supplied that specified in the minimizer header")}
        }
        if err != nil {
            return ParseError{"FullKmers.ReadFullKmers()", filepath, err}
        }

        copy(fullKmer[:], kmerBuf)
        rg.FullKmers = append(rg.FullKmers, fullKmer)
    }

    bytesRead, err = outfile.Read(kmerBuf)
    if err != io.EOF {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("extra data in file")}
    }

    if len(rg.SortedMins) > 0 {
        fmt.Println("WARNING: RepeatGenome.ReadFullKmers() overwriting rg.SortedMins")
        rg.SortedMins = rg.SortedMins[:0]
    }
    if len(rg.OffsetsToMin) > 0 {
        fmt.Println("WARNING: RepeatGenome.ReadFullKmers() overwriting rg.OffsetsToMin")
        rg.OffsetsToMin = rg.OffsetsToMin[:0]
    }

    for i, fullKmer := range rg.FullKmers {
        minInt := *(*MinInt)(unsafe.Pointer(&fullKmer[8]))
        if len(rg.SortedMins) == 0 || minInt != rg.SortedMins[len(rg.SortedMins)-1] {
            rg.SortedMins = append(rg.SortedMins, minInt)
            rg.OffsetsToMin = append(rg.OffsetsToMin, uint64(i))
        }
    }

    if !sort.IsSorted(rg.SortedMins) {
        return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("minimizers not sorted")}
    }

    for i, offset := range rg.OffsetsToMin {
        var end uint64
        if i == len(rg.OffsetsToMin) - 1 {
            end = uint64(len(rg.FullKmers))
        } else {
            end = rg.OffsetsToMin[i+1]
        }
        if !sort.IsSorted(rg.FullKmers[offset : end]) {
            return ParseError{"FullKmers.ReadFullKmers()", filepath, fmt.Errorf("kmers not sorted")}
        }
    }

    fmt.Println("read", comma(uint64(len(rg.FullKmers))), "kmers")
    fmt.Println("read", comma(uint64(len(rg.SortedMins))), "minimizers")

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

func readSimSeqReads(filepath string) (error, Seqs) {
    err, lines := fileLines(filepath)
    if err != nil {
        return err, nil
    }

    simReads := make(Seqs, 0, len(lines))
    for _, line := range lines {
        simReads = append(simReads, TextSeq(line).Seq())
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

func parseReadSAMs(filepath string) (error, []ReadSAM) {
    err, lines := fileLines(filepath)
    if err != nil {
        return err, nil
    }

    if len(lines) < 3 {
        err := ParseError{"parseReadSAMs()", filepath, fmt.Errorf("only %d lines in SAM file - need at least three", len(lines))}
        return err, nil
    }

    // drop header
    lines = lines[3:]
    readSAMs := make([]ReadSAM, len(lines), len(lines))

    for i, line := range lines {
        fields := bytes.Fields(line)
        if len(fields) != 12 {
            err := ParseError{"parseReadSAMs()", filepath, fmt.Errorf("SAM line has %d fields - 12 expected", len(fields))}
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

// passes all file names in the dir to parseReadSAMs and returns the concatenated results
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
        if len(fileinfo.Name()) > 10 && fileinfo.Name()[len(fileinfo.Name())-10 : ] == ".fasta.sam" {
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

func (rg *RepeatGenome) parseLib(filepath string) error {
    err, lines := fileLines(filepath)
    if err != nil {
        return err
    }

    var minimizer MinInt
    var minSet bool = false

    for _, line := range lines {
        if line[0] == '>' {
            minimizer = TextSeq(bytes.TrimSpace(line[1:])).minInt()
            if len(rg.SortedMins) > 0 && minimizer > rg.SortedMins[len(rg.SortedMins)-1] {
                return IOError{"RepeatGenome.parseLib()", fmt.Errorf("minimizers not written in sorted order")}
            }
            rg.SortedMins = append(rg.SortedMins, minimizer)
            rg.OffsetsToMin = append(rg.OffsetsToMin, uint64(len(rg.Kmers)))
            minSet = true
        } else if minSet {
            return IOError{"RepeatGenome.parseLib()}", fmt.Errorf("missing or empty minimizer")}
        } else {
            fields := bytes.Fields(line)
            kmerInt := TextSeq(fields[0]).kmerInt()
            if len(rg.Kmers) > 0 && kmerInt > *(*KmerInt)(unsafe.Pointer(&rg.Kmers[len(rg.Kmers)-1][0])) {
                return IOError{"RepeatGenome.parseLib()", fmt.Errorf("kmers not written in sorted order")}
            }
            lca := rg.ClassTree.ClassNodes[string(fields[1])]
            var kmer Kmer
            *(*KmerInt)(unsafe.Pointer(&kmer[0])) = kmerInt
            *(*uint16)(unsafe.Pointer(&kmer[8])) = lca.ID
            rg.Kmers = append(rg.Kmers, kmer)
        }
    }
    return nil
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
