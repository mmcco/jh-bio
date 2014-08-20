package repeatgenome

/*
   A collection of trivial functions used in other source files of the repeatgenome package.
*/

import (
    "bytes"
    "fmt"
    "io/ioutil"
    "strconv"
    "strings"
    "unsafe"
)

type ParseError struct {
    FuncName string
    Filepath string
    SubError error
}

type IOError struct {
    FuncName string
    SubError error
}

func (parseError ParseError) Error() string {
    if len(parseError.Filepath) == 0 {
        return fmt.Sprintf("%s: parsing error: %s", parseError.FuncName, parseError.SubError.Error())
    } else {
        return fmt.Sprintf("%s: parsing error in file %s: %s", parseError.FuncName, parseError.Filepath, parseError.SubError.Error())
    }
}

func (ioError IOError) Error() string {
    return fmt.Sprintf("IO error in %s: %s", ioError.FuncName, ioError.SubError.Error())
}

/*
   Courtesy of https://github.com/dustin/go-humanize
   Returns a string representing the int, with commas for readability.
*/
func comma(v int) string {
    sign := ""
    if v < 0 {
        sign = "-"
        v = 0 - v
    }

    parts := []string{"", "", "", "", "", "", "", ""}
    j := len(parts) - 1

    for v > 999 {
        parts[j] = strconv.FormatInt(int64(v%1000), 10)
        switch len(parts[j]) {
        case 2:
            parts[j] = "0" + parts[j]
        case 1:
            parts[j] = "00" + parts[j]
        }
        v = v / 1000
        j--
    }
    parts[j] = strconv.Itoa(v)
    return sign + strings.Join(parts[j:], ",")
}

/*
   Splits a byte slice on newlines, returning the resulting 2D slice.
*/
func lines(byteSlice []byte) [][]byte {
    var lines [][]byte = bytes.Split(byteSlice, []byte{'\n'})
    // drop trailing empty lines
    for len(lines) > 0 && len(lines[len(lines)-1]) == 0 {
        lines = lines[:len(lines)-1]
    }
    return lines
}

/*
   Reads the file at the supplied path and returns an error or its lines' bytes.
*/
func fileLines(filepath string) (err error, linesBytes [][]byte) {
    rawBytes, err := ioutil.ReadFile(filepath)
    if err != nil {
        return err, nil
    } else {
        return nil, lines(rawBytes)
    }
}

func max(a int, b int) int {
    if a > b {
        return a
    } else {
        return b
    }
}

func max32(a int32, b int32) int32 {
    if a > b {
        return a
    } else {
        return b
    }
}

func maxU32(a uint32, b uint32) uint32 {
    if a > b {
        return a
    } else {
        return b
    }
}

func max64(a int64, b int64) int64 {
    if a > b {
        return a
    } else {
        return b
    }
}

func min(a, b int) int {
    if a < b {
        return a
    } else {
        return b
    }
}

func minU64(a, b uint64) uint64 {
    if a < b {
        return a
    } else {
        return b
    }
}

// deprecated in favor of MinInt.canonicalRepr()
/*
func minMinInt(a, b MinInt) MinInt {
    if a < b {
        return a
    } else {
        return b
    }
}
*/

// deprecated in favor of KmerInt.canonicalRepr()
func minKmerInt(a, b KmerInt) KmerInt {
    if a < b {
        return a
    } else {
        return b
    }
}

// needed for sort.Interface
func (pkmers PKmers) Len() int {
    return len(pkmers)
}

func (pkmers PKmers) Swap(i, j int) {
    pkmers[i], pkmers[j] = pkmers[j], pkmers[i]
}

func (pkmers PKmers) Less(i, j int) bool {
    iVal := *(*uint64)(unsafe.Pointer(&pkmers[i][0]))
    jVal := *(*uint64)(unsafe.Pointer(&pkmers[j][0]))
    return iVal < jVal
}

func (kmers Kmers) Len() int {
    return len(kmers)
}

func (kmers Kmers) Swap(i, j int) {
    kmers[i], kmers[j] = kmers[j], kmers[i]
}

func (kmers Kmers) Less(i, j int) bool {
    return *(*KmerInt)(unsafe.Pointer(&kmers[i])) < *(*KmerInt)(unsafe.Pointer(&kmers[j]))
}

// sorts by both min and kmer - too slow
/*
func (kmers Kmers) Less(i, j int) bool {
    iVal := *(*KmerInt)(unsafe.Pointer(&kmers[i][0]))
    jVal := *(*KmerInt)(unsafe.Pointer(&kmers[j][0]))
    iMin := iVal.Minimize()
    jMin := jVal.Minimize()

    if iMin < jMin {
        return true
    } else if iMin > jMin {
        return false
    } else {
        return iVal < jVal
    }
}
*/

func (minInts MinInts) Len() int {
    return len(minInts)
}

func (minInts MinInts) Swap(i, j int) {
    minInts[i], minInts[j] = minInts[j], minInts[i]
}

func (minInts MinInts) Less(i, j int) bool {
    return minInts[i] < minInts[j]
}

func boolToInt(a bool) int {
    if a {
        return 1
    } else {
        return 0
    }
}

// has become unnecessarily brief
/*
// Returns the indices of the minimizer's kmer range in RepeatGenome.Kmers
func (rg *RepeatGenome) getMinIndices(minInt MinInt) (uint64, uint64) {
    startInd := rg.MinOffsets[minInt]
    endInd := startInd + rg.MinCounts[minInt]
    return startInd, endInd
}
*/

// Takes as an argument the index of a minimizer in RepeatGenome.SortedMins - NOT the minimizer's actual value.
// Returns a slice of the kmers associated with that minimizer.
// Intended for non-performance-critical sections, like file writing.
func (rg *RepeatGenome) getMinsKmers(minInt MinInt) Kmers {
    startInd := rg.MinOffsets[minInt]
    endInd := startInd + int64(rg.MinCounts[minInt])
    return rg.Kmers[startInd:endInd]
}

func ceilDiv_U64(a, b uint64) uint64 {
    if a == 0 {
        return 0
    } else {
        return ((a - 1) / b) + 1
    }
}

func ceilDiv_64(a, b int64) int64 {
    return ((a - 1) / b) + 1
}

func ceilDiv(a, b int) int {
    return ((a - 1) / b) + 1
}

// Assumes that rg.MinCounts and rg.SortedMins are populated.
func (rg *RepeatGenome) populateMinOffsets() {
    if rg.MinOffsets != nil && len(rg.MinOffsets) > 0 {
        fmt.Println("!!! WARNING !!! RepeatGenome.populateMinOffsets overwriting previous non-empty value")
    }

    rg.MinOffsets = make([]int64, minSliceSize, minSliceSize)
    // initialize to -1 to make explicit which indices are not mins
    for i := 0; i < len(rg.MinOffsets); i++ {
        rg.MinOffsets[i] = -1
    }

    if len(rg.SortedMins) > 0 {
        rg.MinOffsets[rg.SortedMins[0]] = 0
        lastMin := rg.SortedMins[0]
        for i := 1; i < len(rg.SortedMins); i++ {
            thisMin := rg.SortedMins[i]
            rg.MinOffsets[thisMin] = rg.MinOffsets[lastMin] + int64(rg.MinCounts[lastMin])
            lastMin = thisMin
        }
    }
}

func zero(byteSlice []byte) {
    for i := range byteSlice {
        byteSlice[i] = 0
    }
}

func (repeat Repeat) IsRoot() bool {
    return repeat.ID == 0
}
