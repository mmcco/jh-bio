package repeatgenome

/*
   A collection of trivial functions used in other source files of the repeatgenome package.
*/

import (
    "bytes"
    "fmt"
    "io/ioutil"
    //"reflect"
    "strconv"
    "strings"
    "sync"
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

// Courtesy of https://github.com/dustin/go-humanize
// Returns a string representing the int, with commas for readability.
func comma(v uint64) string {
    sign := ""
    if v < 0 {
        sign = "-"
        v = 0 - v
    }

    parts := []string{"", "", "", "", "", "", "", ""}
    j := len(parts) - 1

    for v > 999 {
        parts[j] = strconv.FormatUint(v%1000, 10)
        switch len(parts[j]) {
        case 2:
            parts[j] = "0" + parts[j]
        case 1:
            parts[j] = "00" + parts[j]
        }
        v = v / 1000
        j--
    }
    parts[j] = strconv.Itoa(int(v))
    return sign + strings.Join(parts[j:len(parts)], ",")
}

func mergeThreadResp(cs [](chan ThreadResponse)) <-chan ThreadResponse {
    var wg sync.WaitGroup
    //elemType := reflect.TypeOf(cs).Elem()
    //chanType := reflect.ChanOf(RecvDir, elemType)
    out := make(chan ThreadResponse)

    // Start an output goroutine for each input channel in cs.  output
    // copies values from c to out until c is closed, then calls wg.Done.
    wg.Add(len(cs))
    for _, c := range cs {
        go func(c chan ThreadResponse) {
            for n := range c {
                    out <- n
            }
            wg.Done()
        }(c)
    }

    // Start a goroutine to close out once all the output goroutines are
    // done.  This must start after the wg.Add call.
    go func() {
        wg.Wait()
        close(out)
    }()
    return out
}

// returns the number of lines and a slice of the lines
func lines(byteSlice []byte) [][]byte {
    var lines [][]byte = bytes.Split(byteSlice, []byte{'\n'})
    // drop the trailing newlines
    newline := []byte("\n")
    for lastLine := lines[len(lines)-1]; len(lines) > 0 && (len(lastLine) == 0 || bytes.Equal(lastLine, newline)); lastLine = lines[len(lines)-1] {
        lines = lines[:len(lines)-1]
    }
    return lines
}

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

func minMinInt(a, b MinInt) MinInt {
    if a < b {
        return a
    } else {
        return b
    }
}

func minKmerInt(a, b KmerInt) KmerInt {
    if a < b {
        return a
    } else {
        return b
    }
}

// needed for sort.Interface
func (kmers PKmers) Len() int {
    return len(kmers)
}

func (kmers PKmers) Swap(i, j int) {
    kmers[i], kmers[j] = kmers[j], kmers[i]
}

func (kmers PKmers) Less(i, j int) bool {
    iVal := *(*uint64)(unsafe.Pointer(&kmers[i][0]))
    jVal := *(*uint64)(unsafe.Pointer(&kmers[j][0]))
    return iVal < jVal
}

func (kmers Kmers) Len() int {
    return len(kmers)
}

func (kmers Kmers) Swap(i, j int) {
    kmers[i], kmers[j] = kmers[j], kmers[i]
}

func (kmers Kmers) Less(i, j int) bool {
    iVal := *(*uint64)(unsafe.Pointer(&kmers[i][0]))
    jVal := *(*uint64)(unsafe.Pointer(&kmers[j][0]))
    return iVal < jVal
}

func (minInts MinInts) Len() int {
    return len(minInts)
}

func (minInts MinInts) Swap(i, j int) {
    minInts[i], minInts[j] = minInts[j], minInts[i]
}

func (minInts MinInts) Less(i, j int) bool {
    return minInts[i] < minInts[j]
}

func (minPairs MinPairs) Len() int {
    return len(minPairs)
}

func (minPairs MinPairs) Swap(i, j int) {
    minPairs[i], minPairs[j] = minPairs[j], minPairs[i]
}

func (minPairs MinPairs) Less(i, j int) bool {
    iMin := *(*MinInt)(unsafe.Pointer(&minPairs[i][8]))
    jMin := *(*MinInt)(unsafe.Pointer(&minPairs[j][8]))
    if iMin < jMin {
        return true
    } else if iMin > jMin {
        return false
    } else {
        return *(*KmerInt)(unsafe.Pointer(&minPairs[i])) < *(*KmerInt)(unsafe.Pointer(&minPairs[j]))
    }
}

func (fullKmers FullKmers) Len() int {
    return len(fullKmers)
}

func (fullKmers FullKmers) Swap(i, j int) {
    fullKmers[i], fullKmers[j] = fullKmers[j], fullKmers[i]
}

// could also just loop through bytes and compare, if the MinInt came first
func (fullKmers FullKmers) Less(i, j int) bool {
    iMin := *(*MinInt)(unsafe.Pointer(&fullKmers[i][8]))
    jMin := *(*MinInt)(unsafe.Pointer(&fullKmers[j][8]))
    if iMin < jMin {
        return true
    } else if iMin > jMin {
        return false
    } else {
        return *(*KmerInt)(unsafe.Pointer(&fullKmers[i])) < *(*KmerInt)(unsafe.Pointer(&fullKmers[j]))
    }
}

func boolToInt(a bool) int {
    if a {
        return 1
    } else {
        return 0
    }
}

// Takes as an argument the index of a minimizer in RepeatGenome.SortedMins - NOT the minimizer's actual value.
// Returns the indices of the minimizer's kmer range in RepeatGenome.Kmers
func (rg *RepeatGenome) getMinIndices(minIndex uint64) (uint64, uint64) {
    startInd := rg.OffsetsToMin[minIndex]
    var endInd uint64
    if minIndex == uint64(len(rg.SortedMins)) - 1 {
        endInd = uint64(len(rg.Kmers))
    } else {
        endInd = rg.OffsetsToMin[minIndex+1]
    }
    return startInd, endInd
}

// Takes as an argument the index of a minimizer in RepeatGenome.SortedMins - NOT the minimizer's actual value.
// Returns a slice of the kmers associated with that minimizer.
// Intended for non-performance-critical sections, like file writing.
func (rg *RepeatGenome) getMinsKmers(minIndex uint64) []Kmer {
    startInd, endInd := rg.getMinIndices(minIndex)
    return rg.Kmers[startInd : endInd]
}

func (rg *RepeatGenome) getMinsFullKmers(minIndex uint64) FullKmers {
    startInd, endInd := rg.getMinIndices(minIndex)
    return rg.FullKmers[startInd : endInd]
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
