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
    Filename string
    SubError error
}

type IOError struct {
    FuncName string
    SubError error
}

func (parseError ParseError) Error() string {
    if len(parseError.Filename) == 0 {
        return fmt.Sprintf("%s: parsing error: %s", parseError.FuncName, parseError.SubError.Error())
    } else {
        return fmt.Sprintf("%s: parsing error in file %s: %s", parseError.FuncName, parseError.Filename, parseError.SubError.Error())
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

// needed for sort.Interface
type Uint64Slice []uint64

func (uint64s Uint64Slice) Len() int {
    return len(uint64s)
}

func (uint64s Uint64Slice) Swap(i, j int) {
    uint64s[i], uint64s[j] = uint64s[j], uint64s[i]
}

func (uint64s Uint64Slice) Less(i, j int) bool {
    return uint64s[i] < uint64s[j]
}

func boolToInt(a bool) int {
    if a {
        return 1
    } else {
        return 0
    }
}
