package main

import (
    "bytes"
    "flag"
    "fmt"
    "github.com/plsql/jh-bio/repeatgenome"
    "io/ioutil"
    "os"
    "runtime/pprof"
    "strconv"
    "strings"
    "time"
)

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
    return sign + strings.Join(parts[j:], ",")
}

func main() {

    if len(os.Args) < 2 {
        fmt.Println("arg error - usage: ./minimize <flags> <reference genome dir>")
        os.Exit(1)
    }
    genomeName := os.Args[len(os.Args)-1]

    k_arg := flag.Uint("k", 31, "kmer length")
    m_arg := flag.Uint("m", 15, "minimizer length")
    forceGen := flag.Bool("force_gen", false, "force Kraken database generation, regardless of whether it already exists in stored form")
    writeStats := flag.Bool("write_stats", false, "write various tab-delimited and JSON files representing peripheral Kraken and repeat data")
    dontWriteLib := flag.Bool("no_write_lib", false, "don't write the Kraken library to file")
    verifyClass := flag.Bool("verify_class", false, "run classification a second time, with SAM-formatted reads, to find percent correct classification")
    debug := flag.Bool("debug", false, "run and print debugging tests")
    cpuProfile := flag.Bool("cpuprof", false, "write cpu profile to file <genomeName>.cpuprof")
    memProfile := flag.Bool("memprof", false, "write memory profile to <genomeName>.memprof")
    lcaClassify := flag.Bool("lca_classify", false, "use the LCA of all recognized kmers' classes as a read's classification")
    //useRoot := flag.Bool("use_root", false, "include kmers with root as their LCA in the Kraken DB, and return root read classifications rather than nil")
    flag.Parse()

    if *cpuProfile {
        fmt.Println("CPU profiler enabled")
        fmt.Println()
        os.Mkdir("profiles", os.ModeDir)
        f, err := os.Create("profiles/" + genomeName + ".cpuprof")
        if err != nil {
            panic(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    }

    if *memProfile {
        fmt.Println("memory profiler enabled")
        fmt.Println()
    }

    if *debug {
        fmt.Println("debug tests enabled")
        fmt.Println()
    }

    if *dontWriteLib {
        fmt.Println("Writing Kraken library to file disabled")
        fmt.Println()
    }

    var k, m uint8
    if *k_arg > 255 || *m_arg > 255 {
        panic("k and m must be >= 255")
    } else {
        k = uint8(*k_arg)
        m = uint8(*m_arg)
    }

    err, rg := repeatgenome.New(repeatgenome.Config{
        Name:       genomeName,
        K:          k,
        M:          m,
        Debug:      *debug,
        CPUProfile: *cpuProfile,
        MemProfile: *memProfile,
        WriteLib:   !*dontWriteLib,
        ForceGen:   *forceGen,
        WriteStats: *writeStats,
    })

    if err != nil {
        fmt.Println("./run-rg: RepeatGenome generation failed:")
        fmt.Println(err)
        os.Exit(1)
    }

    fmt.Println(comma(uint64(len(rg.Repeats))), "repeat types")
    fmt.Println(comma(uint64(len(rg.ClassTree.ClassNodes))), "class nodes")
    fmt.Println(comma(uint64(len(rg.Matches))), "matches")
    fmt.Println()

    if *lcaClassify {
        fmt.Println("using LCA logic to classify reads")
        fmt.Println()
    }

    err, reads := rg.GetSAMReads()
    if err != nil {
        fmt.Println(err)
        os.Exit(1)
    }
    respChan := rg.GetClassChan(reads, *lcaClassify)
    startTime := time.Now()
    for range respChan {
    }
    netTime := time.Since(startTime)

    var numReads, numClassifiedReads, rootReads uint64 = 0, 0, 0
    var responses []repeatgenome.ReadResponse
    err, reads = rg.GetSAMReads()
    if err != nil {
        fmt.Println(err)
        os.Exit(1)
    }
    respChan = rg.GetClassChan(reads, *lcaClassify)
    for response := range respChan {
        responses = append(responses, response)
        _, classNode := response.Seq, response.ClassNode
        numReads++
        if classNode != nil {
            numClassifiedReads++
            if classNode == rg.ClassTree.Root {
                rootReads++
            }
        }
    }

    /*  deprecated by changed functions - should eventually find another way of making this check
        if numReads != uint64(len(readsBytes)) {
            panic("not all reads, or too many reads, returned from RepeatGenome.GetReadClassChan()")
        }
    */

    /*
       if rg.Flags.Debug {
           classCount := make(map[*repeatgenome.ClassNode]uint64)
           err, respChan = rg.ProcessReads()
           for response := range respChan {
               if response.ClassNode != nil {
                   classCount[response.ClassNode]++
               }
           }

           for classNode, count := range classCount {
               fmt.Printf("%s\t%d\n", classNode.Name, count)
           }
       }
    */

    var nonRootResps []repeatgenome.ReadResponse
    for _, resp := range responses {
        if resp.ClassNode != nil && resp.ClassNode.Name != "root" {
            nonRootResps = append(nonRootResps, resp)
        }
    }

    fmt.Printf("RepeatGenome.Kmers comprises %.2f GB\n", rg.KmersGBSize())
    fmt.Println()

    fmt.Printf("%.2f million reads processed per minute\n", (float64(numReads)/1000000)/netTime.Minutes())
    fmt.Printf("%.2f%% of the genome consists of repeat sequences\n", rg.PercentRepeats())
    fmt.Printf("%.2f%% of reads were classified with a repeat sequence (%s of %s)\n", 100*(float64(numClassifiedReads)/float64(numReads)), comma(numClassifiedReads), comma(numReads))
    fmt.Printf("%.2f%% of classified reads were classified at the class tree root (%s reads)\n", 100*(float64(rootReads)/float64(numReads)), comma(rootReads))
    fmt.Printf("on average, a classification restricted a read's possible location to %.2f%% of the genome\n", rg.AvgPossPercentGenome(responses, true))
    fmt.Printf("on average, a non-root classification restricted a read's possible location to %.2f%% of the genome\n", rg.AvgPossPercentGenome(nonRootResps, true))
    fmt.Printf("on average, a classification restricted a read's possible location to %.2f%% of the genome (non-strict)\n", rg.AvgPossPercentGenome(responses, false))
    fmt.Printf("on average, a non-root classification restricted a read's possible location to %.2f%% of the genome (non-strict)\n", rg.AvgPossPercentGenome(nonRootResps, false))
    fmt.Println()

    if *verifyClass {
        fmt.Println("...using SAM-formatted reads to check classification correctness...")

        workingDirName, err := os.Getwd()
        if err != nil {
            fmt.Println(err)
            os.Exit(1)
        }
        readsDirName := workingDirName + "/" + rg.Name + "-reads"

        err, readSAMs := repeatgenome.GetReadSAMs(readsDirName)
        if err != nil {
            fmt.Println(err)
            os.Exit(1)
        }

        seqToClass := make(map[string]*repeatgenome.ClassNode, len(responses))
        for _, response := range responses {
            seqToClass[string(response.Seq)] = response.ClassNode
        }

        readSAMResps := []repeatgenome.ReadSAMResponse{}
        for _, readSAM := range readSAMs {
            if _, exists := seqToClass[string(readSAM.Seq)]; !exists {
                fmt.Println(string(readSAM.Seq), "present in ReadSAM but not Read")
                os.Exit(1)
            }
            resp := repeatgenome.ReadSAMResponse{
                ReadSAM:   readSAM,
                ClassNode: seqToClass[string(readSAM.Seq)],
            }
            readSAMResps = append(readSAMResps, resp)
        }

        fmt.Println("parsed", len(readSAMResps), "ReadSAMResponses")

        fmt.Printf("%.2f%% of classified reads overlapped an instance of their assigned repeat class\n", rg.PercentTrueClassifications(readSAMResps, false))
        fmt.Printf("%.2f%% of classified reads overlapped an instance of their assigned repeat class (strict)\n\n", rg.PercentTrueClassifications(readSAMResps, true))
    }
}
