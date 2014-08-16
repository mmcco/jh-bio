package repeatgenome

/*
   Functions used in read classification.
*/

import (
    "bytes"
    "fmt"
    "os"
    "sync"
)

/*
   The type sent back from read-classifying goroutines of RepeatGenome.ClassifyReads()
*/
type ReadResponse struct {
    Seq       TextSeq
    ClassNode *ClassNode
}

/*
   Classifies each read in a slice of reads, stored as type TextSeq.
   The read and its classification are returned through responseChan.
   In the future, the reads may be of type Seq.
   However, this currently seems to be the fastest way of doing things.
   This version uses the first recognized kmer for classification - the Kraken-Q technique.
*/
func (rg *RepeatGenome) QuickClassifyReads(readTextSeqs []TextSeq, responseChan chan ReadResponse) {
ReadLoop:
    for _, read := range readTextSeqs {
        // signed integers prevent underflow in down-counting for loops
        k_ := int64(k)
        numKmers := int64(len(read)) - k_ + 1

        var i int64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            for j := k_ + i - 1; j >= i; j-- {
                if read[j] == byte('n') {
                    i += j - i
                    continue KmerLoop
                }
            }

            kmerBytes := read[i : i+k_]
            kmerInt := kmerBytes.kmerInt().canonicalRepr()
            kmerLCA := rg.getKmerLCA(kmerInt)

            if kmerLCA != nil {
                // only use the first matched kmer
                responseChan <- ReadResponse{read, kmerLCA}
                continue ReadLoop
            }
        }
        responseChan <- ReadResponse{read, nil}
    }
    close(responseChan)
}

/*
   Classifies each read in a slice of reads, stored as type TextSeq.
   The read and its classification are returned through responseChan.
   In the future, the reads may be of type Seq.
   However, this currently seems to be the fastest way of doing things.
   This version returns the LCA of all recognized kmers' classifications.
*/
func (rg *RepeatGenome) LCA_ClassifyReads(readTextSeqs []TextSeq, responseChan chan ReadResponse) {
    for _, read := range readTextSeqs {
        // signed integers prevent underflow in down-counting for loops
        k_ := int64(k)
        numKmers := int64(len(read)) - k_ + 1

        var lca *ClassNode = nil
        var i int64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            for j := k_ + i - 1; j >= i; j-- {
                if read[j] == byte('n') {
                    i += j - i
                    continue KmerLoop
                }
            }

            kmerBytes := read[i : i+k_]
            kmerInt := kmerBytes.kmerInt().canonicalRepr()
            kmerLCA := rg.getKmerLCA(kmerInt)

            if kmerLCA != nil {
               if lca == nil {
                   lca = kmerLCA
               } else {
                   lca = rg.ClassTree.getLCA(lca, kmerLCA)
                }
            }
        }
        responseChan <- ReadResponse{read, lca}
    }
    close(responseChan)
}

/*
   Dispatches as many read-classifying goroutines as there are CPUs, giving each a subslice of the slice of reads provided.
   Each read-classifying goroutine is given a unique response chan.
   These are then merged into a single response chan, which is the return value.
   The useLCA parameter determines whether to use Quick or LCA read classification logic.
*/
func (rg *RepeatGenome) GetClassChan(reads []TextSeq, useLCA bool) chan ReadResponse {
    responseChans := make([]chan ReadResponse, 0, numCPU)

    // dispatch one read-classifying goroutine per CPU
    for i := 0; i < numCPU; i++ {
        responseChans = append(responseChans, make(chan ReadResponse, 50))
        startInd := i * (len(reads) / numCPU)
        endInd := ((i + 1) * len(reads)) / numCPU
        if useLCA {
            go rg.LCA_ClassifyReads(reads[startInd:endInd], responseChans[i])
        } else {
            go rg.QuickClassifyReads(reads[startInd:endInd], responseChans[i])
        }
    }

    // the rest of this function is spent merging the response chans
    var wg sync.WaitGroup
    wg.Add(len(responseChans))
    master := make(chan ReadResponse)

    for _, respChan := range responseChans {
        // need respChan argument to prevent all goroutines from using same instance of respChan - language quirk
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

/*
   A rather hairy function that classifies all reads in ./<genome-name>-reads/*.proc if any exist.
   .proc files are our own creation for ease of parsing and testing.
   They contain one lowercase read sequence per line, and nothing else.
   We have a script that will convert FASTQ files to .proc files: github.com/plsql/bioinformatics/blob/master/scripts/format-FASTA-reads.py
   This is generally really easy to do.
   However, we will used a FASTQ reader when we get past the initial testing phase.
   This could be done concurrently, considering how many disk accesses there are.
*/
func (rg *RepeatGenome) GetProcReads() (error, []TextSeq) {
    workingDirName, err := os.Getwd()
    if err != nil {
        return err, nil
    }
    readsDirName := workingDirName + "/" + rg.Name + "-reads"
    currDir, err := os.Open(readsDirName)
    if err != nil {
        return err, nil
    }
    // get os.FileInfo struct values of all files in the reads dir
    fileinfos, err := currDir.Readdir(-1)
    if err != nil {
        return err, nil
    }
    // populate a list of os.FileInfo values representing .proc files
    var procFiles []os.FileInfo
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 5 && fileinfo.Name()[len(fileinfo.Name())-5:] == ".proc" {
            procFiles = append(procFiles, fileinfo)
        }
    }
    // populate a list of the reads contained in each of these files
    var reads []TextSeq
    for _, fileinfo := range procFiles {
        err, readsBytes := fileLines(readsDirName + "/" + fileinfo.Name())
        if err != nil {
            return err, nil
        }
        for _, lineBytes := range readsBytes {
            reads = append(reads, TextSeq(lineBytes))
        }
    }

    return nil, reads
}

/*
   
*/
func (rg *RepeatGenome) GetSAMReads() (error, []TextSeq) {
    workingDirName, err := os.Getwd()
    if err != nil {
        return err, nil
    }
    readsDirName := workingDirName + "/" + rg.Name + "-reads"
    currDir, err := os.Open(readsDirName)
    if err != nil {
        return err, nil
    }
    // get os.FileInfo struct values of all files in the reads dir
    fileinfos, err := currDir.Readdir(-1)
    if err != nil {
        return err, nil
    }
    // populate a list of os.FileInfo values representing .sam files
    var samFiles []os.FileInfo
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 4 && fileinfo.Name()[len(fileinfo.Name())-4:] == ".sam" {
            samFiles = append(samFiles, fileinfo)
        }
    }
    // populate a list of the reads contained in each of these files
    var reads []TextSeq
    for _, fileinfo := range samFiles {
        err, readsBytes := fileLines(readsDirName + "/" + fileinfo.Name())
        if err != nil {
            return err, nil
        }
        // drop header
        if len(readsBytes) < 3 {
            return ParseError{"RepeatGenome.GetSAMReads", fileinfo.Name(), fmt.Errorf("too few lines in SAM file - header missing")}, nil
        } else {
            readsBytes = readsBytes[3:]
        }

        for _, lineBytes := range readsBytes {
            fields := bytes.Fields(lineBytes)
            if len(fields) != 12 {
                return ParseError{"RepeatGenome.GetSAMReads", fileinfo.Name(), fmt.Errorf("too few fields in SAM file line")}, nil
            } else {
                reads = append(reads, TextSeq(bytes.ToLower(fields[9])))
            }
        }
    }

    return nil, reads
}
