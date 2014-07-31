package repeatgenome

import (
    "os"
    "sync"
    "unsafe"
)

type ReadResponse struct {
    Seq       TextSeq
    ClassNode *ClassNode
}

// This function is a mess and needs cleaning up
func (rg *RepeatGenome) ClassifyReads(readTextSeqs []TextSeq, responseChan chan ReadResponse) {
ReadLoop:
    for _, read := range readTextSeqs {
        // we use sign int64s in this triple-nested loop because the innermost one counts down and would otherwise overflow
        // this isn't a significant limitation because reads are never big enough to overflow one
        k_ := int64(k)
        numKmers := int64(len(read)) - k_ + 1

        var i int64
    KmerLoop:
        for i = 0; i < numKmers; i++ {

            // if this is ever revived, all the loop ints must be made signed again
            for j := k_ + i - 1; j >= i ; j-- {
                if read[j] == byte('n') {
                    i += j - i
                    continue KmerLoop
                }
            }

            kmerBytes := read[i : i+k_]
            kmerInt := kmerBytes.kmerInt().canonicalRepr()
            kmer := rg.getKmer(kmerInt)

            if kmer != nil {
                lcaID := *(*uint16)(unsafe.Pointer(&kmer[8]))
                responseChan <- ReadResponse{read, rg.ClassTree.NodesByID[lcaID]}
                // only use the first matched kmer
                continue ReadLoop
            }
        }
        responseChan <- ReadResponse{read, nil}
    }
    close(responseChan)
}

func (rg *RepeatGenome) GetReadClassChan(reads []TextSeq) chan ReadResponse {
    responseChans := make([]chan ReadResponse, 0, numCPU)

    numReads := uint64(len(reads))
    numCPU_ := uint64(numCPU)
    var i uint64
    for i = 0; i < uint64(numCPU_); i++ {
        responseChans = append(responseChans, make(chan ReadResponse, 50))
        startInd := i * (numReads / numCPU_)
        endInd := ((i + 1) * numReads) / numCPU_
        go rg.ClassifyReads(reads[startInd : endInd], responseChans[i])
    }

    var wg sync.WaitGroup
    wg.Add(len(responseChans))
    master := make(chan ReadResponse)

    for _, respChan := range responseChans {
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

func (rg *RepeatGenome) ProcessReads() (error, chan ReadResponse) {
    workingDirName, err := os.Getwd()
    if err != nil {
        return err, nil
    }
    readsDirName := workingDirName + "/" + rg.Name + "-reads"
    currDir, err := os.Open(readsDirName)
    if err != nil {
        return err, nil
    }
    fileinfos, err := currDir.Readdir(-1)
    if err != nil {
        return err, nil
    }
    processedFiles := []os.FileInfo{}
    for _, fileinfo := range fileinfos {
        if len(fileinfo.Name()) > 5 && fileinfo.Name()[len(fileinfo.Name())-5 : ] == ".proc" {
            processedFiles = append(processedFiles, fileinfo)
        }
    }
    var reads []TextSeq
    for _, fileinfo := range processedFiles {
        _, theseReadsBytes := fileLines(readsDirName + "/" + fileinfo.Name())
        for _, bytesLine := range theseReadsBytes {
            reads = append(reads, TextSeq(bytesLine))
        }
    }

    return nil, rg.GetReadClassChan(reads)
}
