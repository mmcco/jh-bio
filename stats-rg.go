package repeatgenome

import (
    "bytes"
    "unsafe"
)

func (repeatGenome *RepeatGenome) Size() uint64 {
    var numBases uint64 = 0
    for _, seqs := range repeatGenome.chroms {
        for _, seq := range seqs {
            numBases += uint64(len(seq))
        }
    }
    return numBases
}

func (rg *RepeatGenome) KmersGBSize() float64 {
    return (float64(len(rg.Kmers))/1000000000) * float64(unsafe.Sizeof(Kmer{}))
}

func (rg *RepeatGenome) PercentRepeats() float64 {
    var totalBases uint64 = rg.Size()

    var repeatBases uint64 = 0
    for _, match := range rg.Matches {
        repeatBases += match.SeqEnd - match.SeqStart
    }

    return 100 * (float64(repeatBases) / float64(totalBases))
}

func (repeat *Repeat) RepeatSize() uint64 {
    var repeatSize uint64 = 0

    for _, match := range repeat.Instances {
        repeatSize += match.SeqEnd - match.SeqStart
    }

    return repeatSize
}

func (classNode *ClassNode) Size() uint64 {
    if classNode == nil {
        return 0
    }
    
    var classNodeSize uint64 = 0

    if classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            classNodeSize += match.SeqEnd - match.SeqStart
        }
    }

    for _, child := range classNode.Children {
        classNodeSize += child.Size()
    }

    return classNodeSize
}

func (readResp ReadResponse) HangingSize() uint64 {
    classNode := readResp.ClassNode
    if classNode == nil {
        return 0
    }
    
    var classNodeSize uint64 = 0

    if classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            classNodeSize += match.SeqEnd - match.SeqStart
            classNodeSize += 2 * uint64(len(readResp.Seq) - int(k))
        }
    }

    for _, child := range classNode.Children {
        classNodeSize += ReadResponse{readResp.Seq, child}.HangingSize()
    }

    return classNodeSize
}

// returns the average percent of the genome that a classified read could exist in, in regard to the supplied list of classified reads
// uses a cumulative average to prevent overflow
func (rg *RepeatGenome) AvgPossPercentGenome(resps []ReadResponse, strict bool) float64 {
    classNodeSizes := make(map[*ClassNode]float64, len(rg.ClassTree.ClassNodes))
    if strict {
        for _, classNode := range rg.ClassTree.NodesByID {
            classNodeSizes[classNode] = float64(classNode.Size())
        }
    } else {
        for _, resp := range resps {
            if _, exists := classNodeSizes[resp.ClassNode]; !exists {
                classNodeSizes[resp.ClassNode] = float64(resp.HangingSize())
            }
        }
    }

    var classesProcessed, avgClassSize float64 = 0, 0
    for _, resp := range resps {
        avgClassSize += (classNodeSizes[resp.ClassNode] - avgClassSize) / (classesProcessed + 1)
        classesProcessed++
    }

    return 100 * (avgClassSize / float64(rg.Size()))
}

// written for the PercentTrueClassification() below
// determines whether a read overlaps any repeat instances in the given ClassNode's subtree
func (rg *RepeatGenome) recNodeSearch(classNode *ClassNode, readSAM ReadSAM, strict bool) bool {
    if classNode != nil && classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            // must compute where the read ends
            endInd := readSAM.StartInd + uint64(len(readSAM.Seq))
            overlap := readSAM.SeqName == match.SeqName && readSAM.StartInd < match.SeqEnd && endInd > match.SeqStart
            if overlap && !strict {
                return true
            } else if overlap && strict {
                // below logic is for checking for at least rg.K overlap
                var overlap uint64 = endInd - match.SeqStart
                if readSAM.StartInd > match.SeqStart {
                    overlap -= readSAM.StartInd - match.SeqStart
                }
                if overlap >= uint64(k) {
                    return true
                }
            }
        }
    }
    if classNode != nil && classNode.Children != nil {
        for _, child := range classNode.Children {
            if rg.recNodeSearch(child, readSAM, strict) {
                return true
            }
        }
    }
    return false
}

/*
func TestNodeSearch(classNode *ClassNode, readSAM ReadSAM) bool {
    if classNode == nil || (classNode.Name != "root" && classNode.Name != "Satellite" && classNode.Name != "Satellite/HETRP_DM") {
        return false
    }
    fmt.Println("testing", classNode.Name)
    if classNode != nil && classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            // must compute where the read ends
            endInd := readSAM.StartInd + uint64(len(readSAM.Seq))
            fmt.Printf("testing match %s[%d:%d] against %s[%d:%d]\n", match.SeqName, match.SeqStart, match.SeqEnd, readSAM.SeqName, readSAM.StartInd, endInd)
            if readSAM.SeqName == match.SeqName && readSAM.StartInd < match.SeqEnd && endInd > match.SeqStart {
                fmt.Println("true")
                return true
                // below logic is for checking for at least rg.K overlap
                /*
                var overlap uint64 := readSAM.SeqEnd - match.SeqStart
                if readSAM.SeqStart > match.SeqStart {
                    overlap -= readSAM.SeqStart - match.SeqStart
                }
                if overlap >= uint64(rg.K) {
                    return true
                }
            }
        }
    } else {
        fmt.Println("no classNode or no classNode.Repeat")
    }
    if classNode != nil && classNode.Children != nil {
        for _, child := range classNode.Children {
            if TestNodeSearch(child, readSAM) {
                fmt.Println("true")
                return true
            } else {
                fmt.Println("(child false)")
            }
        }
    } else {
        fmt.Println("no classNode or no classNode.Children")
    }
    fmt.Println("false")
    return false
}
*/

// we currently use the simple metric that the read and one of the repeat's instances overlap at all
func (rg *RepeatGenome) PercentTrueClassifications(responses []ReadSAMResponse, strict bool) float64 {
    var classifications, correctClassifications uint64 = 0, 0

    for _, resp := range responses {
        if resp.ClassNode != nil {
            classifications++
        }
        if rg.recNodeSearch(resp.ClassNode, resp.ReadSAM, strict) {
            correctClassifications++
        }
    }
    return 100 * (float64(correctClassifications) / float64(classifications))
}

func (rg *RepeatGenome) numRawKmers() uint64 {
    var k_ = int(k)
    var numRawKmers uint64 = 0

    splitOnN := func(c rune) bool { return c == 'n' }

    for i := range rg.Matches {
        match := &rg.Matches[i]
        seq := rg.chroms[match.SeqName][match.SeqName][match.SeqStart:match.SeqEnd]
        seqs := bytes.FieldsFunc([]byte(seq), splitOnN)
        for j := range seqs {
            if len(seqs[j]) >= k_ {
                numRawKmers += uint64(len(seqs[j]) - k_ + 1)
            }
        }
    }
    return numRawKmers
}
