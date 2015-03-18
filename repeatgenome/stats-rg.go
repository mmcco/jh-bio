/*
   Functions used for (mostly simple) statistical analysis of the repeat sequence and Kraken data.
*/

package repeatgenome

import (
    "bytes"
    "unsafe"
)

/*
   Returns the total number of bases in a RepeatGenome's reference chromosomes.
*/
func (repeatGenome *RepeatGenome) Size() uint64 {
    var numBases uint64 = 0
    for _, seqs := range repeatGenome.chroms {
        for _, seq := range seqs {
            numBases += uint64(len(seq))
        }
    }
    return numBases
}

/*
   Returns the size in gigabytes of the supplied RepeatGenome's Kmers field.
*/
func (rg *RepeatGenome) KmersGBSize() float64 {
    return (float64(len(rg.Kmers)) / 1000000000) * float64(unsafe.Sizeof(Kmer{}))
}

/*
   Returns the percent of a RepeatGenome's reference bases that are contained in a repeat instance.
   It makes the assumption that no base is contained in more than one repeat instance.
*/
func (rg *RepeatGenome) PercentRepeats() float64 {
    var repeatBases uint64 = 0
    for _, match := range rg.Matches {
        repeatBases += match.SeqEnd - match.SeqStart
    }

    return 100 * (float64(repeatBases) / float64(rg.Size()))
}

/*
   Returns the sum of the sizes of all of a repeat sequence type's instances.
*/
func (repeat *Repeat) Size() uint64 {
    var repeatSize uint64 = 0

    for _, match := range repeat.Instances {
        repeatSize += match.SeqEnd - match.SeqStart
    }

    return repeatSize
}

/*
   Returns the sum of the sizes of all repeat instances in the supplied ClassNode's subtree.
*/
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

/*
   Returns the number of base pairs from which the supplied read could have originated, assuming that its classification was correct.
   This is done in terms of Kraken-Q logic, meaning that there is at least one kmer shared between the repeat reference and the read.
   Therefore, the read must overlap a repeat reference from the classified subtree by at least k bases.
   This function is used to calculate the probability of correct classification assuming random selection, and the amount to which a classification narrows a read's potential origin.
*/
func (readResp ReadResponse) HangingSize() uint64 {
    classNode := readResp.ClassNode
    if classNode == nil {
        return 0
    }

    var classNodeSize uint64 = 0

    if classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            classNodeSize += match.SeqEnd - match.SeqStart
            classNodeSize += 2 * uint64(len(readResp.Seq)-int(k))
        }
    }

    for _, child := range classNode.Children {
        classNodeSize += ReadResponse{readResp.Seq, child}.HangingSize()
    }

    return classNodeSize
}

/*
   Returns the average percent of the genome a read from the given set could have originated from, assuming their classification was correct.
   This is used to estimate how much the classification assisted us in locating reads' origins.
   The more specific and helpful the classifications are, the lower the percentage will be.
   Uses a cumulative average to prevent overflow.
*/
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

/*
   A helper function for the PercentTrueClassification() below
   It recursively determines whether a read originated from a reference repeat instance in the subtree indicated by the supplied ClassNode.
*/
func (rg *RepeatGenome) recNodeSearch(classNode *ClassNode, readSAM ReadSAM, strict bool) bool {
    if classNode != nil && classNode.Repeat != nil {
        for _, match := range classNode.Repeat.Instances {
            if match.SeqName != readSAM.SeqName {
                continue
            }
            // must compute where the read ends
            endInd := readSAM.StartInd + uint64(len(readSAM.TextSeq))
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


func (rg *RepeatGenome) RepeatIsCorrect(readSAMRepeat ReadSAMRepeat, strict bool) bool {
    // awkward unpacking - maybe use separate args?
    readSAM, repeat := readSAMRepeat.ReadSAM, readSAMRepeat.Repeat
    read, seqName, startInd := readSAM.TextSeq, readSAM.SeqName, readSAM.StartInd
    if repeat == nil {
        // We will for now use a panic rather than an error.
        // This is for speed and simplicity, and because the caller
        // logically should ensure than the repeat is non-nil.
        panic("RepeatGenome.RepeatIsCorrect(): readSAMRepeat.Repeat is nil")
    }
    for _, match := range repeat.Instances {
        if match.SeqName != seqName {
            continue
        }
        // must compute where the read ends - it isn't stored
        endInd := readSAM.StartInd + uint64(len(read))
        overlap := seqName == match.SeqName && startInd < match.SeqEnd && endInd > match.SeqStart
        if overlap && !strict {
            return true
        } else if overlap && strict {
            // below logic is for checking for at least rg.K overlap
            var overlap uint64 = endInd - match.SeqStart
            if startInd > match.SeqStart {
                overlap -= startInd - match.SeqStart
            }
            if overlap >= uint64(k) {
                return true
            }
        }
    }
    return false
}
/*
func TestNodeSearch(classNode *ClassNode, readSAM ReadSAM) bool {
    if classNode == nil {
        return false
    }

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

/*
   Determines whether a read overlaps any repeat instances in the given ClassNode's subtree.
   If the argument strict is true, the read must be entirely contained in a reference repeat instance (classic Kraken logic).
   Otherwise, the read must overlap a reference repeat instance by at least k bases.
*/
func (rg *RepeatGenome) PercentTrueClassifications(responses []ReadSAMResponse, useStrict bool) float64 {
    var classifications, correctClassifications uint64 = 0, 0

    for _, resp := range responses {
        if resp.ClassNode != nil {
            classifications++
        }
        if rg.recNodeSearch(resp.ClassNode, resp.ReadSAM, useStrict) {
            correctClassifications++
        }
    }
    return 100 * (float64(correctClassifications) / float64(classifications))
}

/*
   Returns the number of non-ambiguous (non-n-containing), non-unique kmers in the reference genome.
   It is used for simple printed statistics, and to determine the amount of memory to allocate for raw-kmer-containing data structures.
   This is different from the count returned by krakenFirstPass() because it does not allow ambiguous kmers.
*/
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
