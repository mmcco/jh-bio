package bioutils

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

/*
   RepeatMasker is a program that takes as input a set of reference repeat
   sequences and a reference genome. It outputs "matches", specific instances
   of supplied reference repeats in the supplied reference genome. These are
   stored in the file <genome-name>.fa.out, and are parsed line-by-line into
   values of this type.

   Match.SW_Score - Smith-Waterman score, describing the likeness of
       this match to the repeat reference sequence.
   Match.PercDiv - "% substitutions in matching region compared to the
       consensus" - RepeatMasker docs
   Match.PercDel - "% of bases opposite a gap in the query sequence
       (deleted bp)" - RepeatMasker docs
   Match.PercIns - "% of bases opposite a gap in the repeat consensus
       (inserted bp)" - RepeatMasker docs
   Match.SeqName - The name (without ".fa") of the reference genome
       FASTA file this match came from. It is typically the chromosome
       name, such as "chr2L". This is inherently unsound, as RepeatMasker
       gives only a 1-dimensional qualification, but FASTA-formatted
       reference genomes are 2-dimensional, using both filename and
       sequence name.
       Reference genome FASTA files generally contain only a single
       sequence, with the same name as the file. If this is not the
       case when parsing a FASTA reference genome, we print an ominous
       warning to stdout and use the sequence name.
   Match.SeqStart - The match's start index (inclusive and
       zero-indexed) in the reference genome. Note that RepeatMasker's
       output is one-indexed.
   Match.SeqEnd - The end index (exclusive and zero-indexed) in the
       reference genome.
   Match.SeqRemains - The number of bases past the end of the match in
       the relevant reference sequence.
   Match.IsRevComp - Whether the match was for the reverse complement
       of the reference repeat sequence. In this case, we manually adjust
       some location fields, as RepeatMasker's output gives indexes for
       the reverse complement of the reference sequence. This allows us to
       treat all matches with the same logic.
   Match.RepeatClass - The repeat's full ancestry in a slice of
       strings. This includes its repeat class and repeat name, which are
       listed separately in the RepeatMasker output file. Root is implicit
       and excluded.
   Match.RepeatStart - The start index (inclusive and zero-indexed) of
       this match in the repeat consensus sequence. A signed integer is
       used because it can be negative in weird cases.
   Match.RepeatEnd - The end sequence (exclusive and zero-indexed) of
       this match in the consensus repeat sequence. A signed integer is
       used, in agreement with Match.RepeatStart.
   Match.RepeatRemains - The number of bases at the end of the
       consensus repeat sequence that this match excludes.
   Match.InsertionID - A numerical ID that is the same only for
       matches of the same long terminal repeat (LTR) instance. The
       sequences classified as <LTR name>_I or <LTR name>_int are the
       internal sequences of LTRs. These are less well defined than the
       core LTR sequence. These IDs begin at 1.

   The below fields are not parsed, but rather calculated:
   Match.RepeatName - Simply Match.RepeatClass's items concatenated.
       It is used for quick printing, and is not necessarily going to
       remain in the long-term.
   Match.ClassNode - A pointer to the match's corresponding ClassNode
   in RepeatGenome.ClassTree.
   Match.ID - A unique ID, used as a quick way of referencing and
       indexing a repeat. Pointers can generally be used, so this isn't
       necessarily going to remain in the long-term.
*/
type Match struct {
	SW_Score      int32
	PercDiv       float64
	PercDel       float64
	PercIns       float64
	SeqName       string
	SeqStart      uint64
	SeqEnd        uint64
	SeqRemains    uint64
	IsRevComp     bool
	RepeatClass   []string
	RepeatStart   int64
	RepeatEnd     int64
	RepeatRemains int64
	InsertionID   uint64

	// these are generated, not parsed
	RepeatName string
	ID         uint64
}

type Matches []Match

func ParseMatches(reader io.Reader) (error, Matches) {

	lineScanner := bufio.NewScanner(reader)
	// drop the header lines
	for i := 0; i < 3; i++ {
		if !lineScanner.Scan() {
			return fmt.Errorf("repeatgenome.parseMatches(): error reading header lines of RepeatMasker file"), nil
		}
	}

	var matches Matches

	for lineScanner.Scan() {
		matchLine := lineScanner.Text()
		rawVals := strings.Fields(matchLine)
		if len(rawVals) != 15 {
			return fmt.Errorf("repeatgenome.parseMatches(): supplied match line in RepeatMasker file is not 15 fields long (has %d fields and length %d)",
					len(rawVals), len(matchLine)),
				nil
		}
		var match Match
		match.IsRevComp = rawVals[8] == "C"

		// TODO: in the future, check to ensure that the parentheses exist
		// should be added
		// TODO: it would also be sensible to check that rawVals[8] is either
		// "C" or "+"
		// remove enclosing parens
		rawVals[7] = rawVals[7][1 : len(rawVals[7])-1]
		if match.IsRevComp {
			rawVals[11] = rawVals[11][1 : len(rawVals[11])-1]
		} else {
			rawVals[13] = rawVals[13][1 : len(rawVals[13])-1]
		}

		/*
		   Everything in this block is just typical trimming, converting, and
		   error checking.
		*/
		sw_Score, err := strconv.ParseInt(rawVals[0], 10, 32)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.SW_Score = int32(sw_Score)
		match.PercDiv, err = strconv.ParseFloat(rawVals[1], 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.PercDel, err = strconv.ParseFloat(rawVals[2], 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.PercIns, err = strconv.ParseFloat(rawVals[3], 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.SeqName = strings.TrimSpace(rawVals[4])
		match.SeqStart, err = strconv.ParseUint(rawVals[5], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.SeqEnd, err = strconv.ParseUint(rawVals[6], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.SeqRemains, err = strconv.ParseUint(rawVals[7], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		// match.IsComplement, rawVals[8], moved above
		match.RepeatClass = append(strings.Split(strings.TrimSpace(rawVals[10]), "/"), strings.TrimSpace(rawVals[9]))
		match.RepeatStart, err = strconv.ParseInt(rawVals[11], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.RepeatEnd, err = strconv.ParseInt(rawVals[12], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.RepeatRemains, err = strconv.ParseInt(rawVals[13], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}
		match.InsertionID, err = strconv.ParseUint(rawVals[14], 10, 64)
		if err != nil {
			return fmt.Errorf("repeatgenome.parseMatches():" + err.Error()), nil
		}

		/*
		   Necessary swaps to convert reverse complement repeat indexes to
		   positive-strand indexes.
		*/
		if match.IsRevComp {
			match.RepeatStart = match.RepeatRemains
			match.RepeatEnd = match.RepeatStart
			match.RepeatRemains = match.RepeatRemains + (match.RepeatEnd - match.RepeatStart)
		}

		/*
		   Decrement match.SeqStart and match.RepeatStart so that they
		   work from a start index of 0 rather than 1.
		   That way, we can use them without modification in slices.
		*/
		match.SeqStart--
		match.RepeatStart--

		/*
		   "Other" and "Unknown" classes are heirarchically meaningless and
		   really just mean "root", so we remove them
		*/
		if match.RepeatClass[0] == "Other" || match.RepeatClass[0] == "Unknown" {
			match.RepeatClass = match.RepeatClass[1:]
		}

		match.RepeatName = strings.Join(match.RepeatClass, "/")

		match.ID = uint64(len(matches))

		matches = append(matches, match)
	}

	err := lineScanner.Err()
	if err != nil {
		return fmt.Errorf("repeatgenome.parseMatches(): error reading RepeatMasker file: %s", err), nil
	}

	return nil, matches
}

// Returns the size in bases of a repeat instance.
func (match Match) Size() uint64 {
	return match.SeqEnd - match.SeqStart
}

// Writes a representation of the Matches to the supplied filename.
func (matches Matches) Write(filename string) error {
	outfile, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("Matches.Write():" + err.Error())
	}
	defer outfile.Close()

	fmt.Fprintln(outfile, "ID\tClassNode_ID\tSize\tChrom\tSW_Score")
	fmt.Fprintln(outfile)
	for _, match := range matches {
		//fmt.Fprintf(outfile, "%d\t%d\t%d\t%s\t%d", match.ID, match.ClassNode.ID, match.SeqEnd-match.SeqStart, match.SeqName, match.SW_Score)
		fmt.Fprintf(outfile, "%d\t%d\t%s\t%d", match.ID, match.SeqEnd-match.SeqStart, match.SeqName, match.SW_Score)
		fmt.Fprintln(outfile)
	}

	return nil
}
