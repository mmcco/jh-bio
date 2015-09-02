package bioutils

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
)

/*
   Ignores semicolon syntax for now.

   Could later use speed hack of accumulating slices and calling bytes.Join() at the end.
*/
func ReadFASTA(reader io.Reader) (error, map[string][]byte) {
	fileBytes, err := ioutil.ReadAll(reader)
	if err != nil {
		return err, nil
	}

	seqMap := make(map[string][]byte)
	var seqName string
	for _, line := range bytes.Split(fileBytes, []byte{'\n'}) {

		line = bytes.TrimSpace(line)
		if len(line) == 0 {
			continue
		}

		// Hit an identifier line
		if line[0] == '>' {
			seqName = string(line[1:])

		} else {
			// Append here since multi-line DNA strings are possible
			if seqName == "" {
				return fmt.Errorf("bioutils.ReadFASTAFile(): missing sequence name in FASTA file"), nil
			}
			for _, byt := range line {
				seqMap[seqName] = append(seqMap[seqName], byt)
			}
		}
	}

	return nil, seqMap
}

/*
func getLine(reader bufio.Reader) []byte {
    line, lineFrag []byte
    err error
    isPrefix := false

    for !isPrefix {
        lineFrag, isPrefix, err = reader.ReadLine()
        if err != nil && err != io.EOF {
            return err, line
        }
        for byt in lineFrag {
            line = append(line, byt)
        }
        if err == io.EOF {
            return err, line
        }
    }

    return line
}
*/
