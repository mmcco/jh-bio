/*
   Somewhat archaic functions that deal solely with TextSeq
   These operations are far slower than the equivalents performed on the corresponding int types.
   They therefore generally shouldn't be used.
*/

package repeatgenome

/*
   A representation of a genetic sequence using one byte letter per base.
   A type synonym is used to differentiate one-byte-per-base sequences from two-bits-per-base sequences (type Seq), which also use the concrete type []byte.
*/
type TextSeq []byte

/*
   Returns the reverse complement of the supplied TextSeq.
   This drains memory and should therefore not be used outside of debugging and printing.
*/
func (seq TextSeq) revComp() TextSeq {
    var revCompSeq = make(TextSeq, 0, len(seq))

    for i := range seq {
        switch seq[len(seq)-i-1] {
        case 'a':
            revCompSeq = append(revCompSeq, 't')
            break
        case 't':
            revCompSeq = append(revCompSeq, 'a')
            break
        case 'c':
            revCompSeq = append(revCompSeq, 'g')
            break
        case 'g':
            revCompSeq = append(revCompSeq, 'c')
            break
        default:
            panic("byte other than 'a', 'c', 'g', or 't' supplied to revComp")
        }
    }
    return revCompSeq
}

/* 
   Returns a bool describing whether the first TextSeq is lexicographically smaller than the second.
*/
func Less(a, b TextSeq) bool {
    size := min(len(a), len(b))

    for i := 0; i < size; i++ {
        if a[i] < b[i] {
            return true
        }
        if a[i] > b[i] {
            return false
        }
    }

    return len(a) < len(b)
}
