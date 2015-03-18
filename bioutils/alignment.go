/*
    Adapted from Ant Zucaro's matchr package:
        https://github.com/antzucaro/matchr

    matchr is GPLv2 licensed.
*/
package bioutils

import (
    "fmt"
)

const GAP_COST float64 = 0.5

func getCost(a []byte, aIndex int, b []byte, bIndex int) float64 {
	if a[aIndex] == b[bIndex] {
		return 1.0
	} else {
		return -2.0
	}
}

func NeedlemanWunsch(a, b []byte) float64 {

	if len(a) == 0 {
		return float64(len(b))
	}

	if len(b) == 0 {
		return float64(len(a))
	}

	d := make([][]float64, len(a)+1)
	for i := range d {
		d[i] = make([]float64, len(b)+1)
	}

	for i := 0; i < len(a)+1; i++ {
        d[i][0] = float64(-i) * GAP_COST
	}

	for j := 0; j < len(b)+1; j++ {
        d[0][j] = float64(-j) * GAP_COST
	}

	for i := 1; i < len(a)+1; i++ {
		for j := 1; j < len(b)+1; j++ {
            cost := getCost(a, i-1, b, j-1)

			// find the lowest cost
			d[i][j] = maxF64(
				d[i-1][j]-GAP_COST,
				maxF64(d[i][j-1]-GAP_COST, d[i-1][j-1]+cost))
		}
	}

    printMatrix(d)

	return d[len(a)][len(b)]
}


func SmithWaterman(a, b []byte) float64 {

	if len(a) == 0 {
		return float64(len(b))
	}

	if len(b) == 0 {
		return float64(len(a))
	}

	d := make([][]float64, len(a)+1)
	for i := range d {
		d[i] = make([]float64, len(b)+1)
	}

	var maxSoFar float64

	for i := 1; i < len(a)+1; i++ {
		for j := 1; j < len(b)+1; j++ {
            cost := getCost(a, i-1, b, j-1)

			// find the lowest cost
			d[i][j] = maxF64(
				maxF64(0, d[i-1][j]-GAP_COST),
				maxF64(d[i][j-1]-GAP_COST, d[i-1][j-1]+cost))

			// save if it is the biggest thus far
			if d[i][j] > maxSoFar {
				maxSoFar = d[i][j]
			}
		}
	}

    printMatrix(d)

	return maxSoFar
}


/*
    Print matrix for debugging.
*/
func printMatrix(matrix [][]float64) {
    for _, row := range matrix {
        for _, score := range row {
            fmt.Printf("%.1f  ", score)
        }
        fmt.Print("\n\n")
    }
}
