package bioutils

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

func ceilDiv_U64(a, b uint64) uint64 {
    if a == 0 {
        return 0
    } else {
        return ((a - 1) / b) + 1
    }
}

func ceilDiv_64(a, b int64) int64 {
    return ((a - 1) / b) + 1
}

func ceilDiv(a, b int) int {
    return ((a - 1) / b) + 1
}
