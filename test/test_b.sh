#!/usr/bin/env bash
set -euo pipefail

FILE="bin/test_b_laplacian.tsv"

awk '
BEGIN {
    expected["A_B"]=1
    expected["A_C"]=1
    expected["B_C"]=1
    expected["D_E"]=1
    expected["D_F"]=1
    expected["E_F"]=1
}

NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "p_CorrelationScore") pcol = i
    }
    next
}

{
    gsub(/\r/, "", $0)

    # split on ANY whitespace (tabs or spaces)
    n = split($0, f, /[ \t]+/)

    pair = f[1]
    pval = f[pcol]

    #  SKIP EMPTY OR BROKEN LINES (IMPORTANT FIX)
    if (pair == "" || pval == "" || n < 2) {
        next
    }

    seen[pair] = 1

    if (pval >= 0.05) {
        print " FAIL p-value too high:", pair, pval
        exit 1
    }

    if (!(pair in expected)) {
        print " UNEXPECTED PAIR:", pair
        exit 1
    }
}

END {
    for (p in expected) {
        if (!(p in seen)) {
            print " MISSING PAIR:", p
            exit 1
        }
    }

    print "ALL CHECKS PASSED"
}
' "$FILE"