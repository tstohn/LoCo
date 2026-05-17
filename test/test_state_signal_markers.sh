#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Run LOCO
# ============================================================

./bin/loco \
    -i ./test/data3.tsv \
    -o bin/ \
    -p data3 \
    -c \
    -n 1 \
    -s 200 \
    -x 0.4 \
    -a 0 \
    -z 1 \
    -t 10 \
    -w test/data3_signalMarkers.txt \
    -v test/data3_stateMarkers.txt

# ============================================================
# Validate correlations file
# ============================================================

FILE="bin/data3_correlations.tsv"

if [[ ! -f "$FILE" ]]; then
    echo "ERROR: correlations file not found: $FILE"
    exit 1
fi

echo "Checking correlation constraints in $FILE"

awk '
BEGIN {
    FS="\t"
    fail=0
}

NR==1 {
    # Store header names
    for(i=1; i<=NF; i++) {
        header[i]=$i
    }
    next
}

{
    rowName=$1

    for(i=2; i<=NF; i++) {

        val=$i
        col=header[i]

        # Columns containing CORR_5
        if(col ~ /CORR_5/) {

            if(val >= -0.5) {
                printf("FAIL: row=%s column=%s value=%f expected < -0.5\n",
                       rowName, col, val)
                fail=1
            }

        } else {

            if(val <= 0.5) {
                printf("FAIL: row=%s column=%s value=%f expected > 0.5\n",
                       rowName, col, val)
                fail=1
            }
        }
    }
}

END {
    if(fail) {
        print "TEST FAILED"
        exit 1
    } else {
        print "TEST PASSED"
    }
}
' "$FILE"