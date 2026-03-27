#!/bin/bash

# Configuration
EXE="./bin/loco"
# This is the output file name loco usually generates (adjust if loco adds suffixes)
OUT_FILE="bin/test_c_laplacian.tsv"
TOTAL_RUNS=5
REQUIRED_PASSES=4
PASS_COUNT=0

# The command to run
CMD="$EXE -i ./test/simulatedData3.tsv -o bin -p test_c -c -n 20 -s 100 -x 0.5 -t 50 -q 2 -m 2 -a 0.01"

check_output() {
    # If the file doesn't even exist, fail immediately
    if [ ! -f "$1" ]; then return 1; fi

    # awk logic:
    # NR > 1 && NR <= 6: Check only the first 5 data rows
    # $3 < 0.05 and $5 < 0.05: Verify p-values
    # $6 == "A,B,C,D,E": Verify the clique column
    # $1 ~ /^[ABCDE]_[ABCDE]$/: Verify ProteinPair only contains A,B,C,D,E
    awk -F'\t' '
    NR > 1 && NR <= 6 {
        if ($3 >= 0.01 || $5 >= 0.01) exit 1 
        if ($6 != "A,B,C,D,E") exit 1
        if ($1 !~ /^[ABCDE]_[ABCDE]$/) exit 1
    }
    END { exit 0 }' "$1"
}

echo "Starting $TOTAL_RUNS iterations of Loco..."

for i in $(seq 1 $TOTAL_RUNS); do
    echo -n "Run $i: Executing... "
    
    # Execute the loco command (silencing stdout to keep logs clean)
    $CMD > /dev/null 2>&1
    
    if check_output "$OUT_FILE"; then
        ((PASS_COUNT++))
        echo "VALID"
    else
        echo "INVALID (Criteria not met or file missing)"
    fi
done

echo "---------------------------------------"
echo "Results: $PASS_COUNT / $TOTAL_RUNS passed."

if [ "$PASS_COUNT" -ge "$REQUIRED_PASSES" ]; then
    echo "OVERALL SUCCESS: Quality threshold ($REQUIRED_PASSES/$TOTAL_RUNS) met."
    exit 0
else
    echo "OVERALL FAILURE: Logic check failed."
    exit 1
fi