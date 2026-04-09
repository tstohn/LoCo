#!/usr/bin/env bash
set -euo pipefail

# Configuration
EXE="./bin/loco"
OUT_FILE="bin/test_c_laplacian.tsv"
TOTAL_RUNS=5
REQUIRED_PASSES=4
PASS_COUNT=0

CMD="$EXE -i ./test/simulatedData3.tsv -o bin -p test_c -c -n 20 -s 100 -x 0.5 -t 50 -q 2 -m 2 -a 0.01"

check_output() {
    if [ ! -f "$1" ]; then
        return 1
    fi

    awk -F'\t' '
    NR > 1 && NR <= 6 {
        gsub(/\r/, "", $0)

        if ($3 >= 0.01 || $5 >= 0.01) exit 1
        if ($6 != "A,B,C,D,E") exit 1
        if ($1 !~ /^[ABCDE]_[ABCDE]$/) exit 1
    }
    END { exit 0 }' "$1"
}

echo "Starting $TOTAL_RUNS iterations of Loco..."

# ✅ POSIX-safe loop (NO seq!)
i=1
while [ "$i" -le "$TOTAL_RUNS" ]; do
    echo -n "Run $i: Executing... "

    # run silently
    $CMD > /dev/null 2>&1 || true

    if check_output "$OUT_FILE"; then
        PASS_COUNT=$((PASS_COUNT + 1))
        echo "VALID"
    else
        echo "INVALID"
    fi

    i=$((i + 1))
done

echo "---------------------------------------"
echo "Results: $PASS_COUNT / $TOTAL_RUNS passed."

if [ "$PASS_COUNT" -ge "$REQUIRED_PASSES" ]; then
    echo "OVERALL SUCCESS: Quality threshold met."
    exit 0
else
    echo "OVERALL FAILURE: Logic check failed."
    exit 1
fi