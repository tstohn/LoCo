#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
EXE="./bin/loco"
OUT_FILE="bin/data_1_laplacian.tsv"
TOTAL_RUNS=5
REQUIRED_PASSES=4
PASS_COUNT=0

CMD="$EXE -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s 50 -x 0.3 -z 1 -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0"

# Store results
RESULTS=()

check_families() {
    if [ ! -f "$1" ]; then return 1; fi

    awk -F'\t' '
    NR > 1 && NR <= 6 {
        gsub(/\r/, "", $0)

        if ($3 >= 0.01 || $5 >= 0.01) exit 1
        if ($1 !~ /^E[0-9]+_E[0-9]+$/ && $1 !~ /^M[0-9]+_M[0-9]+$/) exit 1
    }
    END { exit 0 }' "$1"
}

echo "================================================="
echo "   LOCO CONSISTENCY TEST: 5 RUNS REQUIRED"
echo "================================================="

# ✅ NO seq → portable loop
i=1
while [ "$i" -le "$TOTAL_RUNS" ]; do
    echo ""
    echo ">>> RUN $i of $TOTAL_RUNS..."

    # remove old output safely (Windows-safe)
    rm -f "$OUT_FILE" 2>/dev/null || true

    # run command safely (no eval needed)
    $CMD > /dev/null 2>&1 || true
    STATUS=$?

    if [ $STATUS -ne 0 ]; then
        echo "LOGIC: Binary CRASHED (Exit Code: $STATUS)"
        RESULTS[$i]="CRASH"
    else
        if check_families "$OUT_FILE"; then
            echo "LOGIC: PASS (Top 5 are valid E-E/M-M)"
            RESULTS[$i]="PASS"
            PASS_COUNT=$((PASS_COUNT + 1))
        else
            echo "LOGIC: FAIL (Invalid family mixing or p-values)"
            RESULTS[$i]="FAIL"
        fi
    fi

    i=$((i + 1))
done

# --- Final Summary ---
echo ""
echo "================================================="
echo "                FINAL SUMMARY"
echo "================================================="

i=1
while [ "$i" -le "$TOTAL_RUNS" ]; do
    echo "Run $i: ${RESULTS[$i]}"
    i=$((i + 1))
done

echo "-------------------------------------------------"
echo "Total Successful Runs: $PASS_COUNT / $TOTAL_RUNS"
echo "Target Required:      $REQUIRED_PASSES"

# --- Final Exit ---
if [ "$PASS_COUNT" -ge "$REQUIRED_PASSES" ]; then
    echo "OVERALL RESULT: SUCCESS"
    exit 0
else
    echo "OVERALL RESULT: FAILURE (Did not reach threshold)"
    exit 1
fi