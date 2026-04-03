#!/bin/bash

# --- Configuration ---
EXE="./bin/loco"
OUT_FILE="bin/data_1_laplacian.tsv"
TOTAL_RUNS=5
REQUIRED_PASSES=4
PASS_COUNT=0
# Your specific command
CMD="$EXE -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s 50 -x 0.3 -z 1 -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0"

# Store results for a final summary
RESULTS=()

check_families() {
    # 1. Check if file exists
    if [ ! -f "$1" ]; then return 1; fi

    # 2. Use awk to check top 5 rows (NR 2 to 6)
    # Returns 0 if all 5 rows pass, 1 if any fail
    awk -F'\t' '
    NR > 1 && NR <= 6 {
        # Check P-values
        if ($3 >= 0.01 || $5 >= 0.01) exit 1
        # Check Family Consistency (E_E or M_M only)
        if ($1 !~ /^E[0-9]+_E[0-9]+$/ && $1 !~ /^M[0-9]+_M[0-9]+$/) exit 1
    }
    END { exit 0 }' "$1"
}

echo "================================================="
echo "   LOCO CONSISTENCY TEST: 5 RUNS REQUIRED"
echo "================================================="

for i in $(seq 1 $TOTAL_RUNS); do
    echo -e "\n>>> RUN $i of $TOTAL_RUNS..."
    
    # CRITICAL: Remove old output so we don't test stale data
    rm -f "$OUT_FILE"

    # Execute the command
    eval "$CMD"
    STATUS=$?

    if [ $STATUS -ne 0 ]; then
        echo "LOGIC: Binary CRASHED (Exit Code: $STATUS)"
        RESULTS[$i]="CRASH"
    else
        if check_families "$OUT_FILE"; then
            echo "LOGIC: PASS (Top 5 are valid E-E/M-M)"
            RESULTS[$i]="PASS"
            ((PASS_COUNT++))
        else
            echo "LOGIC: FAIL (Invalid family mixing or p-values)"
            RESULTS[$i]="FAIL"
        fi
    fi
done

# --- Final Summary Table ---
echo -e "\n================================================="
echo "                FINAL SUMMARY"
echo "================================================="
for i in $(seq 1 $TOTAL_RUNS); do
    echo "Run $i: ${RESULTS[$i]}"
done
echo "-------------------------------------------------"
echo "Total Successful Runs: $PASS_COUNT / $TOTAL_RUNS"
echo "Target Required:      $REQUIRED_PASSES"

# --- Final Exit Decision ---
if [ "$PASS_COUNT" -ge "$REQUIRED_PASSES" ]; then
    echo -e "\033[0;32mOVERALL RESULT: SUCCESS\033[0m"
    exit 0
else
    echo -e "\033[0;31mOVERALL RESULT: FAILURE (Did not reach 4/5 goal)\033[0m"
    exit 1
fi