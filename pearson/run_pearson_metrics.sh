#!/bin/bash

set -e

EXECUTABLE="./pearson"
DATA_DIR="./pearson/data"
OUTPUT_DIR="./pearson/data_o"
METRICS_FILE="metrics.csv"
INPUT_SIZES=("128" "256" "512" "1024")

# Use sudo perf if needed
PERF_CMD="perf"
if ! perf stat -e cycles true &>/dev/null; then
    PERF_CMD="sudo perf"
    echo "Using sudo for perf"
fi

echo "Compiling with profiling..."
make clean && make CXXFLAGS="-pg -O2 -Wno-sign-compare"

echo "Writing metrics to: $METRICS_FILE"
echo "input,cycles,instructions,real_time,user_time,sys_time,cpu_percent,max_mem_mb,vol_ctx,invol_ctx,page_faults,read_calls,write_calls,read_time,write_time,read_pct,write_pct,pearson_time_pct" > "$METRICS_FILE"

for size in "${INPUT_SIZES[@]}"; do
  input_file="${DATA_DIR}/${size}.data"
  output_file="${OUTPUT_DIR}/${size}_seq.data"

  echo " Profiling input: $size"

  [[ ! -f "$input_file" ]] && echo "$size,MISSING" >> "$METRICS_FILE" && continue

  # Cleanup old profiling files
  rm -f gmon.out time.log strace.log perf_output.tmp gprof.txt

  # Run with all profiling tools
  /usr/bin/time -v -o time.log \
  $PERF_CMD stat -e cycles,instructions -o perf_output.tmp \
  strace -c -o strace.log \
  "$EXECUTABLE" "$input_file" "$output_file" > /dev/null 2>&1

  # Generate gprof report
  gprof "$EXECUTABLE" gmon.out > gprof.txt

  # Extract perf metrics
  CYCLES=$(grep "cycles" perf_output.tmp | awk '{print $1}' | tr -d ',')
  INSTRS=$(grep "instructions" perf_output.tmp | awk '{print $1}' | tr -d ',')

  # Extract time log metrics
  REAL_TIME=$(grep "Elapsed (wall clock) time" time.log | awk '{print $8}')
  USER_TIME=$(grep "User time (seconds)" time.log | awk '{print $5}')
  SYS_TIME=$(grep "System time (seconds)" time.log | awk '{print $5}')
  CPU_PERCENT=$(grep "Percent of CPU this job got" time.log | awk '{print $7}')
  MAX_MEM_KB=$(grep "Maximum resident set size" time.log | awk '{print $6}')
  MAX_MEM_MB=$(awk "BEGIN {printf \"%.2f\", $MAX_MEM_KB/1024}")

  VOL_CTX=$(grep "Voluntary context switches" time.log | awk '{print $5}')
  INVOL_CTX=$(grep "Involuntary context switches" time.log | awk '{print $5}')
  PAGE_FAULTS=$(grep "Page faults" time.log | awk '{print $3}')

  # Extract strace metrics
  READ_CALLS=$(grep "^read" strace.log | awk '{print $2}')
  WRITE_CALLS=$(grep "^write" strace.log | awk '{print $2}')
  READ_TIME=$(grep "^read" strace.log | awk '{print $4}')
  WRITE_TIME=$(grep "^write" strace.log | awk '{print $4}')
  TOTAL_SYSCALL_TIME=$(awk '/^% time/ {getline; print $2}' strace.log)
  [[ -z "$TOTAL_SYSCALL_TIME" ]] && TOTAL_SYSCALL_TIME=$(awk '/^read/ {sum+=$4} END {print sum}' strace.log)

  if (( $(echo "$TOTAL_SYSCALL_TIME > 0" | bc -l) )); then
    READ_PCT=$(awk "BEGIN {printf \"%.2f\", ($READ_TIME/$TOTAL_SYSCALL_TIME)*100}")
    WRITE_PCT=$(awk "BEGIN {printf \"%.2f\", ($WRITE_TIME/$TOTAL_SYSCALL_TIME)*100}")
  else
    READ_PCT=0
    WRITE_PCT=0
  fi

  # Extract gprof time percentage
  PEARSON_TIME_PCT=$(grep -A1 "pearson" gprof.txt | grep "[0-9]\+" | awk '{print $1}')
  [[ -z "$PEARSON_TIME_PCT" ]] && PEARSON_TIME_PCT="N/A"

  # Write to CSV
  echo "$size,$CYCLES,$INSTRS,$REAL_TIME,$USER_TIME,$SYS_TIME,$CPU_PERCENT,$MAX_MEM_MB,$VOL_CTX,$INVOL_CTX,$PAGE_FAULTS,$READ_CALLS,$WRITE_CALLS,$READ_TIME,$WRITE_TIME,$READ_PCT,$WRITE_PCT,$PEARSON_TIME_PCT" >> "$METRICS_FILE"

done

echo "All metrics collected. See: $METRICS_FILE"