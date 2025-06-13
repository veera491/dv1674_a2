#!/bin/bash

# Paths
PEARSON_EXEC="./pearson"
DATA_DIR="./data"
OUTPUT_DIR="./data_o"
METRICS_FILE="metrics.csv"

# Ensure output dir exists
mkdir -p "$OUTPUT_DIR"

# Compile pearson with profiling enabled (-pg for gprof)
echo "Compiling pearson with profiling enabled..."
g++ -pg -std=c++17 -g -Wall vector.cpp dataset.cpp analysis.cpp pearson.cpp -o pearson

# Clean previous metrics file
rm -f "$METRICS_FILE"

# Print CSV header
echo "Input Size,Real Time (s),User Time (s),Sys Time (s),CPU (%),Max Mem (MiB),Page Faults,Voluntary Ctx Switches,Involuntary Ctx Switches,Read Calls,Write Calls,% Read Time,% Write Time,Pearson Function Time (%)" > "$METRICS_FILE"

# Loop over each data file in data dir
for input_file in "$DATA_DIR"/*.data; do
  filename=$(basename "$input_file")
  input_size=$(echo "$filename" | grep -oP '^\d+')
  output_file="$OUTPUT_DIR/${input_size}_seq.data"

  echo "Profiling input: $filename"

  # Remove old profiling files
  rm -f gmon.out strace.log time.log

  # 1. /usr/bin/time to get resource usage
  /usr/bin/time -v $PEARSON_EXEC "$input_file" "$output_file" 1> /dev/null 2> time.log

  # 2. strace to get syscall info (read/write)
  strace -T -c -e trace=read,write $PEARSON_EXEC "$input_file" "$output_file" 2> strace.log > /dev/null

  # 3. Run again to generate gmon.out for gprof
  $PEARSON_EXEC "$input_file" "$output_file" 1> /dev/null 2> /dev/null

  # Generate gprof report
  gprof $PEARSON_EXEC gmon.out > gprof.txt

  # ---- Parse metrics ----
  REAL_TIME=$(grep "Elapsed (wall clock) time" time.log | awk '{print $8}')
  USER_TIME=$(grep "User time (seconds)" time.log | awk '{print $5}')
  SYS_TIME=$(grep "System time (seconds)" time.log | awk '{print $5}')
  CPU_PERCENT=$(grep "Percent of CPU this job got" time.log | awk '{print $7}')
  MAX_MEM_KB=$(grep "Maximum resident set size" time.log | awk '{print $6}')
  MAX_MEM_MB=$(awk "BEGIN {printf \"%.2f\", $MAX_MEM_KB/1024}")

  VOL_CTX=$(grep "Voluntary context switches" time.log | awk '{print $5}')
  INVOL_CTX=$(grep "Involuntary context switches" time.log | awk '{print $5}')
  PAGE_FAULTS=$(grep "Page faults" time.log | awk '{print $3}')

  READ_CALLS=$(grep "^read" strace.log | awk '{print $2}')
  WRITE_CALLS=$(grep "^write" strace.log | awk '{print $2}')
  READ_TIME=$(grep "^read" strace.log | awk '{print $4}')
  WRITE_TIME=$(grep "^write" strace.log | awk '{print $4}')
  TOTAL_SYSCALL_TIME=$(awk '/^% time/ {getline; print $2}' strace.log)
  [ -z "$TOTAL_SYSCALL_TIME" ] && TOTAL_SYSCALL_TIME=$(awk '/^read/ {sum+=$4} END {print sum}' strace.log)

  if (( $(echo "$TOTAL_SYSCALL_TIME > 0" | bc -l) )); then
    READ_PCT=$(awk "BEGIN {printf \"%.2f\", ($READ_TIME/$TOTAL_SYSCALL_TIME)*100}")
    WRITE_PCT=$(awk "BEGIN {printf \"%.2f\", ($WRITE_TIME/$TOTAL_SYSCALL_TIME)*100}")
  else
    READ_PCT=0
    WRITE_PCT=0
  fi

  PEARSON_TIME_PCT=$(grep -A1 "pearson" gprof.txt | grep "[0-9]\+" | awk '{print $1}')
  [ -z "$PEARSON_TIME_PCT" ] && PEARSON_TIME_PCT="N/A"

  # Append metrics to CSV
  echo "$input_size,$REAL_TIME,$USER_TIME,$SYS_TIME,$CPU_PERCENT,$MAX_MEM_MB,$PAGE_FAULTS,$VOL_CTX,$INVOL_CTX,$READ_CALLS,$WRITE_CALLS,$READ_PCT,$WRITE_PCT,$PEARSON_TIME_PCT" >> "$METRICS_FILE"

done

echo "Profiling completed. Metrics saved in $METRICS_FILE"
#!/bin/bash

# Paths
PEARSON_EXEC="./pearson"
DATA_DIR="./data"
OUTPUT_DIR="./data_o"
METRICS_FILE="metrics.csv"

# Ensure output dir exists
mkdir -p "$OUTPUT_DIR"

# Compile pearson with profiling enabled (-pg for gprof)
echo "Compiling pearson with profiling enabled..."
g++ -pg -std=c++17 -g -Wall vector.cpp dataset.cpp analysis.cpp pearson.cpp -o pearson

# Clean previous metrics file
rm -f "$METRICS_FILE"

# Print CSV header
echo "Input Size,Real Time (s),User Time (s),Sys Time (s),CPU (%),Max Mem (MiB),Page Faults,Voluntary Ctx Switches,Involuntary Ctx Switches,Read Calls,Write Calls,% Read Time,% Write Time,Pearson Function Time (%)" > "$METRICS_FILE"

# Loop over each data file in data dir
for input_file in "$DATA_DIR"/*.data; do
  filename=$(basename "$input_file")
  input_size=$(echo "$filename" | grep -oP '^\d+')
  output_file="$OUTPUT_DIR/${input_size}_seq.data"

  echo "Profiling input: $filename"

  # Remove old profiling files
  rm -f gmon.out strace.log time.log

  # 1. /usr/bin/time to get resource usage
  /usr/bin/time -v $PEARSON_EXEC "$input_file" "$output_file" 1> /dev/null 2> time.log

  # 2. strace to get syscall info (read/write)
  strace -T -c -e trace=read,write $PEARSON_EXEC "$input_file" "$output_file" 2> strace.log > /dev/null

  # 3. Run again to generate gmon.out for gprof
  $PEARSON_EXEC "$input_file" "$output_file" 1> /dev/null 2> /dev/null

  # Generate gprof report
  gprof $PEARSON_EXEC gmon.out > gprof.txt

  # ---- Parse metrics ----
  REAL_TIME=$(grep "Elapsed (wall clock) time" time.log | awk '{print $8}')
  USER_TIME=$(grep "User time (seconds)" time.log | awk '{print $5}')
  SYS_TIME=$(grep "System time (seconds)" time.log | awk '{print $5}')
  CPU_PERCENT=$(grep "Percent of CPU this job got" time.log | awk '{print $7}')
  MAX_MEM_KB=$(grep "Maximum resident set size" time.log | awk '{print $6}')
  MAX_MEM_MB=$(awk "BEGIN {printf \"%.2f\", $MAX_MEM_KB/1024}")

  VOL_CTX=$(grep "Voluntary context switches" time.log | awk '{print $5}')
  INVOL_CTX=$(grep "Involuntary context switches" time.log | awk '{print $5}')
  PAGE_FAULTS=$(grep "Page faults" time.log | awk '{print $3}')

  READ_CALLS=$(grep "^read" strace.log | awk '{print $2}')
  WRITE_CALLS=$(grep "^write" strace.log | awk '{print $2}')
  READ_TIME=$(grep "^read" strace.log | awk '{print $4}')
  WRITE_TIME=$(grep "^write" strace.log | awk '{print $4}')
  TOTAL_SYSCALL_TIME=$(awk '/^% time/ {getline; print $2}' strace.log)
  [ -z "$TOTAL_SYSCALL_TIME" ] && TOTAL_SYSCALL_TIME=$(awk '/^read/ {sum+=$4} END {print sum}' strace.log)

  if (( $(echo "$TOTAL_SYSCALL_TIME > 0" | bc -l) )); then
    READ_PCT=$(awk "BEGIN {printf \"%.2f\", ($READ_TIME/$TOTAL_SYSCALL_TIME)*100}")
    WRITE_PCT=$(awk "BEGIN {printf \"%.2f\", ($WRITE_TIME/$TOTAL_SYSCALL_TIME)*100}")
  else
    READ_PCT=0
    WRITE_PCT=0
  fi

  PEARSON_TIME_PCT=$(grep -A1 "pearson" gprof.txt | grep "[0-9]\+" | awk '{print $1}')
  [ -z "$PEARSON_TIME_PCT" ] && PEARSON_TIME_PCT="N/A"

  # Append metrics to CSV
  echo "$input_size,$REAL_TIME,$USER_TIME,$SYS_TIME,$CPU_PERCENT,$MAX_MEM_MB,$PAGE_FAULTS,$VOL_CTX,$INVOL_CTX,$READ_CALLS,$WRITE_CALLS,$READ_PCT,$WRITE_PCT,$PEARSON_TIME_PCT" >> "$METRICS_FILE"

done

echo "Profiling completed. Metrics saved in $METRICS_FILE"
