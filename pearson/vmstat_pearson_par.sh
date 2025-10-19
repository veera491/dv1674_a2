#!/bin/bash

# Counter for unique log file names
counter=1

# Function to monitor a running process
# Now accepts thread count for logging purposes
monitor_metrics() {
    local pid=$1
    local file=$2    # e.g., "128.data"
    local thread=$3  # e.g., "4"
    local start_time=$4

    echo "Monitoring PID: $pid (File: $file, Threads: $thread)"

    # Base name for log files, includes file and thread count
    local base_log_name="${counter}_${file%.*}_t${thread}"
    ((counter++))

    # Use -t with pidstat to monitor all threads of the process
    pidstat -p $pid -t 1 > "${base_log_name}_pidstat.log" &
    PIDSTAT_PID=$!

    iostat 1 > "${base_log_name}_iostat.log" &
    IOSTAT_PID=$!

    vmstat 1 > "${base_log_name}_vmstat.log" &
    VMSTAT_PID=$!

    # Wait for the main pearson_par process to finish
    wait $pid
    local exit_status=$?

    local end_time=$(date +%s%3N)
    local elapsed_ms=$(( end_time - start_time ))

    # Stop the monitoring processes
    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null

    if [ $exit_status -ne 0 ]; then
        printf "Run FAILED (Exit: %d) for file %s, t:%s in %d.%03d s (%d ms).\n" \
            "$exit_status" "$file" "$thread" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"
    else
        printf "Run completed for file %s, t:%s in %d.%03d s (%d ms).\n" \
            "$file" "$thread" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"
    fi
}

# --- Main Script ---

DATA_SIZES=("128" "256" "512" "1024")
# Added 1 thread for baseline parallel/optimized sequential measurement
THREADS=(1 2 4 8 16 32)

# Ensure output directory exists
mkdir -p ./data_o

# Outer loop for thread counts
for THREAD in "${THREADS[@]}"; do

    # Inner loop for input data sizes
    for SIZE in "${DATA_SIZES[@]}"; do

        FILE_NAME="${SIZE}.data"
        INPUT_FILE="data/$FILE_NAME"
        # Unique output file name for this run
        OUTPUT_FILE="data_o/${SIZE}_par_t${THREAD}.data"

        echo "====================================================================="
        echo "Running vmstat monitor for $FILE_NAME, Threads=$THREAD"
        echo "====================================================================="

        start_time=$(date +%s%3N)

        # Run pearson_par in the background
        ./pearson_par "$INPUT_FILE" "$OUTPUT_FILE" $THREAD &
        PEARSON_PID=$!

        # Monitor the background process
        monitor_metrics $PEARSON_PID $FILE_NAME $THREAD $start_time


        # --- Valgrind Callgrind Run ---

        echo "---------------------------------------------------------------------"
        echo "Running Valgrind Callgrind for $FILE_NAME, Threads=$THREAD"
        echo "---------------------------------------------------------------------"

        # Unique callgrind output log
        callgrind_log="callgrind.${SIZE}_t${THREAD}.out"
        callgrind_start=$(date +%s%3N)

        # Run valgrind in the foreground
        valgrind --tool=callgrind --callgrind-out-file="$callgrind_log" \
            ./pearson_par "$INPUT_FILE" "$OUTPUT_FILE" $THREAD

        callgrind_end=$(date +%s%3N)
        callgrind_elapsed=$(( callgrind_end - callgrind_start ))

        printf "Callgrind completed for %s, t:%s in %d.%03d s (%d ms). Log: %s\n" \
            "$FILE_NAME" "$THREAD" \
            $(( callgrind_elapsed / 1000 )) $(( callgrind_elapsed % 1000 )) "$callgrind_elapsed" \
            "$callgrind_log"

        # Clean up the generated data file (as in verify.sh)
        rm "$OUTPUT_FILE"

        sleep 2 # Give the system 2 seconds to settle
    done
done

echo "All input files have been processed for all thread counts."