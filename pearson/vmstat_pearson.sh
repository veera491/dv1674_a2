#!/bin/bash

# Initialize counter for incremental numerical naming
counter=1

# Function to monitor CPU, I/O, and memory usage for the Pearson process
monitor_metrics() {
    local pid=$1
    local file=$2
    local start_time=$3

    # Start monitoring CPU, I/O, and memory usage for the given process
    echo "Monitoring process ID: $pid for file: $file"

    # Increment the counter for unique filenames
    ((counter++))

    # Start monitoring and log the metrics
    pidstat -p $pid 1 > "${counter}_metrics_${file%.*}_pearson.log" &
    PIDSTAT_PID=$!

    iostat 1 > "${counter}_disk_${file%.*}_pearson.log" &
    IOSTAT_PID=$!

    vmstat 1 > "${counter}_vmstat_${file%.*}_pearson.log" &
    VMSTAT_PID=$!

    # Wait for the Pearson process to complete
    wait $pid

    # Calculate and display elapsed time
    end_time=$(date +%s%3N)
    elapsed_ms=$(( end_time - start_time ))
    printf "Run completed for file %s in %d.%03d seconds (%d ms).\n" \
        "$file" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"

    # Cleanup monitoring processes
    cleanup_monitoring
}

# Function to clean up monitoring processes
cleanup_monitoring() {
    echo "Stopping monitoring..."
    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null
}

# Define an array of dataset files
DATA_FILES=("128.data" "256.data" "512.data" "1024.data")

# Loop through each dataset file
for FILE in "${DATA_FILES[@]}"
do
    # Define the output filename based on the input file
    OUTPUT_FILE="output_${FILE%.*}_pearson.data"

    # Capture start time
    iteration_start=$(date +%s%3N)
    start_time=$(date +%s%3N)

    # Start Pearson process
    echo "Running Pearson for file $FILE..."
    ./pearson "./data/$FILE" "./data/$OUTPUT_FILE" &
    PEARSON_PID=$!

    # Start monitoring the Pearson process
    monitor_metrics $PEARSON_PID $FILE $start_time

    echo "Running Valgrind Callgrind for file $FILE..."
    callgrind_start=$(date +%s%3N)
    valgrind --tool=callgrind ./pearson "./data/$FILE" "./data/$OUTPUT_FILE"
    callgrind_end=$(date +%s%3N)
    callgrind_elapsed=$(( callgrind_end - callgrind_start ))
    printf "Callgrind completed for file %s in %d.%03d seconds (%d ms).\n" \
        "$FILE" $(( callgrind_elapsed / 1000 )) $(( callgrind_elapsed % 1000 )) "$callgrind_elapsed"

    iteration_end=$(date +%s%3N)
    iteration_elapsed=$(( iteration_end - iteration_start ))
    printf "Total processing time for %s: %d.%03d seconds (%d ms).\n" \
        "$FILE" $(( iteration_elapsed / 1000 )) $(( iteration_elapsed % 1000 )) "$iteration_elapsed"

    # Optional: Sleep for a bit to avoid overwhelming the system
    sleep 2
done

echo "All dataset files have been processed and monitored!"
