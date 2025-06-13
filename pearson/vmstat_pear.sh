#!/bin/bash

# Initialize counter for incremental numerical naming
counter=1

# Function to monitor CPU, I/O, and memory usage for the Pearson process
monitor_metrics() {
    local pid=$1
    local file=$2

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
    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))
    echo "Run completed for file $file  $elapsed seconds."


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
        start_time=$(date +%s)

    # Start Pearson process
    echo "Running Pearson for file $FILE..."
    ./pearson "./data/$FILE" "./data/$OUTPUT_FILE" &
    PEARSON_PID=$!

    # Start monitoring the Pearson process
    monitor_metrics $PEARSON_PID $FILE

    echo "Running Valgrind Callgrind for file $FILE..."
    valgrind --tool=callgrind ./pearson "./data/$FILE" "./data/$OUTPUT_FILE"

    # Optional: Sleep for a bit to avoid overwhelming the system
    sleep 2
done

echo "All dataset files have been processed and monitored!"
