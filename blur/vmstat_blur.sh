#!/bin/bash

# Initialize counter for incremental numerical naming
counter=1

# Function to monitor CPU, I/O, and memory usage for the blur process
monitor_metrics() {
    local pid=$1
    local file=$2
    local radius=$3
    local start_time=$4

    # Start monitoring CPU, I/O, and memory usage for the given blur process
    echo "Monitoring process ID: $pid for file: $file with radius: $radius"

    # Increment the counter for the next file (pre-increment)
    ((counter++))

    # Start monitoring in the background and use the incremental counter in the filenames
    pidstat -p $pid 1 > "${counter}_metrics_${file%.*}_r${radius}.log" &
    PIDSTAT_PID=$!

    iostat 1 > "${counter}_disk_${file%.*}_r${radius}.log" &
    IOSTAT_PID=$!

    vmstat 1 > "${counter}_vmstat_${file%.*}_r${radius}.log" &
    VMSTAT_PID=$!

    # Wait for the blur process to complete
    wait $pid

    # Calculate and display elapsed time
    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))
    echo "Run completed for file $file with radius $radius in $elapsed seconds."

    # Cleanup monitoring processes
    cleanup_monitoring
}

# Function to clean up monitoring processes
cleanup_monitoring() {
    echo "Stopping monitoring..."
    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null
}

# Define an array of radii to test
RADIUS_VALUES=(1 100 500 900 1000)

# Define an array of input PPM files
INPUT_FILES=("im1.ppm" "im2.ppm" "im3.ppm" "im4.ppm")

# Loop through each input file
for FILE in "${INPUT_FILES[@]}"
do
    # Loop through each radius
    for RADIUS in "${RADIUS_VALUES[@]}"
    do
        # Define the output filename based on the input file and radius
        OUTPUT_FILE="output_${FILE%.*}_r${RADIUS}.ppm"
        
        # Capture start time
        start_time=$(date +%s)

        # Run the blur program and capture its PID
        echo "Running blur for file $FILE with radius $RADIUS..."
        ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE" &
        BLUR_PID=$!

        # Start monitoring the blur process
        monitor_metrics $BLUR_PID $FILE $RADIUS $start_time

        echo "Running Valgrind Callgrind for file $FILE with radius $RADIUS..."
        valgrind --tool=callgrind ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE"


        # Optional: Sleep for a bit to avoid overwhelming the system
        sleep 2
    done
done


echo "All combinations of input files and radii have been processed!"
