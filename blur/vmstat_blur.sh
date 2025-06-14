#!/bin/bash

counter=1

monitor_metrics() {
    local pid=$1
    local file=$2
    local radius=$3
    local start_time=$4

    echo "Monitoring process ID: $pid for file: $file with radius: $radius"

    ((counter++))

    pidstat -p $pid 1 > "${counter}_metrics_${file%.*}_r${radius}.log" &
    PIDSTAT_PID=$!

    iostat 1 > "${counter}_disk_${file%.*}_r${radius}.log" &
    IOSTAT_PID=$!

    vmstat 1 > "${counter}_vmstat_${file%.*}_r${radius}.log" &
    VMSTAT_PID=$!

    wait $pid

    end_time=$(date +%s)
    elapsed=$(( end_time - start_time ))
    echo "Run completed for file $file with radius $radius in $elapsed seconds."

    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null
}

INPUT_FILES=("im1.ppm" "im2.ppm" "im3.ppm" "im4.ppm")
RADIUS=15

for FILE in "${INPUT_FILES[@]}"; do
    OUTPUT_FILE="output_${FILE%.*}_r${RADIUS}.ppm"
    start_time=$(date +%s)

    echo "Running blur for file $FILE with radius $RADIUS..."
    ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE" &
    BLUR_PID=$!

    monitor_metrics $BLUR_PID $FILE $RADIUS $start_time

    echo "Running Valgrind Callgrind for file $FILE with radius $RADIUS..."
    valgrind --tool=callgrind ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE"

    sleep 2
done

echo "All input files have been processed at radius 15 for baseline measurements."
