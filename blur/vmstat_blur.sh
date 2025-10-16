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

    end_time=$(date +%s%3N)
    elapsed_ms=$(( end_time - start_time ))
    printf "Run completed for file %s with radius %s in %d.%03d seconds (%d ms).\n" \
        "$file" "$radius" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"

    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null
}

INPUT_FILES=("im1.ppm" "im2.ppm" "im3.ppm" "im4.ppm")
RADIUS=15

for FILE in "${INPUT_FILES[@]}"; do
    OUTPUT_FILE="output_${FILE%.*}_r${RADIUS}.ppm"
    iteration_start=$(date +%s%3N)
    start_time=$(date +%s%3N)

    echo "Running blur for file $FILE with radius $RADIUS..."
    ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE" &
    BLUR_PID=$!

    monitor_metrics $BLUR_PID $FILE $RADIUS $start_time

    echo "Running Valgrind Callgrind for file $FILE with radius $RADIUS..."
    callgrind_start=$(date +%s%3N)
    valgrind --tool=callgrind ./blur $RADIUS "./data/$FILE" "./data_o/$OUTPUT_FILE"
    callgrind_end=$(date +%s%3N)
    callgrind_elapsed=$(( callgrind_end - callgrind_start ))
    printf "Callgrind completed for file %s with radius %s in %d.%03d seconds (%d ms).\n" \
        "$file" "$RADIUS" $(( callgrind_elapsed / 1000 )) $(( callgrind_elapsed % 1000 )) "$callgrind_elapsed"

    iteration_end=$(date +%s%3N)
    iteration_elapsed=$(( iteration_end - iteration_start ))
    printf "Total processing time for %s with radius %s: %d.%03d seconds (%d ms).\n" \
        "$file" "$RADIUS" $(( iteration_elapsed / 1000 )) $(( iteration_elapsed % 1000 )) "$iteration_elapsed"

    sleep 2
done

echo "All input files have been processed at radius 15 for baseline measurements."
