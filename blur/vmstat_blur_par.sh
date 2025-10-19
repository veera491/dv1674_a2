#!/bin/bash

# Counter for unique log file names
counter=1

# Function to monitor a running process
# Now accepts thread count for logging purposes
monitor_metrics() {
    local pid=$1
    local file=$2
    local radius=$3
    local thread=$4
    local start_time=$5

    echo "Monitoring PID: $pid (File: $file, Radius: $radius, Threads: $thread)"

    # Base name for log files, includes file, radius, and thread count
    local base_log_name="${counter}_${file%.*}_r${radius}_t${thread}"
    ((counter++))

    # Use -t with pidstat to monitor all threads of the process
    pidstat -p $pid -t 1 > "${base_log_name}_pidstat.log" &
    PIDSTAT_PID=$!

    iostat 1 > "${base_log_name}_iostat.log" &
    IOSTAT_PID=$!

    vmstat 1 > "${base_log_name}_vmstat.log" &
    VMSTAT_PID=$!

    # Wait for the main blur_par process to finish
    wait $pid
    local exit_status=$?

    local end_time=$(date +%s%3N)
    local elapsed_ms=$(( end_time - start_time ))

    # Stop the monitoring processes
    kill $PIDSTAT_PID $IOSTAT_PID $VMSTAT_PID 2>/dev/null

    if [ $exit_status -ne 0 ]; then
        printf "Run FAILED (Exit: %d) for file %s, r:%s, t:%s in %d.%03d s (%d ms).\n" \
            "$exit_status" "$file" "$radius" "$thread" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"
    else
        printf "Run completed for file %s, r:%s, t:%s in %d.%03d s (%d ms).\n" \
            "$file" "$radius" "$thread" $(( elapsed_ms / 1000 )) $(( elapsed_ms % 1000 )) "$elapsed_ms"
    fi
}

# --- Main Script ---

INPUT_FILES=("im1.ppm" "im2.ppm" "im3.ppm" "im4.ppm")
THREADS=(1 2 4 8 16 32)
RADIUS=15

# Ensure output directory exists
mkdir -p ./data_o

# Outer loop for thread counts
for THREAD in "${THREADS[@]}"; do

    # Inner loop for input files
    for FILE in "${INPUT_FILES[@]}"; do

        # REMOVED 'local' from here
        base_name="${FILE%.*}"
        # REMOVED 'local' from here
        output_file_name="blur_${base_name}_par.ppm"

        echo "====================================================================="
        echo "Running vmstat monitor for $FILE, Radius=$RADIUS, Threads=$THREAD"
        echo "====================================================================="

        # REMOVED 'local' from here
        start_time=$(date +%s%3N)

        # Run blur_par in the background
        # This command is now correct because $output_file_name is set properly
        ./blur_par $RADIUS "./data/$FILE" "./data_o/$output_file_name" $THREAD &
        BLUR_PID=$!

        # Monitor the background process
        monitor_metrics $BLUR_PID $FILE $RADIUS $THREAD $start_time


        # --- Valgrind Callgrind Run ---

        echo "---------------------------------------------------------------------"
        echo "Running Valgrind Callgrind for $FILE, Radius=$RADIUS, Threads=$THREAD"
        echo "---------------------------------------------------------------------"

        # REMOVED 'local' from here
        callgrind_log="callgrind.${base_name}_r${RADIUS}_t${THREAD}.out"
        # REMOVED 'local' from here
        callgrind_start=$(date +%s%3N)

        # Run valgrind in the foreground
        # This command is also correct now
        valgrind --tool=callgrind --callgrind-out-file="$callgrind_log" \
            ./blur_par $RADIUS "./data/$FILE" "./data_o/$output_file_name" $THREAD

        # REMOVED 'local' from here
        callgrind_end=$(date +%s%3N)
        # REMOVED 'local' from here
        callgrind_elapsed=$(( callgrind_end - callgrind_start ))

        # This will now print the correct time
        printf "Callgrind completed for %s, r:%s, t:%s in %d.%03d s (%d ms). Log: %s\n" \
            "$FILE" "$RADIUS" "$THREAD" \
            $(( callgrind_elapsed / 1000 )) $(( callgrind_elapsed % 1000 )) "$callgrind_elapsed" \
            "$callgrind_log"

        # Clean up the generated PPM file (as in verify.sh)
        # UNCOMMENTED this line
        rm "./data_o/$output_file_name"

        sleep 2 # Give the system 2 seconds to settle
    done
done

echo "All input files have been processed for all thread counts."