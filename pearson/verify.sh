#!/bin/bash

# Determine script's directory to use absolute paths
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

status=0
if tput setaf 1 &> /dev/null; then
    red=$(tput setaf 1)
    yellow=$(tput setaf 3)
    green=$(tput setaf 2)
    reset=$(tput sgr0)
else
    red=""
    yellow=""
    green=""
    reset=""
fi

errors_found=0
warnings_found=0

data_dir="$ROOT_DIR/data"
out_dir="$ROOT_DIR/data_o"

# Run pearson to generate sequential output
"$ROOT_DIR/./pearson" "$data_dir/128.data" "$out_dir/128_seq.data"
"$ROOT_DIR/./pearson" "$data_dir/256.data" "$out_dir/256_seq.data"
"$ROOT_DIR/./pearson" "$data_dir/512.data" "$out_dir/512_seq.data"
"$ROOT_DIR/./pearson" "$data_dir/1024.data" "$out_dir/1024_seq.data"

for thread in 2 4 8 16 32; do
    for size in 128 256 512 1024; do
        # FIXED: Call pearson_par, NOT the directory!
        "$ROOT_DIR/./pearson_par" "$data_dir/${size}.data" "$out_dir/${size}_par.data" $thread

        # Run the verify binary and capture the return code
        "$ROOT_DIR/./verify" "$out_dir/${size}_seq.data" "$out_dir/${size}_par.data"
        ret=$?

        # Check the return code and print corresponding message
        if [ $ret -eq 2 ]; then
            echo "${red}ERROR: Significant mismatch found in size ${size} with ${thread} thread(s).${reset}"
            errors_found=1
        elif [ $ret -eq 1 ]; then
            echo "${yellow}WARNING: Minor differences found in size ${size} with ${thread} thread(s).${reset}"
            warnings_found=1
        elif [ $ret -eq 0 ]; then
            echo "${green}Success: Files match for size ${size} with ${thread} thread(s).${reset}"
        else
            echo "${red}ERROR: An unexpected error occurred while processing size ${size} with ${thread} thread(s).${reset}"
            errors_found=1
        fi

        # Clean up parallel output
        rm -f "$out_dir/${size}_par.data"
    done
    echo
done

# Final output based on results
if [ $errors_found -eq 1 ]; then
    echo "${red}Errors found during the tests.${reset}"
    status=1
elif [ $warnings_found -eq 1 ]; then
    echo "${yellow}Warnings were found, but no major errors.${reset}"
else
    echo "${green}Success: All tests passed successfully.${reset}"
fi

exit $status
