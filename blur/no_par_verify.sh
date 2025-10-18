#!/bin/bash

echo "NOTE: This script verifies the sequential 'blur' binary."

status=0
red=$(tput setaf 1)
green=$(tput setaf 2)
reset=$(tput sgr0)

# Define the reference output files that are already in your data_o directory
declare -A ref_files
ref_files["im1"]="output_im1_r15.ppm"
ref_files["im2"]="output_im2_r15.ppm"
ref_files["im3"]="output_im3_r15.ppm"
ref_files["im4"]="output_im4_r15.ppm"

# Loop through each image
for image in im1 im2 im3 im4
do
    # Define a temporary output file for the test
    new_output_file="./data_o/${image}_seq_test.ppm"

    # Define the known-good reference file to compare against
    reference_file="./data_o/${ref_files[$image]}"

    echo "--- Testing sequential blur with $image.ppm ---"

    # Run the sequential blur program with its 3 required arguments
    ./blur 15 "data/$image.ppm" "$new_output_file"

    # --- Verification Steps ---

    # 1. Check if the reference file exists
    if [ ! -f "$reference_file" ]; then
        echo "${red}Error: Reference file $reference_file not found!${reset}"
        status=1
    # 2. Check if the blur program actually created an output file
    elif [ ! -f "$new_output_file" ]; then
        echo "${red}Error: Sequential program failed to create $new_output_file ${reset}"
        status=1
    # 3. Compare the new output against the reference output
    elif ! cmp -s "$reference_file" "$new_output_file"
    then
        echo "${red}Error: Incongruent output data detected for image $image.ppm.${reset}"
        echo "       (Reference: $reference_file, Generated: $new_output_file)"
        status=1
    else
        echo "${green}Success: Image $image.ppm output is correct.${reset}"
    fi

    # Clean up the temporary test file
    if [ -f "$new_output_file" ]; then
        rm "$new_output_file"
    fi

    echo "" # Add a blank line for readability
done

exit $status