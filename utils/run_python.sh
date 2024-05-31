#!/bin/bash

# Define the Python script
python_script="/plots/video.py"

# Define the input file containing argument pairs
input_file="pyargs.dat"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file not found."
    exit 1
fi

# Read the arguments from the input file and launch the Python script
while IFS=' ' read -r arg1 arg2; do
    echo "Launching $python_script with arguments $arg1 and $arg2"

    # Run the Python script with the arguments
    python "$python_script" "$arg1" "$arg2" &

done <"$input_file"

# Wait for all background jobs to finish
wait

echo "All instances of the script have been executed."
