#!/bin/bash

# Function to search for a file
search_file() {
    local filename=$1
    local search_dir=${PWD}"/.."
    local location

    location=$(find "${search_dir}" -name "${filename}" 2>/dev/null)

    if [ -e "${location}" ]; then
        echo "${location}"
    fi
}

# Check if an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

echo search_file "$1"

# Define the Python script
python_script="plots/video.py"

# Read the arguments from the input file and launch the Python script
while IFS=' ' read -r arg1 arg2; do
    echo "Launching $python_script with arguments $arg1 and $arg2"

    # Run the Python script with the arguments
    python "$python_script" "$arg1" "$arg2" &

done < <(cat <"$1")

# Wait for all background jobs to finish
wait

echo "All instances of the script have been executed."
