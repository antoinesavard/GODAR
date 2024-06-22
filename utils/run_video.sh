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
    echo "Usage: sh $0 <input_file>"
    exit 1
fi

path_to_file=$(search_file "$1")
if [ -e "$path_to_file" ]; then
    echo "File $1 found at:"
    echo "${path_to_file}"
else
    echo "File $1 not found."
fi

# Define the Python script
python_script="../tools/plotter/video.py"

# Get the absolute path to the parent directory of 'tools'
PROJECT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../"
    pwd
)"

# Read the arguments from the input file and launch the Python script
cat <$1 | {
    while IFS=' ' read -r arg1 arg2; do
        echo "Launching $python_script with arguments $arg1 and $arg2"

        # Run the Python script with the arguments
        PYTHONPATH="$PROJECT_DIR" python "$python_script" "$arg1" "$arg2" &

    done

    # Wait for all background jobs to finish
    wait
}

echo "All instances of the script have been executed."
