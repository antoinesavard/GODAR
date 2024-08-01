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
python_script="../tools/init/gridded_placement.py"

# Read the arguments from the input file and launch the Python script
cat $1 | {
    while IFS=' ' read -r arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10; do
        echo "Launching $python_script with arguments $arg1, $arg2, $arg3, $arg4, $arg5, $arg6, $arg7, $arg8, $arg9, $arg10"

        # Run the Python script with the arguments
        python "$python_script" "$arg1" "$arg2" "$arg3" "$arg4" "$arg5" "$arg6" "$arg7" "$arg8" "$arg9" "$arg10" &

    done < <(tail -n +2 $1)

    # Wait for all background jobs to finish
    wait
}

echo "All instances of the script have been executed."
