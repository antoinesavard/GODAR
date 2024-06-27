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
python_video="../tools/plotter/video.py"
python_void="../tools/analysis/void_ratio.py"

# Get the absolute path to the parent directory of 'tools'
PROJECT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")/../" || exit
    pwd
)"

# read the first 2 lines of the input file
# Initialize variables
video=
analysis=

# Read each line of the file
while IFS=' ' read -r number rest_of_line; do
    # Check if 'number' is not empty (to skip any empty lines)
    if [[ -n $number ]]; then
        if [[ -z $video ]]; then
            video="$number"
        else
            analysis="$number"
            break
        fi
    fi
done <$1

# Print the results (optional)
echo "video.py will be run? $video"
echo "void_ratio.py will be run? $analysis"

# Read the arguments from the input file and launch the Python script
while IFS=' ' read -r arg1 arg2; do
    # store the arguments
    args_list+=("$arg1 $arg2")

    # Run the Python script with the arguments
    if [ "$video" -eq 1 ]; then
        echo "Launching $python_video with arguments $arg1 and $arg2"

        PYTHONPATH="$PROJECT_DIR" python "$python_video" "$arg1" "$arg2" &
    fi
done < <(tail -n +3 $1)

# Wait for all background jobs to finish
wait

# run the second script if prompted to do it
if [ "$analysis" -eq 1 ]; then
    echo "Performing void ratio computation"

    for args in "${args_list[@]}"; do
        # split the arguments back into individual variables
        set -- $args

        arg1=$1
        arg2=$2
        arg3="strip-collision${arg1}.mp4"

        # making sure the video file exists
        path_to_file=$(search_file "$arg3")
        if [ -e "$path_to_file" ]; then
            echo "File $arg3 found at:"
            echo "${path_to_file}"
        else
            echo "File $arg3 not found."
            exit 1
        fi
        
        # Add your additional analysis script here
        PYTHONPATH="$PROJECT_DIR" python "$python_void" "$arg1" "$arg2" "$path_to_file"
    done
fi

echo "All instances of the script have been executed."
