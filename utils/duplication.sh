#!/bin/bash

#---------------------------------------------------------
# Program duplication.sh
#
# This program takes a given namelist template for GoDAR
# and creates copies of it with modified values inside
# which were provided by the user. The values that can be
# modified are: the forcings, the plate forcings, and the
# input file to use.
#
# It does the same with the input files necessary for
# godar to be launched. It creates files with the correct
# experience number and the appropriate number of cores
# and pointing towards the proper namelist file.
#
# There is a companion program that creates the job
# scripts in SLURM for Alliance Can associated with the
# same experiments.
#
# The program can work with a lot of different argument
# types. For example, it absolutely needs 4 arguments
# which corresponds to the basic things needed to
# duplicate a file: the namefile, and the start and end
# of the experiment numbers and the number of cores to
# use in openmp.
#
# You can also provide the program with an forcing
# variable and 3 parameters. E.g. uw 0 20 1, which would
# create 21 namelists, 1 for each case of the 0 to 20 uw
# by increment of 1. And a total of 6 forcings can be
# passed to this program, such that a lot of files can be
# created at once (uw, vw, ua, va, pfn, pfs).
#
# Note however that the forcings are not independant of
# each other and specifying for example uw and vw will
# create an array of namelist for all combinations.

# Written by Antoine Savard, 2024
#---------------------------------------------------------

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

get_cores() {
    if command -v nproc &>/dev/null; then
        # Linux
        num_cores=$(nproc)
    elif [[ "$(uname)" == "Darwin" ]]; then
        # macOS
        num_cores=$(sysctl -n hw.ncpu)
    elif command -v getconf &>/dev/null; then
        # POSIX-compliant
        num_cores=$(getconf _NPROCESSORS_ONLN)
    else
        # Fallback
        num_cores="Unknown"
    fi
    echo "Total number of CPU cores: $num_cores"
}

# function to prompt the user for arguments
prompt_for_inputs() {
    # template of namelist to copy
    while true; do
        read -rp "Enter the name of the template namelist you want to copy: " template
        echo "You provided ${template} as a namefile."

        path_to_file=$(search_file "${template}")
        if [ -e "$path_to_file" ]; then
            echo "File ${template} found at:"
            echo "${path_to_file}"
            break
        else
            echo "File ${template} not found."
        fi
    done
    echo ""

    # restart="USER INPUT"
    while true; do
        read -rp "Should we use restart template (0/1): " restart
	echo "You said ${restart} restart"
        if [ "${restart}" -ne 1 ] && [ "${restart}" -ne 0 ]; then
            echo "This should be 0 or 1"
        else
            break
        fi
    done
    
    # first="USER INPUT"
    while true; do
        read -rp "Enter the first experiment number: " first
	echo "You provided ${first} as a starting point."
        if [ "${first}" -gt 99 ]; then
	    echo "This number is above max value (99)."
        else
            break
        fi
    done

    # last="USER INPUT"
    while true; do
        read -rp "Enter the last experiment number: " last
        echo "You provided ${last} as an end point."
        if [ "${last}" -gt 99 ]; then
            echo "This number is above max value (99)."
        elif [ "${last}" -le "${first}" ]; then
            echo "The last number needs to be bigger than the first number."
        else
            break
        fi
    done
    echo ""

    # cores="USER INPUT"
    while true; do
        read -rp "Enter the number of openmp cores to use: " cores
        echo "You want to use ${cores} cores."
        get_cores
        if [ "${cores}" -gt "${num_cores}" ]; then
            echo "Your machine does not have enough cores (${num_cores})."
            exit 1
        else
            break
        fi
    done
    echo ""
}

forcings() {
    local list
    local input

    # forcings to modify
    while true; do
        read -rp "Enter the name of the variable to modify: " input
        list=($input)
        var="${list[0]}"
        fstart="${list[1]}"
        fend="${list[2]}"
        fstep="${list[3]}"

        if [ -z "${var}" ]; then
            break
        elif [ "${var}" = "/" ]; then
            echo "End of forcings to read."
            echo ""
            continue
        elif [ "${var}" != "uw" ] && [ "${var}" != "vw" ] && [ "${var}" != "ua" ] && [ "${var}" != "va" ] && [ "${var}" != "pfn" ] && [ "${var}" != "pfs" ]; then
            echo "The variable you provided cannot be modified either because it does not exist or there is a typo."
            echo ""
            continue
        fi
        while true; do
            echo "You provided ${var} as variable to modify."
            echo "The inputs are:" "${fstart}" "${fend}" "${fstep}"
            if [ -z "${fstart}" ] || [ -z "${fend}" ]; then
                echo "You need to provide start point, end point and step size."
                echo "Or value and number of copies"
                break
            elif [ "${var}" = "uw" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C uw=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    uw=()
                    for ((i = 0; i < fend; i++)); do
                        uw+=(${fstart})
                    done
                fi
                echo "The list to use is uw =" "${uw[@]}"
                break
            elif [ "${var}" = "vw" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C vw=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    vw=()
                    for ((i = 0; i < fend; i++)); do
                        vw+=(${fstart})
                    done
                fi
                echo "The list to use is vw =" "${vw[@]}"
                break
            elif [ "${var}" = "ua" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C ua=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    ua=()
                    for ((i = 0; i < fend; i++)); do
                        ua+=(${fstart})
                    done
                fi
                echo "The list to use is ua =" "${ua[@]}"
                break
            elif [ "${var}" = "va" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C va=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    va=()
                    for ((i = 0; i < fend; i++)); do
                        va+=(${fstart})
                    done
                fi
                echo "The list to use is va =" "${va[@]}"
                break
            elif [ "${var}" = "pfn" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C pfn=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    pfn=()
                    for ((i = 0; i < fend; i++)); do
                        pfn+=(${fstart})
                    done
                fi
                echo "The list to use is pfn =" "${pfn[@]}"
                break
            elif [ "${var}" = "pfs" ]; then
                if [ -n "${fstep}" ]; then
                    LANG=C pfs=($(seq "${fstart}" "${fstep}" "${fend}"))
                elif [ -z "${fstep}" ]; then
                    pfs=()
                    for ((i = 0; i < fend; i++)); do
                        pfs+=(${fstart})
                    done
                fi
                echo "The list to use is pfs =" "${pfs[@]}"
                break
            else
                break
            fi
        done
        echo ""
    done
}

count() {
    local count
    local array
    total_count=1

    # Iterate through each array name
    for array in "$@"; do
        # verify array is non empty
        if [ -z "${array}" ]; then
            continue
        fi
        # Count the number of elements in the current array
        count=$(eval echo \${#${array}[@]})

        # Add to total count
        if [ "${count}" -ne 0 ]; then
            total_count=$((total_count * count))
        fi

        # Print the count for the current array
        echo "The array ${array} has ${count} elements."
    done
    echo "The total number of namelists needed is" "${total_count}"
    echo ""
}

check_forcing() {
    if [ -z "${uw}" ]; then
        uw=0
    fi
    if [ -z "${vw}" ]; then
        vw=0
    fi
    if [ -z "${ua}" ]; then
        ua=0
    fi
    if [ -z "${va}" ]; then
        va=0
    fi
    if [ -z "${pfn}" ]; then
        pfn=0
    fi
    if [ -z "${pfs}" ]; then
        pfs=0
    fi
}

#---------------------------------------------------------
# From here we just check the inputs and stuff nothing
# else happens of interest
#---------------------------------------------------------
inputs=$#

if [[ "${inputs}" -ne 0 ]]; then
    echo "Number of arguments: $# > 0"
    echo "This function will only work with arguments passed with '< args.txt'"
    exit 1

else
    echo "I will try to read an input file or ask for inputs."
    echo ""
    prompt_for_inputs
    forcings
    check_forcing
fi

#---------------------------------------------------------
# From here we are actually creating the files, and
# for each file we are changing the values with the
# corresponding values given by the user
#
# Things that can and need to change:
#  -forcings
#  -plate forces
#  -initial input files
#---------------------------------------------------------

# create an array of all experiment numbers
exp_num=($(seq "${first}" 1 "${last}"))

# create an array of all arrays for forcings
force=("uw" "vw" "ua" "va" "pfn" "pfs")
count "${force[@]}"

if [ "${#exp_num[@]}" -ne "${total_count}" ]; then
    echo "The number of experiment provided was not accurate."
    last=$((total_count + first - 1))
    exp_num=($(seq "${first}" 1 "${last}"))
    echo "It now has been adjusted to ${#exp_num[@]}"
    echo ""
fi

# loop through the exp_num and create namelist files for each one
for i in "${!exp_num[@]}"; do
    # create the file
    generic="namelist.nml"
    filename=../namelist/"${exp_num[i]}${generic}"
    cp "${path_to_file}" "${filename}"
    echo "Creating namelist file with name ${exp_num[i]}${generic}"

    if [[ "$(uname)" == "Linux" ]]; then
        # Linux
        sed -i "s/^    Xfile.*/    Xfile = \"files\/x${exp_num[i]}.dat\"/" "${filename}"
        sed -i "s/^    Yfile.*/    Yfile = \"files\/y${exp_num[i]}.dat\"/" "${filename}"
        sed -i "s/^    Rfile.*/    Rfile = \"files\/r${exp_num[i]}.dat\"/" "${filename}"
        sed -i "s/^    Hfile.*/    Hfile = \"files\/h${exp_num[i]}.dat\"/" "${filename}"
        sed -i "s/^    Tfile.*/    Tfile = \"files\/theta${exp_num[i]}.dat\"/" "${filename}"
        sed -i "s/^    Ofile.*/    Ofile = \"files\/omega${exp_num[i]}.dat\"/" "${filename}"
    elif [[ "$(uname)" == "Darwin" ]]; then
        # macOS
        sed -i "" "s/^    Xfile.*/    Xfile = \"files\/x${exp_num[i]}.dat\"/" "${filename}"
        sed -i "" "s/^    Yfile.*/    Yfile = \"files\/y${exp_num[i]}.dat\"/" "${filename}"
        sed -i "" "s/^    Rfile.*/    Rfile = \"files\/r${exp_num[i]}.dat\"/" "${filename}"
        sed -i "" "s/^    Hfile.*/    Hfile = \"files\/h${exp_num[i]}.dat\"/" "${filename}"
        sed -i "" "s/^    Tfile.*/    Tfile = \"files\/theta${exp_num[i]}.dat\"/" "${filename}"
        sed -i "" "s/^    Ofile.*/    Ofile = \"files\/omega${exp_num[i]}.dat\"/" "${filename}"
    else
        echo "I don't know what to do with this..."
    fi

    if [ "${restart}" -eq 0 ]; then
	# create the associate input file
	filename=../inputs/"${exp_num[i]}input"

	# Create the file and write the contents
	cat <<EOL >"${filename}"
1                   input namelist
${exp_num[i]}${generic}      namelist name
0                   input restart
${exp_num[i]}                  exp version
${cores}                  number of threads
EOL

	# Print a message indicating the file was created
	echo "Input file ${exp_num[i]}input has been created."

    elif [ "${restart}" -eq 1 ]; then
	# create the associate input file
	filename=../inputs/"${exp_num[i]}input_restart"

	# Create the file and write the contents
	cat <<EOL >"${filename}"
1                   input namelist
${exp_num[i]}${generic}      namelist name
1                   input restart
$((${exp_num[i]} - ${total_count}))                  restart exp
1000                restart time
${exp_num[i]}                  exp version
${cores}                  number of threads
EOL

	# Print a message indicating the file was created
	echo "Input file ${exp_num[i]}input_restart has been created."
    fi
done

# loop through the file and modify the values at the correct places
itnum=${first}
for i in "${!uw[@]}"; do
    for j in "${!vw[@]}"; do
        for k in "${!ua[@]}"; do
            for l in "${!va[@]}"; do
                for m in "${!pfn[@]}"; do
                    for o in "${!pfs[@]}"; do

                        generic="namelist.nml"
                        filename="../namelist/${itnum}${generic}"
                        if [[ "$(uname)" == "Linux" ]]; then
                            # Linux
                            # change the forcings here
                            LANG=C sed -i "s/^    uw.*/    uw = ${uw[i]}/" "${filename}"
                            LANG=C sed -i "s/^    vw.*/    vw = ${vw[j]}/" "${filename}"
                            LANG=C sed -i "s/^    ua.*/    ua = ${ua[k]}/" "${filename}"
                            LANG=C sed -i "s/^    va.*/    va = ${va[l]}/" "${filename}"

                            # change the plate normal and shear forces here
                            LANG=C sed -i "s/^    pfn.*/    pfn = ${pfn[m]}/" "${filename}"
                            LANG=C sed -i "s/^    pfs.*/    pfs = ${pfs[o]}/" "${filename}"

                        elif [[ "$(uname)" == "Darwin" ]]; then
                            # macOS
                            # change the forcings here
                            LANG=C sed -i "" "s/^    uw.*/    uw = ${uw[i]}/" "${filename}"
                            LANG=C sed -i "" "s/^    vw.*/    vw = ${vw[j]}/" "${filename}"
                            LANG=C sed -i "" "s/^    ua.*/    ua = ${ua[k]}/" "${filename}"
                            LANG=C sed -i "" "s/^    va.*/    va = ${va[l]}/" "${filename}"

                            # change the plate normal and shear forces here
                            LANG=C sed -i "" "s/^    pfn.*/    pfn = ${pfn[m]}/" "${filename}"
                            LANG=C sed -i "" "s/^    pfs.*/    pfs = ${pfs[o]}/" "${filename}"

                        else
                            echo "I don't know that to do..."
                        fi

                        itnum=$((itnum + 1))
                    done
                done
            done
        done
    done
done

# creation of the slurm sbatch script
if [ "${restart}" -eq 0 ]; then
    cat <<EOL >"../jobs/${first}${last}init_plate.sh"
#!/bin/bash

#SBATCH --array=${first}-${last}%$((${last} - ${first} + 1))
#SBATCH --time=3-0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${cores}
#SBATCH --mem=4G
#SBATCH --output=../output/out.%a

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

cd ..

srun --cpus-per-task=\$SLURM_CPUS_PER_TASK ./godar < inputs/\${SLURM_ARRAY_TASK_ID}input

EOL

elif [ "${restart}" -eq 1 ]; then
    cat <<EOL >"../jobs/${first}${last}init_plate.sh"
#!/bin/bash

#SBATCH --array=${first}-${last}%$((${last} - ${first} + 1))
#SBATCH --time=3-0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${cores}
#SBATCH --mem=4G
#SBATCH --output=../output/out.%a

export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

cd ..

srun --cpus-per-task=\$SLURM_CPUS_PER_TASK ./godar < inputs/\${SLURM_ARRAY_TASK_ID}input_restart

EOL
    
fi
