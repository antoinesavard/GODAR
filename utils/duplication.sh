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
# There is a companion program that creates the job
# scripts in SLURM for Alliance Can associated with the
# same experiment, as well as the input file required by
# the program to find the appropriate namelist.
#
# The program can work with a lot of different argument
# types. For example, it absolutely needs 3 arguments
# which corresponds to the basic things needed to
# duplicate a file: the namefile, and the start and end
# of the experiment numbers.
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


#---------------------------------------------------------
# From here we just check the inputs and stuff nothing
# else happens of interest
#---------------------------------------------------------
inputs=$#

if [ "${inputs}" -ne 3 ] || [ "${inputs}" -ne 7 ] ||  [ "${inputs}" -ne 10 ]; then
    echo "Number of arguments: $#"
    echo "Not enough arguments provided (need 3 for filename, and 4 per forcing up to 6 forcings)."
    echo "Reading argument file."

    # template of namelist to copy
    while true; do
	read -p "Enter the path to the template of namelist you want to copy: " template

        path_to_file=$(find $PWD -name $template)

        if [ -e "${path_to_file}" ]; then
     	    echo "${path_to_file} exists."
       	    break
	else
       	    echo "${template} does not exist."
	fi
    done
    
    # first="USER INPUT"
    while true; do
        read -p "Enter the first experiment number: " first
        if [ "${first}" -gt 99 ]; then
     	    echo "This number is above max value (99)."
	else
       	    break
	fi
    done
    
    # last="USER INPUT"
    while true; do
        read -p "Enter the last experiment number: " last
        if [ "${last}" -gt 99 ]; then
       	    echo "This number is above max value (99)."
	elif [ "${last}" -le "${first}" ]; then
       	    echo "The last number needs to be bigger than the first number."
	else
	    break
	fi
    done

elif [ "${inputs}" -e 3 ]; then
    echo "Initial namelist file: $1"
    echo "First number: $2"
    echo "Last number: $3"

    template=$1
    first=$2
    last=$3

    path_to_file=$(find $PWD -name $template)

    if [ -e "${path_to_file}" ]; then
        echo "${path_to_file} exists."
        break
    else
        echo "${template} does not exist."
	exit 1
    fi

    if [ "${first}" -gt 99 ]; then
        echo "This number is above max value (99)."
    else
       	break
    fi

    if [ "${last}" -gt 99 ]; then
        echo "This number is above max value (99)."
    elif [ "${last}" -le "${first}" ]; then
        echo "The last number needs to be bigger than the first number."
    else
	break
    fi

elif [ "${inputs}" -e 6 ]; then

    echo "uw  = $1"
    echo "vw  = $2"
    echo "ua  = $3"
    echo "va  = $4"
    echo "pfn = $5"
    echo "pfs = $6"

    uw=$1
    vw=$2
    ua=$3
    va=$4
    pfn=$5
    pfs=$6

    path_to_file=$(find $PWD -name $template)

    if [ -e "${path_to_file}" ]; then
        echo "${path_to_file} exists."
        break
    else
        echo "${template} does not exist."
	exit 1
    fi

    if [ "${first}" -gt 99 ]; then
        echo "This number is above max value (99)."
    else
       	break
    fi

    if [ "${last}" -gt 99 ]; then
        echo "This number is above max value (99)."
    elif [ "${last}" -le "${first}" ]; then
        echo "The last number needs to be bigger than the first number."
    else
	break
    fi
    
elif [ "${inputs}" -e 9 ]; then

    echo "Initial namelist file: $1"
    echo "First number: $2"
    echo "Last number: $3"

    template=$1
    first=$2
    last=$3

    path_to_file=$(find $PWD -name $template)

    if [ -e "${path_to_file}" ]; then
        echo "${path_to_file} exists."
        break
    else
        echo "${template} does not exist."
	exit 1
    fi

    if [ "${first}" -gt 99 ]; then
        echo "This number is above max value (99)."
    else
       	break
    fi

    if [ "${last}" -gt 99 ]; then
        echo "This number is above max value (99)."
    elif [ "${last}" -le "${first}" ]; then
        echo "The last number needs to be bigger than the first number."
    else
	break
    fi
    
else
    echo "Something went wrong."
    exit 1
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
exp_num=($(seq ${first} 1 ${last}))

# loop through them and create namelist files for each one
for i in ${!exp_num[@]}; do
    # create the file
    generic="namelist.nml"
    filename="${exp_num[i]}${generic}"
    cp ${template} ${filename}
    echo "Creating namelist file with name ${exp_num[i]}${generic}"

    # change the forcings here
    uw='0d0'
    vw='0d0'
    ua='0d0'
    va='0d0'

    # change the plate normal and shear forces here
    pfn='5d7'
    pfs='5d7'

    sed -i "s/^uw.*/uw = ${uw}/" ${filename}
    sed -i "s/^vw.*/vw = ${vw}/" ${filename}
    sed -i "s/^ua.*/ua = ${ua}/" ${filename}
    sed -i "s/^va.*/va = ${va}/" ${filename}

    sed -i "s/^pfn.*/pfn = ${pfn}/" ${filename}
    sed -i "s/^pfs.*/pfs = ${pfs}/" ${filename}

    sed -i "s/^Xfile.*/Xfile = x${exp_num[i]}.dat/" ${filename}
    sed -i "s/^Yfile.*/Yfile = y${exp_num[i]}.dat/" ${filename}
    sed -i "s/^Rfile.*/Rfile = r${exp_num[i]}.dat/" ${filename}
    sed -i "s/^Hfile.*/Hfile = h${exp_num[i]}.dat/" ${filename}
    sed -i "s/^Tfile.*/Tfile = t${exp_num[i]}.dat/" ${filename}
    sed -i "s/^Ofile.*/Ofile = o${exp_num[i]}.dat/" ${filename}
done
