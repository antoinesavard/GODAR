#!/bin/bash

mkdir output
mkdir files
mkdir namelist
mkdir inputs
mkdir jobs
mkdir plots
mkdir plots/anim
mkdir plots/plot

cp generic/SConstruct_generic SConstruct
cp generic/namelist_generic.nml namelist/namelist.nml
cp generic/input_generic inputs/input
cp generic/input_restart_generic inputs/input_restart
cp generic/start_generic.sh jobs/start.sh
cp generic/args_generic.dat utils/args.dat
cp generic/args_restart_generic.dat utils/args_restart.dat
cp generic/init_args_generic.dat utils/init_args.dat
cp generic/video_args_generic.dat utils/video_args.dat

# Set executable permissions for all .sh files in the current directory and subdirectories
find . -type f -name "*.sh" -exec chmod +x {} \;

scons=SConstruct
start=jobs/start.sh

# modifies what needs to be modified
if [[ "$(uname)" == "Linux" ]]; then
    # Linux
    echo "This is some form of Linux, we will deal with it accordingly"
    echo "If this is a supercomputer, it probably won't find coretran by itself, and you might want to do this manually."

    # checks possible places where coretran could be installed
    echo "Looking for coretran installation"

    while true; do
        # Prompt the user for the directory name
        read -rp "Enter the directory name to search from: " dir_name

        # Search for the directory coretran
        core_lib=$(find "$dir_name" \( -path /dev -o -path /home -o -path /sys -o -path /etc -o -path /mnt -o -path /root -o -path /tmp -o -path /var \) -prune -o \( -type d \( -path "*/coretran/lib" \) \) -print 2>/dev/null)

        core_inc=$(find "$dir_name" \( -path /dev -o -path /home -o -path /sys -o -path /etc -o -path /mnt -o -path /root -o -path /tmp -o -path /var \) -prune -o \( -type d \( -path "*/coretran/include/coretran" \) \) -print 2>/dev/null)

        # Check if any directories were found
        if [[ -n "$core_lib" && -n "$core_inc" ]]; then
            echo "Found coretran installation directories:"
            echo "$core_lib"
            echo "$core_inc"
            break
        else
            echo "No coretran installation found."
            echo "Please try again"
        fi
    done

    echo "Modifying SConstruct to match coretran install path"

    sed -i "s|^env\.Append(LIBPATH=\"/path.*|env.Append(LIBPATH=\"${core_lib}\")|" "${scons}"

    sed -i "s|^env\.Append(F90PATH=\"/path.*|env.Append(F90PATH=\"${core_inc}\")|" "${scons}"

    echo "Done"

    echo "Adding coretran install path to PATH"

    echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${core_lib}" >>~/.bashrc
    exec bash

elif [[ "$(uname)" == "Darwin" ]]; then
    # macOS
    echo "This is macOS, we will deal with it accordingly"

    # checks possible places where coretran could be installed
    echo "Looking for coretran installation"

    core_lib=$(find / \( -path /System -o -path /Library -o -path /private -o -path $HOME/Library \) -prune -o \( -type d \( -path "*/coretran/lib" \) \) -print 2>/dev/null)

    core_inc=$(find / \( -path /System -o -path /Library -o -path /private -o -path $HOME/Library \) -prune -o \( -type d \( -path "*/coretran/include" \) \) -print 2>/dev/null)

    # Check if any directories were found
    if [[ -n "$core_lib" && -n "$core_inc" ]]; then
        echo "Found coretran installation directories:"
        echo "$core_lib"
        echo "$core_inc"

    else
        echo "No coretran installation found."
        exit 1
    fi

    echo "Modifying SConstruct to match coretran install path"

    sed -i "" "s|^env\.Append(LIBPATH=\"/path.*|env.Append(LIBPATH=\"${core_lib}\")|" "${scons}"

    sed -i "" "s|^env\.Append(F90PATH=\"/path.*|env.Append(F90PATH=\"${core_inc}\")|" "${scons}"

    echo "Done"

    # changing the option that is not the same on macOS
    echo "Changing the -mcmodel option to 'small'"

    sed -i "" "s|^        \"-mcmodel=.*|        \"-mcmodel=small\",|" "${scons}"

    # changing jobs/start.sh to accomodate for the dylib of coretran
    echo "Changing ${start} to accomodate for the dylib of coretran"

    sed -i '' '4i\
install_name_tool -change @rpath/libcoretran.dylib '"${core_lib}"'/libcoretran.dylib godar' "${start}"

else
    # Fallback
    echo "Why are you using Windows? This is a crime."
fi
