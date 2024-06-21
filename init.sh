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
cp generic/init_args_generic utils/init_args.dat
cp generic/video_args_generic.dat utils/video_args.dat
