#!/bin/bash

mkdir output
mkdir files
mkdir namelist
mkdir inputs
mkdir jobs

cp generic/SConstruct_generic SConstruct
cp generic/namelist_generic.nml namelist/namelist.nml
cp generic/input_generic inputs/input
cp generic/input_restart_generic inputs/input_restart
cp generic/start_generic.sh jobs/start.sh
cp generic/args_generic.dat namelist/args.dat
cp generic/input_init_generic plots/input_init.dat
