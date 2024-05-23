mkdir output
mkdir files
mkdir namelist
mkdir inputs

cp generic/SConstruct_generic SConstruct
cp generic/namelist_generic.nml namelist/namelist.nml
cp generic/input_generic inputs/input
cp generic/input_restart_generic inputs/input_restart
cp generic/start_generic.sh start.sh