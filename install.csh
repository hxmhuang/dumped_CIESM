#!/bin/csh -f 
################################################
# Install ciesm on your own path
###############################################
set ciesmroot=$cwd
echo "Your ciesm root is ${ciesmroot}"

set casedir="${ciesmroot}/ciesm.run"
set archdir="${ciesmroot}/ciesm.archive"

if !(-d ${casedir}) then 
	mkdir -p ${casedir} 
endif
if !(-d ${archdir}) then 
	mkdir -p ${archdir}
endif


set RUNDIR="${casedir}/"'$'"CASE/run"
set EXEROOT="${casedir}/"'$'"CASE/bld"
set DOUT_S_ROOT="${archdir}/"'$'"CASE"

echo "RUNDIR:" $RUNDIR
echo "EXEROOT:" $EXEROOT
echo "DOUT_S_ROOT:" $DOUT_S_ROOT

cd ${ciesmroot}/ciesm.model/scripts/ccsm_utils/Machines
sed -i "37c <RUNDIR>${RUNDIR}</RUNDIR><!--- complete path to the run directory -->" config_machines.xml
sed -i "38c <EXEROOT>${EXEROOT}</EXEROOT> <!--- complete path to the build directory -->" config_machines.xml
sed -i "41c <DOUT_S_ROOT>${DOUT_S_ROOT}</DOUT_S_ROOT><!--- complete path to a short term archiving directory -->" config_machines.xml

sed -i "54c <RUNDIR>${RUNDIR}</RUNDIR><!--- complete path to the run directory -->" config_machines.xml
sed -i "55c <EXEROOT>${EXEROOT}</EXEROOT> <!--- complete path to the build directory -->" config_machines.xml
sed -i "58c <DOUT_S_ROOT>${DOUT_S_ROOT}</DOUT_S_ROOT><!--- complete path to a short term archiving directory -->" config_machines.xml

echo "The ciesm have successfully installed on ${ciesmroot}"
