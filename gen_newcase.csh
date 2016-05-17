#!/bin/csh -f 
################################################
# For POP users:
#   Pls define grid, compset, pp, number, casename, SourceMods,
#     stop_option, stop_n by yourself
###############################################


####   USER Define ##############################
# resolution
set grid="ne30_g16_rx1"
#set grid="ne30_g16"
#set grid="f09_g16"
# compset
set compset="G"
# experimental years
set stop_option="nmonths"
set stop_n=2
# process number
set pp=100
# machine name
set machine="cessb"

#set iopp=`expr ${pp} / 8`
set iopp=0
set totalpp=`expr ${pp} + ${iopp}`
set node=`expr ${totalpp} / 20`

# case name
set casename="mycase2_${compset}_${grid}_${totalpp}" #@@@@@PLEASE MODIFY CASE NAME!!!!@@@@@@

#set SourceMods="SourceMods.B2000CAM5.fv"
#set SourceMods="SourceMods.B2000CAM5.ne"
#set SourceMods="SourceMods.F2000CAM5.fv"
#set SourceMods="SourceMods.F2000CAM5.ne"
#set SourceMods="SourceMods.G_NORMAL_YEAR.ne"

# set true to save the restart files,always set it to true
set dout="TRUE"

#######  USER DEFINE End ############################# 
cd CIESM.MODEL

set rootpwd=$cwd/scripts
echo "Your script path is ${rootpwd}"

cd ..
set ciesmroot=$cwd
echo "Your cesm root is ${ciesmroot}"

set casedir="${ciesmroot}/CIESM.RUN"
set archdir="${ciesmroot}/CIESM.ARCHIVE"
if !(-d ${casedir}) then 
	mkdir -p ${casedir} 
	touch .gitignore
endif
if !(-d ${archdir}) then 
	mkdir -p ${archdir}
	touch .gitignore
endif

set case1="${casedir}/${casename}" #@@@@@PLEASE MODIFY CASE NAME!!!!@@@@@@

# create the folder
set RUNDIR="${casedir}/"'$'"CASE/run"
set EXEROOT="${casedir}/"'$'"CASE/bld"
set DOUT_S_ROOT="${ciesmroot}/CIESM.ARCHIVE/"'$'"CASE"

cd $rootpwd
sed -i "500c <RUNDIR>${RUNDIR}</RUNDIR><!--- complete path to the run directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "501c <EXEROOT>${EXEROOT}</EXEROOT> <!--- complete path to the build directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "503c <DOUT_S_ROOT>${DOUT_S_ROOT}</DOUT_S_ROOT><!--- complete path to a short term archiving directory -->" ccsm_utils/Machines/config_machines.xml

sed -i "517c <RUNDIR>${RUNDIR}</RUNDIR><!--- complete path to the run directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "518c <EXEROOT>${EXEROOT}</EXEROOT> <!--- complete path to the build directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "520c <DOUT_S_ROOT>${DOUT_S_ROOT}</DOUT_S_ROOT><!--- complete path to a short term archiving directory -->" ccsm_utils/Machines/config_machines.xml

sed -i "534c <RUNDIR>${RUNDIR}</RUNDIR><!--- complete path to the run directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "535c <EXEROOT>${EXEROOT}</EXEROOT> <!--- complete path to the build directory -->" ccsm_utils/Machines/config_machines.xml
sed -i "537c <DOUT_S_ROOT>${DOUT_S_ROOT}</DOUT_S_ROOT><!--- complete path to a short term archiving directory -->" ccsm_utils/Machines/config_machines.xml


echo "You will create a case in $case1"
./create_newcase -compset $compset -res $grid -mach $machine -case $case1

########  MODIFY CONFIG FILES ######################
cd $case1
source $case1/Tools/ccsm_getenv

# set up the # of processes, don't touch it
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $pp
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $pp

./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val 0

./xmlchange -file env_run.xml -id STOP_OPTION -val $stop_option
./xmlchange -file env_run.xml -id STOP_N -val $stop_n
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val $dout
#./xmlchange -file env_run.xml -id PIO_TYPENAME -val "cfio" 
#./xmlchange -file env_run.xml -id PIO_TYPENAME -val "pnetcdf" 

############CFIO config#####################
#./xmlchange -file env_run.xml -id CFIO_USED -val "TRUE"
#./xmlchange -file env_run.xml -id CFIO_OCNIO_NUM -val $iopp
#./xmlchange -file env_run.xml -id PIO_STRIDE -val 1 
#./xmlchange -file env_run.xml -id PIO_ROOT -val 0 
#./xmlchange -file env_run.xml -id PIO_NUMTASKS -val $pp 
 
########  MODIFY CONFIG FILES END ###########


######## CESM_STEP  #########################
cd $case1
./cesm_setup

#cp $rootpwd/Macros.cess Macros
#cp $rootpwd/Makefile.cess Tools/Makefile

#cp -r $rootpwd/$SourceMods/* SourceMods/
#cp -r $rootpwd/SourceMods.ACC/* SourceMods/
#cp -r $rootpwd/Buildconf/* Buildconf/

#add the config file for models
#cp -r $rootpwd/user_nl_* .


#add the optimization of gather_scatter
#cp -r $rootpwd/SourceMods.gather/* SourceMods/


#=============================================
#add the PCSI solver, to replace the ChronGear solver
cp -r $rootpwd/user_nl_pop2 .
cp -r $rootpwd/SourceMods.POP/CESM-PCSI/* SourceMods/

#add the shr_flux correction
cp -r $rootpwd/SourceMods.POP/CESM-FLUX/* SourceMods/

#add the tke correction of vmix_kpp
#cp -r $rootpwd/SourceMods.POP/CESM-TKE/* SourceMods/
#=============================================


#add the CFIO optimization to the PIO
#cp -r $rootpwd/SourceMods.POP/CESM-CFIO/* SourceMods/

#cp -r $rootpwd/SourceMods.CAM5/CESM-LIN/* SourceMods/src.cam/

######## CESM_STEP END ######################

######## CESM BUILD #########################
./$casename.build
######## CESM BUILD END #####################
