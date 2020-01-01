#!/bin/csh -f 
################################################
# For ciesm user on Sunway supercomputer:
#   Pls define grid, compset, pp, number, casename, SourceMods,
#     stop_option, stop_n by yourself
###############################################

####   USER Define ##############################
# resolution
 set grid="ne30_g16"
# set grid="ne30_s05"
# set grid="f09_s05"
#  set grid="f09_g16"
# set grid="f09_f09"
# set grid="f09_g16_rx1"

# mach
# set mach = "sunway"
set mach = "sunwayx86"

#compset
# set compset = "G"
#set compset="F_AMIP_CAM5"
#set compset="F_2000_CAM5"
set compset="B_1850_CAM5"
#set compset="B20TRC5"
#set compset = "GIAF"

#experimental years
set stop_option="nyears"
set stop_n=150

setenv mycompset $compset

#set resubmit
set resubmit = 0

#process number
 set pp=960
set iopp=0
set totalpp=`expr ${pp} + ${iopp}`
set node=`expr ${totalpp} / 4`

#set CLM_FORCE_COLDSTART
#if s05 grid, set CLM_FORCE_COLDSTART on
#  set clm_force_coldstart = on
  set clm_force_coldstart = off

# case name
set casename="B1850_CAM5" #@@@@@PLEASE MODIFY CASE NAME!!!!@@@@@@

# set true to save the restart files,always set it to true
set dout="TRUE"

#######  USER DEFINE End ############################## 

cd ciesm.model
setenv rootmod $cwd/models
echo "Your models path is ${rootmod}"
setenv rootpwd $cwd/scripts
echo "Your script path is ${rootpwd}"

cd ..
#set ciesmroot=$cwd
setenv  ciesmroot $cwd
echo "Your ciesm root is ${ciesmroot}"

setenv case1 "${ciesmroot}/ciesm.run/${casename}" #@@@@@PLEASE MODIFY CASE NAME!!!!@@@@@@

#echo $RUNDIR
cd $rootpwd
echo "You will create a case in $case1"
./create_newcase -compset $compset -res $grid -mach $mach -case $case1
#./create_newcase -user_compset $compset -res $grid -mach $mach -case $case1

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
./xmlchange -file env_run.xml -id RESUBMIT -val $resubmit
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val $dout
./xmlchange -file env_run.xml -id CLM_FORCE_COLDSTART -val $clm_force_coldstart
./xmlchange -file env_run.xml -id PIO_TYPENAME -val "pnetcdf" 
./xmlchange -file env_run.xml -id BFBFLAG  -val "TRUE"

########  MODIFY CONFIG FILES END ###########


######## CESM_STEP  #########################
cd $case1

cp $rootpwd/Makefile Tools/
 
#add Modified ZM, Single-ice Micro-physics, PDF Macro-physics, Four Stream schemes
#  setenv atm_schemes  "off"
  setenv atm_schemes  "on"

#add the PCSI solver, to replace the ChronGear solver
#  setenv PCSI_SOLVER  "off"
  setenv PCSI_SOLVER  "on"

#add the albedo-wind scheme
#  setenv albedo_wind  "off"
 setenv albedo_wind  "on"

#add the shr_flux correction
#  setenv shr_flux  "off"
 setenv shr_flux  "on"

#add the GSDE land dataset  
#  setenv GSDE "off"
  setenv GSDE "on"

#add the lateral melting sea-ice shceme  
#  setenv LM_ice "off"
  setenv LM_ice "on"

########CESM_SETUP####i######################
cd $case1

echo $case1

./cesm_setup
touch $rootpwd/../models/utils/pio/*.in

########SET USER_NL FILES####################
cd $rootpwd
./user_nl_setups.csh

######## CESM BUILD #########################
cd $case1
./$casename.build
######## CESM BUILD END #####################
