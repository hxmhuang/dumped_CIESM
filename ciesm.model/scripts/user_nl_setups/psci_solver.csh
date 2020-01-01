#! /bin/csh -f

if ($PSCI_SOLVER == on) then

cat >> $case1/user_nl_pop2 << EOF 
!Namelist setups for PSCI SOLVER
solverchoice = 'PCSI'
preconditionerchoice = 'evp'
diag_global_freq_opt = 'nday'
convergencecheckstart = 200
convergencecheckfreq = 20

EOF
    
echo ">>>>PSCI SOLVER                       ON<<<<"

else 
   
echo ">>>>PSCI SOLVER                      OFF<<<<"

endif

