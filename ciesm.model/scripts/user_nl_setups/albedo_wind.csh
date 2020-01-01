#! /bin/csh -f

if ($albedo_wind == on) then

cat >> $case1/user_nl_cpl << EOF 
!The albedo-wind scheme
albedo_wind = .true.
ocean_tight_coupling = .true.

EOF
    
echo ">>>>ALBEDO-WIND SCHEME                ON<<<<"

else 
   
echo ">>>>ALBEDO-WIND SCHEME               OFF<<<<"

endif

