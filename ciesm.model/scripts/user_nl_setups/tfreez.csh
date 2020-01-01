#! /bin/csh -f

if ($tfreez == on) then

cat >> $case1/user_nl_pop2 << EOF 
!ice freeze temperature scheme
tfreez_salt = .true.
EOF
    
cat >> $case1/user_nl_cice << EOF 
!ice freeze temperature scheme
tfreez_salt_ice = .true.
EOF

echo ">>>>FREEZE TEMPERATURE SCHEME         ON<<<<"

else 
   
echo ">>>>FREEZE TEMPERATURE SCHEME        OFF<<<<"

endif

