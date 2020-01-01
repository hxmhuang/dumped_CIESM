#! /bin/csh -f

if ($sub_oro_drag == on) then

cat >> $case1/user_nl_cam << EOF 
!Namelist setups for subgrid orographic drag scheme
oro_drag = .true.
do_tms   = .false.

EOF
    
echo ">>>>SUBGRID OROGRAPHIC DRAG SCHEME    ON<<<<"

else 
   
echo ">>>>SUBGRID OROGRAPHIC DRAG SCHEME   OFF<<<<"

endif

