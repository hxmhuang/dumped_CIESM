#! /bin/csh -f

if ($shr_flux == on) then

cp -r $rootpwd/SourceMods.POP/CESM-FLUX/* $case1/SourceMods/
    
echo ">>>>SHR_FLUX SCHEME                   ON<<<<"

else 
   
echo ">>>>SHR_FLUX SCHEME                  OFF<<<<"

endif

