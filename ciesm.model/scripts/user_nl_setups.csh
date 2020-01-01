#! /bin/csh -f

echo ">>>>====================================<<<<"
cat >> $case1/README.case << EOF

>>>>====================================<<<<
EOF

#############PCSI_SOLVER#################
if ($PCSI_SOLVER == on) then

cat >> $case1/user_nl_pop2 << EOF 
!Namelist setups for PCSI SOLVER
solverchoice = 'PCSI'
preconditionerchoice = 'evp'
diag_global_freq_opt = 'nday'
convergencecheckstart = 200
convergencecheckfreq = 20

EOF
    
echo ">>>>PCSI SOLVER                       ON<<<<"
cat >> $case1/README.case << EOF
>>>>PCSI SOLVER                       ON<<<<
EOF

else 
   
echo ">>>>PCSI SOLVER                      OFF<<<<"
cat >> $case1/README.case << EOF
>>>>PCSI SOLVER                      OFF<<<<
EOF

endif

#############ALBEDO-WIND#################
if ($albedo_wind == on) then

cat >> $case1/user_nl_cpl << EOF 
!The albedo-wind scheme
albedo_wind = .true.

EOF
    
echo ">>>>ALBEDO-WIND SCHEME                ON<<<<"
cat >> $case1/README.case << EOF
>>>>ALBEDO-WIND SCHEME                ON<<<<
EOF

else 
   
echo ">>>>ALBEDO-WIND SCHEME               OFF<<<<"
cat >> $case1/README.case << EOF
>>>>ALBEDO-WIND SCHEME               OFF<<<<
EOF

endif

cd
#############SHR_FLUX#################
if ($shr_flux == on) then

cp -r $rootpwd/SourceMods.POP/CESM-FLUX/* $case1/SourceMods/
    
echo ">>>>SHR_FLUX SCHEME                   ON<<<<"
cat >> $case1/README.case << EOF
>>>>SHR_FLUX SCHEME                   ON<<<<
EOF

else 
   
echo ">>>>SHR_FLUX SCHEME                  OFF<<<<"
cat >> $case1/README.case << EOF
>>>>SHR_FLUX SCHEME                  OFF<<<<
EOF

endif




#########GSDE########

if( $GSDE == on) then
cd $case1

cp -r $rootpwd/SourceMods.CLM/imporved_clm40/*  $case1/SourceMods/src.clm

echo ">>>>GSDE land dataset                 ON<<<<"
cat >> $case1/README.case << EOF
>>>>GSDE land dataset                 ON<<<<
EOF

else

echo ">>>>GSDE land dataset                OFF<<<<"
cat >> $case1/README.case << EOF
>>>>GSDE land dataset                OFF<<<<
EOF
endif

#####lateral melting sea-ice shceme#####
if ($LM_ice == on) then

cp -r $rootpwd/SourceMods.CICE/src.cice/*   $case1/SourceMods/src.cice

echo ">>>>LATERAL MELTING SEA-ICE SCHEME    ON<<<<"
cat >> $case1/README.case << EOF
>>>>LATERAL MELTING SEA-ICE SCHEME    ON<<<<
EOF

else

echo ">>>>LATERAL MELTING SEA-ICE SCHEME   OFF<<<<"
cat >> $case1/README.case << EOF
>>>>LATERAL MELTING SEA-ICE SCHEME   OFF<<<<
EOF

endif


#####Modified ZM, Single-ice Micro-physics, PDF Macro-physics, Four Stream schemes#####
if ($atm_schemes == on) then

cp -r $rootpwd/SourceMods.CAM5/src.cam/*   $case1/SourceMods/src.cam
cd $case1
./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -val "-phys cam5 -chem none"

cp -r ${rootmod}/atm/cam/bld/namelist_files/atm.namelist_definition.xml ${rootmod}/atm/cam/bld/namelist_files/namelist_definition.xml

#Namelist setups
cat >> $case1/user_nl_cam << EOF
zmconv_microp = .true.
zmconv_stc    = .true.
zmconv_ke = 1.0E-5

history_aerosol				=.true.
history_aero_optics			=.true.

cldfrc_dp1 = 0.05D0
cldfrc_dp2 = 200.0D0
cldfrc_rhminh = 0.850D0
cldfrc_rhminl = 0.918D0

micro_mg_sub_version 			= 5

ext_frc_type           			= 'CYCLICAL'
ext_frc_cycle_yr       			= 1850
srf_emis_type          			= 'CYCLICAL'
srf_emis_cycle_yr      			= 1850
tracer_cnst_type       			= 'CYCLICAL'
tracer_cnst_cycle_yr   			= 1850

cmip6_forcing 				= .false.
cmip6_specifier 			= '/home/export/online1/cess/CESM_INPUT/aerosol-cmpi6/MACv2-SP_3Dfields_f09_g16_1850-2016.nc'
cmip6_type 	     			= 'CYCLICAL' !Can be set to 'CYCLICAL', 'SERIAL', 'INTERP_MISSING_MONTHS', or 'FIXED'.
cmip6_cycle_yr 				= 1850        ! used when cmip6_type is 'CYCLICAL'
cmip6_fixed_ymd 			= 0       ! used when cmip6_type is 'FIXED'
cmip6_fixed_tod 			= 0       ! used when cmip6_type is 'FIXED'

prescribed_volcaero_datapath		= '/home/export/online1/cess/CESM_INPUT/aerosol-cmpi6'
prescribed_volcaero_file	    	= 'CMIP_1850_2014_volc_v2.nc'
prescribed_volcaero_type		    = 'CYCLICAL'
prescribed_volcaero_cycle_yr		= 1850

prescribed_volc_rad_datapath		= '/home/export/online1/cess/CESM_INPUT/aerosol-cmpi6'
prescribed_volc_rad_file    		= 'CMIP_DESS_radiation_volc_v3_PI_cyclical.nc' ! use 1850-2014 monthly climatology
prescribed_volc_rad_type    		= 'CYCLICAL' 
prescribed_volc_rad_cycle_yr		= 1850

solar_data_file                		= '/home/export/online1/cess/CESM_INPUT/cmip6-forcing/SOLAR_cmip6_PI_climo.nc' ! 1850-1873 mean
solar_data_type        			    = 'FIXED'
solar_data_ymd         			    = 18500101
solar_htng_spctrl_scl			    = .true.

prescribed_ozone_datapath      		= '/home/export/online1/cess/CESM_INPUT/cmip6-forcing'
prescribed_ozone_file          		= 'new_ozone_cmip6_hist.nc'
prescribed_ozone_cycle_yr      		= 1850
prescribed_ozone_name          		= 'O3'
prescribed_ozone_type          		= 'CYCLICAL'

aerodep_flx_type           = 'CYCLICAL'
aerodep_flx_cycle_yr       = 1850
aerodep_flx_file           = '/home/export/online1/cess/CESM_INPUT/atm/cam/chem/trop_mam/aero/mam3_1.9x2.5_L30_1850clim_c130319.nc'
prescribed_aero_type       = 'CYCLICAL'
prescribed_aero_cycle_yr   = 1850
prescribed_aero_file       = '/home/export/online1/cess/CESM_INPUT/atm/cam/chem/trop_mam/aero/mam3_1.9x2.5_L30_1850clim_c130319.nc'

ch4vmr					= 808.249e-9
co2vmr					= 284.317e-6
f11vmr					= 32.1102e-12
f12vmr					= 0.0
n2ovmr                  = 273.0211e-9

EOF


cat >> $case1/user_nl_pop2 << EOF
 ltidal_mixing = .true.
 tidal_energy_file = 'tidal_energy_Canuto_MC_gx1v6_20170918.nc'
 tidal_energy_file_fmt = 'nc'
 overflows_on = .true.
 overflows_interactive = .true.

 ah&hmix_gm_nml         = 3.0e7
 ah_bkg_srfbl&hmix_gm_nml   = 3.0e7
 ah_bolus&hmix_gm_nml       = 3.0e7
 kappa_max_eg&hmix_gm_nml   = 2.0e7
 kappa_min_eg&hmix_gm_nml   = 0.35e7
 rich_mix&vmix_kpp_nml = 50.0
EOF

echo ">>>>MODIFIED ZM SCHEME                ON<<<<"
echo ">>>>SINGLE ICE MICRO-PHYSICS SCHEME   ON<<<<"
echo ">>>>PDF MACRO-PHYSICS SCHEME          ON<<<<"
echo ">>>>FOUR STREAM SCHEME                ON<<<<"
cat >> $case1/README.case << EOF
>>>>MODIFIED ZM SCHEME                ON<<<<
>>>>SINGLE ICE MICRO-PHYSICS SCHEME   ON<<<<
>>>>PDF MACRO-PHYSICS SCHEME          ON<<<<
>>>>FOUR STREAM SCHEME                ON<<<<
EOF

else

echo ">>>>MODIFIED ZM SCHEME               OFF<<<<"
echo ">>>>SINGLE ICE MICRO-PHYSICS SCHEME  OFF<<<<"
echo ">>>>PDF MACRO-PHYSICS SCHEME         OFF<<<<"
echo ">>>>FOUR STREAM SCHEME               OFF<<<<"
cat >> $case1/README.case << EOF
>>>>MODIFIED ZM SCHEME               OFF<<<<
>>>>SINGLE ICE MICRO-PHYSICS SCHEME  OFF<<<<
>>>>PDF MACRO-PHYSICS SCHEME         OFF<<<<
>>>>FOUR STREAM SCHEME               OFF<<<<
EOF

endif

###############################
echo ">>>>====================================<<<<"
cat >> $case1/README.case << EOF
>>>>====================================<<<<
EOF
