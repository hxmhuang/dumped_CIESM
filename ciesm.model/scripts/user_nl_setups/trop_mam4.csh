#####Aerosol-model-MAM4#####
if ($trop_mam4 == on) then
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/build-namelist $rootmod/atm/cam/bld/build-namelist
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/definition.xml $rootmod/atm/cam/bld/config_files/definition.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/configure $rootmod/atm/cam/bld/configure
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/master_mam_dep_list.xml $rootmod/atm/cam/bld/namelist_files/master_mam_dep_list.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/namelist_defaults_cam.xml $rootmod/atm/cam/bld/namelist_files/namelist_defaults_cam.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/namelist_definition.xml $rootmod/atm/cam/bld/namelist_files/namelist_definition.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/seq_drydep_mod.F90 $rootmod/drv/shr/seq_drydep_mod.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/atm_comp_esmf.F90 $rootmod/atm/cam/src/cpl_esmf/atm_comp_esmf.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/atm_comp_mct.F90 $rootmod/atm/cam/src/cpl_mct/atm_comp_mct.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/1850-2005_cam5.xml $rootmod/atm/cam/bld/namelist_files/use_cases/1850-2005_cam5.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/1850_cam5.xml $rootmod/atm/cam/bld/namelist_files/use_cases/1850_cam5.xml
rm -rf $rootmod/atm/cam/src/chemistry
cp -rf $rootpwd/SourceMods.CAM5/trop_mam4/CESM_MAM4_CODE/chemistry $rootmod/atm/cam/src/
else
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/build-namelist $rootmod/atm/cam/bld/build-namelist
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/config_files/definition.xml $rootmod/atm/cam/bld/config_files/definition.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/configure $rootmod/atm/cam/bld/configure
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/namelist_files/master_mam_dep_list.xml $rootmod/atm/cam/bld/namelist_files/master_mam_dep_list.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/namelist_files/namelist_defaults_cam.xml $rootmod/atm/cam/bld/namelist_files/namelist_defaults_cam.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/namelist_files/namelist_definition.xml $rootmod/atm/cam/bld/namelist_files/namelist_definition.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/drv/shr/seq_drydep_mod.F90 $rootmod/drv/shr/seq_drydep_mod.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/src/cpl_esmf/atm_comp_esmf.F90 $rootmod/atm/cam/src/cpl_esmf/atm_comp_esmf.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/src/cpl_mct/atm_comp_mct.F90 $rootmod/atm/cam/src/cpl_mct/atm_comp_mct.F90
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/namelist_files/use_cases/1850-2005_cam5.xml $rootmod/atm/cam/bld/namelist_files/use_cases/1850-2005_cam5.xml
cp $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/bld/namelist_files/use_cases/1850_cam5.xml $rootmod/atm/cam/bld/namelist_files/use_cases/1850_cam5.xml
rm -rf $rootmod/atm/cam/src/chemistry
cp -rf $rootpwd/SourceMods.CAM5/trop_mam4/models/atm/cam/src/chemistry $rootmod/atm/cam/src/
endif
