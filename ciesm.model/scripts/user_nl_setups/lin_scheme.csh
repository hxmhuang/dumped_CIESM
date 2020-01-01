#! /bin/csh -f

if ($lin_scheme == on) then

cp -r $rootpwd/SourceMods.CAM5/CESM-LIN/* $case1/SourceMods/src.cam

cat >> $case1/user_nl_cam << EOF 
!Namelist setups for LIN scheme
fincl1 = "PRECZ","PRECSH","ZMDT","ZMDQ","ZMDLIQ","CMFDT","CMFDQ","CMFDLIQ","ZMEU","ZMED","ZMMU","ZMMD","ZMDU","CAPE","dCAPE","MPDT","MPDQ","MPDLIQ","MACPDT","MACPDQ","MACPDLIQ","MACPDICE","DPDLFT","SHDLFT","DPDLFLIQ","DPDLFICE","AST","CONCLD","CBMF","tten_PBL","qlten_PBL","sgm_tota","sgm_shal","sgm_turb","qtu_shal","umf_shal","clddep2","wstarPBL","dqwdz","dthldz","KVH","beta","TKE","aa","bb","UW_sh","UW_leng","DTCOND","DTCORE","UW_turbtype","UW_qrl","QQ_out1","QQ_out2","T_before","qv_before","ql_before","T_after","qv_after","ql_after","CMELIQ","ql_st","delta_q"

micro_mg_sub_version   =  5

zmconv_c0_ocn      =  0.00650D0

EOF

echo ">>>>LIN SCHEME                        ON<<<<"

else 
   
echo ">>>>LIN SCHEME                       OFF<<<<"

endif

