#!/bin/bash
cd /home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/atm/obj/
chmod +x a.sh
./a.sh
cd /home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/cesm/obj/
mpif90.new  -o /home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/cesm.exe ccsm_comp_mod.o ccsm_driver.o mrg_mod.o seq_avdata_mod.o seq_diag_mct.o seq_domain_mct.o seq_flux_mct.o seq_frac_mct.o seq_hist_mod.o seq_map_esmf.o seq_map_mod.o seq_mctext_mod.o seq_rest_mod.o  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -latm  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -lice  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -llnd  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -locn  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -lrof  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -lglc  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib/ -lwav -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/lib -lcsm_share  -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/mct/mct -lmct -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/mct/mpeu -lmpeu -L/home/export/online1/qhtest/swqh/CESM/cesm_v1.2.1/scripts/F2000C5_f19g16_qhtest/bld/pio -lpio -lgptl -L/home/export/online1/qhtest/swqh/netcdf/lib -lnetcdf -lm_old  -L/home/export/online1/qhtest/swqh/netcdf/lib -lnetcdf -lm_old