The present CIESM package includes:

(1) ciesm.model: includes earth system model codes developed by DESS of Tsinghua University based on NCAR CESM_v1.2.1. 

(2) install.csh: This script is used to setup model run and archive directories （ciesm.run, ciesm.archive） and setup config_machines.xml. 

(3) gen_newcase.csh: This script is used to setup and create model experiments.

================================================= How to use:

(1) Execute ./install.csh to create ciesm.run, ciesm.archive directories.

(2) Edit config_machines.xml, mkbatch.$machine and config_compilers.xml (in the path .../ciesm.model/scripts/ccsm_utils/Machines/) to setup your machine, job system, initial data directory and compiling environment.

(3) Edit gen_newcase.csh to setup your experiments and execute gen_newcase.csh to generate your experiments.

(4) Goto your experiment directory in ciesm.run, then execute the $CASE.submit script to submit.

(5) All output data will be moved to ciesm.archive/$CASE directory after your job finished.
