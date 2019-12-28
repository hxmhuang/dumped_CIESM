The present CIESM package includes:

(1) SourceMods.tar: SourceMods files of CIESM version 1.1
(2) namelist_definition.tar: Namelist definition files of CIESM version 1.1

=================================================
How to use:
=================================================

 CIESM v1.1 is developed based on NCAR CESM version 1.2.1.

 (1) Copy the namelist definition files in namelist_definition.tar into the corresponding model directories:
     ---cam.namelist_definition.xml     -> .../models/atm/cam/bld/namelist_files/namelist_definition.xml
     ---pop.namelist_definition.xml     -> .../models/ocn/pop2/bld/namelist_files/namelist_definition.xml
     ---cpl.namelist_definition_drv.xml -> .../models/drv/bld/namelist_files/namelist_definition_drv.xml

 (2) Use the same methods to create and configure a case as in CESM v1.2.1.

 (3) Copy the SourceMods folder in SourceMods.tar into your $CASE directory before compiling the case.

 (4) Use the same methods to compile and run a case as in CESM v1.2.1.
