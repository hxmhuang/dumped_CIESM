
program ccsm_driver


!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CCSM4/cpl7. Can have different
!          land, sea-ice, and ocean models plugged in at compile-time.
!          These models can be either: stub, dead, data, or active
!          components or some combination of the above.
!
!               stub -------- Do nothing.
!               dead -------- Send analytic data back.
!               data -------- Send data back interpolated from input files.
!               active ------ Prognostically simulate the given component.
!
! Method: Call appropriate initialization, run (time-stepping), and 
!         finalization routines.
! 
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use shr_sys_mod,       only: shr_sys_abort
   use perf_mod
   use ESMF
   use ccsm_comp_mod
   use seq_comm_mct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   use cfio
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


   implicit none


#include <mpif.h>
   !--------------------------------------------------------------------------
   ! Local Variables
   !--------------------------------------------------------------------------
   integer                    :: localrc
#ifdef USE_ESMF_LIB
   character(len=ESMF_MAXSTR) :: compName
   type(ESMF_CplComp)         :: ccsmComp
#endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   integer                    :: mype, tmppe, tmpcom, tmpionum, ierr

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!#ifdef CFIO_USED
!   write( *, *) "CFIO_CLIENT_NUM ", CFIO_CLIENT_NUM, " CFIO_IO_NUM ", CFIO_IO_NUM 
!#endif
   !--------------------------------------------------------------------------
   ! Setup and initialize the communications and logging.  
   !--------------------------------------------------------------------------
   call ccsm_pre_init()

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   call mpi_comm_rank(MPI_Comm_World, mype, ierr)
   tmpcom = MPI_COMM_NULL
   tmppe = CFIO_CLIENT_NUM
   if(mype >= tmppe .and. mype < tmppe + CFIO_CPLIO_NUM) then
        tmpcom = seq_comms(CPLID)%iocom
        tmpionum = seq_comms(CPLID)%ionum
   endif

   tmppe = tmppe + CFIO_CPLIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_ATMIO_NUM) then
        tmpcom = seq_comms(ATMID(1))%iocom
        tmpionum = seq_comms(ATMID(1))%ionum
   endif

   tmppe = tmppe + CFIO_ATMIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_LNDIO_NUM) then
        tmpcom = seq_comms(LNDID(1))%iocom
        tmpionum = seq_comms(LNDID(1))%ionum
   endif

   tmppe = tmppe + CFIO_LNDIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_OCNIO_NUM) then
        tmpcom = seq_comms(OCNID(1))%iocom
        tmpionum = seq_comms(OCNID(1))%ionum
   endif

   tmppe = tmppe + CFIO_OCNIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_ICEIO_NUM) then
        tmpcom = seq_comms(ICEID(1))%iocom
        tmpionum = seq_comms(ICEID(1))%ionum
   endif

   tmppe = tmppe + CFIO_ICEIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_GLCIO_NUM) then
        tmpcom = seq_comms(GLCID(1))%iocom
        tmpionum = seq_comms(GLCID(1))%ionum
   endif

   tmppe = tmppe + CFIO_GLCIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_ROFIO_NUM) then
        tmpcom = seq_comms(ROFID(1))%iocom
        tmpionum = seq_comms(ROFID(1))%ionum
   endif

   tmppe = tmppe + CFIO_ROFIO_NUM

   if(mype >= tmppe .and. mype < tmppe + CFIO_WAVIO_NUM) then
        tmpcom = seq_comms(WAVID(1))%iocom
        tmpionum = seq_comms(WAVID(1))%ionum
   endif

   tmppe = tmppe + CFIO_WAVIO_NUM
   !print *,"[wwj], iocom", seq_comms(OCNID(1))%iocom
   ierr = cfio_init(seq_comms(OCNID(1))%iocom, POP_NX_BLOCKS, POP_NY_BLOCKS, 8)
    if(tmpcom == MPI_COMM_NULL) then

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !--------------------------------------------------------------------------
   ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
   ! because it is needed for the time manager, even if the ESMF_INTERFACE
   ! is not used.
   !--------------------------------------------------------------------------
   call ESMF_Initialize()

#ifdef USE_ESMF_LIB

   !--------------------------------------------------------------------------
   ! Create the "Cap" component and set the services using the register
   ! routine.  Setting the services is where the initialize, run and 
   ! finalize routines for this component are set.
   !--------------------------------------------------------------------------
   compName = "CESM_Component"
   ccsmComp = ESMF_CplCompCreate(name=compName, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to create CESM Component')

   call ESMF_CplCompSetServices(ccsmComp, userRoutine=ccsm_comp_register, &
                                rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to set services for CESM Comp')

   !--------------------------------------------------------------------------
   ! Call the initialize, run and finalize routines registered with the
   ! cap component.
   !--------------------------------------------------------------------------
   call ESMF_CplCompInitialize(ccsmComp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf initialize')

   call ESMF_CplCompRun(ccsmComp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf run')

   call ESMF_CplCompFinalize(ccsmComp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf finalize')

#else

   !--------------------------------------------------------------------------
   ! If ESMF is not defined, then just call the initialize, run and finalize
   ! routines directly.
   !--------------------------------------------------------------------------
   call ccsm_init()
   call ccsm_run()
   call ccsm_final()

#endif

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   endif
   ierr = cfio_finalize()

   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )


end program ccsm_driver
