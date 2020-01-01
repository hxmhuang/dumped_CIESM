
module mo_chemini

  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmd_utils,   only : masterproc
  use cam_logfile,  only : iulog
  use cam_control_mod, only: cmip6_forcing !wmq

  implicit none

  private
  public :: chemini

contains

!---wmq----
  subroutine chemini &
       ( cmip6_specifier &
       , cmip6_type &
       , cmip6_cycle_yr &
       , cmip6_fixed_ymd &
       , cmip6_fixed_tod &   
       , cmip6_1850 & 
       )
!-------------   

    !-----------------------------------------------------------------------
    ! 	... Chemistry module intialization
    !-----------------------------------------------------------------------

    use spmd_utils,        only : iam
    use MAC_SP_forcing,    only : MAC_SP_inti  !wmq
    implicit none

!---wmq:CMIP6--
    character(len=*), intent(in) :: cmip6_specifier
    character(len=*), intent(in) :: cmip6_type
    integer,          intent(in) :: cmip6_cycle_yr
    integer,          intent(in) :: cmip6_fixed_ymd
    integer,          intent(in) :: cmip6_fixed_tod
    character(len=*), intent(in) :: cmip6_1850
!--------------
    !wmq
    !-----------------------------------------------------------------------
    !   ... initialize CMIP6 forcings module
    !-----------------------------------------------------------------------
    if(cmip6_forcing) then
       call MAC_SP_inti(cmip6_specifier, cmip6_type, cmip6_cycle_yr, cmip6_fixed_ymd, cmip6_fixed_tod)
       if (masterproc) write(iulog,*) 'chemini: after MAC_SP_inti on node ',iam
    end if

  end subroutine chemini

end module mo_chemini
