!================================================================================================
! This is the 'none' chemistry module.
! Most of the routines return without doing anything.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8, shr_kind_cl
  use physics_types,       only: physics_state, physics_ptend
  use ppgrid,              only: begchunk, endchunk, pcols
  use cam_control_mod, only: cmip6_forcing !wmq
  use cam_logfile,  only : iulog
  use spmd_utils,       only : masterproc

  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public :: chem_is                        ! identify which chemistry is being used
  public :: chem_register                  ! register consituents
  public :: chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public :: chem_implements_cnst           ! returns true if consituent is implemented by this package
  public :: chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public :: chem_init                      ! initialize (history) variables
  public :: chem_timestep_init             ! time interpolate chemical loss frequencies
  public :: chem_timestep_tend             ! interface to tendency computation
  public :: chem_final
  public :: chem_write_restart
  public :: chem_read_restart
  public :: chem_init_restart
  public :: chem_readnl                    ! read chem namelist 
  public :: chem_reset_fluxes

  interface chem_write_restart
     module procedure chem_write_restart_bin
     module procedure chem_write_restart_pio
  end interface
  interface chem_read_restart
     module procedure chem_read_restart_bin
     module procedure chem_read_restart_pio
  end interface

!--wmq:cmip6
  character(len=shr_kind_cl) :: cmip6_specifier = ''
  character(len=24)  :: cmip6_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer            :: cmip6_cycle_yr  = 0
  integer            :: cmip6_fixed_ymd = 0
  integer            :: cmip6_fixed_tod = 0
  character(len=shr_kind_cl) :: cmip6_1850 = ''

  character(len=24)  :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer            :: srf_emis_cycle_yr  = 0
  character(len=24)  :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer            :: ext_frc_cycle_yr  = 0


  ! Private data

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'none' ) then
       chem_is = .true.
    end if

  end function chem_is

!================================================================================================

  subroutine chem_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for parameterized greenhouse gas chemistry
    ! 
    !-----------------------------------------------------------------------

  end subroutine chem_register

!================================================================================================

  subroutine chem_readnl(nlfile)

    ! args
    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    integer :: unitn, ierr

    real(r8)           :: sol_facti_cloud_borne
    integer            :: tracer_cnst_cycle_yr
    character(len=shr_kind_cl) :: tracer_cnst_datapath          ! absolute path of prescribed data files
    character(len=shr_kind_cl) :: tracer_cnst_file              ! prescribed data file
    character(len=24)  :: tracer_cnst_type              ! 'INTERP_MISSING_MONTHS' | 'CYCLICAL' | 'SERIAL' (default)

    namelist /chem_inparm/ cmip6_specifier, cmip6_type, cmip6_cycle_yr, &
    cmip6_fixed_ymd, cmip6_fixed_tod, cmip6_forcing, cmip6_1850 !wmq:cmip6

    namelist /chem_inparm/ ext_frc_cycle_yr, ext_frc_type, sol_facti_cloud_borne, &
    srf_emis_cycle_yr, srf_emis_type, tracer_cnst_cycle_yr, tracer_cnst_datapath, &
    tracer_cnst_file, tracer_cnst_type

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'chem_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, chem_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('chem_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

    ! Broadcast namelist variables

    ! control

!wmq:cmip6
#ifdef SPMD
    call mpibcast (cmip6_specifier, len(cmip6_specifier),            mpichar, 0, mpicom)
    call mpibcast (cmip6_type,      len(cmip6_type),                 mpichar, 0, mpicom)
    call mpibcast (cmip6_cycle_yr,  1,                               mpiint,  0, mpicom)
    call mpibcast (cmip6_fixed_ymd, 1,                               mpiint,  0, mpicom)
    call mpibcast (cmip6_fixed_tod, 1,                               mpiint,  0, mpicom)
    call mpibcast (cmip6_forcing,   1,                               mpilog,  0, mpicom)
    call mpibcast (cmip6_1850,      len(cmip6_1850),                 mpichar, 0, mpicom)
#endif
  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .false.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    chem_implements_cnst = .false.

  end function chem_implements_cnst

!===============================================================================

  subroutine chem_init(phys_state, pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize parameterized greenhouse gas chemistry
    !          (declare history variables)
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only : physics_buffer_desc
    use cam_history,    only: addfld, add_default, phys_decomp
    use mo_chemini,          only : chemini !wmq

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
!-----------------------------------------------------------------------
! Initialize chemistry modules
!-----------------------------------------------------------------------

    call chemini &
       ( cmip6_specifier & !wmq
       , cmip6_type &
       , cmip6_cycle_yr &
       , cmip6_fixed_ymd &
       , cmip6_fixed_tod &
       , cmip6_1850 &      !--
       )


  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state, pbuf2d)
    use physics_buffer, only : physics_buffer_desc
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
         is_perpetual
    use MAC_SP_forcing,    only : MAC_SP_timestep_init !wmq

    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    ! wmq: Set CMIP6 forcings timing factors
    !-----------------------------------------------------------------------
    if(cmip6_forcing) then
      call MAC_SP_timestep_init( pbuf2d, phys_state )
    end if

  end subroutine chem_timestep_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf, fh2o, fsds )
    use physics_buffer,           only: physics_buffer_desc
    use cam_history,      only: outfld
    use camsrfexch,       only: cam_in_t, cam_out_t
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)    :: dt              ! time step
    type(physics_state), intent(in)    :: state           ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend           ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    real(r8),            intent(out)   :: fh2o(pcols)     ! h2o flux to balance source from chemistry
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    real(r8),            intent(in)    :: fsds(pcols)     ! longwave down at sfc

    return
  end subroutine chem_timestep_tend

!===============================================================================

  subroutine chem_init_cnst(name, q, gcid)

    character(len=*), intent(in) :: name         ! constituent name
    real(r8), intent(out) :: q(:,:)   !  mass mixing ratio (gcol, plev)
    integer, intent(in) :: gcid(:)    !  global column id

    return
  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_write_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_write_restart_bin
!===============================================================================
  subroutine chem_read_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_read_restart_bin
!===============================================================================
  subroutine chem_write_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_write_restart_pio
!===============================================================================
  subroutine chem_read_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_read_restart_pio
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_init_restart
!================================================================================
  subroutine chem_reset_fluxes( fptr, cam_in )
    use camsrfexch,          only : cam_in_t     

    real(r8), pointer             :: fptr(:,:)        ! pointer into    array data
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)

  end subroutine chem_reset_fluxes

end module chemistry
