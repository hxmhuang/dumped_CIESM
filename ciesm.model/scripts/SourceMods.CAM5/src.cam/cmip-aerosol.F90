module MAC_SP_forcing
!---------------------------------------------------------------------------------
! Purpose:
!   Get inputdata to run CMIP6 historical run (1850-2016)
!   Author: Wang minqi
!   begun: September 2017
!   
!   Data get from MACv2-SP (Bjorn Stevens et al.,2017) 

  !---------------------------------------------------------------
  ! 	... insitu forcing module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ppgrid,       only : pcols, begchunk, endchunk, pver, pverp
  use chem_mods,    only : gas_pcnst, extcnt
  use spmd_utils,   only : masterproc,iam
  use abortutils,   only : endrun
  use cam_history,  only : addfld, outfld, phys_decomp, add_default
  use cam_logfile,  only : iulog
  use tracer_data,  only : trfld,trfile
  use cam_control_mod, only: N_CMIP

  implicit none

  type :: forcing_cmip
     integer           :: frc_ndx
     real(r8)              :: mw
     character(len=265) :: filename
     real(r8), pointer     :: times(:)
     real(r8), pointer     :: levi(:)
     character(len=8)  :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type forcing_cmip

  private
  public  :: MAC_SP_inti
  public  :: MAC_SP_set_cdnc
  public  :: MAC_SP_set_rad
  public  :: MAC_SP_timestep_init

  save

  integer, parameter :: time_span = 1

  character(len=256) ::   filename

  type(forcing_cmip), allocatable  :: forcings_cmip(:)
  integer :: extfrc_cnt = 0

!  integer, parameter :: N_CMIP = 4
  character(len=11)    :: input_names(N_CMIP) = (/'aod','ssa','asy','dNovrN'/)

  integer :: number_flds

  logical :: has_extfrc(N_CMIP)

contains

  subroutine MAC_SP_inti( cmip_specifier, cmip_type, cmip_cycle_yr, cmip_fixed_ymd, cmip_fixed_tod)

    !-----------------------------------------------------------------------
    ! 	... initialize the CMIP6 forcings
    !-----------------------------------------------------------------------
    use cam_pio_utils, only : cam_pio_openfile
    use pio, only : pio_inq_dimid, pio_inquire, pio_inq_varndims, pio_closefile, &
         pio_inq_varname, pio_nowrite, file_desc_t
    use mo_tracname,   only : solsym
    use tracer_data,   only : trcdata_init
    use phys_control,  only : phys_getopts
    use physics_buffer, only : physics_buffer_desc

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: cmip_specifier
    character(len=*), intent(in) :: cmip_type
    integer  , intent(in)        :: cmip_cycle_yr
    integer  , intent(in)        :: cmip_fixed_ymd
    integer  , intent(in)        :: cmip_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer :: astat
    integer :: j, l, m, n, i,mm                          ! Indices
    character(len=16)  :: species
    character(len=16)  :: spc_name
    character(len=256) :: locfn
    character(len=256) :: spc_fnames(N_CMIP)

    integer ::  vid, ndims, nvars, isec, ierr
    type(file_desc_t) :: ncid
    character(len=32)  :: varname

    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies

    !-----------------------------------------------------------------------
    call phys_getopts( history_aerosol_out        = history_aerosol   )

    do i = 1, N_CMIP
       has_extfrc(i) = .false.
       spc_fnames(i) = ''
    enddo

    !-----------------------------------------------------------------------
    ! 	... species has insitu forcing ?
    !-----------------------------------------------------------------------

    !write(iulog,*) 'Species with insitu forcings'

       if ( len_trim(cmip_specifier) == 0 ) then
            write(iulog,*) 'MAC_SP_inti: can not find the path of cmip6_specifier'
            call endrun
       endif

    count_emis: do n=1,N_CMIP

       spc_name = trim(input_names(n))
       filename = trim(cmip_specifier)

       m = get_cmip_ndx( spc_name )

       if ( m < 1 ) then
          call endrun('MAC_SP_inti: '//trim(spc_name)// ' does not have an external source')
       endif

       spc_fnames(n) = filename

       has_extfrc(n) = .true.
       write(iulog,*) '   ',  spc_name ,' : filename = ',trim(spc_fnames(n)),' spc ndx = ',n

    enddo count_emis

    extfrc_cnt = count( has_extfrc(:) )

    if( extfrc_cnt < 1 ) then
       if (masterproc) write(iulog,*) 'There are no species with insitu forcings'
       return
    end if

    if (masterproc) write(iulog,*) ' '

    !-----------------------------------------------------------------------
    ! 	... allocate forcings type array
    !-----------------------------------------------------------------------
    allocate( forcings_cmip(extfrc_cnt), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'MAC_SP_inti: failed to allocate forcings array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the forcing type array
    !-----------------------------------------------------------------------
    n = 0
    species_loop : do m = 1,N_CMIP
       has_forcing : if( has_extfrc(m) ) then

          allocate( forcings_cmip(m)%sectors(1), stat=astat )
          if( astat/= 0 ) then
             write(iulog,*) 'MAC_SP_init: failed to allocate forcings_cmip%sectors array; error = ',astat
             call endrun
          end if

          allocate( forcings_cmip(m)%fields(1), stat=astat )
          if( astat/= 0 ) then
             write(iulog,*) 'MAC_SP_init: failed to allocate forcings_cmip%fields array; error = ',astat
             call endrun
          end if

          spc_name = trim( input_names(m) )
          n        = n + 1
          !-----------------------------------------------------------------------
          ! 	... default settings
          !-----------------------------------------------------------------------
          forcings_cmip(n)%frc_ndx          = get_cmip_ndx( spc_name )
          forcings_cmip(n)%species          = spc_name
          forcings_cmip(n)%sectors          = spc_name
          forcings_cmip(n)%nsectors         = 1
          forcings_cmip(n)%filename         = spc_fnames(m)
         
          if(trim(forcings_cmip(m)%species)==trim('dNovrN')) then
             call addfld( trim(spc_name),  '', 1, 'A', &
                         'cmip6 forcing for '//trim(spc_name),   phys_decomp )
          else
             call addfld( trim(spc_name),  '', pver, 'A', &
                         'cmip6 forcing for '//trim(spc_name),   phys_decomp )
          end if
          if ( history_aerosol ) then 
             call add_default( trim(spc_name), 1, ' ' )
          endif
       end if has_forcing
    end do species_loop

    if (masterproc) then
       !-----------------------------------------------------------------------
       ! 	... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'MAC_SP_inti: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'MAC_SP timing specs'
       write(iulog,*) 'type = ',cmip_type
       if( cmip_type == 'FIXED' ) then
          write(iulog,*) ' fixed date = ', cmip_fixed_ymd
          write(iulog,*) ' fixed time = ', cmip_fixed_tod
       else if( cmip_type == 'CYCLICAL' ) then
          write(iulog,*) ' cycle year = ', cmip_cycle_yr
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ',extfrc_cnt,' species with external forcing files'
       do m = 1,extfrc_cnt
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',m
          write(iulog,*) 'species = ',trim(forcings_cmip(m)%species)
          write(iulog,*) 'frc ndx = ',forcings_cmip(m)%frc_ndx
          write(iulog,*) 'filename= ',trim(forcings_cmip(m)%filename)
       end do
       write(iulog,*) ' '
    endif

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    frcing_loop: do m = 1, extfrc_cnt

       allocate(forcings_cmip(m)%file%in_pbuf(size(forcings_cmip(m)%sectors)))
       forcings_cmip(m)%file%in_pbuf(:) = .false.
       call trcdata_init( forcings_cmip(m)%sectors, &
                          forcings_cmip(m)%filename, filelist, datapath, &
                          forcings_cmip(m)%fields,  &
                          forcings_cmip(m)%file, &
                          rmv_file, cmip_cycle_yr, cmip_fixed_ymd, cmip_fixed_tod, cmip_type)

       number_flds = 0
       if (associated(forcings_cmip(m)%fields)) number_flds = size( forcings_cmip(m)%fields )

       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(iulog,*) 'There are no cimp6 data'
             write(iulog,*) ' '
             call endrun
          endif
       end if

    enddo frcing_loop


  end subroutine MAC_SP_inti

  subroutine MAC_SP_timestep_init( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,extfrc_cnt
       call advance_trcdata( forcings_cmip(m)%fields, forcings_cmip(m)%file, state, pbuf2d  )
    end do

  end subroutine MAC_SP_timestep_init

  subroutine MAC_SP_set_cdnc( lchnk, flux, ncol )
    !--------------------------------------------------------
    !   ... form the external forcing
    !--------------------------------------------------------
    use phys_grid,    only : pcols, begchunk, endchunk

    implicit none

    !--------------------------------------------------------
    !   ... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    :: ncol                  ! columns in chunk
    integer,  intent(in)    :: lchnk                 ! chunk index
    real(r8), intent(inout) :: flux(pcols,N_CMIP)
    !--------------------------------------------------------
    !   ... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n, c
    character(len=16) :: xfcname
    integer  :: k, isec

    if( extfrc_cnt < 1) then
       return
    end if
 
    do m = 1,extfrc_cnt

       n = forcings_cmip(m)%frc_ndx

       flux(:ncol,n) = 0._r8

       do isec = 1,forcings_cmip(m)%nsectors
              if(trim(forcings_cmip(m)%species)==trim('dNovrN')) then
                 flux(:ncol,n) = flux(:ncol,n) + forcings_cmip(m)%fields(isec)%data(:ncol,1,lchnk)
              end if
       enddo

       xfcname = trim(forcings_cmip(m)%species)
       if(trim(forcings_cmip(m)%species)==trim('dNovrN')) then
          call outfld( xfcname, flux(:ncol,n), ncol, lchnk )
       end if

    end do 

  end subroutine MAC_SP_set_cdnc

  subroutine MAC_SP_set_rad( lchnk, frcing, ncol )
    !--------------------------------------------------------
    !	... form the external forcing
    !--------------------------------------------------------
    use phys_grid,    only : pcols, begchunk, endchunk

    implicit none

    !--------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)    :: ncol                  ! columns in chunk
    integer,  intent(in)    :: lchnk                 ! chunk index
    real(r8), intent(inout) :: frcing(pcols,pver,N_CMIP)   ! insitu forcings (molec/cm^3/s)

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n, c
    character(len=16) :: xfcname
    integer  :: k, isec

    if( extfrc_cnt < 1 ) then
       return
    end if

    !--------------------------------------------------------
    !	... set non-zero forcings
    !--------------------------------------------------------
!       frcing(:,:,:) = 0._r8
!       flux(:,:) = 0._r8

    do m = 1,extfrc_cnt

       n = forcings_cmip(m)%frc_ndx

       frcing(:ncol,:,n) = 0._r8

       do isec = 1,forcings_cmip(m)%nsectors
              if(trim(forcings_cmip(m)%species) .ne. trim('dNovrN')) then
                 frcing(:ncol,:,n) = frcing(:ncol,:,n) + forcings_cmip(m)%fields(isec)%data(:ncol,pver:1:-1,lchnk)
              end if
       enddo

       xfcname = trim(forcings_cmip(m)%species)
       if(trim(forcings_cmip(m)%species) .ne. trim('dNovrN')) then
          call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )
       end if
 
    end do 

  end subroutine MAC_SP_set_rad

 integer function get_cmip_ndx( name )

    implicit none
    character(len=*), intent(in) :: name

    integer :: i

    get_cmip_ndx = 0
    do i = 1,N_CMIP
      if ( trim(name) == trim(input_names(i)) ) then
        get_cmip_ndx = i
        return
      endif
    enddo

  end function get_cmip_ndx

end module MAC_SP_forcing
