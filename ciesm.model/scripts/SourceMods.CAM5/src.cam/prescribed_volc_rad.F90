!-------------------------------------------------------------------
! manages reading and interpolation of volcanic aerosol radiative properties
! Created by: Wang Minqi
!-------------------------------------------------------------------
module prescribed_volc_rad

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog
  use ppgrid,       only : pcols, begchunk, endchunk, pver, pverp

  implicit none

  type :: forcing_volc
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
  end type forcing_volc

  private
  public :: prescribed_volc_rad_readnl
  public :: prescribed_volc_rad_init
  public :: prescribed_volc_rad_adv
  public :: has_prescribed_volc_rad
  public :: prescribed_volc_rad_set

  save

  logical :: has_prescribed_volc_rad = .false.

  ! These variables are settable via the namelist (with longer names)
  integer, parameter :: N_rad = 6
  character(len=16)  :: input_names(N_rad) = (/'ext_earth','omega_earth','g_earth','ext_sun','omega_sun','g_sun'/)
  character(len=256) :: spc_fnames(N_rad)
  character(len=16)  :: spc_name
  integer :: number_flds

  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  type(forcing_volc), allocatable  :: forcings_rad(:)

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_volc_rad_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_volc_rad_readnl'

   character(len=256) :: prescribed_volc_rad_file
   character(len=256) :: prescribed_volc_rad_filelist
   character(len=256) :: prescribed_volc_rad_datapath
   character(len=32)  :: prescribed_volc_rad_type
   logical            :: prescribed_volc_rad_rmfile
   integer            :: prescribed_volc_rad_cycle_yr
   integer            :: prescribed_volc_rad_fixed_ymd
   integer            :: prescribed_volc_rad_fixed_tod

   namelist /prescribed_volc_rad_nl/ &
      prescribed_volc_rad_file,      &
      prescribed_volc_rad_filelist,  &
      prescribed_volc_rad_datapath,  &
      prescribed_volc_rad_type,      &
      prescribed_volc_rad_rmfile,    &
      prescribed_volc_rad_cycle_yr,  &
      prescribed_volc_rad_fixed_ymd, &
      prescribed_volc_rad_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_volc_rad_file     = filename
   prescribed_volc_rad_filelist = filelist
   prescribed_volc_rad_datapath = datapath
   prescribed_volc_rad_type     = data_type
   prescribed_volc_rad_rmfile   = rmv_file
   prescribed_volc_rad_cycle_yr = cycle_yr
   prescribed_volc_rad_fixed_ymd= fixed_ymd
   prescribed_volc_rad_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_volc_rad_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_volc_rad_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_volc_rad_file,     len(prescribed_volc_rad_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_rad_filelist, len(prescribed_volc_rad_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_rad_datapath, len(prescribed_volc_rad_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_rad_type,     len(prescribed_volc_rad_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volc_rad_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_volc_rad_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volc_rad_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volc_rad_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   filename   = prescribed_volc_rad_file
   filelist   = prescribed_volc_rad_filelist
   datapath   = prescribed_volc_rad_datapath
   data_type  = prescribed_volc_rad_type
   rmv_file   = prescribed_volc_rad_rmfile
   cycle_yr   = prescribed_volc_rad_cycle_yr
   fixed_ymd  = prescribed_volc_rad_fixed_ymd
   fixed_tod  = prescribed_volc_rad_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_volc_rad = .true.

end subroutine prescribed_volc_rad_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volc_rad_init()

    use tracer_data,    only : trcdata_init
    use cam_history,    only : addfld, phys_decomp
    use ppgrid,         only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc, pbuf_get_index
    use phys_control,   only : phys_getopts

    implicit none

    integer :: ndx, astat, n
    integer :: errcode
    logical  :: history_aerosol      ! Output the MAM aerosol tendencies

    call phys_getopts( history_aerosol_out        = history_aerosol   )
    
    if ( has_prescribed_volc_rad ) then
       if ( masterproc ) then
          write(iulog,*) 'volcanic aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate( forcings_rad(N_rad), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'prescribed_volc_rad_init: failed to allocate forcings array; error = ',astat
       call endrun
    end if

do n=1,N_rad
    spc_fnames(n) = filename
    spc_name = trim( input_names(n) )
    write(iulog,*) '   ',  spc_name ,' : filename = ',trim(spc_fnames(n)),' spc ndx = ',n

    allocate( forcings_rad(n)%sectors(1), stat=astat )
    if( astat/= 0 ) then
        write(iulog,*) 'prescribed_volc_rad_init: failed to allocate forcings_rad%sectors array; error = ',astat
        call endrun
    end if

    allocate( forcings_rad(n)%fields(1), stat=astat )
    if( astat/= 0 ) then
        write(iulog,*) 'MAC_SP_init: failed to allocate forcings_rad%fields array; error = ',astat
        call endrun
    end if

    !-----------------------------------------------------------------------
    !     ... default settings
    !-----------------------------------------------------------------------
          forcings_rad(n)%frc_ndx          = n
          forcings_rad(n)%species          = spc_name
          forcings_rad(n)%sectors          = spc_name
          forcings_rad(n)%nsectors         = 1
          forcings_rad(n)%filename         = spc_fnames(n)

!          if(trim(forcings_cmip(n)%species)==trim('dNovrN')) then
!             call addfld( trim(spc_name),  '', 1, 'A', &
!                         'cmip6 stratospheric aerosols for '//trim(spc_name),   phys_decomp )
!          else
!             call addfld( trim(spc_name),  '', pver, 'A', &
!                         'cmip6 stratospheric aerosols for '//trim(spc_name),   phys_decomp )
!          end if
!          if ( history_aerosol ) then
!             call add_default( trim(spc_name), 1, ' ')
!          endif
end do

    if (masterproc) then
       !-----------------------------------------------------------------------
       !        ... diagnostics
       !-----------------------------------------------------------------------
       write(iulog,*) ' '
       write(iulog,*) 'prescribed_volc_rad_init: diagnostics'
       write(iulog,*) ' '
       write(iulog,*) 'volc_rad timing specs'
       write(iulog,*) 'type = ',data_type
       if( data_type == 'FIXED' ) then
          write(iulog,*) ' fixed date = ', fixed_ymd
          write(iulog,*) ' fixed time = ', fixed_tod
       else if( data_type == 'CYCLICAL' ) then
          write(iulog,*) ' cycle year = ', cycle_yr
       end if
       write(iulog,*) ' '
       write(iulog,*) 'there are ', N_rad,' species with external forcing files'
       do n = 1,N_rad
          write(iulog,*) ' '
          write(iulog,*) 'forcing type ',n
          write(iulog,*) 'species = ',trim(forcings_rad(n)%species)
          write(iulog,*) 'frc ndx = ',forcings_rad(n)%frc_ndx
          write(iulog,*) 'filename= ',trim(forcings_rad(n)%filename)
       end do
       write(iulog,*) ' '
    endif

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    do n = 1, N_rad

       allocate(forcings_rad(n)%file%in_pbuf(size(forcings_rad(n)%sectors)))
       forcings_rad(n)%file%in_pbuf(:) = .false.
       call trcdata_init( forcings_rad(n)%sectors, &
                          forcings_rad(n)%filename, filelist, datapath, &
                          forcings_rad(n)%fields,  &
                          forcings_rad(n)%file, &
                          rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

       number_flds = 0
       if (associated(forcings_rad(n)%fields)) number_flds = size( forcings_rad(n)%fields )

       if( number_flds < 1 ) then
          if ( masterproc ) then
             write(iulog,*) 'There are no cimp6 data'
             write(iulog,*) ' '
             call endrun
          endif
       end if

    enddo

  end subroutine prescribed_volc_rad_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volc_rad_adv( state, pbuf2d)
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
    integer :: n

    if( .not. has_prescribed_volc_rad ) return

    do n = 1,N_rad
       call advance_trcdata( forcings_rad(n)%fields, forcings_rad(n)%file, state, pbuf2d  )
    end do

  end subroutine prescribed_volc_rad_adv

  subroutine prescribed_volc_rad_set( lchnk, frcing_sw, frcing_lw, ncol )
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
    real(r8), intent(inout) :: frcing_sw(pcols,pver,14,3)  
    real(r8), intent(inout) :: frcing_lw(pcols,pver,16,3)
    !--------------------------------------------------------
    !   ... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n, c
    character(len=16) :: xfcname
    integer  :: k, isec

    if( N_rad < 1 ) then
       return
    end if

    !--------------------------------------------------------
    !   ... set non-zero forcings
    !--------------------------------------------------------
    n=0
    src_loop1 : do m = 1, 3
       n = n+1
       frcing_lw(:ncol,:,:,n) = 0._r8
       do isec = 1,forcings_rad(m)%nsectors
          frcing_lw(:ncol,:,:,n) = frcing_lw(:ncol,:,:,n) + forcings_rad(m)%fields(isec)%data_volc(:ncol,pver:1:-1,lchnk,:)
       enddo
!       xfcname = trim(forcings_cmip(m)%species)
!       call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )
   end do src_loop1

    n=0
    src_loop2 : do m = 4, 6
       n = n+1
       frcing_sw(:ncol,:,:,n) = 0._r8
       do isec = 1,forcings_rad(m)%nsectors
          frcing_sw(:ncol,:,:,n) = frcing_sw(:ncol,:,:,n) + forcings_rad(m)%fields(isec)%data_volc(:ncol,pver:1:-1,lchnk,:)
       enddo

!      write(iulog,*) 'frcing_sw(5,:pver,9,1)=',frcing_sw(5,:pver,9,1)
!      call endrun
!       xfcname = trim(forcings_cmip(m)%species)
!       call outfld( xfcname, frcing(:ncol,:,n), ncol, lchnk )
   end do src_loop2

  end subroutine prescribed_volc_rad_set

end module prescribed_volc_rad
