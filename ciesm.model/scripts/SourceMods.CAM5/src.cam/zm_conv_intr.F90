

module zm_conv_intr
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the Zhang-McFarlane deep convection scheme
!
! Author: D.B. Coleman
! January 2010 modified by J. Kay to add COSP simulator fields to physics buffer
!---------------------------------------------------------------------------------
   !++wy 
   !use shr_kind_mod, only: r8=>shr_kind_r8
   use shr_kind_mod, only: r8=>shr_kind_r8, i8=>shr_kind_i8
   !--wy
   use physconst,    only: cpair                              
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use zm_conv,      only: zm_conv_evap, zm_convr, convtran, momtran, zm_microp, zm_stc
!<songxl 2014-11-20----------
   use zm_microphysics,  only: zm_aero_t
   use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                               rad_cnst_get_aer_props, rad_cnst_get_mode_props  !, &
!                               rad_cnst_get_mode_num_idx, rad_cnst_get_mam_mmr_idx
   use ndrop_bam,        only: ndrop_bam_init
   use abortutils,       only: endrun
   use physconst,        only: pi
   use spmd_utils,       only: masterproc
!>songxl 2014-11-20----------
   use cam_history,  only: outfld, addfld, add_default, phys_decomp
   use perf_mod
   use cam_logfile,  only: iulog
   
   implicit none
   private
   save

   ! Public methods

   public ::&
      zm_conv_register,           &! register fields in physics buffer
      zm_conv_init,               &! initialize donner_deep module
      zm_conv_tend,               &! return tendencies
      zm_conv_tend_2               ! return tendencies

   ! Private module data

   real(r8), allocatable, dimension(:,:,:) :: mu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: eu  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: du  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: md  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: ed  !(pcols,pver,begchunk:endchunk)
   real(r8), allocatable, dimension(:,:,:) :: dp  !(pcols,pver,begchunk:endchunk) 
	! wg layer thickness in mbs (between upper/lower interface).
   real(r8), allocatable, dimension(:,:)   :: dsubcld  !(pcols,begchunk:endchunk)
	! wg layer thickness in mbs between lcl and maxi.

   integer, allocatable, dimension(:,:) :: jt   !(pcols,begchunk:endchunk)
        ! wg top  level index of deep cumulus convection.
   integer, allocatable, dimension(:,:) :: maxg !(pcols,begchunk:endchunk)
        ! wg gathered values of maxi.
   integer, allocatable, dimension(:,:) :: ideep !(pcols,begchunk:endchunk)               
	! w holds position of gathered points vs longitude index

   integer, allocatable, dimension(:) :: lengath !(begchunk:endchunk)

   integer ::& ! indices for fields in the physics buffer
      dp_flxprc_idx, &
      dp_flxsnw_idx, &
      dp_cldliq_idx, &
      dp_cldice_idx, &
!<songxl 2014-11-20------------
      dlfzm_idx,     &     ! detrained convective cloud water mixing ratio.
      difzm_idx,     &     ! detrained convective cloud ice mixing ratio.
      dnlfzm_idx,    &     ! detrained convective cloud water num concen.
      dnifzm_idx,    &     ! detrained convective cloud ice num concen.
!>songxl 2014-11-20------------
      prec_dp_idx,   &
      snow_dp_idx
!++wy
   integer :: t_temp_idx
   integer :: q_temp_idx
   integer :: u_temp_idx
   integer :: v_temp_idx

   integer :: t_stemp_idx
   integer :: q_stemp_idx
   integer :: u_stemp_idx
   integer :: v_stemp_idx
!--wy      



!  indices for fields in the physics buffer
   integer  ::    cld_idx          = 0    
   integer  ::    icwmrdp_idx      = 0     
   integer  ::    rprddp_idx       = 0    
   integer  ::    fracis_idx       = 0   
   integer  ::    nevapr_dpcu_idx  = 0    
!<songxl 2014-11-20-----------
   integer  ::    dgnum_idx        = 0
   integer  ::    lambdadpcu_idx   = 0
   integer  ::    mudpcu_idx       = 0
   integer  ::    icimrdp_idx      = 0


   logical :: old_snow  = .true.   ! set true to use old estimate of snow production in zm_conv_evap
                                   ! set false to use snow production from zm
                                   ! microphysics
   integer :: nmodes
   integer :: nbulk

   type(zm_aero_t), allocatable :: aero(:)   ! object contains information about the aerosols
!>songxl 2014-11-20-----------

!=========================================================================================
contains
!=========================================================================================

!================================================================================================

subroutine zm_conv_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  use physics_buffer, only : pbuf_add_field, dtype_r8

  implicit none

  integer idx

! Flux of precipitation from deep convection (kg/m2/s)
   call pbuf_add_field('DP_FLXPRC','global',dtype_r8,(/pcols,pverp/),dp_flxprc_idx) 

! Flux of snow from deep convection (kg/m2/s) 
   call pbuf_add_field('DP_FLXSNW','global',dtype_r8,(/pcols,pverp/),dp_flxsnw_idx) 

! deep gbm cloud liquid water (kg/kg)
   call pbuf_add_field('DP_CLDLIQ','global',dtype_r8,(/pcols,pver/), dp_cldliq_idx)  

! deep gbm cloud liquid water (kg/kg)    
   call pbuf_add_field('DP_CLDICE','global',dtype_r8,(/pcols,pver/), dp_cldice_idx)  

!<songxl 2014-11-20---------
   ! detrained convective cloud water mixing ratio.
   call pbuf_add_field('DLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dlfzm_idx)
   ! detrained convective cloud ice mixing ratio.
   call pbuf_add_field('DIFZM', 'physpkg', dtype_r8, (/pcols,pver/), difzm_idx)

   if (zm_microp) then
      ! Only add the number conc fields if the microphysics is active.

      ! detrained convective cloud water num concen.
      call pbuf_add_field('DNLFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnlfzm_idx)
      ! detrained convective cloud ice num concen.
      call pbuf_add_field('DNIFZM', 'physpkg', dtype_r8, (/pcols,pver/), dnifzm_idx)
   end if

!>songxl 2014-11-20---------

!++wy
   !if (zm_stc) then
   !    call pbuf_add_field('T_TEMP','global',dtype_r8,(/pcols,pver,6/), t_temp_idx)
   !    call pbuf_add_field('Q_TEMP','global',dtype_r8,(/pcols,pver,6/), q_temp_idx)
   !    call pbuf_add_field('U_TEMP','global',dtype_r8,(/pcols,pver,6/), u_temp_idx)
   !    call pbuf_add_field('V_TEMP','global',dtype_r8,(/pcols,pver,6/), v_temp_idx)

   !    call pbuf_add_field('T_STEMP','global',dtype_r8,(/pcols,pver,6/), t_stemp_idx)
   !    call pbuf_add_field('Q_STEMP','global',dtype_r8,(/pcols,pver,6/), q_stemp_idx)
   !    call pbuf_add_field('U_STEMP','global',dtype_r8,(/pcols,pver,6/), u_stemp_idx)
   !    call pbuf_add_field('V_STEMP','global',dtype_r8,(/pcols,pver,6/), v_stemp_idx)
   !end if
!--wy

end subroutine zm_conv_register

!=========================================================================================

subroutine zm_conv_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: outfld, addfld, add_default, phys_decomp
  use ppgrid,         only: pcols, pver
  use zm_conv,        only: zm_convi
  use pmgrid,         only: plev,plevp
  use spmd_utils,     only: masterproc
  use error_messages, only: alloc_err	
  use phys_control,   only: phys_deepconv_pbl, phys_getopts, cam_physpkg_is
  use physics_buffer, only: pbuf_get_index
!<songxl 2014-11-20------
  use zm_microphysics, only: zm_mphyi
!>songxl 2014-11-20------

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces


  logical :: no_deep_pbl    ! if true, no deep convection in PBL
  integer  limcnv           ! top interface level limit for convection
  integer k, istat
  logical :: history_budget ! output tendencies and state variables for CAM4
                            ! temperature, water vapor, cloud ice and cloud
                            ! liquid budgets.
  integer :: history_budget_histfile_num ! output history file number for budget fields

!<songxl 2014-11-20----------
  ! Aerosols
  integer :: i
  character(len=*), parameter :: routine = 'zm_conv_init' 
!  integer :: iaer, l, m
!  integer :: nspecmx   ! max number of species in a mode

!  character(len=20), allocatable :: aername(:)
!  character(len=32) :: str32

!  real(r8) :: sigmag, dgnumlo, dgnumhi
!  real(r8) :: alnsg
 
 
!>songxl 2014-11-20----------

!
! Allocate space for arrays private to this module
!
     allocate( mu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'mu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( eu(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'eu', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( du(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'du', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( md(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'md', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( ed(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'ed', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dp(pcols,pver,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'dp', &
                      pcols*pver*((endchunk-begchunk)+1) )
     allocate( dsubcld(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'dsubcld', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( jt(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'jt', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( maxg(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'maxg', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( ideep(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'ideep', &
                      pcols*((endchunk-begchunk)+1) )
     allocate( lengath(begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'zm_conv_tend', 'lengath', &
                      ((endchunk-begchunk)+1) )


! 
! Register fields with the output buffer
!


    call addfld ('PRECZ   ','m/s     ',1,    'A','total precipitation from ZM convection',        phys_decomp)
    call addfld ('ZMDT    ','K/s     ',pver, 'A','T tendency - Zhang-McFarlane moist convection', phys_decomp)
    call addfld ('ZMDQ    ','kg/kg/s ',pver, 'A','Q tendency - Zhang-McFarlane moist convection', phys_decomp)
    call addfld ('ZMDICE ','kg/kg/s ',pver, 'A','Cloud ice tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('ZMDLIQ ','kg/kg/s ',pver, 'A','Cloud liq tendency - Zhang-McFarlane convection',phys_decomp)
    call addfld ('EVAPTZM ','K/s     ',pver, 'A','T tendency - Evaporation/snow prod from Zhang convection',phys_decomp)
    call addfld ('FZSNTZM ','K/s     ',pver, 'A','T tendency - Rain to snow conversion from Zhang convection',phys_decomp)
    call addfld ('EVSNTZM ','K/s     ',pver, 'A','T tendency - Snow to rain prod from Zhang convection',phys_decomp)
    call addfld ('EVAPQZM ','kg/kg/s ',pver, 'A','Q tendency - Evaporation from Zhang-McFarlane moist convection',phys_decomp)
    
    call addfld ('ZMFLXPRC','kg/m2/s ',pverp, 'A','Flux of precipitation from ZM convection'       ,phys_decomp)
    call addfld ('ZMFLXSNW','kg/m2/s ',pverp, 'A','Flux of snow from ZM convection'                ,phys_decomp)
    call addfld ('ZMNTPRPD','kg/kg/s ',pver , 'A','Net precipitation production from ZM convection',phys_decomp)
    call addfld ('ZMNTSNPD','kg/kg/s ',pver , 'A','Net snow production from ZM convection'         ,phys_decomp)
    call addfld ('ZMEIHEAT','W/kg'    ,pver , 'A','Heating by ice and evaporation in ZM convection',phys_decomp)
    
    call addfld ('CMFMCDZM','kg/m2/s ',pverp,'A','Convection mass flux from ZM deep ',phys_decomp)
    call addfld ('PRECCDZM','m/s     ',1,    'A','Convective precipitation rate from ZM deep',phys_decomp)

    call addfld ('PCONVB','Pa'    ,1 , 'A','convection base pressure',phys_decomp)
    call addfld ('PCONVT','Pa'    ,1 , 'A','convection top  pressure',phys_decomp)

    call addfld ('CAPE',   'J/kg',       1, 'A', 'Convectively available potential energy', phys_decomp)
!<songxl 2011-09-20----------------------------------
    call addfld ('ZMDU    ','m/s2    ',pver , 'A','U tendency - zhang CMT scheme', phys_decomp)
    call addfld ('ZMDV    ','m/s2    ',pver , 'A','V tendency - zhang CMT scheme', phys_decomp)
    call addfld ('ZMPU    ','m/s2    ',pver , 'A','dpdx - zhang CMT scheme', phys_decomp)
    call addfld ('ZMPV    ','m/s2    ',pver , 'A','dpdy - zhang CMT scheme', phys_decomp)
!>songxl 2011-09-20----------------------------------
    !++wy
    if (zm_stc) then
        call addfld ('CAPE_AVG',   'J/kg',       1, 'A', 'Convectively available potential energy', phys_decomp)
    end if
    !--wy
    call addfld ('FREQZM ','fraction  ',1  ,'A', 'Fractional occurance of ZM convection',phys_decomp) 

    call addfld ('ZMMTT ', 'K/s',     pver, 'A', 'T tendency - ZM convective momentum transport',phys_decomp)
    call addfld ('ZMMTU',  'm/s2',    pver, 'A', 'U tendency - ZM convective momentum transport',  phys_decomp)
    call addfld ('ZMMTV',  'm/s2',    pver, 'A', 'V tendency - ZM convective momentum transport',  phys_decomp)

    call addfld ('ZMMU',   'kg/m2/s', pver, 'A', 'ZM convection updraft mass flux',   phys_decomp)
    call addfld ('ZMMD',   'kg/m2/s', pver, 'A', 'ZM convection downdraft mass flux', phys_decomp)

    call addfld ('ZMUPGU', 'm/s2',    pver, 'A', 'zonal force from ZM updraft pressure gradient term',       phys_decomp)
    call addfld ('ZMUPGD', 'm/s2',    pver, 'A', 'zonal force from ZM downdraft pressure gradient term',     phys_decomp)
    call addfld ('ZMVPGU', 'm/s2',    pver, 'A', 'meridional force from ZM updraft pressure gradient term',  phys_decomp)
    call addfld ('ZMVPGD', 'm/s2',    pver, 'A', 'merdional force from ZM downdraft pressure gradient term', phys_decomp)

    call addfld ('ZMICUU', 'm/s',     pver, 'A', 'ZM in-cloud U updrafts',      phys_decomp)
    call addfld ('ZMICUD', 'm/s',     pver, 'A', 'ZM in-cloud U downdrafts',    phys_decomp)
    call addfld ('ZMICVU', 'm/s',     pver, 'A', 'ZM in-cloud V updrafts',      phys_decomp)
    call addfld ('ZMICVD', 'm/s',     pver, 'A', 'ZM in-cloud V downdrafts',    phys_decomp)

!<songxl 2014-11-20---------

    if (zm_microp) then
       call addfld ('DCAPE', 'J/kg', 1, 'A', 'CAPE change due to freezing heating', phys_decomp)
       call add_default ('DCAPE',1, ' ')

       call addfld ('CLDLIQZM','g/m3'    ,pver, 'A','Cloud liquid water - ZM convection',phys_decomp)
       call addfld ('CLDICEZM','g/m3'    ,pver, 'A','Cloud ice water - ZM convection',phys_decomp)
       call addfld ('CLIQSNUM','1'       ,pver, 'A','Cloud liquid water sample number - ZM convection',phys_decomp)
       call addfld ('CICESNUM','1'       ,pver, 'A','Cloud ice water sample number - ZM convection',phys_decomp)
       call addfld ('QRAINZM' ,'g/m3'    ,pver, 'A','rain water - ZM convection',phys_decomp)
       call addfld ('QSNOWZM' ,'g/m3'    ,pver, 'A','snow - ZM convection',phys_decomp)
       call addfld ('CRAINNUM','1'       ,pver, 'A','Cloud rain water sample number - ZM convection',phys_decomp)
       call addfld ('CSNOWNUM','1'       ,pver, 'A','Cloud snow sample number - ZM convection',phys_decomp)

       call addfld ('DIFZM','kg/kg/s '   ,pver, 'A','Detrained ice water from ZM convection',phys_decomp)
       call addfld ('DLFZM','kg/kg/s '   ,pver, 'A','Detrained liquid water from ZM convection',phys_decomp)
       call addfld ('DNIFZM','1/kg/s '   ,pver, 'A','Detrained ice water num concen from ZM convection',phys_decomp)
       call addfld ('DNLFZM','1/kg/s '   ,pver, 'A','Detrained liquid water num concen from ZM convection',phys_decomp)
       call addfld ('WUZM','m/s'         ,pver, 'A','vertical velocity - ZM convection',phys_decomp)
       call addfld ('WUZMSNUM','1'       ,pver, 'A','vertical velocity sample number - ZM convection',phys_decomp)

       call addfld ('QNLZM','1/m3'       ,pver, 'A','Cloud liquid water number concen - ZM convection',phys_decomp)
       call addfld ('QNIZM','1/m3'       ,pver, 'A','Cloud ice number concen - ZM convection',phys_decomp)
       call addfld ('QNRZM','1/m3'       ,pver, 'A','Cloud rain water number concen - ZM convection',phys_decomp)
       call addfld ('QNSZM','1/m3'       ,pver, 'A','Cloud snow number concen - ZM convection',phys_decomp)

       call addfld ('FRZZM','1/s'       ,pver, 'A','mass tendency due to freezing - ZM convection',phys_decomp)    

       call addfld ('PSZM','Pa'          ,1,    'A','surface pressure on convective grid point', phys_decomp)
       call addfld ('AUTOL_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to autoconversion of droplets to rain',phys_decomp)
       call addfld ('ACCRL_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to accretion of droplets by rain',phys_decomp)
       call addfld ('BERGN_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to Bergeron process',phys_decomp)
       call addfld ('FHTIM_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to immersion freezing',phys_decomp)
       call addfld ('FHTCT_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to contact freezing',phys_decomp)
       call addfld ('FHML_M'  ,'kg/kg/m' ,pver, 'A','mass tendency due to homogeneous freezing of droplet',phys_decomp)
       call addfld ('HMPI_M'  ,'kg/kg/m' ,pver, 'A','mass tendency due to HM process',phys_decomp)
       call addfld ('ACCSL_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to accretion of droplet by snow',phys_decomp)
       call addfld ('DLF_M'   ,'kg/kg/m' ,pver, 'A','mass tendency due to detrainment of droplet',phys_decomp)
       call addfld ('COND_M'  ,'kg/kg/m' ,pver, 'A','mass tendency due to condensation',phys_decomp)

       call addfld ('AUTOL_N' ,'1/kg/m' ,pver, 'A','num tendency due to autoconversion of droplets to rain',phys_decomp)
       call addfld ('ACCRL_N' ,'1/kg/m' ,pver, 'A','num tendency due to accretion of droplets by rain',phys_decomp)
       call addfld ('BERGN_N' ,'1/kg/m' ,pver, 'A','num tendency due to Bergeron process',phys_decomp)
       call addfld ('FHTIM_N' ,'1/kg/m' ,pver, 'A','num tendency due to immersion freezing',phys_decomp)
       call addfld ('FHTCT_N' ,'1/kg/m' ,pver, 'A','num tendency due to contact freezing',phys_decomp)
       call addfld ('FHML_N'  ,'1/kg/m' ,pver, 'A','num tendency due to homogeneous freezing of droplet',phys_decomp)
       call addfld ('ACCSL_N' ,'1/kg/m' ,pver, 'A','num tendency due to accretion of droplet by snow',phys_decomp)
       call addfld ('ACTIV_N' ,'1/kg/m' ,pver, 'A','num tendency due to droplets activation',phys_decomp)
       call addfld ('DLF_N'   ,'1/kg/m' ,pver, 'A','num tendency due to detrainment of droplet',phys_decomp)

       call addfld ('AUTOI_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to autoconversion of ice to snow',phys_decomp)
       call addfld ('ACCSI_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to accretion of ice by snow',phys_decomp)
       call addfld ('DIF_M'   ,'kg/kg/m' ,pver, 'A','mass tendency due to detrainment of cloud ice',phys_decomp)
       call addfld ('DEPOS_M' ,'kg/kg/m' ,pver, 'A','mass tendency due to deposition',phys_decomp)

       call addfld ('NUCLI_N' ,'1/kg/m' ,pver, 'A','num tendency due to ice nucleation',phys_decomp)
       call addfld ('AUTOI_N' ,'1/kg/m' ,pver, 'A','num tendency due to autoconversion of ice to snow',phys_decomp)
       call addfld ('ACCSI_N' ,'1/kg/m' ,pver, 'A','num tendency due to accretion of ice by snow',phys_decomp)
       call addfld ('HMPI_N'  ,'1/kg/s' ,pver, 'A','num tendency due to HM process',phys_decomp)
       call addfld ('DIF_N'   ,'1/kg/m' ,pver, 'A','num tendency due to detrainment of cloud ice',phys_decomp)

       call addfld ('TRSPC_M' ,'kg/kg/m',pver, 'A','mass tendency of droplets due to convective transport',phys_decomp)
       call addfld ('TRSPC_N' ,'1/kg/m' ,pver, 'A','num tendency of droplets due to convective transport',phys_decomp)
       call addfld ('TRSPI_M' ,'kg/kg/m',pver, 'A','mass tendency of ice crystal due to convective transport',phys_decomp)
       call addfld ('TRSPI_N' ,'1/kg/m' ,pver, 'A','num tendency of ice crystal due to convective transport',phys_decomp)

    call addfld ('DUSNUM'  ,'#'      ,pver, 'A','ZM convection mass detrainment rate sample number',phys_decomp)

    call addfld ('AUTOL_SN' ,'#' ,pver, 'A','smaple num of tendency due to autoconversion of droplets to rain',phys_decomp)
    call addfld ('ACCRL_SN' ,'#' ,pver, 'A','smaple num of tendency due to accretion of droplets by rain',phys_decomp)
    call addfld ('BERGN_SN' ,'#' ,pver, 'A','sample num of tendency due to Bergeron process',phys_decomp)
    call addfld ('FHTIM_SN' ,'#' ,pver, 'A','smaple num of tendency due to immersion freezing',phys_decomp)
    call addfld ('FHTCT_SN' ,'#' ,pver, 'A','sample num of tendency due to contact freezing',phys_decomp)
    call addfld ('FHML_SN'  ,'#' ,pver, 'A','sample num of tendency due to homogeneous freezing of droplet',phys_decomp)
    call addfld ('ACCSL_SN' ,'#' ,pver, 'A','smaple num of tendency due to accretion of droplet by snow',phys_decomp)
    call addfld ('ACTIV_SN' ,'#' ,pver, 'A','sample num of tendency due to droplets activation',phys_decomp)
    call addfld ('DLF_SN'   ,'#' ,pver, 'A','sample num of tendency due to detrainment of droplet',phys_decomp)
    call addfld ('NUCLI_SN' ,'#' ,pver, 'A','sample num of tendency due to ice nucleation',phys_decomp)
    call addfld ('AUTOI_SN' ,'#' ,pver, 'A','sample num of tendency due to autoconversion of ice to snow',phys_decomp)
    call addfld ('ACCSI_SN' ,'#' ,pver, 'A','sample num of tendency due to accretion of ice by snow',phys_decomp)
    call addfld ('HMPI_SN'  ,'#' ,pver, 'A','sample num of tendency due to HM process',phys_decomp)
    call addfld ('DIF_SN'   ,'#' ,pver, 'A','sample num of tendency due to detrainment of cloud ice',phys_decomp)
    call addfld ('TRSPC_SN' ,'#' ,pver, 'A','smaple num of droplet num tendency due to convective transport',phys_decomp)
    call addfld ('TRSPI_SN' ,'#' ,pver, 'A','smaple num of ice crystal num tendency due to convective transport',phys_decomp)

    call addfld ('COND_SN'  ,'#' ,pver, 'A','sample num of mass tendency due to condensation',phys_decomp)
    call addfld ('DEPOS_SN' ,'#' ,pver, 'A','sample num of mass tendency due to deposition',phys_decomp)

    call addfld ('FRZZM_SN','#' ,pver, 'A','sample num of mass tendency due to freezing -ZM convection',phys_decomp)
    call addfld ('ZMDT_SN'  ,'#' ,pver, 'A','sample num of ZMDT',phys_decomp)
    call addfld ('PRECZ_SN' ,'#' ,1,    'A','sample num of convective precipitation rate from ZM deep',phys_decomp)
    call addfld ('PCT_SN'   ,'#' ,1,    'A','sample num of convective cld top pressure',phys_decomp)

       call add_default ('CLDLIQZM', 1, ' ')
       call add_default ('CLDICEZM', 1, ' ')
       call add_default ('CLIQSNUM', 1, ' ')
       call add_default ('CICESNUM', 1, ' ')
       call add_default ('DIFZM',    1, ' ')
       call add_default ('DLFZM',    1, ' ')
       call add_default ('DNIFZM',   1, ' ')
       call add_default ('DNLFZM',   1, ' ')
       call add_default ('WUZM',     1, ' ')
       call add_default ('QRAINZM',  1, ' ')
       call add_default ('QSNOWZM',  1, ' ')
       call add_default ('CRAINNUM', 1, ' ')
       call add_default ('CSNOWNUM', 1, ' ')
       call add_default ('QNLZM',    1, ' ')
       call add_default ('QNIZM',    1, ' ')
       call add_default ('QNRZM',    1, ' ')
       call add_default ('QNSZM',    1, ' ')
       call add_default ('FRZZM',   1, ' ')
       call add_default ('PSZM',   1, ' ')

!       call add_default ('AUTOL_M',  1, ' ')
!       call add_default ('ACCRL_M',  1, ' ')
!       call add_default ('BERGN_M',  1, ' ')
!       call add_default ('FHTIM_M',  1, ' ')
!       call add_default ('FHTCT_M',  1, ' ')
!       call add_default ('FHML_M',   1, ' ')
!       call add_default ('HMPI_M',   1, ' ')
!       call add_default ('ACCSL_M',  1, ' ')
!       call add_default ('DLF_M',    1, ' ')
!       call add_default ('COND_M',   1, ' ')

!       call add_default ('AUTOL_N',  1, ' ')
!       call add_default ('ACCRL_N',  1, ' ')
!       call add_default ('BERGN_N',  1, ' ')
!       call add_default ('FHTIM_N',  1, ' ')
!       call add_default ('FHTCT_N',  1, ' ')
!       call add_default ('FHML_N',   1, ' ')
!       call add_default ('ACCSL_N',  1, ' ')
!       call add_default ('ACTIV_N',  1, ' ')
!       call add_default ('DLF_N',    1, ' ')
!       call add_default ('AUTOI_M',  1, ' ')
!       call add_default ('ACCSI_M',  1, ' ')
!       call add_default ('DIF_M',    1, ' ')
!       call add_default ('DEPOS_M',  1, ' ')

!       call add_default ('NUCLI_N',  1, ' ')
!       call add_default ('AUTOI_N',  1, ' ')
!       call add_default ('ACCSI_N',  1, ' ')
!       call add_default ('HMPI_N',   1, ' ')
!       call add_default ('DIF_N',    1, ' ')

!       call add_default ('TRSPC_M',  1, ' ')
!       call add_default ('TRSPC_N',  1, ' ')
!       call add_default ('TRSPI_M',  1, ' ')
!       call add_default ('TRSPI_N',  1, ' ')

!    call add_default ('DUSNUM',  1, ' ')

!    call add_default ('AUTOL_SN',  1, ' ')
!    call add_default ('ACCRL_SN',  1, ' ')
!    call add_default ('BERGN_SN',  1, ' ')
!    call add_default ('FHTIM_SN',  1, ' ')
!    call add_default ('FHTCT_SN',  1, ' ')
!    call add_default ('FHML_SN',   1, ' ')
!    call add_default ('ACCSL_SN',  1, ' ')
!    call add_default ('ACTIV_SN',  1, ' ')
!    call add_default ('DLF_SN',    1, ' ')
!    call add_default ('NUCLI_SN',  1, ' ')
!    call add_default ('AUTOI_SN',  1, ' ')
!    call add_default ('ACCSI_SN',  1, ' ')
!    call add_default ('HMPI_SN',   1, ' ')
!    call add_default ('DIF_SN',    1, ' ')
!    call add_default ('TRSPC_SN',  1, ' ')
!    call add_default ('TRSPI_SN',  1, ' ')

!    call add_default ('COND_SN',   1, ' ')
!    call add_default ('DEPOS_SN',  1, ' ')

!    call add_default ('FRZZM_SN', 1, ' ')
!    call add_default ('ZMDT_SN',   1, ' ')
!    call add_default ('PRECZ_SN',  1, ' ')
!    call add_default ('PCT_SN',    1, ' ')

    end if
  
!>songxl 2014-11-20---------    

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num)

    if ( history_budget ) then
       call add_default('EVAPTZM  ', history_budget_histfile_num, ' ')
       call add_default('EVAPQZM  ', history_budget_histfile_num, ' ')
       call add_default('ZMDT     ', history_budget_histfile_num, ' ')
       call add_default('ZMDQ     ', history_budget_histfile_num, ' ')
       call add_default('ZMDLIQ   ', history_budget_histfile_num, ' ')
       call add_default('ZMDICE   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam4') .or. cam_physpkg_is('cam5') ) then
          call add_default('ZMMTT    ', history_budget_histfile_num, ' ')
       end if

    end if
!
! Limit deep convection to regions below 40 mb
! Note this calculation is repeated in the shallow convection interface
!
    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,plev
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = plevp
    end if
    
    if (masterproc) then
       write(iulog,*)'ZM_CONV_INIT: Deep convection will be capped at intfc ',limcnv, &
            ' which is ',pref_edge(limcnv),' pascals'
       !++wy
       if (zm_microp) then
           write(iulog,*) 'ZM_CONV_INIT: Convective microphysics is turned on'
       end if
       if (zm_stc) then
           write(iulog,*) 'ZM_CONV_INIT: Stochastic convection is turned on'
       end if
       !--wy
    end if

    no_deep_pbl = phys_deepconv_pbl()

!<songxl 2014-11-20--------------
    call zm_convi(limcnv,no_deep_pbl_in = no_deep_pbl)

    lambdadpcu_idx  = pbuf_get_index('LAMBDADPCU')
    mudpcu_idx      = pbuf_get_index('MUDPCU')   
    icimrdp_idx     = pbuf_get_index('ICIMRDP') 
!>songxl 2014-11-20--------------
    cld_idx         = pbuf_get_index('CLD')
    icwmrdp_idx     = pbuf_get_index('ICWMRDP')
    rprddp_idx      = pbuf_get_index('RPRDDP')
    fracis_idx      = pbuf_get_index('FRACIS')
    nevapr_dpcu_idx = pbuf_get_index('NEVAPR_DPCU')
    prec_dp_idx     = pbuf_get_index('PREC_DP')
    snow_dp_idx     = pbuf_get_index('SNOW_DP')

!<songxl 2014-11-20----------
    ! Initialization for the microphysics
    if (zm_microp) then
       old_snow = (.not. zm_microp)
       call zm_mphyi()

       ! Initialize the aerosol object with data from the modes/species
       ! affecting climate,
       ! i.e., the list index is hardcoded to 0.

       call rad_cnst_get_info(0, nmodes=nmodes, naero=nbulk)

       allocate(aero(begchunk:endchunk))

       do i = begchunk, endchunk
          call zm_aero_init(nmodes, nbulk, aero(i))
       end do

       if (nmodes > 0) then

          dgnum_idx = pbuf_get_index('DGNUM')

       else if (nbulk > 0) then

          ! This call is needed to allow running the ZM microphysics with the
          ! cam4 physics package.
          call ndrop_bam_init()

       end if

     end if ! zmconv_microp

    !-------------------------------------------------------------------------------------
    contains
    !-------------------------------------------------------------------------------------

    subroutine zm_aero_init(nmodes, nbulk, aero)

       ! Initialize the zm_aero_t object for modal aerosols

       integer,         intent(in)  :: nmodes
       integer,         intent(in)  :: nbulk
       type(zm_aero_t), intent(out) :: aero

       integer :: iaer, l, m
       integer :: nspecmx   ! max number of species in a mode

       character(len=20), allocatable :: aername(:)
       character(len=32) :: str32

       real(r8) :: sigmag, dgnumlo, dgnumhi
       real(r8) :: alnsg
       !----------------------------------------------------------------------------------

       aero%nmodes = nmodes
       aero%nbulk  = nbulk

       if (nmodes > 0) then

          ! Initialize the modal aerosol information

          aero%scheme = 'modal'

          ! Get number of species in each mode, and find max.
          allocate(aero%nspec(aero%nmodes))
          nspecmx = 0
          do m = 1, aero%nmodes

             call rad_cnst_get_info(0, m, nspec=aero%nspec(m), mode_type=str32)

             nspecmx = max(nspecmx, aero%nspec(m))

             ! save mode index for specified mode types
             select case (trim(str32))
             case ('accum')
                aero%mode_accum_idx = m
             case ('aitken')
                aero%mode_aitken_idx = m
             case ('coarse')
                aero%mode_coarse_idx = m
             end select

          end do

          ! Check that required mode types were found
          if (aero%mode_accum_idx == -1 .or. aero%mode_aitken_idx == -1 .or. aero%mode_coarse_idx == -1) then
             write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
                aero%mode_accum_idx, aero%mode_aitken_idx, aero%mode_coarse_idx
             call endrun(routine//': ERROR required mode type not found')
          end if

          ! find indices for the dust and seasalt species in the coarse mode
          do l = 1, aero%nspec(aero%mode_coarse_idx)
             call rad_cnst_get_info(0, aero%mode_coarse_idx, l, spec_type=str32)
             select case (trim(str32))
             case ('dust')
                aero%coarse_dust_idx = l
             case ('seasalt')
                aero%coarse_nacl_idx = l
             end select
          end do
          ! Check that required modal specie types were found
          if (aero%coarse_dust_idx == -1 .or. aero%coarse_nacl_idx == -1) then
             write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
                aero%coarse_dust_idx, aero%coarse_nacl_idx
             call endrun(routine//': ERROR required mode-species type not found')
          end if

          allocate( &
             aero%num_a(nmodes), &
             aero%mmr_a(nspecmx,nmodes), &
             aero%numg_a(pcols,pver,nmodes), &
             aero%mmrg_a(pcols,pver,nspecmx,nmodes), &
             aero%voltonumblo(nmodes), &
             aero%voltonumbhi(nmodes), &
             aero%specdens(nspecmx,nmodes), &
             aero%spechygro(nspecmx,nmodes), &
             aero%dgnum(nmodes), &
             aero%dgnumg(pcols,pver,nmodes) )

          do m = 1, nmodes

             ! Properties of modes
             call rad_cnst_get_mode_props(0, m, &
                sigmag=sigmag, dgnumlo=dgnumlo, dgnumhi=dgnumhi)

             alnsg               = log(sigmag)
             aero%voltonumblo(m) = 1._r8 / ( (pi/6._r8)*(dgnumlo**3._r8)*exp(4.5_r8*alnsg**2._r8) )
             aero%voltonumbhi(m) = 1._r8 / ( (pi/6._r8)*(dgnumhi**3._r8)*exp(4.5_r8*alnsg**2._r8) )

             ! save sigmag of aitken mode
             if (m == aero%mode_aitken_idx) aero%sigmag_aitken = sigmag

             ! Properties of modal species
             do l = 1, aero%nspec(m)
                call rad_cnst_get_aer_props(0, m, l, density_aer=aero%specdens(l,m), &
                   hygro_aer=aero%spechygro(l,m))
             end do
          end do

       else if (nbulk > 0) then

          aero%scheme = 'bulk'

          ! Props needed for BAM number concentration calcs.
          allocate( &
             aername(nbulk),                   &
             aero%num_to_mass_aer(nbulk),      &
             aero%mmr_bulk(nbulk),             &
             aero%mmrg_bulk(pcols,plev,nbulk)  )

          do iaer = 1, aero%nbulk
             call rad_cnst_get_aer_props(0, iaer, &
                aername         = aername(iaer), &
                num_to_mass_aer = aero%num_to_mass_aer(iaer) )

             ! Look for sulfate aerosol in this list (Bulk aerosol only)
             if (trim(aername(iaer)) == 'SULFATE') aero%idxsul = iaer
             if (trim(aername(iaer)) == 'DUST1')   aero%idxdst1 = iaer
             if (trim(aername(iaer)) == 'DUST2')   aero%idxdst2 = iaer
             if (trim(aername(iaer)) == 'DUST3')   aero%idxdst3 = iaer
             if (trim(aername(iaer)) == 'DUST4')   aero%idxdst4 = iaer
             if (trim(aername(iaer)) == 'BCPHI')   aero%idxbcphi = iaer
          end do

       end if

    end subroutine zm_aero_init
!>songxl 2014-11-20-------------

end subroutine zm_conv_init
!=========================================================================================
!subroutine zm_conv_tend(state, ptend, tdt)

subroutine zm_conv_tend(pblh    ,mcon    ,cme     , &
!<songxl 2014-11-20---
!     tpert   ,dlf     ,pflx    ,zdu      , &
     tpert   ,pflx    ,zdu      , &
!>songxl 2014-11-20---
     rliq    ,rice    ,  &
     ztodt   , &
     jctop   ,jcbot , &
     state   ,ptend_all   ,landfrac,  pbuf)
  

   use cam_history,   only: outfld
   use physics_types, only: physics_state, physics_ptend
   use physics_types, only: physics_ptend_init, physics_update
   use physics_types, only: physics_state_copy, physics_state_dealloc
   use physics_types, only: physics_ptend_sum, physics_ptend_dealloc
    
   !++wy
   !use phys_grid,     only: get_lat_p, get_lon_p
   use phys_grid,     only: get_lat_p, get_lon_p, get_area_all_p    
   !--wy

   use time_manager,  only: get_nstep, is_first_step
   use physics_buffer, only : pbuf_get_field, physics_buffer_desc, pbuf_old_tim_idx
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use check_energy,  only: check_energy_chng
   use physconst,     only: gravit
   use phys_control,  only: cam_physpkg_is
   !++wy
   !use spatial_avg_tquv, only: t_spatial_avg, q_spatial_avg, u_spatial_avg, v_spatial_avg
   !--wy    

   ! Arguments
!<songxl 2014-11-20--------
!   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_state), target, intent(in ) :: state      ! Physics state variables
!>songxl 2014-11-20---------
   type(physics_ptend), intent(out)   :: ptend_all      ! individual parameterization tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                       ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblh(pcols)                 ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: landfrac(pcols)             ! RBN - Landfrac 

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
!songxl 2014-11-20   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: rice(pcols) ! reserved ice (not yet in cldice) for energy integrals

   ! Local variables
!<songxl 11-20---------
!   integer :: i,k,m
   integer :: i,k,l,m
!>songxl 11-20--------
   integer :: ilon                      ! global longitude index of a column
   integer :: ilat                      ! global latitude index of a column
   integer :: nstep
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.
   integer :: lchnk                   ! chunk identifier
   integer :: ncol                    ! number of atmospheric columns
   integer :: itim                    ! for physics buffer fields

   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: ntprprd(pcols,pver)    ! evap outfld: net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)    ! evap outfld: net snow production in layer
   real(r8) :: tend_s_snwprd  (pcols,pver) ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver) ! Heating rate of evap/melting of snow
   real(r8) :: fake_dpdry(pcols,pver) ! used in convtran call

   ! physics types
   type(physics_state) :: state1        ! locally modify for evaporation to use, not returned
!   type(physics_ptend) :: ptend_loc     ! package tendencies
   type(physics_ptend),target :: ptend_loc     ! package tendencies

   ! physics buffer fields
   real(r8), pointer, dimension(:)   :: prec         ! total precipitation
   real(r8), pointer, dimension(:)   :: snow         ! snow from ZM convection 
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: ql           ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd         ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   real(r8), pointer, dimension(:,:) :: evapcdp      ! Evaporation of deep convective precipitation
   real(r8), pointer, dimension(:,:) :: flxprec      ! Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: flxsnow      ! Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8), pointer, dimension(:,:) :: dp_cldliq
   real(r8), pointer, dimension(:,:) :: dp_cldice

!<songxl 2011-09-20------------------------------
   real(r8) :: zmpu(pcols,pver)                ! Zhang convective dp/dx
   real(r8) :: zmpv(pcols,pver)                ! Zhang convective dp/dy
!>songxl 2011-09-20------------------------------

   !++wy
   real(r8), pointer, dimension(:,:,:) :: t_temp
   real(r8), pointer, dimension(:,:,:) :: q_temp
   real(r8), pointer, dimension(:,:,:) :: u_temp
   real(r8), pointer, dimension(:,:,:) :: v_temp

   real(r8), pointer, dimension(:,:,:) :: t_stemp
   real(r8), pointer, dimension(:,:,:) :: q_stemp
   real(r8), pointer, dimension(:,:,:) :: u_stemp
   real(r8), pointer, dimension(:,:,:) :: v_stemp
   !--wy

!<songxl 2014-11-20-----
   real(r8), pointer :: dlf(:,:)    ! detrained convective cloud water mixing ratio.
   real(r8), pointer :: dif(:,:)    ! detrained convective cloud ice mixing ratio.
   real(r8), pointer :: dnlf(:,:)   ! detrained convective cloud water num concen.
   real(r8), pointer :: dnif(:,:)   ! detrained convective cloud ice num concen.
   real(r8), pointer :: lambdadpcu(:,:) ! slope of cloud liquid size distr
   real(r8), pointer :: mudpcu(:,:)     ! width parameter of droplet size distr
   real(r8), pointer :: qi(:,:)         ! wg grid slice of cloud ice.
!>songxl 2014-11-20-----

   real(r8) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8) :: jcbot(pcols)  ! o row of base of cloud indices passed out.

   real(r8) :: pcont(pcols), pconb(pcols), freqzm(pcols)

   ! history output fields
   real(r8) :: cape(pcols)        ! w  convective available potential energy.
   !++wy
   real(r8) :: cape_avgd(pcols)   ! w  convective available potential energy.
   !--wy
   real(r8) :: mu_out(pcols,pver)
   real(r8) :: md_out(pcols,pver)

   ! used in momentum transport calculation
   real(r8) :: winds(pcols, pver, 2)
   real(r8) :: wind_tends(pcols, pver, 2)
   real(r8) :: pguall(pcols, pver, 2)
   real(r8) :: pgdall(pcols, pver, 2)
   real(r8) :: icwu(pcols,pver, 2)
   real(r8) :: icwd(pcols,pver, 2)
   real(r8) :: seten(pcols, pver)
   logical  :: l_windt(2)
   real(r8) :: tfinal1, tfinal2
   integer  :: ii

   logical  :: lq(pcnst)

   !++wy
   integer  :: timesteps
   real(r8) :: u_avgd(pcols, pver)                       ! spatial average of U
   real(r8) :: v_avgd(pcols, pver)                       ! spatial average of V
   real(r8) :: t_avgd(pcols, pver)                       ! spatial average of T
   real(r8) :: q_avgd(pcols, pver)                       ! spatial average of Q
   real(r8) :: area(pcols)                               ! column surface area

   integer :: kinv                                       ! to invert the vertical grid labelling
   integer :: n
   !integer :: a_flag = 0
   !integer :: a_flag = 3
   integer :: a_flag = 2
   !--wy

!<songxl 2014-11-20-----------
   real(r8) :: qice(pcols,pver)           ! convective cloud ice.
   real(r8) :: qliq(pcols,pver)           ! convective cloud liquid water.
   real(r8) :: cice_snum(pcols,pver)      ! convective cloud ice sample number.
   real(r8) :: cliq_snum(pcols,pver)      ! convective cloud liquid sample number.
   real(r8) :: qrain(pcols,pver)          ! convective rain water.
   real(r8) :: qsnow(pcols,pver)          ! convective snow.
   real(r8) :: crain_snum(pcols,pver)     ! convective rain water sample number.
   real(r8) :: csnow_snum(pcols,pver)     ! convective snow sample number.
   real(r8) :: wu(pcols,pver)             ! vertical velocity
   real(r8) :: wu_snum(pcols,pver)        ! vertical velocity sample number

   real(r8) :: sprd(pcols,pver)
   real(r8) :: frz(pcols,pver)            ! heating rate due to freezing
   real(r8) :: qni(pcols,pver)            ! convective cloud ice num concen.
   real(r8) :: qnl(pcols,pver)            ! convective cloud liquid water num concen.
   real(r8) :: qni_snum(pcols,pver)       ! convective cloud ice number sample number.
   real(r8) :: qnl_snum(pcols,pver)       ! convective cloud liquid number sample number.
   real(r8) :: qnr(pcols,pver)            ! convective rain water num concen.
   real(r8) :: qns(pcols,pver)            ! convective snow num concen.

   real(r8) dcape(pcols)           ! CAPE change due to freezing heating

   real(r8) autolm(pcols,pver)    !mass tendency due to autoconversion of droplets to rain
   real(r8) accrlm(pcols,pver)    !mass tendency due to accretion of droplets by rain
   real(r8) bergnm(pcols,pver)    !mass tendency due to Bergeron process
   real(r8) fhtimm(pcols,pver)    !mass tendency due to immersion freezing
   real(r8) fhtctm(pcols,pver)    !mass tendency due to contact freezing
   real(r8) fhmlm (pcols,pver)    !mass tendency due to homogeneous freezing
   real(r8) hmpim (pcols,pver)    !mass tendency due to HM process
   real(r8) accslm(pcols,pver)    !mass tendency due to accretion of droplets by snow
   real(r8) dlfm  (pcols,pver)    !mass tendency due to detrainment of droplet
   real(r8) cmel  (pcols,pver)    !mass tendency due to condensation

   real(r8) autoln(pcols,pver)    !num tendency due to autoconversion of droplets to rain
   real(r8) accrln(pcols,pver)    !num tendency due to accretion of droplets by rain
   real(r8) bergnn(pcols,pver)    !num tendency due to Bergeron process
   real(r8) fhtimn(pcols,pver)    !num tendency due to immersion freezing
   real(r8) fhtctn(pcols,pver)    !num tendency due to contact freezing
   real(r8) fhmln (pcols,pver)    !num tendency due to homogeneous freezing
   real(r8) accsln(pcols,pver)    !num tendency due to accretion of droplets by snow
   real(r8) activn(pcols,pver)    !num tendency due to droplets activation
   real(r8) dlfn  (pcols,pver)    !num tendency due to detrainment of droplet

   real(r8) autoim(pcols,pver)    !mass tendency due to autoconversion of cloud ice to snow
   real(r8) accsim(pcols,pver)    !mass tendency due to accretion of cloud ice by snow
   real(r8) difm  (pcols,pver)    !mass tendency due to detrainment of cloud ice
   real(r8) cmei  (pcols,pver)    !mass tendency due to deposition

   real(r8) nuclin(pcols,pver)    !num tendency due to ice nucleation
   real(r8) autoin(pcols,pver)    !num tendency due to autoconversion of cloud ice to snow
   real(r8) accsin(pcols,pver)    !num tendency due to accretion of cloud ice by snow
   real(r8) hmpin (pcols,pver)    !num tendency due to HM process
   real(r8) difn  (pcols,pver)    !num tendency due to detrainment of cloud ice

   real(r8) trspcm(pcols,pver)    !LWC tendency due to convective transport
   real(r8) trspcn(pcols,pver)    !droplet num tendency due to convective transport
   real(r8) trspim(pcols,pver)    !IWC tendency due to convective transport
   real(r8) trspin(pcols,pver)    !ice crystal num tendency due to convective transport

   real(r8) du_snum(pcols,pver)

   real(r8)  autoln_snum(pcols,pver)
   real(r8)  accrln_snum(pcols,pver)
   real(r8)  bergnn_snum(pcols,pver)
   real(r8)  fhtimn_snum(pcols,pver)
   real(r8)  fhtctn_snum(pcols,pver)
   real(r8)  fhmln_snum(pcols,pver)
   real(r8)  accsln_snum(pcols,pver)
   real(r8)  activn_snum(pcols,pver)
   real(r8)  dlfn_snum(pcols,pver)

   real(r8)  nuclin_snum(pcols,pver)
   real(r8)  autoin_snum(pcols,pver)
   real(r8)  accsin_snum(pcols,pver)
   real(r8)  hmpin_snum(pcols,pver)
   real(r8)  difn_snum(pcols,pver)

   real(r8)  cmel_snum(pcols,pver)
   real(r8)  cmei_snum(pcols,pver)

   real(r8)  trspcn_snum(pcols,pver)
   real(r8)  trspin_snum(pcols,pver)

   real(r8)  frz_snum(pcols,pver)
   real(r8)  zmdt_snum(pcols,pver)

   real(r8)  precz_snum(pcols)

   real(r8) pszm(pcols)
   real(r8) pcont_snum(pcols)
!>songxl 2014-11-20-----------

   !----------------------------------------------------------------------

   ! initialize

   !++wy
   u_avgd = 0._r8
   v_avgd = 0._r8
   t_avgd = 0._r8
   q_avgd = 0._r8
   !--wy

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ftem = 0._r8   
   mu_out(:,:) = 0._r8
   md_out(:,:) = 0._r8
   wind_tends(:ncol,:pver,:) = 0.0_r8

   call physics_state_copy(state,state1)             ! copy state to local state1.

   lq(:) = .FALSE.
   lq(1) = .TRUE.
!<songxl 2011-09-20--------------------
   !call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq)! initialize local ptend type
    call physics_ptend_init(ptend_loc, state%psetcols, 'zm_convr', ls=.true., lq=lq, lu=.true., lv=.true.) 
!>songxl 2011-09-20-------------------- 
!
! Associate pointers with physics buffer fields
!
   itim = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim/), kount=(/pcols,pver,1/) )

   call pbuf_get_field(pbuf, icwmrdp_idx,     ql )
   call pbuf_get_field(pbuf, rprddp_idx,      rprd )
   call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/),    kount=(/pcols, pver, pcnst/) )
   call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
   call pbuf_get_field(pbuf, prec_dp_idx,     prec )
   call pbuf_get_field(pbuf, snow_dp_idx,     snow )
!<songxl 2011-09-20---------------------------
   zmpu(:ncol,:pver)  = 0._r8
   zmpv(:ncol,:pver)  = 0._r8
!>songxl 2011-09-20---------------------------
!<songxl 2014-11-20------------------
   call pbuf_get_field(pbuf, icimrdp_idx,     qi )
   call pbuf_get_field(pbuf, dlfzm_idx,  dlf)
   call pbuf_get_field(pbuf, difzm_idx,  dif)

   if (zm_microp) then
      call pbuf_get_field(pbuf, dnlfzm_idx, dnlf)
      call pbuf_get_field(pbuf, dnifzm_idx, dnif)
   else
      allocate(dnlf(pcols,pver), dnif(pcols,pver))
   end if

   call pbuf_get_field(pbuf, lambdadpcu_idx, lambdadpcu)
   call pbuf_get_field(pbuf, mudpcu_idx,     mudpcu)

   if (zm_microp) then

      if (nmodes > 0) then

         ! Associate pointers with the modes and species that affect the climate (list 0)

         do m = 1, nmodes
            call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, aero(lchnk)%num_a(m)%val)
            call pbuf_get_field(pbuf, dgnum_idx, aero(lchnk)%dgnum(m)%val, start=(/1,1,m/), kount=(/pcols,pver,1/))

            do l = 1, aero(lchnk)%nspec(m)
               call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, aero(lchnk)%mmr_a(l,m)%val)
            end do
         end do

      else if (nbulk > 0) then

         ! Associate pointers with the bulk aerosols that affect the climate (list 0)

         do m = 1, nbulk
            call rad_cnst_get_aer_mmr(0, m, state, pbuf, aero(lchnk)%mmr_bulk(m)%val)
         end do

      end if
   end if


!>songxl 2014-11-20------------------
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')

   !++wy
   !if (zm_stc) then
   !    call pbuf_get_field(pbuf, t_temp_idx,  t_temp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, q_temp_idx,  q_temp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, u_temp_idx,  u_temp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, v_temp_idx,  v_temp, start=(/1,1,1/), kount=(/pcols,pver,6/) )

   !    call pbuf_get_field(pbuf, t_stemp_idx,  t_stemp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, q_stemp_idx,  q_stemp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, u_stemp_idx,  u_stemp, start=(/1,1,1/), kount=(/pcols,pver,6/) )
   !    call pbuf_get_field(pbuf, v_stemp_idx,  v_stemp, start=(/1,1,1/), kount=(/pcols,pver,6/) )

   !    call avg_tquv(lchnk, ncol, state%lat, state%lon, a_flag,   &
   !                  state%t, state%q(:,:,1), state%u, state%v,  &
   !                  t_temp, q_temp, u_temp, v_temp,  &
   !                  t_stemp, q_stemp, u_stemp, v_stemp,  &
   !                  t_spatial_avg(:,lchnk,:), q_spatial_avg(:,lchnk,:), &
   !                  u_spatial_avg(:,lchnk,:), v_spatial_avg(:,lchnk,:), &
   !                  t_avgd, q_avgd, u_avgd, v_avgd)

   !    do timesteps = 2, 6
   !        t_temp(:,:,timesteps-1) = t_temp(:,:,timesteps)
   !        q_temp(:,:,timesteps-1) = q_temp(:,:,timesteps)
   !        u_temp(:,:,timesteps-1) = u_temp(:,:,timesteps)
   !        v_temp(:,:,timesteps-1) = v_temp(:,:,timesteps)
   !    end do
   !    t_temp(:,:,6) = state%t(:,:)
   !    q_temp(:,:,6) = state%q(:,:,1)
   !    u_temp(:,:,6) = state%u(:,:)
   !    v_temp(:,:,6) = state%v(:,:)

   !    do timesteps = 2, 6
   !        t_stemp(:,:,timesteps-1) = t_stemp(:,:,timesteps)
   !        q_stemp(:,:,timesteps-1) = q_stemp(:,:,timesteps)
   !        u_stemp(:,:,timesteps-1) = u_stemp(:,:,timesteps)
   !        v_stemp(:,:,timesteps-1) = v_stemp(:,:,timesteps)
   !    end do
   !    t_stemp(:,:,6) = t_spatial_avg(:,lchnk,:)
   !    q_stemp(:,:,6) = q_spatial_avg(:,lchnk,:)
   !    u_stemp(:,:,6) = u_spatial_avg(:,lchnk,:)
   !    v_stemp(:,:,6) = v_spatial_avg(:,lchnk,:)

   !    call get_area_all_p(lchnk, ncol, area)

   !    do i = 1, pcols
   !        do k = 1, pver
   !            if (t_avgd(i,k) .eq. 0._r8) then
   !                t_avgd(i,k) = state%t(i,k)
   !            end if

   !            if (q_avgd(i,k) .eq. 0._r8) then
   !                q_avgd(i,k) = state%q(i,k,1)
   !            end if
   !        end do
   !    end do
   !end if
   !--wy

   call zm_convr(   lchnk   ,ncol    , &
                    state%t       ,state%q(:,:,1)     ,prec    ,jctop   ,jcbot   , &
                    pblh    ,state%zm      ,state%phis    ,state%zi      ,ptend_loc%q(:,:,1)    , &
                    ptend_loc%s    ,state%pmid     ,state%pint    ,state%pdel     , &
                    .5_r8*ztodt    ,mcon    ,cme     , cape,      &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu(:,:,lchnk),md(:,:,lchnk),du(:,:,lchnk),eu(:,:,lchnk),ed(:,:,lchnk)      , &
                    dp(:,:,lchnk) ,dsubcld(:,lchnk) ,jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk)   , &
!<songxl 2014-11-20-----------
!                   lengath(lchnk) ,ql      ,rliq  ,landfrac   ) 
                    lengath(lchnk) ,ql      ,rliq  ,landfrac, qi, qliq, &
                    qice, dif, dnlf, dnif, wu, &
                    sprd, qrain, qsnow, qnl, qni, &
                    qnr, qns, frz, aero(lchnk),  &
                    autolm, accrlm, bergnm, fhtimm, fhtctm, fhmlm, hmpim, accslm, dlfm,       &
                    autoln, accrln, bergnn, fhtimn, fhtctn, fhmln, accsln, activn, dlfn,      &
                    autoim, accsim, difm, nuclin, autoin, accsin, hmpin, difn, cmel, cmei,    &
                    trspcm, trspcn, trspim, trspin, lambdadpcu, mudpcu, dcape, rice ,         &
                    !++wy
                    !t_avgd, q_avgd, u_avgd, v_avgd, area, &
                    state%t, state%q(:,:,1), state%u, state%v, area, &
                    cape_avgd,                            &
                    !songxl 2011-09-20
                    state%u, state%v  ,ptend_loc%u,  ptend_loc%v, zmpu,  zmpv)

   if (.not. zm_microp) then
      deallocate(dnlf, dnif)
   end if
!>songxl 2014-11-20-------------

   call outfld('CAPE', cape, pcols, lchnk)        ! RBN - CAPE output
   !++wy
   if (zm_stc) then
       call outfld('CAPE_AVG', cape_avgd, pcols, lchnk)        ! RBN - CAPE output
   end if
   !--wy
!
! Output fractional occurance of ZM convection
!
   freqzm(:) = 0._r8
   do i = 1,lengath(lchnk)
      freqzm(ideep(i,lchnk)) = 1.0_r8
   end do
   call outfld('FREQZM  ',freqzm          ,pcols   ,lchnk   )
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit
!<songxl 2014-11-20-------
   call outfld('CMFMCDZM', mcon, pcols, lchnk)
!>songxl 2014-11-20-------

   ! Store upward and downward mass fluxes in un-gathered arrays
   ! + convert from mb/s to kg/m^2/s
   do i=1,lengath(lchnk) 
      do k=1,pver
         ii = ideep(i,lchnk)
         mu_out(ii,k) = mu(i,k,lchnk) * 100._r8/gravit
         md_out(ii,k) = md(i,k,lchnk) * 100._r8/gravit
      end do
   end do

   call outfld('ZMMU', mu_out(1,1), pcols, lchnk)
   call outfld('ZMMD', md_out(1,1), pcols, lchnk)

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
!<songxl 2011-09-20----------------------------
   call outfld('ZMDU    ',ptend_loc%u    ,pcols   ,lchnk   )
   call outfld('ZMDV    ',ptend_loc%v    ,pcols   ,lchnk   )
   call outfld('ZMPU    ',zmpu           ,pcols   ,lchnk   )
   call outfld('ZMPV    ',zmpv           ,pcols   ,lchnk   )
!>songxl 2011-09-20---------------------------
   call t_stopf ('zm_convr')

!<songxl 2014-11-20--------
    
   if (zm_microp) then !++wy fix
       call outfld('DCAPE',dcape, pcols, lchnk)
   end if                   
   do k = 1,pver
      do i = 1,ncol
         if (qice(i,k) .gt. 0.0_r8) then
            cice_snum(i,k) = 1.0_r8
         else
            cice_snum(i,k) = 0.0_r8
         end if
         if (qliq(i,k) .gt. 0.0_r8) then
            cliq_snum(i,k) = 1.0_r8
         else
            cliq_snum(i,k) = 0.0_r8
         end if
         if (qsnow(i,k) .gt. 0.0_r8) then
            csnow_snum(i,k) = 1.0_r8
         else
            csnow_snum(i,k) = 0.0_r8
         end if
         if (qrain(i,k) .gt. 0.0_r8) then
            crain_snum(i,k) = 1.0_r8
         else
            crain_snum(i,k) = 0.0_r8
         end if

         if (qnl(i,k) .gt. 0.0_r8) then
            qnl_snum(i,k) = 1.0_r8
         else
            qnl_snum(i,k) = 0.0_r8
         end if
         if (qni(i,k) .gt. 0.0_r8) then
            qni_snum(i,k) = 1.0_r8
         else
            qni_snum(i,k) = 0.0_r8
         end if
         if (wu(i,k) .gt. 0.0_r8) then
            wu_snum(i,k) = 1.0_r8
         else
            wu_snum(i,k) = 0.0_r8
         end if

         if (zdu(i,k) .gt. 0.0_r8) then
            du_snum(i,k) = 1.0_r8
         else
            du_snum(i,k) = 0.0_r8
         end if


         if (autoln(i,k) .ne. 0.0_r8) then
           autoln_snum(i,k) = 1.0_r8
         else
           autoln_snum(i,k) = 0.0_r8
         end if

         if (accrln(i,k) .ne. 0.0_r8) then
           accrln_snum(i,k) = 1.0_r8
         else
           accrln_snum(i,k) = 0.0_r8
         end if

         if (bergnn(i,k) .ne. 0.0_r8) then
           bergnn_snum(i,k) = 1.0_r8
         else
           bergnn_snum(i,k) = 0.0_r8
         end if

         if (fhtimn(i,k) .ne. 0.0_r8) then
           fhtimn_snum(i,k) = 1.0_r8
         else
           fhtimn_snum(i,k) = 0.0_r8
         end if

         if (fhtctn(i,k) .ne. 0.0_r8) then
           fhtctn_snum(i,k) = 1.0_r8
         else
           fhtctn_snum(i,k) = 0.0_r8
         end if

         if (fhmln(i,k) .ne. 0.0_r8) then
           fhmln_snum(i,k) = 1.0_r8
         else
           fhmln_snum(i,k) = 0.0_r8
         end if

         if (accsln(i,k) .ne. 0.0_r8) then
           accsln_snum(i,k) = 1.0_r8
         else
           accsln_snum(i,k) = 0.0_r8
         end if

         if (activn(i,k) .ne. 0.0_r8) then
           activn_snum(i,k) = 1.0_r8
         else
           activn_snum(i,k) = 0.0_r8
         end if

         if (dlfn(i,k) .ne. 0.0_r8) then
           dlfn_snum(i,k) = 1.0_r8
         else
           dlfn_snum(i,k) = 0.0_r8
         end if

         if (nuclin(i,k) .ne. 0.0_r8) then
           nuclin_snum(i,k) = 1.0_r8
         else
           nuclin_snum(i,k) = 0.0_r8
         end if

         if (autoin(i,k) .ne. 0.0_r8) then
           autoin_snum(i,k) = 1.0_r8
         else
           autoin_snum(i,k) = 0.0_r8
         end if
         if (accsin(i,k) .ne. 0.0_r8) then
           accsin_snum(i,k) = 1.0_r8
         else
           accsin_snum(i,k) = 0.0_r8
         end if

         if (hmpin(i,k) .ne. 0.0_r8) then
           hmpin_snum(i,k) = 1.0_r8
         else
           hmpin_snum(i,k) = 0.0_r8
         end if

         if (difn(i,k) .ne. 0.0_r8) then
           difn_snum(i,k) = 1.0_r8
         else
           difn_snum(i,k) = 0.0_r8
         end if

         if (cmel(i,k) .ne. 0.0_r8) then
           cmel_snum(i,k) = 1.0_r8
         else
           cmel_snum(i,k) = 0.0_r8
         end if

         if (cmei(i,k) .ne. 0.0_r8) then
           cmei_snum(i,k) = 1.0_r8
         else
           cmei_snum(i,k) = 0.0_r8
         end if

         if (trspcn(i,k) .ne. 0.0_r8) then
           trspcn_snum(i,k) = 1.0_r8
         else
           trspcn_snum(i,k) = 0.0_r8
         end if

         if (trspin(i,k) .ne. 0.0_r8) then
           trspin_snum(i,k) = 1.0_r8
         else
           trspin_snum(i,k) = 0.0_r8
         end if

         if (frz(i,k) .ne. 0.0_r8) then
           frz_snum(i,k) = 1.0_r8
         else
           frz_snum(i,k) = 0.0_r8
         end if
      end do
   end do
   
   if (zm_microp) then
      call outfld('CLDLIQZM',qliq           ,pcols, lchnk)
      call outfld('CLDICEZM',qice           ,pcols, lchnk)
      call outfld('CLIQSNUM',cliq_snum      ,pcols, lchnk)
      call outfld('CICESNUM',cice_snum      ,pcols, lchnk)
      call outfld('QRAINZM' ,qrain          ,pcols, lchnk)
      call outfld('QSNOWZM' ,qsnow          ,pcols, lchnk)
      call outfld('CRAINNUM',crain_snum     ,pcols, lchnk)
      call outfld('CSNOWNUM',csnow_snum     ,pcols, lchnk)

      call outfld('DIFZM'   ,dif            ,pcols, lchnk)
      call outfld('DLFZM'   ,dlf            ,pcols, lchnk)
      call outfld('DNIFZM'  ,dnif           ,pcols, lchnk)
      call outfld('DNLFZM'  ,dnlf           ,pcols, lchnk)
      call outfld('WUZM'    ,wu             ,pcols, lchnk)
      call outfld('WUZMSNUM',wu_snum        ,pcols, lchnk)
      call outfld('QNLZM'   ,qnl            ,pcols, lchnk)
      call outfld('QNIZM'   ,qni            ,pcols, lchnk)
      call outfld('QNRZM'   ,qnr            ,pcols, lchnk)
      call outfld('QNSZM'   ,qns            ,pcols, lchnk)
      call outfld('FRZZM'  ,frz            ,pcols, lchnk)

      call outfld('AUTOL_M' ,autolm         ,pcols, lchnk)
      call outfld('ACCRL_M' ,accrlm         ,pcols, lchnk)
      call outfld('BERGN_M' ,bergnm         ,pcols, lchnk)
      call outfld('FHTIM_M' ,fhtimm         ,pcols, lchnk)
      call outfld('FHTCT_M' ,fhtctm         ,pcols, lchnk)
      call outfld('FHML_M'  ,fhmlm          ,pcols, lchnk)
      call outfld('HMPI_M'  ,hmpim          ,pcols, lchnk)
      call outfld('ACCSL_M' ,accslm         ,pcols, lchnk)
      call outfld('DLF_M'   ,dlfm           ,pcols, lchnk)

      call outfld('AUTOL_N' ,autoln         ,pcols, lchnk)
      call outfld('ACCRL_N' ,accrln         ,pcols, lchnk)
      call outfld('BERGN_N' ,bergnn         ,pcols, lchnk)
      call outfld('FHTIM_N' ,fhtimn         ,pcols, lchnk)
      call outfld('FHTCT_N' ,fhtctn         ,pcols, lchnk)
      call outfld('FHML_N'  ,fhmln          ,pcols, lchnk)
      call outfld('ACCSL_N' ,accsln         ,pcols, lchnk)
      call outfld('ACTIV_N' ,activn         ,pcols, lchnk)
      call outfld('DLF_N'   ,dlfn           ,pcols, lchnk)
      call outfld('AUTOI_M' ,autoim         ,pcols, lchnk)
      call outfld('ACCSI_M' ,accsim         ,pcols, lchnk)
      call outfld('DIF_M'   ,difm           ,pcols, lchnk)
      call outfld('NUCLI_N' ,nuclin         ,pcols, lchnk)
      call outfld('AUTOI_N' ,autoin         ,pcols, lchnk)
      call outfld('ACCSI_N' ,accsin         ,pcols, lchnk)
      call outfld('HMPI_N'  ,hmpin          ,pcols, lchnk)
      call outfld('DIF_N'   ,difn           ,pcols, lchnk)
      call outfld('COND_M'  ,cmel           ,pcols, lchnk)
      call outfld('DEPOS_M' ,cmei           ,pcols, lchnk)

      call outfld('TRSPC_M' ,trspcm         ,pcols, lchnk)
      call outfld('TRSPC_N' ,trspcn         ,pcols, lchnk)
      call outfld('TRSPI_M' ,trspim         ,pcols, lchnk)
      call outfld('TRSPI_N' ,trspin         ,pcols, lchnk)

      call outfld('DUSNUM'    ,du_snum       ,pcols   ,lchnk   )

      call outfld('AUTOL_SN'  ,autoln_snum   ,pcols   ,lchnk   )
      call outfld('ACCRL_SN'  ,accrln_snum   ,pcols   ,lchnk   )
      call outfld('BERGN_SN'  ,bergnn_snum   ,pcols   ,lchnk   )
      call outfld('FHTIM_SN'  ,fhtimn_snum   ,pcols   ,lchnk   )
      call outfld('FHTCT_SN'  ,fhtctn_snum   ,pcols   ,lchnk   )
      call outfld('FHML_SN'   ,fhmln_snum    ,pcols   ,lchnk   )
      call outfld('ACCSL_SN'  ,accsln_snum   ,pcols   ,lchnk   )
      call outfld('ACTIV_SN'  ,activn_snum   ,pcols   ,lchnk   )
      call outfld('DLF_SN'    ,dlfn_snum     ,pcols   ,lchnk   )
      call outfld('NUCLI_SN'  ,nuclin_snum   ,pcols   ,lchnk   )
      call outfld('AUTOI_SN'  ,autoin_snum   ,pcols   ,lchnk   )
      call outfld('ACCSI_SN'  ,accsin_snum   ,pcols   ,lchnk   )
      call outfld('HMPI_SN'   ,hmpin_snum    ,pcols   ,lchnk   )
      call outfld('DIF_SN'    ,difn_snum     ,pcols   ,lchnk   )
      call outfld('COND_SN'   ,cmel_snum     ,pcols   ,lchnk   )
      call outfld('DEPOS_SN'  ,cmei_snum     ,pcols   ,lchnk   )
      call outfld('TRSPC_SN'  ,trspcn_snum   ,pcols   ,lchnk   )
      call outfld('TRSPI_SN'  ,trspin_snum   ,pcols   ,lchnk   )
      call outfld('FRZZM_SN' ,frz_snum      ,pcols   ,lchnk   )
      call outfld('ZMDT_SN'   ,zmdt_snum     ,pcols   ,lchnk   )
   end if

!>songxl 2014-11-20--------


!    do i = 1,pcols
!    do i = 1,nco
   pcont(:ncol) = state%ps(:ncol)
   pconb(:ncol) = state%ps(:ncol)
   do i = 1,lengath(lchnk)
       if (maxg(i,lchnk).gt.jt(i,lchnk)) then
          pcont(ideep(i,lchnk)) = state%pmid(ideep(i,lchnk),jt(i,lchnk))  ! gathered array (or jctop ungathered)
          pconb(ideep(i,lchnk)) = state%pmid(ideep(i,lchnk),maxg(i,lchnk))! gathered array
       endif
       !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
    call outfld('PCONVT  ',pcont          ,pcols   ,lchnk   )
    call outfld('PCONVB  ',pconb          ,pcols   ,lchnk   )

    if (zm_microp) then
       do i = 1,ncol
          if ( pcont(i).eq.state%ps(i)) then
             pcont_snum(i) =0.0_r8
             pszm(i) = state%ps(i)
          else
             pcont_snum(i)= 1.0_r8
             pszm(i) = 0.0_r8
          end if
       end do

       call outfld('PCT_SN'   ,pcont_snum     ,pcols   ,lchnk   )
       call outfld('PSZM'     ,pszm           ,pcols   ,lchnk   )
    end if 

  ! This name triggers a special case in physics_types.F90:physics_update()
  call physics_ptend_init(ptend_all, state%psetcols, 'convect_deep')

  ! add tendency from this process to tendencies from other processes
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

  ! initialize ptend for next process
  lq(:) = .FALSE.
  lq(1) = .TRUE.
  call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

   call t_startf ('zm_conv_evap')
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
!

    call pbuf_get_field(pbuf, dp_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, dp_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, dp_cldliq_idx, dp_cldliq  )
    call pbuf_get_field(pbuf, dp_cldice_idx, dp_cldice  )
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8

    call zm_conv_evap(state1%ncol,state1%lchnk, &
         state1%t,state1%pmid,state1%pdel,state1%q(:pcols,:pver,1), &
         ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt, ptend_loc%q(:pcols,:pver,1), &
         rprd, cld, ztodt, &
!<songxl 2014-11-20---------
!         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow)
         prec, snow, ntprprd, ntsnprd , flxprec, flxsnow, sprd, old_snow)
!>songxl 2014-11-20----------
 
    evapcdp(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)
!
! Write out variables from zm_conv_evap
!
   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   call outfld('EVAPTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwprd  (:ncol,:pver)/cpair
   call outfld('FZSNTZM ',ftem           ,pcols   ,lchnk   )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver)/cpair
   call outfld('EVSNTZM ',ftem           ,pcols   ,lchnk   )
   call outfld('EVAPQZM ',ptend_loc%q(1,1,1) ,pcols   ,lchnk   )
   call outfld('ZMFLXPRC', flxprec, pcols, lchnk)
   call outfld('ZMFLXSNW', flxsnow, pcols, lchnk)
   call outfld('ZMNTPRPD', ntprprd, pcols, lchnk)
   call outfld('ZMNTSNPD', ntsnprd, pcols, lchnk)
   call outfld('ZMEIHEAT', ptend_loc%s, pcols, lchnk)
!songxl 2014-11-20   call outfld('CMFMCDZM   ',mcon ,  pcols   ,lchnk   )
   call outfld('PRECCDZM   ',prec,  pcols   ,lchnk   )


   call t_stopf ('zm_conv_evap')

   call outfld('PRECZ   ', prec   , pcols, lchnk)

   if (zm_microp) then
      do i = 1,ncol
         if (prec(i) .gt. 0.0_r8) then
            precz_snum(i) = 1.0_r8
         else
            precz_snum(i) = 0.0_r8
         end if
      end do
      call outfld('PRECZ_SN', precz_snum, pcols, lchnk)
   end if 

  ! add tendency from this process to tend from other processes here
  call physics_ptend_sum(ptend_loc,ptend_all, ncol)

  ! update physics state type state1 with ptend_loc 
  call physics_update(state1, ptend_loc, ztodt)

!<songxl 2011-09-20---------------------------------
!#if 0
#if 1
  !write(6,*) "NCAR CMT scheme"

  ! Momentum Transport (non-cam3 physics)

  if ( .not. cam_physpkg_is('cam3')) then

     call physics_ptend_init(ptend_loc, state1%psetcols, 'momtran', ls=.true., lu=.true., lv=.true.)

     winds(:ncol,:pver,1) = state1%u(:ncol,:pver)
     winds(:ncol,:pver,2) = state1%v(:ncol,:pver)
   
     l_windt(1) = .true.
     l_windt(2) = .true.

     call t_startf ('momtran')
     call momtran (lchnk, ncol,                                        &
                   l_windt,winds, 2,  mu(1,1,lchnk), md(1,1,lchnk),   &
                   du(1,1,lchnk), eu(1,1,lchnk), ed(1,1,lchnk), dp(1,1,lchnk), dsubcld(1,lchnk),  &
                   jt(1,lchnk),maxg(1,lchnk), ideep(1,lchnk), 1, lengath(lchnk),  &
                   nstep,  wind_tends, pguall, pgdall, icwu, icwd, ztodt, seten )  
     call t_stopf ('momtran')

     ptend_loc%u(:ncol,:pver) = wind_tends(:ncol,:pver,1)
     ptend_loc%v(:ncol,:pver) = wind_tends(:ncol,:pver,2)
     ptend_loc%s(:ncol,:pver) = seten(:ncol,:pver)  

     call physics_ptend_sum(ptend_loc,ptend_all, ncol)

     ! update physics state type state1 with ptend_loc 
     call physics_update(state1, ptend_loc, ztodt)

     ftem(:ncol,:pver) = seten(:ncol,:pver)/cpair
     call outfld('ZMMTT', ftem             , pcols, lchnk)
     call outfld('ZMMTU', wind_tends(1,1,1), pcols, lchnk)
     call outfld('ZMMTV', wind_tends(1,1,2), pcols, lchnk)
   
     ! Output apparent force from  pressure gradient
     call outfld('ZMUPGU', pguall(1,1,1), pcols, lchnk)
     call outfld('ZMUPGD', pgdall(1,1,1), pcols, lchnk)
     call outfld('ZMVPGU', pguall(1,1,2), pcols, lchnk)
     call outfld('ZMVPGD', pgdall(1,1,2), pcols, lchnk)

     ! Output in-cloud winds
     call outfld('ZMICUU', icwu(1,1,1), pcols, lchnk)
     call outfld('ZMICUD', icwd(1,1,1), pcols, lchnk)
     call outfld('ZMICVU', icwu(1,1,2), pcols, lchnk)
     call outfld('ZMICVD', icwd(1,1,2), pcols, lchnk)

   end if
#endif
!>songxl 2011-09-20-------------------

   ! Transport cloud water and ice only
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   lq(:)  = .FALSE.
   lq(2:) = cnst_is_convtran1(2:)
   call physics_ptend_init(ptend_loc, state1%psetcols, 'convtran1', lq=lq)


   ! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
   ! ratios are moist
   fake_dpdry(:,:) = 0._r8

   call t_startf ('convtran1')
   call convtran (lchnk,                                        &
                  ptend_loc%lq,state1%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                  du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dp(:,:,lchnk), dsubcld(:,lchnk),  &
                  jt(:,lchnk),maxg(:,lchnk), ideep(:,lchnk), 1, lengath(lchnk),  &
                  nstep,   fracis,  ptend_loc%q, fake_dpdry, ztodt)
   call t_stopf ('convtran1')

   call outfld('ZMDICE ',ptend_loc%q(1,1,ixcldice) ,pcols   ,lchnk   )
   call outfld('ZMDLIQ ',ptend_loc%q(1,1,ixcldliq) ,pcols   ,lchnk   )

   ! add tendency from this process to tend from other processes here
   call physics_ptend_sum(ptend_loc,ptend_all, ncol)

   call physics_state_dealloc(state1)
   call physics_ptend_dealloc(ptend_loc)

end subroutine zm_conv_tend
!=========================================================================================


subroutine zm_conv_tend_2( state,  ptend,  ztodt, pbuf)

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   use time_manager,  only: get_nstep
   use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   use constituents,  only: pcnst, cnst_get_ind, cnst_is_convtran1
   use error_messages, only: alloc_err	
 
! Arguments
   type(physics_state), intent(in )   :: state          ! Physics state variables
   type(physics_ptend), intent(out)   :: ptend          ! indivdual parameterization tendencies
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

! Local variables
   integer :: i, lchnk, istat
   integer :: nstep
   real(r8), dimension(pcols,pver) :: dpdry

! physics buffer fields 
   integer itim, ifld
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble
   logical   :: lq(pcnst)

!
! Initialize
!
  lq(:) = .FALSE.
  lq(:) = .not. cnst_is_convtran1(:)
  call physics_ptend_init(ptend, state%psetcols, 'convtran2', lq=lq )

!
! Associate pointers with physics buffer fields
!
   ifld = pbuf_get_index('FRACIS')
   call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

!
! Transport all constituents except cloud water and ice
!

  lchnk = state%lchnk

   nstep = get_nstep()

   if (any(ptend%lq(:))) then
      ! initialize dpdry for call to convtran
      ! it is used for tracers of dry mixing ratio type
      dpdry = 0._r8
      do i = 1,lengath(lchnk)
         dpdry(i,:) = state%pdeldry(ideep(i,lchnk),:)/100._r8
      end do

      call t_startf ('convtran2')
      call convtran (lchnk,                                        &
                     ptend%lq,state%q, pcnst,  mu(:,:,lchnk), md(:,:,lchnk),   &
                     du(:,:,lchnk), eu(:,:,lchnk), ed(:,:,lchnk), dp(:,:,lchnk), dsubcld(:,lchnk),  &
                     jt(:,lchnk),maxg(:,lchnk),ideep(:,lchnk), 1, lengath(lchnk),  &
                     nstep,   fracis,  ptend%q, dpdry, ztodt)
      call t_stopf ('convtran2')
   end if

end subroutine zm_conv_tend_2

!=========================================================================================

!++wy
subroutine avg_tquv(lchnk , ncol  , lat_c, lon_c, a_flag ,          &
                    t     , q     , u     , v     ,                 &
                    t_temp, q_temp, u_temp, v_temp,                 &
                    t_stemp, q_stemp, u_stemp, v_stemp,             &
                    t_savg, q_savg, u_savg, v_savg,                 &
                    t_avgd, q_avgd, u_avgd, v_avgd)

    use ppgrid,           only: pcols, pver
    use time_manager,     only: get_nstep
    use abortutils,       only: endrun
    use phys_grid,        only: get_lon_p

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    real(r8), intent(in) :: lat_c(pcols)
    real(r8), intent(in) :: lon_c(pcols)
    !average flag: 0, without performing average
    !average flag: 1, performing spatial average
    !average flag: 2, performing temporal average 
    !average flag: 3, performing spatial and temporal average 
    integer, intent(in) :: a_flag                  ! average flag

    real(r8), intent(in) :: t_temp(pcols,pver,6)
    real(r8), intent(in) :: q_temp(pcols,pver,6)
    real(r8), intent(in) :: u_temp(pcols,pver,6)
    real(r8), intent(in) :: v_temp(pcols,pver,6)

    real(r8), intent(in) :: t_stemp(pcols,pver,6)
    real(r8), intent(in) :: q_stemp(pcols,pver,6)
    real(r8), intent(in) :: u_stemp(pcols,pver,6)
    real(r8), intent(in) :: v_stemp(pcols,pver,6)
    real(r8), intent(in) :: t_savg(pcols,pver)
    real(r8), intent(in) :: q_savg(pcols,pver)
    real(r8), intent(in) :: u_savg(pcols,pver)
    real(r8), intent(in) :: v_savg(pcols,pver)

    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: q(pcols,pver)
    real(r8), intent(in) :: u(pcols,pver)
    real(r8), intent(in) :: v(pcols,pver)

    real(r8), intent(out) :: t_avgd(pcols,pver)
    real(r8), intent(out) :: q_avgd(pcols,pver)
    real(r8), intent(out) :: u_avgd(pcols,pver)
    real(r8), intent(out) :: v_avgd(pcols,pver)

    integer i, k
    integer kinv
    integer :: nstep
    integer :: tt


    if (a_flag .ne. 0 .and. a_flag .ne. 1 .and. &
        a_flag .ne. 2 .and. a_flag .ne. 3) then
        print *, "a_flag: ", a_flag
        call endrun('ERROR: average flag should be 0, 1, 2, or 3')
    end if

    if (a_flag == 0) then
        t_avgd = t
        q_avgd = q
        u_avgd = u
        v_avgd = v
    end if

    if (a_flag == 1) then
        ! Perform spatial averaging
        !**************************************************************
        !NOTE: When performing average, please go to phys_grid.F90
        !to set the option def_lbal_opt to 4 or 5 
        !**************************************************************
        t_avgd = t_savg
        q_avgd = q_savg
        u_avgd = u_savg
        v_avgd = v_savg
    end if

    if (a_flag == 2) then
        nstep = get_nstep()

        t_avgd = t
        q_avgd = q
        u_avgd = u
        v_avgd = v
        if (nstep >= 6) then
            do tt = 2, 6
                t_avgd(:,:) = t_avgd(:,:)+t_temp(:,:,tt)
                q_avgd(:,:) = q_avgd(:,:)+q_temp(:,:,tt)
                u_avgd(:,:) = u_avgd(:,:)+u_temp(:,:,tt)
                v_avgd(:,:) = v_avgd(:,:)+v_temp(:,:,tt)
            end do
            t_avgd = t_avgd/6._r8
            q_avgd = q_avgd/6._r8
            u_avgd = u_avgd/6._r8
            v_avgd = v_avgd/6._r8
        else
            do tt = 6, 6-nstep+2, -1
                t_avgd(:,:) = t_avgd(:,:)+t_temp(:,:,tt)
                q_avgd(:,:) = q_avgd(:,:)+q_temp(:,:,tt)
                u_avgd(:,:) = u_avgd(:,:)+u_temp(:,:,tt)
                v_avgd(:,:) = v_avgd(:,:)+v_temp(:,:,tt)
            end do
            if (nstep /= 0) then
                t_avgd = t_avgd/nstep
                q_avgd = q_avgd/nstep
                u_avgd = u_avgd/nstep
                v_avgd = v_avgd/nstep
            else
                t_avgd = t
                q_avgd = q
                u_avgd = u
                v_avgd = v
            end if
        end if
    end if

    if (a_flag == 3) then
        nstep = get_nstep()

        t_avgd = t_savg
        q_avgd = q_savg
        u_avgd = u_savg
        v_avgd = v_savg
        if (nstep >= 6) then
            do tt = 2, 6
                t_avgd(:,:) = t_avgd(:,:)+t_stemp(:,:,tt)
                q_avgd(:,:) = q_avgd(:,:)+q_stemp(:,:,tt)
                u_avgd(:,:) = u_avgd(:,:)+u_stemp(:,:,tt)
                v_avgd(:,:) = v_avgd(:,:)+v_stemp(:,:,tt)
            end do
            t_avgd = t_avgd/6._r8
            q_avgd = q_avgd/6._r8
            u_avgd = u_avgd/6._r8
            v_avgd = v_avgd/6._r8
        else
            do tt = 6, 6-nstep+2, -1
                t_avgd(:,:) = t_avgd(:,:)+t_stemp(:,:,tt)
                q_avgd(:,:) = q_avgd(:,:)+q_stemp(:,:,tt)
                u_avgd(:,:) = u_avgd(:,:)+u_stemp(:,:,tt)
                v_avgd(:,:) = v_avgd(:,:)+v_stemp(:,:,tt)
            end do
            if (nstep /= 0) then
                t_avgd = t_avgd/nstep
                q_avgd = q_avgd/nstep
                u_avgd = u_avgd/nstep
                v_avgd = v_avgd/nstep
            else
                t_avgd = t_savg
                q_avgd = q_savg
                u_avgd = u_savg
                v_avgd = v_savg
            end if
        end if
    end if

end subroutine avg_tquv
!--wy
end module zm_conv_intr
