
module zm_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective 
! precip from the ZM scheme
!
! Apr 2006: RBN: Code added to perform a dilute ascent for closure of the CM mass flux
!                based on an entraining plume a la Raymond and Blythe (1992)
!
! Author: Byron Boville, from code in tphysbc
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8, i8 => shr_kind_i8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp
  use cloud_fraction,  only: cldfrc_fice
  use physconst,       only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpwv, cpliq, rh2o
  use abortutils,      only: endrun
  use cam_logfile,     only: iulog
!<songxl 2014-11-20--------
  use zm_microphysics, only: zm_mphyi, zm_mphy, zm_aero_t
!>songxl 2014-11-20--------

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zmconv_readnl            ! read zmconv_nl namelist
  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! evaporation of precip from ZM schemea
  public convtran                 ! convective transport
  public momtran                  ! convective momentum transport
  public zm_microp                ! true if convective microphysics 
  !++wy
  public zm_stc                   ! true if stochastic convection
  !--wy


!
! Private data
!
!<songxl 2014-11-20--------
   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   real(r8) :: zmconv_c0_lnd = unset_r8    
   real(r8) :: zmconv_c0_ocn = unset_r8    
   real(r8) :: zmconv_ke     = unset_r8    
   logical  :: zmconv_microp = .false.
!>songxl 2014-11-20--------
   !++wy
   logical  :: zmconv_stc    = .false.
   !--wy

   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.
   real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   real(r8) tau   ! convective time scale
   real(r8),parameter :: c1 = 6.112_r8
   real(r8),parameter :: c2 = 17.67_r8
   real(r8),parameter :: c3 = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1
      

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 

!<songxl 2014-11-20-----------   
!   logical :: zm_microp   ! switch for ZM microphysics
    logical :: zm_microp = .false.
!>songxl 2014-11-20-----------

!++wy
    logical :: zm_stc    = .false.
!--wy

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair
   
   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8   

contains

!<songxl 2014-11-20-----------
subroutine zmconv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zmconv_readnl'

   namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_microp, zmconv_stc
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zmconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zmconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      c0_lnd = zmconv_c0_lnd
      c0_ocn = zmconv_c0_ocn
      ke = zmconv_ke
      zm_microp = zmconv_microp
      zm_stc    = zmconv_stc

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(ke,                1, mpir8,  0, mpicom)
   call mpibcast(zm_microp,         1, mpilog, 0, mpicom)
   call mpibcast(zm_stc,            1, mpilog, 0, mpicom)
#endif

end subroutine zmconv_readnl
!>songxl 2014-11-20--------------

!<songxl 2014-11-20--------
subroutine zm_convi(limcnv_in, no_deep_pbl_in)
 
!>songxl 2014-11-20--------
   use dycore,       only: dycore_is, get_resolution

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL 

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   if ( present(no_deep_pbl_in) )  then
      no_deep_pbl = no_deep_pbl_in
   else
      no_deep_pbl = .false.
   endif

   ! tau=4800. were used in canadian climate center. however, in echam3 t42, 
   ! convection is too weak, thus adjusted to 2400.

   hgrid = get_resolution()
   tau = 3600._r8

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convi: tau',tau
      write(iulog,*) 'tuning parameters zm_convi: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn 
      write(iulog,*) 'tuning parameters zm_convi: ke',ke
      write(iulog,*) 'tuning parameters zm_convi: no_deep_pbl',no_deep_pbl
      !++wy
      write(iulog,*) 'logical parameters zm_convi: zm_microp',zm_microp
      write(iulog,*) 'logical parameters zm_convi: zm_stc',zm_stc
      !--wy
   endif

   if (masterproc) write(iulog,*)'**** ZM: DILUTE Buoyancy Calculation ****'

end subroutine zm_convi



subroutine zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
!<songxl 2014-11-20---------
!                    lengath ,ql      ,rliq    ,landfrac)
                    lengath ,ql      ,rliq    ,landfrac,qi      , &   
                    qliq    ,qice    ,dif     ,dnlf    ,dnif    , &
                    wu      ,sprd    ,qrain   ,qsnow   ,qnl     , &
                    qni     ,qnr     ,qns     ,frz     ,aero    , &
                    autolm  ,accrlm  ,bergnm  ,fhtimm  ,fhtctm  , &
                    fhmlm   ,hmpim   ,accslm  ,dlfm    ,autoln  , &
                    accrln  ,bergnn  ,fhtimn  ,fhtctn  ,fhmln   , &
                    accsln  ,activn  ,dlfn    ,autoim  ,accsim  , &
                    difm    ,nuclin  ,autoin  ,accsin  ,hmpin   , &
                    difn    ,cmel    ,cmei    ,trspcm  ,trspcn  , &
                    trspim  ,trspin  ,lambdadpcu, mudpcu, dcape , &
                    rice    ,                                     &
                    !++wy
                    t_avgd  ,q_avgd  ,u_avgd  ,v_avgd, area, &
                    cape_avgd,                               &
                    !songxl 2011-09-20
                    u       ,v       ,utg     , &
                    vtg     ,putg    ,pvtg )


!>songxl 2014-11-20----------
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
! 
!-----------------------------------------------------------------------
   use phys_control, only: cam_physpkg_is
   !++wy
   use time_manager, only : is_first_step, get_nstep
   use shr_const_mod,  only: SHR_CONST_PI, SHR_CONST_REARTH
   use time_manager,     only: get_nstep
   use random_number_generator, only: RAN1
   !--wy

!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qh(pcols,pver)   ! grid slice of specific humidity.
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
!<songxl 2011-09-20------------------------------
   real(r8), intent(in) :: u(pcols,pver)           ! grid slice of U-wind at mid-layer.
   real(r8), intent(in) :: v(pcols,pver)           ! grid slice of V-wind at mid-layer.
!>songxl 2011-09-20------------------------------

!++wy
   real(r8), intent(in) ::  t_avgd(pcols,pver)          ! spatial average of T
   real(r8), intent(in) ::  q_avgd(pcols,pver)          ! spatial average of Q
   real(r8), intent(in) ::  u_avgd(pcols,pver)          ! spatial average of U
   real(r8), intent(in) ::  v_avgd(pcols,pver)          ! spatial average of V
   real(r8), intent(in) ::  area(pcols)                 ! column surface area
!--wy

!<songxl 2014-11-20---------
   type(zm_aero_t), intent(inout) :: aero         ! aerosol object. intent(inout) because the
                                                  ! gathered arrays are set here
                                                  ! before passing object
                                                  ! to microphysics

!>songxl 20140-11-20--------

!
! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   !++wy
   real(r8), intent(out) :: cape_avgd(pcols)   ! w  convective available potential energy.
   !--wy
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
!<songxl 2014-11-20---------
   real(r8), intent(out) :: sprd(pcols,pver)       ! snow production rate
   real(r8), intent(out) :: qice(pcols,pver)       ! convective cloud ice.
   real(r8), intent(out) :: qliq(pcols,pver)       ! convective cloud liquid water.
   real(r8), intent(out) :: qni(pcols,pver)        ! convective cloud ice num concen.
   real(r8), intent(out) :: qnl(pcols,pver)        ! convective cloud liquid water num concen.
   real(r8), intent(out) :: qnr(pcols,pver)        ! convective rain water num concen.
   real(r8), intent(out) :: qns(pcols,pver)        ! convective snow num concen. 
   real(r8), intent(out) :: qrain(pcols,pver)      ! rain water mixing ratio (g/m3)
   real(r8), intent(out) :: qsnow(pcols,pver)      ! snow mixing ratio (g/m3)
   real(r8), intent(out) :: dif(pcols,pver)        ! detrained convective cloud ice mixing ratio.
   real(r8), intent(out) :: dnlf(pcols,pver)       ! detrained convective cloud water num concen.
   real(r8), intent(out) :: dnif(pcols,pver)       ! detrained convective cloud ice num concen.
   real(r8), intent(out) :: wu(pcols,pver)         ! vertical velocity
   real(r8), intent(out) :: frz(pcols,pver)        ! freezing heating
   real(r8), intent(out) :: dcape(pcols)           ! CAPE change due to freezing heating

! tendency for output

   real(r8), intent(out) :: autolm(pcols,pver)    !mass tendency due to autoconversion of droplets to rain
   real(r8), intent(out) :: accrlm(pcols,pver)    !mass tendency due to accretion of droplets by rain
   real(r8), intent(out) :: bergnm(pcols,pver)    !mass tendency due to Bergeron process
   real(r8), intent(out) :: fhtimm(pcols,pver)    !mass tendency due to immersion freezing
   real(r8), intent(out) :: fhtctm(pcols,pver)    !mass tendency due to contact freezing
   real(r8), intent(out) :: fhmlm (pcols,pver)    !mass tendency due to homogeneous freezing
   real(r8), intent(out) :: hmpim (pcols,pver)    !mass tendency due to HM process
   real(r8), intent(out) :: accslm(pcols,pver)    !mass tendency due to accretion of droplets by snow
   real(r8), intent(out) :: dlfm  (pcols,pver)    !mass tendency due to detrainment of droplet
   real(r8), intent(out) :: cmel  (pcols,pver)    !mass tendency due to condensation

   real(r8), intent(out) :: autoln(pcols,pver)    !num tendency due to autoconversion of droplets to rain
   real(r8), intent(out) :: accrln(pcols,pver)    !num tendency due to accretion of droplets by rain
   real(r8), intent(out) :: bergnn(pcols,pver)    !num tendency due to Bergeron process
   real(r8), intent(out) :: fhtimn(pcols,pver)    !num tendency due to immersion freezing
   real(r8), intent(out) :: fhtctn(pcols,pver)    !num tendency due to contact freezing
   real(r8), intent(out) :: fhmln (pcols,pver)    !num tendency due to homogeneous freezing
   real(r8), intent(out) :: accsln(pcols,pver)    !num tendency due to accretion of droplets by snow
   real(r8), intent(out) :: activn(pcols,pver)    !num tendency due to droplets activation
   real(r8), intent(out) :: dlfn  (pcols,pver)    !num tendency due to detrainment of droplet

   real(r8), intent(out) :: autoim(pcols,pver)    !mass tendency due to autoconversion of cloud ice to snow
   real(r8), intent(out) :: accsim(pcols,pver)    !mass tendency due to accretion of cloud ice by snow
   real(r8), intent(out) :: difm  (pcols,pver)    !mass tendency due to detrainment of cloud ice
   real(r8), intent(out) :: cmei  (pcols,pver)    !mass tendency due to deposition

   real(r8), intent(out) :: nuclin(pcols,pver)    !num tendency due to ice nucleation
   real(r8), intent(out) :: autoin(pcols,pver)    !num tendency due to autoconversion of cloud ice to snow
   real(r8), intent(out) :: accsin(pcols,pver)    !num tendency due to accretion of cloud ice by snow
   real(r8), intent(out) :: hmpin (pcols,pver)    !num tendency due to HM process
   real(r8), intent(out) :: difn  (pcols,pver)    !num tendency due to detrainment of cloud ice

   real(r8), intent(out) :: trspcm(pcols,pver)    !mass tendency of droplets due to convective transport
   real(r8), intent(out) :: trspcn(pcols,pver)    !num tendency of droplets due to convective transport
   real(r8), intent(out) :: trspim(pcols,pver)    !mass tendency of ice crystal due to convective transport
   real(r8), intent(out) :: trspin(pcols,pver)    !num tendency of ice crystal due to convective transport

   real(r8), intent(out) :: lambdadpcu(pcols,pver) ! slope of cloud liquid size distr
   real(r8), intent(out) :: mudpcu(pcols,pver)     ! width parameter of droplet size distr

!>songxl 2014-11-20---------

! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8), intent(out) :: eu(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8), intent(out) :: ed(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
   real(r8), intent(out) :: rice(pcols) ! reserved ice (not yet in cldce) for energy integrals

!<songxl 2011-09-20-------------------------
   real(r8), intent(out) :: utg(pcols,pver)      !grid slice of u-wind tendency (real).m/s2
   real(r8), intent(out) :: vtg(pcols,pver)      !grid slice of v-wind tendency (real).m/s2
   real(r8), intent(out) :: putg(pcols,pver)
   real(r8), intent(out) :: pvtg(pcols,pver)
!>songxl 2011-09-20-------------------------

   real(r8) zs(pcols)
   real(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) pflxg(pcols,pverp) ! gather precip flux at each level
   real(r8) cug(pcols,pver)    ! gathered condensation rate
!<songxl 2014-11-20---------
   real(r8) frzg(pcols,pver)   ! gathered freezing rate
!>songxl 2014-11-20---------
   real(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) mumax(pcols)
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.
!<songxl 2014-11-20--------
   real(r8) qi(pcols,pver)                    ! wg grid slice of cloud ice.
!>songxl 2014-11-20-------
!
   real(r8) pblt(pcols)           ! i row of pbl top indices.




!
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!

   !++wy
   integer :: nstep
   real(r8) tmp_tp(pcols,pver)
   real(r8) tmp_qstp(pcols,pver)
   real(r8) tmp_tl(pcols)
   integer  tmp_lcl(pcols)
   integer  tmp_lel(pcols)
   integer  tmp_lon(pcols)
   integer  tmp_maxi(pcols)
   integer  tmp_msg
   !--wy

   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
   real(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
   real(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.


   integer lcl(pcols)                  ! w  base level index of deep cumulus convection.
   integer lel(pcols)                  ! w  index of highest theoretical convective plume.
   integer lon(pcols)                  ! w  index of onset level for deep convection.
   integer maxi(pcols)                 ! w  index of level with largest moist static energy.
   integer index(pcols)
   real(r8) precip
!
! gathered work fields:
!
   real(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
   real(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
   real(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
   real(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
   real(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
   real(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
   real(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
   real(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
   real(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
   real(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
   real(r8) cmeg(pcols,pver)

   real(r8) rprdg(pcols,pver)           ! wg gathered rain production rate

   real(r8) capeg(pcols)               ! wg gathered convective available potential energy.
   !++wy
   real(r8) capeg_avgd(pcols)          ! wg gathered convective available potential energy.
   !--wy
   real(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
   real(r8) landfracg(pcols)            ! wg grid slice of landfrac  

!++wy
   real(r8) areag(pcols)                  ! wg grid area
!--wy

   integer lclg(pcols)       ! wg gathered values of lcl.
   integer lelg(pcols)
!
! work fields arising from gathered calculations.
!
   real(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
   real(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.

   real(r8) tg_avgd(pcols, pver)
   real(r8) qg_avgd(pcols, pver)
   real(r8) ug_avgd(pcols, pver)
   real(r8) vg_avgd(pcols, pver)

   real(r8) s_avgd(pcols, pver)
   real(r8) sg_avgd(pcols, pver)
!--wy

!      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
   real(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
   !++wy
   real(r8) qhat_avgd(pcols,pver)      ! wg grid slice of upper interface mixing ratio.
   !--wy
   real(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
   !++wy
   real(r8) qs_avgd(pcols,pver)        ! wg grid slice of saturation mixing ratio.
   real(r8) es_avgd(pcols)
   !--wy
   real(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
   !++wy
   real(r8) shat_avgd(pcols,pver)      ! wg grid slice of upper interface dry static energy.
   !--wy
   real(r8) hmn(pcols,pver)            ! wg moist static energy.
   real(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
   real(r8) qlg(pcols,pver)
   real(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
   real(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.
!<songxl 2011-09-20----------------------------
   real(r8) dpdx(pcols,pver)
   real(r8) dpdy(pcols,pver)
!>songxl 2011-09-20----------------------------
!      real(r8) ud(pcols,pver)
!      real(r8) vd(pcols,pver)

!<songxl 2014-11-20--------

   real(r8) sprdg(pcols,pver)          ! wg gathered snow production rate
   real(r8) wug(pcols,pver)            ! wg vertical velocity
   real(r8) qliqg(pcols,pver)
   real(r8) qiceg(pcols,pver)
   real(r8) qraing(pcols,pver)         ! rain water mixing ratio (kg/kg)
   real(r8) qsnowg(pcols,pver)         ! snow mixing ratio (kg/kg)
   real(r8) dig(pcols,pver)
   real(r8) dnlg(pcols,pver)
   real(r8) dnig(pcols,pver)
   real(r8) qnlg(pcols,pver)
   real(r8) qnig(pcols,pver)
   real(r8) qnrg(pcols,pver)
   real(r8) qnsg(pcols,pver)

   real(r8) qldeg(pcols,pver)     ! cloud water mixing ratio for detrainment (kg/kg)
   real(r8) qideg(pcols,pver)     ! cloud ice mixing ratio for detrainment (kg/kg)
   real(r8) qncdeg(pcols,pver)    ! cloud water number concentration for detrainment (1/kg)
   real(r8) qnideg(pcols,pver)    ! cloud ice number concentration for detrainment (1/kg)

   real(r8) autolmg(pcols,pver)    !mass tendency due to autoconversion of droplets to rain
   real(r8) accrlmg(pcols,pver)    !mass tendency due to accretion of droplets by rain
   real(r8) bergnmg(pcols,pver)    !mass tendency due to Bergeron process 
   real(r8) fhtimmg(pcols,pver)    !mass tendency due to immersion freezing
   real(r8) fhtctmg(pcols,pver)    !mass tendency due to contact freezing
   real(r8) fhmlmg (pcols,pver)    !mass tendency due to homogeneous freezing
   real(r8) hmpimg (pcols,pver)    !mass tendency due to HM process
   real(r8) accslmg(pcols,pver)    !mass tendency due to accretion of droplets by snow
   real(r8) dlfmg  (pcols,pver)    !mass tendency due to detrainment of droplet
   real(r8) cmelg  (pcols,pver)    !mass tendency due to condensation

   real(r8) autolng(pcols,pver)    !num tendency due to autoconversion of droplets to rain
   real(r8) accrlng(pcols,pver)    !num tendency due to accretion of droplets by rain
   real(r8) bergnng(pcols,pver)    !num tendency due to Bergeron process
   real(r8) fhtimng(pcols,pver)    !num tendency due to immersion freezing
   real(r8) fhtctng(pcols,pver)    !num tendency due to contact freezing
   real(r8) fhmlng (pcols,pver)    !num tendency due to homogeneous freezing
   real(r8) accslng(pcols,pver)    !num tendency due to accretion of droplets by snow
   real(r8) activng(pcols,pver)    !num tendency due to droplets activation
   real(r8) dlfng  (pcols,pver)    !num tendency due to detrainment of droplet

   real(r8) autoimg(pcols,pver)    !mass tendency due to autoconversion of cloud ice to snow
   real(r8) accsimg(pcols,pver)    !mass tendency due to accretion of cloud ice by snow
   real(r8) difmg  (pcols,pver)    !mass tendency due to detrainment of cloud ice
   real(r8) cmeig  (pcols,pver)    !mass tendency due to deposition

   real(r8) nucling(pcols,pver)    !num tendency due to ice nucleation
   real(r8) autoing(pcols,pver)    !num tendency due to autoconversion of cloud ice to snow
   real(r8) accsing(pcols,pver)    !num tendency due to accretion of cloud ice by snow
   real(r8) hmping (pcols,pver)    !num tendency due to HM process
   real(r8) difng  (pcols,pver)    !num tendency due to detrainment of cloud ice

   real(r8) trspcmg(pcols,pver)    !mass tendency of droplets due to convective transport
   real(r8) trspcng(pcols,pver)    !num tendency of droplets due to convective transport
   real(r8) trspimg(pcols,pver)    !mass tendency of ice crystal due to convective transport
   real(r8) trsping(pcols,pver)    !num tendency of ice crystal due to convective transport

   real(r8) :: lambdadpcug(pcols,pver) ! slope of cloud liquid size distr
   real(r8) :: mudpcug(pcols,pver)     ! width parameter of droplet size distr

   real(r8) dcapeg(pcols)           ! CAPE change due to freezing heating
!>songxl 2014-11-20--------

   real(r8) mb(pcols)                  ! wg cloud base mass flux.
    
   !++wy
   real(r8) sum_mb(pcols)              ! wg cloud base mass flux.
   !--wy

   integer jlcl(pcols)
   integer j0(pcols)                 ! wg detrainment initiation level index.
   integer jd(pcols)                 ! wg downdraft initiation level index.

   real(r8) delt                     ! length of model time-step in seconds.

   integer i
   integer ii
   integer k, kk
   integer msg                      !  ic number of missing moisture levels at the top of model.
   real(r8) qdifr
   real(r8) sdifr
!<songxl 2011-09-20---------------------------
   logical momentm
!   data momentm/.TRUE./
  data momentm/.FALSE./
!>songxl 2011-09-20---------------------------
!++wy
   real(r8) qdifr_avgd
   real(r8) sdifr_avgd
!--------------------For Stochastic Parameterization---------------------------------
! Namelist parameters
   !real(r8), parameter :: MeanCMF = 1.0E7                                    ! mean mass flux <m> for a single cloud
   real(r8), parameter :: MeanCMF = 8.0E7                                    ! mean mass flux <m> for a single cloud
   integer,  parameter :: nbins  = 100
   real(r8) :: CLMASSRANGE(2)
   integer seed
   integer(i8) :: ran_seed
   real(r8) rannumbers(pcols,nbins)         ! wg random numbers
   real(r8) rannumbersg(pcols,nbins)        ! wg random numbers
!----------------------------------------------------------------------------------
   real(r8), parameter :: rad2deg = 180._r8/SHR_CONST_PI
   real(r8), parameter :: pi = 3.1415926
   real(r8) :: DELTAR, LOOPRADIUS
   real(r8) :: DELTAM, LOOPMASS
   real(r8) :: LOOPMASS1                            ! units in mb/s for a grid area 
   real(r8) :: MASS
   real(r8) :: MASS1
   real(r8) :: REQMF
   real(r8) :: PROB
   real(r8) :: ENSNC
   logical :: launchcl
   integer :: irannum
   integer  :: n                                    ! counter for convective clouds 
!--wy

!<eongxl 2014-11-20-----------
   integer l, m
   real(r8), parameter :: dcon  = 25.e-6_r8
   real(r8), parameter :: mucon = 5.3_r8
   real(r8) negadq
   logical doliq
!esongxl 2014-11-20-----------
!
!++wy
   CLMASSRANGE = (/0.0_r8, 9.0_r8*MeanCMF/)
!--wy
!--------------------------Data statements------------------------------
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! initialize necessary arrays.
! zero out variables not used in cam
!
   qtnd(:,:) = 0._r8
   heat(:,:) = 0._r8
   mcon(:,:) = 0._r8
   rliq(:ncol)   = 0._r8
   rice(:ncol)   = 0._r8
!<songxl 2011-09-20--------------------
   utg(:,:) = 0._r8
   vtg(:,:) = 0._r8
   putg(:,:) = 0._r8
   pvtg(:,:) = 0._r8
!>songxl 2011-09-20--------------------
!
! initialize convective tendencies
!
   prec(:ncol) = 0._r8
   !++wy
   sum_mb(:ncol) = 0._r8
   !--wy
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k)  = 0._r8
         dsdt(i,k)  = 0._r8
         dudt(i,k)  = 0._r8
         dvdt(i,k)  = 0._r8
         dpdx(i,k)  = 0._r8     !songxl 2011-09-20
         dpdy(i,k)  = 0._r8     !songxl 2011-09-20
         pflx(i,k)  = 0._r8
         pflxg(i,k) = 0._r8
         cme(i,k)   = 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         qlg(i,k)   = 0._r8
         dlf(i,k)   = 0._r8
         dlg(i,k)   = 0._r8
!<songxl 2014-11-20-------
         qi(i,k)    = 0._r8
         frz(i,k)   = 0._r8
         sprd(i,k)  = 0._r8
         qliq(i,k)  = 0._r8
         qice(i,k)  = 0._r8
         qnl(i,k)  = 0._r8
         qni(i,k)  = 0._r8
         qnr(i,k)  = 0._r8
         qns(i,k)  = 0._r8
         qrain(i,k) = 0._r8
         qsnow(i,k) = 0._r8
         dif(i,k)   = 0._r8
         dnlf(i,k)  = 0._r8
         dnif(i,k)  = 0._r8
         wu(i,k)    = 0._r8

         qliqg(i,k) = 0._r8
         qiceg(i,k) = 0._r8
         dig(i,k)   = 0._r8
         qraing(i,k)= 0._r8
         qsnowg(i,k)= 0._r8
         dnlg(i,k)  = 0._r8
         dnig(i,k)  = 0._r8
         wug(i,k)   = 0._r8
         qnlg(i,k)  = 0._r8
         qnig(i,k)  = 0._r8
         qnrg(i,k)  = 0._r8
         qnsg(i,k)  = 0._r8
         frzg(i,k)  = 0._r8
         sprdg(i,k) = 0._r8
         qldeg(i,k)  = 0._r8
         qideg(i,k)  = 0._r8
         qncdeg(i,k) = 0._r8
         qnideg(i,k) = 0._r8

         autolm(i,k) = 0._r8
         accrlm(i,k) = 0._r8
         bergnm(i,k) = 0._r8
         fhtimm(i,k) = 0._r8
         fhtctm(i,k) = 0._r8
         fhmlm (i,k) = 0._r8
         hmpim (i,k) = 0._r8
         accslm(i,k) = 0._r8
         dlfm  (i,k) = 0._r8

         autoln(i,k) = 0._r8
         accrln(i,k) = 0._r8
         bergnn(i,k) = 0._r8
         fhtimn(i,k) = 0._r8
         fhtctn(i,k) = 0._r8
         fhmln (i,k) = 0._r8
         accsln(i,k) = 0._r8
         activn(i,k) = 0._r8
         dlfn  (i,k) = 0._r8
         cmel  (i,k) = 0._r8

         autoim(i,k) = 0._r8
         accsim(i,k) = 0._r8
         difm  (i,k) = 0._r8
         cmei  (i,k) = 0._r8

         nuclin(i,k) = 0._r8
         autoin(i,k) = 0._r8
         accsin(i,k) = 0._r8
         hmpin (i,k) = 0._r8
         difn  (i,k) = 0._r8

         trspcm(i,k) = 0._r8
         trspcn(i,k) = 0._r8
         trspim(i,k) = 0._r8
         trspin(i,k) = 0._r8

         autolmg(i,k) = 0._r8
         accrlmg(i,k) = 0._r8
         bergnmg(i,k) = 0._r8
         fhtimmg(i,k) = 0._r8
         fhtctmg(i,k) = 0._r8
         fhmlmg (i,k) = 0._r8
         hmpimg (i,k) = 0._r8
         accslmg(i,k) = 0._r8
         dlfmg  (i,k) = 0._r8

         autolng(i,k) = 0._r8
         accrlng(i,k) = 0._r8
         bergnng(i,k) = 0._r8
         fhtimng(i,k) = 0._r8
         fhtctng(i,k) = 0._r8
         fhmlng (i,k) = 0._r8
         accslng(i,k) = 0._r8
         activng(i,k) = 0._r8
         dlfng  (i,k) = 0._r8
         cmelg  (i,k) = 0._r8

         autoimg(i,k) = 0._r8
         accsimg(i,k) = 0._r8
         difmg  (i,k) = 0._r8
         cmeig  (i,k) = 0._r8

         nucling(i,k) = 0._r8
         autoing(i,k) = 0._r8
         accsing(i,k) = 0._r8
         hmping (i,k) = 0._r8
         difng  (i,k) = 0._r8

         trspcmg(i,k) = 0._r8
         trspcng(i,k) = 0._r8
         trspimg(i,k) = 0._r8
         trsping(i,k) = 0._r8
!>songxl 2014-11-20-------
      end do
   end do


!++wy
   if (zm_stc) then
   ! generate random numbers in each column
     nstep = get_nstep()
     do i = 1,ncol
        seed = i*lchnk+nstep
        seed = mod(seed,2147483646_i8)+1
        ran_seed = -1*seed
        do k = 1, nbins 
            rannumbers(i,k) = RAN1(ran_seed)
        end do
     end do
   end if
!--wy

!<songxl 2014-11-20----------
   lambdadpcu  = (mucon + 1._r8)/dcon
   mudpcu      = mucon
   lambdadpcug = lambdadpcu
   mudpcug     = mudpcu   
!>songxl 2014-11-20----------

   do i = 1,ncol
      pflx(i,pverp) = 0
      pflxg(i,pverp) = 0
   end do
!
   do i = 1,ncol
      pblt(i) = pver
      dsubcld(i) = 0._r8

      jctop(i) = pver
      jcbot(i) = 1

   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01_r8
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01_r8
         pf(i,k) = paph(i,k)*0.01_r8
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
   do k = pver - 1,msg + 1,-1
      do i = 1,ncol
         if (abs(z(i,k)-zs(i)-pblh(i)) < (zf(i,k)-zf(i,k+1))*0.5_r8) pblt(i) = k
      end do
   end do
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qh(i,k)
         s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
         !++wy
         s_avgd(i,k) = t_avgd(i,k) + (grav/cpres)*z(i,k)
         !--wy
         tp(i,k)=0.0_r8
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
      end do
   end do

   do i = 1,ncol
      capeg(i) = 0._r8
      !++wy
      capeg_avgd(i) = 0._r8
      !--wy
      lclg(i) = 1
      lelg(i) = pver
      maxg(i) = 1
      tlg(i) = 400._r8
      dsubcld(i) = 0._r8
      !++wy
      rannumbersg(i,:) = 0._r8
      areag(i) = area(i)
      !--wy
!<songxl 2014-10-10----------
      dcape(i) = 0._r8
      dcapeg(i) = 0._r8
!>songxl 2014-10-10----------
   end do

   if( cam_physpkg_is('cam3')) then

      !  For cam3 physics package, call non-dilute

      !++wy
      tmp_tp   = tp
      tmp_qstp = qstp
      tmp_tl   = tl
      tmp_lcl  = lcl
      tmp_lel  = lel
      tmp_lon  = lon
      tmp_maxi = maxi
      tmp_msg  = msg
      !--wy

      call buoyan(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf       , &
                  tp      ,qstp    ,tl      ,rl      ,cape     , &
                  pblt    ,lcl     ,lel     ,lon     ,maxi     , &
                  rgas    ,grav    ,cpres   ,msg     , &
                  tpert   )
      !++wy calculate the averaged cape with the smoothed fields
      if (zm_stc)  then
          call buoyan(lchnk   ,ncol    , &
                      q_avgd       ,t_avgd       ,p       ,z       ,pf       , &
                      tmp_tp       ,tmp_qstp     ,tmp_tl  ,rl      ,cape_avgd          , &
                      pblt    ,tmp_lcl     ,tmp_lel     ,tmp_lon     ,tmp_maxi     , &
                      rgas    ,grav    ,cpres   ,tmp_msg     , &
                      tpert   )
      end if
      !--wy
   else

      !  Evaluate Tparcel, qs(Tparcel), buoyancy and CAPE, 
      !     lcl, lel, parcel launch level at index maxi()=hmax

      !++wy
      tmp_tp   = tp
      tmp_qstp = qstp
      tmp_tl   = tl
      tmp_lcl  = lcl
      tmp_lel  = lel
      tmp_lon  = lon
      tmp_maxi = maxi
      tmp_msg  = msg
      !--wy

!<songxl 2015-01-10------------------------------
      call buoyan_dilute(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf       , &
                  tp      ,qstp    ,tl      ,rl      ,cape     , &
                  pblt    ,lcl     ,lel     ,lon     ,maxi     , &
                  rgas    ,grav    ,cpres   ,msg     , &
                  tpert   )



!>songxl 2015-01-10--------------

      !++wy calculate the averaged cape with the smoothed fields
      if (zm_stc)  then
          call buoyan_dilute(lchnk   ,ncol    , &
                      q_avgd       ,t_avgd       ,p       ,z       ,pf       , &
                      tmp_tp       ,tmp_qstp    ,tmp_tl      ,rl      ,cape_avgd     , &
                      pblt    ,tmp_lcl     ,tmp_lel     ,tmp_lon     ,tmp_maxi     , &
                      rgas    ,grav    ,cpres   ,tmp_msg     , &
                      tpert   )
      end if
      !--wy
   end if

!
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).
!
   lengath = 0
   do i=1,ncol
!<songxl 2015-01-10----------------
      if (cape(i) > capelmt) then
         lengath = lengath + 1
         index(lengath) = i
      end if
!>songxl 2015-01-10----------------
   end do

   if (lengath.eq.0) return
   do ii=1,lengath
      i=index(ii)
      ideep(ii)=i
   end do
!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,pver
      do i = 1,lengath
         dp(i,k) = 0.01_r8*dpp(ideep(i),k)
         !++wy
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         qg_avgd(i,k) = q_avgd(ideep(i),k)
         tg_avgd(i,k) = t_avgd(ideep(i),k)
         !--wy
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         pg(i,k) = p(ideep(i),k)
         zg(i,k) = z(ideep(i),k)
         sg(i,k) = s(ideep(i),k)
         !++wy
         sg_avgd(i,k) = s_avgd(ideep(i),k)
         !--wy
         tpg(i,k) = tp(ideep(i),k)
         zfg(i,k) = zf(ideep(i),k)
         qstpg(i,k) = qstp(ideep(i),k)
         !<songxl 2011-09-20-----------------------------
!         ug(i,k) = 0._r8
!         vg(i,k) = 0._r8
         ug(i,k) = u(ideep(i),k)
         vg(i,k) = v(ideep(i),k)
         !++wy
         ug_avgd(i,k) = u_avgd(ideep(i),k)
         vg_avgd(i,k) = v_avgd(ideep(i),k)
         !--wy
      end do
   end do
!
!<songxl 2014-11-20--------
   if (zm_microp) then

      if (aero%scheme == 'modal') then

         do m = 1, aero%nmodes

            do k = 1,pver
               do i = 1,lengath
!                   lnum = aero%numidx_a(m)
!                   aero%numg_a(i,k,m) = qh(ideep(i),k,lnum)
                  aero%numg_a(i,k,m) = aero%num_a(m)%val(ideep(i),k)
                  aero%dgnumg(i,k,m) = aero%dgnum(m)%val(ideep(i),k)
               end do
            end do

            do l = 1, aero%nspec(m)
               do k = 1,pver
                  do i = 1,lengath
!                     lmass = aero%mmridx_a(l,m) ! interstitial
!                     aero%mmrg_a(i,k,l,m) = qh(ideep(i),k,lmass)
                     aero%mmrg_a(i,k,l,m) = aero%mmr_a(l,m)%val(ideep(i),k)
                  end do
               end do
            end do

         end do

      else if (aero%scheme == 'bulk') then

         do m = 1, aero%nbulk
            do k = 1,pver
               do i = 1,lengath
                  aero%mmrg_bulk(i,k,m) = aero%mmr_bulk(m)%val(ideep(i),k)
               end do
            end do
         end do

      end if

   end if
!>songxl 2014-11-20---------

   do i = 1,lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
   end do
   do i = 1,lengath
      capeg(i) = cape(ideep(i))
      !++wy
      capeg_avgd(i) = cape_avgd(ideep(i))
      !--wy
      lclg(i) = lcl(ideep(i))
      lelg(i) = lel(ideep(i))
      maxg(i) = maxi(ideep(i))
      tlg(i) = tl(ideep(i))
      landfracg(i) = landfrac(ideep(i))
      !++wy
      areag(i) = area(ideep(i))
      if (zm_stc) then
          rannumbersg(i,:) = rannumbers(ideep(i),:)
      end if
      !--wy
   end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = msg + 1,pver
      do i = 1,lengath
         if (k >= maxg(i)) then
            dsubcld(i) = dsubcld(i) + dp(i,k)
         end if
      end do
   end do
!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
!            alpha(i,k) = 0.5
         sdifr = 0._r8
         qdifr = 0._r8
         !++wy
         sdifr_avgd = 0._r8
         qdifr_avgd = 0._r8
         !--wy
         if (sg(i,k) > 0._r8 .or. sg(i,k-1) > 0._r8) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._r8 .or. qg(i,k-1) > 0._r8) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         !++wy
         if (sg_avgd(i,k) > 0._r8 .or. sg_avgd(i,k-1) > 0._r8) &
            sdifr_avgd = abs((sg_avgd(i,k)-sg_avgd(i,k-1))/max(sg_avgd(i,k-1),sg_avgd(i,k)))
         if (qg_avgd(i,k) > 0._r8 .or. qg_avgd(i,k-1) > 0._r8) &
            qdifr_avgd = abs((qg_avgd(i,k)-qg_avgd(i,k-1))/max(qg_avgd(i,k-1),qg_avgd(i,k)))
         !--wy
         if (sdifr > 1.E-6_r8) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
         end if
         !++wy
         if (sdifr_avgd > 1.E-6_r8) then
            shat_avgd(i,k) = log(sg_avgd(i,k-1)/sg_avgd(i,k))*sg_avgd(i,k-1)*sg_avgd(i,k)/(sg_avgd(i,k-1)-sg_avgd(i,k))
         else
            shat_avgd(i,k) = 0.5_r8* (sg_avgd(i,k)+sg_avgd(i,k-1))
         end if
         !--wy
         if (qdifr > 1.E-6_r8) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
         end if
         !++wy
         if (qdifr_avgd > 1.E-6_r8) then
            qhat_avgd(i,k) = log(qg_avgd(i,k-1)/qg_avgd(i,k))*qg_avgd(i,k-1)*qg_avgd(i,k)/(qg_avgd(i,k-1)-qg_avgd(i,k))
         else
            qhat_avgd(i,k) = 0.5_r8* (qg_avgd(i,k)+qg_avgd(i,k-1))
         end if
         !--wy
      end do
   end do
!
! obtain cloud properties.
!

   call cldprp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg     , &
               cmeg    ,maxg    ,lelg    ,jt      ,jlcl    , &
               maxg    ,j0      ,jd      ,rl      ,lengath , &
               rgas    ,grav    ,cpres   ,msg     , &
!<songxl 2014-11-20------
!               pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  ,landfracg)
               pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  ,landfracg, &
               qnlg    ,qnig    ,qldeg   ,qideg   ,qncdeg  , qnideg   , &
               qraing  ,qsnowg  ,qnrg    ,qnsg    ,wug     , sprdg    , &
               aero    ,autolmg ,accrlmg ,bergnmg ,fhtimmg , fhtctmg  , &
               fhmlmg  ,hmpimg  ,accslmg ,dlfmg   ,autolng , accrlng  , &
               bergnng ,fhtimng ,fhtctng ,fhmlng  ,accslng , activng  , &
               dlfng   ,autoimg ,accsimg ,difmg   ,nucling , autoing  , &
               accsing ,hmping  ,difng   ,cmelg   ,cmeig   , trspcmg  , &
               trspcng ,trspimg ,trsping ,lambdadpcug, mudpcug, frzg  , &
               qliqg   ,qiceg   ,dcapeg  ,qhat ) 

   do i = 1,lengath
      capeg(i) = capeg(i)+ dcapeg(i)
      !++wy
      capeg_avgd(i) = capeg_avgd(i)+dcapeg(i)
      !--wy 
   end do

!>songxl 2014-11-20------

!++wy
    if (zm_stc) then
        do k = 1,pver
            do i = 1,lengath
                call qmmr_hPa(t_avgd(i,k), p(i,k), es_avgd(i), qs_avgd(i,k))
                if ( p(i,k)-es_avgd(i) <= 0._r8 ) then
                    qs_avgd(i,k) = 1.0_r8
                end if
            end do
        end do
    end if
!--wy

!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du   (i,k) = du   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu   (i,k) = eu   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed   (i,k) = ed   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug  (i,k) = cug  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg (i,k) = cmeg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg(i,k) = rprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
!<songxl 2014-11-20----------
         sprdg(i,k) = sprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         frzg (i,k) = frzg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
!>songxl 2014-11-20----------
         evpg (i,k) = evpg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   !++wy
   if (.not.zm_stc) then
       call closure(lchnk   , &
                    qg      ,tg      ,pg      ,zg      ,sg      , &
                    tpg     ,qs      ,qu      ,su      ,mc      , &
                    du      ,mu      ,md      ,qd      ,sd      , &
                    qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                    qlg     ,dsubcld ,mb      ,capeg   ,tlg     , &
                    lclg    ,lelg    ,jt      ,maxg    ,1       , &
                    lengath ,rgas    ,grav    ,cpres   ,rl      , &
                    msg     ,capelmt    )
   else
       call closure(lchnk     , &
                    qg_avgd   ,tg_avgd   ,pg      ,zg           ,sg_avgd , &
                    tpg       ,qs_avgd   ,qu      ,su           ,mc      , &
                    du        ,mu        ,md      ,qd           ,sd      , &
                    qhat_avgd ,shat_avgd ,dp      ,qstpg        ,zfg     , &
                    qlg       ,dsubcld   ,mb      ,capeg_avgd   ,tlg     , &
                    lclg      ,lelg      ,jt      ,maxg         ,1       , &
                    lengath   ,rgas      ,grav    ,cpres        ,rl      , &
                    msg       ,capelmt   )
   end if
   !--wy

!
! limit cloud base mass flux to theoretical upper bound.
!
   do i=1,lengath
      mumax(i) = 0
   end do
   do k=msg + 2,pver
      do i=1,lengath
        mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
   end do

   do i=1,lengath
      if (mumax(i) > 0._r8) then
         mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
      else
         mb(i) = 0._r8
      endif
   end do
   ! If no_deep_pbl = .true., don't allow convection entirely 
   ! within PBL (suggestion of Bjorn Stevens, 8-2000)

   if (no_deep_pbl) then
      do i=1,lengath
         if (zm(ideep(i),jt(i)) < pblh(ideep(i))) mb(i) = 0
      end do
   end if

!++wy 
!for stochastic deep convection scheme 
   if (zm_stc) then
       do i = 1, lengath
            ENSNC = mb(i)/MeanCMF
    
            DELTAR = (sqrt(CLMASSRANGE(2)/10.0/pi)-sqrt(CLMASSRANGE(1)/10.0/pi))/nbins
    
   !Initialize random number index
            irannum=0
             
            !Cloud radius loop
            LOOPRADIUS = 0._r8
            LOOPMASS   = CLMASSRANGE(1)
    
            do while (LOOPMASS < CLMASSRANGE(2))
                !Increment loop cloud mass flux
                LOOPRADIUS = LOOPRADIUS + DELTAR  ! m
                LOOPMASS   = LOOPRADIUS**2*10._r8*pi
                LOOPMASS1  = (LOOPRADIUS-DELTAR)**2*10._r8*pi
                DELTAM     = LOOPMASS-LOOPMASS1
    
                !Propability of getting a cloud for the current bin
                PROB = EXP(-1.0_r8*LOOPMASS/MeanCMF)*DELTAM/MeanCMF
    
                !In stochastic mode: Draw random number and check cloud launch
                !Increment random number index
                irannum=irannum+1
                !Is a cloud of this size required?
    
                launchcl=rannumbersg(i,irannum).LE.PROB
    
                !Called only at points where a new cloud is to be launched
                if (launchcl) then!Set Cloud radius to bin center
    
                    !Assign mass flux
                    MASS = LOOPMASS-0.5*DELTAM
                    REQMF = MASS*ENSNC
                    sum_mb(i)      = sum_mb(i)+REQMF
                end if
            end do ! radius loop
        end do !i loop
    end if
!--wy

   do k=msg+1,pver
      do i=1,lengath
         
         !++wy mp
         if (zm_stc) then
             mb(i) = sum_mb(i)
             if (mumax(i) > 0._r8) then
                 mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
             else
                 mb(i) = 0._r8
             end if
         end if
         !--wy mp

         mu   (i,k)  = mu   (i,k)*mb(i)
         md   (i,k)  = md   (i,k)*mb(i)
         mc   (i,k)  = mc   (i,k)*mb(i)
         du   (i,k)  = du   (i,k)*mb(i)
         eu   (i,k)  = eu   (i,k)*mb(i)
         ed   (i,k)  = ed   (i,k)*mb(i)
         cmeg (i,k)  = cmeg (i,k)*mb(i)
         rprdg(i,k)  = rprdg(i,k)*mb(i)
         cug  (i,k)  = cug  (i,k)*mb(i)
         evpg (i,k)  = evpg (i,k)*mb(i)
         pflxg(i,k+1)= pflxg(i,k+1)*mb(i)*100._r8/grav
         !<songxl 2014-11-20----------
         sprdg(i,k)  = sprdg(i,k)*mb(i)
         frzg (i,k)  = frzg (i,k)*mb(i)
         if ( zm_microp .and. mb(i).eq.0._r8) then
            qlg (i,k) = 0._r8
            qliqg (i,k) = 0._r8
            qiceg (i,k) = 0._r8
            qraing(i,k) = 0._r8
            qsnowg(i,k) = 0._r8
            wug(i,k) = 0._r8
            qnlg (i,k) = 0._r8
            qnig (i,k) = 0._r8
            qnrg (i,k) = 0._r8
            qnsg (i,k) = 0._r8

            autolmg(i,k) = 0._r8
            accrlmg(i,k) = 0._r8
            bergnmg(i,k) = 0._r8
            fhtimmg(i,k) = 0._r8
            fhtctmg(i,k) = 0._r8
            fhmlmg (i,k) = 0._r8
            hmpimg (i,k) = 0._r8
            accslmg(i,k) = 0._r8
            dlfmg  (i,k) = 0._r8

            autolng(i,k) = 0._r8
            accrlng(i,k) = 0._r8
            bergnng(i,k) = 0._r8
            fhtimng(i,k) = 0._r8
            fhtctng(i,k) = 0._r8
            fhmlng (i,k) = 0._r8
            accslng(i,k) = 0._r8
            activng(i,k) = 0._r8
            dlfng  (i,k) = 0._r8
            cmelg  (i,k) = 0._r8


            autoimg(i,k) = 0._r8
            accsimg(i,k) = 0._r8
            difmg  (i,k) = 0._r8
            cmeig  (i,k) = 0._r8

            nucling(i,k) = 0._r8
            autoing(i,k) = 0._r8
            accsing(i,k) = 0._r8
            hmping (i,k) = 0._r8
            difng  (i,k) = 0._r8

            trspcmg(i,k) = 0._r8
            trspcng(i,k) = 0._r8
            trspimg(i,k) = 0._r8
            trsping(i,k) = 0._r8
         end if
!>songxl 2014-11-20----------

      end do
   end do

   !
   ! compute temperature and moisture changes due to convection.
   !
      call q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
   !<songxl 2014-11-20---------
   !                 mu      ,md      ,sd      ,qd      ,qlg     , &
                    mu      ,md      ,sd      ,qd      ,qldeg   , &
   !>songxl 2014-11-20--------
                    dsubcld ,jt      ,maxg    ,1       ,lengath , &
                    cpres   ,rl      ,msg     ,          &
   !<songxl 2014-11-20------------
   !                 dlg     ,evpg    ,cug     )
                    dlg     ,evpg    ,cug     ,qideg   ,dig     , &
                    qncdeg  ,qnideg  ,dnlg    ,dnig    ,frzg    ) 
   !>songxl 2014-11-20------------

   !++wy
   sum_mb = 0._r8
   !--wy
       
!<songxl 2011-09-20-------------------------------------
    if(momentm) then
      call moment(dudt,dvdt,dpdx,dpdy,                         &
                 dp,ed,eu,mc,md,mu,mb,                         &
                 pg,qd,qu,qhat,sd,su,shat,tg,ug,vg,zg,zfg,     &
                 dsubcld,maxg,jd,jt,                           &
                 msg,2._r8*delt,grav,cpres,rgas,pver,1,lengath,pcols)
      end if
!-songxl 2011-09-20-------------------------------------

!
! gather back temperature and mixing ratio.
!
  if (zm_microp) then
    do k = msg + 1,pver
!DIR$ CONCURRENT
      do i = 1,lengath
         if (dqdt(i,k)*2._r8*delt+qg(i,k)<0._r8) then
            negadq = (dqdt(i,k)+0.5_r8*qg(i,k)/delt)/0.9999_r8
            dqdt(i,k) = dqdt(i,k)-negadq

            do kk=k,jt(i),-1
              if (negadq<0._r8) then
                 if (rprdg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                    dsdt(i,k) = dsdt(i,k) + negadq*rl/cpres
                    if (rprdg(i,kk)>sprdg(i,kk)) then
                       if(rprdg(i,kk)-sprdg(i,kk)<-negadq*dp(i,k)/dp(i,kk)) then
                         dsdt(i,k) = dsdt(i,k) + (negadq+ (rprdg(i,kk)-sprdg(i,kk))*dp(i,kk)/dp(i,k))*latice/cpres
                         sprdg(i,kk) = negadq*dp(i,k)/dp(i,kk)+rprdg(i,kk)
                       end if
                    else
                       sprdg(i,kk) = sprdg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                       dsdt(i,k) = dsdt(i,k) + negadq*latice/cpres
                    end if
                    rprdg(i,kk) = rprdg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                    negadq = 0._r8
                 else
                    negadq = rprdg(i,kk)*dp(i,kk)/dp(i,k)+negadq
                    dsdt(i,k) = dsdt(i,k) - rprdg(i,kk)*rl/cpres*dp(i,kk)/dp(i,k)
                    if (rprdg(i,kk)>sprdg(i,kk)) then
                       dsdt(i,k) = dsdt(i,k) - sprdg(i,kk)*latice/cpres*dp(i,kk)/dp(i,k)
                       sprdg(i,kk) = 0._r8
                    else
                       dsdt(i,k) = dsdt(i,k) -rprdg(i,kk)*latice/cpres*dp(i,kk)/dp(i,k)
                       sprdg(i,kk)= sprdg(i,kk)- rprdg(i,kk)
                    end if
                    rprdg(i,kk) = 0._r8
                 end if

                 if (dlg(i,kk)>dig(i,kk)) then
                    doliq= .true.
                 else
                    doliq= .false.
                 end if

                 if (negadq<0._r8) then
                   if (doliq) then
                     if (dlg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                        dsdt(i,k) = dsdt(i,k) + negadq*rl/cpres
                        dnlg(i,kk) = dnlg(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dlg(i,kk))
                        dlg(i,kk)  = dlg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                        negadq = 0._r8
                     else
                        negadq = negadq + dlg(i,kk)*dp(i,kk)/dp(i,k)
                        dsdt(i,k) = dsdt(i,k) - dlg(i,kk)*dp(i,kk)/dp(i,k)*rl/cpres
                        dlg(i,kk) = 0._r8
                        dnlg(i,kk) = 0._r8
                     end if
                   else
                     if (dig(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                       dsdt(i,k) = dsdt(i,k) + negadq*(rl+latice)/cpres
                       dnig(i,kk) = dnig(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dig(i,kk))
                       dig(i,kk)  = dig(i,kk)+negadq*dp(i,k)/dp(i,kk)
                       negadq = 0._r8
                     else
                       negadq = negadq + dig(i,kk)*dp(i,kk)/dp(i,k)
                       dsdt(i,k) = dsdt(i,k) - dig(i,kk)*dp(i,kk)/dp(i,k)*(rl+latice)/cpres
                       dig(i,kk) = 0._r8
                       dnig(i,kk) = 0._r8
                     end if
                     doliq= .false.
                   end if
                 end if
                 if (negadq<0._r8 .and. doliq ) then
                   if (dlg(i,kk)> -negadq*dp(i,k)/dp(i,kk)) then
                      dsdt(i,k) = dsdt(i,k) + negadq*rl/cpres
                      dnlg(i,kk) = dnlg(i,kk)*(1._r8+negadq*dp(i,k)/dp(i,kk)/dlg(i,kk))
                      dlg(i,kk)  = dlg(i,kk)+negadq*dp(i,k)/dp(i,kk)
                      negadq = 0._r8
                   else
                      negadq = negadq + dlg(i,kk)*dp(i,kk)/dp(i,k)
                      dsdt(i,k) = dsdt(i,k) - dlg(i,kk)*dp(i,kk)/dp(i,k)*rl/cpres
                      dlg(i,kk) = 0._r8
                      dnlg(i,kk) = 0._r8
                   end if
                 end if

              end if
            end do

            if (negadq<0._r8) then
!               write(*,*) "can not fix negadq=",negadq
               dqdt(i,k) = dqdt(i,k) + negadq
            end if

          end if
      end do
   end do
  end if


   do k = msg + 1,pver
!DIR$ CONCURRENT
      do i = 1,lengath
!
! q is updated to compute net precip.
!
         q(ideep(i),k) = qh(ideep(i),k) + 2._r8*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         cme (ideep(i),k) = cmeg (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*cpres
         dlf (ideep(i),k) = dlg  (i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql  (ideep(i),k) = qlg  (i,k)
!<songxl 2011-09-20-------------------------------
         utg(ideep(i),k) = dudt(i,k)
         vtg(ideep(i),k) = dvdt(i,k)
         putg(ideep(i),k) = dpdx(i,k)
         pvtg(ideep(i),k) = dpdy(i,k)
!>songxl 2011-09-20-------------------------------
!<songxl 2014-11-20---------
         qi  (ideep(i),k) = qiceg  (i,k)
         frz(ideep(i),k)  = frzg(i,k)*latice/cpres
         sprd(ideep(i),k) = sprdg(i,k)
         wu  (ideep(i),k) = wug  (i,k)
         qliq(ideep(i),k) = qliqg (i,k)
         qice(ideep(i),k) = qiceg (i,k)
         qrain(ideep(i),k) = qraing (i,k)
         qsnow(ideep(i),k) = qsnowg (i,k)
         dif (ideep(i),k) = dig  (i,k)
         dnlf(ideep(i),k) = dnlg  (i,k)
         dnif(ideep(i),k) = dnig  (i,k)
         qnl(ideep(i),k)  = qnlg(i,k)
         qni(ideep(i),k)  = qnig(i,k)
         qnr(ideep(i),k)  = qnrg(i,k)
         qns(ideep(i),k)  = qnsg(i,k)

         autolm(ideep(i),k) = autolmg(i,k)
         accrlm(ideep(i),k) = accrlmg(i,k)
         bergnm(ideep(i),k) = bergnmg(i,k)
         fhtimm(ideep(i),k) = fhtimmg(i,k)
         fhtctm(ideep(i),k) = fhtctmg(i,k)
         fhmlm (ideep(i),k) = fhmlmg (i,k)
         hmpim (ideep(i),k) = hmpimg (i,k)
         accslm(ideep(i),k) = accslmg(i,k)
         dlfm  (ideep(i),k) = dlfmg  (i,k)

         autoln(ideep(i),k) = autolng(i,k)
         accrln(ideep(i),k) = accrlng(i,k)
         bergnn(ideep(i),k) = bergnng(i,k)
         fhtimn(ideep(i),k) = fhtimng(i,k)
         fhtctn(ideep(i),k) = fhtctng(i,k)
         fhmln (ideep(i),k) = fhmlng (i,k)
         accsln(ideep(i),k) = accslng(i,k)
         activn(ideep(i),k) = activng(i,k)
         dlfn  (ideep(i),k) = dlfng  (i,k)
         cmel  (ideep(i),k) = cmelg  (i,k)

         autoim(ideep(i),k) = autoimg(i,k)
         accsim(ideep(i),k) = accsimg(i,k)
         difm  (ideep(i),k) = difmg  (i,k)
         cmei  (ideep(i),k) = cmeig  (i,k)

         nuclin(ideep(i),k) = nucling(i,k)
         autoin(ideep(i),k) = autoing(i,k)
         accsin(ideep(i),k) = accsing(i,k)
         hmpin (ideep(i),k) = hmping (i,k)
         difn  (ideep(i),k) = difng  (i,k)

         trspcm(ideep(i),k) = trspcmg(i,k)
         trspcn(ideep(i),k) = trspcng(i,k)
         trspim(ideep(i),k) = trspimg(i,k)
         trspin(ideep(i),k) = trsping(i,k)

         lambdadpcu(ideep(i),k) = lambdadpcug(i,k)
         mudpcu(ideep(i),k)     = mudpcug(i,k)
!>songxl 2014-11-20---------
      end do
   end do

!<songxl 2014-11-20------------
   do k = msg + 1,pver
      do i = 1,ncol

         !convert it from units of "kg/kg" to "g/m3"

         if(k.lt.pver) then
            qice (i,k) = 0.5_r8*(qice(i,k)+qice(i,k+1))
            qliq (i,k) = 0.5_r8*(qliq(i,k)+qliq(i,k+1))
            qrain (i,k) = 0.5_r8*(qrain(i,k)+qrain(i,k+1))
            qsnow (i,k) = 0.5_r8*(qsnow(i,k)+qsnow(i,k+1))
            qni (i,k) = 0.5_r8*(qni(i,k)+qni(i,k+1))
            qnl (i,k) = 0.5_r8*(qnl(i,k)+qnl(i,k+1))
            qnr (i,k) = 0.5_r8*(qnr(i,k)+qnr(i,k+1))
            qns (i,k) = 0.5_r8*(qns(i,k)+qns(i,k+1))
            wu(i,k)   = 0.5_r8*(wu(i,k)+wu(i,k+1))
         end if

         if (t(i,k).gt. 273.15_r8 .and. t(i,k-1).le.273.15_r8) then
             qice (i,k-1) = qice (i,k-1) + qice (i,k)
             qice (i,k) = 0._r8
             qni (i,k-1) = qni (i,k-1) + qni (i,k)
             qni (i,k) = 0._r8
             qsnow (i,k-1) = qsnow (i,k-1) + qsnow (i,k)
             qsnow (i,k) = 0._r8
             qns (i,k-1) = qns (i,k-1) + qns (i,k)
             qns (i,k) = 0._r8
         end if

         qice (i,k) = qice(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         qliq (i,k) = qliq(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         qrain (i,k) = qrain(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         qsnow (i,k) = qsnow(i,k) * pap(i,k)/t(i,k)/rgas *1000._r8
         qni (i,k) = qni(i,k) * pap(i,k)/t(i,k)/rgas
         qnl (i,k) = qnl(i,k) * pap(i,k)/t(i,k)/rgas
         qnr (i,k) = qnr(i,k) * pap(i,k)/t(i,k)/rgas
         qns (i,k) = qns(i,k) * pap(i,k)/t(i,k)/rgas
      end do
   end do

!>songxl 2014-11-20------------

!
!DIR$ CONCURRENT
   do i = 1,lengath
      jctop(ideep(i)) = jt(i)
!++bee
      jcbot(ideep(i)) = maxg(i)
!--bee
      pflx(ideep(i),pverp) = pflxg(i,pverp)
      dcape(ideep(i)) =  dcapeg(i)            !songxl 2014-10-10
   end do

! Compute precip by integrating change in water vapor minus detrained cloud water
   do k = pver,msg + 1,-1
      do i = 1,ncol
!<songxl 2014-11-20---------
!         prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k)) - dpp(i,k)*dlf(i,k)*2*delt
          prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k)) - dpp(i,k)*(dlf(i,k)+dif(i,k))*2._r8*delt
!>songxl 2014-11-20---------
      end do
   end do

! obtain final precipitation rate in m/s.
   do i = 1,ncol
      prec(i) = rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
   end do

! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
   do k = 1, pver
      do i = 1, ncol
!<songxl 2014-11-20------
!         rliq(i) = rliq(i) + dlf(i,k)*dpp(i,k)/gravit
          rliq(i) = rliq(i) + (dlf(i,k)+dif(i,k))*dpp(i,k)/gravit
          rice(i) = rice(i) + dif(i,k)*dpp(i,k)/gravit
!>songxl 2014-11-20------
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8
   rice(:ncol) = rice(:ncol) /1000._r8

   return

end subroutine zm_convr

!===============================================================================
subroutine zm_conv_evap(ncol,lchnk, &
     t,pmid,pdel,q, &
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat,  &
!<songxl 2014-11-20--------
!     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow )
     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow, prdsnow, old_snow )
!>songxl 2014-11-20--------      

!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

    use wv_saturation,  only: qsat
    use phys_grid, only: get_rlat_all_p

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol, lchnk             ! number of columns and chunk index
    real(r8),intent(in), dimension(pcols,pver) :: t          ! temperature (K)
    real(r8),intent(in), dimension(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    real(r8),intent(in), dimension(pcols,pver) :: pdel       ! layer thickness (Pa)
    real(r8),intent(in), dimension(pcols,pver) :: q          ! water vapor (kg/kg)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwprd ! Heating rate of snow production
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    


    real(r8), intent(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step

    real(r8), intent(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    real(r8), intent(out)   :: snow(pcols)        ! Convective-scale snowfall rate
!<songxl 2014-11-20-------------
    real(r8), intent(in   ) :: prdsnow(pcols,pver)! snow production (kg/ks/s)
    logical,  intent(in)    :: old_snow           ! true for old estimate of snow production
!>songxl 2014-11-20-------------
!
!---------------------------Local storage-------------------------------

    real(r8) :: es    (pcols,pver)    ! Saturation vapor pressure
    real(r8) :: fice   (pcols,pver)    ! ice fraction in precip production
    real(r8) :: fsnow_conv(pcols,pver) ! snow fraction in precip production
    real(r8) :: qs   (pcols,pver)    ! saturation specific humidity
    real(r8),intent(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),intent(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8),intent(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    real(r8),intent(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
    real(r8) :: work1                  ! temp variable (pjr)
    real(r8) :: work2                  ! temp variable (pjr)

    real(r8) :: evpvint(pcols)         ! vertical integral of evaporation
    real(r8) :: evpprec(pcols)         ! evaporation of precipitation (kg/kg/s)
    real(r8) :: evpsnow(pcols)         ! evaporation of snowfall (kg/kg/s)
    real(r8) :: snowmlt(pcols)         ! snow melt tendency in layer
    real(r8) :: flxsntm(pcols)         ! flux of snow into layer, after melting

    real(r8) :: evplimit               ! temp variable for evaporation limits
    real(r8) :: rlat(pcols)
!<songxl 2014-11-20----
    real(r8) :: dum 
    real(r8) :: omsm
!>songxl 2014-11-20---

    integer :: i,k                     ! longitude,level indices


!-----------------------------------------------------------------------

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._r8

! determine saturation vapor pressure
    call qsat(t(1:ncol, 1:pver), pmid(1:ncol, 1:pver), &
         es(1:ncol, 1:pver), qs(1:ncol, 1:pver))

! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldfrc_fice(ncol, t, fice, fsnow_conv)

! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0._r8
    flxsnow(:ncol,1) = 0._r8
    evpvint(:ncol)   = 0._r8
    omsm=0.9999_r8

    do k = 1, pver
       do i = 1, ncol

! Melt snow falling into layer, if necessary. 
!<songxl 2014-11-20------------
        if( old_snow ) then
          if (t(i,k) > tmelt) then
             flxsntm(i) = 0._r8
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
          end if
        else
! make sure melting snow doesn't reduce temperature below threshold
          if (t(i,k) > tmelt) then
              dum = -latice/cpres*flxsnow(i,k)*gravit/pdel(i,k)*deltat
              if (t(i,k) + dum .le. tmelt) then
                dum = (t(i,k)-tmelt)*cpres/latice/deltat
                dum = dum/(flxsnow(i,k)*gravit/pdel(i,k))
                dum = max(0._r8,dum)
                dum = min(1._r8,dum)
              else
                dum = 1._r8
              end if
              dum = dum*omsm
              flxsntm(i) = flxsnow(i,k)*(1.0_r8-dum)
              snowmlt(i) = dum*flxsnow(i,k)*gravit/ pdel(i,k)
           else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
           end if
         end if
!>songxl 2014-10-10--------------------

! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._r8 - q(i,k)/qs(i,k), 0._r8)

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))
!**********************************************************
!!          evpprec(i) = 0.    ! turn off evaporation for now
!**********************************************************

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qs
          evplimit   = max(0._r8, (qs(i,k)-q(i,k)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
          evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = min(evplimit, evpprec(i))

          if( .not.old_snow ) then
            evpprec(i) = max(0._r8, evpprec(i))
            evpprec(i) = evpprec(i)*omsm
          end if

! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(i,k) > 0._r8) then
!            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
!            prevent roundoff problems
             work1 = min(max(0._r8,flxsntm(i)/flxprec(i,k)),1._r8)
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._r8
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
!pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
!pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
! the small amount added to flxprec in the work1 expression has been increased from 
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.

!<songxl 2014-11-20--------
      if( old_snow ) then  
!>songxl 2014-11-20--------

#ifdef PERGRO
          work1 = min(max(0._r8,flxsnow(i,k)/(flxprec(i,k)+8.64e-11_r8)),1._r8)
#else
          if (flxprec(i,k).gt.0._r8) then
             work1 = min(max(0._r8,flxsnow(i,k)/flxprec(i,k)),1._r8)
          else
             work1 = 0._r8
          endif
#endif
          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i).gt.0._r8) work2 = 0._r8
!         work2 = fsnow_conv(i,k)
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k)*work2*latice
!<songxl 2014-11-20----------
          tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*latice
       else
          ntsnprd(i,k) = prdsnow(i,k) - min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i))
          tend_s_snwprd  (i,k) = prdsnow(i,k)*latice
          tend_s_snwevmlt(i,k) = -min(flxsnow(i,k)*gravit/pdel(i,k), evpsnow(i)+snowmlt(i) )*latice
       end if        
!>songxl 2014-11-20----------

! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

! protect against rounding error
          flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
          flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)
! more protection (pjr)
!         flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

! heating (cooling) and moistening due to evaporation 
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
!<songxl 2014-11-20--------
!          tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice        
          if( old_snow ) then
             tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          else
             tend_s(i,k)   =-evpprec(i)*latvap + tend_s_snwevmlt(i,k)
          end if 
!>songxl 2014-11-20--------
          tend_q(i,k) = evpprec(i)
       end do
    end do

! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000._r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000._r8

!**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
!**********************************************************

  end subroutine zm_conv_evap



subroutine convtran(lchnk   , &
                    doconvtran,q       ,ncnst   ,mu      ,md      , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry   ,dt )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index

   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dt               ! 2 delta t (model time increment)

! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered tracer array
   real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
   real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
   real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
   real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces

   real(r8) total(pcols)
   real(r8) negadt,qtmp

!-----------------------------------------------------------------------
!
   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

!<songxl test----------
       if (zm_microp) then
         do i = il1g,il2g
           do k = jt(i),mx(i)
             if (dcondt(i,k)*dt+const(i,k)<0._r8) then
                negadt = dcondt(i,k)+const(i,k)/dt
                dcondt(i,k) = -const(i,k)/dt
                do kk= k+1, mx(i)
                  if (negadt<0._r8 .and. dcondt(i,kk)*dt+const(i,kk)>0._r8 ) then
                    qtmp = dcondt(i,kk)+negadt*dptmp(i,k)/dptmp(i,kk)
                    if (qtmp*dt+const(i,kk)>0._r8) then
                      dcondt(i,kk)= qtmp
                      negadt=0._r8
                    else
                      negadt= negadt+(const(i,kk)/dt+dcondt(i,kk))*dptmp(i,kk)/dptmp(i,k)
                      dcondt(i,kk)= -const(i,kk)/dt
                    end if

                  end if
                end do
                do kk= k-1, jt(i), -1
                  if (negadt<0._r8 .and. dcondt(i,kk)*dt+const(i,kk)>0._r8 ) then
                    qtmp = dcondt(i,kk)+negadt*dptmp(i,k)/dptmp(i,kk)
                    if (qtmp*dt+const(i,kk)>0._r8) then
                      dcondt(i,kk)= qtmp
                      negadt=0._r8
                    else
                      negadt= negadt+(const(i,kk)/dt+dcondt(i,kk))*dptmp(i,kk)/dptmp(i,k)
                      dcondt(i,kk)= -const(i,kk)/dt
                    end if
                  end if
                end do

                if (negadt<0._r8) then
!                   write(*,*) "can not fix negadt=",negadt
                   dcondt(i,k) = dcondt(i,k) + negadt
                end if
             end if
           end do
         end do
       end if
!>songxl test----------


! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   return
end subroutine convtran

!=========================================================================================

subroutine momtran(lchnk, ncol, &
                    domomtran,q       ,ncnst   ,mu      ,md    , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,dqdt    ,pguall     ,pgdall, icwu, icwd, dt, seten    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of momentum
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! Based on the convtran subroutine by P. Rasch
! <Also include any applicable external references.> 
! 
! Author: J. Richter and P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: domomtran(ncnst)      ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Wind array
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in)    :: dt                       !  time step in seconds : 2*delta_t

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index



! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer kkm1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index
   integer ii                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered wind array
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) momcu                ! constant for updraft pressure gradient term
   real(r8) momcd                ! constant for downdraft pressure gradient term
   real(r8) sum                  ! sum
   real(r8) sum2                  ! sum2
 
   real(r8) mududp(pcols,pver) ! working variable
   real(r8) mddudp(pcols,pver)     ! working variable

   real(r8) pgu(pcols,pver)      ! Pressure gradient term for updraft
   real(r8) pgd(pcols,pver)      ! Pressure gradient term for downdraft

   real(r8),intent(out) ::  pguall(pcols,pver,ncnst)      ! Apparent force from  updraft PG
   real(r8),intent(out) ::  pgdall(pcols,pver,ncnst)      ! Apparent force from  downdraft PG

   real(r8),intent(out) ::  icwu(pcols,pver,ncnst)      ! In-cloud winds in updraft
   real(r8),intent(out) ::  icwd(pcols,pver,ncnst)      ! In-cloud winds in downdraft

   real(r8),intent(out) ::  seten(pcols,pver) ! Dry static energy tendency
   real(r8)                 gseten(pcols,pver) ! Gathered dry static energy tendency

   real(r8)  mflux(pcols,pverp,ncnst)   ! Gathered momentum flux

   real(r8)  wind0(pcols,pver,ncnst)       !  gathered  wind before time step
   real(r8)  windf(pcols,pver,ncnst)       !  gathered  wind after time step
   real(r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2
   

!-----------------------------------------------------------------------
!

! Initialize outgoing fields
   pguall(:,:,:)     = 0.0_r8
   pgdall(:,:,:)     = 0.0_r8
! Initialize in-cloud winds to environmental wind
   icwu(:ncol,:,:)       = q(:ncol,:,:)
   icwd(:ncol,:,:)       = q(:ncol,:,:)

! Initialize momentum flux and  final winds
   mflux(:,:,:)       = 0.0_r8
   wind0(:,:,:)         = 0.0_r8
   windf(:,:,:)         = 0.0_r8

! Initialize dry static energy

   seten(:,:)         = 0.0_r8
   gseten(:,:)         = 0.0_r8

! Define constants for parameterization

   momcu = 0.4_r8
   momcd = 0.4_r8

   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each wind component
   do m = 1, ncnst                    !start at m = 1 to transport momentum
      if (domomtran(m)) then

! Gather up the winds and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
                wind0(i,k,m) = const(i,k)
            end do
         end do


! From now on work only with gathered data

! Interpolate winds to interfaces

         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g

               ! use arithmetic mean
               chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do


!
! Pressure Perturbation Term
! 

      !Top boundary:  assume mu is zero 

         k=1
         pgu(:il2g,k) = 0.0_r8
         pgd(:il2g,k) = 0.0_r8

         do k=2,pver-1
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
            
               !interior points

               mududp(i,k) =  ( mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  mu(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgu(i,k) = - momcu * 0.5_r8 * mududp(i,k)
                           

               mddudp(i,k) =  ( md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  md(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgd(i,k) = - momcd * 0.5_r8 * mddudp(i,k)


            end do
         end do

       ! bottom boundary 
       k = pver
       km1 = max(1,k-1)
       do i=il1g,il2g

          mududp(i,k) =   mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)
          pgu(i,k) = - momcu *  mududp(i,k)
          
          mddudp(i,k) =   md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) 

          pgd(i,k) = - momcd * mddudp(i,k)
          
       end do
       

!
! In-cloud velocity calculations
!

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         kkm1 = max(1,kk-1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then
                 
               conu(i,kk) = (+eu(i,kk)*const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-ed(i,km1)*const(i,km1)*dp(i,km1))-pgd(i,km1)*dp(i,km1)/md(i,k)
            endif

                        
         end do



! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkm1 = max(1,kk-1)
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp > mbsth) then
            
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eu(i,kk)* &
                                  const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
               endif
            end do

         end do


! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                            
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-ed(i,km1)*const(i,km1) &
                                  *dp(i,km1)-pgd(i,km1)*dp(i,km1) )/md(i,k)

               endif
            end do
         end do


         sum = 0._r8
         sum2 = 0._r8


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               ii = ideep(i)
	
! version 1 hard to check for roundoff errors
               dcondt(i,k) =  &
                           +(mu(i,kp1)* (conu(i,kp1)-chat(i,kp1)) &
                           -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                           +md(i,kp1)* (cond(i,kp1)-chat(i,kp1)) &
                           -md(i,k)*   (cond(i,k)-chat(i,k)) &
                          )/dp(i,k)

            end do
         end do

  ! dcont for bottom layer
          !
          !DIR$ NOINTERCHANGE
          do k = kbm,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
                if (k == mx(i)) then

                   ! version 1
                   dcondt(i,k) = (1._r8/dp(i,k))*   &  
                        (-mu(i,k)*(conu(i,k)-chat(i,k)) &
                        -md(i,k)*(cond(i,k)-chat(i,k)) &
                        )
                end if
             end do
          end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8

         do k = 1,pver
            do i = il1g,il2g
               ii = ideep(i)
               dqdt(ii,k,m) = dcondt(i,k)
    ! Output apparent force on the mean flow from pressure gradient
               pguall(ii,k,m) = -pgu(i,k)
               pgdall(ii,k,m) = -pgd(i,k)
               icwu(ii,k,m)   =  conu(i,k)
               icwd(ii,k,m)   =  cond(i,k)
            end do
         end do

          ! Calculate momentum flux in units of mb*m/s2 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                mflux(i,k,m) = &
                     -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                     -md(i,k)*   (cond(i,k)-chat(i,k))
             end do
          end do


          ! Calculate winds at the end of the time step 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                km1 = max(1,k-1)
                kp1 = k+1
                windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)

             end do
          end do

       end if      ! for domomtran
   end do

 ! Need to add an energy fix to account for the dissipation of kinetic energy
    ! Formulation follows from Boville and Bretherton (2003)
    ! formulation by PJR

    do k = ktm,pver
       km1 = max(1,k-1)
       kp1 = min(pver,k+1)
       do i = il1g,il2g

          ii = ideep(i)

          ! calculate the KE fluxes at top and bot of layer 
          ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
          utop = (wind0(i,k,1)+wind0(i,km1,1))/2._r8
          vtop = (wind0(i,k,2)+wind0(i,km1,2))/2._r8
          ubot = (wind0(i,kp1,1)+wind0(i,k,1))/2._r8
          vbot = (wind0(i,kp1,2)+wind0(i,k,2))/2._r8
          fket = utop*mflux(i,k,1)   + vtop*mflux(i,k,2)    ! top of layer
          fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2)  ! bot of layer

          ! divergence of these fluxes should give a conservative redistribution of KE
          ketend_cons = (fket-fkeb)/dp(i,k)

          ! tendency in kinetic energy resulting from the momentum transport
          ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))*0.5_r8/dt

          ! the difference should be the dissipation
          gset2 = ketend_cons - ketend
          gseten(i,k) = gset2

       end do

    end do

    ! Scatter dry static energy to full array
    do k = 1,pver
       do i = il1g,il2g
          ii = ideep(i)
          seten(ii,k) = gseten(i,k)

       end do
    end do

   return
end subroutine momtran

!=========================================================================================

subroutine buoyan(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author:
! This is contributed code not fully standardized by the CCM core group.
! The documentation has been enhanced to the degree that we are able.
! Reviewed:          P. Rasch, April 1996
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)         ! perturbation temperature by pbl processes

!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)       ! parcel temperature
   real(r8), intent(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel
   real(r8), intent(out) :: tl(pcols)            ! parcel temperature at lcl
   real(r8), intent(out) :: cape(pcols)          ! convective aval. pot. energy.
   integer lcl(pcols)        !
   integer lel(pcols)        !
   integer lon(pcols)        ! level of onset of deep convection
   integer mx(pcols)         ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
   real(r8) capeten(pcols,5)     ! provisional value of cape
   real(r8) tv(pcols,pver)       !
   real(r8) tpv(pcols,pver)      !
   real(r8) buoy(pcols,pver)

   real(r8) a1(pcols)
   real(r8) a2(pcols)
   real(r8) estp(pcols)
   real(r8) pl(pcols)
   real(r8) plexp(pcols)
   real(r8) hmax(pcols)
   real(r8) hmn(pcols)
   real(r8) y(pcols)

   logical plge600(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)

   real(r8) cp
   real(r8) e
   real(r8) grav

   integer i
   integer k
   integer msg
   integer n

   real(r8) rd
   real(r8) rl
#ifdef PERGRO
   real(r8) rhd
#endif
!
!-----------------------------------------------------------------------
!
   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = pver
      knt(i) = 0
      lel(i) = pver
      mx(i) = lon(i)
      cape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:ncol,:) = t(:ncol,:) *(1._r8+1.608_r8*q(:ncol,:))/ (1._r8+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._r8

!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
#ifdef PERGRO
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
!
! Reset max moist static energy level when relative difference exceeds 1.e-4
!
         rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#else
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#endif
!
   do i = 1,ncol
      lcl(i) = mx(i)
      e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      tl(i) = 2840._r8/ (3.5_r8*log(t(i,mx(i)))-log(e)-4.805_r8) + 55._r8
      if (tl(i) < t(i,mx(i))) then
         plexp(i) = (1._r8/ (0.2854_r8* (1._r8-0.28_r8*q(i,mx(i)))))
         pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)
      else
         tl(i) = t(i,mx(i))
         pl(i) = p(i,mx(i))
      end if
   end do

!
! calculate lifting condensation level (lcl).
!
   do k = pver,msg + 2,-1
      do i = 1,ncol
         if (k <= mx(i) .and. (p(i,k) > pl(i) .and. p(i,k-1) <= pl(i))) then
            lcl(i) = k - 1
         end if
      end do
   end do
!
! if lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8
   end do
!
! initialize parcel properties in sub-cloud layer below lcl.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k > lcl(i) .and. k <= mx(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**(0.2854_r8* (1._r8-0.28_r8*q(i,mx(i))))
!
! buoyancy is increased by 0.5 k as in tiedtke
!
!-jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
!-jjh     1                     (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))*(1._r8+1.608_r8*q(i,mx(i)))/ (1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do

!
! define parcel properties at lcl (i.e. level immediately above pl).
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k == lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = tl(i)* (p(i,k)/pl(i))**(0.2854_r8* (1._r8-0.28_r8*qstp(i,k)))
!              estp(i)  =exp(21.656_r8 - 5418._r8/tp(i,k))
! use of different formulas for es has about 1 g/kg difference
! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
! above giving larger qs.
            call qmmr_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
            a1(i) = cp / rl + qstp(i,k) * (1._r8+ qstp(i,k) / eps1) * rl * eps1 / &
                    (rd * tp(i,k) ** 2)
            a2(i) = .5_r8* (qstp(i,k)* (1._r8+2._r8/eps1*qstp(i,k))* &
                    (1._r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                    (rd**2*tp(i,k)**4)-qstp(i,k)* &
                    (1._r8+qstp(i,k)/eps1)*2._r8*eps1*rl/ &
                    (rd*tp(i,k)**3))
            a1(i) = 1._r8/a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = q(i,mx(i)) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
            call qmmr_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
!
! buoyancy is increased by 0.5 k in cape calculation.
! dec. 9, 1994
!-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
!
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+1.608_r8*qstp(i,k)) / (1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do
!
! main buoyancy calculation.
!
   do k = pver - 1,msg + 1,-1
      do i=1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**(0.2854_r8* (1._r8-0.28_r8*qstp(i,k)))
            call qmmr_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
            a1(i) = cp/rl + qstp(i,k)* (1._r8+qstp(i,k)/eps1)*rl*eps1/ (rd*tp(i,k)**2)
            a2(i) = .5_r8* (qstp(i,k)* (1._r8+2._r8/eps1*qstp(i,k))* &
                    (1._r8+qstp(i,k)/eps1)*eps1**2*rl*rl/ &
                    (rd**2*tp(i,k)**4)-qstp(i,k)* &
                    (1._r8+qstp(i,k)/eps1)*2._r8*eps1*rl/ &
                    (rd*tp(i,k)**3))
            a1(i) = 1._r8/a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
            call qmmr_hPa(tp(i,k), p(i,k), estp(i), qstp(i,k))
!-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
!jt            (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+1.608_r8*qstp(i,k))/(1._r8+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add
         end if
      end do
   end do

!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0._r8 .and. buoy(i,k) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
!
   return
end subroutine buoyan

subroutine cldprp(lchnk   , &
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      , &
                  cmeg    ,jb      ,lel     ,jt      ,jlcl    , &
                  mx      ,j0      ,jd      ,rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     , &
!<songxl 2014-11-20---------
!                  pflx    ,evp     ,cu      ,rprd    ,limcnv  ,landfrac)
                  pflx    ,evp     ,cu      ,rprd    ,limcnv  ,landfrac, &
                  qnc     ,qni     ,qcde    ,qide    ,qncde   ,qnide   , &
                  qrain   ,qsnow   ,qnr     ,qns     ,wu      ,sprd    , &
                  aero    ,autolm  ,accrlm  ,bergnm  ,fhtimm  ,fhtctm  , &
                  fhmlm   ,hmpim   ,accslm  ,dlfm    ,autoln  ,accrln  , &
                  bergnn  ,fhtimn  ,fhtctn  ,fhmln   ,accsln  ,activn  , &
                  dlfn    ,autoim  ,accsim  ,difm    ,nuclin  ,autoin  , &
                  accsin  ,hmpin   ,difn    ,cmel    ,cmei    ,trspcm  , &
                  trspcn  ,trspim  ,trspin  ,lambdadpcu, mudpcu, frz1  , &
                  qc      ,qi      ,dcape   ,qhat    )
 
!>songxl 2014-11-20---------
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
! 
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
!
!-----------------------------------------------------------------------

   implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier

   real(r8), intent(in) :: q(pcols,pver)         ! spec. humidity of env
   real(r8), intent(in) :: t(pcols,pver)         ! temp of env
   real(r8), intent(in) :: p(pcols,pver)         ! pressure of env
   real(r8), intent(in) :: z(pcols,pver)         ! height of env
   real(r8), intent(in) :: s(pcols,pver)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
   real(r8), intent(in) :: u(pcols,pver)         ! zonal velocity of env
   real(r8), intent(in) :: v(pcols,pver)         ! merid. velocity of env

   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac

   integer, intent(in) :: jb(pcols)              ! updraft base level
   integer, intent(in) :: lel(pcols)             ! updraft launch level
   integer, intent(out) :: jt(pcols)              ! updraft plume top
   integer, intent(out) :: jlcl(pcols)            ! updraft lifting cond level
   integer, intent(in) :: mx(pcols)              ! updraft base level (same is jb)
   integer, intent(out) :: j0(pcols)              ! level where updraft begins detraining
   integer, intent(out) :: jd(pcols)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pcols,pver)      ! interface values of dry stat energy
!<songxl 2014-11-20-----------
   real(r8), intent(in) :: qhat(pcols,pver)      ! wg grid slice of upper interface mixing ratio.
   type(zm_aero_t), intent(in) :: aero           ! aerosol object
!>songxl 2014-11-20-----------
!
! output
!
   real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
   real(r8), intent(out) :: qst(pcols,pver)      ! saturation mixing ratio of env.
   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft

!<songxl 2014-11-20------------
   real(r8), intent(out) :: wu(pcols,pver)       ! vertical velocity 
   real(r8), intent(out) :: qc(pcols,pver)       ! cloud water mixing ratio (kg/kg)
   real(r8), intent(out) :: qi(pcols,pver)       ! cloud ice mixing ratio (kg/kg)
   real(r8), intent(out) :: qnc(pcols,pver)      ! cloud water number concentration(1/kg)
   real(r8), intent(out) :: qni(pcols,pver)      ! cloud ice number concentration(1/kg)
   real(r8), intent(out) :: qcde(pcols,pver)     ! cloud water mixing ratio for detrainment (kg/kg)
   real(r8), intent(out) :: qide(pcols,pver)     ! cloud ice mixing ratio for detrainment (kg/kg)
   real(r8), intent(out) :: qncde(pcols,pver)    ! cloud water number concentration for detrainment (1/kg)
   real(r8), intent(out) :: qnide(pcols,pver)    ! cloud ice number concentration for detrainment (1/kg)
   real(r8), intent(out) :: qnr(pcols,pver)      ! rain water number concentration(1/kg)
   real(r8), intent(out) :: qns(pcols,pver)      ! snow number concentration(1/kg)
   real(r8), intent(out) :: qrain(pcols,pver)    ! rain water mixing ratio (kg/kg)
   real(r8), intent(out) :: qsnow(pcols,pver)    ! snow mixing ratio (kg/kg) 
   real(r8), intent(out) :: sprd(pcols,pver)     ! rate of production of snow at that layer
!   real(r8), intent(out) :: frz(pcols,pver)      ! rate of freezing

   ! tendency for output

   real(r8), intent(out) :: autolm(pcols,pver)   !mass tendency due to autoconversion of droplets to rain
   real(r8), intent(out) :: accrlm(pcols,pver)   !mass tendency due to accretion of droplets by rain
   real(r8), intent(out) :: bergnm(pcols,pver)   !mass tendency due to Bergeron process
   real(r8), intent(out) :: fhtimm(pcols,pver)   !mass tendency due to immersion freezing
   real(r8), intent(out) :: fhtctm(pcols,pver)   !mass tendency due to contact freezing
   real(r8), intent(out) :: fhmlm (pcols,pver)   !mass tendency due to homogeneous freezing
   real(r8), intent(out) :: hmpim (pcols,pver)   !mass tendency due to HM process
   real(r8), intent(out) :: accslm(pcols,pver)   !mass tendency due to accretion of droplets by snow
   real(r8), intent(out) :: dlfm  (pcols,pver)   !mass tendency due to detrainment of droplet

   real(r8), intent(out) :: autoln(pcols,pver)   !num tendency due to autoconversion of droplets to rain
   real(r8), intent(out) :: accrln(pcols,pver)   !num tendency due to accretion of droplets by rain
   real(r8), intent(out) :: bergnn(pcols,pver)   !num tendency due to Bergeron process
   real(r8), intent(out) :: fhtimn(pcols,pver)   !num tendency due to immersion freezing
   real(r8), intent(out) :: fhtctn(pcols,pver)   !num tendency due to contact freezing
   real(r8), intent(out) :: fhmln (pcols,pver)   !num tendency due to homogeneous freezing
   real(r8), intent(out) :: accsln(pcols,pver)   !num tendency due to accretion of droplets by snow
   real(r8), intent(out) :: activn(pcols,pver)   !num tendency due to droplets activation
   real(r8), intent(out) :: dlfn  (pcols,pver)   !num tendency due to detrainment of droplet

   real(r8), intent(out) :: autoim(pcols,pver)   !mass tendency due to autoconversion of cloud ice to snow
   real(r8), intent(out) :: accsim(pcols,pver)   !mass tendency due to accretion of cloud ice by snow
   real(r8), intent(out) :: difm  (pcols,pver)   !mass tendency due to detrainment of cloud ice

   real(r8), intent(out) :: nuclin(pcols,pver)   !num tendency due to ice nucleation
   real(r8), intent(out) :: autoin(pcols,pver)   !num tendency due to autoconversion of cloud ice to snow
   real(r8), intent(out) :: accsin(pcols,pver)   !num tendency due to accretion of cloud ice by snow
   real(r8), intent(out) :: hmpin (pcols,pver)   !num tendency due to HM process 
   real(r8), intent(out) :: difn  (pcols,pver)   !num tendency due to detrainment of cloud ice

   real(r8), intent(out) :: trspcm(pcols,pver)   !mass tendency of droplets due to convective transport
   real(r8), intent(out) :: trspcn(pcols,pver)   !num tendency of droplets due to convective transport
   real(r8), intent(out) :: trspim(pcols,pver)   !mass tendency of ice crystal due to convective transport
   real(r8), intent(out) :: trspin(pcols,pver)   !num tendency of ice crystal due to convective transport

   real(r8), intent(inout) :: lambdadpcu(pcols,pver) ! slope of cloud liquid size distr
   real(r8), intent(inout) :: mudpcu(pcols,pver)     ! width parameter of droplet size distr

   real(r8) :: dcape(pcols)          ! CAPE change due to freezing heating (J/kg)
!>songxl 2014-11-20------------

   real(r8) rd                   ! gas constant for dry air
   real(r8) grav                 ! gravity
   real(r8) cp                   ! heat capacity of dry air

!
! Local workspace
!
   real(r8) gamma(pcols,pver)
   real(r8) dz(pcols,pver)
   real(r8) iprm(pcols,pver)
   real(r8) hu(pcols,pver)
   real(r8) hd(pcols,pver)
   real(r8) eps(pcols,pver)
   real(r8) f(pcols,pver)
   real(r8) k1(pcols,pver)
   real(r8) i2(pcols,pver)
   real(r8) ihat(pcols,pver)
   real(r8) i3(pcols,pver)
   real(r8) idag(pcols,pver)
   real(r8) i4(pcols,pver)
   real(r8) qsthat(pcols,pver)
   real(r8) hsthat(pcols,pver)
   real(r8) gamhat(pcols,pver)
   real(r8) cu(pcols,pver)
   real(r8) evp(pcols,pver)
   real(r8) cmeg(pcols,pver)
   real(r8) qds(pcols,pver)
! RBN For c0mask
   real(r8) c0mask(pcols)

   real(r8) hmin(pcols)
   real(r8) expdif(pcols)
   real(r8) expnum(pcols)
   real(r8) ftemp(pcols)
   real(r8) eps0(pcols)
   real(r8) rmue(pcols)
   real(r8) zuef(pcols)
   real(r8) zdef(pcols)
   real(r8) epsm(pcols)
   real(r8) ratmjb(pcols)
   real(r8) est(pcols)
   real(r8) totpcp(pcols)
   real(r8) totevp(pcols)
   real(r8) alfa(pcols)
   real(r8) ql1
   real(r8) tu
   real(r8) estu
   real(r8) qstu

   real(r8) small
   real(r8) mdt

!<songxl 2014-11-20---------
   real(r8) fice(pcols,pver)        ! ice fraction in precip production
   real(r8) tug(pcols,pver)
   real(r8) cmel(pcols,pver)        ! condensation rate of updraft
   real(r8) cmei(pcols,pver)        ! condensation rate of updraft

   real(r8) tvuo(pcols,pver)        ! updraft virtual T w/o freezing heating
   real(r8) tvu(pcols,pver)         ! updraft virtual T with freezing heating
   real(r8) totfrz(pcols)
   real(r8) frz1(pcols,pver)        ! rate of freezing
   real(r8) frz (pcols,pver)        ! rate of freezing
   integer  jto(pcols)              ! updraft plume old top
   integer  tmplel(pcols)

   integer  iter, itnum
   integer  m 
!>songxl 2014-11-20---------

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical doit(pcols)
   logical done(pcols)
!
!------------------------------------------------------------------------------
!

!<songxl 2014-11-20---------
      autolm(:il2g,:) = 0._r8
      accrlm(:il2g,:) = 0._r8
      bergnm(:il2g,:) = 0._r8
      fhtimm(:il2g,:) = 0._r8
      fhtctm(:il2g,:) = 0._r8
      fhmlm (:il2g,:) = 0._r8
      hmpim (:il2g,:) = 0._r8
      accslm(:il2g,:) = 0._r8
      dlfm  (:il2g,:) = 0._r8

      autoln(:il2g,:) = 0._r8
      accrln(:il2g,:) = 0._r8
      bergnn(:il2g,:) = 0._r8
      fhtimn(:il2g,:) = 0._r8
      fhtctn(:il2g,:) = 0._r8
      fhmln (:il2g,:) = 0._r8
      accsln(:il2g,:) = 0._r8
      activn(:il2g,:) = 0._r8
      dlfn  (:il2g,:) = 0._r8

      autoim(:il2g,:) = 0._r8
      accsim(:il2g,:) = 0._r8
      difm  (:il2g,:) = 0._r8

      nuclin(:il2g,:) = 0._r8
      autoin(:il2g,:) = 0._r8
      accsin(:il2g,:) = 0._r8
      hmpin (:il2g,:) = 0._r8
      difn  (:il2g,:) = 0._r8

      trspcm(:il2g,:) = 0._r8
      trspcn(:il2g,:) = 0._r8
      trspim(:il2g,:) = 0._r8
      trspin(:il2g,:) = 0._r8
!>songxl 2014-11-20---------

   do i = 1,il2g
      ftemp(i) = 0._r8
      expnum(i) = 0._r8
      expdif(i) = 0._r8
      c0mask(i)  = c0_ocn * (1._r8-landfrac(i)) +   c0_lnd * landfrac(i) 
      dcape(i) = 0._r8                   !songxl2014-11-20
   end do
!
!jr Change from msg+1 to 1 to prevent blowup
!
   do k = 1,pver
      do i = 1,il2g
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!
! initialize many output and work variables to zero
!
   pflx(:il2g,1) = 0

   do k = 1,pver
      do i = 1,il2g
         k1(i,k) = 0._r8
         i2(i,k) = 0._r8
         i3(i,k) = 0._r8
         i4(i,k) = 0._r8
         mu(i,k) = 0._r8
         f(i,k) = 0._r8
         eps(i,k) = 0._r8
         eu(i,k) = 0._r8
         du(i,k) = 0._r8
         ql(i,k) = 0._r8
         cu(i,k) = 0._r8
         evp(i,k) = 0._r8
         cmeg(i,k) = 0._r8
         qds(i,k) = q(i,k)
         md(i,k) = 0._r8
         ed(i,k) = 0._r8
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0._r8
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
         call qmmr_hPa(t(i,k), p(i,k), est(i), qst(i,k))
!++bee
         if ( p(i,k)-est(i) <= 0._r8 ) then
            qst(i,k) = 1.0_r8
         end if
!--bee
         gamma(i,k) = qst(i,k)*(1._r8 + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         rprd(i,k) = 0._r8
!<songxl 2014-11-20----------------
         sprd(i,k) = 0._r8
         wu(i,k)   = 0._r8
         fice(i,k) = 0._r8
         tug(i,k)  = 0._r8
         cmel(i,k) = 0._r8
         cmei(i,k) = 0._r8
         qc(i,k)   = 0._r8
         qi(i,k)   = 0._r8
         qnc(i,k)  = 0._r8
         qni(i,k)  = 0._r8
         qcde(i,k)   = 0._r8
         qide(i,k)   = 0._r8
         qncde(i,k)  = 0._r8
         qnide(i,k)  = 0._r8
         qnr(i,k)  = 0._r8
         qns(i,k)  = 0._r8
         qrain(i,k)= 0._r8
         qsnow(i,k)= 0._r8
         frz(i,k)  = 0._r8
         frz1(i,k) = 0._r8
         tvuo(i,k) = (shat(i,k) - grav/cp*zf(i,k))*(1._r8 + 0.608_r8*qhat(i,k))
         tvu(i,k) = tvuo(i,k)
!>songxl 2014-11-20----------------
      end do
   end do
!
!jr Set to zero things which make this routine blow up
!
   do k=1,msg
      do i=1,il2g
         rprd(i,k) = 0._r8
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do k = 1, msg+1
      do i = 1,il2g
         hsthat(i,k) = hsat(i,k)
         qsthat(i,k) = qst(i,k)
         gamhat(i,k) = gamma(i,k)
      end do
   end do
   do i = 1,il2g
      totpcp(i) = 0._r8
      totevp(i) = 0._r8
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6_r8) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_r8) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                                (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
      end do
   end do
!
! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
!
   jt(:) = pver
   jto(:) = pver              !songxl 2014-11-20
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.E6_r8
   end do
!
! find the level of minimum hsat, where detrainment starts
!

   do k = msg + 1,pver
      do i = 1,il2g
         if (hsat(i,k) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
            hmin(i) = hsat(i,k)
            j0(i) = k
         end if
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)
!
! Fix from Guang Zhang to address out of bounds array reference
!
      j0(i) = min(j0(i),pver)
   end do
!
! Initialize certain arrays inside cloud
!
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
            hu(i,k) = hmn(i,mx(i)) + cp*tiedke_add
            su(i,k) = s(i,mx(i)) + tiedke_add
         end if
      end do
   end do
!
! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
!
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5_r8* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5_r8* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5_r8* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
         end if
      end do
   end do
!
! re-initialize hmin array for ensuing calculation.
!
   do i = 1,il2g
      hmin(i) = 1.E6_r8
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(i,k) <= hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
         end if
      end do
   end do
!
! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************
!
   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0._r8
         ftemp(i) = 0._r8
         if (k < jt(i) .or. k >= jb(i)) then
            k1(i,k) = 0._r8
            expnum(i) = 0._r8
         else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) + &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
         end if
         if ((expdif(i) > 100._r8 .and. expnum(i) > 0._r8) .and. &
	     k1(i,k) > expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 + &
                     (2._r8*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2* &
                     ftemp(i)**3 + (-5._r8*k1(i,k)*i2(i,k)*i3(i,k)+ &
                     5._r8*i2(i,k)**3+k1(i,k)**2*i4(i,k))/ &
                     k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._r8)
            f(i,k) = min(f(i,k),0.0002_r8)
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(i,j0(i)) < 1.E-6_r8 .and. f(i,j0(i)+1) > f(i,j0(i))) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(i,j0(i))
      eps(i,jb(i)) = eps0(i)
   end do
!
! This is set to match the Rasch and Kristjansson paper
!
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
            eps(i,k) = f(i,j0(i))
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(i,k) = f(i,k)
      end do
   end do

!<songxl 2014-11-20---------------------------
   itnum = 1
   if (zm_microp) itnum = 2  

   do iter=1, itnum

      do k = pver,msg + 1,-1
        do i = 1,il2g
           cu(i,k) = 0._r8
           qc(i,k) = 0._r8
           qi(i,k) = 0._r8
           ql(i,k) = 0._r8
           frz1(i,k) = 0._r8
        end do
      end do
      do i = 1,il2g
          totpcp(i) = 0._r8
          hu(i,jb(i)) = hmn(i,jb(i)) + cp*tiedke_add
      end do

!
! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
!
   do i = 1,il2g
      if (eps0(i) > 0._r8) then
         mu(i,jb(i)) = 1._r8
         eu(i,jb(i)) = mu(i,jb(i))/dz(i,jb(i))
      end if
      if (zm_microp) then
        tmplel(i) = lel(i)
      else
        tmplel(i) = jt(i)
      end if
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
!<songxl 2014-11-20---------
!         if (eps0(i) > 0._r8 .and. (k >= jt(i) .and. k < jb(i))) then
!         if (eps0(i) > 0._r8 .and. (k >= lel(i) .and. k < jb(i))) then
          if (eps0(i) > 0._r8 .and. (k >= tmplel(i) .and. k < jb(i))) then
!>songxl 2014-11-20---------
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1._r8/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._r8)/zuef(i)
            mu(i,k) = (1._r8/eps0(i))* (exp(eps(i,k  )*zuef(i))-1._r8)/zuef(i)
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
         end if
      end do
   end do
!
   khighest = pverp
   klowest = 1
   do i=1,il2g
      khighest = min(khighest,lel(i))
      klowest = max(klowest,jb(i))
   end do
   do k = klowest-1,khighest,-1
      do i = 1,il2g
         if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0._r8) then
            if (mu(i,k) < 0.02_r8) then
               hu(i,k) = hmn(i,k)
               mu(i,k) = 0._r8
               eu(i,k) = 0._r8
               du(i,k) = mu(i,k+1)/dz(i,k)
            else
!<songxl 2014-11-20--------
           if (zm_microp) then
               hu(i,k) = (mu(i,k+1)*hu(i,k+1) + dz(i,k)*(eu(i,k)*hmn(i,k) +   &
                            latice*frz(i,k)))/(mu(i,k)+ dz(i,k)*du(i,k))
            else
               hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) + &
                         dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)- du(i,k)*hsat(i,k))
            end if
               

!>songxl 2014-11-20--------
            end if
         end if
      end do
   end do
!
! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
!
   do i=1,il2g
      doit(i) = .true.
!<songxl 2014-10-10-------------
      totfrz(i)= 0._r8
      do k = pver,msg + 1,-1
         totfrz(i)= totfrz(i)+ frz(i,k)*dz(i,k)
      end do
      if (iter == 1 .and. totfrz(i) > 0._r8) write(*,*) "totfrz=",totfrz(i)
!>songxl 2014-10-10--------------
   end do
   do k=klowest-2,khighest-1,-1
      do i=1,il2g
         if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
  	   if (hu(i,k) <= hsthat(i,k) .and. hu(i,k+1) > hsthat(i,k+1) &
	       .and. mu(i,k) >= 0.02_r8) then
               if (hu(i,k)-hsthat(i,k) < -2000._r8) then
                  jt(i) = k + 1
                  doit(i) = .false.
               else
                  jt(i) = k
                  doit(i) = .false.
               end if
!<songxl 2014-11-20----------
!            else if (hu(i,k) > hu(i,jb(i)) .or. mu(i,k) < 0.02_r8) then
            else if ( (hu(i,k) > hu(i,jb(i)) .and. totfrz(i)<=0._r8) .or. mu(i,k) < 0.02_r8) then
!>songxl 2014-11-20----------
               jt(i) = k + 1
               doit(i) = .false.
            end if
         end if
      end do
   end do

!<songxl 2014-11-20------
   do i = 1,il2g
      if (iter == 1)  jto(i) = jt(i)
   end do
!>songxl 2014-11-20------

   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0._r8) then
            mu(i,k) = 0._r8
            eu(i,k) = 0._r8
            du(i,k) = 0._r8
            hu(i,k) = hmn(i,k)
         end if
         if (k == jt(i) .and. eps0(i) > 0._r8) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            eu(i,k) = 0._r8
            mu(i,k) = 0._r8
         end if
      end do
   end do

!>songxl 2014-11-20------------

!
   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k == jb(i) .and. eps0(i) > 0._r8) then
            qu(i,k) = q(i,mx(i))
            su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
         end if
         if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0._r8) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                      dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                            du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cp*zf(i,k)
            call qmmr_hPa(tu, (p(i,k)+p(i,k-1))/2._r8, estu, qstu)
            if (qu(i,k) >= qstu) then
               jlcl(i) = k
               kount = kount + 1
               done(i) = .true.
            end if
         end if
      end do
      if (kount >= il2g) goto 690
   end do
690 continue
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1._r8+gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                     (rl* (1._r8+gamhat(i,k)))
         end if
      end do
   end do

! compute condensation in updraft

   do i = 1,il2g
     if (zm_microp) then
!     if(.false.) then
        tmplel(i) = jlcl(i)+1
     else
        tmplel(i) = jb(i)
     end if
   end do

! compute condensation in updraft
   do k = pver,msg + 2,-1
      do i = 1,il2g
!<songxl 2014-11-20------------------
!         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
          if (k >= jt(i) .and. k < tmplel(i) .and. eps0(i) > 0._r8) then
            if (zm_microp) then
               cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                      dz(i,k)- eu(i,k)*s(i,k)+du(i,k)*su(i,k))/(rl/cp)  &
                       - latice*frz(i,k)/rl
            else

               cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                      dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/(rl/cp)
            end if

!>songxl 2014-11-20--------------------

            if (k == jt(i)) cu(i,k) = 0._r8
            cu(i,k) = max(0._r8,cu(i,k))
         end if
      end do
   end do

!<songxl 2014-11-20----------

   if (zm_microp) then

      tug(:il2g,:) = t(:il2g,:)
      fice(:,:)    = 0._r8

      do k = pver, msg+2, -1
         do i = 1, il2g
            tug(i,k) = su(i,k) - grav/cp*zf(i,k)
         end do
      end do

      do k = 1, pver-1
         do i = 1, il2g

            if (tug(i,k+1) > 273.15_r8) then
               ! If warmer than tmax then water phase
               fice(i,k) = 0._r8

            else if (tug(i,k+1) < 233.15_r8) then
               ! If colder than tmin then ice phase
               fice(i,k) = 1._r8

            else
               ! Otherwise mixed phase, with ice fraction decreasing linearly
               ! from tmin to tmax
               fice(i,k) =(273.15_r8 - tug(i,k+1)) / 40._r8
            end if
         end do
      end do

      do k = 1, pver
         do i = 1,il2g
            cmei(i,k) = cu(i,k)* fice(i,k)
            cmel(i,k) = cu(i,k) * (1._r8-fice(i,k))
         end do
      end do

      call  zm_mphy(su,   qu,    mu,    du,   eu,    cmel,  cmei,   zf,    p,    t,    q,         &
                    eps0, jb,    jt,    jlcl, msg,   il2g,  grav,   cp,    rd,   aero, gamhat,    &
                    qc,   qi,    qnc,   qni,  qcde,  qide,  qncde,  qnide, rprd, sprd, frz,       &
                    wu,   qrain, qsnow, qnr,  qns,   autolm, accrlm, bergnm, fhtimm, fhtctm,      &
                    fhmlm,  hmpim,  accslm, dlfm,   autoln, accrln, bergnn, fhtimn, fhtctn,       &
                    fhmln,  accsln, activn, dlfn,   autoim, accsim, difm,   nuclin, autoin,       &
                    accsin, hmpin,  difn,   trspcm, trspcn, trspim, trspin, lambdadpcu, mudpcu  )


      do k = pver,msg + 2,-1
         do i = 1,il2g
            ql(i,k) = qc(i,k)+ qi(i,k)
!            if (iter == 1) frz1(i,k) = frz(i,k)         !songxl2014-11-20
            frz1(i,k) = frz(i,k)         !songxl2014-11-20 
         end do
      end do

      do i = 1,il2g
        if (iter == 2 .and. jt(i)> jto(i)) then
          do k = jt(i), jto(i), -1
             frz1(i,k) = 0.0_r8
             cu(i,k)=0.0_r8
          end do
        end if
      end do

      do k = pver,msg + 2,-1
         do i = 1,il2g
            if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
               totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*(qcde(i,k+1)+qide(i,k+1) ))
            end if
         end do
      end do

      do k = msg + 2,pver
        do i = 1,il2g
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
            if (iter == 1) tvuo(i,k)= (su(i,k) - grav/cp*zf(i,k))*(1._r8+0.608_r8*qu(i,k)  &
                                       -qc(i,k)-qi(i,k)-qrain(i,k)-qsnow(i,k) )
!            if (iter == 2 .and. k > max(jt(i),jto(i)) ) then
            if (iter == 2) then
               tvu(i,k) = (su(i,k) - grav/cp*zf(i,k))*(1._r8 +0.608_r8*qu(i,k)   &
                          -qc(i,k)-qi(i,k)-qrain(i,k)-qsnow(i,k))
               dcape(i) = dcape(i)+ rd*(tvu(i,k)-tvuo(i,k))*log(p(i,k)/p(i,k-1))
            end if
         end if
        end do
      end do

!>songxl 2014-10-10--------------

   else  ! no microphysics
!>songxl 2014-11-20----------

        ! compute condensed liquid, rain production rate
        ! accumulate total precipitation (condensation - detrainment of liquid)
        ! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
        ! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
        ! consistently applied.
        !    mu, ql are interface quantities
        !    cu, du, eu, rprd are midpoint quantites
      do k = pver,msg + 2,-1
         do i = 1,il2g
            rprd(i,k) = 0._r8
            if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(i,k) >= 0.0_r8) then
               if (mu(i,k) > 0._r8) then
                  ql1 = 1._r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                     dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
                  ql(i,k) = ql1/ (1._r8+dz(i,k)*c0mask(i))
               else
                  ql(i,k) = 0._r8
               end if
               totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*ql(i,k+1))
               rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
!<songxl 2014-11-20------
               qcde(i,k) = ql(i,k)
               qide(i,k) = 0._r8
               qncde(i,k) = 0._r8
               qnide(i,k) = 0._r8
               sprd(i,k) = 0._r8
!>songxl 2014-11-20------

             end if
         end do
      end do
!<songxl 2014-11-20-------
    end if  ! zm_microp

 end do   !iter   
!
! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
!
   do i = 1,il2g
!
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
!
      alfa(i) = 0.1_r8
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(i,jd(i)) = hmn(i,jd(i)-1)
      if (jd(i) < jb(i) .and. eps0(i) > 0._r8) then
         epsm(i) = eps0(i)
         md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
      end if
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i)
         end if
      end do
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)
         end if
      end do
   end do

   small = 1.e-20_r8
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0._r8) then
            ed(i,k-1) = (md(i,k-1)-md(i,k))/dz(i,k-1)
            mdt = min(md(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) - dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
         end if
      end do
   end do
!
! calculate updraft and downdraft properties.
!
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/ &
               (rl*(1._r8 + gamhat(i,k)))
         end if
      end do
   end do

!>songxl 2014-11-20------
!
   do i = 1,il2g
      qd(i,jd(i)) = qds(i,jd(i))
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            mdt = min(md(i,k+1),-small)
            if (zm_microp) then
              evp(i,k) = min(evp(i,k),rprd(i,k))
            end if
            sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) + md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   end do
   do i = 1,il2g
!*guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
      totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) - md(i,jb(i))*qd(i,jb(i))
   end do
!!$   if (.true.) then
   if (.false.) then
      do i = 1,il2g
         k = jb(i)
         if (eps0(i) > 0._r8) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._r8)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
         end if
      end do
   endif

   do i = 1,il2g
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do
!
   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0._r8 .and. totpcp(i) > 0._r8) then
            md(i,k)  = md (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(i,k)  = ed (i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(i,k) = 0._r8
            ed(i,k) = 0._r8
            evp(i,k) = 0._r8
         end if
! cmeg is the cloud water condensed - rain water evaporated
! rprd is the cloud water converted to rain - (rain evaporated)
         cmeg(i,k) = cu(i,k) - evp(i,k)
         rprd(i,k) = rprd(i,k)-evp(i,k)
      end do
   end do

! compute the net precipitation flux across interfaces
   pflx(:il2g,1) = 0._r8
   do k = 2,pverp
      do i = 1,il2g
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
      end do
   end do
!
   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
      end do
   end do
!
   return
end subroutine cldprp

subroutine closure(lchnk   , &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,capelmt )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
   use dycore,    only: dycore_is, get_resolution

   implicit none

!
!-----------------------------Arguments---------------------------------
!
   integer, intent(in) :: lchnk                 ! chunk identifier

   real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
   real(r8), intent(inout) :: t(pcols,pver)        ! temperature
   real(r8), intent(inout) :: p(pcols,pver)        ! pressure (mb)
   real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
   real(r8), intent(in) :: z(pcols,pver)        ! height (m)
   real(r8), intent(in) :: s(pcols,pver)        ! normalized dry static energy
   real(r8), intent(in) :: tp(pcols,pver)       ! parcel temp
   real(r8), intent(in) :: qs(pcols,pver)       ! sat spec humidity
   real(r8), intent(in) :: qu(pcols,pver)       ! updraft spec. humidity
   real(r8), intent(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   real(r8), intent(in) :: mc(pcols,pver)       ! net convective mass flux
   real(r8), intent(in) :: du(pcols,pver)       ! detrainment from updraft
   real(r8), intent(in) :: mu(pcols,pver)       ! mass flux of updraft
   real(r8), intent(in) :: md(pcols,pver)       ! mass flux of downdraft
   real(r8), intent(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
   real(r8), intent(in) :: sd(pcols,pver)       ! dry static energy of downdraft
   real(r8), intent(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
   real(r8), intent(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
   real(r8), intent(in) :: dp(pcols,pver)       ! pressure thickness of layers
   real(r8), intent(in) :: qstp(pcols,pver)     ! spec humidity of parcel
   real(r8), intent(in) :: zf(pcols,pver+1)     ! height of interface levels
   real(r8), intent(in) :: ql(pcols,pver)       ! liquid water mixing ratio
   real(r8), intent(in) :: cape(pcols)          ! available pot. energy of column
   real(r8), intent(in) :: tl(pcols)
   real(r8), intent(in) :: dsubcld(pcols)       ! thickness of subcloud layer

   integer, intent(in) :: lcl(pcols)        ! index of lcl
   integer, intent(in) :: lel(pcols)        ! index of launch leve
   integer, intent(in) :: jt(pcols)         ! top of updraft
   integer, intent(in) :: mx(pcols)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real(r8) dtpdt(pcols,pver)
   real(r8) dqsdtp(pcols,pver)
   real(r8) dtmdt(pcols,pver)
   real(r8) dqmdt(pcols,pver)
   real(r8) dboydt(pcols,pver)
   real(r8) thetavp(pcols,pver)
   real(r8) thetavm(pcols,pver)
   real(r8) dtbdt(pcols),dqbdt(pcols),dtldt(pcols)
   real(r8) beta
   real(r8) capelmt
   real(r8) cp
   real(r8) dadt(pcols)
   real(r8) debdt
   real(r8) dltaa
   real(r8) eb
   real(r8) grav

   integer i
   integer il1g
   integer il2g
   integer k, kmin, kmax
   integer msg


   real(r8) rd
   real(r8) rl
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
      mb(i) = 0._r8
      eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      dtbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
      dqbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
                 md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
      debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
      dtldt(i) = -2840._r8* (3.5_r8/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t(i,mx(i)))-log(eb)-4.805_r8)**2
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(i,k) = 0._r8
         dqmdt(i,k) = 0._r8
      end do
   end do
!
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
            dqmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                         qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
         end if
      end do
   end do
!
   beta = 0._r8
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                         dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
            dqsdtp(i,k) = qstp(i,k)* (1._r8+qstp(i,k)/eps1)*eps1*rl/(rd*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(i,k) = tp(i,k)/ (1._r8+rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                        (dtbdt(i)/t(i,mx(i))+rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                         tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._r8/(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))* &
                          (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_r8/ &
                          (1._r8+0.608_r8*q(i,k))*dqmdt(i,k)))*grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+0.608_r8/ (1._r8+0.608_r8*q(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t(i,k)-0.608_r8/ (1._r8+0.608_r8*q(i,k))*dqmdt(i,k))* &
                          grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   dadt(il1g:il2g)  = 0._r8
   kmin = minval(lel(il1g:il2g))
   kmax = maxval(mx(il1g:il2g)) - 1
   do k = kmin, kmax
      do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         endif
      end do
   end do

   do i = il1g,il2g
      dltaa = -1._r8* (cape(i)-capelmt)
      if (dadt(i) /= 0._r8) mb(i) = max(dltaa/tau/dadt(i),0._r8)
   end do
!
   return
end subroutine closure

subroutine q1q2_pjr(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat    ,shat    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
!<songxl 2014-11-20----------
!                    dl      ,evp     ,cu      )
                    dl      ,evp     ,cu      ,qice    ,dice    , &
                    qnl     ,qni     ,dnl     ,dni     ,frz  )
!>songxl 2014-11-20----------

   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
! 
!-----------------------------------------------------------------------


   real(r8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   real(r8), intent(in) :: q(pcols,pver)
   real(r8), intent(in) :: qs(pcols,pver)
   real(r8), intent(in) :: qu(pcols,pver)
   real(r8), intent(in) :: su(pcols,pver)
   real(r8), intent(in) :: du(pcols,pver)
   real(r8), intent(in) :: qhat(pcols,pver)
   real(r8), intent(in) :: shat(pcols,pver)
   real(r8), intent(in) :: dp(pcols,pver)
   real(r8), intent(in) :: mu(pcols,pver)
   real(r8), intent(in) :: md(pcols,pver)
   real(r8), intent(in) :: sd(pcols,pver)
   real(r8), intent(in) :: qd(pcols,pver)
   real(r8), intent(in) :: ql(pcols,pver)
   real(r8), intent(in) :: evp(pcols,pver)
   real(r8), intent(in) :: cu(pcols,pver)
   real(r8), intent(in) :: dsubcld(pcols)
!<songxl 2014-11-20------------
   real(r8), intent(in) :: frz(pcols,pver)
   real(r8), intent(in) :: qice(pcols,pver)
   real(r8), intent(in) :: qnl(pcols,pver)
   real(r8), intent(in) :: qni(pcols,pver)
!>songxl 2014-11-20------------

   real(r8),intent(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
   real(r8),intent(out) :: dl(pcols,pver)
!<songxl 2014-11-20------------
   real(r8),intent(out) :: dice(pcols,pver)
   real(r8),intent(out) :: dnl(pcols,pver)
   real(r8),intent(out) :: dni(pcols,pver)  
!>songxl 2014-11-20------------

   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
!
! work fields:
!
   integer i
   integer k

   real(r8) emc
   real(r8) rl
!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0._r8
         dqdt(i,k) = 0._r8
         dl(i,k) = 0._r8
!<songxl 2014-11-20------
         dice(i,k) = 0._r8
         dnl(i,k)  = 0._r8
         dni(i,k)  = 0._r8
!>songxl 2014-11-20------
      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

!>songxl 2014-11-20-------------
         if (zm_microp) dsdt(i,k) = dsdt(i,k) + latice/cp*frz(i,k) 
!<songxl 2014-11-20--------------

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = du(i,k)*ql(i,k+1)
!<songxl 2014-11-20------------
         if (zm_microp) then
           dice(i,k) = du(i,k)*qice(i,k+1)
           dnl(i,k)  = du(i,k)*qnl(i,k+1)
           dni(i,k)  = du(i,k)*qni(i,k+1)
         end if
!>songxl 2014-11-20------------

      end do
   end do

!
!DIR$ NOINTERCHANGE!
   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        )
            dqdt(i,k) = (1._r8/dsubcld(i))* &
                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        )
         else if (k > mx(i)) then
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
         end if
      end do
   end do
!
   return
end subroutine q1q2_pjr


subroutine buoyan_dilute(lchnk   ,ncol    , &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,lon     ,mx      , &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)         ! perturbation temperature by pbl processes

!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)       ! parcel temperature
   real(r8), intent(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(r8), intent(out) :: tl(pcols)            ! parcel temperature at lcl
   real(r8), intent(out) :: cape(pcols)          ! convective aval. pot. energy.
   integer lcl(pcols)        !
   integer lel(pcols)        !
   integer lon(pcols)        ! level of onset of deep convection
   integer mx(pcols)         ! level of max moist static energy
!
!--------------------------Local Variables------------------------------
!
   real(r8) capeten(pcols,5)     ! provisional value of cape
   real(r8) tv(pcols,pver)       !
   real(r8) tpv(pcols,pver)      !
   real(r8) buoy(pcols,pver)

   real(r8) a1(pcols)
   real(r8) a2(pcols)
   real(r8) estp(pcols)
   real(r8) pl(pcols)
   real(r8) plexp(pcols)
   real(r8) hmax(pcols)
   real(r8) hmn(pcols)
   real(r8) y(pcols)

   logical plge600(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)

   real(r8) cp
   real(r8) e
   real(r8) grav

   integer i
   integer k
   integer msg
   integer n

   real(r8) rd
   real(r8) rl
#ifdef PERGRO
   real(r8) rhd
#endif
!
!-----------------------------------------------------------------------
!
   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = pver
      knt(i) = 0
      lel(i) = pver
      mx(i) = lon(i)
      cape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:ncol,:) = t(:ncol,:) *(1._r8+1.608_r8*q(:ncol,:))/ (1._r8+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._r8

!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
#ifdef PERGRO
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
!
! Reset max moist static energy level when relative difference exceeds 1.e-4
!
         rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#else
   do k = pver,msg + 1,-1
      do i = 1,ncol
         hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
#endif


! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

   do i = 1,ncol ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl(i) = t(i,mx(i))
      pl(i) = p(i,mx(i))
   end do

!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!

   call parcel_dilute(lchnk, ncol, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)


! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
   end do

!
! Main buoyancy calculation.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            buoy(i,k) = tpv(i,k) - tv(i,k) + tiedke_add  ! +0.5K or not?
         else
            qstp(i,k) = q(i,k)
            tp(i,k)   = t(i,k)            
            tpv(i,k)  = tv(i,k)
         endif
      end do
   end do



!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0._r8 .and. buoy(i,k) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do

!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
!
   return
end subroutine buoyan_dilute

subroutine parcel_dilute (lchnk, ncol, msg, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: lchnk
integer, intent(in) :: ncol
integer, intent(in) :: msg

integer, intent(in), dimension(pcols) :: klaunch(pcols)

real(r8), intent(in), dimension(pcols,pver) :: p
real(r8), intent(in), dimension(pcols,pver) :: t
real(r8), intent(in), dimension(pcols,pver) :: q
real(r8), intent(in), dimension(pcols) :: tpert ! PBL temperature perturbation.

real(r8), intent(inout), dimension(pcols,pver) :: tp    ! Parcel temp.
real(r8), intent(inout), dimension(pcols,pver) :: qstp  ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(pcols) :: tl         ! Actual temp of LCL.
real(r8), intent(inout), dimension(pcols) :: pl          ! Actual pressure of LCL. 

integer, intent(inout), dimension(pcols) :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(pcols,pver) :: tpv   ! Define tpv within this routine.

!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(pcols,pver)        ! Tempertaure of the entraining parcel.
real(r8) qtmix(pcols,pver)       ! Total water of the entraining parcel.
real(r8) qsmix(pcols,pver)       ! Saturated mixing ratio at the tmix.
real(r8) smix(pcols,pver)        ! Entropy of the entraining parcel.
real(r8) xsh2o(pcols,pver)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(pcols,pver)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(pcols,pver)   ! Entropy change sue to freezing of precip.

real(r8) mp(pcols)    ! Parcel mass flux.
real(r8) qtp(pcols)   ! Parcel total water.
real(r8) sp(pcols)    ! Parcel entropy.

real(r8) sp0(pcols)    ! Parcel launch entropy.
real(r8) qtp0(pcols)   ! Parcel launch total water.
real(r8) mp0(pcols)    ! Parcel launch relative mass flux.

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat     ! Number of iterations for condensation/freezing loop.
integer i,k,ii   ! Loop counters.

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
!dmpdz=-1.e-3_r8        ! Entrainment rate. (-ve for /m)
!++wy tuning
!lys dmpdz=-1.5e-3_r8        ! Entrainment rate. (-ve for /m)
dmpdz=-2.0e-3_r8        ! Entrainment rate. (-ve for /m)
!--wy tuning
!dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
lwmax = 1.e-3_r8    ! Need to put formula in for this.
tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

! **** Begin loops ****

do k = pver, msg+1, -1
   do i=1,ncol 

! Initialize parcel values at launch level.

      if (k == klaunch(i)) then 
         qtp0(i) = q(i,k)   ! Parcel launch total water (assuming subsaturated) - OK????.
         sp0(i)  = entropy(t(i,k),p(i,k),qtp0(i))  ! Parcel launch entropy.
         mp0(i)  = 1._r8       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
         smix(i,k)  = sp0(i)
         qtmix(i,k) = qtp0(i)
         tfguess = t(i,k)
         rcall = 1
         call ientropy (rcall,i,lchnk,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess)
      end if

! Entraining levels
      
      if (k < klaunch(i)) then 

! Set environmental values for this level.                 
         
         dp = (p(i,k)-p(i,k+1)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5_r8*(q(i,k)+q(i,k+1))         ! Total water of environment.
         tenv  = 0.5_r8*(t(i,k)+t(i,k+1)) 
         penv  = 0.5_r8*(p(i,k)+p(i,k+1))

         senv  = entropy(tenv,penv,qtenv)  ! Entropy of environment.   

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1._r8/dpdz                  ! in m/mb
         dmpdp = dmpdz*dzdp              ! /mb Fractional entrainment

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv 
         qtp(i) = qtp(i) - dmpdp*dp*qtenv 
         mp(i)  = mp(i)  - dmpdp*dp
            
! Entrain s and qt to next level.

         smix(i,k)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(i,k) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

         tfguess = tmix(i,k+1)
         rcall = 2
         call ientropy(rcall,i,lchnk,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess)   

!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.

         if (qsmix(i,k) <= qtmix(i,k) .and. qsmix(i,k+1) > qtmix(i,k+1)) then
            lcl(i) = k
            qxsk   = qtmix(i,k) - qsmix(i,k)
            qxskp1 = qtmix(i,k+1) - qsmix(i,k+1)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(i,k+1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp   = (smix(i,k)  - smix(i,k+1))/dp
            dqtdp  = (qtmix(i,k) - qtmix(i,k+1))/dp
            slcl   = smix(i,k+1)  +  dsdp* (pl(i)-p(i,k+1))  
            qtlcl  = qtmix(i,k+1) +  dqtdp*(pl(i)-p(i,k+1))

            tfguess = tmix(i,k)
            rcall = 3
            call ientropy (rcall,i,lchnk,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)

!            write(iulog,*)' '
!            write(iulog,*)' p',p(i,k+1),pl(i),p(i,lcl(i))
!            write(iulog,*)' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
!            write(iulog,*)' s',smix(i,k+1),slcl,smix(i,lcl(i))
!            write(iulog,*)'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
!            write(iulog,*)'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

         endif
!         
      end if !  k < klaunch

 
   end do ! Levels loop
end do ! Columns loop

!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

xsh2o = 0._r8
ds_xsh2o = 0._r8
ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do k = pver, msg+1, -1
   do i=1,ncol    
      
! Initialize variables at k=klaunch
      
      if (k == klaunch(i)) then

! Set parcel values at launch level assume no liquid water.            

         tp(i,k)    = tmix(i,k)
         qstp(i,k)  = q(i,k) 
         tpv(i,k)   =  (tp(i,k) + tpert(i)) * (1._r8+1.608_r8*qstp(i,k)) / (1._r8+qstp(i,k))
         
      end if

      if (k < klaunch(i)) then
            
! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(i,k) = max (0._r8, qtmix(i,k) - qsmix(i,k) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     
                     
            ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - cpliq * log (tmix(i,k)/tfreez) * max(0._r8,(xsh2o(i,k)-xsh2o(i,k+1)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
!
 
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) == 0._r8) then ! One off freezing of condensate. 
               ds_freeze(i,k) = (latice/tmix(i,k)) * max(0._r8,qtmix(i,k)-qsmix(i,k)-xsh2o(i,k)) ! Gain of LH
            end if
            
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) /= 0._r8) then ! Continual freezing of additional condensate.
               ds_freeze(i,k) = ds_freeze(i,k+1)+(latice/tmix(i,k)) * max(0._r8,(qsmix(i,k+1)-qsmix(i,k)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k) 

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(i,k) - xsh2o(i,k)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(i,k)
            rcall =4
            call ientropy (rcall,i,lchnk,new_s, p(i,k), new_q, tmix(i,k), qsmix(i,k), tfguess)
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(i,k)    = tmix(i,k)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(i,k)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(i,k) = qsmix(i,k)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(i,k) = new_q
         end if

         tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+1.608_r8*qstp(i,k)) / (1._r8+ new_q) 

      end if ! k < klaunch
      
   end do ! Loop for columns
   
end do  ! Loop for vertical levels.


return
end subroutine parcel_dilute

!-----------------------------------------------------------------------------------------
real(r8) function entropy(TK,p,qtot)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(r8), intent(in) :: p,qtot,TK
     real(r8) :: qv,qst,e,est,L,eref,pref

pref = 1000.0_r8           ! mb
eref = 6.106_r8            ! sat p at tfreez (mb)

L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE

! Replace call to satmixutils.

call qmmr_hPa(TK, p, est, qst)

qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
e = qv*p / (eps1 +qv)

entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qst)
! 
return
end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
   SUBROUTINE ientropy (rcall,icol,lchnk,s,p,qt,T,qst,Tfg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

     use phys_grid, only: get_rlon_p, get_rlat_p

     integer, intent(in) :: icol, lchnk, rcall
     real(r8), intent(in)  :: s, p, Tfg, qt
     real(r8), intent(out) :: qst, T
     real(r8) :: qv,Ts,dTs,fs1,fs2,est
     real(r8) :: pref,eref,L,e
     real(r8) :: this_lat,this_lon
     integer :: LOOPMAX,i

LOOPMAX = 100                   !* max number of iteration loops 

! Values for entropy
pref = 1000.0_r8           ! mb ref pressure.
eref = 6.106_r8           ! sat p at tfreez (mb)

! Invert the entropy equation -- use Newton's method

Ts = Tfg                  ! Better first guess based on Tprofile from conv.

converge: do i=0, LOOPMAX

   L = rl - (cpliq - cpwv)*(Ts-tfreez) 

   call qmmr_hPa(Ts, p, est, qst)
   qv = min(qt,qst) 
   e = qv*p / (eps1 +qv)  ! Bolton (eq. 16)
   fs1 = (cpres + qt*cpliq)*log( Ts/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/Ts - qv*rh2o*log(qv/qst) - s
   
   L = rl - (cpliq - cpwv)*(Ts-1._r8-tfreez)         

   call qmmr_hPa(Ts-1._r8, p, est, qst)
   qv = min(qt,qst) 
   e = qv*p / (eps1 +qv)
   fs2 = (cpres + qt*cpliq)*log( (Ts-1._r8)/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/(Ts-1._r8) - qv*rh2o*log(qv/qst) - s 
   
   dTs = fs1/(fs2 - fs1)
   Ts  = Ts+dTs
   if (abs(dTs).lt.0.001_r8) exit converge
   if (i .eq. LOOPMAX - 1) then
      this_lat = get_rlat_p(lchnk, icol)*57.296_r8
      this_lon = get_rlon_p(lchnk, icol)*57.296_r8
      write(iulog,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(iulog,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
       ' lat: ',this_lat,' lon: ',this_lon, &
       ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
       ' qst(g/kg) = ', 1000._r8*qst,', s(J/kg) = ',s
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
   end if
enddo converge

! Replace call to satmixutils.

call qmmr_hPa(Ts, p, est, qst)

qv = min(qt,qst)                             !       /* check for saturation */
T = Ts 

 100    format (A,I1,I4,I4,7(A,F6.2))

return
end SUBROUTINE ientropy

! Wrapper for qmmr that does translation between Pa and hPa
! qmmr uses Pa internally, so get qmmr right, need to pass in Pa.
! Afterward, set es back to hPa.
elemental subroutine qmmr_hPa(t, p, es, qm)
  use wv_saturation, only: qmmr

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature (K)
  real(r8), intent(in) :: p    ! Pressure (hPa)
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure (hPa)
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

  call qmmr(t, p*100._r8, es, qm)

  es = es*0.01_r8

end subroutine qmmr_hPa

end module zm_conv
