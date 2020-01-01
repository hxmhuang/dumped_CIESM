
  module cldwat2m_macro

  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !--------------------------------------------------- !

   use shr_kind_mod,     only: r8=>shr_kind_r8
   use spmd_utils,       only: masterproc
   use ppgrid,           only: pcols, pver, pverp
   use abortutils,       only: endrun
   use physconst,        only: cpair, latvap, latice, rh2o
   use wv_saturation,    only: qsat_water, svp_water, svp_ice
   use cam_history,      only: addfld, add_default, phys_decomp, outfld 
   use cam_logfile,      only: iulog
   use ref_pres,         only: top_lev=>trop_cloud_top_lev

! qinyi
  use cldfrac_conden_mod, only: subgrid_var,cldfrac_conden_single,cldfrac_conden, cldfrac_conden_single_only


   implicit none
   private
   save

   public :: ini_macro, mmacro_pcond, aist_vector, aist_single

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! -------------------------- !
   ! Parameters for Ice Stratus !
   ! -------------------------- !
   real(r8), public,  parameter :: rhmini       = 0.80_r8       ! Minimum rh for ice cloud fraction > 0.
   real(r8), public,  parameter :: rhmaxi       = 1.1_r8        ! rhi at which ice cloud fraction = 1.
   real(r8), parameter          :: qist_min     = 1.e-7_r8      ! Minimum in-stratus ice IWC constraint [ kg/kg ]
   real(r8), parameter          :: qist_max     = 5.e-3_r8      ! Maximum in-stratus ice IWC constraint [ kg/kg ]

   ! ----------------------------- !
   ! Parameters for Liquid Stratus !
   ! ----------------------------- !

                                                                ! use Slingo (triangular PDF-based) liquid stratus fraction
   logical,  parameter          :: freeze_dry   = .false.       ! If .true., use 'freeze dry' in liquid stratus fraction formula
   real(r8), parameter          :: qlst_min     = 2.e-5_r8      ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: qlst_max     = 3.e-3_r8      ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter          :: cc           = 0.1_r8        ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )

   real(r8), parameter          :: ramda        = 0.5_r8        ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(r8), private            :: rhminl                       ! Critical RH for low-level  liquid stratus clouds
   real(r8), private            :: rhminl_adj_land              ! rhminl adjustment for snowfree land
   real(r8), private            :: rhminh                       ! Critical RH for high-level liquid stratus clouds
   real(r8), private            :: premit                       ! Top    height for mid-level liquid stratus fraction
   real(r8), private            :: premib                       ! Bottom height for mid-level liquid stratus fraction
   integer,  private            :: iceopt                       ! option for ice cloud closure 
                                                                ! 1=wang & sassen 2=schiller (iciwc)  
                                                                ! 3=wood & field, 4=Wilson (based on smith)
                                                                ! 5=modified slingo (ssat & empyt cloud)        
   real(r8), private            :: icecrit                      ! Critical RH for ice clouds in Wilson & Ballard closure
                                                                ! ( smaller = more ice clouds )

   contains

   ! -------------- !
   ! Initialization !
   ! -------------- !

   subroutine ini_macro

   !--------------------------------------------------------------------- !
   !                                                                      ! 
   ! Purpose: Initialize constants for the liquid stratiform macrophysics !
   !          called from stratiform.F90.                                 !  
   ! Author:  Sungsu Park, Dec.01.2009.                                   !
   !                                                                      !
   !--------------------------------------------------------------------- !

   use cloud_fraction, only: cldfrc_getparams

   call cldfrc_getparams(rhminl_out = rhminl, rhminl_adj_land_out = rhminl_adj_land,    &
                         rhminh_out = rhminh, premit_out = premit,   premib_out = premib, &
                         iceopt_out = iceopt, icecrit_out = icecrit)

   if( masterproc ) then
       write(iulog,*) 'ini_macro: tuning parameters : rhminl = ', rhminl, &
                                                     'rhminl_adj_land = ', rhminl_adj_land, & 
                                                     'rhminh = ', rhminh, & 
                                                     'premit = ', premit, & 
                                                     'premib = ',  premib,  &
                                                     'iceopt = ',  iceopt,  &
                                                     'icecrit = ', icecrit
   end if

   return
   end subroutine ini_macro

   ! ------------------------------ !
   ! Stratiform Liquid Macrophysics !
   ! ------------------------------ !

   ! In the version, 'macro --> micro --> advective forcing --> macro...'
   ! So, 'A' and 'C' are exclusive. 

   subroutine mmacro_pcond( lchnk      , ncol       , dt         , p            , dp         ,              &
                            T0         , qv0        , ql0        , qi0          , nl0        , ni0        , &
                            a_cud      , a_cu0      , landfrac   , snowh        ,                           & 
                            s_tendout  , qv_tendout , ql_tendout , qi_tendout   , nl_tendout , ni_tendout , &
                            qme        , qvadj      , qladj      , qiadj        , qllim      , qilim      , &
                            cld        , al_st_star , ai_st_star , ql_st_star   , qi_st_star , do_cldice  , &
                            rpdel, zm,  &
                            tke, kvh, kvm, lengi, shi, smi, wstarPBL, pblh, &
                            qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                            alst_gp, cond_gp, alst_def, N1, &
                            sgm1, sgm, sgm_dprime,&
                            delta_q, deltaq_sat, deltaq_uns, Q1_sat, Q1_uns,&
                            adjust_factor, & ! qinyi 2017-1-1 22:03:09
                            explwi2qtu) ! qinyi 2017-11-30 21:47:10

   use constituents,     only : qmin, cnst_get_ind
   use time_manager,     only : is_first_step, get_nstep
   use wv_saturation,    only : findsp_vc

   implicit none

   integer   icol
   integer,  intent(in)    :: lchnk                        ! Chunk number
   integer,  intent(in)    :: ncol                         ! Number of active columns

   ! qinyi
   real(r8), intent(in) :: rpdel(pcols,pver)
   real(r8), intent(in) :: zm(pcols,pver)

   real(r8), intent(in) :: tke(pcols,pverp)
   real(r8), intent(in) :: kvh(pcols,pverp)
   real(r8), intent(in) :: kvm(pcols,pverp)
   real(r8), intent(in) :: lengi(pcols,pverp)
   real(r8), intent(in) :: shi(pcols,pverp)
   real(r8), intent(in) :: smi(pcols,pverp)
   real(r8), intent(in) :: wstarPBL(pcols)
   real(r8), intent(in) :: pblh(pcols)

   ! shallow
   real(r8), intent(in) :: qtu_shal(pcols,pverp)
   real(r8), intent(in) :: umf_shal(pcols,pverp)
   real(r8), intent(in) :: cnt_shal(pcols)
   real(r8), intent(in) :: cnb_shal(pcols)
   real(r8), intent(in) :: thlu_shal(pcols,pverp)
   ! qinyi 2017-11-30 21:47:27
   real(r8), intent(in) :: explwi2qtu(pcols,pver)
   ! qinyi 

   ! Input-Output variables

   real(r8), intent(inout) :: T0(pcols,pver)               ! Temperature [K]
   real(r8), intent(inout) :: qv0(pcols,pver)              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(inout) :: ql0(pcols,pver)              ! Grid-mean liquid water content [kg/kg]
   real(r8), intent(inout) :: qi0(pcols,pver)              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(inout) :: nl0(pcols,pver)              ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(inout) :: ni0(pcols,pver)              ! Grid-mean number concentration of cloud ice    droplet [#/kg]

   ! Input variables

   real(r8), intent(in)    :: dt                           ! Model integration time step [s]
   real(r8), intent(in)    :: p(pcols,pver)                ! Pressure at the layer mid-point [Pa]
   real(r8), intent(in)    :: dp(pcols,pver)               ! Pressure thickness [Pa] > 0
                                                          ! within liquid stratus [kg/kg/s]
   real(r8), intent(in)    :: a_cud(pcols,pver)            ! Old cumulus fraction before update
   real(r8), intent(in)    :: a_cu0(pcols,pver)            ! New cumulus fraction after update

   real(r8), intent(in)    :: landfrac(pcols)              ! Land fraction
   real(r8), intent(in)    :: snowh(pcols)                 ! Snow depth (liquid water equivalent)
   logical, intent(in)     :: do_cldice                    ! Whether or not cldice should be prognosed

   ! Output variables

   real(r8), intent(out)   :: s_tendout(pcols,pver)        ! Net tendency of grid-mean s  from 'Micro+Macro' processes [J/kg/s]
   real(r8), intent(out)   :: qv_tendout(pcols,pver)       ! Net tendency of grid-mean qv from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: ql_tendout(pcols,pver)       ! Net tendency of grid-mean ql from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: qi_tendout(pcols,pver)       ! Net tendency of grid-mean qi from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: nl_tendout(pcols,pver)       ! Net tendency of grid-mean nl from 'Micro+Macro' processes [#/kg/s]
   real(r8), intent(out)   :: ni_tendout(pcols,pver)       ! Net tendency of grid-mean ni from 'Micro+Macro' processes [#/kg/s]

   real(r8), intent(out)   :: qme  (pcols,pver)            ! Net condensation rate [kg/kg/s]
   real(r8), intent(out)   :: qvadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (vapor)
   real(r8), intent(out)   :: qladj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (liquid)
   real(r8), intent(out)   :: qiadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (ice)
   real(r8), intent(out)   :: qllim(pcols,pver)            ! tendency from "instratus_condensate" call (liquid)
   real(r8), intent(out)   :: qilim(pcols,pver)            ! tendency from "instratus_condensate" call (ice)

   real(r8), intent(out)   :: cld(pcols,pver)              ! Net cloud fraction ( 0 <= cld <= 1 )
   real(r8), intent(out)   :: al_st_star(pcols,pver)       ! Physical liquid stratus fraction
   real(r8), intent(out)   :: ai_st_star(pcols,pver)       ! Physical ice stratus fraction
   real(r8), intent(out)   :: ql_st_star(pcols,pver)       ! In-stratus LWC [kg/kg] 
   real(r8), intent(out)   :: qi_st_star(pcols,pver)       ! In-stratus IWC [kg/kg] 

! qinyi
   real(r8), intent(out) :: alst_gp(pcols,pver)            ! cloud fraction from Gauss-PDF
   real(r8), intent(out) :: cond_gp(pcols,pver)            ! cloud condensate from Gauss-PDF
   real(r8), intent(out) :: alst_def(pcols,pver)           ! cloud fraction calculated by the default cloud scheme
   real(r8), intent(out) :: N1(pcols,pver)                 ! use to count the occurrence of over-saturated situation

   real(r8), intent(out) :: sgm1(pcols,pver)
   real(r8), intent(out) :: sgm(pcols,pver)
   real(r8), intent(out) :: sgm_dprime(pcols,pver)

   real(r8), intent(out) :: delta_q(pcols,pver)
   real(r8), intent(out) :: deltaq_sat(pcols,pver)
   real(r8), intent(out) :: deltaq_uns(pcols,pver)
   real(r8), intent(out) :: Q1_sat(pcols,pver)
   real(r8), intent(out) :: Q1_uns(pcols,pver)

   real(r8), intent(out) :: adjust_factor(pcols,pver)

   ! --------------- !
   ! Local variables !
   ! --------------- !
   integer :: ixcldliq, ixcldice
 
   integer :: i, j, k, ii, jj                              ! Loop indexes

   ! Thermodynamic state variables

   real(r8) T(pcols,pver)                                  ! Temperature of equilibrium reference state
                                                           ! from which 'Micro & Macro' are computed [K]
   real(r8) T1(pcols,pver)                                 ! Temperature after 'fice_force' on T01  
   real(r8) T_0(pcols,pver)                                ! Temperature after 'instratus_condensate' on T1
   real(r8) T_05(pcols,pver)                               ! Temperature after 'advection' on T_0 
   real(r8) T_prime0(pcols,pver)                           ! Temperature after 'Macrophysics (QQ)' on T_05star
   real(r8) T_dprime(pcols,pver)                           ! Temperature after 'fice_force' on T_prime
   real(r8) T_star(pcols,pver)                             ! Temperature after 'instratus_condensate' on T_dprime

   real(r8) qv(pcols,pver)                                 ! Grid-mean qv of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qv1(pcols,pver)                                ! Grid-mean qv after 'fice_force' on qv01  
   real(r8) qv_0(pcols,pver)                               ! Grid-mean qv after 'instratus_condensate' on qv1
   real(r8) qv_05(pcols,pver)                              ! Grid-mean qv after 'advection' on qv_0 
   real(r8) qv_prime0(pcols,pver)                          ! Grid-mean qv after 'Macrophysics (QQ)' on qv_05star
   real(r8) qv_dprime(pcols,pver)                          ! Grid-mean qv after 'fice_force' on qv_prime
   real(r8) qv_star(pcols,pver)                            ! Grid-mean qv after 'instratus_condensate' on qv_dprime

   real(r8) ql(pcols,pver)                                 ! Grid-mean ql of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ql1(pcols,pver)                                ! Grid-mean ql after 'fice_force' on ql01  
   real(r8) ql_0(pcols,pver)                               ! Grid-mean ql after 'instratus_condensate' on ql1
   real(r8) ql_05(pcols,pver)                              ! Grid-mean ql after 'advection' on ql_0 
   real(r8) ql_prime0(pcols,pver)                          ! Grid-mean ql after 'Macrophysics (QQ)' on ql_05star
   real(r8) ql_dprime(pcols,pver)                          ! Grid-mean ql after 'fice_force' on ql_prime
   real(r8) ql_star(pcols,pver)                            ! Grid-mean ql after 'instratus_condensate' on ql_dprime

   real(r8) qi(pcols,pver)                                 ! Grid-mean qi of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qi1(pcols,pver)                                ! Grid-mean qi after 'fice_force' on qi01  
   real(r8) qi_0(pcols,pver)                               ! Grid-mean qi after 'instratus_condensate' on qi1
   real(r8) qi_05(pcols,pver)                              ! Grid-mean qi after 'advection' on qi_0 
   real(r8) qi_prime0(pcols,pver)                          ! Grid-mean qi after 'Macrophysics (QQ)' on qi_05star
   real(r8) qi_dprime(pcols,pver)                          ! Grid-mean qi after 'fice_force' on qi_prime
   real(r8) qi_star(pcols,pver)                            ! Grid-mean qi after 'instratus_condensate' on qi_dprime

   real(r8) nl(pcols,pver)                                 ! Grid-mean nl of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) nl1(pcols,pver)                                ! Grid-mean nl after 'fice_force' on nl01  
   real(r8) nl_0(pcols,pver)                               ! Grid-mean nl after 'instratus_condensate' on nl1
   real(r8) nl_05(pcols,pver)                              ! Grid-mean nl after 'advection' on nl_0 
   real(r8) nl_prime0(pcols,pver)                          ! Grid-mean nl after 'Macrophysics (QQ)' on nl_05star
   real(r8) nl_dprime(pcols,pver)                          ! Grid-mean nl after 'fice_force' on nl_prime
   real(r8) nl_star(pcols,pver)                            ! Grid-mean nl after 'instratus_condensate' on nl_dprime

   real(r8) ni(pcols,pver)                                 ! Grid-mean ni of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ni1(pcols,pver)                                ! Grid-mean ni after 'fice_force' on ni01  
   real(r8) ni_0(pcols,pver)                               ! Grid-mean ni after 'instratus_condensate' on ni1
   real(r8) ni_05(pcols,pver)                              ! Grid-mean ni after 'advection' on ni_0 
   real(r8) ni_prime0(pcols,pver)                          ! Grid-mean ni after 'Macrophysics (QQ)' on ni_05star
   real(r8) ni_dprime(pcols,pver)                          ! Grid-mean ni after 'fice_force' on ni_prime
   real(r8) ni_star(pcols,pver)                            ! Grid-mean ni after 'instratus_condensate' on ni_dprime

   real(r8) a_st(pcols,pver)                               ! Stratus fraction of equilibrium reference state 
   real(r8) a_st_0(pcols,pver)                             ! Stratus fraction at '_0' state
   real(r8) a_st_star(pcols,pver)                          ! Stratus fraction at '_star' state

   real(r8) al_st(pcols,pver)                              ! Liquid stratus fraction of equilibrium reference state 
   real(r8) al_st_0(pcols,pver)                            ! Liquid stratus fraction at '_0' state
   real(r8) al_st_nc(pcols,pver)                           ! Non-physical liquid stratus fraction in the non-cumulus pixels

   real(r8) ai_st(pcols,pver)                              ! Ice stratus fraction of equilibrium reference state 
   real(r8) ai_st_0(pcols,pver)                            ! Ice stratus fraction at '_0' state
   real(r8) ai_st_nc(pcols,pver)                           ! Non-physical ice stratus fraction in the non-cumulus pixels

   real(r8) ql_st(pcols,pver)                              ! In-stratus LWC of equilibrium reference state [kg/kg] 
   real(r8) ql_st_0(pcols,pver)                            ! In-stratus LWC at '_0' state

   real(r8) qi_st(pcols,pver)                              ! In-stratus IWC of equilibrium reference state [kg/kg] 
   real(r8) qi_st_0(pcols,pver)                            ! In-stratus IWC at '_0' state

 ! Cumulus properties 

   real(r8) a_cu(pcols,pver)

 ! Adjustment tendency in association with 'positive_moisture'

   real(r8) Tten_pwi1(pcols,pver)                          ! Pre-process T  tendency of input equilibrium state [K/s] 
   real(r8) qvten_pwi1(pcols,pver)                         ! Pre-process qv tendency of input equilibrium state [kg/kg/s]
   real(r8) qlten_pwi1(pcols,pver)                         ! Pre-process ql tendency of input equilibrium state [kg/kg/s]
   real(r8) qiten_pwi1(pcols,pver)                         ! Pre-process qi tendency of input equilibrium state [kg/kg/s]
   real(r8) nlten_pwi1(pcols,pver)                         ! Pre-process nl tendency of input equilibrium state [#/kg/s]
   real(r8) niten_pwi1(pcols,pver)                         ! Pre-process ni tendency of input equilibrium state [#/kg/s] 

   real(r8) Tten_pwi2(pcols,pver)                          ! Post-process T  tendency of provisional equilibrium state [K/s] 
   real(r8) qvten_pwi2(pcols,pver)                         ! Post-process qv tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qlten_pwi2(pcols,pver)                         ! Post-process ql tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qiten_pwi2(pcols,pver)                         ! Post-process qi tendency of provisional equilibrium state [kg/kg/s]
   real(r8) nlten_pwi2(pcols,pver)                         ! Post-process nl tendency of provisoonal equilibrium state [#/kg/s]
   real(r8) niten_pwi2(pcols,pver)                         ! Post-process ni tendency of provisional equilibrium state [#/kg/s] 

 ! Adjustment tendency in association with 'instratus_condensate'

   real(r8) QQw1(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi1(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQw2(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi2(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]

   real(r8) QQnl1(pcols,pver)          ! Tendency of nl associated with QQw1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQni1(pcols,pver)          ! Tendency of ni associated with QQi1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQnl2(pcols,pver)          ! Tendency of nl associated with QQw2 only when QQw2<0 (net evaporation) [#/kg/s]
   real(r8) QQni2(pcols,pver)          ! Tendency of ni associated with QQi2 only when QQw2<0 (net evaporation) [#/kg/s]

 ! Macrophysical process tendency variables

   real(r8) QQ(pcols,pver)             ! Net condensation rate into water+ice           [kg/kg/s] 
   real(r8) QQw(pcols,pver)            ! Net condensation rate into water               [kg/kg/s] 
   real(r8) QQi(pcols,pver)            ! Net condensation rate into ice                 [kg/kg/s]
   real(r8) QQnl(pcols,pver)           ! Tendency of nl associated with QQw both for condensation and evaporation [#/kg/s]
   real(r8) QQni(pcols,pver)           ! Tendency of ni associated with QQi both for condensation and evaporation [#/kg/s]
   real(r8) ACnl(pcols,pver)           ! Cloud liquid droplet (nl) activation tendency [#/kg/s]
   real(r8) ACni(pcols,pver)           ! Cloud ice    droplet (ni) activation tendency [#/kg/s]

   real(r8) QQw_prev(pcols,pver)   
   real(r8) QQi_prev(pcols,pver)   
   real(r8) QQnl_prev(pcols,pver)  
   real(r8) QQni_prev(pcols,pver)  

   real(r8) QQw_prog(pcols,pver)   
   real(r8) QQi_prog(pcols,pver)   
   real(r8) QQnl_prog(pcols,pver)  
   real(r8) QQni_prog(pcols,pver)  

   real(r8) QQ_final(pcols,pver)                           
   real(r8) QQw_final(pcols,pver)                           
   real(r8) QQi_final(pcols,pver)                           
   real(r8) QQn_final(pcols,pver)                           
   real(r8) QQnl_final(pcols,pver)                          
   real(r8) QQni_final(pcols,pver)                          

   real(r8) QQ_all(pcols,pver)         ! QQw_all    + QQi_all
   real(r8) QQw_all(pcols,pver)        ! QQw_final  + QQw1  + QQw2  + qlten_pwi1 + qlten_pwi2 + A_ql_adj [kg/kg/s]
   real(r8) QQi_all(pcols,pver)        ! QQi_final  + QQi1  + QQi2  + qiten_pwi1 + qiten_pwi2 + A_qi_adj [kg/kg/s]
   real(r8) QQn_all(pcols,pver)        ! QQnl_all   + QQni_all
   real(r8) QQnl_all(pcols,pver)       ! QQnl_final + QQnl1 + QQnl2 + nlten_pwi1 + nlten_pwi2 + ACnl [#/kg/s]
   real(r8) QQni_all(pcols,pver)       ! QQni_final + QQni1 + QQni2 + niten_pwi1 + niten_pwi2 + ACni [#/kg/s]

 ! Coefficient for computing QQ and related processes

   real(r8) U(pcols,pver)                                  ! Grid-mean RH
   real(r8) U_nc(pcols,pver)                               ! Mean RH of non-cumulus pixels
   real(r8) G_nc(pcols,pver)                               ! d(U_nc)/d(a_st_nc)


   real(r8) zeros(pcols,pver)

   real(r8) qmin1(pcols,pver)
   real(r8) qmin2(pcols,pver)
   real(r8) qmin3(pcols,pver)

   real(r8) esat_a(pcols)                                  ! Saturation water vapor pressure [Pa]
   real(r8) qsat_a(pcols)                                  ! Saturation water vapor specific humidity [kg/kg]
   real(r8) dqsdT_a(pcols)                                 ! dqsat/dT [kg/kg/K]
   real(r8) Twb_aw(pcols,pver)                             ! Wet-bulb temperature [K]
   real(r8) qvwb_aw(pcols,pver)                            ! Wet-bulb water vapor specific humidity [kg/kg]

   real(r8) esat_b(pcols)                                 
   real(r8) qsat_b(pcols)                                 
   real(r8) dqsdT_b(pcols)                                 

   real(r8) QQmax,QQmin,QQwmin,QQimin                      ! For limiting QQ
   real(r8) cone                                           ! Number close to but smaller than 1
   real(r8) qsmall                                         ! Smallest mixing ratio considered in the macrophysics

   ! qinyi
   real(r8) conden_st_nc(pcols,pver)

   real(r8) G_nc_def(pcols,pver)

   conden_st_nc(:ncol,:) = 0._r8

   sgm1(:ncol,:) = 0._r8
   sgm(:ncol,:) = 0._r8
   sgm_dprime(:ncol,:) = 0._r8

   G_nc_def(:ncol,:) = 0._r8
   
   deltaq_sat(:ncol,:) = 0._r8
   deltaq_uns(:ncol,:) = 0._r8
   Q1_sat(:ncol,:) = 0._r8
   Q1_uns(:ncol,:) = 0._r8
   N1(:ncol,:) = 0._r8
   ! qinyi 2017-1-1 22:02:12
   adjust_factor(:ncol,:) = 0._r8

   alst_def(:ncol,:) = 0._r8


   cone            = 0.999_r8
   qsmall          = 1.e-18_r8
   zeros(:ncol,:)  = 0._r8

   ! ------------------------------------ !
   ! Global initialization of main output !
   ! ------------------------------------ !

     s_tendout(:ncol,:)     = 0._r8
     qv_tendout(:ncol,:)    = 0._r8
     ql_tendout(:ncol,:)    = 0._r8
     qi_tendout(:ncol,:)    = 0._r8
     nl_tendout(:ncol,:)    = 0._r8
     ni_tendout(:ncol,:)    = 0._r8

     qme(:ncol,:)           = 0._r8

     cld(:ncol,:)           = 0._r8
     al_st_star(:ncol,:)    = 0._r8
     ai_st_star(:ncol,:)    = 0._r8
     ql_st_star(:ncol,:)    = 0._r8
     qi_st_star(:ncol,:)    = 0._r8

! qinyi
     alst_gp(:ncol,:)       = 0._r8
     cond_gp(:ncol,:)       = 0._r8
     alst_def(:ncol,:)      = 0._r8
     N1(:ncol,:)            = 0._r8

     sgm1(:ncol,:)          = 0._r8
     sgm(:ncol,:)           = 0._r8
     sgm_dprime(:ncol,:)    = 0._r8

     delta_q(:ncol,:)       = 0._r8
     deltaq_sat(:ncol,:)    = 0._r8
     deltaq_uns(:ncol,:)    = 0._r8
     Q1_sat(:ncol,:)        = 0._r8
     Q1_uns(:ncol,:)        = 0._r8
     adjust_factor(:ncol,:) = 0._r8


   ! --------------------------------------- !
   ! Initialization of internal 2D variables !
   ! --------------------------------------- !

     T(:ncol,:)             = 0._r8
     T1(:ncol,:)            = 0._r8
     T_0(:ncol,:)           = 0._r8
     T_05(:ncol,:)          = 0._r8
     T_prime0(:ncol,:)      = 0._r8
     T_dprime(:ncol,:)      = 0._r8
     T_star(:ncol,:)        = 0._r8

     qv(:ncol,:)            = 0._r8
     qv1(:ncol,:)           = 0._r8
     qv_0(:ncol,:)          = 0._r8
     qv_05(:ncol,:)         = 0._r8
     qv_prime0(:ncol,:)     = 0._r8
     qv_dprime(:ncol,:)     = 0._r8
     qv_star(:ncol,:)       = 0._r8

     ql(:ncol,:)            = 0._r8
     ql1(:ncol,:)           = 0._r8
     ql_0(:ncol,:)          = 0._r8
     ql_05(:ncol,:)         = 0._r8
     ql_prime0(:ncol,:)     = 0._r8
     ql_dprime(:ncol,:)     = 0._r8
     ql_star(:ncol,:)       = 0._r8

     qi(:ncol,:)            = 0._r8
     qi1(:ncol,:)           = 0._r8
     qi_0(:ncol,:)          = 0._r8
     qi_05(:ncol,:)         = 0._r8
     qi_prime0(:ncol,:)     = 0._r8
     qi_dprime(:ncol,:)     = 0._r8
     qi_star(:ncol,:)       = 0._r8

     nl(:ncol,:)            = 0._r8
     nl1(:ncol,:)           = 0._r8
     nl_0(:ncol,:)          = 0._r8
     nl_05(:ncol,:)         = 0._r8
     nl_prime0(:ncol,:)     = 0._r8
     nl_dprime(:ncol,:)     = 0._r8
     nl_star(:ncol,:)       = 0._r8

     ni(:ncol,:)            = 0._r8
     ni1(:ncol,:)           = 0._r8
     ni_0(:ncol,:)          = 0._r8
     ni_05(:ncol,:)         = 0._r8
     ni_prime0(:ncol,:)     = 0._r8
     ni_dprime(:ncol,:)     = 0._r8
     ni_star(:ncol,:)       = 0._r8

     a_st(:ncol,:)          = 0._r8
     a_st_0(:ncol,:)        = 0._r8
     a_st_star(:ncol,:)     = 0._r8

     al_st(:ncol,:)         = 0._r8
     al_st_0(:ncol,:)       = 0._r8
     al_st_nc(:ncol,:)      = 0._r8

     ai_st(:ncol,:)         = 0._r8
     ai_st_0(:ncol,:)       = 0._r8
     ai_st_nc(:ncol,:)      = 0._r8

     ql_st(:ncol,:)         = 0._r8
     ql_st_0(:ncol,:)       = 0._r8

     qi_st(:ncol,:)         = 0._r8
     qi_st_0(:ncol,:)       = 0._r8

 ! Cumulus properties 

     a_cu(:ncol,:)          = 0._r8

 ! Adjustment tendency in association with 'positive_moisture'

     Tten_pwi1(:ncol,:)     = 0._r8
     qvten_pwi1(:ncol,:)    = 0._r8
     qlten_pwi1(:ncol,:)    = 0._r8
     qiten_pwi1(:ncol,:)    = 0._r8
     nlten_pwi1(:ncol,:)    = 0._r8
     niten_pwi1(:ncol,:)    = 0._r8

     Tten_pwi2(:ncol,:)     = 0._r8
     qvten_pwi2(:ncol,:)    = 0._r8
     qlten_pwi2(:ncol,:)    = 0._r8
     qiten_pwi2(:ncol,:)    = 0._r8
     nlten_pwi2(:ncol,:)    = 0._r8
     niten_pwi2(:ncol,:)    = 0._r8

     qvadj   (:ncol,:)      = 0._r8
     qladj   (:ncol,:)      = 0._r8
     qiadj   (:ncol,:)      = 0._r8

 ! Adjustment tendency in association with 'instratus_condensate'

     QQw1(:ncol,:)          = 0._r8
     QQi1(:ncol,:)          = 0._r8
     QQw2(:ncol,:)          = 0._r8
     QQi2(:ncol,:)          = 0._r8

     QQnl1(:ncol,:)         = 0._r8
     QQni1(:ncol,:)         = 0._r8
     QQnl2(:ncol,:)         = 0._r8
     QQni2(:ncol,:)         = 0._r8

     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8

 ! Macrophysical process tendency variables

     QQ(:ncol,:)            = 0._r8
     QQw(:ncol,:)           = 0._r8
     QQi(:ncol,:)           = 0._r8
     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8
     ACnl(:ncol,:)          = 0._r8
     ACni(:ncol,:)          = 0._r8

     QQw_prev(:ncol,:)      = 0._r8
     QQi_prev(:ncol,:)      = 0._r8
     QQnl_prev(:ncol,:)     = 0._r8
     QQni_prev(:ncol,:)     = 0._r8

     QQw_prog(:ncol,:)      = 0._r8
     QQi_prog(:ncol,:)      = 0._r8
     QQnl_prog(:ncol,:)     = 0._r8
     QQni_prog(:ncol,:)     = 0._r8

     QQ_final(:ncol,:)      = 0._r8                        
     QQw_final(:ncol,:)     = 0._r8                  
     QQi_final(:ncol,:)     = 0._r8           
     QQn_final(:ncol,:)     = 0._r8    
     QQnl_final(:ncol,:)    = 0._r8
     QQni_final(:ncol,:)    = 0._r8

     QQ_all(:ncol,:)        = 0._r8
     QQw_all(:ncol,:)       = 0._r8
     QQi_all(:ncol,:)       = 0._r8
     QQn_all(:ncol,:)       = 0._r8
     QQnl_all(:ncol,:)      = 0._r8
     QQni_all(:ncol,:)      = 0._r8

 ! Coefficient for computing QQ and related processes

     U(:ncol,:)             = 0._r8
     U_nc(:ncol,:)          = 0._r8
     G_nc(:ncol,:)          = 0._r8

 ! Other

     qmin1(:ncol,:)         = 0._r8
     qmin2(:ncol,:)         = 0._r8
     qmin3(:ncol,:)         = 0._r8

     esat_b(:ncol)          = 0._r8     
     qsat_b(:ncol)          = 0._r8    
     dqsdT_b(:ncol)         = 0._r8 

     esat_a(:ncol)          = 0._r8     
     qsat_a(:ncol)          = 0._r8    
     dqsdT_a(:ncol)         = 0._r8 
     Twb_aw(:ncol,:)        = 0._r8
     qvwb_aw(:ncol,:)       = 0._r8

   ! ---------------- !
   ! Main computation ! 
   ! ---------------- !

   ! ---------------------------------------------------------------------- !
   ! set to zero for levels above
   ! ---------------------------------------------------------------------- !
   ql0(:ncol,:top_lev-1) = 0._r8
   qi0(:ncol,:top_lev-1) = 0._r8
   nl0(:ncol,:top_lev-1) = 0._r8
   ni0(:ncol,:top_lev-1) = 0._r8
   
   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !                
   ! ---------------------------------------------------------------------- !

   T1(:ncol,:)    =  T0(:ncol,:) 
   qv1(:ncol,:)   = qv0(:ncol,:) 
   ql1(:ncol,:)   = ql0(:ncol,:) 
   qi1(:ncol,:)   = qi0(:ncol,:) 
   nl1(:ncol,:)   = nl0(:ncol,:) 
   ni1(:ncol,:)   = ni0(:ncol,:) 

   
   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )


   qmin1(:ncol,:) = qmin(1)
   qmin2(:ncol,:) = qmin(ixcldliq)
   qmin3(:ncol,:) = qmin(ixcldice)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv1, ql1, qi1, T1, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql1(i,k) .lt. qsmall ) then
          nlten_pwi1(i,k) = -nl1(i,k)/dt
          nl1(i,k)        = 0._r8
      endif 
      if( qi1(i,k) .lt. qsmall ) then
          niten_pwi1(i,k) = -ni1(i,k)/dt
          ni1(i,k)        = 0._r8
      endif 
   enddo ! do i=1,ncol
   enddo ! do k=top_lev,pver

   ! ------------------------------------------------------------- !
   ! Impose 'in-stratus condensate amount constraint'              !
   ! such that it is bounded by two limiting values.               !      
   ! This should also use 'old' cumulus properties since it is     !
   ! before applying external forcings.                            ! 
   ! Below 'QQw1,QQi1' are effective adjustive condensation        ! 
   ! Although this process also involves freezing of cloud         !
   ! liquid into ice, they can be and only can be expressed        !
   ! in terms of effective condensation.                           !
   ! ------------------------------------------------------------- !

! qinyi
   call subgrid_var(lchnk, ncol, rpdel, p, zm, T1, qv1, ql1, &
                    tke, kvh, kvm, lengi, shi, smi, wstarPBL, pblh, &
                    qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal,&
                    explwi2qtu, & ! qinyi 2017-11-30 21:48:04
                    sgm1)

   do k = top_lev, pver
! qinyi
! although using this subroutine, it does not have any impact on the
! temperature, specific humidity et al. just for keeping this structure as
! before.
     call instratus_condensate1( lchnk, ncol, k,                                   &
                                 p(:,k), T1(:,k), qv1(:,k), ql1(:,k), qi1(:,k),    &
                                 a_cud(:,k), zeros(:,k), zeros(:,k),               &
                                 zeros(:,k), zeros(:,k), zeros(:,k),               &
                                 landfrac, snowh,                                  &
                                 T_0(:,k), qv_0(:,k), ql_0(:,k), qi_0(:,k),        & 
                                 al_st_0(:,k), ai_st_0(:,k), ql_st_0(:,k), qi_st_0(:,k), sgm1(:,k))

      a_st_0(:ncol,k) = max(al_st_0(:ncol,k),ai_st_0(:ncol,k))
      QQw1(:ncol,k)   = (ql_0(:ncol,k) - ql1(:ncol,k))/dt
      QQi1(:ncol,k)   = (qi_0(:ncol,k) - qi1(:ncol,k))/dt

      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! Set a limit such that negative state not happens.  ! 
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw1(i,k) .le. 0._r8 ) then
             if( ql1(i,k) .gt. qsmall ) then
                 QQnl1(i,k) = QQw1(i,k)*nl1(i,k)/ql1(i,k)
                 QQnl1(i,k) = min(0._r8,cone*max(QQnl1(i,k),-nl1(i,k)/dt))
             else
                 QQnl1(i,k) = 0._r8
             endif  
         endif 
         if( QQi1(i,k) .le. 0._r8 ) then
             if( qi1(i,k) .gt. qsmall ) then
                 QQni1(i,k) = QQi1(i,k)*ni1(i,k)/qi1(i,k)
                 QQni1(i,k) = min(0._r8,cone*max(QQni1(i,k),-ni1(i,k)/dt))
             else
                 QQni1(i,k) = 0._r8
             endif  
         endif 
      enddo ! do i=1,ncol

   enddo !do k=top_lev,pver

   nl_0(:ncol,top_lev:) = max(0._r8,nl1(:ncol,top_lev:)+QQnl1(:ncol,top_lev:)*dt) 
   ni_0(:ncol,top_lev:) = max(0._r8,ni1(:ncol,top_lev:)+QQni1(:ncol,top_lev:)*dt)

   ! -------------------------------------------------------------- !
   ! the main part starts.    
   ! -------------------------------------------------------------- !

   T(:ncol,top_lev:)     = T_0(:ncol,top_lev:)
   qv(:ncol,top_lev:)    = qv_0(:ncol,top_lev:)
   ql(:ncol,top_lev:)    = ql_0(:ncol,top_lev:)
   qi(:ncol,top_lev:)    = qi_0(:ncol,top_lev:)
   al_st(:ncol,top_lev:) = al_st_0(:ncol,top_lev:)
   ai_st(:ncol,top_lev:) = ai_st_0(:ncol,top_lev:)
   a_st(:ncol,top_lev:)  = a_st_0(:ncol,top_lev:)
   ql_st(:ncol,top_lev:) = ql_st_0(:ncol,top_lev:)
   qi_st(:ncol,top_lev:) = qi_st_0(:ncol,top_lev:)
   nl(:ncol,top_lev:)    = nl_0(:ncol,top_lev:)
   ni(:ncol,top_lev:)    = ni_0(:ncol,top_lev:)

! qinyi

call subgrid_var(lchnk, ncol, rpdel, p, zm, T, qv, ql, &
                 tke, kvh, kvm, lengi, shi, smi, wstarPBL, pblh, &
                 qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                 explwi2qtu, & ! qinyi 2017-11-30 21:48:23
                 sgm)

do k = top_lev, pver
      ! "False" means that ice will not be taken into account.
      call findsp_vc(qv(:ncol,k),T(:ncol,k),p(:ncol,k), .false., &
           Twb_aw(:ncol,k),qvwb_aw(:ncol,k))

      call qsat_water(T(1:ncol,k), p(1:ncol,k), &
           esat_a(1:ncol), qsat_a(1:ncol), dqsdt=dqsdT_a(1:ncol))
      call qsat_water(T(1:ncol,k), p(1:ncol,k), &
           esat_b(1:ncol), qsat_b(1:ncol), dqsdt=dqsdT_b(1:ncol))

          a_cu(:ncol,k) = a_cud(:ncol,k)

      do i = 1, ncol
         U(i,k)    =  qv(i,k)/qsat_b(i)
         U_nc(i,k) =  U(i,k)
      enddo

      call cldfrac_conden(p(:,k),T(:,k),qv(:,k),ql(:,k),sgm(:,k),al_st_nc(:,k),conden_st_nc(:,k),G_nc(:,k),ncol)
      call aist_vector(qv(:,k),T(:,k),p(:,k),qi(:,k),landfrac(:),snowh(:),ai_st_nc(:,k),ncol)

      ai_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*ai_st_nc(:ncol,k)
      al_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*al_st_nc(:ncol,k)
      a_st(:ncol,k)   =  max(al_st(:ncol,k),ai_st(:ncol,k))  

      ! use default cloud fraction scheme
      call astG_PDF(U_nc(:,k),p(:,k),qv(:,k),landfrac(:),snowh(:),alst_def(:,k),G_nc_def(:,k),ncol)


      do i = 1, ncol
        ! qinyi
         call cldfrac_conden_single_only(p(i,k),T(i,k),qv(i,k),ql(i,k),sgm(i,k),al_st_nc(i,k),conden_st_nc(i,k),G_nc(i,k),&
                                         deltaq_sat(i,k),deltaq_uns(i,k),Q1_sat(i,k),Q1_uns(i,k),&
                                         adjust_factor(i,k))
        ! qinyi 2017-1-1 22:01:46
        cond_gp(i,k) = adjust_factor(i,k)*(conden_st_nc(i,k)-ql(i,k))/dt

        QQ(i,k) = cond_gp(i,k)

        ! in order to diagnose the output cloud fraction and condensate
        alst_gp(i,k) = al_st_nc(i,k)
        cond_gp(i,k) = conden_st_nc(i,k) 
        delta_q(i,k) = G_nc(i,k)

        if(delta_q(i,k).ge.0._r8)then
            N1(i,k) = 1._r8
        else
            N1(i,k) = 0._r8
        endif

        ! adjust tendency
        if( QQ(i,k) .ge. 0._r8 ) then
            QQmax    = (qv(i,k) - qmin(1))/dt ! For ghost cumulus & semi-ghost ice stratus
            QQmax    = max(0._r8,QQmax) 
            QQ(i,k)  = min(QQ(i,k),QQmax)
            QQw(i,k) = QQ(i,k)
            QQi(i,k) = 0._r8 
        else
            QQmin  = 0._r8
            if( qv(i,k) .lt. qsat_a(i) ) QQmin = min(0._r8,cone*(qv(i,k)-qvwb_aw(i,k))/dt)
            QQ(i,k)  = max(QQ(i,k),QQmin)
            QQw(i,k) = QQ(i,k)
            QQi(i,k) = 0._r8
            QQwmin   = min(0._r8,-cone*ql(i,k)/dt)
            QQimin   = min(0._r8,-cone*qi(i,k)/dt)
            QQw(i,k) = min(0._r8,max(QQw(i,k),QQwmin))
            QQi(i,k) = min(0._r8,max(QQi(i,k),QQimin))
        endif

       ! -------------------------------------------------- !
       ! Reduce droplet concentration if evaporation occurs !
       ! -------------------------------------------------- !

         if( QQw(i,k) .lt. 0._r8 ) then
             if( ql(i,k) .gt. qsmall ) then
                 QQnl(i,k) = QQw(i,k)*nl(i,k)/ql(i,k)
                 QQnl(i,k) = min(0._r8,cone*max(QQnl(i,k),-nl(i,k)/dt))
             else
                 QQnl(i,k) = 0._r8
             endif  
         endif 

         if( QQi(i,k) .lt. 0._r8 ) then
             if( qi(i,k) .gt. qsmall ) then
                 QQni(i,k) = QQi(i,k)*ni(i,k)/qi(i,k)
                 QQni(i,k) = min(0._r8,cone*max(QQni(i,k),-ni(i,k)/dt))
             else
                 QQni(i,k) = 0._r8
             endif  
         endif 

      enddo ! do i=1,ncol
enddo !do k=top_lev,pver

    ! -------------------------------------------------------------------- !
    ! Until now, we have finished computing all necessary tendencies       ! 
    ! from the equilibrium input state (T_0).                              !
    ! -------------------------------------------------------------------- !

      QQw_prog(:ncol,top_lev:)   = QQw(:ncol,top_lev:) 
      QQi_prog(:ncol,top_lev:)   = QQi(:ncol,top_lev:) 
      QQnl_prog(:ncol,top_lev:)  = QQnl(:ncol,top_lev:)
      QQni_prog(:ncol,top_lev:)  = QQni(:ncol,top_lev:)

    ! -------------------------------------------------------- !
    ! Compute final prognostic state on which final diagnosts  !
    ! -------------------------------------------------------- !

    do k = top_lev, pver
    do i = 1, ncol
       T_prime0(i,k)  = T_0(i,k)  + dt*(  &
            (latvap*QQw_prog(i,k)+(latvap+latice)*QQi_prog(i,k))/cpair )
       qv_prime0(i,k) = qv_0(i,k) + dt*( 0._r8 - QQw_prog(i,k) - QQi_prog(i,k) )
       ql_prime0(i,k) = ql_0(i,k) + dt*( 0._r8 + QQw_prog(i,k) )
       qi_prime0(i,k) = qi_0(i,k) + dt*( 0._r8 + QQi_prog(i,k) )
       nl_prime0(i,k) = max(0._r8,nl_0(i,k) + dt*( 0._r8 + QQnl_prog(i,k) ))
       ni_prime0(i,k) = max(0._r8,ni_0(i,k) + dt*( 0._r8 + QQni_prog(i,k) ))
       if( ql_prime0(i,k) .lt. qsmall ) nl_prime0(i,k) = 0._r8
       if( qi_prime0(i,k) .lt. qsmall ) ni_prime0(i,k) = 0._r8
    enddo
    enddo

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !

   T_dprime(:ncol,top_lev:)  =  T_prime0(:ncol,top_lev:) 
   qv_dprime(:ncol,top_lev:) = qv_prime0(:ncol,top_lev:) 
   ql_dprime(:ncol,top_lev:) = ql_prime0(:ncol,top_lev:) 
   qi_dprime(:ncol,top_lev:) = qi_prime0(:ncol,top_lev:) 
   nl_dprime(:ncol,top_lev:) = nl_prime0(:ncol,top_lev:) 
   ni_dprime(:ncol,top_lev:) = ni_prime0(:ncol,top_lev:) 

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp,          & 
                           qv_dprime, ql_dprime, qi_dprime, T_dprime,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice)

   do k = top_lev, pver
   do i = 1, ncol
      if( ql_dprime(i,k) .lt. qsmall ) then
          nlten_pwi2(i,k) = -nl_dprime(i,k)/dt
          nl_dprime(i,k)   = 0._r8
      endif 
      if( qi_dprime(i,k) .lt. qsmall ) then
          niten_pwi2(i,k) = -ni_dprime(i,k)/dt
          ni_dprime(i,k)   = 0._r8
      endif 
   enddo
   enddo

   ! ---------------------------------------------------------- !
   ! Impose diagnostic upper and lower limits on the in-stratus !
   ! condensate amount. This produces a final equilibrium state !
   ! at the end of each iterative process.                      !
   ! ---------------------------------------------------------- !

   ! qinyi
   call subgrid_var(lchnk,ncol,rpdel,p,zm,T_0,qv_0,ql_0, &
                    tke,kvh,kvm,lengi,shi,smi,wstarPBL,pblh, &
                    qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                    explwi2qtu, & ! qinyi 2017-11-30 21:48:47
                    sgm_dprime)

   do k = top_lev, pver

! qinyi
! currently this subroutine is just used to calculate the cloud fraction one more time (actually it is same to the before:al_st_nc)

     call instratus_condensate1( lchnk          , ncol           , k              , p(:,k)        , &
                                 T_0(:,k)       , qv_0(:,k)      , ql_0(:,k)      , qi_0(:,k),      &
                                 a_cu0(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 & 
                                 zeros(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 &
                                 landfrac       , snowh          ,                                  &
                                 T_star(:,k)    , qv_star(:,k)   , ql_star(:,k)   , qi_star(:,k)  , & 
                                 al_st_star(:,k), ai_st_star(:,k), ql_st_star(:,k), qi_st_star(:,k) , sgm_dprime(:,k))
 
      a_st_star(:ncol,k)  = max(al_st_star(:ncol,k),ai_st_star(:ncol,k))

      ! qinyi added in order to fully eliminate tendency QQw2
      T_star(:ncol,k)  = T_dprime(:ncol,k)
      qv_star(:ncol,k) = qv_dprime(:ncol,k)
      ql_star(:ncol,k) = ql_dprime(:ncol,k)
      qi_star(:ncol,k) = qi_dprime(:ncol,k)

      QQw2(:ncol,k) = (ql_star(:ncol,k) - ql_dprime(:ncol,k))/dt
      QQi2(:ncol,k) = (qi_star(:ncol,k) - qi_dprime(:ncol,k))/dt


      ! -------------------------------------------------- !
      ! Reduce droplet concentration if evaporation occurs !
      ! -------------------------------------------------- !
      do i = 1, ncol
         if( QQw2(i,k) .le. 0._r8 ) then
             if( ql_dprime(i,k) .ge. qsmall ) then
                 QQnl2(i,k) = QQw2(i,k)*nl_dprime(i,k)/ql_dprime(i,k)
                 QQnl2(i,k) = min(0._r8,cone*max(QQnl2(i,k),-nl_dprime(i,k)/dt))
             else
                 QQnl2(i,k) = 0._r8
             endif  
         endif 
         if( QQi2(i,k) .le. 0._r8 ) then
             if( qi_dprime(i,k) .gt. qsmall ) then
                 QQni2(i,k) = QQi2(i,k)*ni_dprime(i,k)/qi_dprime(i,k)
                 QQni2(i,k) = min(0._r8,cone*max(QQni2(i,k),-ni_dprime(i,k)/dt))
             else
                 QQni2(i,k) = 0._r8
             endif  
         endif 
      enddo
   enddo
   nl_star(:ncol,top_lev:) = max(0._r8,nl_dprime(:ncol,top_lev:)+QQnl2(:ncol,top_lev:)*dt) 
   ni_star(:ncol,top_lev:) = max(0._r8,ni_dprime(:ncol,top_lev:)+QQni2(:ncol,top_lev:)*dt)

   ! ------------------------------------------ !
   ! Final adjustment of droplet concentration. !
   ! Set # to zero if there is no cloud.        !
   ! ------------------------------------------ !

   do k = top_lev, pver
   do i = 1, ncol 
      if( ql_star(i,k) .lt. qsmall ) then
          ACnl(i,k) = - nl_star(i,k)/dt
          nl_star(i,k) = 0._r8
      endif
      if( qi_star(i,k) .lt. qsmall ) then
          ACni(i,k) = - ni_star(i,k)/dt
          ni_star(i,k) = 0._r8
      endif
   enddo
   enddo

   ! ------------------------------------------------------------------------ !
   ! Compute final tendencies of main output variables and diagnostic outputs !
   ! ------------------------------------------------------------------------ !

   ! ------------------ !
   ! Process tendencies !
   ! ------------------ !

   QQw_final(:ncol,top_lev:)  = QQw_prog(:ncol,top_lev:)
   QQi_final(:ncol,top_lev:)  = QQi_prog(:ncol,top_lev:)
   QQ_final(:ncol,top_lev:)   = QQw_final(:ncol,top_lev:) + QQi_final(:ncol,top_lev:)
   QQw_all(:ncol,top_lev:)    = QQw_prog(:ncol,top_lev:)  + QQw1(:ncol,top_lev:) + QQw2(:ncol,top_lev:) + &
        qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:)
   QQi_all(:ncol,top_lev:)    = QQi_prog(:ncol,top_lev:)  + QQi1(:ncol,top_lev:) + QQi2(:ncol,top_lev:) + &
        qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:)
   QQ_all(:ncol,top_lev:)     = QQw_all(:ncol,top_lev:)   + QQi_all(:ncol,top_lev:)
   QQnl_final(:ncol,top_lev:) = QQnl_prog(:ncol,top_lev:)
   QQni_final(:ncol,top_lev:) = QQni_prog(:ncol,top_lev:)
   QQn_final(:ncol,top_lev:)  = QQnl_final(:ncol,top_lev:) + QQni_final(:ncol,top_lev:)
   QQnl_all(:ncol,top_lev:)   = QQnl_prog(:ncol,top_lev:)  + QQnl1(:ncol,top_lev:) + QQnl2(:ncol,top_lev:) + &
        nlten_pwi1(:ncol,top_lev:) + nlten_pwi2(:ncol,top_lev:) + ACnl(:ncol,top_lev:)
   QQni_all(:ncol,top_lev:)   = QQni_prog(:ncol,top_lev:)  + QQni1(:ncol,top_lev:) + QQni2(:ncol,top_lev:) + &
        niten_pwi1(:ncol,top_lev:) + niten_pwi2(:ncol,top_lev:) + ACni(:ncol,top_lev:)
   QQn_all(:ncol,top_lev:)    = QQnl_all(:ncol,top_lev:)   + QQni_all(:ncol,top_lev:)
   qme(:ncol,top_lev:)        = QQ_final(:ncol,top_lev:)   
   qvadj(:ncol,top_lev:)      = qvten_pwi1(:ncol,top_lev:) + qvten_pwi2(:ncol,top_lev:)
   qladj(:ncol,top_lev:)      = qlten_pwi1(:ncol,top_lev:) + qlten_pwi2(:ncol,top_lev:)
   qiadj(:ncol,top_lev:)      = qiten_pwi1(:ncol,top_lev:) + qiten_pwi2(:ncol,top_lev:)
   qllim(:ncol,top_lev:)      = QQw1      (:ncol,top_lev:) + QQw2      (:ncol,top_lev:)
   qilim(:ncol,top_lev:)      = QQi1      (:ncol,top_lev:) + QQi2      (:ncol,top_lev:)

   ! ----------------- !
   ! Output tendencies !
   ! ----------------- !

   s_tendout(:ncol,top_lev:)  = cpair*( T_dprime(:ncol,top_lev:)  -  T0(:ncol,top_lev:) )/dt
   qv_tendout(:ncol,top_lev:) =       ( qv_dprime(:ncol,top_lev:) - qv0(:ncol,top_lev:) )/dt
   ql_tendout(:ncol,top_lev:) =       ( ql_dprime(:ncol,top_lev:) - ql0(:ncol,top_lev:) )/dt
   qi_tendout(:ncol,top_lev:) =       ( qi_dprime(:ncol,top_lev:) - qi0(:ncol,top_lev:) )/dt

   nl_tendout(:ncol,top_lev:) =       ( nl_dprime(:ncol,top_lev:) - nl0(:ncol,top_lev:) )/dt
   ni_tendout(:ncol,top_lev:) =       ( ni_dprime(:ncol,top_lev:) - ni0(:ncol,top_lev:) )/dt

   if (.not. do_cldice) then
      do k = top_lev, pver
         do i = 1, ncol
            ! Don't want either qi or ni tendencies, but the code above is somewhat convoluted and
            ! is trying to adjust both (small numbers). Just force it to zero here.
            qi_tendout(i,k) = 0._r8
            ni_tendout(i,k) = 0._r8
          end do
      end do
   end if

   ! ------------------ !
   ! Net cloud fraction !
   ! ------------------ !

   cld(:ncol,top_lev:) = a_st_star(:ncol,top_lev:) + a_cu0(:ncol,top_lev:)

   ! --------------------------------- !
   ! Updated grid-mean state variables !
   ! --------------------------------- !

   T0(:ncol,top_lev:)  = T_dprime(:ncol,top_lev:)
   qv0(:ncol,top_lev:) = qv_dprime(:ncol,top_lev:)
   ql0(:ncol,top_lev:) = ql_dprime(:ncol,top_lev:)
   qi0(:ncol,top_lev:) = qi_dprime(:ncol,top_lev:)
   nl0(:ncol,top_lev:) = nl_dprime(:ncol,top_lev:)
   ni0(:ncol,top_lev:) = ni_dprime(:ncol,top_lev:)


   return
   end subroutine mmacro_pcond

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_condensate1( lchnk, ncol, k,                     &  
                                    p_in, T0_in, qv0_in, ql0_in, qi0_in, & 
                                    a_dc_in, ql_dc_in, qi_dc_in,         &
                                    a_sc_in, ql_sc_in, qi_sc_in,         & 
                                    landfrac, snowh,                     &
                                    T_out, qv_out, ql_out, qi_out,       &
                                    al_st_out, ai_st_out, ql_st_out, qi_st_out, sgm_in)

   ! ------------------------------------------------------- !
   ! Diagnostically force in-stratus condensate to be        ! 
   ! in the range of 'qlst_min < qc_st < qlst_max'           !
   ! whenever stratus exists in the equilibrium state        !
   ! ------------------------------------------------------- !

   use time_manager,  only: is_first_step, get_nstep

   implicit none

   integer,  intent(in)  :: lchnk                ! Chunk identifier
   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index

   real(r8), intent(in)  :: p_in(pcols)          ! Pressure [Pa]
   real(r8), intent(in)  :: T0_in(pcols)         ! Temperature [K]
   real(r8), intent(in)  :: qv0_in(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0_in(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0_in(pcols)        ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)  :: a_dc_in(pcols)       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc_in(pcols)      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc_in(pcols)      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc_in(pcols)       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc_in(pcols)      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc_in(pcols)      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: landfrac(pcols)      ! Land fraction
   real(r8), intent(in)  :: snowh(pcols)         ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: T_out(pcols)         ! Temperature [K]
   real(r8), intent(out) :: qv_out(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(out) :: ql_out(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(out) :: qi_out(pcols)        ! Grid-mean IWC [kg/kg]

   real(r8), intent(out) :: al_st_out(pcols)     ! Liquid stratus fraction
   real(r8), intent(out) :: ai_st_out(pcols)     ! Ice stratus fraction
   real(r8), intent(out) :: ql_st_out(pcols)     ! In-stratus LWC [kg/kg]
   real(r8), intent(out) :: qi_st_out(pcols)     ! In-stratus IWC [kg/kg]
   ! qinyi
   real(r8), intent(in)  :: sgm_in(pcols)

   ! Local variables

   integer i                                     ! Column    index

   real(r8) p    
   real(r8) T0   
   real(r8) qv0    
   real(r8) ql0    
   real(r8) qi0    
   real(r8) a_dc   
   real(r8) ql_dc  
   real(r8) qi_dc  
   real(r8) a_sc   
   real(r8) ql_sc  
   real(r8) qi_sc  
   real(r8) esat0  
   real(r8) qsat0  
   real(r8) U0     
   real(r8) U0_nc  
   real(r8) G0_nc
   real(r8) al0_st_nc            
   real(r8) al0_st
   real(r8) ai0_st_nc            
   real(r8) ai0_st               
   real(r8) a0_st               
   real(r8) ql0_nc
   real(r8) qi0_nc
   real(r8) qc0_nc
   real(r8) ql0_st
   real(r8) qi0_st
   real(r8) qc0_st
   real(r8) T   
   real(r8) qv    
   real(r8) ql    
   real(r8) qi
   real(r8) ql_st
   real(r8) qi_st
   real(r8) es  
   real(r8) qs  

   real(r8) esat_in(pcols)  
   real(r8) qsat_in(pcols)  
   real(r8) U0_in(pcols)  
   real(r8) al0_st_nc_in(pcols)
   real(r8) ai0_st_nc_in(pcols)
   real(r8) G0_nc_in(pcols)

   real(r8) U
   real(r8) U_nc
   real(r8) al_st_nc
   real(r8) ai_st_nc
   real(r8) G_nc
   real(r8) a_st
   real(r8) al_st
   real(r8) ai_st
   real(r8) Tmin0
   real(r8) Tmax0
   real(r8) Tmin
   real(r8) Tmax

   ! qinyi
   real(r8) conden0_st_nc
   real(r8) conden0_st_nc_in(pcols)
   real(r8) conden_st_nc
   real(r8) sgm

   ! ---------------- !
   ! Main Computation ! 
   ! ---------------- !

   call qsat_water(T0_in(1:ncol), p_in(1:ncol), &
        esat_in(1:ncol), qsat_in(1:ncol))
   U0_in(:ncol) = qv0_in(:ncol)/qsat_in(:ncol)
   
   call cldfrac_conden(p_in(:),T0_in(:),qv0_in(:),ql0_in(:),sgm_in(:),al0_st_nc_in(:),conden0_st_nc_in(:),G0_nc_in(:),ncol)

   call aist_vector(qv0_in(:),T0_in(:),p_in(:),qi0_in(:),landfrac(:),snowh(:),ai0_st_nc_in(:),ncol)

   do i = 1, ncol

      ! ---------------------- !
      ! Define local variables !
      ! ---------------------- !

      p   = p_in(i)

      T0  = T0_in(i)
      qv0 = qv0_in(i)
      ql0 = ql0_in(i)
      qi0 = qi0_in(i)

      a_dc  = a_dc_in(i)
      ql_dc = ql_dc_in(i)
      qi_dc = qi_dc_in(i)

      a_sc  = a_sc_in(i)
      ql_sc = ql_sc_in(i)
      qi_sc = qi_sc_in(i)

      ql_dc = 0._r8
      qi_dc = 0._r8
      ql_sc = 0._r8
      qi_sc = 0._r8

      es  = esat_in(i) 
      qs  = qsat_in(i) 

      ! qinyi
      sgm           = sgm_in(i)
      conden0_st_nc = conden0_st_nc_in(i)

   ! -------------------------------------------------- !
   ! Force final energy-moisture conserving consistency !
   ! -------------------------------------------------- !

     qi = qi0

     ai_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
     al_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)

! qinyi 2016-05-21 16:43:52
     if( al_st .eq. 0._r8 ) then
         ql_st = 0._r8
     else
         ql_st = ql0/al_st
!         ql_st = min(qlst_max,max(qlst_min,ql_st)) ! PJR
     endif
     if( ai_st .eq. 0._r8 ) then
         qi_st = 0._r8
     else
         qi_st = qi0/ai_st
     endif

   ! -------------- !
   ! Send to output !
   ! -------------- !

   T_out(i)  = T0
   qv_out(i) = qv0
   ql_out(i) = ql0
   qi_out(i) = qi0
   al_st_out(i) = al_st
   ai_st_out(i) = ai_st
   ql_st_out(i) = ql_st
   qi_st_out(i) = qi_st

   enddo 

   return
   end subroutine instratus_condensate1

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice)

   ! ------------------------------------------------------------------------------- !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': top layer, 'pver' : near-surface layer         ! 
   ! ------------------------------------------------------------------------------- !

   implicit none
   integer,  intent(in)     :: ncol
   real(r8), intent(in)     :: dt
   real(r8), intent(in)     :: dp(pcols,pver), qvmin(pcols,pver), qlmin(pcols,pver), qimin(pcols,pver)
   real(r8), intent(inout)  :: qv(pcols,pver), ql(pcols,pver), qi(pcols,pver), t(pcols,pver)
   real(r8), intent(out)    :: qvten(pcols,pver), qlten(pcols,pver), qiten(pcols,pver), tten(pcols,pver)
   logical, intent(in)      :: do_cldice
   integer   i, k
   real(r8)  dql, dqi, dqv, sum, aa, dum 

   tten(:ncol,:pver)  = 0._r8
   qvten(:ncol,:pver) = 0._r8
   qlten(:ncol,:pver) = 0._r8
   qiten(:ncol,:pver) = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         if( qv(i,k) .lt. qvmin(i,k) .or. ql(i,k) .lt. qlmin(i,k) .or. qi(i,k) .lt. qimin(i,k) ) then
             goto 10
         endif
      enddo
      goto 11
   10 continue
      do k = top_lev, pver    ! From the top to the 1st (lowest) layer from the surface
         dql = max(0._r8,1._r8*qlmin(i,k)-ql(i,k))

         if (do_cldice) then
         dqi = max(0._r8,1._r8*qimin(i,k)-qi(i,k))
         else
           dqi = 0._r8
         end if

         qlten(i,k) = qlten(i,k) +  dql/dt
         qiten(i,k) = qiten(i,k) +  dqi/dt
         qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
         tten(i,k)  = tten(i,k)  + (latvap/cpair)*(dql/dt) + ((latvap+latice)/cpair)*(dqi/dt)
         ql(i,k)    = ql(i,k) + dql
         qi(i,k)    = qi(i,k) + dqi
         qv(i,k)    = qv(i,k) - dql - dqi
         t(i,k)     = t(i,k)  + (latvap * dql + (latvap+latice) * dqi)/cpair
         dqv        = max(0._r8,1._r8*qvmin(i,k)-qv(i,k))
         qvten(i,k) = qvten(i,k) + dqv/dt
         qv(i,k)    = qv(i,k)    + dqv
         if( k .ne. pver ) then 
             qv(i,k+1)    = qv(i,k+1)    - dqv*dp(i,k)/dp(i,k+1)
             qvten(i,k+1) = qvten(i,k+1) - dqv*dp(i,k)/dp(i,k+1)/dt
         endif
         qv(i,k) = max(qv(i,k),qvmin(i,k))
         ql(i,k) = max(ql(i,k),qlmin(i,k))
         qi(i,k) = max(qi(i,k),qimin(i,k))
      end do
      ! Extra moisture used to satisfy 'qv(i,pver)=qvmin' is proportionally 
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture. 
      if( dqv .gt. 1.e-20_r8 ) then
          sum = 0._r8
          do k = top_lev, pver
             if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) sum = sum + qv(i,k)*dp(i,k)
          enddo
          aa = dqv*dp(i,pver)/max(1.e-20_r8,sum)
          if( aa .lt. 0.5_r8 ) then
              do k = top_lev, pver
                 if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) then
                     dum        = aa*qv(i,k)
                     qv(i,k)    = qv(i,k) - dum
                     qvten(i,k) = qvten(i,k) - dum/dt
                 endif
              enddo 
          else 
              write(iulog,*) 'Full positive_moisture is impossible in Park Macro'
          endif
      endif 
11 continue
   enddo
   return

   end subroutine positive_moisture

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine astG_PDF( U_in, p_in, qv_in, landfrac_in, snowh_in, a_out, Ga_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute 'stratus fraction(a)' and Gs=(dU/da) from the     !
   ! analytical formulation of triangular PDF.                 !
   ! Here, 'dV' is the ratio of 'half-width of PDF / qs(p,T)', !
   ! so using constant 'dV' assume that width is proportional  !
   ! to the saturation specific humidity.                      !
   !    dV ~ 0.1.                                              !
   !    cldrh : RH of in-stratus( = 1 if no supersaturation)   !
   ! Note that if U > 1, Ga = 1.e10 instead of Ga = 0, that is !
   ! G is discontinuous across U = 1.  In fact, it does not    !
   ! matter whether Ga = 1.e10 or 0 at a = 1: I derived that   !
   ! they will produce the same results.                       !
   ! --------------------------------------------------------- !

   implicit none

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: U_in(pcols)           ! Relative humidity
   real(r8), intent(in)  :: p_in(pcols)           ! Pressure [Pa]
   real(r8), intent(in)  :: qv_in(pcols)          ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols)    ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)       ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: a_out(pcols)          ! Stratus fraction
   real(r8), intent(out) :: Ga_out(pcols)         ! dU/da

   real(r8)              :: U                     ! Relative humidity
   real(r8)              :: p                     ! Pressure [Pa]
   real(r8)              :: qv                    ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8)              :: landfrac              ! Land fraction
   real(r8)              :: snowh                 ! Snow depth (liquid water equivalent)

   real(r8)              :: a                     ! Stratus fraction
   real(r8)              :: Ga                    ! dU/da

   ! Local variables
   integer :: i                                   ! Loop indexes
   real(r8) dV                                    ! Width of triangular PDF
   real(r8) cldrh                                 ! RH of stratus cloud
   real(r8) rhmin                                 ! Critical RH
   real(r8) rhwght
                            
   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! ---------- !
   ! Parameters !
   ! ---------- !

   cldrh  = 1.0_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   a_out(:)  = 0._r8
   Ga_out(:) = 0._r8

   do i = 1, ncol

   U        = U_in(i)      
   p        = p_in(i)        
   qv       = qv_in(i)       
   landfrac = landfrac_in(i) 
   snowh    = snowh_in(i)    

   if( p .ge. premib ) then

       if( land(i) .and. (snowh.le.0.000001_r8) ) then
           rhmin = rhminl - rhminl_adj_land
       else
           rhmin = rhminl
       endif

       dV = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

       if( freeze_dry ) then
           a  = a *max(0.15_r8,min(1.0_r8,qv/0.0030_r8))
           Ga = Ga/max(0.15_r8,min(1.0_r8,qv/0.0030_r8)) 
       endif

   elseif( p .lt. premit ) then

       rhmin = rhminh
       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                      (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   else

       rhwght = (premib-(max(p,premit)))/(premib-premit)

     ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
     !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
     ! else
           rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
     ! endif

       dV    = cldrh - rhmin

       if( U .ge. 1._r8 ) then
           a  = 1._r8
           Ga = 1.e10_r8
       elseif( U .gt. (cldrh-dV/6._r8) .and. U .lt. 1._r8 ) then
           a  = 1._r8 - (-3._r8/sqrt(2._r8)*(U-cldrh)/dV)**(2._r8/3._r8)
           Ga = dV/sqrt(2._r8)*sqrt(1._r8-a)
       elseif( U .gt. (cldrh-dV) .and. U .le. (cldrh-dV/6._r8) ) then
           a  = 4._r8*(cos((1._r8/3._r8)*(acos((3._r8/2._r8/sqrt(2._r8))* & 
                         (1._r8+(U-cldrh)/dV))-2._r8*3.141592_r8)))**2._r8
           Ga = dV/sqrt(2._r8)*(1._r8/sqrt(a)-sqrt(a))
       elseif( U .le. (cldrh-dV) ) then
           a  = 0._r8
           Ga = 1.e10_r8
       endif

   endif

   a_out(i)  = a
   Ga_out(i) = Ga 

   enddo

   return
   end subroutine astG_PDF

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine aist_single( qv, T, p, qi, landfrac, snowh, aist )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   use physconst,     only: rair

   implicit none
  
   real(r8), intent(in)  :: qv              ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T               ! Temperature
   real(r8), intent(in)  :: p               ! Pressure [Pa]
   real(r8), intent(in)  :: qi              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac        ! Land fraction
   real(r8), intent(in)  :: snowh           ! Snow depth (liquid water equivalent)

   real(r8), intent(out) :: aist            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables
   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) es, qs

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme


   ! Statement functions
   logical land
   land = nint(landfrac) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8 
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     call qsat_water(T, p, es, qs)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

   return
   end subroutine aist_single

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine aist_vector( qv_in, T_in, p_in, qi_in, landfrac_in, snowh_in, aist_out, ncol )

   ! --------------------------------------------------------- !
   ! Compute non-physical ice stratus fraction                 ! 
   ! --------------------------------------------------------- !

   use physconst,     only: rair

   implicit none
  
   integer,  intent(in)  :: ncol 
   real(r8), intent(in)  :: qv_in(pcols)       ! Grid-mean water vapor[kg/kg]
   real(r8), intent(in)  :: T_in(pcols)        ! Temperature
   real(r8), intent(in)  :: p_in(pcols)        ! Pressure [Pa]
   real(r8), intent(in)  :: qi_in(pcols)       ! Grid-mean ice water content [kg/kg]
   real(r8), intent(in)  :: landfrac_in(pcols) ! Land fraction
   real(r8), intent(in)  :: snowh_in(pcols)    ! Snow depth (liquid water equivalent)
   real(r8), intent(out) :: aist_out(pcols)    ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   ! Local variables

   real(r8) qv                              ! Grid-mean water vapor[kg/kg]
   real(r8) T                               ! Temperature
   real(r8) p                               ! Pressure [Pa]
   real(r8) qi                              ! Grid-mean ice water content [kg/kg]
   real(r8) landfrac                        ! Land fraction
   real(r8) snowh                           ! Snow depth (liquid water equivalent)
   real(r8) aist                            ! Non-physical ice stratus fraction ( 0<= aist <= 1 )

   real(r8) rhmin                           ! Critical RH
   real(r8) rhwght

   real(r8) a,b,c,as,bs,cs                  ! Fit parameters
   real(r8) Kc                              ! Constant for ice cloud calc (wood & field)
   real(r8) ttmp                            ! Limited temperature
   real(r8) icicval                         ! Empirical IWC value [ kg/kg ]
   real(r8) rho                             ! Local air density
   real(r8) esl                             ! Liq sat vapor pressure
   real(r8) esi                             ! Ice sat vapor pressure
   real(r8) ncf,phi                         ! Wilson and Ballard parameters
   real(r8) qs
   real(r8) esat_in(pcols)
   real(r8) qsat_in(pcols)

   real(r8) rhi                             ! grid box averaged relative humidity over ice
   real(r8) minice                          ! minimum grid box avg ice for having a 'cloud'
   real(r8) mincld                          ! minimum ice cloud fraction threshold
   real(r8) icimr                           ! in cloud ice mixing ratio
 ! real(r8) qist_min                        ! minimum in cloud ice mixing ratio
 ! real(r8) qist_max                        ! maximum in cloud ice mixing ratio                
   real(r8) rhdif                           ! working variable for slingo scheme

   integer i


   ! Statement functions
   logical land
   land(i) = nint(landfrac_in(i)) == 1

   ! --------- !
   ! Constants !
   ! --------- !

   ! Wang and Sassen IWC paramters ( Option.1 )
     a = 26.87_r8
     b = 0.569_r8
     c = 0.002892_r8
   ! Schiller parameters ( Option.2 )
     as = -68.4202_r8
     bs = 0.983917_r8
     cs = 2.81795_r8
   ! Wood and Field parameters ( Option.3 )
     Kc = 75._r8
   ! Wilson & Ballard closure ( Option.4. smaller = more ice clouds)
   ! Slingo modified (option 5)
     minice = 1.e-12_r8
     mincld = 1.e-4_r8
   ! qist_min = 1.e-7_r8
   ! qist_max = 5.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

     aist_out(:) = 0._r8
     esat_in(:)  = 0._r8
     qsat_in(:)  = 0._r8

     call qsat_water(T_in(1:ncol), p_in(1:ncol), &
          esat_in(1:ncol), qsat_in(1:ncol))
     
     do i = 1, ncol

     landfrac = landfrac_in(i)     
     snowh = snowh_in(i)   
     T = T_in(i)
     qv = qv_in(i)
     p = p_in(i)
     qi = qi_in(i)
     qs = qsat_in(i)
     esl = svp_water(T)
     esi = svp_ice(T)
          
     if( iceopt.lt.3 ) then
         if( iceopt.eq.1 ) then
             ttmp = max(195._r8,min(T,253._r8)) - 273.16_r8
             icicval = a + b * ttmp + c * ttmp**2._r8
             rho = p/(rair*T)
             icicval = icicval * 1.e-6_r8 / rho 
         else
             ttmp = max(190._r8,min(T,273.16_r8))
             icicval = 10._r8 **(as * bs**ttmp + cs)
             icicval = icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
         endif
         aist =  max(0._r8,min(qi/icicval,1._r8)) 
     elseif( iceopt.eq.3 ) then
         aist = 1._r8 - exp(-Kc*qi/(qs*(esi/esl)))
         aist = max(0._r8,min(aist,1._r8))
     elseif( iceopt.eq.4) then
         if( p .ge. premib ) then
             if( land(i) .and. (snowh.le.0.000001_r8) ) then
                 rhmin = rhminl - rhminl_adj_land
             else
                 rhmin = rhminl
             endif
         elseif( p .lt. premit ) then
             rhmin = rhminh
         else
             rhwght = (premib-(max(p,premit)))/(premib-premit)
           ! if( land(i) .and. (snowh.le.0.000001_r8) ) then
           !     rhmin = rhminh*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
           ! else
                 rhmin = rhminh*rhwght + rhminl*(1.0_r8-rhwght)
           ! endif
         endif
         ncf = qi/((1._r8 - icecrit)*qs)
         if( ncf.le.0._r8 ) then 
             aist = 0._r8
         elseif( ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8 ) then 
             aist = 0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
         elseif( ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8 ) then
             phi = (acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             aist = (1._r8 - 4._r8 * cos(phi) * cos(phi))
         else
             aist = 1._r8
         endif
             aist = max(0._r8,min(aist,1._r8))
     elseif (iceopt.eq.5) then 
! set rh ice cloud fraction
             rhi= (qv+qi)/qs * (esl/esi)
             rhdif= (rhi-rhmini) / (rhmaxi - rhmini)
             aist = min(1.0_r8, max(rhdif,0._r8)**2)

! limiter to remove empty cloud and ice with no cloud
! and set icecld fraction to mincld if ice exists

             if (qi.lt.minice) then
                aist=0._r8
             else
                aist=max(mincld,aist)
             endif

! enforce limits on icimr
             if (qi.ge.minice) then
                icimr=qi/aist

!minimum
                if (icimr.lt.qist_min) then
                   aist = max(0._r8,min(1._r8,qi/qist_min))
                endif
!maximum
                if (icimr.gt.qist_max) then
                   aist = max(0._r8,min(1._r8,qi/qist_max))
                endif

             endif
     endif 

   ! 0.999_r8 is added to prevent infinite 'ql_st' at the end of instratus_condensate
   ! computed after updating 'qi_st'.  

     aist = max(0._r8,min(aist,0.999_r8))

     aist_out(i) = aist

     enddo

   return
   end subroutine aist_vector

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module cldwat2m_macro
