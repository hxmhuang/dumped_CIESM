module cldfrac_conden_mod

use wv_saturation, only: qsat_water
use physconst, only: cpair, latvap, latice, rh2o, gravit,rair
use shr_kind_mod, only: r8=>shr_kind_r8
use ppgrid, only: pcols,pverp,pver
use cam_history, only: outfld,addfld,add_default,phys_decomp
use spmd_utils,      only: masterproc
use abortutils,      only: endrun

! the unit of mass flux from deep convection is mb/s

implicit none
private
save

public :: cldfrac_conden_init,subgrid_var,cldfrac_conden,cldfrac_conden_single,cldfrac_conden_single_only



! qinyi 2017-03-02 16:00:28
real(r8), parameter :: unset_r8 = huge(1.0_r8)

real(r8) :: cldfrac_conden_rcon2 = unset_r8
real(r8) :: cldfrac_conden_rcon3 = unset_r8

real(r8) :: rcon2
real(r8) :: rcon3


contains

subroutine cldfrac_conden_init()

use cam_history, only: outfld,addfld,add_default,phys_decomp
use ppgrid, only: pver, pverp,pcols

!shallow
call addfld('qtu_shal','kg/kg',pverp,'A','updraft specific humidity of shallow conv',phys_decomp)
call addfld('thlu_shal','K',pverp,'A','updraft liquid potential temperature  of shallow conv',phys_decomp)
call addfld('umf_shal','(kg/kg)*(m/s)',pverp,'A','shal convective mass flux',phys_decomp)
call addfld('cnt_shal','no',1,'A','shallow convection top index',phys_decomp)
call addfld('cnb_shal','no',1,'A','shallow convection bot index',phys_decomp)
call addfld('clddep2','m',1,'A','shal cloud depth',phys_decomp)
!turbulent 
call addfld('wstarPBL','m/s',1,'A','convective velocity scale',phys_decomp)
call addfld('lengi','m',pverp,'A','turbulent length scale',phys_decomp)
call addfld('shi','no',pverp,'A','stability function',phys_decomp)

call addfld('beta_gp','no',pver,'A','coefficient beta',phys_decomp)
call addfld('aa_gp','no',pver,'A','coefficient a',phys_decomp)
call addfld('bb_gp','no',pver,'A','coefficient b',phys_decomp)
call addfld('dqwdz','kg/kg/m',pver,'A','vertial gradient of qw',phys_decomp)
call addfld('dthldz','K/m',pver,'A','vertical gradient of thl',phys_decomp)

call addfld('qw2','(kg/kg)^2',pver,'A','qwprime^2',phys_decomp)
call addfld('thl2','(K)^2',pver,'A','thlprime^2',phys_decomp)
call addfld('qwthl','(kg/kg)*(K)',pver,'A','qwprime*thlprime',phys_decomp)

call addfld('sgm_turb','no',pver,'A','standard deviation by turbulence',phys_decomp)
call addfld('sgm_shal','no',pver,'A','standard deviation by shallow convection',phys_decomp)
call addfld('sgm_tota','no',pver,'A','standard deviation by all processes',phys_decomp)


end subroutine cldfrac_conden_init

subroutine subgrid_var(lchnk,ncol,rpdel_in,p_in,z_in,T_in,qv_in,ql_in, &
                       tke_in,kvh_in,kvm_in,lengi_in,shi_in,smi_in,wstarPBL_in,pblh_in,&
                       qtu_shal, umf_shal, cnt_shal, cnb_shal, thlu_shal, &
                       explwi2qtu, & ! qinyi 2017-11-30 21:49:36
                       sgm_out)

use cam_history, only: outfld,addfld,add_default,phys_decomp
use ppgrid, only: pcols,pverp,pver
use physconst, only: cpair, latvap, latice, rh2o, gravit,rair

implicit none

integer,intent(in) :: lchnk
integer,intent(in) :: ncol
real(r8),intent(in) :: p_in(pcols,pver) ! pressure
real(r8),intent(in) :: z_in(pcols,pver) ! height
real(r8),intent(in) :: rpdel_in(pcols,pver) ! pressure difference
real(r8),intent(in) :: T_in(pcols,pver) ! temperature
real(r8),intent(in) :: qv_in(pcols,pver)! water vapor mixing ratio
real(r8),intent(in) :: ql_in(pcols,pver) ! liquid water mixing ratio
real(r8),intent(in) :: tke_in(pcols,pverp) ! turbulent kinetic energy
real(r8),intent(in) :: kvh_in(pcols,pverp) ! turbulent diffusivity -- heat
real(r8),intent(in) :: kvm_in(pcols,pverp) ! turbulent diffusivity -- momentum
real(r8),intent(in) :: lengi_in(pcols,pverp) ! master turbulent length scale
real(r8),intent(in) :: shi_in(pcols,pverp) ! stability function -- heat
real(r8),intent(in) :: smi_in(pcols,pverp) ! stability function -- momentum
real(r8),intent(in) :: wstarPBL_in(pcols)  ! convective vertical velocity scale
! 2017-01-14 19:16:01
real(r8),intent(in) :: pblh_in(pcols) ! PBL height

!shallow
real(r8),intent(in) :: qtu_shal(pcols,pverp) ! total specific humidity in updraft
real(r8),intent(in) :: thlu_shal(pcols,pverp) ! total liquid poten. temp. in updraft
real(r8),intent(in) :: umf_shal(pcols,pverp) ! shallow mass flux 
real(r8),intent(in) :: cnt_shal(pcols) ! cloud top
real(r8),intent(in) :: cnb_shal(pcols) ! cloud base 
! qinyi 2017-11-30 21:49:49
real(r8),intent(in) :: explwi2qtu(pcols,pver) ! PE from shallow convection
! qinyi

real(r8),intent(out) :: sgm_out(pcols,pver) ! total subgrid scale variance


real(r8),parameter :: c1_gp = 0.32_r8 ! c1=2*ck/cab
real(r8),parameter :: kafa_gp = 7.44 ! consant for calculation of ld
real(r8),parameter :: p00_gp = 100000._r8 ! reference surface pressure
real(r8),parameter :: lkfix_gp=50._r8 ! assumed fixed turbulent length
real(r8),parameter :: rcon22 = 7.14_r8 !coefficient
real(r8),parameter :: rcon33 = 1.0_r8 ! coefficient

! qinyi 2017-12-11 15:34:39
! introduce PE factor to further amplify the effect from shallow precipitation's
! effect.
!qy real(r8),parameter :: pe_factor = 0.4_r8
real(r8),parameter :: pe_factor = 1.0_r8

!local variables
integer :: i,j,k

real(r8) :: exner_gp(pcols,pver)
real(r8) :: tl_gp(pcols,pver) ! liquid water temperature
real(r8) :: thl_gp(pcols,pver) ! liquid potential temp.
real(r8) :: qw_gp(pcols,pver) ! total specific humidity
real(r8) :: rho_gp(pcols,pver) ! density


real(r8) :: a_gp(pcols,pver) ! coefficent "a"
real(r8) :: b_gp(pcols,pver)! coefficient "b"
real(r8) :: beta_gp(pcols,pver) ! coefficient "beta"
real(r8) :: qs_gp(pcols,pver) ! saturated specific humidity
real(r8) :: es_gp(pcols,pver) ! satureated water vapor pressure
real(r8) :: qstl_gp(pcols,pver) ! saturation specific humidity for liquid water temp.
real(r8) :: estl_gp(pcols,pver) ! saturation pressure for liquid water temp.

real(r8) :: dqwdz(pcols,pver) ! vertical gradient of total water
real(r8) :: dthldz(pcols,pver) ! vertical gradient of liquid potential temperature

real(r8) :: cld_dep2(pcols) ! cloud depth from shallow convection

real(r8) :: qw2(pcols,pver) ! qw*qw
real(r8) :: thl2(pcols,pver) ! thl*thl
real(r8) :: qwthl(pcols,pver) ! qu*thl

integer :: cnt_int(pcols) ! cloud top
integer :: cnb_int(pcols) ! cloud base

real(r8) :: sgm_s1(pcols,pver) ! variance from turbulent
real(r8) :: sgm_s2(pcols,pver) ! variance from shallow convection

!!!!!!!!!!!!
! intialize
!!!!!!!!!!!!

! output variable
sgm_out(:ncol,:pver) = 0._r8

! local variables
exner_gp(:ncol,:pver) = 0._r8
tl_gp(:ncol,:pver) = 0._r8
thl_gp(:ncol,:pver) = 0._r8
qw_gp(:ncol,:pver) = 0._r8
rho_gp(:ncol,:pver) = 0._r8

a_gp(:ncol,:pver) = 0._r8
b_gp(:ncol,:pver)= 0._r8
beta_gp(:ncol,:pver)= 0._r8
qs_gp(:ncol,:pver)= 0._r8
es_gp(:ncol,:pver)= 0._r8
qstl_gp(:ncol,:pver)= 0._r8
estl_gp(:ncol,:pver)= 0._r8

dqwdz(:ncol,:pver) = 0._r8
dthldz(:ncol,:pver) = 0._r8

qw2(:ncol,:pver) = 0._r8
thl2(:ncol,:pver) = 0._r8
qwthl(:ncol,:pver) = 0._r8

cnt_int(:ncol) = int(cnt_shal(:ncol))
cnb_int(:ncol) = int(cnb_shal(:ncol))

sgm_s1(:ncol,:pver) = 0._r8
sgm_s2(:ncol,:pver) = 0._r8

!!!!!!!!!!!!!!!!!!!!!!
! calculation starts.
!!!!!!!!!!!!!!!!!!!!!!

do i=1,ncol
do k=1,pver

exner_gp(i,k) = (p_in(i,k)/p00_gp)**0.286_r8
qw_gp(i,k) = qv_in(i,k)+ql_in(i,k)
tl_gp(i,k) = T_in(i,k)-latvap/cpair*ql_in(i,k)
thl_gp(i,k) = tl_gp(i,k)/exner_gp(i,k)
rho_gp(i,k) = p_in(i,k)/(rair*T_in(i,k))


call qsat_water(tl_gp(i,k),p_in(i,k),estl_gp(i,k),qstl_gp(i,k))
call qsat_water(t_in(i,k),p_in(i,k),es_gp(i,k),qs_gp(i,k))

beta_gp(i,k) = 0.622_r8*latvap*qstl_gp(i,k)/rair/tl_gp(i,k)**2

a_gp(i,k) = 1._r8/(1._r8+beta_gp(i,k)*latvap/cpair)
b_gp(i,k) = exner_gp(i,k)*beta_gp(i,k)/(1._r8+beta_gp(i,k)*latvap/cpair)

enddo !i=1,ncol
enddo !k=1,pver

!================================================================!
!=========================turbulent effect=======================!
do i=1,ncol
        do k = 1,pver-1
        ! part I: caused by turbulent process
        !dqwdz(i,k) = -rho(i,k)*gravit*(qw(i,k+1)-qw(i,k))*rpdel_in(i,k)
        !dthldz(i,k) =-rho(i,k)*gravit*(thl(i,k+1)-thl(i,k))*rpdel_in(i,k) 
        
        dqwdz(i,k) = (qw_gp(i,k+1)-qw_gp(i,k))/(z_in(i,k+1)-z_in(i,k))
        dthldz(i,k) = (thl_gp(i,k+1)-thl_gp(i,k))/(z_in(i,k+1)-z_in(i,k))
        enddo

        dqwdz(i,pver) = 0._r8
        dthldz(i,pver) = 0._r8
enddo


do i=1,ncol
        do k=1,pver
        qw2(i,k) = 2.0_r8*rcon22*lengi_in(i,k+1)*lengi_in(i,k+1)*shi_in(i,k+1)*dqwdz(i,k)**2_r8
        thl2(i,k) = 2.0_r8*rcon22*lengi_in(i,k+1)*lengi_in(i,k+1)*shi_in(i,k+1)*dthldz(i,k)**2_r8
        qwthl(i,k) = 2.0_r8*rcon22*lengi_in(i,k+1)*lengi_in(i,k+1)*shi_in(i,k+1)*dqwdz(i,k)*dthldz(i,k)
        
        sgm_s1(i,k) = 0.25_r8*2.0_r8*rcon22*lengi_in(i,k+1)*lengi_in(i,k+1)*shi_in(i,k+1)*(a_gp(i,k)**2*dqwdz(i,k)**2 + &
                      b_gp(i,k)**2*dthldz(i,k)**2 - 2._r8*a_gp(i,k)*b_gp(i,k)*dqwdz(i,k)*dthldz(i,k))
        enddo
enddo 
!=============================================================!
!=====================shallow convection effect===============!

do i=1,ncol
do k=1,pver

! cloud depth
cld_dep2(i) = z_in(i,cnt_int(i))-z_in(i,cnb_int(i)) 

if(cld_dep2(i).le.0)then
        cld_dep2(i) =0._r8
endif

if((wstarPBL_in(i).ne.0._r8).and.(qtu_shal(i,k).ne.0._r8).and.(thlu_shal(i,k).ne.0._r8))then
    sgm_s2(i,k) = 0.5_r8*rcon33*umf_shal(i,k)/rho_gp(i,k)*cld_dep2(i)/wstarPBL_in(i)&
                 *(a_gp(i,k)**2*(qtu_shal(i,k)-qw_gp(i,k))*dqwdz(i,k)&
                  +b_gp(i,k)**2*(thlu_shal(i,k)-thl_gp(i,k))*dthldz(i,k) &
                  -a_gp(i,k)*b_gp(i,k)*((thlu_shal(i,k)-thl_gp(i,k))*dqwdz(i,k)+(qtu_shal(i,k)-qw_gp(i,k))*dthldz(i,k)))

! qinyi 2017-11-30 21:52:23
! consider precipitation's effect on subgrid scale variance
    if(explwi2qtu(i,k).gt.1e-4_r8)then
!      sgm_s2(i,k) = sgm_s2(i,k)*explwi2qtu(i,k)
        ! qinyi 2017-12-11 15:35:57
      sgm_s2(i,k) = sgm_s2(i,k)*explwi2qtu(i,k)*pe_factor
    else
      sgm_s2(i,k) = sgm_s2(i,k)
    endif
endif
! qinyi

enddo
enddo

! calcuate the total subgrid scale variance
do k=1,pver
        do i=1,ncol
        sgm_out(i,k) = abs(sgm_s1(i,k))+abs(sgm_s2(i,k))
        enddo
enddo

!shallow
call outfld('qtu_shal',qtu_shal,pcols,lchnk)
call outfld('thlu_shal',thlu_shal,pcols,lchnk)
call outfld('umf_shal',umf_shal,pcols,lchnk)
call outfld('cnt_shal',cnt_shal,pcols,lchnk)
call outfld('cnb_shal',cnb_shal,pcols,lchnk)
call outfld('clddep2',cld_dep2,pcols,lchnk)
!turbulent
call outfld('wstarPBL',wstarPBL_in,pcols,lchnk)
call outfld('lengi',lengi_in,pcols,lchnk)
call outfld('shi',shi_in,pcols,lchnk)

call outfld('beta_gp',beta_gp,pcols,lchnk)
call outfld('aa_gp',a_gp,pcols,lchnk)
call outfld('bb_gp',b_gp,pcols,lchnk)
call outfld('dqwdz',dqwdz,pcols,lchnk)
call outfld('dthldz',dthldz,pcols,lchnk)

call outfld('qw2',qw2,pcols,lchnk)
call outfld('thl2',thl2,pcols,lchnk)
call outfld('qwthl',qwthl,pcols,lchnk)

call outfld('sgm_turb',sgm_s1,pcols,lchnk)
call outfld('sgm_shal',sgm_s2,pcols,lchnk)
call outfld('sgm_tota',sgm_out,pcols,lchnk)

end subroutine subgrid_var


!===================================================
subroutine cldfrac_conden(p_in,T_in,qv_in,ql_in,sgm_in,cld_frac_out,conden_out,G,ncol)

implicit none

integer, intent(in) :: ncol
real(r8),intent(in) :: p_in(pcols) ! pressure
real(r8),intent(in) :: T_in(pcols) ! temperature
real(r8),intent(in) :: qv_in(pcols)! water vapor specific humidity
real(r8),intent(in) :: ql_in(pcols) ! liquid water specific humidity
real(r8),intent(in) :: sgm_in(pcols) ! total subgrid scale variance

real(r8),intent(out) :: cld_frac_out(pcols) ! diagnosed cloud fraction from Gauss-PDF
real(r8),intent(out) :: conden_out(pcols) ! diagnosed cloud condensate
real(r8),intent(out) :: G(pcols) ! nothing.

real(r8),parameter:: pi = 3.14159265358_r8
real(r8),parameter:: sgm_min = 1.e-8_r8
real(r8),parameter:: p00 = 100000._r8

!local variables
integer :: i,j,k

real(r8) :: conden(pcols) !ql/(2*sigma(s))
real(r8) :: cld_frac(pcols)

real(r8) :: exner_gp(pcols)
real(r8) :: qw_gp(pcols)
real(r8) :: tl_gp(pcols)
real(r8) :: thl_gp(pcols)
real(r8) :: rho_gp(pcols)

real(r8) :: qs_gp
real(r8) :: es_gp
real(r8) :: qstl_gp ! saturation specific humidity for liquid water temp.
real(r8) :: estl_gp ! saturation pressure for liquid water temp.

real(r8) :: a_gp
real(r8) :: b_gp
real(r8) :: beta_gp
real(r8) :: deltaq_gp
real(r8) :: Q1_gp


exner_gp(:ncol) = 0._r8
qw_gp(:ncol) = 0._r8
tl_gp(:ncol)= 0._r8
thl_gp(:ncol) = 0._r8
rho_gp(:ncol) = 0._r8

conden(:ncol) = 0._r8
cld_frac(:ncol) = 0._r8


do k=1,ncol

qs_gp = 0._r8
es_gp = 0._r8
estl_gp = 0._r8
qstl_gp = 0._r8

a_gp = 0._r8
b_gp = 0._r8
beta_gp = 0._r8
deltaq_gp = 0._r8
Q1_gp = 0._r8


exner_gp(k) = (p_in(k)/p00)**0.286_r8
qw_gp(k) = qv_in(k)+ql_in(k)
tl_gp(k) = T_in(k)-latvap/cpair*ql_in(k)
thl_gp(k) = tl_gp(k)/exner_gp(k)
rho_gp(k) = p_in(k)/(rair*T_in(k))

call qsat_water(tl_gp(k),p_in(k),estl_gp,qstl_gp)
call qsat_water(t_in(k),p_in(k),es_gp,qs_gp)

beta_gp = 0.622_r8*latvap*qstl_gp/rair/tl_gp(k)**2
a_gp = 1._r8/(1._r8+beta_gp*latvap/cpair)
b_gp = exner_gp(k)*beta_gp/(1._r8+beta_gp*latvap/cpair)

deltaq_gp = qw_gp(k)-qstl_gp

if(sgm_in(k).ne.0._r8)then
        Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in(k)))
else
        Q1_gp = 0._r8
endif

! calculate cloud fraction and condensate
if(sqrt(sgm_in(k))>sgm_min)then
    cld_frac(k) = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
    conden(k) = cld_frac(k)*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)

    cld_frac_out(k) = min(1._r8,max(0._r8,cld_frac(k)))
    conden_out(k) = conden(k)*2._r8*sqrt(sgm_in(k))
    conden_out(k) = max(0._r8,conden_out(k))

else
    if(deltaq_gp.le.0._r8)then
      cld_frac(k) = 0._r8
      conden(k) = 0._r8

      cld_frac_out(k) = min(1._r8,max(0._r8,cld_frac(k)))
      conden_out(k) = conden(k)*2._r8*sqrt(sgm_in(k))
      conden_out(k) = max(0._r8,conden_out(k))

    else
      cld_frac(k) = 1._r8
      conden_out(k) = a_gp*deltaq_gp

      cld_frac_out(k) = min(1._r8,max(0._r8,cld_frac(k)))
      conden_out(k) = max(0._r8,conden_out(k))

    endif
endif

! qinyi 2016-11-18 21:59:05
! test by simple calculation indicates that there is unstable oscillation while
! cldfrac is less than 1e-4.
! the exact reason for this is the too small deltaq could cause this.
! so in order to eliminate this inapproximate value, I set the threshold for
! cloud fraction.

if((conden_out(k).eq.0._r8).or.(cld_frac_out(k).lt.1.e-4_r8))then
        conden_out(k) = 0._r8
        cld_frac_out(k) = 0._r8
endif


G(k) = 0._r8

enddo !k=1,pcols

end subroutine cldfrac_conden

!===================================================
subroutine cldfrac_conden_single_only(p_in,T_in,qv_in,ql_in,sgm_in,cld_frac_out,conden_out,G,deltaq_sat,deltaq_uns,Q1_sat,Q1_uns,adjust_factor)

real(r8),intent(in) :: p_in
real(r8),intent(in) :: T_in
real(r8),intent(in) :: qv_in
real(r8),intent(in) :: ql_in
real(r8),intent(in) :: sgm_in

real(r8),intent(out) :: cld_frac_out
real(r8),intent(out) :: conden_out
real(r8),intent(out) :: G
real(r8),intent(out) :: deltaq_sat
real(r8),intent(out) :: deltaq_uns
real(r8),intent(out) :: Q1_sat
real(r8),intent(out) :: Q1_uns
real(r8),intent(out) :: adjust_factor

real(r8),parameter:: pi = 3.14159265358_r8
real(r8),parameter:: sgm_min = 1.e-8_r8
real(r8),parameter:: p00 = 100000._r8

!local variables
integer :: i,j,k

real(r8) :: conden !ql/(2*sigma(s))
real(r8) :: cld_frac

real(r8) :: exner_gp
real(r8) :: tl_gp
real(r8) :: thl_gp
real(r8) :: rho_gp
real(r8) :: qw_gp
real(r8) :: theta_gp

real(r8) :: qs_gp
real(r8) :: es_gp
real(r8) :: qstl_gp
real(r8) :: estl_gp

real(r8) :: a_gp
real(r8) :: b_gp
real(r8) :: beta_gp
real(r8) :: deltaq_gp
real(r8) :: Q1_gp

real(r8) :: adjust_mid

! initialized
! output variables
cld_frac_out = 0._r8
conden_out = 0._r8
G = 0._r8
deltaq_sat = 0._r8
deltaq_uns = 0._r8
Q1_sat = 0._r8
Q1_uns = 0._r8
adjust_factor = 0._r8

! local variables
conden = 0._r8
cld_frac = 0._r8

exner_gp = 0._r8
tl_gp = 0._r8
thl_gp = 0._r8
rho_gp = 0._r8
qw_gp = 0._r8
theta_gp = 0._r8

qs_gp = 0._r8
es_gp = 0._r8
qstl_gp = 0._r8
estl_gp = 0._r8

a_gp = 0._r8
b_gp = 0._r8
beta_gp = 0._r8
deltaq_gp = 0._r8
Q1_gp = 0._r8

! qinyi 2017-1-1 21:53:57
adjust_mid = 0._r8

! calculation starts

exner_gp = (p_in/p00)**0.286_r8
tl_gp = T_in-latvap/cpair*ql_in
thl_gp = tl_gp/exner_gp
theta_gp = T_in/exner_gp
rho_gp = p_in/(rair*T_in)
qw_gp = qv_in+ql_in

call qsat_water(tl_gp,p_in,estl_gp,qstl_gp)
call qsat_water(t_in,p_in,es_gp,qs_gp)

beta_gp = 0.622_r8*latvap*qstl_gp/rair/tl_gp**2

a_gp = 1._r8/(1._r8+beta_gp*latvap/cpair)
b_gp = exner_gp*beta_gp/(1._r8+beta_gp*latvap/cpair)


deltaq_gp = qw_gp-qstl_gp

if(sgm_in.ne.0_r8)then
        Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in))
else
        Q1_gp = 0._r8
endif

! calculate cloud fraction and condensate
if(sqrt(sgm_in)>sgm_min)then
    cld_frac = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
    conden = cld_frac*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)
    
    cld_frac_out = min(1._r8,max(0._r8,cld_frac))
    conden_out = conden*2._r8*sqrt(sgm_in)
    conden_out = max(0._r8,conden_out)

else
    if(deltaq_gp.le.0._r8)then
      cld_frac = 0._r8
      conden = 0._r8
      cld_frac_out = min(1._r8,max(0._r8,cld_frac))
      conden_out = conden*2._r8*sqrt(sgm_in)
      conden_out = max(0._r8,conden_out)

    else
      cld_frac = 1._r8
      conden_out = a_gp*deltaq_gp
      cld_frac_out = min(1._r8,max(0._r8,cld_frac))
      conden_out = max(0._r8,conden_out)

    endif
endif
! qinyi 2016-11-18 22:01:15
if((conden_out.eq.0._r8).or.(cld_frac_out.lt.1.e-4_r8))then
        conden_out = 0._r8
        cld_frac_out = 0._r8
endif

G = deltaq_gp

! qinyi 2017-1-1 21:55:23
! this part is for the consistency b/t liquid water latent heating and
! temperature change
adjust_mid = -1._r8*a_gp*(cld_frac_out+1._r8/sqrt(2._r8*pi)*exp(-0.5_r8*Q1_gp**2)*(-1._r8)*Q1_gp)*0.622_r8*latvap*qstl_gp/(rair*tl_gp**2*thl_gp**2)
adjust_factor = cpair*T_in/(latvap*theta_gp)/(cpair*T_in/(latvap*theta_gp)-adjust_mid)


if(deltaq_gp.ge.0._r8) then
    deltaq_sat = deltaq_gp

    if(sgm_in.ne.0_r8)then
    Q1_sat = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in))
    else
    Q1_sat = 0._r8
    endif

else
    deltaq_uns = deltaq_gp

    if(sgm_in.ne.0_r8)then
    Q1_uns = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in))
    else
    Q1_uns = 0._r8
    endif


endif

end subroutine cldfrac_conden_single_only


subroutine cldfrac_conden_single(p_in,T_in,qv_in,ql_in,sgm_in,cld_frac_out,conden_out,G)

real(r8),intent(in) :: p_in
real(r8),intent(in) :: T_in
real(r8),intent(in) :: qv_in
real(r8),intent(in) :: ql_in
real(r8),intent(in) :: sgm_in

real(r8),intent(out) :: cld_frac_out
real(r8),intent(out) :: conden_out
real(r8),intent(out) :: G

real(r8),parameter:: pi = 3.14159265358_r8
real(r8),parameter:: sgm_min = 1.e-8_r8 ! minimum value of SGV
real(r8),parameter:: p00 = 100000._r8

!local variables
integer :: i,j,k

real(r8) :: conden !ql/(2*sigma(s))
real(r8) :: cld_frac

real(r8) :: exner_gp
real(r8) :: tl_gp
real(r8) :: thl_gp
real(r8) :: rho_gp
real(r8) :: qw_gp

real(r8) :: qs_gp
real(r8) :: es_gp
real(r8) :: qstl_gp
real(r8) :: estl_gp
real(r8) :: Q1_gp
real(r8) :: deltaq_gp
real(r8) :: a_gp
real(r8) :: b_gp
real(r8) :: beta_gp

!!!!!!!!!!!!!!!
! initialized 
!!!!!!!!!!!!!!!

! output variables
cld_frac_out = 0._r8
conden_out = 0._r8
G = 0._r8

! local variables

conden = 0._r8
cld_frac = 0._r8

exner_gp = 0._r8
tl_gp = 0._r8
thl_gp = 0._r8
rho_gp = 0._r8
qw_gp = 0._r8

qs_gp = 0._r8
es_gp = 0._r8
qstl_gp = 0._r8
estl_gp = 0._r8
Q1_gp = 0._r8
deltaq_gp = 0._r8
a_gp = 0._r8
b_gp = 0._r8
beta_gp = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!
!calculation starts here
!!!!!!!!!!!!!!!!!!!!!!!!!

exner_gp = (p_in/p00)**0.286_r8
tl_gp = T_in-latvap/cpair*ql_in
thl_gp = tl_gp/exner_gp
rho_gp = p_in/(rair*T_in)
qw_gp = qv_in+ql_in

call qsat_water(tl_gp,p_in,estl_gp,qstl_gp)
call qsat_water(t_in,p_in,es_gp,qs_gp)

beta_gp = 0.622_r8*latvap*qstl_gp/rair/tl_gp**2

a_gp = 1._r8/(1._r8+beta_gp*latvap/cpair)
b_gp = exner_gp*beta_gp/(1._r8+beta_gp*latvap/cpair)

deltaq_gp = qw_gp-qstl_gp

if(sgm_in.ne.0_r8)then
        Q1_gp = a_gp*deltaq_gp/(2._r8*sqrt(sgm_in))
else
        Q1_gp = 0._r8
endif

! calculate cloud fraction and condensate
if(sqrt(sgm_in)>sgm_min)then
    cld_frac = 0.5_r8*(1._r8+erf(Q1_gp/sqrt(2._r8)))
    conden = cld_frac*Q1_gp+exp(-0.5_r8*Q1_gp**2._r8)/sqrt(2._r8*pi)
    
    cld_frac_out = min(1._r8,max(0._r8,cld_frac))
    conden_out = conden*2._r8*sqrt(sgm_in)
    conden_out = max(0._r8,conden_out)

else
    if(deltaq_gp.le.0._r8)then
      cld_frac = 0._r8
      conden = 0._r8
      cld_frac_out = min(1._r8,max(0._r8,cld_frac))
      conden_out = conden*2._r8*sqrt(sgm_in)
      conden_out = max(0._r8,conden_out)

    else
      cld_frac = 1._r8
      conden_out = a_gp*deltaq_gp
      cld_frac_out = min(1._r8,max(0._r8,cld_frac))
      conden_out = max(0._r8,conden_out)

    endif
endif

if((conden_out.eq.0._r8).or.(cld_frac_out.lt.1.e-4_r8))then
conden_out = 0._r8
cld_frac_out = 0._r8
endif

G = deltaq_gp

end subroutine cldfrac_conden_single


!===================================================

end module cldfrac_conden_mod
