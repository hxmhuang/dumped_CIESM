!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_reftra.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_reftra

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrrft, hnamrft

      implicit none

      contains

! --------------------------------------------------------------------
!      subroutine reftra_sw(nlayers, lrtchk, pgg, prmuz, ptau, pw, &
!                           pref, prefd, ptra, ptrad)
! --------------------------------------------------------------------
! --------------------------ZWJ---------------------------------------
       subroutine reftra_sw(kmodts,nlayers, omega1,omega2,omega3, &
                           lrtchk, pgg, prmuz, ptau, pw, &
                           pref, prefd, ptra, ptrad, &
                           prdf, ptdf, prdr, ptdr)
! --------------------------ZWJ---------------------------------------
! --------------------------------------------------------------------

! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer, intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                           ! Logical flag for reflectivity and
                                                                 ! and transmissivity calculation; 
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ptau(:)                     ! optical depth
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pw(:)                       ! single scattering albedo 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: prmuz                       ! cosine of solar zenith angle

! ------- Output -------

      real(kind=r8), intent(inout) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
! ------------------------------------------------------------------------------------------
! --------------------------------ZWJ-------------------------------------------------------
! --------Input---------
      real(kind=r8), intent(in) :: omega1(:),omega2(:),omega3(:)
      integer, intent(in) :: kmodts
! --------Output--------
!  layer diffuse reflection    dimensions: (nlayers,2,2)
      real(kind=r8), intent(inout) :: prdf(:,:,:)

!  layer diffuse transmission    dimensions: (nlayers,2,2)
      real(kind=r8), intent(inout) :: ptdf(:,:,:)

!  layer direct reflection   dimensions: (nlayers,2)
      real(kind=r8), intent(inout) :: prdr(:,:)

!  layer direct transmission   dimensions: (nlayers,2)
      real(kind=r8), intent(inout) :: ptdr(:,:)

! --------Local---------
      real(kind=r8) :: ome1,ome2,ome3
! --------------------------------ZWJ-------------------------------------------------------
! ------------------------------------------------------------------------------------------

! ------- Local -------
 
      integer :: jk, jl
      integer :: itind

      real(kind=r8) :: tblind, ze0
      real(kind=r8) :: za, za1, za2
      real(kind=r8) :: zbeta, zdend, zdenr, zdent
      real(kind=r8) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=r8) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=r8) :: zr1, zr2, zr3, zr4, zr5
      real(kind=r8) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=r8) :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=r8) :: zw, zwcrit, zwo

      real(kind=r8), parameter :: eps = 1.e-08_r8

! ----------------------------------------------------------------------------------
! --------------------------------ZWJ-----------------------------------------------
! ---------local---------
      real(kind=r8) :: d9,dmu,rmu2,rmu3,dmu2,dmu4
      real(kind=r8) :: a0,a1,a2,a3,a01,a23,a03
      real(kind=r8) :: b0,b1,b2,b3,a3b2,a0b1
      real(kind=r8) :: dmub0,dmub3
      real(kind=r8) :: uu,uu2
      real(kind=r8) :: alams1,alams2,alam1,alam2
      real(kind=r8) :: tr1,tr2,delta
      real(kind=r8) :: eta0,eta1,eta2,eta3
      real(kind=r8) :: p1,q1,r1,p2,q2,r2,h1,h2,h3,h4
      real(kind=r8) :: w11,w12,w13,w14,w21,w22,w23,w24
      real(kind=r8) :: wa,wb,wc,wd,we,wf
      real(kind=r8) :: det,ya,yb,yc,yd,ye,yf,yg,yh
      real(kind=r8) :: c1,d1,c2,d2,c1t,c2t,d1t,d2t
      real(kind=r8) :: teu0,teu1,teu2,teu3
      real(kind=r8) :: ted0,ted1,ted2,ted3,vv,pdtr(nlayers+1)
! --------------------------------ZWJ----------------------------------------------
! ---------------------------------------------------------------------------------

!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision: 1.2 $'

      zsr3=sqrt(3._r8)
      zwcrit=0.9999995_r8
!      kmodts=2
     
      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._r8
            ptra(jk) =1._r8
            prefd(jk)=0._r8
            ptrad(jk)=1._r8
! -----------------------------------------------------------------------------------
! -----------------------------------ZWJ---------------------------------------------
            ptdf(jk,1,1) = 1._r8
            ptdf(jk,2,1:2) = 0._r8
            ptdf(jk,1,2) = 0._r8
            prdf(jk,1,1:2) = 0._r8
            prdf(jk,2,1:2) = 0._r8
            ptdr(jk,1) = 1._r8
            ptdr(jk,2) = 0._r8
            prdr(jk,1) = 0._r8
            prdr(jk,2) = 0._r8   
! -----------------------------------ZWJ---------------------------------------------
! -----------------------------------------------------------------------------------
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)

! -----------------------------------------------------------------------------------
! ------------------------------------ZWJ--------------------------------------------
       if (kmodts/=4) then
! ------------------------------------ZWJ--------------------------------------------  
! -----------------------------------------------------------------------------------
! General two-stream expressions

            zg3= 3._r8 * zg
            if (kmodts == 1) then
               zgamma1= (7._r8 - zw * (4._r8 + zg3)) * 0.25_r8
               zgamma2=-(1._r8 - zw * (4._r8 - zg3)) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 2) then  
               zgamma1= (8._r8 - zw * (5._r8 + zg3)) * 0.25_r8
               zgamma2=  3._r8 *(zw * (1._r8 - zg )) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 3) then  
               zgamma1= zsr3 * (2._r8 - zw * (1._r8 + zg)) * 0.5_r8
               zgamma2= zsr3 * zw * (1._r8 - zg ) * 0.5_r8
               zgamma3= (1._r8 - zsr3 * zg * prmuz ) * 0.5_r8
            end if
            zgamma4= 1._r8 - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            zwo= zw / (1._r8 - (1._r8 - zw) * (zg / (1._r8 - zg))**2)
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500._r8)
!               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  ze2 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  ze2 = exp_tbl(itind)
               endif
!

               pref(jk) = (zgt - za1 * (1._r8 - ze2)) / (1._r8 + zgt)
               ptra(jk) = 1._r8 - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._r8 + zgt)
               ptrad(jk) = 1._r8 - prefd(jk)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0_r8) then 
                  pref(jk) = 0.0_r8
                  ptra(jk) = 1.0_r8
                  prefd(jk) = 0.0_r8
                  ptrad(jk) = 1.0_r8
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._r8 + zrp
               zrm1 = 1._r8 - zrp
               zrk2 = 2._r8 * zrk
               zrpp = 1._r8 - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5
               zbeta = (zgamma1 - zrk) / zrkg !- zr5 / zr4
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._r8)
               ze2 = min ( zto1 / prmuz , 500._r8)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
!              zep1 = exp( ze1 )
!              zem1 = 1._r8 / zep1
!              zep2 = exp( ze2 )
!              zem2 = 1._r8 / zep2
!
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  zem1 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  zep1 = 1._r8 / zem1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zem1 = exp_tbl(itind)
                  zep1 = 1._r8 / zem1
               endif

               if (ze2 .le. od_lo) then 
                  zem2 = 1._r8 - ze2 + 0.5_r8 * ze2 * ze2
                  zep2 = 1._r8 / zem2
               else
                  tblind = ze2 / (bpade + ze2)
                  itind = tblint * tblind + 0.5_r8
                  zem2 = exp_tbl(itind)
                  zep2 = 1._r8 / zem2
               endif

! collimated beam

               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(jk) = eps
                  ptra(jk) = zem2
               else
                  pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif 

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._r8 / ( (1._r8 - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._r8 - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

! --------------------------------------------------------------------------------------
! ------------------------------------ZWJ-----------------------------------------------
         else 
! Initialize
             d9 = 1.0_r8 / 9.0_r8
             dmu = 1.0_r8 / prmuz 
             rmu2 = prmuz * prmuz
             rmu3 = prmuz * rmu2
             dmu2 = dmu * dmu
             dmu4 = dmu2 * dmu2
! four-stream
             zw = min(zw,0.999999995_r8)

             ome1 = omega1(jk) * zw 
             ome2 = omega2(jk) * zw 
             ome3 = omega3(jk) * zw
  
             a0 = 1.0_r8 - zw 
             a1 = 3.0_r8 - ome1
             a2 = 5.0_r8 - ome2 
             a3 = 7.0_r8 - ome3 
   
             a01 = a0 * a1
             a23 = a2 * a3
             a03 = a0 * a3

             b0 = 0.25_r8 * zw
             b1 = -0.25_r8 * ome1 * prmuz
             b2 = 0.125_r8 * ome2 * (3.0_r8 * rmu2 -1.0_r8)
             b3 = -0.125_r8 * ome3 * (5.0_r8 * rmu3 - 3.0_r8 * prmuz)
             a3b2 = a3 * b2
             a0b1 = a0 * b1
             dmub0 = dmu * b0
             dmub3 = dmu * b3

             uu = a01 + d9 * (a23 + 4.0_r8 * a03)
             uu2 = uu * uu
             vv = 4.0_r8 * d9 * a01 * a23

             alams1 = 0.5_r8 * (uu + sqrt(uu2 - vv))
             alams2 = 0.5_r8 * (uu - sqrt(uu2 - vv))
             alam1 = sqrt(alams1)
             alam2 = sqrt(alams2)
 
             ze1 = zto1 * dmu
             ze2 = alam1 * zto1
               if (ze1 .le. 1.e-8_r8) then
                  zem1 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  zep1 = 1._r8 / zem1
                  pdtr(jk) = zep1
               else
                  pdtr(jk) = exp((-1.0_r8)* zto1 * dmu)
               endif
               if (ze2 .le. 1.e-8_r8) then
                  zem2 = 1._r8 - ze2 + 0.5_r8 * ze2 * ze2
                  zep2 = 1._r8 / zem2
                  tr1 = zep2
               else
                  tr1 = exp((-1.0_r8)* alam1 * zto1)
               endif
               
             ze2 = alam2 * zto1
               if (ze2 .le. 1.e-8_r8) then
                  zem2 = 1._r8 - ze2 + 0.5_r8 * ze2 * ze2
                  zep2 = 1._r8 / zem2
                  tr2 = zep2
               else
                  tr2 = exp((-1.0_r8)* alam2 * zto1)
               endif

             delta = 1.0_r8/(9.0_r8*dmu4-dmu2*(9.0*a01+a23+4.0*a03)+a01*a23)
             eta0 = ((a1*b0-dmu*b1)*(a23-9.0*dmu2)+2.0*dmu2* &
                    (a3b2-2.0*a3*b0-3.0*dmub3))*delta
             eta1 = ((a0b1-dmub0)*(a23-9.0*dmu2)-2.0*dmu*  &
                    a0*(a3b2-3.0*dmub3))*delta
             eta2 = 0.625*((a3b2-3.0*dmub3)*(a01-dmu2)-   &
                    2.0_r8*dmu*a3*(a0b1-dmub0))*delta
             eta3 = ((a2*b3-3.0*dmu*b2)*(a01-dmu2)+   &
                    dmu2*(6.0*a0b1-4.0*a0*b3-6.0*dmub0))*delta

             p1 = -a0 / alam1
             q1 =  0.3125 * (a01 / alams1 - 1.0)

             r1 = -1.5 * (a01 / alam1 - alam1) / a3
             p2 = -a0 / alam2
             q2 =  0.3125 * (a01 / alams2 - 1.0)

             r2 = -1.5 * (a01 / alam2 - alam2) / a3

             h1 = -(0.5 * eta0 - eta1 + eta2)
             h2 = -(-0.125 * eta0 + eta2 - eta3)
             h3 = -(0.5 * eta0 + eta1 + eta2) * pdtr(jk)
             h4 = -(-0.125 * eta0 + eta2 + eta3) * pdtr(jk)

             w11 = 0.5 - p1 + q1
             w12 = (0.5 + p1 + q1) * tr1
             w13 = 0.5 - p2 + q2
             w14 = (0.5 + p2 + q2) * tr2
             w21 = -0.125 + q1 - r1
             w22 = (-0.125 + q1 + r1) * tr1
             w23 = -0.125 + q2 - r2
             w24 = (-0.125 + q2 + r2) * tr2

             wa =  w11 * w22 - w21 * w12
             wb =  w14 * w23 - w24 * w13
             wc =  w11 * w23 - w21 * w13
             wd =  w11 * w24 - w21 * w14
             we =  w12 * w23 - w22 * w13
             wf =  w12 * w24 - w22 * w14

             det = 1.0/(2.0*wa*wb-wc*wc+wd*wd+we*we-wf*wf)
             ya = ( w22 * wb - w23 * wc + w24 * wd) * det
             yb = (-w12 * wb + w13 * wc - w14 * wd) * det
             yc = ( w23 * we - w24 * wf - w21 * wb) * det
             yd = (-w13 * we + w14 * wf + w11 * wb) * det
             ye = ( w21 * wc - w22 * we - w24 * wa) * det
             yf = (-w11 * wc + w12 * we + w14 * wa) * det
             yg = ( w23 * wa - w21 * wd + w22 * wf) * det
             yh = (-w13 * wa + w11 * wd - w12 * wf) * det

             c1 =  ya * h1 + yb * h2 + yc * h3 + yd * h4
             d1 =  yc * h1 + yd * h2 + ya * h3 + yb * h4
             c2 =  ye * h1 + yf * h2 + yg * h3 + yh * h4
             d2 =  yg * h1 + yh * h2 + ye * h3 + yf * h4

             c1t =  c1 * tr1
             c2t =  c2 * tr2
             d1t =  d1 * tr1
             d2t =  d2 * tr2

             teu0 =  c1 + d1t + c2 + d2t + eta0
             teu1 =  p1 * (c1 - d1t) + p2 * (c2 - d2t) + eta1
             teu2 =  q1 * (c1 + d1t) + q2 * (c2 + d2t) + eta2
             teu3 =  r1 * (c1 - d1t) + r2 * (c2 - d2t) + eta3
  
! --------direct incident--------

             prdr(jk,1) = (teu0 + 2.0 * (teu1 + teu2)) * dmu
             prdr(jk,2) =  2.0 * (-0.125 * teu0 + teu2 + teu3) * dmu
 
             ted0 =  c1t + d1 + c2t + d2 + eta0 * pdtr(jk)
             ted1 =  p1 * (c1t - d1) + p2 * (c2t - d2) + eta1 * pdtr(jk)
             ted2 =  q1 * (c1t + d1) + q2 * (c2t + d2) + eta2 * pdtr(jk)
             ted3 =  r1 * (c1t - d1) + r2 * (c2t - d2) + eta3 * pdtr(jk)

             ptdr(jk,1) = (ted0 - 2.0 * (ted1 - ted2))  * dmu
             ptdr(jk,2) =  2.0 * (-0.125 * ted0 + ted2 - ted3) * dmu

! --------diffuse radiation-------

             c1t =  ya * tr1
             d1t =  yc * tr1
             c2t =  ye * tr2
             d2t =  yg * tr2

             teu0 =  ya + d1t + ye + d2t
             teu1 =  p1 * (ya - d1t) + p2 * (ye - d2t)
             teu2 =  q1 * (ya + d1t) + q2 * (ye + d2t)
             teu3 =  r1 * (ya - d1t) + r2 * (ye - d2t)
             prdf(jk,1,1) =  0.5_r8 * teu0 + teu1 + teu2
             prdf(jk,2,1) = -0.125_r8 * teu0 + teu2 + teu3
 
             ted0 =  c1t + yc + c2t + yg
             ted1 =  p1 * (c1t - yc) + p2 * (c2t - yg)
             ted2 =  q1 * (c1t + yc) + q2 * (c2t + yg)
             ted3 =  r1 * (c1t - yc) + r2 * (c2t - yg)
             ptdf(jk,1,1) =  0.5_r8 * ted0 - ted1 + ted2
             ptdf(jk,2,1) = -0.125_r8 * ted0 + ted2 - ted3

             c1t =  yb * tr1
             d1t =  yd * tr1
             c2t =  yf * tr2
             d2t =  yh * tr2

             teu0 =  yb + d1t + yf + d2t
             teu1 =  p1 * (yb - d1t) + p2 * (yf - d2t)
             teu2 =  q1 * (yb + d1t) + q2 * (yf + d2t)
             teu3 =  r1 * (yb - d1t) + r2 * (yf - d2t)
             prdf(jk,1,2) =  0.5_r8 * teu0 + teu1 + teu2
             prdf(jk,2,2) = -0.125_r8 * teu0 + teu2 + teu3

             ted0 =  c1t + yd + c2t + yh
             ted1 =  p1 * (c1t - yd) + p2 * (c2t - yh)
             ted2 =  q1 * (c1t + yd) + q2 * (c2t + yh)
             ted3 =  r1 * (c1t - yd) + r2 * (c2t - yh)
             ptdf(jk,1,2) =  0.5_r8 * ted0 - ted1 + ted2
             ptdf(jk,2,2) = -0.125_r8 * ted0 + ted2 - ted3
 
         endif
! ------------------------------------ZWJ-----------------------------------------------
! --------------------------------------------------------------------------------------
         endif         

      end do  

      end subroutine reftra_sw

      end module rrtmg_sw_reftra

