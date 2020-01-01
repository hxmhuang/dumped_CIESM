!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:15 $
!
      module rrtmg_sw_vrtqdr

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

!      use parkind, only: jpim, jprb
!      use parrrsw, only: ngptsw

      implicit none

      contains

! --------------------------------------------------------------------------
!      subroutine vrtqdr_sw(klev, kw, &
!                           pref, prefd, ptra, ptrad, &
!                           pdbt, prdnd, prup, prupd, ptdbt, &
!                           pfd, pfu)
! --------------------------------------------------------------------------
! ----------------------------ZWJ-------------------------------------------
       subroutine vrtqdr_sw(kmodts, klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu, prdf, ptdf, prdr, ptdr)
! ----------------------------ZWJ-------------------------------------------
! --------------------------------------------------------------------------

! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer, intent (in) :: klev                   ! number of model layers
      integer, intent (in) :: kw                     ! g-point index

      real(kind=r8), intent(in) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(in) :: pdbt(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptdbt(:)
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(inout) :: prdnd(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prup(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prupd(:)
                                                                 !   Dimensions: (nlayers+1)

! Output
      real(kind=r8), intent(out) :: pfd(:,:)                   ! downwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle
      real(kind=r8), intent(out) :: pfu(:,:)                   ! upwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle

! Local

      integer :: ikp, ikx, jk

      real(kind=r8) :: zreflect
      real(kind=r8) :: ztdn(klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! -------------------------------ZWJ------------------------------------------
! --------Input--------

!  layer diffuse reflection    dimensions: (nlayers,2,2)
      real(kind=r8), intent(in) :: prdf(:,:,:)

!  layer diffuse transmission    dimensions: (nlayers,2,2)
      real(kind=r8), intent(in) :: ptdf(:,:,:)

!  layer direct reflection   dimensions: (nlayers,2)
      real(kind=r8), intent(in) :: prdr(:,:)

!  layer direct transmission   dimensions: (nlayers,2)
      real(kind=r8), intent(in) :: ptdr(:,:)

      integer, intent(in) :: kmodts

! ----------------local----------------

      real(kind=r8) :: rmur(klev+1,2),rmuf(klev+1,2,2)
      real(kind=r8) :: c11, c12, c21, c22
      real(kind=r8) :: rmdf(klev+1,2,2)
      real(kind=r8) :: pmod, t11, t12, t21, t22, b11, b12, b21, b22
      real(kind=r8) :: d1, d2, uu1, uu2, du1, du2, a11, a12, a21, a22
      real(kind=r8) :: r11, r12, r21, r22
      real(kind=r8) :: tmdr(klev+1,2)

      integer :: k, l, km1, lp1
! -------------------------------ZWJ------------------------------------------
! ----------------------------------------------------------------------------           

! ----------------------------------------------------------------------------
! -------------------------------ZWJ------------------------------------------
      if (kmodts /= 4) then
! -------------------------------ZWJ------------------------------------------
! ----------------------------------------------------------------------------        
! Link lowest layer with surface
             
      zreflect = 1._r8 / (1._r8 - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._r8 / (1._r8 -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._r8
      prdnd(1) = 0._r8
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._r8 / (1._r8 - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._r8 / (1._r8 - prdnd(jk) * prupd(jk))
         pfu(jk,kw) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(jk,kw) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

! -------------------------------------------------------------------------
! -------------------------ZWJ---------------------------------------------
      else
        tmdr(1,1:2) = 0.0_r8
        rmdf(1,1:2,1:2) = 0._r8
        rmur(klev+1,1)        = prup(klev+1)
        rmur(klev+1,2)        = -0.25 * prup(klev+1)
        rmuf(klev+1,1,1)      = prupd(klev+1)
        rmuf(klev+1,2,1)      = -0.25 * prupd(klev+1)
        rmuf(klev+1,1:2,2)    = 0.0

! ... add the layers upward from one layer above surface to the klev. 

        do  k = 2, klev+1
            km1 = k - 1
            l = klev+1 - k + 1
            lp1 = l + 1

            c11            =  1.0 - prdf(km1,1,1) * rmdf(km1,1,1) - &
                                prdf(km1,1,2) * rmdf(km1,2,1)
            c12            =      - prdf(km1,1,1) * rmdf(km1,1,2) - &
                                prdf(km1,1,2) * rmdf(km1,2,2)
            c21            =      - prdf(km1,2,1) * rmdf(km1,1,1) - &
                                prdf(km1,2,2) * rmdf(km1,2,1)
            c22            =  1.0 - prdf(km1,2,1) * rmdf(km1,1,2) - &
                                prdf(km1,2,2) * rmdf(km1,2,2)
            pmod           =  c11 * c22 - c12 * c21 
            t11            =  c22 / pmod
            t12            = -c12 / pmod
            t21            = -c21 / pmod
            t22            =  c11 / pmod

            b11            =  t11 * ptdf(km1,1,1) + &
                              t12 * ptdf(km1,2,1)
            b12            =  t11 * ptdf(km1,1,2) + &
                              t12 * ptdf(km1,2,2)
            b21            =  t21 * ptdf(km1,1,1) + &
                              t22 * ptdf(km1,2,1)
            b22            =  t21 * ptdf(km1,1,2) + &
                              t22 * ptdf(km1,2,2)

            d1             =  prdr(km1,1) * ptdbt(km1) + &
                              prdf(km1,1,1) * tmdr(km1,1) + &
                              prdf(km1,1,2) * tmdr(km1,2)
            d2             =  prdr(km1,2) * ptdbt(km1) + &
                              prdf(km1,2,1) * tmdr(km1,1) + &
                              prdf(km1,2,2) * tmdr(km1,2)
            uu1            =  d1 * t11 + d2 * t12
            uu2            =  d1 * t21 + d2 * t22
            du1            =  tmdr(km1,1) + uu1 * rmdf(km1,1,1) + &
                                 uu2 * rmdf(km1,1,2)
            du2            =  tmdr(km1,2) + uu1 * rmdf(km1,2,1) + &
                                 uu2 * rmdf(km1,2,2)
            tmdr(k,1)   =  ptdr(km1,1) * ptdbt(km1) + du1 * &
                              ptdf(km1,1,1) + du2 * ptdf(km1,1,2)
            tmdr(k,2)   =  ptdr(km1,2) * ptdbt(km1) + du1 * &
                              ptdf(km1,2,1) + du2 * ptdf(km1,2,2)

            a11            =  ptdf(km1,1,1) * rmdf(km1,1,1) + &
                              ptdf(km1,1,2) * rmdf(km1,2,1)
            a12            =  ptdf(km1,1,1) * rmdf(km1,1,2) + &
                              ptdf(km1,1,2) * rmdf(km1,2,2)
            a21            =  ptdf(km1,2,1) * rmdf(km1,1,1) + &
                              ptdf(km1,2,2) * rmdf(km1,2,1)
            a22            =  ptdf(km1,2,1) * rmdf(km1,1,2) + &
                              ptdf(km1,2,2) * rmdf(km1,2,2)
            rmdf(k,1,1) =  prdf(km1,1,1) + a11 * b11 + a12 * b21
            rmdf(k,1,2) =  prdf(km1,1,2) + a11 * b12 + a12 * b22
            rmdf(k,2,1) =  prdf(km1,2,1) + a21 * b11 + a22 * b21
            rmdf(k,2,2) =  prdf(km1,2,2) + a21 * b12 + a22 * b22

! ... add the layerS upward from layer above surface to the klev.

            c11            =  1.0 - rmuf(lp1,1,1) * prdf(l,1,1) - &
                               rmuf(lp1,1,2) * prdf(l,2,1)
            c12            =      - rmuf(lp1,1,1) * prdf(l,1,2) - &
                               rmuf(lp1,1,2) * prdf(l,2,2)
            c21            =      - rmuf(lp1,2,1) * prdf(l,1,1) - &
                               rmuf(lp1,2,2) * prdf(l,2,1)
            c22            =  1.0 - rmuf(lp1,2,1) * prdf(l,1,2) - &
                               rmuf(lp1,2,2) * prdf(l,2,2)
            pmod           =  c11 * c22 - c12 * c21
            t11            =  c22 / pmod
            t12            = -c12 / pmod
            t21            = -c21 / pmod
            t22            =  c11 / pmod
            d1             =  rmur(lp1,1) * pdbt(l) + &
                              rmuf(lp1,1,1) * ptdr(l,1) + &
                              rmuf(lp1,1,2) * ptdr(l,2)
            d2             =  rmur(lp1,2) * pdbt(l) + &
                              rmuf(lp1,2,1) * ptdr(l,1) + &
                              rmuf(lp1,2,2) * ptdr(l,2)
            uu1            =  d1 * t11 + d2 * t12
            uu2            =  d1 * t21 + d2 * t22
            rmur(l,1)   =  prdr(l,1) + uu1 * ptdf(l,1,1) + &
                                   uu2 * ptdf(l,1,2)
            rmur(l,2)   =  prdr(l,2) + uu1 * ptdf(l,2,1) + &
                                   uu2 * ptdf(l,2,2)

            c11            =  1.0 - prdf(l,1,1) * rmuf(lp1,1,1) - &
                               prdf(l,1,2) * rmuf(lp1,2,1)
            c12            =      - prdf(l,1,1) * rmuf(lp1,1,2) - &
                               prdf(l,1,2) * rmuf(lp1,2,2)
            c21            =      - prdf(l,2,1) * rmuf(lp1,1,1) - &
                               prdf(l,2,2) * rmuf(lp1,2,1)
            c22            =  1.0 - prdf(l,2,1) * rmuf(lp1,1,2) - &
                               prdf(l,2,2) * rmuf(lp1,2,2)
            pmod           =  c11 * c22 - c12 * c21
            t11            =  c22 / pmod
            t12            = -c12 / pmod
            t21            = -c21 / pmod
            t22            =  c11 / pmod

            r11            =  ptdf(l,1,1) * (t11 * rmuf(lp1,1,1) + &
                                      t21 * rmuf(lp1,1,2)) + &
                              ptdf(l,1,2) * (t11 * rmuf(lp1,2,1) + &
                                      t21 * rmuf(lp1,2,2))
            r12            =  ptdf(l,1,1) * (t12 * rmuf(lp1,1,1) + &
                                      t22 * rmuf(lp1,1,2)) + &
                              ptdf(l,1,2) * (t12 * rmuf(lp1,2,1) + &
                                      t22 * rmuf(lp1,2,2))
            r21            =  ptdf(l,2,1) * (t11 * rmuf(lp1,1,1) + &
                                      t21 * rmuf(lp1,1,2)) + &
                              ptdf(l,2,2) * (t11 * rmuf(lp1,2,1) + &
                                      t21 * rmuf(lp1,2,2))
            r22            =  ptdf(l,2,1) * (t12 * rmuf(lp1,1,1) + &
                                      t22 * rmuf(lp1,1,2)) + &
                              ptdf(l,2,2) * (t12 * rmuf(lp1,2,1) + &
                                      t22 * rmuf(lp1,2,2))

            rmuf(l,1,1) =  prdf(l,1,1) + r11 * ptdf(l,1,1) + &
                                           r12 * ptdf(l,2,1) 
            rmuf(l,1,2) =  prdf(l,1,2) + r11 * ptdf(l,1,2) + &
                                           r12 * ptdf(l,2,2)
            rmuf(l,2,1) =  prdf(l,2,1) + r21 * ptdf(l,1,1) + &
                                           r22 * ptdf(l,2,1)
            rmuf(l,2,2) =  prdf(l,2,2) + r21 * ptdf(l,1,2) + &
                                           r22 * ptdf(l,2,2)

        end do

        do k = 1, klev+1

            c11            =  1.0 - rmuf(k,1,1) * rmdf(k,1,1) - &
                               rmuf(k,1,2) * rmdf(k,2,1)
            c12            =      - rmuf(k,1,1) * rmdf(k,1,2) - &
                               rmuf(k,1,2) * rmdf(k,2,2)
            c21            =      - rmuf(k,2,1) * rmdf(k,1,1) - &
                               rmuf(k,2,2) * rmdf(k,2,1)
            c22            =  1.0 - rmuf(k,2,1) * rmdf(k,1,2) - &
                               rmuf(k,2,2) * rmdf(k,2,2)
            pmod           =  c11 * c22 - c12 * c21 
            t11            =  c22 / pmod
            t12            = -c12 / pmod
            t21            = -c21 / pmod
            t22            =  c11 / pmod

            d1             =  rmur(k,1) * ptdbt(k) + &
                              rmuf(k,1,1) * tmdr(k,1) + &
                              rmuf(k,1,2) * tmdr(k,2)
            d2             =  rmur(k,2) * ptdbt(k) + &
                              rmuf(k,2,1) * tmdr(k,1) + &
                              rmuf(k,2,2) * tmdr(k,2)
            uu1            =  d1 * t11 + d2 * t12
            du1            =  tmdr(k,1) + rmdf(k,1,1) * uu1 + &
                              rmdf(k,1,2) * (d1 * t21 + d2 * t22)
            pfu(k,kw)    =  uu1
            pfd(k,kw)    =  du1 + ptdbt(k)
 
        end do

      endif

! -------------------------ZWJ---------------------------------------------
! -------------------------------------------------------------------------

      end subroutine vrtqdr_sw

      end module rrtmg_sw_vrtqdr
