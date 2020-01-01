!-----------------------------------------------
!--Subroutines for oro_drag scheme--------------
!--Added by Liang Yishuang on 2016.08.29--------
!-----------------------------------------------

 subroutine zm_calculate( nx, ib, ie, jb, je, kb, ke, pe, pt, zm)

   use shr_kind_mod,  only : r8 => shr_kind_r8


   implicit none

   real(r8), parameter :: cappa  = 0.2857
   real(r8), parameter :: rog    = 29.2692

   integer , intent(in   ) :: nx
   integer , intent(in   ) :: ib, ie, jb, je, kb, ke
   real(r8), intent(in   ) :: pe(ib:ie, kb:ke+1, jb:je) 
   real(r8), intent(in   ) :: pt(ib:ie, jb:je,   kb:ke) 
   real(r8), intent(inout) :: zm(ib:ie, jb:je,   kb:ke) 

   !--------------------local variables--------------------
   integer  :: i, j, k, m, n
   real(r8) :: pk  (ib:ie, jb:je,  kb:ke+1) 
   real(r8) :: pkz (ib:ie, jb:je,  kb:ke  )
   real(r8) :: peln(ib:ie, kb:ke+1,jb:je)
   real(r8) :: zi  (ib:ie, jb:je,  kb:ke+1)
   real(r8) :: pk2 (ib:ie, kb:ke+1)
   real(r8) :: pek, lnp
   real(r8) :: rpdel, hkl, hkk, tv 
   integer  :: ixj, jp, it, i1, i2, itot, nxu

   !------------------------------------------------------
   nxu  = nx
   itot = ie -ib + 1
   it = itot / nxu
   jp = nxu * ( je -jb + 1 )

   do k = kb, ke+1
      do j = jb, je
         do i = ib, ie
            pk(i,j,k) = pe(i,k,j) ** cappa
         enddo
      enddo
   enddo
   !write(6,*) 'itot & jp=', itot, jp , grid%ilastxy, grid%ifirstxy, grid%jlastxy, grid%jfirstxy
  !---------------------Calculate virtual Temperature-----------
   do ixj = 1, jp
      j  = jb + ( ixj -1 ) /nxu
      i1 = ib + it * mod( ixj - 1, nxu )
      i2 = i1 + it - 1
     ! write(6,*) 'j i1 & i2=', j, i1, i2

      if( kb == 1 ) then
         pek = pk(i1,j,1)
         lnp = log ( pe(i1,1,j) )
      endif
      do i = i1, i2
         pk2 (i,1) = pek
         peln(i,1,j) = lnp
      enddo
     ! write(6,*) '---------------pk2, pek  OK-----------'
      do k = max( 2, kb ), ke + 1
         do i = i1, i2
            peln(i,k,j) = log( pe(i,k,j) )
            pk2 (i,k)   = pk(i,j,k)
         enddo
      enddo
     ! write(6,*) '---------------peln pk2  OK-----------'
      do k = kb, ke
         do i = i1, i2
            pkz(i,j,k) = ( pk2(i,k+1) - pk2(i,k) ) / ( cappa * ( peln(i,k+1,j) - peln(i,k,j) ) )
         enddo
      enddo

   enddo  ! ixj  

      !-------------------------------Caculated ZM-------------------------------
   !zi = 0.0
   zi = 0.0
   do k = ke, kb, -1
      do j = jb, je
         do i = ib, ie
            rpdel = 1. / ( pe(i,k+1,j) - pe(i,k,j) )
            !hkl = log( pe(i,k+1,j) ) - log( pe(i,k,j) )
            hkl = peln(i,k+1,j) - peln(i,k,j)
            hkk = 1 - pe(i,k,j) * hkl * rpdel
            tv  = pt(i,j,k) * pkz(i,j,k)
            zm(i,j,k) = zi(i,j,k+1) + rog * tv * hkk 
            zi(i,j,k) = zi(i,j,k+1) + rog * tv * hkl
!            write(6,*) ' I, J, K=', i,j,k
!            write(6,*) ' RPDEL=', rpdel
!            write(6,*) 'hkk & hkl=', hkk, hkl
!            write(6,*) ' Pt & PKZ=', pt(i,j,k), pkz(i,k)
!            write(6,*) ' Tv & ZM=', tv, zm(i,j,k)
         enddo
      enddo 
   enddo ! Count ixj 
   !write(6,*) '------------------Finished the PKZ Calculated'

!   write(6,*) ' Vertical Coor=', kb, ke
!   do k = kb, ke
!      write(6,*) 'Max & Min ZM=', k, maxval(zm(:,:,k)), minval(zm(:,:,k))
!   enddo
!   write(6,*) '------------------Finished the ZM_WIND  Calculated'
  
 
 end subroutine zm_calculate
 
!================================================================================================= 
 subroutine orography_drag(  ib, ie, jb, je, kb, ke, zm, sgh, U, V, dt, pe, t)

   use shr_kind_mod,  only : r8 => shr_kind_r8
   use shr_const_mod, only : shr_const_rdair
   use cam_history,  only: dyn_decomp, addfld, add_default,outfld
   use spmd_utils, only: npes, masterproc
   use cam_logfile,        only: iulog


   implicit none

   real(r8), parameter :: h1500 = 1. / 1500.
   real(r8), parameter :: a     = 12.
   real(r8), parameter :: b     = 1.
   real(r8), parameter :: Cmd   = 0.005
   real(r8), parameter :: Ccorr = 0.6
   real(r8), parameter :: Ih    = 0.00102
   real(r8), parameter :: Kflt  = 0.00035
   real(r8), parameter :: n1    = -1.9
   real(r8), parameter :: n2    = -2.8
   real(r8), parameter :: z1    = 1.5
   real(r8), parameter :: z2    = -1.2
   real(r8), parameter :: k1    = 0.003
   real(r8), parameter :: c     = 2.109
   real(r8), parameter :: aa    = a * b * Cmd * Ccorr * c
!---------------liangys 2016-08-29 for exponent type integer-------------
!   real(r8), parameter :: a1    = 1.0 / ( Ih * (Kflt ** n1) )
!   real(r8), parameter :: a2    = a1 * ( k1 ** (n1-n2) )
!   real(r8), parameter :: ax    = aa * a2
   real(r8) :: a1
   real(r8) :: a2
   real(r8) :: ax
!---------------liangys 2016-08-29 for exponent type integer-------------
   real(r8), parameter :: rair  = shr_const_rdair


   integer,  intent(in   ) :: ib, ie
   integer,  intent(in   ) :: jb, je
   integer,  intent(in   ) :: kb, ke
   real(r8), intent(in   ) :: zm (ib:ie, jb:je, kb:ke)
   real(r8), intent(in   ) :: sgh(ib:ie, jb:je)
   real(r8), intent(in   ) :: pe (ib:ie, jb:je, kb:ke)     !pressure
   real(r8), intent(in   ) :: t  (ib:ie, jb:je, kb:ke)     !tempurature
   real(r8), intent(inout) :: u (ib:ie, jb:je, kb:ke)
   real(r8), intent(inout) :: v (ib:ie, jb:je, kb:ke)
   real(r8), intent(in   ) :: dt

   !--------------------local variables--------------------
   integer  :: i, j, k, m, n
   real(r8) :: uv, z12, expz, z1500
   real(r8) :: ux(ib:ie, jb:je, kb:ke)
   real(r8) :: vx(ib:ie, jb:je, kb:ke)
!  real(r8) :: rho(ib:ie, jb:je, kb:ke)    !density of air
   real(r8) :: rho    !density of air
   real(r8) :: ftaux(ib:ie, jb:je, kb:ke)   !mid var
   real(r8) :: ftauy(ib:ie, jb:je, kb:ke)   !mid var
   real(r8) :: taubbwx(ib:ie, jb:je)       !taux
   real(r8) :: taubbwy(ib:ie, jb:je)       !tauy
   real(r8) :: staux(ib:ie, jb:je)        !mid var
   real(r8) :: stauy(ib:ie, jb:je)        !mid var
   integer  :: idim
   logical, save :: first_call = .true.

!---------------liangys 2017-11-29 for 0.7lev-------------
   real(r8) :: rdc !reduce coff from z
   real(r8) :: u_save !
   real(r8) :: v_save !
!---------------liangys 2017-11-29 for 0.7lev-------------
   ftaux = u
   ftauy = v

!---------------liangys 2016-08-29 for exponent type integer-------------
   a1    = 1.0 / ( Ih * (Kflt ** n1) )
   a2    = a1 * ( k1 ** (n1-n2) )
   ax    = aa * a2
!---------------liangys 2016-08-29 for exponent type integer-------------

!   write(iulog,*) '------------------LIANGYS KB KE-------------------',kb,ke
   do j = jb, je
       do i = ib, ie 
         do k = kb, ke
            uv    = sqrt( u(i,j,k) **2 + v(i,j,k) ** 2 )
            z1500 = ( zm(i,j,k) * h1500 ) ** z1
            expz  = exp ( - z1500 )
!---------------liangys 2017-11-29 for 0.7lev-------------
            rdc   = 1. + 0.55  * expz
!---------------liangys 2017-11-29 for 0.7lev-------------
            z12   = zm(i,j,k) ** z2
            rho    = - pe(i,j,k) / ( rair * t(i,j,k) )
            ftaux(i,j,k) = rho * ax * uv * expz * z12 * sgh(i,j) * u(i,j,k) 
            ftauy(i,j,k) = rho * ax * uv * expz * z12 * sgh(i,j) * v(i,j,k)
!---------------liangys 2017-11-29 for 0.7lev-------------
!            u(i,j,k) = u(i,j,k) - dt * ax * uv * expz * z12 * sgh(i,j) * u(i,j,k)            
!            v(i,j,k) = v(i,j,k) - dt * ax * uv * expz * z12 * sgh(i,j) * v(i,j,k)            
            u_save   = u(i,j,k)
            v_save   = v(i,j,k)
            u(i,j,k) =  u(i,j,k) - rdc * ( dt * ax * uv * expz * z12 * sgh(i,j) * u(i,j,k) )            
            v(i,j,k) =  v(i,j,k) - rdc * ( dt * ax * uv * expz * z12 * sgh(i,j) * v(i,j,k) )           
            
!            if ((u_save * u(i,j,k)) .lt. 0.) then
!               write(iulog,*) 'u wind zhuanxiang, u_before = ',u_save,', u_after = ',u(i,j,k)
!               u(i,j,k) = 0.
!            end if
!            if ((v_save * v(i,j,k)) .lt. 0.) then
!               write(iulog,*) 'v wind zhuanxiang, v_before = ',v_save,', v_after = ',v(i,j,k)
!               v(i,j,k) = 0.
!            end if
!---------------liangys 2017-11-29 for 0.7lev-------------
        enddo 
      enddo   
   enddo 

   taubbwx(:,:) = 0._r8
   taubbwy(:,:) = 0._r8
   staux(:,:)   = 0._r8
   stauy(:,:)   = 0._r8
   idim = ie - ib + 1
   do j = jb, je
       do i = ib, ie 
         do k = kb, ke
         
         if (k .lt. ke) then
            staux(i,j) = ( ftaux(i,j,k) + ftaux(i,j,k+1) ) * ( zm(i,j,k) - zm(i,j,k+1) ) * 0.5_r8
            stauy(i,j) = ( ftauy(i,j,k) + ftauy(i,j,k+1) ) * ( zm(i,j,k) - zm(i,j,k+1) ) * 0.5_r8
         else
            staux(i,j) =  ftaux(i,j,k)  *  zm(i,j,k)
            stauy(i,j) =  ftauy(i,j,k)  *  zm(i,j,k)
         end if

            taubbwx(i,j) = taubbwx(i,j) + staux(i,j)
            taubbwy(i,j) = taubbwy(i,j) + stauy(i,j)

        enddo 
      enddo   
   if( .not. first_call ) then
!     write(iulog,*) '------------------LIANGYS CALL OUTFLD-------------------',first_call
     call outfld( 'TAUBBWX'  , taubbwx, idim, j )
     call outfld( 'TAUBBWY'  , taubbwy, idim, j )   
   end if 
     
   enddo 

   first_call = .false.
   
 
 end subroutine orography_drag
 
 
