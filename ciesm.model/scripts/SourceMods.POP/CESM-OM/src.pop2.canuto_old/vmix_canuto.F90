

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     beginning of turbulence models
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!030424Z****INTRODUCE THE CORIOLIS PARAMETER "2 \Omega sin \phi" FOR ROTATION MODELS****
!030401Y****INTRODUCE THE SQUARE OF THE BRUNT VAISALA FREQUENCY FOR ZERO SHEAR MODEL****
!020912X******INTRODUCE SURFACE FLUX INPUTS******


module arguments
        
		integer,parameter::isurfuse=1,ifextermld=0,ifoutput=0
		integer,parameter::ifsmooth_ri=0,if_smooth_amld=0        !choice for smooth ri and amld: 0 NO 1 YES
		integer,parameter::nmodel=1,ntbl= 501
		integer,parameter::ifback=5,ifsali=1
		integer,parameter::icondear=0,ifepson2=2
		integer,parameter::ifbotenhance=0,ifdeeplat=1
		integer,parameter::ilomega=0,ifcheckbottomeps=0
		integer,parameter::ifrafgmax=1,ifsalback=5
	  
		integer,parameter::ifexpabstable=1,ifast=1
		integer,parameter::ifastexpabs=ifast*ifexpabstable

	 
		integer,parameter::ifchengcon=0,idefmld=0
		integer,parameter::ifpolartablewrite=0
		integer,parameter::ifbg_theta_interp=1
	  
		logical,parameter::ifshearmin=.TRUE.,ifzeroshear=.TRUE.
        logical,parameter :: if_convect = .true.
end module arguments
	
module constant
		use arguments
        use kinds_mod
		implicit none
		

		
         
         
		integer,parameter::nbig=100
		integer,parameter::nposapprox=101
		integer,parameter::nextrtbl1=1000,nextrtbl0= 62
		integer,parameter::mt_ra_r=nposapprox-1
		integer,parameter::nextrtbl=nextrtbl0+ifexpabstable*nextrtbl1	
		integer,parameter::mt0=ntbl-nposapprox
		integer,parameter::mt=mt0+nextrtbl
  
		real(kind=r8),parameter::e=2.71828182845904509D0
		real(kind=r8),parameter::pi=3.14159265358979312D0
		real(kind=r8),parameter::s2min=1.D-14
		real(kind=r8),parameter::eps_bot0=2.D-5
		real(kind=r8),parameter::epson2__=.288D0
		real(kind=r8),parameter::scale_bot=5.D4
		real(kind=r8),parameter::eplatidependmin=7.D-2
		real(kind=r8),parameter::ri0=- 4.D0
		real(kind=r8),parameter::deltemld=0.5D0
		real(kind=r8),parameter::delrhmld=0.125D-3
		real(kind=r8),parameter::amldminlom=5.D4


		
		real(kind=r8),parameter::back_ph_0=(6.D-5)*(1.D2/(2.D0*pi))
		real(kind=r8),parameter::adjust_gargett=1.D0
		real(kind=r8),parameter::back_k_0=(0.1D0)*(2.D0)*pi*(1.D-2)*adjust_gargett
		real(kind=r8),parameter::back_del_0=pi/back_k_0
		real(kind=r8),parameter::back_s2=back_ph_0*back_k_0
		real(kind=r8),parameter::back_sm2=1.D0/back_s2
		real(kind=r8),parameter::ri_internal=1.D0
		real(kind=r8),parameter::backfrac = 85.D-2,backfact = e**(-1)
		  
		real(kind=r8),parameter::v_back0=0.01,t_back0=0.01,s_back0=0.01
	  
		
  
		integer,parameter::n_theta_r_oct=INT(((pi/4.D0)*mt_ra_r)/15.D0)*15
end module constant

module param_mod
	    integer,parameter::mytid=3
	    integer,parameter::my_control=0,num_v_smooth_Ri=1    
end module param_mod

module global_g
  
      use kinds_mod
	 
  
      real(kind=r8)::g1,g2,g3,g4,g5,g6,g7,g8
      real(kind=r8)::d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
	 
	
	  
	  
      common/const0/g1,g2,g3,g4,g5,g6,g7,g8
      common/const1/d0,d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9
     

end module global_g
	  
module global_t

      use kinds_mod
	  
	  real(kind=r8)::ttot,tcot,tctot,tptot,tpcot,tpvot
	  real(kind=r8)::rit,ric
	  
	  
	  common /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot
      common /bb/rit,ric
end module global_t


module vmix_canuto
	 
    use kinds_mod            
    use blocks
    use distribution
    use domain_size
    use domain
    use constants
    use grid
    use broadcast
    use io
    use state_mod
    use exit_mod
    use sw_absorption
    use tavg, only: define_tavg_field, accumulate_tavg_field, accumulate_tavg_now
    use io_types, only: stdout
    use communicate, only: my_task, master_task
    use niw_mixing
    use tidal_mixing, only: TIDAL_COEF, tidal_mix_max, ltidal_mixing
    use registry
    use prognostic
    use time_management
	use param_mod
	use arguments 
#ifdef NO_R16
    integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
    integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif


  
	  
	
    ! !PUBLIC MEMBER FUNCTIONS:
      public :: vmix_coeffs_canuto, &
                fct_sal, &
                fctysal,     &
                eplatidepend_,interp_mix,calRid,calRi,calAN2,&
			    calculate_s2,nasagiss_canuto,formld,ccoeff,ourl2,fct,mksmsh,&
				oursal2_2a,smshsc_a3,rtwi,oursal2_zeroshear,quad_sal_pureconv,&
				mksmshss1,interp1d_expabs,interp2d_expabs
	  	   

	  
	  
	  
!***********************************************************************

    contains
	
!***********************************************************************

	
	subroutine vmix_coeffs_canuto(VDC, VVC,&
	                                TMIX,UMIX, VMIX,RHOMIX,&
									STF,SHF_QSW, &
									this_block,&
	                                SMFT,SMF)


	  
	  
	  implicit none 

	  
	  real(kind=r8),dimension(nx_block,ny_block,km),intent(in) :: &
      UMIX,           &
	  VMIX,           &! U,V     at mix time level
      RHOMIX         ! density at mix time level
	  
	  real(kind=r8),dimension(nx_block,ny_block,km,2),intent(in):: &
	  TMIX         	  ! tracers at mix time level, 1 for T,2 for S
	   
	  real(kind=r8):: TMIX_(nx_block,ny_block,km,2) ! For Salt Consistent with licom
	   
	   
	  real(kind=r8),dimension(nx_block,ny_block,2), intent(in) :: &
      STF                 ! surface forcing for all tracers

      real(kind=r8),dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW             ! short-wave forcing 
	   
	  real(kind=r8),dimension(nx_block,ny_block,2),intent(in), optional :: &
      SMF,               &! surface momentum forcing at U points
      SMFT                ! surface momentum forcing at T points
	   
	   
	  type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

      real(kind=r8), dimension(nx_block,ny_block,km),intent(inout) ::      &
      VVC      ! viscosity for momentum diffusion

      real(kind=r8), dimension(nx_block,ny_block,km,2),intent(inout) :: &
      VDC      ! diffusivity for tracer diffusion 1 for t ;2 for s
	  
!local variables

      integer:: ii,jj,kk,na,nmax,n,bid,k,nml_error
	  real(kind=r8)::Coriol,ustar_,buoysol,buoytur,amld0
	  real(kind=r8)::fricmx,wndmix!increase 
      real(kind=r8)::akm(km),akh(km),aks(km),z(km),an2(km),s2(km),ri(km),rid(km),&
	                                            t0(km),s0(km), rh0(km),&
												t(km),s(km),rh(km)
      real(kind=r8)::v_back(km),t_back(km),s_back(km)
	  real(kind=r8),dimension(nx_block,ny_block):: &
	  USTAR,             &! surface friction velocity
	  work,              &! temp array
	  BO,                &! surface buoyancy forcing = buoytur
      BOSOL,             &! radiative buoyancy forcing =buoysol
      TALPHA,            &! temperature expansion coeff 
      SBETA,             &! salinity    expansion coeff
      RHO1,              &! density at the surface
	  amld                !  mixing layer depth
	  real(kind=r8),dimension(nx_block,ny_block,km):: &
	  VSHEAR ,           & ! s~2
	  RID_XY,            &
	  AN2_XY
	  
	  real(kind=r8),parameter:: RH00=1.026    
	  real (kind=r8), parameter :: convect_visc = 1.0d+4 
      real (kind=r8), parameter :: convect_diff = 1.0d+4 
      character (char_len) :: bid_string
      character (char_len) :: bid_fi_string
      character (char_len) :: bid_fi_r_string

      integer :: i_min, j_min, k_min
      integer :: i_max, j_max, k_max
      integer :: k_111, k_222
	  
	  namelist /vmix_canuto_nml/fricmx, wndmix          ! &
                          ! bckgrnd_vdc_eq, bckgrnd_vdc_psim,     &
                          ! bckgrnd_vdc_ban,                      &
                          ! bckgrnd_vdc_dpth, bckgrnd_vdc_linv,   &
                          ! Prandtl, rich_mix,                    &
                          ! num_v_smooth_Ri, lrich, ldbl_diff,    &
                          ! lshort_wave, lcheckekmo,              &
                          ! lhoriz_varying_bckgrnd, llangmuir,    &
                          ! linertial
	  
	  
!***********************************************************************  
	write(bid_string,*) my_task 
    bid_fi_string    =   "diag_fi_"//trim(adjustl(bid_string))//".txt"
    bid_fi_r_string  =   "turb_nd2on2_"//trim(adjustl(bid_string))//".txt"
    bid = this_block%local_id

!*********************************************************************** 	
!                            initialization
!*********************************************************************** 
	TMIX_(:,:,:,1)=TMIX(:,:,:,1)
	TMIX_(:,:,:,2)=(TMIX(:,:,:,2)*1000.D0-35.D0)/1000.D0  !
    v_back=0.D0
	t_back=0.D0
	s_back=0.D0
	vvc   =0.D0
    vdc   =0.D0
	  
	  
	open (nml_in, file=nml_filename, status='old',iostat=nml_error)
    if (nml_error /= 0) then
         nml_error = -1
    else
         nml_error =  1
    end if
    do while (nml_error > 0)
        read(nml_in, nml=vmix_canuto_nml, iostat=nml_error)
    end do
    if (nml_error == 0) close(nml_in)
	  
!*********************************************************************** 
	
	call calculate_s2(VSHEAR,bid,UMIX,VMIX)
	
	
	 
	if (present(SMFT)) then
		USTAR = sqrt(sqrt(SMFT(:,:,1)**2 + SMFT(:,:,2)**2)/RH00)
	else
		WORK = sqrt(sqrt(SMF(:,:,1)**2 + SMF(:,:,2)**2)/RH00)
		call ugrid_to_tgrid(USTAR,WORK,bid)
	endif
	!-----------------------------------------------------------------------
!
!  compute density and expansion coefficients at surface:369

!
!-----------------------------------------------------------------------

    WORK = merge(-c2,TMIX(:,:,1,1),TMIX(:,:,1,1) < -c2)

	call state(1,1,WORK,TMIX(:,:,1,2),this_block, &
                  RHOFULL=RHO1, DRHODT=TALPHA, DRHODS=SBETA)

!-----------------------------------------------------------------------
!
!  compute turbulent and radiative sfc buoyancy forcing
!
!-----------------------------------------------------------------------

	do jj=1,ny_block
		do ii=1,nx_block
			if (RHO1(ii,jj) /= c0) then
				BO   (ii,jj) = grav*(-TALPHA(ii,jj)*STF(ii,jj,1) - &
									SBETA (ii,jj)*STF(ii,jj,2))/RHO1(ii,jj)

				BOSOL(ii,jj) = -grav*TALPHA(ii,jj)*SHF_QSW(ii,jj)/RHO1(ii,jj)
			else
				BO   (ii,jj) = c0
				BOSOL(ii,jj) = c0
			endif
		end do
	end do
   
    call calAN2(an2_xy,TMIX(:,:,:,1),TMIX(:,:,:,2),this_block)
    call calRid(rid_xy,TMIX(:,:,:,1),TMIX(:,:,:,2),this_block) 
	
  
	 

	
	
	do ii=1,nx_block
		do jj=1,ny_block
	  	  
			if(kmt(ii,jj,bid).gt.1)then
	   
				na=kmt(ii,jj,bid)-1
				n=na
		 
				t0=TMIX_(ii,jj,:,1)
				call interp_mix(t,t0,na)
		 
				do k=1,na+1
					z(k) = zw(k) 
				end do
		 
		 
		 
	    
				if(n.eq.0) then
					amld(ii,jj)=z(1)
		 
				else
					call formld(z,t,amld(ii,jj),n)
				endif
	   
	   
	   
			else
				amld(ii,jj)=0.D0
	   
			end if
	  
		end do
	end do
	
	
	if(if_smooth_amld.EQ.1)then
		call smooth_amld(amld,bid)
	end if
	
	
	
	
	
	
	
	
	do ii=1,nx_block
	   do jj=1,ny_block
			if (kmt(ii,jj,bid) .gt. 1) then
				na=kmt(ii,jj,bid)-1
				Coriol=fcort(ii,jj,bid)
				nmax=km
				n=na
				t0=TMIX_(ii,jj,:,1)
				s0=TMIX_(ii,jj,:,2)
				rh0=RHOMIX(ii,jj,:)
		 
				buoytur=BO(ii,jj)
				buoysol=BOSOL(ii,jj)
				ustar_=USTAR(ii,jj)
				s2=VSHEAR(ii,jj,:)		 
				an2=AN2_XY(ii,jj,:)
		 
				do k=1,km
					if(abs(VSHEAR(ii,jj,k)) .GT. 1.0d-25 ) then
						rid(k)=RID_XY(ii,jj,k)/VSHEAR(ii,jj,k)
					else
						rid(k)=RID_XY(ii,jj,k)/1.0d-25
					end if
				end do


		 
				do k=1,na+1
					z(k) = zw(k) 
				end do
		 
		 
				call calRi(ri,an2,s2)
		 
		 
		 
				if(ifsmooth_ri .EQ.1)then
		 
					call smooth_ri(ri,num_v_smooth_Ri,na)
					call smooth_ri(rid,num_v_smooth_Ri,na)
				end if
		 
				call interp_mix(t,t0,na)
				call interp_mix(s,s0,na)
				call interp_mix(rh,rh0,na)
	   

		
				call nasagiss_canuto(z(1:na+1),t(1:na),s(1:na),rh(1:na),ri(1:na),rid(1:na),s2(1:na), &
					v_back,t_back,s_back, &
					an2(1:na),&					!030401 N^2
					ustar_,buoytur,buoysol,&		!020912 surface fluxes
                    Coriol,		&		!030424 f_coriolis
					amld(ii,jj),akm,akh,aks,n,na,nmax,fricmx,wndmix,bid_fi_string,bid_fi_r_string)
                    
		
				if (if_convect) then 				
					do k = 1, na
						if ((an2(k) .lt. 0.0d0 ).and.(z(k).GT.amld(ii,jj))) then 
							akm(k) = akm(k)+convect_visc
							akh(k) = akh(k)+convect_diff
							aks(k) = aks(k)+convect_diff
						end if
					end do
				end if


		
					
				IF(my_control.EQ.0) THEN
		
					vvc(ii,jj,:)=akm
					vdc(ii,jj,:,1)=akh
					vdc(ii,jj,:,2)=aks
		
				else if (my_control.EQ.1) then
		
					do k=1,na
		
						vvc(ii,jj,k)=akm(k)+1.D-1
						vdc(ii,jj,k,1)=akh(k)+2.D-2
						vdc(ii,jj,k,2)=aks(k)+2.D-2
			
					end do
				ELSE
					vvc(ii,jj,:)=akm
					vdc(ii,jj,:,1)=akh
					vdc(ii,jj,:,2)=aks
		
					do k=1,na
		
						vvc(ii,jj,k)=max(vvc(ii,jj,k),1.D-1)
						vdc(ii,jj,k,1)=max(vdc(ii,jj,k,1),2.D-2)
						vdc(ii,jj,k,2)=max(vdc(ii,jj,k,2),2.D-2)
		 
					end do
			
				ENDIF
		
		
				call area_smooth(60.0_r8,87.0_r8,300.0_r8,360.0_r8,10.0d0,10.0d0,na,bid,umix(ii,jj,:),vmix(ii,jj,:),vvc(ii,jj,:),vdc(ii,jj,:,:),ii,jj)
				call area_smooth(60.0_r8,87.0_r8,0.0_r8,30.0_r8,10.0d0,10.0d0,na,bid,umix(ii,jj,:),vmix(ii,jj,:),vvc(ii,jj,:),vdc(ii,jj,:,:),ii,jj)
				call area_smooth(-80.0_r8,-65.0_r8,-1.0_r8,361.0_r8,10.0d0,10.0d0,na,bid,umix(ii,jj,:),vmix(ii,jj,:),vvc(ii,jj,:),vdc(ii,jj,:,:),ii,jj)
				call area_smooth(70.0_r8,85.0_r8,230.0_r8,360.0_r8,10.0d0,10.0d0,na,bid,umix(ii,jj,:),vmix(ii,jj,:),vvc(ii,jj,:),vdc(ii,jj,:,:),ii,jj)
					
		
	  
			end if
	   end do
    end do
	

    end subroutine vmix_coeffs_canuto
	
	subroutine area_smooth(LAT_MIN,LAT_MAX,LON_MIN,LON_MAX,UMAX,VMAX,NA,BID,UMIX,VMIX,VVC,VDC,II,JJ)
	implicit none
	integer::K
	integer,intent(in)::NA,BID,II,JJ
	real(kind=r8),intent(in)::LAT_MAX,LAT_MIN,LON_MAX,LON_MIN,UMAX,VMAX
	real(kind=r8), dimension(km),intent(in)::UMIX,VMIX
	real(kind=r8), dimension(km),intent(inout) :: VVC
	real(kind=r8), dimension(km,2),intent(inout) :: VDC
	real (kind=r8), parameter :: convect_visc = 1.0d+4
	  
	if ((TLATD(ii,jj,bid) .lt. LAT_MAX)  .and. (TLATD(ii,jj,bid) .gt. LAT_MIN)  .and.  & 
           (TLOND(ii,jj,bid) .ge. LON_MIN) .and. (TLOND(ii,jj,bid) .le. LON_MAX)) then
        do K = 1, NA
            if(UMIX(k) .ge. UMAX .or. VMIX(k) .ge. VMAX) then 
                if (k .eq. 1) then
                    vvc(k)  =   vvc(k)+convect_visc
                    vdc(k,1)=   vdc(k,1)+convect_visc
                    vdc(k,2)=   vdc(k,2)+convect_visc
                end if
                if (k .le. na .and. k .gt. 1) then
                    vvc(k-1:k)  =   vvc(k-1:k)+convect_visc
                    vdc(k-1:k,1)=   vdc(k-1:k,1)+convect_visc
                    vdc(k-1:k,2)=   vdc(k-1:k,2)+convect_visc
                end if
            end if
        end do
    end if
	
	
	
	end subroutine area_smooth
	
	
	
	
	
	subroutine smooth_amld(amld,bid)
	implicit none
	integer,intent(in)::bid
	real(kind=r8)::amld(nx_block,ny_block),work1(nx_block,ny_block)
	real(kind=r8)::cw,ce,cn,cc,cs,c0
	integer::ii,jj
	
	work1=amld
    do jj=2,ny_block-1
      do ii=2,nx_block-1
        if ( KMT(ii,jj,bid) /= 0 ) then
         cw = 1.25D-1
         ce = 1.25D-1
         cn = 1.25D-1
         cs = 1.25D-1
         cc = 5.0D-1
		 c0 = 0.D0
         if ( KMT(ii-1,jj,bid) == 0 ) then
           cc = cc + cw
           cw = c0
         endif
         if ( KMT(ii+1,jj,bid) == 0 ) then
           cc = cc + ce
           ce = c0
         endif
         if ( KMT(ii,jj-1,bid) == 0 ) then
           cc = cc + cs
           cs = c0
         endif
         if ( KMT(ii,jj+1,bid) == 0 ) then
           cc = cc + cn
           cn = c0
         endif
         amld(ii,jj) =  cw * work1(ii-1,jj)   &
                     + ce * work1(ii+1,jj)   &
                     + cs * work1(ii,jj-1)   &
                     + cn * work1(ii,jj+1)   &
                     + cc * work1(ii,jj)
       endif
     enddo
   enddo
	
	
	
	
	end subroutine smooth_amld
	
	
	
	subroutine smooth_ri(rii,num_v_smooth_Ri,na)
	
	implicit none
	
	integer,intent(in)::num_v_smooth_Ri,na
	real(kind=r8)::rii(na)
	integer::n,k
	real(kind=r8)::FRI1
	
	
	do n = 1,num_v_smooth_Ri
 
      FRI1            =  0.25D0 * rii(1)
     
 
     
            if (na >= 3) then
            			   
                do k=1,na-1
      
			      rii(k) = FRI1 + 0.5D0*rii(k) + 0.25D0*rii(k+1)
           
                  FRI1 = 0.25D0*rii(k)
                end do
				
			    rii(na)=0.25D0*rii(na-1)+0.75D0*rii(na) 
			 
		    endif
       

     enddo

	
	
	
	
	
	end subroutine smooth_ri
	
	
	
	
	subroutine interp_mix(var_new,var_old,na)
	implicit none
	real(kind=r8),intent(in)::var_old(km)
	integer,intent(in)::na
	   
	real(kind=r8),intent(out)::var_new(km)
	   
	integer::k
	   
	   
	do k=1,na
	   var_new(k)=var_old(k)+(var_old(k+1)-var_old(k))*(dz(k)/(dz(k+1)+dz(k)))
   end do 
	   
	end subroutine interp_mix
	
	
	subroutine calRid(rid_xy,t,s,this_block)
	implicit none
	
    real(kind=r8),dimension(nx_block,ny_block,km),intent(in) :: &
	t,&
	s

    type (block), intent(in) :: this_block    

! !OUTPUT PARAMETERS:  
    real (kind=r8), intent(out) :: rid_xy(nx_block,ny_block,km)
	 
	 
	 
	integer :: k,na,ii,jj             
  

    real(kind=r8),dimension(nx_block,ny_block) :: &
     
      TALPHAS,            &! temperature expansion coefficient on up level
      SBETAS,             &! salinity    expansion coefficient on up level
	  TALPHAD,            &! temperature expansion coefficient on next level
      SBETAD              ! salinity    expansion coefficient on next level

	  
    real(kind=r8),parameter:: RH0=1.026    
	
	
	
		    
	call state(1,1,t(:,:,1),s(:,:,1),&
							this_block,DRHODT=TALPHAS,DRHODS=SBETAS )	    
			
	do k=2,km							
		call state(k,k,t(:,:,k),s(:,:,k),&
							this_block,DRHODT=TALPHAD,DRHODS=SBETAD )					
							
	    do ii=1,nx_block
			do jj= 1,ny_block	
				rid_xy(ii,jj,k-1)= (p5*(TALPHAS(ii,jj)+TALPHAD(ii,jj))*(t(ii,jj,k)-t(ii,jj,k-1))/dzw(k-1)+ &
				                   p5*(SBETAS(ii,jj)+SBETAD(ii,jj))*(s(ii,jj,k-1)-s(ii,jj,k))/dzw(k-1))/RH0
			end do
		end do 
		
		TALPHAS=TALPHAD
		SBETAS=SBETAD                !conserve  alpha,beta into up level for next loop 
	end do
	
	
	
	
	
	
	end subroutine calRid
	
	  
	subroutine calRi(ri,an2,s2)
    implicit none
    real(kind=r8),intent(in)::an2(km),s2(km)
	real(kind=r8),intent(out)::ri(km)
	integer::k
	   
	do k=1,km
		if(abs(s2(k)) .ge. 1.0d-25) then
			ri(k)=an2(k)/s2(k)
		else
			ri(k)=an2(k)/1.0d-25
		end if
	end do
	
	end subroutine calRi	 
	  
	  
	  
	  
	  
	subroutine calAN2(an2_xy,t,s,this_block)

! !DESCRIPTION:
!  This routine calculates the buoyancy differences at model levels.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
     implicit none
	 


 !    integer,intent(in):: na
     real (kind=r8), dimension(nx_block,ny_block,km),intent(in) :: t,s

	 type (block), intent(in) :: &
     this_block          ! block information for current block

! !OUTPUT PARAMETERS:  
     real (r8), dimension(nx_block,ny_block,km), intent(out) :: an2_xy


!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) :: &
      k,                 &! vertical level index
      ii,jj             ! horizontal indices

    integer :: bid
  

    real (r8), dimension(nx_block,ny_block) :: &
      RHOD,              &! density of sfc t,s displaced to k
      RHOS,              &! density of t(k-1),s(k-1) displaced to k
    !  RHOK,              &! density at level k
   !   TEMPSFC,           &! adjusted temperature at surface
      TALPHA,            &! temperature expansion coefficient
      SBETA               ! salinity    expansion coefficient


	  
    real(kind=r8),parameter:: RH0=1.026    !






!-----------------------------------------------------------------------
!
!  calculate DBLOC and DBSFC for all other levels
!
!-----------------------------------------------------------------------
    bid = this_block%local_id
    an2_xy(:,:,:) = 0.0d0

	do k = 2,km

	!			TEMPK(:,:,klvl) = merge(-c2,t(k),t(k) < -c2)

		call state(k-1,k,t(:,:,k-1),s(:,:,k-1),&
	                    this_block,RHOFULL=RHOS )
		call state(k,k,t(:,:,k),s(:,:,k),&
	                    this_block,RHOFULL=RHOD)	
        do ii=1,nx_block
			do jj=1,ny_block
                if (k .le. kmt(ii,jj,bid)) then
				    an2_xy(ii,jj,k-1)=grav*(RHOD(ii,jj)-RHOS(ii,jj))/(DZW(k-1)*RH0)
                end if
				
		    end do
		end do
		
	end do
	

 end subroutine calAN2 
	  
     subroutine calculate_s2(VSHEAR,bid,UMIX,VMIX)
!	    use domain_size only:km
!		use grid only:ugrid_to_tgrid
		
		
	    integer,intent(in)::bid
        real(kind=r8),intent(in)::UMIX(nx_block,ny_block,km),VMIX(nx_block,ny_block,km)
!		type(block),intent(in)::this_block
		real(kind=r8),intent(out)::VSHEAR(nx_block,ny_block,km)
		
		integer::k
		real(kind=r8)::FRI(nx_block,ny_block,km)
	  
	  
	  
	  
	  
	     do k = 1,km

!-----------------------------------------------------------------------
!
!     compute velocity shear squared and average to T points:
!     VSHEAR = (UUU(k)-UUU(k+1))**2+(VVV(k)-VVV(k+1))**2
!     Use FRI here as a temporary.
!
!-----------------------------------------------------------------------

			if (k < km) then

				FRI(:,:,k) = ((UMIX(:,:,k)-UMIX(:,:,k+1))**2 + &
						(VMIX(:,:,k)-VMIX(:,:,k+1))**2)/(dzw(k)**2)

     

				call ugrid_to_tgrid(VSHEAR,FRI,bid)

			else

				VSHEAR = c0

			endif
	  
	   end do
	  
 end subroutine calculate_s2


      ! subroutine turb_2(z,t,s,rh,ri,rid,s2, &
					! fricmx,wndmix,v_back,t_back,s_back, &
					! an2,&					!030401 N^2
					! ustar_,buoytur,buoysol,&		!020912 surface fluxes
                    ! Coriol,		&		!030424 f_coriolis
                    ! amld,akm,akh,aks,n,na,nmax,&
					! isurfuse,			&	!020912 choice of use of surface fluxes
      		        ! ifextermld,ifoutput,ii,jj)
    subroutine nasagiss_canuto(z,t,s,rh,ri,rid,s2, &
					v_back,t_back,s_back, &
					an2,&					!030401 N^2
					ustar_,buoytur,buoysol,&		!020912 surface fluxes
                    Coriol,		&		!030424 f_coriolis
                    amld,akm,akh,aks,n,na,nmax,fricmx,wndmix,bid_fi_string,bid_fi_r_string)
      !use precision_mod
	!  use precision_
      use param_mod, only: mytid
	  use arguments
	  use constant
!	  use global
	  
      implicit real(kind=r8) (a-h,o-z)
!	  integer,intent(in)::ifextermld,ifoutput,isurfuse
	  integer,intent(in)::n,na,nmax
	  real(kind=r8),intent(in)::ustar_,buoytur,buoysol,Coriol,amld
	  real(kind=r8),intent(out)::akm(nmax),akh(nmax),aks(nmax)
	  
	  
	  
	  real(kind=r8)::fricmx,wndmix
	  real(kind=r8)::z(na+1),t(na),s(na),rh(na),ri(na),rid(na),s2(na),v_back(nmax),t_back(nmax),s_back(nmax),an2(na),epsy(na)
	 
!	  integer::ii,jj

      real(kind=r8),save::epson2,epson2_,epson2_bot,eps_bot,eps_bot__under
	  real(kind=r8),save::and2on2a1(-mt:mt),amtaun2a1(-mt:mt),dand2on2,sma1(-mt:mt),sha1(-mt:mt),ssa1(-mt:mt),rri,rnd2on2
	  real(kind=r8),save::aldeep(nbig),ria(ntbl),slq2a(ntbl),sma(ntbl),sha(ntbl)
	  real(kind=r8),save::rib(-mt:mt),ridb(-mt:mt), slq2b(-mt:mt,-mt:mt),smb(-mt:mt,-mt:mt),shb(-mt:mt,-mt:mt),ssb(-mt:mt,-mt:mt)
	 
	  
	  
	  
	  
	  
	  real(kind=r8),save::sisamax(-n_theta_r_oct:3*n_theta_r_oct)
	  real(kind=r8),save::ra_rmax(-n_theta_r_oct:3*n_theta_r_oct)
	  real(kind=r8),save::ss_r1(-n_theta_r_oct:3*n_theta_r_oct)
	  real(kind=r8),save:: c_y_r0(-n_theta_r_oct:3*n_theta_r_oct) 
	  real(kind=r8),save::back_ra_r(-n_theta_r_oct:3*n_theta_r_oct)
      real(kind=r8),save::sm_r1(-n_theta_r_oct:3*n_theta_r_oct)
      real(kind=r8),save::sh_r1(-n_theta_r_oct:3*n_theta_r_oct)
      real(kind=r8),save::slq2_r1(-n_theta_r_oct:3*n_theta_r_oct)
	  
	  
	  
	  integer::ifirst=0
	  integer,save::irimax(-mt:mt)
	  logical::lifupper,lifepsy(na)
	  complex(kind=r16)::zlomega

	  integer::irisign,iridstep,irid,iri,iristep,idifs,ind2on2sign,ind2on2step,ind2on2,idif
	  integer::nb,mt0s,mtm1s
	  
	  real(kind=r8),save::visc_cbu_limit,diff_cbt_limit,ako,b1,slq2b_00,smb_00,shb_00,ssb_00, c_y0,c_y00
	  real(kind=r8),save::back_l_0,rimax,dri,deltanum,deltaden,delta,rrcrn,rrcrp,theta_rcrn,theta_rcrp,theta_rcrn_deg,theta_rcrp_deg,deltheta_r,delra_r
	  character (char_len),intent(in) :: bid_fi_string,bid_fi_r_string

 !     data ifirst/0/



     
    IF(nmax.GT.nbig) THEN
		if (mytid.eq.1) then
			WRITE(10,*) " "
			WRITE(10,*) " "
			WRITE(10,*) "****************************"
			WRITE(10,*) "*PROBLEM IN TURB_2 ROUTINE.*"
			WRITE(10,*) "Number of model levels exceeds nbig."
			WRITE(10,*) "nbig=",nbig,"	nmax=",nmax
			WRITE(10,*) "If want to use this many levels"// &
             " increase nbig to be bigger than nmax."
			WRITE(10,*) "Program is stopping now."
		endif
        STOP
    END IF
!*****C
    nb=MAX(0,n)
if(ifirst.eq.0) then
!000107	Add nmodel to writeout. nmodel=1,2 for 1,2 pt. closure.
	if (mytid.eq.0) then
		WRITE(10,*) " "
		WRITE(10,*) " "
		WRITE(10,*) "************************************"
		WRITE(10,*) "************************************"
		WRITE(10,*) &
				   "Turbulence calculated by turb_2gi1b 040217 version."
		WRITE(10,*) "************************************"
		WRITE(10,*) "nmodel=",nmodel
		WRITE(10,*) "************************************"
		WRITE(10,*) "************************************"
		WRITE(10,*) " "
		WRITE(10,*) " "
	endif
!*****C
!030401Y Take zero shear unstable case *foreground* calculation from my code for HYCOM
!Y	 in inigiss_fixed2.fs0 . [See NBp.030401-2to3.] 
!Y	 When the Richardson number is more negative than the most negative table value
!Y	 it is more accurate to use the zero shear approximation derived from Canuto's
!Y	 analytic pure convection formula [See NBp.030326-3to9.].
!030324-27AH Amend to make 1D table vs. (N_d^2/N^2) for zero shear unstable case.
!            Choose (N_d^2/N^2) table values to be the same as Ri_d table values.

!980501 Maximum diffusivity, fricmx, and surface minimum, wndmix, set externally
!       fricmx=80.
    visc_cbu_limit=fricmx
    diff_cbt_limit=fricmx
!       wndmix=10.
!990202 *SET THE VALUE OF THE KOLMOGOROV CONSTANT, K_0, HERE NAMED "ako".*
	ako = 1.6D0

!*****C
!991107 START OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
	IF(ifsali.EQ.1) THEN                                                   !go to 10 
!990202-26 Calculate constant lengthscale for the background for ifsalback=3,4,5
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
        dri = -ri0/DFLOAT(mt0)
		IF(nmodel.EQ.1) THEN  
!981104	Need to pass back the value of B_1 from oursal2 for use here. 
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (Ri=0,Rid=0).
			CALL OURSAL2_2A(b1,0.D0,0.D0,slq2b_00, &
                     smb_00,shb_00,ssb_00, &
                     c_y0,c_y00,0,0)
!			write(10,*) "OURSAL2_2A FINISHED AT 818"

	    END IF
		back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0	
		IF(ifsalback.EQ.3) THEN
			if (mytid.eq.0) then
				WRITE(10,*) " "
				WRITE(10,*) "************************************"
				WRITE(10,*) "Internal wave constants for background."
!990303 Add write-out of residual constant backgrounds.   
				WRITE(10,*) "Residual Constant Background Diffusivities:"
				WRITE(10,*) "K_M /(cm^2 sec^{-1})",v_back0
				WRITE(10,*) "K_H /(cm^2 sec^{-1})",t_back0
				WRITE(10,*) "K_S /(cm^2 sec^{-1})",s_back0
				WRITE(10,*) "."
				WRITE(10,*) " "
!*****C
				WRITE(10,*) "Shear^2/(sec^{-2}) =",back_s2
				WRITE(10,*) "Lengthscale, del_0/(cm) =",back_del_0
				WRITE(10,*) "Lengthscale, l_0/(cm) =",back_l_0
				WRITE(10,*) '"adjust_gargett="',adjust_gargett
				WRITE(10,*) "************************************"
				WRITE(10,*) " "
			endif
		ELSE IF(ifsalback.GE.4) THEN
			if (mytid.eq.0) then
				WRITE(10,*) " "
				WRITE(10,*) "************************************"
				WRITE(10,*) "Dubovikov Internal wave constants for background."
			end if
			IF(ifsalback.EQ.4) THEN	
				if (mytid.eq.0) then
					WRITE(10,*) "Internal wave Richardson number=",ri_internal
				end if
			ELSE IF(ifsalback.EQ.5) THEN
!990301
				if (mytid.eq.0) then
					WRITE(10,*) "Ratio of Background to Critical ra_r"// &
                 " [\\equiv ({Ri_T}^2 + {Ri_C}^2)^(1/2)]",backfrac
				end if
			ELSE IF(ifsalback.EQ.6) THEN
!990303-04
				if (mytid.eq.0) then
					WRITE(10,*) "Ratio of Background dimensionless K_M ="// &
					" S_M/(S l/q) to its value at Ri_T=Ri_C=0", &
					backfact
				endif
			END IF
			if (mytid.eq.0) then
				WRITE(10,*) "Lengthscale, del_0/(cm) =",back_del_0
				WRITE(10,*) "Lengthscale, l_0/(cm) =",back_l_0
				WRITE(10,*) "************************************"
				WRITE(10,*) " "
			endif
		END IF
!*****C
!		IF(ifsali.EQ.1) dri = -ri0/DFLOAT(mt0)
!*****C
!** BUILD SALINITY MODEL TABLES VS. "Ri = Ri_T + Ri_C" AND "Ri_d = Ri_T - Ri_C".
!981027 **Use separate loops for calculation of independent table variables.**
		DO iridsign=0,1
			iridstep=(-1)**iridsign 
			DO irid= 0,mt*iridstep,iridstep
!981019 Set Ri_d table values. (See NBP59,63=p#A27,30.)
				IF(ABS(irid).LE.mt0) THEN
					ridb(irid) = DFLOAT(irid)*dri
				ELSE
					mt0s = mt0*iridstep
					mtm1s = (mt0-1)*iridstep
!000320 1st day of spring introduction of exponential absolute val table option.
					IF(ifexpabstable.EQ.0) THEN
						idifs = (ABS(irid)-mt0)*iridstep
						ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
						ridb(mtm1s))**(idifs**2))
					ELSE IF(ifexpabstable.EQ.1) THEN
						idif = ABS(irid)-mt0
						ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
						ridb(mtm1s))**(idif))
					END IF
				END IF
!*****C
			END DO 
		END DO
		DO irisign=0,1
			iristep=(-1)**irisign 
			DO iri= 0,mt*iristep,iristep
!981019 Set Ri table values. (See NBP59,63=p#A27,30.)
				IF(ABS(iri).LE.mt0) THEN
					rib(iri) = DFLOAT(iri)*dri
				ELSE
					mt0s = mt0*iristep
					mtm1s = (mt0-1)*iristep
!000320 1	st day of spring introduction of exponential absolute val table option.
					IF(ifexpabstable.EQ.0) THEN
						idifs = (ABS(iri)-mt0)*iristep
						rib(iri) = rib(mt0s)*((rib(mt0s)/ &
						rib(mtm1s))**(idifs**2))
					ELSE IF(ifexpabstable.EQ.1) THEN
						idif = ABS(iri)-mt0
						rib(iri) = rib(mt0s)*((rib(mt0s)/ &
						rib(mtm1s))**(idif))
					END IF
				END IF
!*****C
			END DO 
		END DO
!*****C
!011107yXI ***If using interp2d_expabs introduce ratio between adjacent Richardson numbers in nonlinear part of table.***
!030327 rnd2on2 is the ratio between adjacent N_d^2/N^2 in the nonlinear part of the zero shear 1D table.

!yXI

!        write(10,*) " FINISHED at 940"
		DO iridsign=0,1
			iridstep=(-1)**iridsign 
			DO irid= 0,mt*iridstep,iridstep
			loop1:DO irisign=0,1
					iristep=(-1)**irisign
					DO iri= 0,mt*iristep,iristep
						IF(nmodel.EQ.1) THEN
!981104	Need to pass back the value of B_1 from oursal2 for use here. 
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
							CALL OURSAL2_2A(b1,rib(iri),ridb(irid),slq2b(iri,irid), &
										smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
										c_y0,c_y00,iri,irid)							
			
						END IF
						IF(slq2b(iri,irid).LT.0) THEN
							irimax(irid) = iri - 1 
							cycle loop1
						END IF
					END DO
!   15     	    CONTINUE
				END DO loop1
			END DO	
		END DO
		
!		write(10,*) "OURSAL2_2A FINISHED AT 953"
!**
!030421Z Only calculate the zero shear table when the zero shear option is enabled.
      	IF(ifastexpabs.EQ.1) THEN
			rri = rib(mt0)/rib(mt0-1)
!030424 Only calculate rnd2on2 when zero shear parameterization is enabled.
			IF(ifzeroshear) rnd2on2 = rri 	
		END IF
      


		IF(ifzeroshear) THEN
!030401-07Y Calculation of table for zero shear approximation from my HYCOM inigiss_fixed2.fs0 .
!030324-26AH Make 1 dimensional table of turbulence functions of N_d^2/N^2
!         to be used for the unstable case with negligible shear.
!         N_d^2 \equiv N_Heat^2 - N_Salt^2 . N_d^2/N^2 ranges from - to + infinity.
!         N_d^2 is analogous to Ri_d^2, so I try making its table values exactly the same.
!         oursal2_zeroshear assumes Shear^2=0 and N^2 < 0.
	
			dand2on2 = dri
! --- Set N_d^2/N^2 table values.
!07Y	Calculate moving out from zero index first postive indices and then negative ones.[See NBp.030407-08.]
			DO ind2on2sign=0,1
				ind2on2step=(-1)**ind2on2sign
				DO ind2on2 = 0,mt*ind2on2step,ind2on2step
					IF(ABS(ind2on2).LE.mt0) THEN
						and2on2a1(ind2on2) = DFLOAT(ind2on2)*dand2on2
					ELSE
						mt0s  = SIGN(mt0,ind2on2)
						mtm1s = SIGN(mt0-1,ind2on2)
! --- introduction of exponential absolute val table option.
						IF(ifexpabstable.EQ.0) THEN
							idifs = SIGN((ABS(ind2on2)-mt0),ind2on2)
							and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
														and2on2a1(mtm1s))**(idifs**2))
						ELSE IF(ifexpabstable.EQ.1) THEN
							idif = ABS(ind2on2)-mt0
							and2on2a1(ind2on2) = and2on2a1(mt0s)*((and2on2a1(mt0s)/ &
												and2on2a1(mtm1s))**(idif))
						END IF
					END IF
				END DO
			END DO
			DO ind2on2 = -mt,mt
				CALL oursal2_zeroshear &
					(and2on2a1(ind2on2),amtaun2a1(ind2on2) &
					,sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2))
			END DO
!******Y
           
		END IF
!******Z
!        write(10,*) "OURSAL2_zeroshear FINISHED at 1002"

!981215-31 Add writes in salinity model case.
		if (mytid.eq.0) then
			WRITE(10,*) " "
			WRITE(10,*) " "
			WRITE(10,*) " "
			WRITE(10,*) "************************************************"
			WRITE(10,*) " "
			WRITE(10,*) "New Temperature-Salinity Model"
			WRITE(10,*) " "
!030722Z1a Add writeout of isurfuse here (See NBp.030722-2 Vol.XIX.)
!020916X Switch for use of surface fluxes to dimensionalize turbulence model.
			WRITE(10,*) "************************************"
			WRITE(10,*) "isurfuse=",isurfuse
			WRITE(10,*) "************************************"
!******X
!*****CZ1a
			WRITE(10,*) "ifsali=",ifsali
			WRITE(10,*) "ifsalback=",ifsalback
!040217Zi1b Write out switch for use of minimum shear and amount if used.
			IF(.NOT.ifzeroshear) THEN
				WRITE(10,*) "ifshearmin=",ifshearmin
				IF(ifshearmin) WRITE(10,*) "s2min=",s2min
			END IF
!030424Z Write switch for whether zero shear model is enabled.
			WRITE(10,*) "ifzeroshear=",ifzeroshear
!030424Z Write parameters for rotation's effect on lengthscale.
			WRITE(10,*) "ilomega=",ilomega
			WRITE(10,*) "amldminlom=",amldminlom
!000215
!030722Z1a Add writeout of icondear here (See NBp.030722-2 Vol.XIX.)
!030721Z1a Write switch for Deardorff treatment.
			WRITE(10,*) " "
			WRITE(10,*) "icondear=",icondear
			IF(icondear.EQ.-1) THEN
				WRITE(10,*) "Do *not* use Deardorff lengthscale modification."
			ELSE IF(icondear.EQ.0) THEN
				WRITE(10,*) "Ye Cheng's old Deardorff:"// &
     "modify l and \tau N leaving S_X unmodified."
			ELSE IF(icondear.EQ.1) THEN
				WRITE(10,*) "Ye Cheng's new Deardorff:"// &
     "modify l but leave *both* \tau N and S_X unmodified."
			END IF
			WRITE(10,*) " "
!*****CZ1a
!*****CZ1a
			WRITE(10,*) "ifepson2=",ifepson2
!020219D Bottom enhancemant writes included.
			IF(ifepson2.GT.0) THEN 
				IF(ifepson2.EQ.2) WRITE(10,*) &
       "epsilon/(N^2) used even for strong mixing beneath weak mixing" 
!030722Z1a Fix writeout of ideeplat here (See NBp.030722-2 Vol.XIX.)
!030429-0502Z1 Write switch for latitude dependence of background mixing.
				WRITE(10,*) "ifdeeplat=",ifdeeplat
				IF(ifdeeplat.GT.0) THEN
					WRITE(10,*) &
        "Use latitude dependence of interior ocean value of \\epsilon/N^2"
!Zi1a',040422Zi1bj Floor on latitude dependent factor in background mixing
					WRITE(10,*) "eplatidependmin=",eplatidependmin
				END IF
!*****CZ1a
!020404D,030502Z1 Must write out the parameter epson2__ for reference \epsilon/N^2 
!Z1	since epson2_ and epson2 are no longer constants.
				WRITE(10,*) "epson2__=",epson2__
				WRITE(10,*) "ifbotenhance=",ifbotenhance
				IF(ifbotenhance.EQ.1) THEN
					WRITE(10,*) "eps_bot0=",eps_bot0
					WRITE(10,*) "scale_bot=",scale_bot
				END IF
			END IF
!*****CD
!000317
			WRITE(10,*)"ifrafgmax=",ifrafgmax
			WRITE(10,*) " "
			WRITE(10,*)"ifbg_theta_interp=",ifbg_theta_interp
			WRITE(10,*) " "
!011108yXI Write out whether using exponential absolute value spacing for nonlinear part of table
!          and whether using new interpolation tailored to exponential absolute value in that part.
			WRITE(10,*) " "
			WRITE(10,*) "***********************************************"
			WRITE(10,*) "ifexpabstable=",ifexpabstable
			WRITE(10,*) "ifastexpabs=",ifastexpabs
			WRITE(10,*) "***********************************************"
			WRITE(10,*) " "
!yXI
			WRITE(10,*) &
			"    i      ", &
			"    rib(i)      ","    ridb(i)     ", &
			"irimax(i)  "
			DO i= -mt,mt
				WRITE(10,9050) i,rib(i),ridb(i),irimax(i)
			END DO
!*****C
			WRITE(10,*) " "
			WRITE(10,*) "irid       Ri_d        Ri(irimax)  " &
					// "S_M        S_H        S_S        " &
					// "S_M/S_H    S_S/S_H    "
			DO irid= -mt,mt
				WRITE(10,9100) irid,ridb(irid),rib(irimax(irid)), &
					smb(irimax(irid),irid), &
					shb(irimax(irid),irid), &
					ssb(irimax(irid),irid), &
					smb(irimax(irid),irid)/shb(irimax(irid),irid), &
					ssb(irimax(irid),irid)/shb(irimax(irid),irid)
			END DO
		endif
!000316 CALCULATE "R_r_Critical" USING 'S 000228 ANALYTIC FORMULA
!	FOR "R_rho_Critical". See NBp.000229-3 and 000316-4.
!	R_rho_Canuto \equiv -Ri_C/Ri_T \equiv -R_r .
!	In a sheet dated 000228 Canuto gave me:
!	"R_\rho^{cr} = {1 \over \Deta} [1 {+\over-} \sqrt{1 - \Delta^2}] 
!	\Delta \equiv {{\pi_2(1 + {15 \over 7} \pi_3)} \over
!	{\pi_3 - \pi_2 + (15 \over 14} \pi_3^2}} ".
!	Note that the + and - choices are reciprocals so this covers
!	both the Salt Fingering and Double Diffusive Critical R_\rho's.
!	From Ocean Turbulence III paper have: 
!	\pi_{1,2,3,4,5} = 
!	(\tau_pc,\tau_c\theta,\tau_c,\tau_p\theta,\tau_\theta)/\tau 
!	R_r_Crit = [-1 -/+ \sqrt{1 - \Delta^2}]/Delta
!	\Delta = {{{\tau_c\theta \over \tau} ( 1 + (15/7)*{\tau_c \over \tau})}
!	        \over {{\tau_c \over \tau} - {\tau_c\theta \over \tau} + 
!	                (15/14) {\tau_c \over \tau}^2}}
		deltanum = tctot*(1.D0 + ((15.D0/7.D0)*tcot))
		deltaden = tcot - tctot + ((15.D0/14.D0)*(tcot**2))
		delta = deltanum/deltaden
		rrcrn = (-1.D0 - SQRT(1.D0 - (delta**2)))/delta
		rrcrp = (-1.D0 + SQRT(1.D0 - (delta**2)))/delta
		theta_rcrn = ATAN(rrcrn)
		theta_rcrp = ATAN(rrcrp)
!990111-000316 Make sure the right choice of arctan(R_r)=[\theta_r] is made.
!	Arctan covers the range (-pi/2,pi/2) while 
!	\theta_r_Crit must be in the range (-pi/4,3pi/4) (The range of Ri>0.)
		IF(theta_rcrn.LT.-pi/4.D0) theta_rcrn_deg = (theta_rcrn + pi)*(180.D0/pi)
		IF(theta_rcrp.LT.-pi/4.D0) theta_rcrp_deg = (theta_rcrp + pi)*(180.D0/pi)
!	theta_rcrn_deg = theta_rcrn*(180.D0/pi)
!	theta_rcrp_deg = theta_rcrp*(180.D0/pi)
		if (mytid.eq.0) then
			WRITE(10,*) "   "
			WRITE(10,*) "   "
			WRITE(10,*) "   "
			WRITE(10,*) "   "
			WRITE(10,*) "R_r_Crit+ =",rrcrp
			WRITE(10,*) "R_r_Crit- =",rrcrn
			WRITE(10,*) "\\theta_r_Crit+ =",theta_rcrp
			WRITE(10,*) "\\theta_r_Crit- =",theta_rcrn
			WRITE(10,*) "\\theta_r_Crit+ in degrees =",theta_rcrp_deg
			WRITE(10,*) "\\theta_r_Crit- in degrees =",theta_rcrn_deg
			WRITE(10,*) "   "
			WRITE(10,*) "   "
!*****C
			WRITE(10,*) " "
			WRITE(10,*) " "
		endif
!	Increments in radial and angular coordinates in (Ri_T,Ri_C) plane.
		delra_r = 1.D0/DFLOAT(mt_ra_r)
		deltheta_r = (pi/4.D0)/DFLOAT(n_theta_r_oct)
!	Calculate the ratio \sigma_sa_max \equiv S_S/S_H as a function 
!	of the angle \theta_r in Ri_T,Ri_C space,
!       \theta_r \equiv arctan(Ri_C/Ri_T). 
!981218 The range of angles where unrealizability occurs is 
!	a subset of theta_r = -pi/4 to 3pi/4.
		if (mytid.eq.0) then
			WRITE(10,*) "S_S/S_H at pre-maximum Ri as a function of" &
             // "\\theta_r \\equiv Arctan(Ri_C/Ri_T)" 
!981220    Absurd default on sisamax \equiv S_S/S_H.
			WRITE(10,*) "Arbitrarily show the absurd value -99.999"
			WRITE(10,*) "at angles where do not have "// &
			"a maximum Ri (or radius ra_r)."
			WRITE(10,*) " "
			WRITE(10,*) "  \\th_r ^o  ra_r      " &
			// "  Ri_T        Ri_C        Ri         Ri_d       " &
			// "  S_M       S_H       S_S      S_S/S_H  "
		endif
!*
!       For Ri_T and Ri_C positive find the realizability limits  
!	in polar coordinates in the (Ri_T,Ri_C) plane : (ra_r,theta_r).
		IF(ifpolartablewrite.EQ.1) then
			if (mytid.eq.0) then
				OPEN(UNIT=68,FILE="turb_ra_th",STATUS="NEW")
			endif
		endif
!		write(10,*)"FINISHED AT 1186"
		DO itheta_r = -n_theta_r_oct,3*n_theta_r_oct
			theta_r = DFLOAT(itheta_r)*deltheta_r
			theta_r_deg = theta_r*(180.D0/pi)
!980119	Introduce jtheta_r, an angle index that begins at zero   
!	for the purposes of letting OURSAL2 know it starts at the origin.
			jtheta_r = itheta_r + n_theta_r_oct
!981220-1 Initialize sisamax to the impossible negative value of -99.999 to 
!	let places where the realizability limit is not reached stand out.
			sisamax(itheta_r) = -99.999
!981221  Initialize sm_r0,sh_r0,ss_r0 to the INCONSISTENT absurd value -9.999999.
			sm_r0 = -9.999999
			sh_r0 = -9.999999
			ss_r0 = -9.999999
!990303 Flag ibg determines if the background value of ra_r has been calculated.
!030116X1a Option had never worked. George Halliwell points out I erroneously had "isalback" here before.
			IF(ifsalback.EQ.6) ibg=0
!*****C
!990315 Flag ifunreal determines if realizability limit has been found.
			ifunreal=0
!000315	 Write to file "back_ra_r_values", the allowed values of background ra_r.
			IF(itheta_r.EQ.-n_theta_r_oct) then
				if (mytid.eq.0) then
!lhl0711	   OPEN(UNIT=69,FILE="back_ra_r_values",STATUS="NEW")
				!	OPEN(UNIT=69,FILE="back_ra_r_values.txt")
				 	OPEN(UNIT=69,FILE=bid_fi_r_string)
				endif
			endif
!*****C
!000318 Make the ra_r max value not too large to try to avoid numerical trouble.
!            write(10,*)"FINISHED AT 1216"
	  loop2:DO ira_r = 0,(mt_ra_r**2)/4
				IF(ira_r.LE.mt_ra_r) THEN
					ra_r = DFLOAT(ira_r)*delra_r
				ELSE
					ra_r = ((1.D0+delra_r)**(ira_r - mt_ra_r)) &
							*(DFLOAT(mt_ra_r)*delra_r)
				END IF
		
!981216-990119	Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
				rit = ra_r*COS(theta_r)
				ric = ra_r*SIN(theta_r)
				ri_r  = rit + ric
				rid_r = rit - ric
!981216 Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
				IF(nmodel.EQ.1) THEN
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
					CALL OURSAL2_2A(b1,ri_r,rid_r,slq2_r, &
							sm_r,sh_r,ss_r, &
							c_y0,c_y00,ira_r,jtheta_r)
					
!			         write(10,*)"FINISHED OURSAL2_2A AT 1237"
				END IF
				
				if (mytid.eq.0) then
!					IF(itheta_r.EQ.-n_theta_r_oct) WRITE(69,9200) ira_r,ra_r
					IF(itheta_r.LE.3*n_theta_r_oct) then
						WRITE(69,*) "itheta_r =",itheta_r," ira= ",ira_r," ra_r= ",ra_r
						WRITE(69,*) "sm =",sm_r," sh= ",sh_r," sc= ",ss_r," slq2=",slq2_r
					else
						close(69)
					end if
				endif
				
				
				
				
				if (mytid.eq.0) then
					IF(ifpolartablewrite.EQ.1) WRITE(68,9001) &
					itheta_r,theta_r_deg,ira_r,ra_r,slq2_r,sm_r,sh_r,ss_r
				endif
!990303 Calculate S_M/(S l/q) and find where it's backfact of its origin value.
				IF(ifsalback.EQ.6) THEN
					smosloq_r = sm_r/SQRT(slq2_r)
					IF(ira_r.EQ.0) smosloq_0r = smosloq_r
!	Use radius where dimensionless K_M falls below backfact*origin value.
					IF((smosloq_r.LE.backfact*smosloq_0r).AND.(ibg.EQ.0)) THEN
						ra_r1            	= ra_r
						rit1    	   	= rit
						ric1  	 	= ric
						ri_r1   		= ri_r
						rid_r1  		= rid_r
						slq2_r1(itheta_r)	= slq2_r
						sm_r1(itheta_r)  	= sm_r
						sh_r1(itheta_r)	= sh_r
						ss_r1(itheta_r)	= ss_r
						ibg=1
					END IF
!*****C
				END IF
				IF(slq2_r.LE.0.D0) THEN 
!981216	Use value of last lattice point on this radius with "slq2" positive.
!	Calculate the ratio of the salt and heat diffusivities there.
					sisamax(itheta_r) = ss_r0/sh_r0 
!990301 Store in an array the maximum radius, ra_r, at this angle, theta_r,
!	in the polar (Ri_T,Ri_C) [that is the (theta_r,ra_r)] plane.
					ra_rmax(itheta_r) = ra_r0
!	Determine the background radius, ra_r, at this \theta_r.
					IF(ifsalback.EQ.5) THEN
!       Use a constant fraction of the maximum radius before model breakdown.
						back_ra_r(itheta_r) = backfrac*ra_rmax(itheta_r)
!990303-15
					ELSE IF(ifsalback.EQ.6) THEN
						back_ra_r(itheta_r) = ra_r1
					END IF
					ifunreal = 1 
!*****C             
         !           write(10,*)"FINISHED AT 1302"
!981230 Skip straight to write out when last point reached.
					EXIT LOOP2
!					GO TO 16
				END IF
				ra_r0   = ra_r
				rit0    = rit
				ric0    = ric
				ri_r0   = ri_r
				rid_r0  = rid_r
				slq2_r0 = slq2_r
				sm_r0   = sm_r
				sh_r0   = sh_r
				ss_r0   = ss_r
!000319 Store c_y as c_y_0 for possible use as a  guess in background calc. .
				c_y_r0(itheta_r) = c_y0
				
			! IF(ira_r.EQ.(mt_ra_r**2)/4) then
                ! write(10,*)"FINISHED AT 1319"
			! end if	
				
			
			END DO LOOP2
!000315 Close file with background ra_r values.

			if (mytid.eq.0) then
				IF(itheta_r.EQ.-n_theta_r_oct) CLOSE(69)
			endif
!981216-30 Write out stability functions, the S's and sisamax.
!   16  continue
			if (mytid.eq.0) then
				WRITE(10,9150) theta_r_deg,ra_r0,rit0,ric0,ri_r0,rid_r0, &
					sm_r0,sh_r0,ss_r0,sisamax(itheta_r)
			endif
!990315	 Set background ra_r large at angles where unrealizability doesn't occur.
!000318 Make the ra_r max value not too large to try to avoid numerical trouble.
			IF(ifunreal.EQ.0) THEN
				ipenra_r = (mt_ra_r**2)/4-1
				back_ra_r(itheta_r) = ((1.D0+delra_r)**(ipenra_r - mt_ra_r)) &
                         *(DFLOAT(mt_ra_r)*delra_r) 
			END IF
!*****C
!990315 For ifsalback=5 case get value for initialization of c_y calculation. 
			IF(ifsalback.EQ.5) THEN
				IF(jtheta_r.EQ.0) THEN
					c_y001 = c_y0
				END IF
			END IF
!*****C

            ! IF(itheta_r.EQ.3*n_theta_r_oct) then
                ! write(10,*)"FINISHED AT 1354"
			! end if

		END DO
		
!		write(10,*)"FINISHED AT 1358"
		
		if (mytid.eq.0) then
			IF(ifpolartablewrite.EQ.1) CLOSE(68)
		endif
!990303-16 Write out stability functions at background ra_r .
		IF(ifsalback.GT.4) THEN
			DO itheta_r = -n_theta_r_oct,3*n_theta_r_oct
				theta_r = DFLOAT(itheta_r)*deltheta_r
				theta_r_deg = theta_r*(180.D0/pi)
!981216-990119	Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
				rit1 = back_ra_r(itheta_r)*COS(theta_r)
				ric1 = back_ra_r(itheta_r)*SIN(theta_r)
				ri_r1  = rit1 + ric1
				rid_r1 = rit1 - ric1
!990315-16 Calculation of turbulence functions for ifsalback=5 case.
				IF(ifsalback.EQ.5) THEN
!981216 Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
					jtheta_r = itheta_r + n_theta_r_oct
!990315 Set second table index to 1 to use last step's value except at start.
!000319 Transform that "last step" value from the most recent angle step to the
!	final realizable ra_r step at {\it this} angle in hope of more accuracy.
					IF(nmodel.EQ.1) THEN
!021210X1 Call version of submodule oursal2 which has an option to use B1={\tau S}^{3/2} (0,0).
						CALL OURSAL2_2A(b1,ri_r1,rid_r1,slq2_r1(itheta_r), &
							sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
							c_y_r0(itheta_r),c_y001,jtheta_r,1)
	    
					END IF
				END IF
				if (mytid.eq.0) then
					IF(itheta_r.EQ.-n_theta_r_oct) THEN
						WRITE(10,*) " "
						WRITE(10,*) &
					"Values at background ra_r=(Ri_T^2 + Ri_C^2)^(1/2)"
						WRITE(10,*) "\\th_r ^o   ra_r       " &
					// "Ri_T       Ri_C       Ri         Ri_d       " &
					// "(Sl/q)^2   S_M       S_H       S_S       S_S/S_H  "
						WRITE(10,*) " "
					END IF
				endif
				sisa1 = ss_r1(itheta_r)/sh_r1(itheta_r)
				if (mytid.eq.0) then
					WRITE(10,9160) theta_r_deg,back_ra_r(itheta_r), &
							rit1,ric1,ri_r1,rid_r1,slq2_r1(itheta_r), &
							sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
							sisa1
				endif
				IF(slq2_r1(itheta_r).LT.0.D0) THEN
					if (mytid.eq.1) then
						WRITE(10,*) &
						"Negative (Sl/q)^2 in table of Background vs. \\theta_r."
						WRITE(10,*) "itheta_r=",itheta_r, &
						"   slq2_r1(itheta_r)=",slq2_r1(itheta_r)
						WRITE(10,*) "Program is stopping in turb_2."
					endif
					STOP
				END IF
			END DO
		END IF
!*****C
!*
 !       write(10,*)"FINISHED AT 1383"
		if (mytid.eq.0) then
			! WRITE(10,*) " "
			! WRITE(10,*) "************************************************"
			! WRITE(10,*) " "
!030404Y Write out to standard output and a file a table of dimensionless turbulence 
!Y	functions of N_d^2/N^2 for the zero shear unstable case.
!030424Z Only write in case zeroshear model enabled.
			IF(ifzeroshear) THEN 
				! CLOSE(68)
				! WRITE(10,*) "************************************************"
				! WRITE(10,*)          "index (N_d^2)/(N^2)  -(\\tau N)^2    "// &
					! "S_M            S_H            S_S            "
				! WRITE(10,*) " "
!lhl0711	OPEN(UNIT=68,FILE="turb_nd2on2",STATUS="NEW")
!				OPEN(UNIT=68,FILE="turb_nd2on2")
				OPEN(UNIT=68,FILE=bid_fi_r_string)
					WRITE(68,*) "************************************************"
					WRITE(68,*)          "index (N_d^2)/(N^2)  -(\\tau N)^2    "// &
					"S_M            S_H            S_S            "
					WRITE(68,*) " "
				DO ind2on2 = -mt,mt
	!				WRITE(10,9268) ind2on2,and2on2a1(ind2on2), &
	!				amtaun2a1(ind2on2),sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2)
					WRITE(68,9268) ind2on2,and2on2a1(ind2on2), &
					amtaun2a1(ind2on2),sma1(ind2on2),sha1(ind2on2),ssa1(ind2on2)
				END DO
	!			WRITE(10,*) " "
				CLOSE(68)
	!			WRITE(10,*) "************************************************"
!******Y
			END IF

		endif
		ifirst=1
!981015-990702 Calculate constants inside subroutine oursal2 
!              or mikesal2 in sali-temp model case.
	else 
!*****C END OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
		if(nmodel.eq.1) then
!980912 Choose which set of model constants to use for Cheng model.
			IF(ifchengcon.EQ.0) THEN
				call ccoeff(b1,rimax)
			END IF
     
		endif
!991107 Calculate constant lengthscale for the background for ifback > 2
!  	\Delta_0 = {B_1 pi \over (3 Ko)^{3/2}} l_0
!	l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
!	"back_l_0" is the constant background l_0 in centimeters.
		back_l_0 = (((3.D0*ako)**(3.D0/2.D0))/(b1*pi))*back_del_0
!*****C
!991107-08C Temperature=Salinity Diffusivity Models Writeouts for background
		IF(ifback.GE.4) THEN
			if (mytid.eq.0) then
				WRITE(10,*) " "
				WRITE(10,*) "************************************"
				WRITE(10,*) "Dubovikov Internal wave constants for background."
			endif
			IF(ifback.EQ.4) THEN	
				if (mytid.eq.0) then
				WRITE(10,*) "Internal wave Richardson number=",ri_internal
				endif
			ELSE IF(ifback.EQ.5) THEN
				if (mytid.eq.0) then
					WRITE(10,*) "Ratio of Background to Critical Ri = ",backfrac
				endif
			END IF
			if (mytid.eq.0) then
				WRITE(10,*) "************************************"
				WRITE(10,*) " "
			endif
		END IF
!******C
!020916X Switch for use of surface fluxes to dimensionalize turbulence model.
		if (mytid.eq.0) then
			WRITE(10,*) "************************************"
			WRITE(10,*) "isurfuse=",isurfuse
			WRITE(10,*) "************************************"
		endif
!******X
!       building the look-up tables of slq2,sm and sh vs. ri
		dri=(rimax-ri0)/DFLOAT(ntbl-1)
		do k=1,ntbl
			ria(k)=ri0+DFLOAT(k-1)*dri
			if(k.eq.ntbl) ria(k)=rimax
			if(nmodel.eq.1) then
				call ourl2(b1,ria(k),slq2a(k),sma(k),sha(k))
			endif
		end do
!		write(10,*)"FINISHED AT 1474"
!       end of building look-up tables
		if (mytid.eq.0) then
			WRITE(10,*) "nmodel=",nmodel," the table's ntbl=",ntbl
!991107
			WRITE(10,*) "Temperature=Salinity diffusivity model"
			WRITE(10,*) "rimax  =",rimax
			WRITE(10,*) "ifback =",ifback
!000215
!030721Z1a Write switch for Deardorff treatment.
			WRITE(10,*) " "
			WRITE(10,*) "icondear=",icondear
			IF(icondear.EQ.-1) THEN
				WRITE(10,*) "Do *not* use Deardorff lengthscale modification."
			ELSE IF(icondear.EQ.0) THEN
				WRITE(10,*) "Ye Cheng's old Deardorff:"// &
		"modify l and \tau N leaving S_X unmodified."
			ELSE IF(icondear.EQ.1) THEN
				WRITE(10,*) "Ye Cheng's new Deardorff:"// &
		"modify l but leave *both* \tau N and S_X unmodified."
			END IF
			WRITE(10,*) " "
!*****CZ1a
			WRITE(10,*) "ifepson2=",ifepson2
			IF(ifepson2.EQ.2) WRITE(10,*) & 
		"epsilon/(N^2) used even for strong mixing beneath weak mixing" 
!030429-0502Z1 Write switch for latitude dependence of background mixing.
			WRITE(10,*) "ifdeeplat=",ifdeeplat
			IF(ifdeeplat.GT.0) THEN
				WRITE(10,*) &
			"Use latitude dependence of interior ocean value of \\epsilon/N^2"
!Zi1a' Floor on latitude dependent factor in background mixing
				WRITE(10,*) "eplatidependmin=",eplatidependmin
			END IF
!020404D,030429Z1 Must write out the parameter epson2__ for reference \epsilon/N^2 
!Z1	since epson2_ and epson2 are no longer constants.
			IF(ifepson2.GT.0) WRITE(10,*) "epson2__=",epson2__
			WRITE(10,*) " "
				
		end if
		ifirst=1
	end if
	 
end if  !  ifirst = 1
!******
    ! write(10,*) "FINISHED AT 1513"
    ! write(10,*)"dri =",dri
!020219D REFERENCE NOTE: END OF INITIALIZATION.
!020219-21[Mummy's 81st Birthday]D Bottom level for use with bottom enhancement.
!	                               (Must be set outside initialization if block.)
    kbot=n+1
!*****CD
!020912-17X ** Surface Buoyancy Flux (can be used for dimensionalization of turbulence model) ** 
!X	**Total Buoyancy Flux = Sum of "turbulent" and "solar" contributions**
	buoytot = buoytur + buoysol
!******X
!980501-0716 Choose the definition of MLD. 
!	     If ifextermld=1, keep the one that was defined externally.
    ! IF(ifextermld.EQ.0) THEN
! !980527,050208 Use Mixed-Layer Routine only when there are at least two levels of sea.
		! IF(n.GT.0) THEN
			! IF(idefmld.EQ.0) THEN
				! call formld(z,t,amld,n)
	
			! END IF
        ! ELSE IF(n.EQ.0) THEN
			! amld = z(1)
        ! ELSE 
			! amld = 0.D0
        ! END IF
    ! END IF
!
    al0=0.17*amld
!980717  Write internal turbulence quantities to fort.91 when writing enabled.
!	 Headers for each outputstep.
!981104 Add S_M,H,S to outputs.
!020924X Add (epsilon (\tau S}^2) to outputs in the case where it is calculated.
    if (mytid.eq.0) then
        IF(ifoutput.EQ.1) THEN
			WRITE(91,*) " "
			IF(isurfuse.EQ.0) THEN
				WRITE(91,*) "z          al         slq2       "// &
					"ri1        rid1       "// &
					"sm         sh         ss         "// &
					"v_back     t_back     s_back     "
			ELSE IF(isurfuse.EQ.1) THEN
				WRITE(91,*) "z          al         slq2       "// &
					"ri1        rid1       "// &
					"sm         sh         ss         "// &
					"epsy       "// &
					"v_back     t_back     s_back     "
			END IF
		END IF
	endif
!******
!  START OF FIRST LOOP THROUGH LEVELS
    IF(ifepson2.EQ.2) THEN
!000302 Initialize switch for sub(background-only) depth. 
		ifbelow=0      
    END IF
!020219D REFERENCE NOTE: START OF PRIMARY LOOP OVER DEPTH LEVELS.
 !   write(10,*) "FINISHED AT 1569"

    do  k=1,n
! !030504-0803Zi1a Section for Gregg et al. parameterization case.
		IF((ifepson2.EQ.2).AND.(ifdeeplat.GT.0)) THEN
!030502Z1 Initialize switch for reversion to deep lengthscale for (N/f)<1 for foreground
!Z1	  in ifepson2=2,ifdeeplat=1 case because Gregg et al. formula is then unusable.
			ifnofsmall=0
!030504-0803Zi1a Calculate N from N^2 and set ifnofsmall flag when (N/f)<1
!Zi1a		 Note that when Gregg et al. use "f", it should really be interpreted as '|f|'.
!Zi1a	         See NBp.030803-3,5 .
			IF(an2(k).GE.0.D0) an = SQRT(an2(k))
			IF((an/ABS(Coriol)).LT.1.D0) THEN	!arccosh(N/f) is undefined can't use Gregg et al.
				ifnofsmall=1
			END IF
		END IF
!*****CZi1a
        ri1=ri(k)
		IF(ifsali.EQ.0) THEN
			if(ri1.ge.rimax) then
				ri(k)=rimax    ! ad hoc
!991107 Traditional received background case for Temperature model.
				IF(ifback.EQ.0) THEN
					akm(k)=v_back(k)
					akh(k)=t_back(k)
!981102 Set background salt=heat diffusivity.
					aks(k)=akh(k)
					exit
!991109C Set turbulence functions to values at rimax for ifback>0.
				ELSE
					sm   = sma(ntbl)
					sh   = sha(ntbl)
					slq2 = slq2a(ntbl)	
					ss = sh
!******C
				END IF
!*****C
			elseif(ri1.le.ri0) then !asymptotically
				sm=sma(1)
				sh=sha(1)
!981102 Set background S_Salt = S_Heat
				ss=sh
				slq2=slq2a(1)*ria(1)/ri1
			else ! linearly interpolate the look-up tables
				m=int((ri1-ri0)/dri)+1
				tmp=(ri1-ria(m))/(ria(m+1)-ria(m))
				sm=sma(m)+(sma(m+1)-sma(m))*tmp
				sh=sha(m)+(sha(m+1)-sha(m))*tmp
				slq2=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
				ss=sh
			endif
		ELSE IF(ifsali.EQ.1) THEN
!981104 Use Ri_d = Ri_C - Ri_T in salinity-temperature turbulence model.
				rid1=rid(k)
!*****C

!030424Z ONLY ENABLE ZERO SHEAR APPROXIMATION WHEN "ifzeroshear" TRUE.[See NBp.030424-12.]
!030401Y Code for Ri => - \infinity approximation taken from my HYCOM turb_2_fixed2.fs0
!Y	 [See NBp.030401-1to4.]For HYCOM N^2 was reconstructed now it's an input argument.
!
!030324-28 For unstable case with shear too small for 2D table of (Ri,Ri_d) to handle,
!       treat as if Ri were minus infinity, the unstable zero shear case, and
!       interpolate 1D table of (N_d^2/N^2): N^2 = N_H^2 + N_S^2 , N_d^2 = N_H^2 - N_S^2 .
!       Estimate N^2 and N_d^2 from the Ri and Ri_d . Note that the minimum, "epsi",
!       on shear in the calling routine can cause Ri and Ri_d to be underestimated.
!       (N^d)^2/N^2 is a more accurate quantity.
				and2 = (rid(k)/ri(k))*an2(k)
				and2on2 = and2/an2(k)
!       Consider shear as being too small when N^2 < Ri_table_minimum * S^2 .
!       Note N^2/S^2, which can go to infinity, is the real Ri, as opposed to the
!       sometimes smaller N^2/(MAX(S^2,epsil)) used in the model which is always finite.
				IF((ifzeroshear).AND.(an2(k).LT.rib(-mt)*s2(k))) THEN
					ifpureshear=1
				ELSE
					ifpureshear=0
				END IF
!******Z
				IF(ifpureshear.EQ.1) THEN
!030326 Introduce a modular 1D table interpolation derived
!       by stripping down the 2D table interpolation.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
!       imax is the maximum positive value of the table index
!       (equal to +mt when there is no realizability limit).
					imax = mt
					IF(ifastexpabs.EQ.1) THEN

						CALL INTERP1D_EXPABS(and2on2, &
									and2on2a1,amtaun2a1,sma1,sha1,ssa1, &
									amtaun2,sm,sh,ss, &
									imax,mt,mt0,dand2on2,rnd2on2)
					END IF
  !                
 9928  format('interp1dtest',3i5,1p,7e12.4)

!yXI

!       Reconstructed slq2 for the outputs - somewhat bogus.
					slq2 = (-amtaun2)/((b1**2)*ri(k))
!					GO TO 5               !Skip 2D interpolation.
					
				ELSE
					IF(ifastexpabs.EQ.1) THEN
				
					
					
!030424 Only calculate rnd2on2 when zero shear parameterization is enabled.
				
						CALL INTERP2D_EXPABS(ri1,rid1, &
							rib,ridb,slq2b,smb,shb,ssb, &
							slq2,sm,sh,ss, &
							irimax,mt,mt0,dri,rri)
					END IF
					
					
				END IF
!******
!******Y

!981015   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
				
!yXI
		END IF
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
		IF(slq2.LT.0.D0) THEN
			if (mytid.eq.1) then
				WRITE(10,*) "************************************************"
				WRITE(10,*) "Error detected in turbulence module." 
				WRITE(10,*) "'slq2' negative in turb_2 subroutine" &
					//" after interpolation."
				WRITE(10,*) "k=",k,"     slq2=",slq2
				WRITE(10,*) "sm=",sm,"   sh=",sh,"   ss=",ss
				WRITE(10,*) "ri1=",ri1,"    rid1=",rid1
				WRITE(10,*) "dri=",dri
				WRITE(10,*) "Program will stop."
				WRITE(10,*) "************************************************"
			endif
			STOP
		END IF
!*****C
!000302 Assume region contiguous with surface where foreground model is
!	realizable has ended when get 0 "slq2".
		IF(slq2.EQ.0.D0) ifbelow = 1
!*****C
!030401 Skipped from 1D table interpolation to here for unstable zero shear approximation.
!    5   CONTINUE

!020912-24,25X ******Calculate epsy(k) \equiv epsilon * (\tau S)^2 for surface forcing dimensionalization.******
!X	       Take epsilon = -buoytot/(1 - ((1/2)(\tau S)^2)S_M))
!X	       epsy = -buoytot/((b1^2 (Sl/q)^2)^{-1} - (1/2)S_M) [See NBp020912-8,12 
!X	       AND NBp020925-2{extension on cover}.(error of (1/4) for (1/2) had to be corrected.) ]
!020919X Keep track of total number of points to compare with number of problem points.
!020923X Define problem points as points where epsy<0 *AND* l^2 S dimensionalization is used 
!X       in isurfuse=0 case AND l is based on MLD. 
!020923X Define negative problem points: problem points with Ri<0 ,
!X	 Keep track of number of points where l^2 S dimensionalization is used AND l is based on MLD.
!020924,26X,030504Z1 Introduce lifupper, logical which is 1 only where l_MLD^2 S would be used for isurfuse=0 .
!26X-030803Zi1a      set lifupper false when unrealizable
		lifupper= (((ifepson2.LT.2).OR.(ifbelow.EQ.0)   &
              .OR.((ifdeeplat.GT.0).AND.(ifnofsmall.EQ.1)).OR.&	!030504 Revert to l^2 S for N/f<1 in Gregg et. al. case.  
            ((ri1.LT.0.D0).AND.(k.LE.2)))   &
            .AND.(slq2.GT.0.D0)) 
		IF(lifupper) THEN
			ilmldpoint=ilmldpoint+1
			IF(ri1.LT.0.D0) ilmldpointneg=ilmldpointneg+1
		END IF
		ipoint=ipoint+1
        IF(k.EQ.1) icall=icall+1
!020924,26X Initialize epsy and lifepsy for this point.
		epsy(k)=0.D0
		lifepsy(k)=.FALSE.
		IF((isurfuse.EQ.1).AND.(n.GT.0)) THEN		
!020924X Although it seems to be being calculated to zero anyway in unrealizable Ri cases,
!X        set epsy to zero  when slq2 is zero for safety's sake.
			IF(slq2.EQ.0.D0) THEN
				epsy(k)=0.D0
			ELSE
				epsy(k) = -buoytot/((1.D0/((b1**2)*(slq2))) - 0.5D0*sm)
			END IF
!020924,26X Introduce lifepsy(k), which is 1 only when epsy dimensionalization used.
			lifepsy(k)= ((epsy(k).GE.0.D0).AND.lifupper)
!020924X Comment out write outs of epsy diagnosis.
!C020918-26X Write negative (epsilon y)'s to fort.67  and all values to fort.68
			IF((epsy(k).LT.0.D0).AND.lifupper) THEN
				iproblem=iproblem+1
				IF(ri1.LT.0.D0) inegproblem=inegproblem+1

			END IF
! 
		END IF
! !******X
        akz=0.4*z(k)
        al=akz*al0/(al0+akz)
!030425Z **MODIFICATION OF LENGTHSCALE BY ROTATION OPTION**
!Z	 Use l=(l_Blackadar^{-1} + l_\Omega^{-1})^{-1} with l_\Omega \equiv \sqrt{-B*/f^3}
!Z	 *when* B*<0 and MLD extends to a set minimum [See NBp.030424-5to8&13&14,25-2&3.].
		IF(ilomega.EQ.1) THEN
			IF((buoytot.LT.0.D0).AND.(amld.GE.amldminlom)) THEN
				rlomega = SQRT((Coriol**3)/(-buoytot))
			ELSE
				rlomega = 0.D0
			END IF
			rlblackadar = 1.D0/al  
			rl = rlblackadar + rlomega
			al = 1.D0/rl
		END IF
!******Z
        al2=al*al
! !000302 Do not use Deardorff limitation when use (\epsilon/N^2) dimensionalization.
! !020924,26X FOR CONSISTENCY OF OUTPUTS WITH SURFACE BUOYANCY FORCING EPSILON DIMENSIONALIZATION, 
! !X       ALTHOUGH "slq2" IS NOT USED, DO NOT APPLY DEARDORFF LIMITATION ON "(Sl/q)^2" .
		IF(.NOT.(((ifepson2.EQ.2).AND.(ifbelow.EQ.1)).OR.lifepsy(k))) THEN
!         length scale reduction by buoyancy
!030717Z1a,050208 ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
			if(ri1.gt.0.D0) then
				anlq2=slq2*ri1
				if((anlq2.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
					al2=0.281D0/anlq2*al2
					IF(icondear.EQ.0) slq2=0.281D0/(ri1+1.D-20)
				endif
			endif
!*****CZ1a
        END IF
! !*****C

		IF(an2(k).LT.0.D0) THEN
			epson2_ = epson2__	 	!Just to "hold the place". Value irrelevent here.
		ELSE 
			IF(ifdeeplat.EQ.0) THEN
				epson2_ = epson2__
			ELSE IF((ifdeeplat.EQ.1).OR.(ifdeeplat.EQ.2)) THEN !030803Zi1a (See NBp.030803-4,5.)
!0504Z1  Since (N/f)<1 => arccosh(N/f) undefined, can't use Gregg et al. formula.
						IF(ifnofsmall.EQ.1) THEN	
							eplatidepend = 0.D0
						ELSE
							eplatidepend = EPLATIDEPEND_(ABS(Coriol),an)
						END IF
						eplatidepend = MAX(eplatidepend,eplatidependmin)
						epson2_ = epson2__*eplatidepend
			END IF
		END IF
! !*****CZ1
		IF(ifepson2.GE.1) THEN
!020214D Option for enhanced mixing near the bottom. 
!	  epson2_ is the "pelagic value" of epsilon/N^2 , that in the interior of the ocean.
!	  eps_bot is the enhanced bottom TKE dissipation, epsilon.
!	  epson2_bot = eps_bot / N^2 
			IF(ifbotenhance.EQ.0) THEN 
				epson2 = epson2_ 
			ELSE IF(ifbotenhance.EQ.1) THEN
				eps_bot = eps_bot0 * EXP((z(k) - z(kbot))/scale_bot)
				epson2_bot = eps_bot/(ri(k)*s2(k))
				epson2 = MAX(epson2_,epson2_bot)
			END IF
        END IF
! !*****CD
! !*****CD
! !************************************************************************
! !************************************************************************
! !991107 Change background diffusivities to value using ocean model N and
! !	background internal wave model Ri for ifback >= 4 .
! !************************************************************************
! !	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
! !	diffusivities calculated using the turbulence model
! !	with Ri and l_0 replaced by constants 'ri_internal' and 'back_l_0' 
! !	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
! !	to represent a modified Dubovikov internal wave generated turbulence
! !	with universal constant Richardson number for ifback=4 case.
		IF((ifback.GE.4).AND.(ifsali.EQ.0)) THEN
!990205 Use a constant background Ri estimate. 
			IF(ifback.EQ.4) THEN
				back_ri1  = ri_internal
			ELSE
!************************************************************************
!      	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with l_0 replaced by a constant 'back_l_0' and 
!	Ri by a constant 
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with constant Richardson number for ifback>4 case.
!
!  Calculate the background Ri .
				back_ri1 = backfrac*rimax
!*****C	    
			END IF
!       Use the calculated background Ri in the turbulence model.
!         linearly interpolate the look-up tables
			m=int((back_ri1-ri0)/dri)+1
			tmp=(back_ri1-ria(m))/(ria(m+1)-ria(m))
			sm_back=sma(m)+(sma(m+1)-sma(m))*tmp
			sh_back=sha(m)+(sha(m+1)-sha(m))*tmp
			slq2_back=slq2a(m)+(slq2a(m+1)-slq2a(m))*tmp
			ss_back=sh_back
			if (mytid.eq.0) then
!				WRITE(10,*) 'i=',ii,'j=',jj,'m',m
				WRITE(10,*)
				WRITE(10,*) 'tmp',tmp
				WRITE(10,*)
				WRITE(10,*) 'sm_back', sm_back,'sma(m)',sma(m),'sma(m+1)',sma(m+1)
				WRITE(10,*)
				WRITE(10,*) 'sh_back', sh_back,'sha(m)',sha(m),'sha(m+1)',sha(m+1)
				WRITE(10,*)
				WRITE(10,*) 'slq2_back', slq2_back,'slq2a(m)',slq2a(m),'slq2a(m+1)',slq2a(m+1)
				WRITE(10,*)
			endif
! Check that "slq2" has been set to 0 where it might have been negative.
			IF(slq2_back.LT.0) THEN
				if (mytid.eq.1) then
					WRITE(10,*) & 
							"************************************************"
					WRITE(10,*) "Error detected in turbulence module."
					WRITE(10,*) "'slq2_back' negative in turb_2 subroutine" &
								//" after 1D interpolation of background."
					WRITE(10,*) "k=",k,"     slq2_back=",slq2_back
					WRITE(10,*) &
								"sm_back=",sm_back,"   sh_back=",sh_back
					WRITE(10,*) "back_ri1=",back_ri1
					WRITE(10,*) "dri=",dri
					WRITE(10,*) "Program will stop."
					WRITE(10,*) & 
								"************************************************"
				endif
				STOP
			END IF
! !*****C
! !990205 Calculate the square of the shear from the background Richardson number.
! !	s2_back   = N^2 / ri_internal = (N^2 / S_ext^2) (S_ext^2 /ri_internal) 
! !	          = (Ri_ext / ri_internal) S_ext^2
			s2_back = (ri1/back_ri1)*s2(k)
!990208-0301 Set square of shear to zero for unstable density stratification.
			IF(ri1.LE.0.D0) s2_back = 0.D0
!*****C
!990316 Set ill-defined S_M,H,S for unstable density stratification to zero.
			IF(ri1.LT.0.D0) THEN
				sm_back = 0.D0
				sh_back = 0.D0
				ss_back = 0.D0
			END IF
! !*****C
! !000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
			IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
				al0_back = back_l_0
				akz=0.4D0*z(k)
				al_back=akz*al0_back/(al0_back+akz)
				al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
				if(back_ri1.gt.0.D0) then
					anlq2_back=slq2_back*back_ri1
					if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
						al2_back=0.281D0/anlq2_back*al2_back
						IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
					endif
				endif
!*****CZ1a
!990203-05,050208 Calculate the background diffusivities.
				tmp_back=0.5D0*b1*al2_back*sqrt(s2_back/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
			ELSE IF(ifepson2.GE.1) THEN
!040422Zi1bj-AH:MIT Dirty act as stopgap response to warning message produced by
!	  MIT ocean model. **Zero** al_back which has not been defined.
!	  [See NBp.040422-7&8.] Think should really calculated an equivalent
!	  lengthscale for the (l^2 S) background dimensionalization,
!	  but don't have time to work this out now. 
!	  al_back is *not* used to calculate the diffusivity in the 
!	  (\epsilon/N^2) dimensionalization case, but could be a diagnostic.
!	  Previously it had been printing out zero anyway because undefined.
					al_back=0.D0
					tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
			END IF
			v_back(k)=tmp_back*sm_back
			t_back(k)=tmp_back*sh_back
			s_back(k)=tmp_back*ss_back
!************************************************************************
		END IF
! !************************************************************************
! !************************************************************************
! !************************************************************************
! !991108 BEGIN SECTION FOR SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
		IF(ifsali.EQ.1) THEN
!*****C
!************************************************************************
!************************************************************************
!981125 Change background salt diffusivity from input value to 
!	(S_S/S_H)*(background heat diffusivity) for ifsalback=1 case.
            IF((ifsalback.EQ.1) .and. (slq2 .EQ. 0.D0) .and. (k .eq. 1))then
			    s_back(k) = t_back(k)
            
			ELSE IF(ifsalback.EQ.1) THEN
					IF(slq2.NE.0.D0) THEN
						sm_last = sm
						sh_last = sh
						ss_last = ss
					! ELSE IF(k.eq.1) THEN
! !       Set salt background = heat background if turbulenceless 1st layer.
					! s_back(k) = t_back(k)
					! GO TO 20 
				! END IF
						s_back(k) = (ss_last/sh_last)*t_back(k)
					end if
!************************************************************
!981216 Change background salt diffusivity to
!	(S_S/S_H)*(background heat diffusivity), 
!	but with S_S,H taken at Ri a little less than Ri_max
!	at the given Ri_T/Ri_C value for ifsalback=2 case.
!990412 When Ri_T = 0, 
!       correctly set the angle 'theta_r' in the (Ri_T,Ri_C) plane to 'pi'/2 , 
!	unless Ri_C = 0, in which case **arbitrarily** set 'theta_r' to 'pi'/4 .
			ELSE IF(ifsalback.EQ.2) THEN
					IF(slq2.NE.0.D0) THEN
						sisa = ss/sh
					ELSE 
!  	Linearly interpolate sisamax array to the angle in (Ri_C,Ri_T) space .
!	Ri \equiv Ri_T + Ri_C 	; Ri_d \equiv Ri_T - Ri_C .
						rit = (ri(k) + rid(k))/2.D0
						ric = (ri(k) - rid(k))/2.D0
!990412-13 Find \theta_r for the Ri_T = 0 case. Treat "0/0 = 1".
						IF(rit.EQ.0.D0) THEN
							IF(ric.EQ.0.D0) THEN
								theta_r = ATAN(1.D0)
							ELSE
								theta_r = pi/2.D0	! Arctangent of infinity.
							END IF
						ELSE
							theta_r = ATAN(ric/rit)
						END IF
!*****C
!990111 Make sure the right choice of arctan(Ri_C/Ri_T) [\theta_r] is made.
!	Arctan only covers the range (-pi/2,pi/2) which theta_r may be outside.
!000323 Choose to have theta_r in range (-pi/4,3pi/4), stable case only.
						IF(ABS(theta_r).GT.(pi/2.D0)) STOP
						IF(theta_r.LT.(-pi)/4.D0) theta_r = theta_r + pi
!000309 MAKE 'jtheta' A NON-NEGATIVE INDEX - ZERO AT THETA = -PI/4 .
!	The fortran function "INT" rounds to the integer *NEAREST TO ZERO*
!	**I.E. ROUNDS **UP** FOR NEGATIVE NUMBERS**, DOWN ONLY FOR POSITIVES.
						jtheta_r0 = INT((theta_r + (pi/4.D0))/deltheta_r)
!000309 INTRODUCE 'itheta' HERE FOR THE INDEX THAT IS ZERO AT THETA=0.
						itheta_r0 = jtheta_r0 - n_theta_r_oct
						itheta_r1 = itheta_r0+1
						theta_r0 = itheta_r0*deltheta_r
						theta_r1 = itheta_r1*deltheta_r
!000314 Angle in degrees.
						theta_r_deg = theta_r*180.D0/pi
!*****C
!	Sound the alarm if have unrealizability outside expected range in angle.
						IF((itheta_r1.GT.3*n_theta_r_oct).OR.  &
							(itheta_r0.LT.-n_theta_r_oct)) THEN
							if (mytid.eq.1) then
								WRITE(10,*) &
											"************************************************"
								WRITE(10,*) "Problem in turbulence module!"
								WRITE(10,*) "Unrealizability outside Ri>0 region. "
								WRITE(10,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
								WRITE(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
								WRITE(10,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
								WRITE(10,*) "theta_r_deg=",theta_r_deg
								WRITE(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
								WRITE(10,*) "n_theta_r_oct=",n_theta_r_oct 
								WRITE(10,*) " "
								WRITE(10,*) "Program will stop."
								WRITE(10,*) &
											"************************************************"
							endif
							STOP
						END IF
						deltheta_r1 = theta_r - theta_r0
						delsisa_r = sisamax(itheta_r1) - sisamax(itheta_r0)
						dsisa_o_dtheta = delsisa_r/deltheta_r
						sisa = sisamax(itheta_r0)+deltheta_r1*dsisa_o_dtheta
					END IF
					s_back(k) = sisa*t_back(k)
!************************************************************
!************************************************************************
!990202-03	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to a
!	small constant plus diffusivities calculated using the turbulence model
!	with l_0 and S^2 replaced by constants intended to represent 
!	an internal-wave-generated turbulence for the ifsalback=3 case.
			ELSE IF(ifsalback.EQ.3) THEN
!990202 Use a constant background shear estimate to calculate background Ri,Ri_d
					back_ri1  = ri1*s2(k)*back_sm2
					back_rid1 = rid1*s2(k)*back_sm2
!990203 Use the calculated background Ri and Ri_d in the turbulence model.
!981015-990203   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
					IF(ifastexpabs.EQ.1) THEN
							CALL INTERP2D_EXPABS(back_ri1,back_rid1, &
							rib,ridb,slq2b,smb,shb,ssb, &
							slq2_back,sm_back,sh_back,ss_back, &
							irimax,mt,mt0,dri,rri)
					END IF
!yXI
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
					IF(slq2_back.LT.0) THEN
						if (mytid.eq.1) then
							WRITE(10,*) "************************************************"
							WRITE(10,*) "Error detected in turbulence module."
							WRITE(10,*) "'slq2_back' negative in turb_2 subroutine" &
							//" after interpolation of background."
							WRITE(10,*) "k=",k,"     slq2_back=",slq2_back
							WRITE(10,*) & 
							"sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
							WRITE(10,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
							WRITE(10,*) "dri=",dri
							WRITE(10,*) "Program will stop."
							WRITE(10,*) "************************************************"
						endif
						STOP
					END IF
!*****C
!000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
					IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
						al0_back = back_l_0
						akz=0.4D0*z(k)
						al_back=akz*al0_back/(al0_back+akz)
						al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
						if(back_ri1.gt.0.D0) then
							anlq2_back=slq2_back*back_ri1
							if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
								al2_back=0.281D0/anlq2_back*al2_back
								IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
							endif
						endif
!*****CZ1a
!990203,050208 Calculate the background diffusivities.
						tmp_back=0.5D0*b1*al2_back*sqrt(back_s2/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
					ELSE IF(ifepson2.GT.0) THEN
						tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
					END IF
					v_back(k)=tmp_back*sm_back+v_back0
					t_back(k)=tmp_back*sh_back+t_back0
					s_back(k)=tmp_back*ss_back+s_back0
!************************************************************************
!************************************************************************
!990205-08	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with Ri and l_0 replaced by constants 'ri_internal' and 'back_l_0' 
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with constant Richardson number for ifsalback=4 case.
			ELSE IF(ifsalback.GE.4) THEN
!990205 Use a constant background Ri estimate. 
				IF(ifsalback.EQ.4) THEN
					back_ri1  = ri_internal
					back_rid1 = (rid1/ri1)*ri_internal
				ELSE
!************************************************************************
!990301	Change ALL THREE BACKGROUND DIFFUSIVITIES from input values to 
!	diffusivities calculated using the turbulence model
!	with l_0 replaced by a constant 'back_l_0' and 
!	Ri by a function of Ri_d
!	and S^2 replaced by  (N^2 / Ri_internal) for N^2>=0 and 0 for N^2 <0
!	to represent a modified Dubovikov internal wave generated turbulence
!	with stability-ratio dependent Richardson number for ifsalback>4 case.
!
!990412 When Ri_T = 0 and Ri_C \ne 0,
!       correctly set the angle 'theta_r' in the (Ri_T,Ri_C) plane to 'pi'/2 . 
!*****C
!990413	Skip background ra_r calculation in unstable *OR NEUTRAL* case.
!000323 Set background ra_r arbitrarily to zero in these cases.
					IF(ri(k).LE.0.D0) THEN
						back_ra_r1 = 0.D0
						back_rit1  = 0.D0
						back_ric1  = 0.D0
						back_ri1   = 0.D0
						back_rid1  = 0.D0
						
!						GO TO 19						
					ELSE

					
!*****C
!  	Linearly interpolate back_ra_r array to this angle in (Ri_C,Ri_T) space.
!	Ri \equiv Ri_T + Ri_C 	; Ri_d \equiv Ri_T - Ri_C .
						rit = (ri(k) + rid(k))/2.D0
						ric = (ri(k) - rid(k))/2.D0
						ra_r = SQRT((rit**2) + (ric**2))
!030403	Use same zero temp. gradient treatment as for isalback=2.[See NBp.030403-13.]
!990412-13 Find \theta_r for the Ri_T = 0 case. Treat "0/0 = 1".
						IF(rit.EQ.0.D0) THEN
							IF(ric.EQ.0.D0) THEN
								theta_r = ATAN(1.D0)
							ELSE
								theta_r = pi/2.D0       ! Arctangent of infinity.
							END IF
						ELSE
							theta_r = ATAN(ric/rit)
						END IF
!*****C
!						theta_r = ATAN(ric/rit)
!990111 Make sure the right choice of arctan(Ri_C/Ri_T) [\theta_r] is made.
!	Arctan only covers the range (-pi/2,pi/2) which theta_r may be outside.
!000323 Want to consider statically stable case only: Ri > 0.
						IF(ABS(theta_r).GT.(pi/2.D0)) STOP
						IF(theta_r.LT.(-pi)/4.D0) theta_r = theta_r + pi
!000309	 MAKE 'jtheta' A NON-NEGATIVE INDEX - ZERO AT THETA = -PI/4 .
!	The fortran function "INT" rounds to the integer *NEAREST TO ZERO*
!	**I.E. ROUNDS **UP** FOR NEGATIVE NUMBERS**, DOWN ONLY FOR POSITIVES.
						jtheta_r0 = INT((theta_r + (pi/4.D0))/deltheta_r)
						jtheta_r1 = jtheta_r0+1
!000309 INTRODUCE 'itheta' HERE FOR THE INDEX THAT IS ZERO AT THETA=0.
						itheta_r0 = jtheta_r0 - n_theta_r_oct
						itheta_r1 = itheta_r0+1
!000330   ***WHEN THE ANGLE IS BETWEEN THE ANGLE FOR REALIZABILITY AT INFINITY***
!	  ***AND THE LAST TABLE ANGLE BEFORE THAT CRITICAL ANGLE, *** 
!	  ***SET IT TO THE LAST TABLE ANGLE BEFORE THE CRITICAL ANGLE.****
						theta_r0 = itheta_r0*deltheta_r
						theta_r1 = itheta_r1*deltheta_r
						IF((theta_r0.LE.theta_rcrp).AND.(theta_r.GT.theta_rcrp)) THEN
							theta_r = theta_r1
							theta_r0 = theta_r1
							itheta_r0 = itheta_r1 
							itheta_r1 = itheta_r1+1
							theta_r1 = theta_r1 + deltheta_r
						ELSE IF((theta_r1.GE.theta_rcrn).AND. (theta_r.LT.theta_rcrn)) THEN
							theta_r = theta_r0
							theta_r1 = theta_r0
							itheta_r1 = itheta_r0 
							itheta_r0 = itheta_r0-1
							theta_r0 = theta_r0 - deltheta_r
						END IF
!*****C
!000314 Angle in degrees.
!modified at 2017-01-17 1:12 by Liu

						theta_r_deg = theta_r*180.D0/pi
!*****C
!	Sound the alarm if have unrealizability outside expected range in angle.
						IF((itheta_r1.GT.3*n_theta_r_oct).OR.(itheta_r0.LT.-n_theta_r_oct)) THEN
							if (mytid.eq.3) then
								WRITE(10,*) & 
									"************************************************"
								WRITE(10,*) "Problem in turbulence module!"
								WRITE(10,*) "Unrealizability outside Ri>0 region. "
								WRITE(10,*) "rit=",rit,"ric=",ric,"    theta_r=",theta_r
								WRITE(10,*) "theta_r_deg =",theta_r_deg," deltheta_r=",deltheta_r
								WRITE(10,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
								WRITE(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
								WRITE(10,*) "n_theta_r_oct=",n_theta_r_oct 
								WRITE(10,*) " "
								WRITE(10,*) "Program will stop."
								WRITE(10,*) & 
						"************************************************"
								write(10,*), "z(1:na+1): ", z(1:na+1)
								write(10,*), "t(1:na): ", t(1:na)
								write(10,*), "s(1:na): ", s(1:na)
								write(10,*), "rh(1:na): ", rh(1:na)
								write(10,*), "ri(1:na): ", ri(1:na)
								write(10,*), "rid(1:na): ", rid(1:na)
								write(10,*), "s2(1:na): ", s2(1:na)
								write(10,*), "an2(1:na): ", an2(1:na)
				
							endif
							STOP
						END IF
						deltheta_r1 = theta_r - theta_r0
						delback_ra_r = back_ra_r(itheta_r1) - back_ra_r(itheta_r0)
						dback_ra_r_o_dtheta = delback_ra_r/deltheta_r
						back_ra_r1 = back_ra_r(itheta_r0) + deltheta_r1*dback_ra_r_o_dtheta
!000316-17 In case choose ifrafgmax=1, ra_r is at maximum the ForeGround ra_r
!	at the "strong" double diffusive \theta_r's 
!       where have turbulence as Ri+> infinity. 
						ifrafglt=0
						IF(ifrafgmax.EQ.1) THEN
							IF((theta_r.LE.theta_rcrp).OR.(theta_r.GE.theta_rcrn)) THEN
								IF(back_ra_r1.GT.ra_r) THEN
									ifrafglt=1
									back_ra_r1=ra_r
								END IF
							END IF
						END IF
!   18   CONTINUE 
						IF(back_ra_r1.LT.0.D0) THEN
							if (mytid.eq.1) then
								WRITE(10,*) & 
											"************************************************"
								WRITE(10,*) "Problem in turbulence module!"
								WRITE(10,*) "Negative bg ra_r \\equiv (Ri_T^2+Ri_C^2)^(1/2)"
								WRITE(10,*) "back_ra_r1 =", back_ra_r1
								WRITE(10,*) "theta_r =", theta_r
								WRITE(10,*) " "
								WRITE(10,*) "slq2=",slq2,"    sm=",sm," sh=",sh," ss=",ss
								WRITE(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
								WRITE(10,*) "rit=",rit,"ric=",ric
								WRITE(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
								WRITE(10,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
								WRITE(10,*) "theta_r_deg =",theta_r_deg
								WRITE(10,*) "n_theta_r_oct=",n_theta_r_oct 
								WRITE(10,*) " "
								WRITE(10,*) " "
								WRITE(10,*) "Program will stop."
								WRITE(10,*) &
												"************************************************"
							endif
							STOP
						END IF 
!990301-0323 Calculate the background Ri and Ri_d .
						back_rit1 = COS(theta_r)*back_ra_r1
						back_ric1 = SIN(theta_r)*back_ra_r1
						back_ri1  = back_rit1 + back_ric1
						back_rid1 = back_rit1 - back_ric1
					END IF
!*****C	    
				END IF
!000315 CALCULATE THE BACKGROUND DIMENSIONLESS TURBULENCE FUNCTIONS
!	USING TABLE OF VALUES FOR BACKGROUND "\theta_r"'S 
!	FOR "ifbg_theta_interp"=1.
!000317 Can only use theta_r table when do *not* reduce ra_r_BackGround
!	to a smaller ra_r_ForeGround.
                IF(ri(k).GT.0.D0) THEN
					IF((ifbg_theta_interp.EQ.0).OR.(ifrafglt.EQ.1)) THEN
!990203 Use the calculated background Ri and Ri_d in the turbulence model.
!981015-990203   Interpolate 2D table for salinity-temperature model case.
!011107yXI ***Option to call instead of generic interpolation one tailored to exponential absolute nonlinear part.***
						IF(ifastexpabs.EQ.1) THEN
							CALL INTERP2D_EXPABS(back_ri1,back_rid1, &
									rib,ridb,slq2b,smb,shb,ssb, &
									slq2_back,sm_back,sh_back,ss_back, &
									irimax,mt,mt0,dri,rri)
						END IF
!yXI
					ELSE IF(ifbg_theta_interp.EQ.1) THEN
!000315 Interpolate 1D table of background vs. theta_r instead.
							deltheta_r1 = theta_r - itheta_r0*deltheta_r
							delsm_back = sm_r1(itheta_r1) - sm_r1(itheta_r0)
							dsm_back_o_dtheta = delsm_back/deltheta_r
							sm_back = sm_r1(itheta_r0) + & 
									  deltheta_r1*dsm_back_o_dtheta
							delsh_back = sh_r1(itheta_r1) - sh_r1(itheta_r0)
							dsh_back_o_dtheta = delsh_back/deltheta_r
							sh_back = sh_r1(itheta_r0) + & 								
							deltheta_r1*dsh_back_o_dtheta
							delss_back = ss_r1(itheta_r1) - ss_r1(itheta_r0)
							dss_back_o_dtheta = delss_back/deltheta_r 
							ss_back = ss_r1(itheta_r0) + & 
							deltheta_r1*dss_back_o_dtheta
							delslq2_back = slq2_r1(itheta_r1) - slq2_r1(itheta_r0)
							dslq2_back_o_dtheta = delslq2_back/deltheta_r
							slq2_back = slq2_r1(itheta_r0) + &
									   deltheta_r1*dslq2_back_o_dtheta
!    	if(ii.eq.14 .and. jj.eq.4) then
!        WRITE(10,*) 'i=',ii,'j=',jj
!        WRITE(10,*)
!        WRITE(10,*) 'itheta_r0', itheta_r0
!        WRITE(10,*)
!        WRITE(10,*) 'deltheta_r1', deltheta_r1
!        WRITE(10,*)
!        WRITE(10,*) 'dsm_back_o_dtheta,sm_r1(itheta_r0)',dsm_back_o_dtheta,sm_r1(itheta_r0)
!        WRITE(10,*)
!        WRITE(10,*) 'dsh_back_o_dtheta,sh_r1(itheta_r0)',dsh_back_o_dtheta,sh_r1(itheta_r0)
!        WRITE(10,*)
!        WRITE(10,*)'dss_back_o_dtheta,ss_r1(itheta_r0)', dss_back_o_dtheta,ss_r1(itheta_r0)
!        WRITE(10,*)
!        WRITE(10,*)'dslq2_back_o_dtheta,slq2_r1(itheta_r0)', dslq2_back_o_dtheta,slq2_r1(itheta_r0)
!        WRITE(10,*)
!	endif

					ELSE
						if (mytid.eq.1) then
							WRITE(10,*) "Problem with choice of background interpolation."
							WRITE(10,*) "ifbg_theta_interp=",ifbg_theta_interp
							WRITE(10,*) "ifrafglt=",ifrafglt
							WRITE(10,*) "Program is stopping."
						endif
						STOP
					END IF
!*****C
!981216-990108 Check that "slq2" has been set to 0 where it might have been negative.
					IF(slq2_back.LT.0) THEN
						if (mytid.eq.1) then
							WRITE(10,*) "************************************************"
							WRITE(10,*) "Error detected in turbulence module."
							WRITE(10,*) "'slq2_back' negative in turb_2 subroutine" &
										//" after interpolation of background."
							WRITE(10,*) "k=",k,"     slq2_back=",slq2_back
							WRITE(10,*) &
								  "sm_back=",sm_back,"   sh_back=",sh_back,"   ss_back=",ss_back
							WRITE(10,*) "back_ri1=",back_ri1,"   back_rid1=",back_rid1
							WRITE(10,*) "dri=",dri
							WRITE(10,*) "Program will stop."
							WRITE(10,*) "************************************************"
						endif
						STOP
					END IF
!*****C
!990205 Calculate the square of the shear from the background Richardson number.
!	s2_back   = N^2 / ri_internal = (N^2 / S_ext^2) (S_ext^2 /ri_internal) 
!	          = (Ri_ext / ri_internal) S_ext^2
					s2_back = (ri1/back_ri1)*s2(k)
				END IF
!990208-0301 Set square of shear to zero for unstable density stratification.
!   19    		IF(ri1.LE.0.D0) s2_back = 0.D0
				IF(ri1.LE.0.D0) s2_back = 0.D0
!*****C
!990316,050208 Set ill-defined S_M,H,S for unstable density stratification to zero.
				IF(ri1.LT.0.D0) THEN
					sm_back = 0.D0
					sh_back = 0.D0
					ss_back = 0.D0
				END IF
!*****C
				IF((sm_back.LT.0.D0).OR.(sh_back.LT.0.D0).OR. (ss_back.LT.0.D0)) THEN
					if (mytid.eq.1) then
						WRITE(10,*) & 
									"************************************************"
						WRITE(10,*) "Problem in turbulence module!"
						WRITE(10,*) "Negative Structure Function in Background."
!						WRITE(10,*) 'i=',ii,'j=',jj
						WRITE(10,*) 'temperature'
						WRITE(10,'(5d15.5)') t
						WRITE(10,*)
						WRITE(10,*) 'salinity'
						WRITE(10,'(5d15.5)') s
						WRITE(10,*)
						WRITE(10,*) 'density'
						WRITE(10,'(5d15.5)') rh
						WRITE(10,*)
						WRITE(10,*) 'richardson'
						WRITE(10,'(5d15.5)') ri
						WRITE(10,*)
						WRITE(10,*) 'rid'
						WRITE(10,'(5d15.5)') rid
						WRITE(10,*)
						WRITE(10,*) 'shear'
						WRITE(10,'(5d15.5)') s2
						WRITE(10,*)
						WRITE(10,*) 'BV frequency'
						WRITE(10,'(5d15.5)') an2			
						WRITE(10,*)
						WRITE(10,*) 'ustart=',ustar_
						WRITE(10,*)
						WRITE(10,*) 'buoytur=',buoytur
						WRITE(10,*)
						WRITE(10,*) 'buoysol=',buoysol
						WRITE(10,*) "slq2_back=",slq2_back
						WRITE(10,*) "sm_back=",sm_back, &
									   " sh_back=",sh_back, &
									   " ss_back=",ss_back
						WRITE(10,*) " "
						WRITE(10,*) "back_ra_r1 =", back_ra_r1
						WRITE(10,*) "theta_r =", theta_r
						WRITE(10,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
						WRITE(10,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
						WRITE(10,*) " "
						WRITE(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
						WRITE(10,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
						WRITE(10,*) "theta_r_deg=",theta_r_deg
						WRITE(10,*) "n_theta_r_oct=",n_theta_r_oct 
						WRITE(10,*) " "
						WRITE(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
						WRITE(10,*) "rit=",rit,"ric=",ric
						WRITE(10,*) " "
						WRITE(10,*) " "
						WRITE(10,*) "Program will stop."
						WRITE(10,*) & 
        "************************************************"
					endif
					STOP
				END IF
!000215 Skip background lengthscale calculation when using K_X/(epsilon/N^2) .
				IF(ifepson2.EQ.0) THEN
!990203,050208 Use the constant background l_0 lengthscale in the turbulence model.
					al0_back = back_l_0
					akz=0.4D0*z(k)
					al_back=akz*al0_back/(al0_back+akz)
					al2_back=al_back*al_back
!         length scale reduction by buoyancy
!030717Z1a ***REMOVE DEARDORFF LIMITATION ON LENGTHSCALE IN "icondear=-1" CASE,***
!	   ***APPLY DEARDORFF LIMITATION ONLY TO LENGTHSCALE *NOT* TO {\tau N} IN "icondear=+1" CASE.***
!	   ***TRADITIONAL METHOD OF LIMITING {\tau N} AFTER S_X SET IS INCONSISTENT, "icondear=0" CASE.***
!	   DEARDORFF LIMITED {\tau N} BUT THIS IS NOT REALLY CONSISTENT WITH PRODUCTION=DISSIPATION.
					if(back_ri1.gt.0.D0) then
						anlq2_back=slq2_back*back_ri1
						if((anlq2_back.gt.0.281D0).AND.(icondear.GE.0)) then  !0.281=0.53**2
							al2_back=0.281D0/anlq2_back*al2_back
							IF(icondear.EQ.0) slq2_back=0.281D0/(back_ri1+1.D-20)
						endif
					endif
!*****CZ1a
!990203-05 Calculate the background diffusivities.
					tmp_back=0.5D0*b1*al2_back*sqrt(s2_back/(slq2_back+1.D-40))
!000215 Use K_X = K_X/(\epsilon/N^2) * (\epsilon/N^2)    
!		From NBp.000215-5, Volume IX : 
!        	 K_X/(\epsilon/N^2) = (1/2) B_1 Ri (S l/q)^2 S_X  .
!	    K_X = (((1/2) B_1^2 Ri (S l/q)^2)* (\epsilon/N^2)) * S_X 
				ELSE IF(ifepson2.GT.0) THEN
							tmp_back=0.5D0*b1**2*back_ri1*slq2_back*epson2
				END IF
				v_back(k)=tmp_back*sm_back
				t_back(k)=tmp_back*sh_back
				s_back(k)=tmp_back*ss_back
!000309 Stop if background diffusivities are negative.
				IF((v_back(k).LT.0.D0).OR.(t_back(k).LT.0.D0).OR. (s_back(k).LT.0.D0)) THEN
					if (mytid.eq.1) then
						WRITE(10,*) &
								"************************************************"
						WRITE(10,*) "Problem in turbulence module!"
						WRITE(10,*) "Negative Background Diffusivity."
!						WRITE(10,*) 'i=',ii,'j=',jj
						WRITE(10,*) 'temperature'
						WRITE(10,'(5d15.5)') t
						WRITE(10,*)
						WRITE(10,*) 'salinity'
						WRITE(10,'(5d15.5)') s
						WRITE(10,*)
						WRITE(10,*) 'density'
						WRITE(10,'(5d15.5)') rh
						WRITE(10,*)
						WRITE(10,*) 'richardson'
						WRITE(10,'(5d15.5)') ri
						WRITE(10,*)
						WRITE(10,*) 'rid'
						WRITE(10,'(5d15.5)') rid
						WRITE(10,*)
						WRITE(10,*) 'shear'
						WRITE(10,'(5d15.5)') s2
						WRITE(10,*)
						WRITE(10,*) 'BV frequency'
						WRITE(10,'(5d15.5)') an2
						WRITE(10,*)
						WRITE(10,*) 'ustart=',ustar_
						WRITE(10,*)
						WRITE(10,*) 'buoytur=',buoytur
						WRITE(10,*)
						WRITE(10,*) 'buoysol=',buoysol
						WRITE(10,*) "v_back=",v_back, &
									   " t_back=",t_back, &
									   " s_back=",s_back
						WRITE(10,*) " "
						WRITE(10,*) "slq2_back=",slq2_back
						WRITE(10,*) "sm_back=",sm_back, &
									   " sh_back=",sh_back, &
									   " ss_back=",ss_back
						WRITE(10,*) " "
						WRITE(10,*) "back_ra_r1 =", back_ra_r1
						WRITE(10,*) "theta_r =", theta_r, &
									   "   theta_r_deg=",theta_r_deg
						WRITE(10,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
						WRITE(10,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
						WRITE(10,*) " "
						WRITE(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
						WRITE(10,*) "rit=",rit,"ric=",ric
						WRITE(10,*) " "
						WRITE(10,*) " "
						WRITE(10,*) "Program will stop."
						WRITE(10,*) &
							 "************************************************"
					endif
					STOP
				END IF
!000314 Stop if background diffusivities are zero at positive Ri.
!                   v_back(k)=0.01D0
!                   t_back(k)=0.01D0
!                   s_back(k)=0.01D0
				IF((ri(k).GT.0.D0).AND.((v_back(k).EQ.0.D0).OR. (t_back(k).EQ.0.D0).OR. (s_back(k).EQ.0.D0))) THEN
					if (mytid.eq.1) then
						WRITE(10,*) & 
								"************************************************"
						WRITE(10,*) "Problem in turbulence module!"
						WRITE(10,*) "Zero Background Diffusivity in stable case."
!						WRITE(10,*) 'i=',ii,'j=',jj
						WRITE(10,*) 'temperature'
						WRITE(10,'(5d15.5)') t
						WRITE(10,*)
						WRITE(10,*) 'salinity'
						WRITE(10,'(5d15.5)') s
						WRITE(10,*)
						WRITE(10,*) 'density'
						WRITE(10,'(5d15.5)') rh
						WRITE(10,*)
						WRITE(10,*) 'richardson'
						WRITE(10,'(5d15.5)') ri
						WRITE(10,*)
						WRITE(10,*) 'rid'
						WRITE(10,'(5d15.5)') rid
						WRITE(10,*)
						WRITE(10,*) 'shear'
						WRITE(10,'(5d15.5)') s2
						WRITE(10,*)
						WRITE(10,*) 'BV frequency'
						WRITE(10,'(5d15.5)') an2
						WRITE(10,*)
						WRITE(10,*) 'ustart=',ustar_
						WRITE(10,*)
						WRITE(10,*) 'buoytur=',buoytur
						WRITE(10,*)
						WRITE(10,*) 'buoysol=',buoysol
						WRITE(10,*) "v_back=",v_back(k), &
									   " t_back=",t_back(k), &
										   " s_back=",s_back(k)
						WRITE(10,*) " "
						WRITE(10,*) "slq2_back=",slq2_back
						WRITE(10,*) "sm_back=",sm_back, &
										   " sh_back=",sh_back, &
										   " ss_back=",ss_back
						WRITE(10,*) " "
						WRITE(10,*) "slq2_r1(itheta_r0)=",slq2_r1(itheta_r0), &
									   " slq2_r1(itheta_r1)=",slq2_r1(itheta_r1)
						WRITE(10,*) "sm_r1(itheta_r0)=",sm_r1(itheta_r0), &
										   " sm_r1(itheta_r1)=",sm_r1(itheta_r1)
						WRITE(10,*) "sh_r1(itheta_r0)=",sh_r1(itheta_r0), &
										   " sh_r1(itheta_r1)=",sh_r1(itheta_r1)
						WRITE(10,*) "ss_r1(itheta_r0)=",ss_r1(itheta_r0), &
										   " ss_r1(itheta_r1)=",ss_r1(itheta_r1)
						WRITE(10,*) " "
						WRITE(10,*) "back_ra_r1 =", back_ra_r1
						WRITE(10,*) "theta_r =", theta_r
						WRITE(10,*) "theta_r_deg =", theta_r_deg
						WRITE(10,*) "back_rit1=",back_rit1,"back_ric1=",back_ric1
						WRITE(10,*) "back_ri1=",back_ri1,"back_rid1=",back_rid1
						WRITE(10,*) "itheta_r0=",itheta_r0," itheta_r1=",itheta_r1
						WRITE(10,*) "jtheta_r0=",jtheta_r0," jtheta_r1=",jtheta_r1
						WRITE(10,*) "n_theta_r_oct=",n_theta_r_oct 
						WRITE(10,*) "deltheta_r=",deltheta_r
						WRITE(10,*) " "
						WRITE(10,*) " "
						WRITE(10,*) "k=",k,"  ri(k)=",ri(k),"  rid(k)=",rid(k)
						WRITE(10,*) "rit=",rit,"ric=",ric
						WRITE(10,*) " "
						WRITE(10,*) " "
						WRITE(10,*) "Program will stop."
						WRITE(10,*) &
									"************************************************"
					endif
					STOP
				END IF
!************************************************************************
			END IF
! 20     			CONTINUE
!*****C
		END IF
!991108 END SECTION FOR SALINITY MODEL BACKGROUND DIFFUSIVITY CALCULATION.
!************************************************************************
!************************************************************************
!************************************************************************
!980717  Write internal turbulence quantities to fort.91 when writing active.
!981104 Add S_M,H,S to outputs.
!020924X,030404Y Add epsy= dissipation times (\tau Shear)^2 to outputs.
!X       [See NBp020912 Vol.V for epsy formula.]
		IF(ifoutput.EQ.1) THEN
			if (mytid.eq.0) then
				IF(isurfuse.EQ.0) THEN
					WRITE(91,9000) z(k),al,slq2,ri1,rid1,sm,sh,ss, &
									v_back(k),t_back(k),s_back(k)
				ELSE IF(isurfuse.EQ.1) THEN
					WRITE(91,9000) z(k),al,slq2,ri1,rid1,sm,sh,ss, &
								epsy(k), &
								v_back(k),t_back(k),s_back(k)
				END IF
			endif
!990208-0322
!991107	Include background turbulence function writeouts for Temp. model bg.
			IF(ifback.EQ.4.OR.ifsalback.EQ.4) THEN
!*****C
				if (mytid.eq.0) then
					IF(k.EQ.1) THEN
						WRITE(94,*) " "
						WRITE(94,*) &
         "z[cm]      l_back[cm] Ri-table   Ri_d-table Ri_d_back  " &
       //"s2_back    slq2_back  sm_back    sh_back    ss_back    "
					END IF
					WRITE(94,9000) z(k),al_back,ri1,rid1,back_rid1, &
										s2_back,slq2_back,sm_back,sh_back,ss_back 
				endif
!991107	Include background turbulence function writeouts for Temp. model bg.
			ELSE IF((ifback.GT.4).OR.(ifsalback.GT.4)) THEN
!*****C
				if (mytid.eq.0) then
					IF(k.EQ.1) THEN
						WRITE(94,*) " "
!020221[Mummy's 81st Birthday]-22D Introduce additional writeouts to fort.94:
!	                        epsilon/N^2 when used for dimensionalization,
!	                        eps_bot when enhanced bottom mixing is enabled
!	                        and eps_bot(z_bottom) only for special check case.
!22		                Also N^2 and S^2 in all cases.
!22				Note that (N^2)/(S^2) is the true Richardson number,
!22				while ri1, like rid1 can be modified by INTERP2D
!22			        to fit within the table. ri(k) and rid(k) should be pristine.

						IF(ifepson2.EQ.0) THEN
							WRITE(94,*) &
											 "z[cm]      l_back[cm] Ri-table   Ri_d-table " &
										   //"Ri_back    Ri_d_back  ra_r_back  " &
										   //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
										   //"N^2        S^2        "
						ELSE
							IF(ifbotenhance.EQ.0) THEN
									WRITE(94,*) &
												 "z[cm]      l_back[cm] Ri-table   Ri_d-table " &
											   //"Ri_back    Ri_d_back  ra_r_back  " &
											   //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
											   //"N^2        S^2        " &
											   //"epson2     "
							ELSE
								WRITE(94,*) & 
											"z[cm]      l_back[cm] Ri-table   Ri_d-table " &
											   //"Ri_back    Ri_d_back  ra_r_back  " &
											   //"s2_back    slq2_back  sm_back    sh_back    ss_back    " &
											   //"N^2        S^2        " &
											   //"epson2     eps_bot    "
							END IF
						END IF
					END IF
				endif

				IF(ifepson2.EQ.0) THEN
					if (mytid.eq.0) then
						WRITE(94,9010) z(k),al_back,ri1,rid1, &
											back_ri1,back_rid1,back_ra_r1, &
											s2_back,slq2_back,sm_back,sh_back,ss_back, &
											ri(k)*s2(k),s2(k)
					endif
				ELSE
					IF(ifbotenhance.EQ.0) THEN
						if (mytid.eq.0) then
							WRITE(94,9020) z(k),al_back,ri1,rid1, &
											back_ri1,back_rid1,back_ra_r1, &
										  s2_back,slq2_back,sm_back,sh_back,ss_back, &
										  ri(k)*s2(k),s2(k), &
											epson2
						endif
					ELSE
						if (mytid.eq.0) then
							 WRITE(94,9030) z(k),al_back,ri1,rid1, &
											back_ri1,back_rid1,back_ra_r1, &
										  s2_back,slq2_back,sm_back,sh_back,ss_back, &
										  ri(k)*s2(k),s2(k), &
										epson2,eps_bot
						endif
!	Special check of bottom epsilon value.
						IF((ifcheckbottomeps.EQ.1).AND.(k.EQ.n)) THEN
								eps_bot__under = &
								eps_bot0 * EXP((z(k+1) - z(kbot))/scale_bot)
							if (mytid.eq.0) then
								WRITE(94,9040) z(k+1),eps_bot__under
							endif
						END IF
					END IF
				END IF
!*****CD	     
			END IF
!*****
		END IF
!******
!030401 Apply zeroing of l_deep taken from HYCOM turb_2_fixed2.fs0
!030328,050208 To avoid confusion set l_deep to zero when it is not calculated.
		aldeep(k)=0.D0

!020924X ******WHEN HAVE SET "isurfuse=1" DO *NOT* DIMENSIONALIZE WITH THE "MLD" ******
!X       ******BASED "Blackadar" LENGTHSCALE. IN CASES WHERE "isurfuse=0" USES******
!X	 ******"l_MLD \equiv MLD based Blackadar lengthscale" USE INSTEAD******
!X       ******'S SURFACE BUOYANCY BASED DIMENSIONALIZATION:******
!X       ******"K_X_foreground = ((epsilon*(\tau S)^2)/(2 S^2))S_X"(SEE NBp020912-12)******
!X       ******EXCEPT WHERE "epsilon*(\tau S)^2 \equiv epsy <0" IN WHICH CASE ******
!X       ******REVERT TO "al^2 S". LEAVE "l_deep^2 S" CASES ALONE.(SEE NBp020920-1,2)******
!000229-0302 In the case where the model is realizable at 
!       the Ri obtained from the external Shear,
!	*but* there is a level above where it is NOT thus realizable, 
!	USE THE "epsilon/(N^2)" DIMENSIONALIZATION FOR "ifepson=2".
!000311 EXCEPT IF  "Ri<0" DO *NOT* USE "epsilon/(N^2)" DIMENSIONALIZATION
!	BECAUSE IT PRODUCES NEGATIVE DIFFUSIVITIES IN THIS CASE.
!	*INSTEAD USE "l_deep^2 S", WHERE "l_deep" IS DERIVED FROM "rho" PROFILE.
!	"|{{d \rho / dz} \over {d2 \rho / dz^2}}| takes place of MLD in l_deep".
!000314 BUT *REVERT* TO "MLD" IN CASES OF FIRST TWO LEVELS !020924 Use epsy if isurfuse=1.
		IF((ifepson2.EQ.2).AND.(ifbelow.EQ.1)) THEN
!030504Z1 **ALSO USE "l_deep^2 S" IN Gregg et al. LATITUDE DEPENDENCE CASE WHEN (N/f)<1.**
			IF(ri1.GE.0.D0.OR.((ifdeeplat.EQ.2).AND.(ifnofsmall.EQ.1))) THEN
				tmp=0.5D0*b1**2*ri1*slq2*epson2
			ELSE IF(k.GT.2) THEN
				IF(k.EQ.n) THEN
					delz = z(k) - z(k-1)
					delrh = rh(k) - rh(k-1)
					del2rh = rh(k) - 2.D0*rh(k-1) + rh(k-2)
				ELSE
					delz = z(k+1) - z(k-1)
					delrh = rh(k+1) - rh(k-1)
					del2rh = rh(k+1) - 2.D0*rh(k) + rh(k-1)
				END IF
				dzrh = delrh/delz
				d2zrh = 4.D0*del2rh/(delz**2)
!000323 rdzlndzrh = *Reciprocal* of Dz_{ln(Dz_{rh})} = Dz_{rh}/Dz2_{rh} .
				rdzlndzrh = dzrh/d2zrh
				al0deep=0.17D0*ABS(rdzlndzrh)
				akz=0.4D0*z(k)
				aldeep(k)=akz*al0deep/(al0deep+akz)
				al2=aldeep(k)*aldeep(k)
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
!040217Zi1b When ifshearmin=.TRUE. introduce a minimum foreground shear to avoid
!	    the singularity that occurs in the (lengthscale^2 Shear) dimensionalization
!	    of the foreground turbulence model at zero shear, or else switch to
!	    (lengthscale^2 Brunt Vaisala frequency) dimensionalization. This deals with the
!	    problem that shear=0  implies (l^2 S)/(slq2+1.D-40) =0, which in turn makes the
!	    diffusivities zero, which is unphysical in the unstable case.
				IF(ifpureshear.EQ.1) THEN 
				   tmp=0.5D0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
!					GO TO 21
				ELSE IF(ifshearmin) THEN
					s2(k) = MAX(s2(k),s2min)
				END IF
				tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
!X In isurfuse=1 case, when buoyancy forcing gives positive dissipation epsilon, 
!X and hence epsy, use epsilon instead of l^2 S as basis for dimensionalization.
			ELSE
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
				IF(ifpureshear.EQ.1) THEN 
!					GO TO 21
					tmp=0.5D0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
					
				ELSE IF(ifshearmin) THEN
					s2(k) = MAX(s2(k),s2min)
				END IF
				IF(lifepsy(k)) THEN
					tmp=0.5D0*epsy(k)/(s2(k)+1.D-40)
				ELSE 
					tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
				END IF
			END IF
		ELSE
!030401Y For case where use pure convection model taken from my HYCOM turb_2_fixed2.
			IF(ifpureshear.EQ.1) THEN 
			    tmp=0.5D0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
!				GO TO 21
			ELSE IF(ifshearmin) THEN
				s2(k) = MAX(s2(k),s2min)
			END IF
!******Zi1b
			IF(lifepsy(k)) THEN
				tmp=0.5D0*epsy(k)/(s2(k)+1.D-40)
			ELSE 
				tmp=0.5D0*b1*al2*sqrt(s2(k)/(slq2+1.D-40))
			END IF
		END IF
!******X
!030401Y,050208 Calculate \epsilon \tau for the pure convection model from HYCOM turb_2_fixed2.
!030324-28 For very negative Ri use zero shear unstable case approximation.
!       K_X = e \tau S_X = (l^2 \sqrt(-N^2)) (B_1^2/2) (-(\tau N)^2)^{-1/2} S_X
!  21   IF(ifpureshear.EQ.1) tmp=0.5D0*(b1**2)*al2*sqrt(-an2(k)/amtaun2)
!******Y
		akm(k)=min(tmp*sm+v_back(k),visc_cbu_limit)
		akh(k)=min(tmp*sh+t_back(k),diff_cbt_limit)
		aks(k)=min(tmp*ss+s_back(k),diff_cbt_limit)
 END DO
 !   write(10,*) "FINISHED AT 2818"
!020219 REFERENCE NOTE: END OF PRIMARY LOOP OVER DEPTH LEVELS.
!980527 Zero diffusivity arrays at points involving rock or land.

    do k=nb+1,nmax 
		akm(k)=0.D0
		akh(k)=0.D0
		aks(k)=0.D0
    end do 
!980527 Introduce minimum windmixing only where there are at least two ocean levs
	IF(n.GT.0) THEN
		if(akm(1).lt.wndmix) akm(1)=wndmix
		if(akh(1).lt.wndmix) akh(1)=wndmix
		if(aks(1).lt.wndmix) aks(1)=wndmix
	END IF
				
!modified by Liu at 2017-07-17 09:34
!980716-19 Write turbulence inputs and outputs to fort.92 when writing active.
!	 Headers for each outputstep.
!981102    Add Ri_d \equiv Ri_T - Ri_C to the outputs.
!000311 Add foreground lengthscale to the outputs. 
!030404Y Add new turb input the square of the Brunt Vaisala frequency to the outputs.
!030425Z Add the Coriolis parameter and rotational lengthscale l_Omega to the outputs.
!25Z	 Also add the surface buoyancy forcing and friction velocity to the outputs.
!25Z	 Calculate the complex diagnostic l_Omega = \sqrt{-B* f_coriolis^{-3}} .
!     write(10,*) "FINISHED AT 2842"
	zlomega = SQRT(CMPLX(-buoytot/(Coriol**3)))
	IF(ifoutput.EQ.1) THEN
		if (mytid.eq.0) then
			WRITE(92,*) " "
			WRITE(92,*) "f_coriolis=",Coriol
			WRITE(92,*) "MLD[cm] = ",amld
			WRITE(92,*) "buoytot[cm^2/s^3] =",buoytot, &
								 "    ustar_[cm/s] =",ustar_
			WRITE(92,*) "l_\\Omega[cm] =",zlomega
			WRITE(92,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "// &
								 "Ri         Ri_d	   S^2[/s2]   "// &
								 "N^2[/s2]   "// &
								 "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "// &
								 "l_deep[cm] "
		endif
		DO k =1,n
			if (mytid.eq.0) then
				WRITE(92,9000) z(k),t(k),s(k),rh(k), &
								ri(k),rid(k),s2(k), &
								an2(k), &
								akm(k),akh(k),aks(k), &
								aldeep(k)
			endif
!000309-12 STOP IF DIFFUSIVITY IS NEGATIVE.
			IF((akm(k).LT.0.D0).OR.(akh(k).LT.0.D0).OR.(akm(k).LT.0.D0)) THEN
				if (mytid.eq.1) then
					WRITE(10,*) "Diffusivity is negative."
					WRITE(10,*) "k=",k
					WRITE(10,*) "z[cm]      tem[C]     sal[ppt]   rho[g/cm3] "// &
									 "Ri         Ri_d	   S^2[/s2]   "// &
									 "K_M[cm2/s] K_H[cm2/s] K_S[cm2/s] "
					WRITE(10,9000) z(k),t(k),s(k),rh(k), &
										ri(k),rid(k),s2(k), &
										akm(k),akh(k),aks(k)
					WRITE(10,*) " "
					WRITE(10,*) "Program will stop."
				endif
				STOP
			END IF
		END DO
	END IF	
!******
     CLOSE(10)
 9000 FORMAT(12(1pe11.3))
!020221[Mummy's 81st Birthday]-22D
 9010 FORMAT(14(1pe11.3))
 9020 FORMAT(15(1pe11.3))
 9030 FORMAT(16(1pe11.3))
 9040 FORMAT(1pe11.3,154X,1pe11.3)
!*****D
 9001 FORMAT(2(I8,'   ',1pe11.3),8(1pe11.3))
 9050 FORMAT(I8,'  ',2E16.4,I8,'  ')
 9100 FORMAT(I8,'  ',2E12.4,3F11.6,2F11.4)
 9150 FORMAT(F11.3,5E12.4,3F10.6,F9.3)
 9160 FORMAT(F11.3,6E10.4,3F10.6,F9.3)
 9200 FORMAT(I12,'    ',5E16.6)
 9268 FORMAT(I6,5F15.9)
      return
      end subroutine nasagiss_canuto



	function eplatidepend_(f,an) !checked

	!use precision_
	!use precision_mod
	implicit none
!	
    real(kind=r8),intent(in)::f,an
    real(kind=r8)			::eplatidepend_
	
	real(kind=r8),parameter	::an0=5.24D-3
	
	real(kind=r8),save		::den
	integer,save			::icalled=0
	
	
    real(kind=r8)			::x,xf,yn
    real(kind=r8)			::ACOSH1,WAVELAT
    real(kind=r8)			::f_30,anum,pi,omega
	

    acosh1(x)       =   log(x + sqrt((x**2) - 1.d0))
    wavelat(xf,yn)  =   xf*acosh1(yn/xf) 	

    !!calculate f_30
    pi     = 4.D0*atan(1.D0)
    omega  = pi/43082.0D0
	f_30=omega

	if(icalled.EQ.0) then
	  den=wavelat(f_30,an0)
	  icalled=icalled+1
	end if

	anum = WAVELAT(f,an)
	eplatidepend_ = anum/den


	end function eplatidepend_
	
!-----------------------------------------------------------------------
!     finds mixed layer depth
!-----------------------------------------------------------------------
!980501 Choice of definitions.
    subroutine formld(z,t,amld,n) !checked 

    use param_mod, only: mytid
	implicit none

    integer,intent(in)::n
	real(kind=r8),intent(in)::t(n),z(n)
	real(kind=r8),intent(out)::amld
    integer::k
    real(kind=r8)::tm



    loop1:do k=1,n
			if (abs(t(k) - t(1)).gt.0.1) then

				tm = t(1) - sign(0.1D0,t(1) - t(k))
				amld = z(k) + (z(k-1) - z(k))* &
					(tm - t(k))/(t(k-1) - t(k) + 1.e-20)
				exit loop1
			else
				amld=z(n)

			endif
		end do loop1

      end subroutine formld

!0501 Temperature difference criterion used, but chose outside \Delta T = delte.

!0501 Density difference criterion used, but chose outside \Delta \rho = delrh.
     
!******


!-----------------------------------------------------------------------
!     beginning of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
    subroutine ccoeff(b1,rimax) !checked
!980501 Make double precision to conform to calling cctmix routine.

      use param_mod, only: mytid
	  use global_g
	  implicit none
	  real(kind=r8)::prt0,qus,beta5,qty,temp1,temp2,cthem1,aa,bb,cc,c7

	  real(kind=r8),intent(out)::rimax,b1

!----------------------------------------------------------------------
!     defines constants and model coefficients
!     b1 == (q/ustar)**3,  qus == q/ustar
!     g2=(1.-3.*v2dq2)/2.
!     g3=(3.*u2dq2-1.-g2)/3.
!     w2dq2=(1+g2-3.*g3)/3.
!     shih & shabbir: g2=0.00336, g3=0.0906
!     apbl: g2=0.06, g3=0.16
!     beta5=0.6*(1-c5), if c5=.3 then beta5=.42
!     qty == us**2*theta2/(w*theta)**2
!     cthem1 == 1/(c sub theta)
!-----------------------------------------------------------------------
      prt0=1.0
      b1=16.6d0
      qus=b1**(1./3.)
      g2 = 0.00336
      g3 = 0.0906
      g6 = 0.4
      g7 = 0.
      beta5=0.42
      qty=3.1
      g1=4./3.*(3.*g3**2-g2**2)+4.*qus**(-4)
      g4=15./8.*beta5*g1
      temp1=1./6.*prt0*(1.+g2-3.*g3)*qus**4
      temp2=9.*(g6-g7)*(g6+g7+2.*prt0)/ &
      ( (1.+g2-3.*g3)*prt0 )**2 *qus**(-4)
      g5=temp1*(1.+sqrt(1.+temp2))
      cthem1=1./(prt0*qus**2)* qty  ! 0.952766
      c7=1./3.
      g8 = (1.-c7)*cthem1
!
      s0=3./2.*g1*g5*g5
      s1=-(g6+g7)*g4 + 2.*(g1-g3-1./3.*g2)*g4*g5 + 3./2.*g1*g5*g8
      s2=-3./8.*(g6*g6-g7*g7)*g1
      s3=3./2.*(g6 + g7)*g4 + (3.*g3 + g2)*g4*g5
      s4=2.*g5
      s5=2.*g4
      s6=-2./3.*g5*(g2*g2-3.*g3*g3)+1./2.*g1*g5*(g2-3.*g3)+3./4.*g1*(g6-g7)
      s7=3.*g5
      s8=3.*g4
      s9=g5*(3.*g3*g3-g2*g2)
!
      d0=3.*g5*g5
      d1=g5*(7.*g4 + 3.*g8)
      d2=g5*g5*(3.*g3*g3 - g2*g2) - 3./4.*(g6*g6 - g7*g7)
      d3=g4*(4.*g4 + 3.*g8)
      d4=g4*(g2*g6 - 3.*g3*g7 - g5*(g2*g2 - g3*g3)) + g5*g8*(3.*g3*g3 - g2*g2)
      d5=1./4.*(g2*g2 - 3.*g3*g3)*(g6*g6 - g7*g7)
!     find rimax:
      aa=(d3+g4)
      bb=(d4-s1/2.+s6/2.)
      cc=d5-s2/2.
      rimax=(-bb+sqrt(bb**2-4.*aa*cc))/(2*aa)
      rimax=rimax*0.999  !ad hoc
      return
      end subroutine ccoeff

     

    subroutine ourl2(b1,ri,slq2,sm,sh)  !checked
!980501 Make double precision to conform to calling cctmix routine.
      use param_mod, only: mytid
	  use global_g
	  implicit none
	  
	  
	  real(kind=r8),intent(in)::b1,ri
	  real(kind=r8),intent(out)::slq2,sm,sh
	  real(kind=r8)::aa,bb,cc,temp,q,y,dd
	  

      aa = (d3+g4)*ri*ri+(d4-s1/2.+s6/2.)*ri+d5-s2/2.
      bb = (d1+g5)*ri+d2-s0/2.
      cc = d0
      temp = bb**2-4.*aa*cc
      if(temp.le.0.) then
         WRITE(10,*) 'in turb_2, temp <= 0, stop, ri =',ri
         stop
      endif

      q = -1./2.*(bb+sign(1.D0,bb)*sqrt(temp))



      if(bb.lt.0.) then
          y=cc/q
      else
          y=q/aa
      endif
      if(y.lt.0.) then
       if (mytid.eq.1) then
          WRITE(10,*) 'in turb_2, y < 0, stop; ri =', ri
       endif
          stop
      endif
      dd=d0+(d1*ri+d2)*y+(d3*ri*ri+d4*ri+d5)*y*y
      sm=(s0+(s1*ri+s2)*y)/dd
      sh=(s4+(s5*ri+s6)*y)/dd
      slq2=y/(b1*b1)
      return
      end subroutine ourl2

!---------------------------------------------------------------------
!     end of improved turbulence model subroutines (nmodel=1)
!-----------------------------------------------------------------------
    function fct(y)                 !checked             
!980501 Make double precision to conform to calling cctmix routine.

       implicit none
	   real(kind=r8),intent(in)::y
	   real(kind=r8)::ri,x,sm,sh
	   real(kind=r8)::fct
!      implicit real*8 (a-h,o-z)
	   common /mikebb/ri
	   x=ri*y
      call mksmsh(x,y,sm,sh)
      fct=2./(sm-ri*sh)
      return                                          
      end function fct        
                                    
    subroutine mksmsh(x,y,sm,sh)    !checked
!980717 Make double precision to conform to calling cctmix routine.

      use param_mod, only: mytid
	  implicit none
	  
	  real(kind=r8),intent(in)::x
	  real(kind=r8),parameter::x1min=15.d0/8.d0
	  real(kind=r8)::x1,fc,y,ri,a,eta,xi
	  real(kind=r8),intent(out)::sm,sh


      x1=0.5d0*sqrt(y)
      if(x1.lt.x1min) then
         fc=5.d0/9.d0
      else
         fc=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
      endif
      y=(2.d0*x1)**2
      ri=x/y
      a=1.d0+1.d0/(1.+0.16d0*x1**2*ri)
      eta=0.05d0*x1**2*ri*a
      sm=1.d0/25.d0*fc/(1.d0+5.d0/(9.d0*fc)*eta)
      xi=x/60.d0
      sh=0.056d0/(1.d0+2.4d0*xi)
      return
      end subroutine mksmsh
	  
	 real(kind=r8) function fct_sal(c_y)         !checked                       

      use param_mod, only: mytid
	  use global_t
	  
	  implicit none
	  real(kind=r8),intent(in)::c_y
	  real(kind=r8)::c_n,c_c,sm,sh,sc

!000210 Decide to use the parameter "tpvot" instead of its value 2/5 \tau .
      c_n = -((tcot*tctot)/(tpvot**2))*c_y*rit
      c_c = -((tcot**2)/(tpvot**2))*c_y*ric
      call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)
!980609 y(S_\nu - Ri_T S_h - Ri_C S_c) = 8/25 . 8/25 = 0.32 . S_\nu = sm.
!	y = 0.32/(S_\nu - Ri_T S_h - Ri_C S_c). 
      fct_sal=(2.D0*(tpvot**2))/(sm-rit*sh-ric*sc)
      
	  return
      end function fct_sal     
	  
	  

	SUBROUTINE oursal2_2a(b1_arg, &                      !checked
                        ri,rid,slq2,sm,sh,sc,c_y0,c_y00,iri,irid)
!****************************************



      use param_mod, only: mytid
	  use global_t
	  implicit none
	  
	  real(kind=r8),parameter::b1_0=16.6D0, c_yst0 = 8.527882D0
	  logical,parameter::ifchengb1=.FALSE.
	  
	  
	  
	  integer,intent(in)::iri,irid
	  real(kind=r8),intent(in)::ri,rid
	  real(kind=r8),intent(out)::sm,sh,sc,slq2,b1_arg
	  
	  
	  integer::ib1set0,iend,ib1set,ier
	  real(kind=r8)::b1, eps,c_yst,c_n,c_c,c_y,val,c_y0,c_y00
	  real(kind=r8)::tau2,tau3,taus2,taus3 

	  real(kind=r8),intrinsic::SQRT 


     ib1set=0
	 ib1set0=ib1set
	 eps=1.D-6                                              
     iend=300  
	 rit = (ri + rid)/2.D0
	 ric = (ri - rid)/2.D0
	 
	 
	 
!021210A Calculate B1 later if use Cheng's formula instead of an a priori constant. 
	IF(.NOT.ifchengb1) THEN
!030124B Only set B1 when it has not been set before. [See NBp.030124-18 Vol.XVII .] 
		IF(ib1set0.EQ.0) THEN
!030128B Set variable "b1" to parameter "b1_0". {See NBp.030128-4 .]
			b1 = b1_0
			b1_arg=b1
			ib1set=ib1set+1
		END IF
	ELSE
		IF(ib1set0.EQ.0) THEN
			taus2 = c_y/(tpvot**2)
			taus3 = taus2*SQRT(taus2)
			b1 = SQRT(taus3)
			b1_arg = b1
			ib1set=ib1set+1
		END IF
	END IF
	
	IF(b1_arg.LE.0.D0) THEN
		if (mytid.eq.1) then
		WRITE(10,*) "B1 <= ZERO."
		WRITE(10,*) "B1=",b1_arg
		WRITE(10,*) "Something must be wrong."
		WRITE(10,*) "Program is stopping in oursal2."
		endif
		STOP
	END IF
	
	IF(ib1set.NE.1) THEN
		if (mytid.eq.1) then
			WRITE(10,*) "Problem in oursal2; B1 not properly set."
			WRITE(10,*) "Number of times B1 set=",ib1set
			WRITE(10,*) "b1_arg=",b1_arg
			WRITE(10,*) "Program is stopping."
		endif
		STOP
	END IF
	
	
		
    call smshsc_a3(0.D0,0.D0,0.D0,sm,sh,sc)
!
                                                 
!     rimax= ?
!-----rtwi finds the root of x=fct(x)                     
!980615  Need a guess at the root, c_yst. Use neighboring solution.
!981007 Initial guess for c_yst for this value of Ri_d.
    IF(iri.EQ.0.AND.irid.EQ.0) THEN
		c_yst = c_yst0
    ELSE IF(iri.EQ.0) THEN
		c_yst = c_y00
    ELSE 
		c_yst = c_y0
    END IF
	
	
!981007 Calculate Ri_T =(Ri + Ri_d)/2 and Ri_C =(Ri - Ri_d)/2.  
	
    call rtwi(c_y,val,fct_sal,c_yst,eps,iend,ier) 
	

    if(ier.ne.0) then
!981022-23 Make error message more specific.
		if (mytid.eq.1) then
			WRITE(10,*) "In oursal2 subroutine"
			WRITE(10,*) "c_y00=",c_y00,"	c_y0=",c_y0
			WRITE(10,*) "ri=",ri,"	rid=",rid
			WRITE(10,*) "rit=",rit,"	ric=",ric
			WRITE(10,*) "Initial guess for rtwi c_yst=",c_yst
!*****C
            WRITE(10,*) "rtwi call problem, ier=",ier
       endif
	   stop
    endif
!981016 Store value of c_y for future guesses.
	IF(c_y.GE.0) THEN
	    c_y0=c_y
	ELSE 
!981022 Turbulence model becomes unphysical for c_y negative.
!981014-16  Realizability for negative Ri
		IF(ri.LT.0) THEN
			if (mytid.eq.1) then
			WRITE(10,*) "c_y negative at negative Ri"
			WRITE(10,*) "Ri=",ri," 	c_y=",c_y
			WRITE(10,*) "Unstable realizability limit unexpected:" 
			WRITE(10,*) "stopping in oursal2."
			endif
			STOP
        END IF
	END IF
	
	IF((iri.EQ.0).AND.(irid.EQ.0).AND.(ABS(c_y - c_yst0).GT.1.D-6)) THEN
		if (mytid.eq.1) then
			WRITE(10,*) "Inconsistency in neutral value of c_y"
			WRITE(10,*) "Value used =",c_yst0
			WRITE(10,*) "Value calculated =",c_y
			WRITE(10,*) "Program stopping in oursal2"
		endif
		STOP
	END IF
    IF(iri.EQ.0) c_y00=c_y

	

	IF(ib1set0.EQ.0) THEN
		if (mytid.eq.0) then
			WRITE(10,*) " "
			WRITE(10,*) "ifchengb1=",ifchengb1
			WRITE(10,*) "B_1=",b1_arg
			WRITE(10,*) " "
		endif
	END IF
	ib1set0=ib1set

	slq2 = c_y/((b1*tpvot)**2)

    c_n = -(tcot*tctot/(tpvot**2))*c_y*rit
    c_c = -((tcot**2)/(tpvot**2))*c_y*ric
    call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)


 1003 format(12(I8))
 1004 format(12(1pe14.5))
      end SUBROUTINE oursal2_2a




                                        
                                                    


    subroutine smshsc_a3(y,n,c,sm,sh,sc)  !checked
!030403Y Include "ttot'=`{\tau_\theta \over \tau}' in the common block with timescale 
!Y       ratios, /bb0/. See NBp030403-8to11.
!000125-27 NEW SUBROUTINE WHICH calculates the "p's" from the timescale ratios.
!	BASED on "smshsc2":
!990513 SUBROUTINE WHICH CALCULATES "p's" from "sgmt". BASED ON "smshsc1":
!990513 NEW SUBROUTINE WHICH USES YE CHENG'S FORTRAN CODE TO CALCULATE CONSTANTS
!	FROM THE "p's" SENT TO ME BY HIM TODAY. BASED ON "smshsc0". 
!980728 **CORRECT THE VALUE OF "p10".**
!	p_10 = {\tau_{p \theta} \tau_{c \theta}} \over {\tau_c ^ 2}
!******
!980609-15 Replace Cheng's smsh with  smshsc, which includes concentration.
!980610 The y,n,c used here are Canuto's "y,n,c" called c_y,c_n,c_c 
!	elsewhere in this program.

      use param_mod, only: mytid
	  use global_t
	  
	  implicit none
	  integer,parameter::ifmodelconstout=0
	  real(kind=r8),parameter::tpvot0=0.4D0,sgmt=0.72D0
	  real(kind=r8),parameter::tptot0=(1.D0/5.D0)*(1.D0/(1.D0+(1.D0/sgmt)))
	  real(kind=r8),parameter::tpcot0=tptot0,ttot0=sgmt,tcot0=ttot0,tctot0=1.D0/3.D0
	  
	  real(kind=r8),intent(in)::y,n,c
	  real(kind=r8),intent(out)::sm,sh,sc
	  real(kind=r8)::Nm,Nh,Nc,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p1m,p2m
	  real(kind=r8)::A0,A1,A2,A3,A4,A5
	  real(kind=r8)::D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D
	  integer::ifrecall
	  
	  
	  

!990513 Switch for whether(1) or not(0) to output p's a's and d's to a file.
 !     PARAMETER(ifmodelconstout=0)
!*****C
!000127	Add `\tau_pv \over \tau' to the common block with timescale ratios.
!      COMMON /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot
!980615	From the printed notes Canuto gave me on 980601 have:
!	 \tau_pv = {2 \over 5} \tau			(B.1)
 !        PARAMETER(tpvot0 = 0.4D0)
!****** "tpv/tau" = 2/5
!990513 Make "sgmt" a parameter. Standard value was 0.72.
! 	PARAMETER(sgmt=0.72D0)
!000125	Determine {\tau_p\theta,tau_pc,tau_\theta,tau_c,tau_c\theta }/tau
!	PARAMETER(tptot0=(1.D0/5.D0)*(1.D0/(1.D0+(1.D0/sgmt))))
!	PARAMETER(tpcot0=tptot0)
!	PARAMETER(ttot0=sgmt)
!	PARAMETER(tcot0=ttot0)
!	PARAMETER(tctot0=1.D0/3.D0)
!*****C

!****************************************
!000125 Set common block passable timescale ratios to parameter statement values
!000127 Include \tau_pv \over \tau .
	tpvot = tpvot0
	tptot = tptot0
	tpcot = tpcot0
	ttot  = ttot0
	tcot  = tcot0
	tctot = tctot0


 
!000125 Calculate the p's.
      p1  = 0.832D0
      p2  = 0.545D0
      p3  = (5.D0/2.D0)*tpcot
      p4  = (1.D0/5.D0)*tpcot*(tcot**(-2))
      p5  = tpcot*tctot*(tcot**(-2))
      p6  = (1.D0/5.D0)*(tcot**(-1))*(tctot**(-1))*tptot
      p7  = 5.D0*tctot
      p8  = (5.D0/2.D0)*tptot
      p9  = ttot*tptot*((tcot*tctot)**(-1))
      p10 = tctot*tptot*(tcot**(-2))
      p11 = tpcot*(tcot**(-1))
      p1m = 1.D0 - p1
      p2m = 1.D0 - p2

!##########################
      A0 = 12
      A1 = p11*(12*p9+8*p6-30*p6*p8-5*p6*(p1m+3*p2m))
      A2 = 5*(2*p4*p6*p7-p4*p9-p6*p11)*(p1m+3*p2m)+8*p6*p11+8*p4*p9-16*p4&
       *p6*p7+12*p11*p9+12*p11*p10-12*p4*p7**2*p6-30*p6*p11*p8+30*p4*p6*&
        p7*p8+30*p6*p4*p7*p3-30*p4*p9*p3
      A3 = p10*(12*p11+8*p4-30*p3*p4-5*p4*(p1m+3*p2m))
      A4 = -p6*(8-30*p8-5*p1m-15*p2m)-12*p9-12*p11
      A5 = -p4*(8-30*p3-5*p1m-15*p2m)-12*p10-12*p11
      D0 = 24
      D1 = p11*((-p6-2*p9)*p1m**2+(p6+6*p9)*p2m**2+2*p6*p8*(p1m-3*p2m))
      D2 = (2*p4*p6*p7-p4*p9-p6*p11)*(p1m**2-p2m**2)+2*(-p11*p10-p11*p9+&
        p4*p7**2*p6)*(p1m**2-3*p2m**2)+2*(-p6*p4*p7*p3-p4*p6*p7*p8+p4*p9*p3&
        +p6*p11*p8)*(p1m-3*p2m)
      D3 = p10*((-p4-2*p11)*p1m**2+2*p4*p3*(p1m-3*p2m)+(6*p11+p4)*p2m**2)
      D4 = -4*p6*p11*(3*p9+2*p6)
      D5 = 4*p4*p6**2*p7*(4+3*p7)-4*p4*p9*(3*p11+2*p6)-4*p6*p11*(3*p9+3*&
          p10+2*p4+2*p6)
      D6 = 4*p4**2*p6*p7*(4+3*p7)-4*p4*p9*(2*p4+3*p11)-8*p4*p6*(p11+p10)&
          -12*p10*p11*(p4+p6)
      D7 = -4*p4*p10*(2*p4+3*p11)
      D8 = (2*p9+2*p11+p6)*p1m**2-2*p6*p8*(p1m-3*p2m)-(p6+6*p9+6*p11)*p2m**2
      D9 = (2*p10+p4+2*p11)*p1m**2-2*p4*p3*(p1m-3*p2m)-(p4+6*p10+6*p11)*p2m**2
      D10 = 8*p6**2+4*(7*p11+3*p9)*p6+24*p11*p9
      D11 = -8*(4+3*p7)*p4*p6*p7+4*p4*(4*p6+7*p9+3*p11)+4*p6*(3*p10+7*p11)&
             +24*p11*(p10+p9)
      D12 = 4*p10*(7*p4+6*p11)+4*p4*(2*p4+3*p11)
      D13 = 6*p2m**2-2*p1m**2
      D14 = -28*p6-24*p9-24*p11
      D15 = -24*p10-28*p4-24*p11
!results.2_1
!****************************************

!980728	Write out the p's.
!000125 Writeout the timescale ratios as well.
!lhl
	ifrecall=1
!lhl
	IF(ifrecall.EQ.0) THEN
       if (mytid.eq.0) then
	  WRITE(10,*) "tau_pv/tau     =",tpvot 
	  WRITE(10,*) "tau_ptheta/tau =",tptot
	  WRITE(10,*) "tau_pc/tau =",tpcot
	  WRITE(10,*) "tau_theta/tau  =",ttot
	  WRITE(10,*) "tau_c/tau  =",tcot
	  WRITE(10,*) "tau_ctheta/tau  =",tctot
	  WRITE(10,*) " "
	  WRITE(10,*) "p1 =",p1
	  WRITE(10,*) "p2 =",p2
	  WRITE(10,*) "p3 =",p3
	  WRITE(10,*) "p4 =",p4
	  WRITE(10,*) "p5 =",p5
	  WRITE(10,*) "p6 =",p6
	  WRITE(10,*) "p7 =",p7
	  WRITE(10,*) "p8 =",p8
	  WRITE(10,*) "p9 =",p9
	  WRITE(10,*) "p10=",p10
	  WRITE(10,*) "p11=",p11
!990513 Write out the a's and d's as well.
	  WRITE(10,*) "a0=",a0
	  WRITE(10,*) "a1=",a1
	  WRITE(10,*) "a2=",a2
	  WRITE(10,*) "a3=",a3
	  WRITE(10,*) "a4=",a4
	  WRITE(10,*) "a5=",a5
	  WRITE(10,*) "d0=",d0
	  WRITE(10,*) "d1=",d1
	  WRITE(10,*) "d2=",d2
	  WRITE(10,*) "d3=",d3
	  WRITE(10,*) "d4=",d4
	  WRITE(10,*) "d5=",d5
	  WRITE(10,*) "d6=",d6
	  WRITE(10,*) "d7=",d7
	  WRITE(10,*) "d8=",d8
	  WRITE(10,*) "d9=",d9
	  WRITE(10,*) "d10=",d10
	  WRITE(10,*) "d11=",d11
	  WRITE(10,*) "d12=",d12
	  WRITE(10,*) "d13=",d13
	  WRITE(10,*) "d14=",d14
	  WRITE(10,*) "d15=",d15
!990513 Output p#, a# and d# to the file model_constants if the switch is set.
!000125 Writeout the timescale ratios as well.
          IF(ifmodelconstout.EQ.1) THEN
            OPEN(UNIT=66,FILE='model_constants',STATUS='UNKNOWN')
	    WRITE(10,*) "tau_pv/tau     =",tpvot 
	    WRITE(10,*) "tau_ptheta/tau =",tptot
	    WRITE(10,*) "tau_pc/tau =",tpcot
	    WRITE(10,*) "tau_theta/tau  =",ttot
	    WRITE(10,*) "tau_c/tau  =",tcot
	    WRITE(10,*) "tau_ctheta/tau  =",tctot
	    WRITE(10,*) " "
	    WRITE(66,*) "p1 =",p1
	    WRITE(66,*) "p2 =",p2
	    WRITE(66,*) "p3 =",p3
	    WRITE(66,*) "p4 =",p4
	    WRITE(66,*) "p5 =",p5
	    WRITE(66,*) "p6 =",p6
	    WRITE(66,*) "p7 =",p7
	    WRITE(66,*) "p8 =",p8
	    WRITE(66,*) "p9 =",p9
	    WRITE(66,*) "p10=",p10
	    WRITE(66,*) "p11=",p11
	    WRITE(66,*) "a0 =",a0
	    WRITE(66,*) "a1 =",a1
	    WRITE(66,*) "a2 =",a2
	    WRITE(66,*) "a3 =",a3
	    WRITE(66,*) "a4 =",a4
	    WRITE(66,*) "a5 =",a5
	    WRITE(66,*) "d0 =",d0
	    WRITE(66,*) "d1 =",d1
	    WRITE(66,*) "d2 =",d2
	    WRITE(66,*) "d3 =",d3
	    WRITE(66,*) "d4 =",d4
	    WRITE(66,*) "d5 =",d5
	    WRITE(66,*) "d6 =",d6
	    WRITE(66,*) "d7 =",d7
	    WRITE(66,*) "d8 =",d8
	    WRITE(66,*) "d9 =",d9
	    WRITE(66,*) "d10=",d10
	    WRITE(66,*) "d11=",d11
	    WRITE(66,*) "d12=",d12
	    WRITE(66,*) "d13=",d13
	    WRITE(66,*) "d14=",d14
	    WRITE(66,*) "d15=",d15
            CLOSE(66)
	  END IF
!*****C
       endif
	END IF
	ifrecall = 1
!******

!980610 Modification of section of "sx" containing the den and nums of the "S"'s

!###############################################

         D = d0 + d1*y*n**2 + d2*y*n*c + d3*y*c**2 + d4*n**3 + d5*n**2*c &
      + d6*n*c**2 + d7*c**3 &
      + d8*y*n + d9*y*c + d10*n**2 + d11*n*c + d12*c**2 + d13*y &
      + d14*n + d15*c

!########################################################################
         Nm = a0 + a1*n**2 + a2*n*c + a3*c**2 + a4*n + a5*c

!###########################################################################


         Nh  = - (30.D0*n*p6 + 30.D0*c*p4 - 60.D0 &
         - ( 2.D0*p1m  + 15.D0*p2m**2 - 6.D0*p2m  - 5.D0*p1m**2 ) &
         * y ) &
         * (c*p4*p7 - c*p11 - n*p11 + 1.D0)


         Nc  =   (30.D0*n*p6 + 30.D0*c*p4 - 60.D0 &
         - ( 2.D0*p1m  + 15.D0*p2m**2 - 6.D0*p2m  - 5.D0*p1m**2 ) &
         * y ) &
         * (c*p10 - 1.D0 - n*p6*p7 + n*p9)


!980610-15 Modification of section of "sx" containing Sm, Sh, Sc

!*******************************************************************************
         Sm = (4.D0/15.D0) * tpvot * Nm/D

         Sh = (4.D0/15.D0) * tptot * Nh/D

         Sc = (4.D0/15.D0) * tpcot * Nc/D
!*******************************************************************************


      return
 1004 format(12(1pe14.5))
      end subroutine smshsc_a3

   


                                                                        
    subroutine rtwi(x,val,fct,xst,eps,iend,ier)              !checked        
!     to solve general nonlinear equations of the form x=fct(x)       
!     by means of wegsteins iteration method                         
!     prepare iteration                                             

      use param_mod, only: mytid
	  implicit none
	  
	  integer,intent(in)::iend
	  real(kind=r8),intent(in)::xst,eps
	  integer,intent(out)::ier
	  real(kind=r8),intent(out)::x
	  
	  
	 real(kind=r8),external::fct
	  
	  integer::i
	  real(kind=r8)::val
	  real(kind=r8)::a,b,tol,d
	  
	  
	  

      ier=0                                                        
      tol=xst                                                     
      x=fct(tol)                                                 
      a=x-xst                                                   
      b=-a                                                     
      tol=x                                                   
      val=x-fct(tol)                                         
!     start iteration loop                                 
      do 6 i=1,iend                                       
!981103 Crude fix to avoid mysterious problem which occurred with a close but not too close guess.

      IF(DABS(val).LT.1.D-12) val =0.D0



      if(val) 1,7,1                                      
!     equation is not satisfied by x                    
 1    b=b/val-1.D0                                       
      if(b) 2,8,2                                     
!     iteration is possible                          
 2    a=a/b                                         
      x=x+a                                        
      b=val                                       
      tol=x                                      
      val=x-fct(tol)                            
!     test on satisfactory accuracy            
      tol=eps                                 
      d=abs(x)                               
      if(d-1.d0) 4,4,3                        
 3    tol=tol*d                            
 4    if(abs(a)-tol) 5,5,6                
 5    if(abs(val)-10.D0*tol) 7,7,6         
 6    continue                          
!     end of iteration loop                                           
!     no convergence after iend iteration steps. error return.       
      ier=1                                          
 7    return                                        
!     error return in case of zero divisor         
 8    ier=2                                       
      return                                     
      end   subroutine rtwi                                     



	subroutine oursal2_zeroshear(and2on2,amtaun2,sm,sh,sc)  !checked

      use param_mod, only: mytid
	  use global_t
	  
	  implicit none
	  real(kind=r8),intent(in)::and2on2
	  real(kind=r8),intent(out)::sm,sh,sc,amtaun2
	  
	  real(kind=r8)::eps,anh2on2,ans2on2,taunh2,tauns2,c_y,c_n,c_c,ss
	  integer::iend,ier
	  

!
!030321-27AH SUBROUTINE TO CALCULATE THE TURBULENCE MODEL FOR THE ZERO SHEAR 
!AH	  *UNSTABLE CASE*. TO DIMENSIONALIZE DIFFUSIVITIES IN THAT CASE NEED
!AH	  THE POSITIVE QUANTITY -(\tau N)^2 WHICH I NAME "amtaun2".
! --- Submodule to calculate turbulence functions, -(\tau N)^2 and S_M,S_H,S_S ,
! --- of the difference of the ratios of the heat and salt contributions 
! --- to the square of the Brunt Vaisala frequency, N_d^2/N^2 , named "and2on2".
! --- Calls quad_sal_pureconv.
! --- Stripped down and altered from the winter 2003 hycom version of oursal2_1a.
! --- Version in which following OTsalche/plot000127 
! --- the timescale ratios are calculated in the 'smshsc' routine 
! --- and passed back through the common block bb0/
! --- to simplify the process of adjustment of timescale ratios.
! --- Stripped and adapted from plot981007.f.
! --- Program to generate contour and 1 variable plots vs. Ri,Ri_d based on
! --- Program to generate contour plots vs. Ri_T and Ri_C based on
! --- Program to generate plots vs. Ri_T at different Ri_C values based on
! --- CORRECTED PROGRAM WITH NEW VALUE OF "p10". 'p10 = tpt*tct/(tc**2)'
! --- Program to generate K_X/((l^2) S) for Canuto based on plot980609.f:
! --- Program to generate data for plots of turbulence functions including
! --- S_{M,H,C} and Canuto's new y = (\tau_pv S)^2
! --- and n,c as functions of stability parameters in the concentration theory
! --- (structure is a 1 point closure like the generalized Mellor-Yamada, 
! --- but the constants are derived based on Dubovikov's model according
! --- to Ye Cheng). The concentration theory dimensionless parameters
! --- associated with the squares of shear, temperature contribution to 
! --- Brunt Vaisala frequency and concentration contribution to it,
! --- the new y,n,c are represented in this program by the variables
! --- c_y,c_n,c_c. 
! --- Adapted from Cheng's program mike_12.f_980528 for the Dubovikov model.
!-----------------------------------------------------------------------
! --- Inputs:
! --- and2on2  	N_d^2/N^2   Difference of heat and salt fractions of N^2
! --- 
! --- Internal quantities:
! --- anh2on2 	N_H^2/N^2   Fraction of N^2 contributed by temperature gradient
! --- ans2on2 	N_S^2/N^2   Fraction of N^2 contributed by salinity gradient
! --- taunh2    (\tau N_H)^2 dimensionless variable for temperature gradient
! --- tauns2    (\tau N_S)^2 dimensionless variable for salinity gradient
! ---
! --- Outputs:
! --- amtaun2			-(\tau N)^2 
! --- sm			S_M
! --- sh			S_H
! --- ss			S_S
!-----------------------------------------------------------------------
!
! --- kx=e*tau*sx=1/2*(b1*l)**2*sqrt(-n^2)/(-(\tau n^2))**(1/2)*sx
!
! --- X = {M,H,C} . 
! --- The program variable "amtaun2" is -(\tau N^2). 
! --- Since \tau=B_1 l/q, (N l/q)^2 = (\tau N^2)/(B1^2) .
! --- Interested in the unstable case so (\tau N^2) is negative.
!
!     parameter(sgmt=0.72)    !Make "sgmt" a parameter.
!    .                          !Standard value was 0.72.
!     parameter(tptot0=(1./5.)*(1./(1.+(1./sgmt)))) ! \tau_p\theta over \tau
!     parameter(tpcot0=tptot0)                                !tau_pc over \tau
!     parameter( ttot0=sgmt)                                  !tau_\theta over \tau
!     parameter( tcot0=ttot0)                                 !tau_c over \tau
!     parameter(tctot0=1./3.)                             ! tau_c\theta } over \tau
!     parameter(tptot = tptot0)
!     parameter(tpcot = tpcot0)
!     parameter(ttot  = ttot0)
!     parameter(tcot  = tcot0)
!     parameter(tctot = tctot0)
!
!
!
! --- Commented excerpt from the file "sx"
!
!
! --- sgmt := 0.72;
!
! --- tpt := 1/(5*(1+1/sgmt))*tau;
! --- tpt  = .08372093019*tau
!
! --- tpc := 1/(5*(1+1/sgmt))*tau;
! --- tpc  = .08372093019*tau
!
! --- tt := sgmt*tau;
! --- tt   = .72*tau
!
! --- tc := sgmt*tau;
! --- tc   = .72*tau
!
! --- tct := 2/15*sgmt*tau;
! --- tct  = .09599999998*tau
!
! --- Calculate the timescale ratios in the 'smshsc' routine instead of here.
! --- Set \sigma_t0. sgmt = .72
!
! --- Calculate {\tau_C \over \tau} and {\tau_{C\theta} \over \tau}.
! --- tcot  = sgmt
! --- tctot = (2./15.)*sgmt
! --- "tpt/tau" and "tpc/tau" from the "sx" excerpt
! --- tptot = 1./(5.*(1+1/sgmt))
! --- tpcot = 1./(5.*(1+1/sgmt))
!
!980610-030403 Common block with ratios of timescales
 !     COMMON /bb0/ ttot,tcot,tctot,tptot,tpcot,tpvot
! --- Timescale ratios are now calculated in the 'smshsc' subroutine.
! --- Make dummy call with c_y=c_n=c_c=0 to get their values for initial use.

      call smshsc_a3(0.D0,0.D0,0.D0,sm,sh,sc)



!
      eps=1.E-6                                              
      iend=300                                              

	

	 anh2on2 = (1. + and2on2)/2.
	 ans2on2 = (1. - and2on2)/2.
!030321AH **SOLVE A QUADRATIC EQUATION TO FIND -(\tau N)^2 **
! ---	  **FOR PURE TURBULENT CONVECTION INCLUDING SALINITY (SHEAR=0).**
	 CALL QUAD_SAL_PURECONV(anh2on2,ans2on2,amtaun2,ier)
!
         if(ier.NE.0) then
       if (mytid.eq.1) then
	    WRITE(10,*) "In oursal2_zeroshear subroutine"
	    WRITE(10,*) "anh2on2",anh2on2,"	ans2on2=",ans2on2
            WRITE(10,*) "Error returned by quad_sal_pureconv ier=",ier
	    WRITE(10,*) "ier=1 means choice of root must be reconsidered" &
               //"for these model constants."
       endif
            STOP
         endif
!
!030322AH-0404Y CALCULATE TEMPERATURE AND SALINITY GRADIENT DIMENSIONLESS TURBULENCE FUNCTIONS.
	 taunh2 = -amtaun2*anh2on2
	 tauns2 = -amtaun2*ans2on2

!
	 c_y = 0.
         c_n = -(tcot*tctot)*taunh2
         c_c = -(tcot**2)*tauns2
         call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)
!030304Y Stop in case any of the S_X is negative.

 1003 format(12(I8))
 1004 format(12(1pe14.5))
      end subroutine oursal2_zeroshear





    SUBROUTINE quad_sal_pureconv(anh2on2,ans2on2,amtaun2,ier)         !checked                     

      use param_mod, only: mytid
	  use global_t
	  implicit none
	  
      real(kind=r8),intent(in)::anh2on2,ans2on2
	  integer,intent(out)::ier
      real(kind=r8),intent(out)::amtaun2
      real(kind=r8)::rootplus,rootminus
	  
	  real(kind=r8),parameter::errorbound=1.D-12
	  real(kind=r8),parameter::amtaun2_negri_0rid= (16.6**2)* &
          (-(6.5489662174209907D-04)*(-56.94668746578522))
		  
	!  integer::ier
      real(kind=r8)::eta3p,amu3p,braasal,braatem,para,a3p,brabsal,brabtem,brab,b3p,anum,bnum,cnum,braa,parb
	  real(kind=r8)::radical,qnum,errorplus,errorminus,errorproplus,errorprominus
	  



	ier=0
	amtaun2=0. 	!Initialize the output variable.
	


      eta3p = tpcot*(tctot*anh2on2 + tcot*ans2on2)
      amu3p = tptot*(tctot*ans2on2 + ttot*anh2on2)


      braasal = tpcot*(amu3p  - tptot*tctot*anh2on2)*ans2on2
      braatem = tptot*(eta3p - tpcot*tctot*ans2on2)*anh2on2
      braa = braasal+braatem
      para = (eta3p*amu3p - tptot*tpcot*(tctot**2)*anh2on2*ans2on2)
      a3p   = -(7./15.)*braa - para

      brabsal = tpcot*ans2on2
      brabtem = tptot*anh2on2
      brab = brabsal+brabtem
      parb = (eta3p + amu3p)
      b3p   = -(7./15.)*brab - parb

!-25     Solution of quadratic for -(\tau N)^2
	
      anum = a3p
      bnum = -b3p
      cnum = -1.
      radical = SQRT(bnum**2 - 4*anum*cnum)
      qnum = (-1./2.)*(bnum + SIGN(radical,bnum))
      rootplus  = qnum/anum
      rootminus = cnum/qnum
!030325-27 Check that calculated solutions actually satisfy the quadratic equation.
!030327 Make acceptable error a small fraction of the root size.
      errorplus  = a3p*(rootplus**2)-b3p*(rootplus)-1.
      errorminus = a3p*(rootminus**2)-b3p*(rootminus)-1.
      errorproplus  = errorplus/rootplus
      errorprominus = errorminus/rootminus
      IF(MAX(ABS(errorproplus),ABS(errorprominus)).GT.errorbound) THEN
       if (mytid.eq.1) then
	WRITE(10,*) "PROBLEM: Error too great in calculating amtaun2 ."
	WRITE(10,*) "In quad_sal_pureconv:"
	WRITE(10,*) "-(\\tau N)^2=q/a \\equiv first root:",rootplus
	WRITE(10,*) "A''' (- \\tau N)^2)^2 - B'''((-\\tau N)^2) - 1 =", &
             errorplus
	WRITE(10,*) "proportional error= ",errorproplus
	WRITE(10,*) "-(\\tau N)^2=c/q \\equiv second root:",rootminus
	WRITE(10,*) "A''' (- \\tau N)^2)^2 - B'''((-\\tau N)^2) - 1 =", &
             errorminus
	WRITE(10,*) "proportional error= ",errorprominus
	WRITE(10,*) "Theoretically both roots should give zero."
	WRITE(10,*) "Maximum proportional error we permit =",errorbound
	WRITE(10,*) "a=",anum," b=",bnum," c=",cnum
	WRITE(10,*) "q=",qnum
	WRITE(10,*) "Program is stopping."
       endif
	STOP
      END IF
!******
!030326-0402Y     I choose "rootminus" because it approximates the solution with shear
!     at very negative Ri and because it bends downward as (N_H^2/N_d^2)
!     departs from zero which is consistent with an additional driving of
!     mixing due to double-diffusivity even in the unstable case.
!     IF CHANGE MODEL CONSTANTS *MUST* REVISIT THE CHOICE OF ROOT.
      amtaun2 = rootminus
      IF(ABS(anh2on2-ans2on2).LT.errorbound) THEN
        IF(ABS(rootminus-amtaun2_negri_0rid).GT.1.) THEN 	
       if (mytid.eq.0) then
	  WRITE(10,*) "Check if have the right root."
	  WRITE(10,*) "rootminus=",rootminus
	  WRITE(10,*) "amtaun2_negri_0rid=",amtaun2_negri_0rid
       endif
	  ier=1
        END IF
      END IF
!
      return                                          
      end    SUBROUTINE quad_sal_pureconv                                        
                                               
      

    function fctysal(y)               !checked               
!990614-15 Function for y \equiv (\tau S)^2 from the Production=Dissipation equation.
!980501 Make double precision to conform to calling cctmix routine.

      use param_mod, only: mytid
	  implicit none
	  real(kind=r8),intent(in)::y
	  real(kind=r8)::x,z,sm,sh,ss
	  
	  real(kind=r8)::ri,rid,ric,rit,rr
      real(kind=r8)::fctysal
      common /mikebb/ri
      common /mikebbs/rid,ric,rit,rr
	  
	  
      x=rit*y
      z=ric*y
      call mksmshss1(x,y,z,sm,sh,ss)
!990614 Use the traditional Production=Dissipation equation:
!	y(S_M - Ri_T S_H - Ri_C S_C) = 2            
      fctysal=2.D0/(sm-rit*sh-ric*ss)
      return                                          
      end     function fctysal                                       


                                                                        
    subroutine mksmshss1(x,y,z,sm,sh,ss)  !checked

      use param_mod, only: mytid
      implicit none
	  
	  integer,parameter::ifnew=0
	  real(kind=r8),parameter::x1min=15.d0/8.d0,calculation_epsilon=1.D-60
	  
	  real(kind=r8)::phi,x,y,z,rit,ric,x1,fc_twid,ri,fc,psi,eta_t,eta_s,eta_sum,sigma_m,st
	  
	  

	  real(kind=r8),intent(out)::sm,sh,ss
	  integer::ifrecall=0
	  
	  



	IF(ifrecall.EQ.0) THEN
       if (mytid.eq.0) then
	  WRITE(10,*) " "
	  WRITE(10,*)"************************************************"
	  WRITE(10,*) " "
	  WRITE(10,*)"Regularization used for Dubovikov salinity model"
	  WRITE(10,*) "ifnew=",ifnew
	  IF (ifnew.EQ.0) THEN
	    WRITE(10,*) "1/1+x regularization"
	  ELSE IF(ifnew.EQ.1) THEN
	    WRITE(10,*) "e^-x regularization"
	  END IF  
	  WRITE(10,*) " "
	  WRITE(10,*)"************************************************"
	  WRITE(10,*) " "
        endif
	END IF
!*****************************************************************************C

      rit=x/y
      ric=z/y
      y = MAX(y,calculation_epsilon)
      x = rit*y
      z = ric*y
!*****C
      x1=0.5d0*sqrt(y)
      if(x1.lt.x1min) then
         fc_twid=5.d0/9.d0
      else
         fc_twid=1.d0-5.d0/(3.d0*x1)+25.d0/(16.d0*x1**2)
      endif
!990615 Total Richardson number, Ri = Ri_T + Ri_C .
      ri = rit + ric

      fc = 1.8D0*fc_twid

	IF(ifnew.EQ.0) THEN
	  phi = 1.D0/(1.D0 + 0.08*x1**2*ri)
        psi = 1.D0/ &
            (  1.D0 & 
            + 0.16D0*x1**2*ri &
             + 0.013D0*x1**4*rit*ric*phi )
	ELSE IF(ifnew.EQ.1) THEN 	!990923 version
	  phi =  EXP(-0.08*x1**2*ri)
          psi = EXP (-0.16D0*x1**2*ri  - 0.013D0*x1**4*rit*ric*phi)
	END IF
!*****C
!990621 Alternate equation 
      eta_t = 0.05D0*x1**2*rit*(1.D0 + (1.D0 +0.16D0*ric*x1**2)*psi) 
      eta_s = 0.05D0*x1**2*ric*(1.D0 + (1.D0 +0.16D0*rit*x1**2)*psi) 
      eta_sum = eta_t + eta_s
      IF(ifnew.EQ.0) THEN
	sigma_m = 1.D0/(1.D0 + eta_sum/fc)
      ELSE IF(ifnew.EQ.1) THEN		!990923 version
	sigma_m = EXP(-eta_sum/fc)
      END IF
!*****C
!990614 Equation (65a) of 990610 paper excerpt: 
!       S_M = {1 \over 25} F_1 F\twiddle_C \sigma_m
!	Set F_1 = 1 .
      sm=1.d0/25.d0*fc_twid*sigma_m
!990621 S_H = 0.056[1 + 0.08 x1**2 Ri_C - 0.006 x1**4 Ri_C Ri] \psi	
        sh=0.056d0*(1.D0 + 0.08D0*x1**2*ric*phi)*psi
        ss=0.056d0*(1.D0 + 0.08D0*x1**2*rit*phi)*psi
!*****C
!990924 Add Concentration Smagorinsky-like function for 990923 version.
!	'S_\tau = 0.056[1 + 0.08 {x_DB}^2 {Ri_DB} (1 - R_\rho) \Phi] \psi_1'
!       S_\tau = 0.056 (1 + 0.08 {x_1}^2 Ri_T (1 + Ri_C/Ri_T) \Phi] \psi_1
!       S_\tau = 0.056 (1 + 0.08 {x_1}^2 Ri \Phi] \psi_1
      IF(ifnew.EQ.1) THEN
        st=0.056d0*(1.D0 + 0.08D0*x1**2*ri*phi)*psi
      END IF
	

      return
      end subroutine mksmshss1

	

	SUBROUTINE interp1d_expabs(x, &        !checked
                    x_1,slq2_1,sm_1,sh_1,ss_1, &
                    slq2,sm,sh,ss, &
                    ixmax,m,m0,delta,rat)
!030129-030328 Subroutine for a faster interpolation calculation in the ifexpabstable=1 case.
! --- 6-990112-C030326	Subroutine for a ONE VARIABLE modular interpolation calculation.
!	x is the independent variable in the turbulence model calculation.
!
! --- 1D input array with table spacing: 	x_1
! --- table limit value:			ixmax
! --- 1D input arrays with table values:	slq2_1,sm_1,sh_1,ss_1
! --- Output interpolated values:		slq2,sm,sh,ss
!01107yXI Input value of ratio between entries  rat


      use param_mod, only: mytid
	  implicit none
	  
	  integer::ixmax,m,m0
	  real(kind=r8)::x
	  real(kind=r8),intent(in)::delta,rat
	  real(kind=r8)::x_1(-m:m),sm_1(-m:m),ss_1(-m:m),sh_1(-m:m),slq2_1(-m:m)
	  
	  
	  real(kind=r8)::tabindx,deltaxta,dslq2_x,dsm_x,dsh_x, dss_x,deltax
      real,intrinsic::SIGN
      integer::lx1,lx0  
      real(kind=r8),intent(out)::sm,sh,ss,slq2

!030326 Take values off the edge of the table back to the table.
!     print *,'start interp1d_expabs'
	IF(x.GT.x_1(ixmax)) THEN
	    x = x_1(m)
	ELSE IF(x.LT.x_1(-m)) THEN
	    x = x_1(-m)
	END IF
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(x).LT.x_1(m0)) THEN
!981019-27 Find Interpolation points in the equally spaced x part of the table.

	  lx1 = INT(x/delta)+NINT(SIGN(DFLOAT(1),x))



!011107-08yXI Find Interpolation points in exponential absolute value spaced x part of the table.
	ELSE IF((ABS(x)).GE.(x_1(m))) &
   THEN
!981103 Special case where have a value which falls at the limit of the table.

	  lx0 = NINT(SIGN(DFLOAT(m),x))



	  lx1 = lx0
	  GO TO 252
!*****C
	ELSE

	  tabindx = SIGN( DFLOAT(m0) + ((LOG(ABS(x)) - LOG(x_1(m0)))/LOG(rat)), x)




	  lx1 = INT(tabindx)+NINT(SIGN(DFLOAT(1),x))



!yXI
	END IF
!011108yXI It is conceivable that rounding errors may in borderline cases 
!	   throw the calculated table indices for x off by one.
!	   Check and allow moving one to either side to take care of this.
 	IF((ABS(x_1(lx1))).LT.(ABS(x))) THEN
	  lx1 = lx1 + SIGN(1,lx1)
	ELSE IF((ABS(x_1(lx1-SIGN(1,lx1)))).GT.(ABS(x))) THEN
	  lx1 = lx1 - SIGN(1,lx1)
	END IF
!yXI
  250	CONTINUE
!981019-27 Make lx0 one less or greater than lx1 according to sgn(x).

        lx0 = lx1 - NINT(SIGN(DFLOAT(1),x)) 



!*****C
        IF(x.EQ.0.D0) lx1 = 1
  252   CONTINUE
!981019-28 Check that the x value falls within the interpolation interval.
	IF((x.GT.0.D0.AND.  &
      (x.LT.x_1(lx0).OR.x.GT.x_1(lx1))).OR.  &
	   (x.LT.0.D0.AND.  &
      (x.GT.x_1(lx0).OR.x.LT.x_1(lx1)))) THEN
       if (mytid.eq.1) then
	   WRITE(10,*) &
      "x is outside interpolation range in interp1d_expabs."
	   WRITE(10,*) "delta=",delta
	   WRITE(10,*) "m0=",m0," m=",m," rat=",rat
	   WRITE(10,*) "x=  ",x," lx0= ",lx0," lx1= ",lx1
	   WRITE(10,*) "x_1(lx0)=  ",x_1(lx0), &
                "   x_1(lx1)= ",x_1(lx1)
	   WRITE(10,*) "Program is stopping."
        endif
	   STOP
	END IF
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
	deltaxta = x_1(lx1) - x_1(lx0)
	deltax = x - x_1(lx0)
!	slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
	IF(lx1.EQ.lx0) THEN
	  dslq2_x = 0.D0
	ELSE
	  dslq2_x = (slq2_1(lx1) - slq2_1(lx0))/ &
                deltaxta
	END IF
	slq2	 = slq2_1(lx0) + &
             dslq2_x*deltax
!	sm
	IF(lx1.EQ.lx0) THEN
	  dsm_x   = 0.D0
	ELSE
	  dsm_x = (sm_1(lx1) - sm_1(lx0))/ &
              deltaxta
	END IF
	sm     = sm_1(lx0) + &
           dsm_x*deltax
!	sh
	IF(lx1.EQ.lx0) THEN
	  dsh_x   = 0.D0
	ELSE
	  dsh_x = (sh_1(lx1) - sh_1(lx0))/ &
              deltaxta
	END IF
	sh     = sh_1(lx0) + &
           dsh_x*deltax
!	ss
	IF(lx1.EQ.lx0) THEN
	  dss_x   = 0.D0
	ELSE
	  dss_x = (ss_1(lx1) - ss_1(lx0))/ &
              deltaxta
	END IF
	ss     = ss_1(lx0) + &
           dss_x*deltax
!*****C


!     print *,'finish interp1d_expabs'
	RETURN
	END SUBROUTINE interp1d_expabs


	

	SUBROUTINE interp2d_expabs(ri,rid, &   !checked
                    ri_1,rid_1,slq2_2,sm_2,sh_2,ss_2, &
                    slq2,sm,sh,ss, &
                    irimax,m,m0,delta,rat)


!01107yXI Input value of ratio between entries  rat
!	1D input arrays with table spacing: 	ri_1,rid_1
!	1D input array with table limits:	irimax
!	2D input arrays with table values:	slq2_2,sm_2,sh_2,ss_2
!	Output interpolated values:		slq2,sm,sh,ss


      use param_mod, only: mytid
	  implicit none
	  
	  
	  integer::m,m0,irimax(-m:m),lrid1,lrid0,lri1,lri0
	  real(kind=r8),intent(in)::delta,rat
	  real(kind=r8)::ri,rid,tabindrid, tabindri
	  real(kind=r8):: ri_1(-m:m),rid_1(-m:m),slq2_2(-m:m,-m:m),sm_2(-m:m,-m:m),sh_2(-m:m,-m:m),ss_2(-m:m,-m:m)
	  real(kind=r8)::deltaridta,deltarita,deltari,deltarid,dslq2_rid,dslq2_ri,dsh_rid,dsh_ri,dss_ri,dss_rid,dsm_ri,dsm_rid
	  real(kind=r8),intent(out)::ss,sm,sh,slq2
	  




! 	Take values off the edge of the table back to the table on radial lines.
!981103 Must use ratio of Ri's taken *before* the cut-off has taken place.
	IF(ri.GT.ri_1(m)) THEN
	  IF(ABS(rid).LE.ri) THEN
	    rid = ri_1(m)*(rid/ri)
	    ri  = ri_1(m)
	  ELSE IF(rid.GT.ri) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	  ELSE IF(rid.LT.-ri) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	  END IF
	ELSE IF(ri.LT.ri_1(-m)) THEN
	  IF(ABS(rid).LE.-ri) THEN
	    rid = ri_1(-m)*(rid/ri)
	    ri  = ri_1(-m)
	  ELSE IF(rid.GT.-ri) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	  ELSE IF(rid.LT.ri) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	  END IF
	ELSE IF(rid.GT.rid_1(m)) THEN
	    ri  = rid_1(m)*(ri/rid)
	    rid = rid_1(m)
	ELSE IF(rid.LT.rid_1(-m)) THEN
	    ri  = rid_1(-m)*(ri/rid)
	    rid = rid_1(-m)
	END IF
!*****C
!981019 Interpolate points within the table range.
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(rid).LT.rid_1(m0)) THEN
	  lrid1 = INT(rid/delta)+NINT(SIGN(DFLOAT(1),rid))
	ELSE IF((ABS(rid)).GE.(rid_1(m)))  THEN
	  lrid0 = NINT(SIGN(DFLOAT(m),rid))
	  lrid1 = lrid0
	  GO TO 252
!*****C
	ELSE
	  tabindrid = SIGN(DFLOAT(m0) + ((LOG(ABS(rid)) - LOG(rid_1(m0)))/LOG(rat)), rid)
	  lrid1 = INT(tabindrid)+NINT(SIGN(DFLOAT(1),rid))

	END IF

 	IF((ABS(rid_1(lrid1))).LT.(ABS(rid))) THEN
	  lrid1 = lrid1 + SIGN(1,lrid1)
	ELSE IF((ABS(rid_1(lrid1-SIGN(1,lrid1)))).GT.(ABS(rid))) THEN
	  lrid1 = lrid1 - SIGN(1,lrid1)
	END IF

  250	CONTINUE
!981019-27 Make lrid0 one less or greater than lrid1 according to sgn(rid).

        lrid0 = lrid1 - NINT(SIGN(DFLOAT(1),rid)) 



!*****C
        IF(rid.EQ.0.D0) lrid1 = 1
  252   CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((rid.GT.0.D0.AND.  &
      (rid.LT.rid_1(lrid0).OR.rid.GT.rid_1(lrid1))).OR.  &
	   (rid.LT.0.D0.AND.  &
      (rid.GT.rid_1(lrid0).OR.rid.LT.rid_1(lrid1)))) THEN
       if (mytid.eq.1) then
	   WRITE(10,*) "Ri_d is outside interpolation range in interp2d."
	   WRITE(10,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
	   WRITE(10,*) "rid_1(lrid0)=  ",rid_1(lrid0), &
                "   rid_1(lrid1)= ",rid_1(lrid1)
	   WRITE(10,*) "Program is stopping."
        endif
	   STOP
	END IF
!*****C
!C981022	Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
!C	This is to conform to the 1D table's realizability limit treatment. 
!	IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
!	  ri = MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))
!	END IF
!C*****C
!981110 Set turbulence to zero if Ri threatens to surpass the realizability limit.
        IF(ri.GT.MIN(ri_1(irimax(lrid0)),ri_1(irimax(lrid1)))) THEN
          slq2=0.D0
          sm = 0.D0
          sh = 0.D0
          ss = 0.D0
          RETURN
        END IF
!*****C
!981022 Table index ranges from -m to m with equal spacing for -m0 to m0.
	IF(ABS(ri).LT.ri_1(m0)) THEN
!981019-27 Find Interpolation points in the equally spaced Ri part of the table.

	  lri1 = INT(ri/delta)+NINT(SIGN(DFLOAT(1),ri)) 



!011107-08yXI Find Interpolation points in exponential absolute value spaced Ri part of the table.
	ELSE IF((ABS(ri)).GE.(ri_1(m))) &
        THEN
!981103 Special case where have a value which falls at the limit of the table.

	  lri0 = NINT(SIGN(DFLOAT(m),ri))



	  lri1 = lri0
	  GO TO 272
!*****C
	ELSE

	  tabindri = SIGN( DFLOAT(m0) + ((LOG(ABS(ri)) - LOG(ri_1(m0)))/LOG(rat)), ri)




	  lri1 = INT(tabindri)+NINT(SIGN(DFLOAT(1),ri))



!yXI
  270	CONTINUE
	END IF
!011108yXI It is conceivable that rounding errors will in borderline cases 
!	   throw the calculated table indices for Ri off by one.
!	   Check and allow moving one to either side to take care of this.
 	IF((ABS(ri_1(lri1))).LT.(ABS(ri))) THEN
	  lri1 = lri1 + SIGN(1,lri1)
	ELSE IF((ABS(ri_1(lri1-SIGN(1,lri1)))).GT.(ABS(ri))) THEN
	  lri1 = lri1 - SIGN(1,lri1)
	END IF
!yXI
!981019-27 Make lri0 one less or greater than lri1 according to sgn(ri).

        lri0 = lri1 - NINT(SIGN(DFLOAT(1),ri)) 



!*****C
        IF(ri.EQ.0.D0) lri1 = 1
  272	CONTINUE
!981019-27 Check that the Ri_d value falls within the interpolation interval.
	IF((ri.GT.0.D0.AND.(ri.LT.ri_1(lri0).OR.ri.GT.ri_1(lri1))).OR.  &
	   (ri.LT.0.D0.AND.(ri.GT.ri_1(lri0).OR.ri.LT.ri_1(lri1)))) THEN
       if (mytid.eq.1) then
	   WRITE(10,*) "Ri is outside interpolation range in interp2d."
	   WRITE(10,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
	   WRITE(10,*) "ri_1(lri0)=  ",ri_1(lri0), &
                "   ri_1(lri1)= ",ri_1(lri1)
	   WRITE(10,*) "Program is stopping."
       endif
	   STOP
	END IF
!*****C
!981019-27 Interpolate turbulence fields.
!981027-990112 Introduce table spacing variables.
	deltaridta = rid_1(lrid1) - rid_1(lrid0)
	deltarita  = ri_1(lri1)  - ri_1(lri0)
	deltarid = rid - rid_1(lrid0)
	deltari  = ri - ri_1(lri0)
!	slq2
!981103 Set delta field to zero in special cases falling at limit of the table. 
	IF(lrid1.EQ.lrid0) THEN
	  dslq2_rid = 0.D0
	ELSE
	  dslq2_rid = (slq2_2(lri0,lrid1) - slq2_2(lri0,lrid0))/ &
                deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dslq2_ri  = 0.D0
	ELSE
	  dslq2_ri = (slq2_2(lri1,lrid0) - slq2_2(lri0,lrid0))/ &
               deltarita
	END IF
	slq2	 = slq2_2(lri0,lrid0) + &
            dslq2_ri*deltari + dslq2_rid*deltarid
!	sm
	IF(lrid1.EQ.lrid0) THEN
	  dsm_rid   = 0.D0
	ELSE
	  dsm_rid = (sm_2(lri0,lrid1) - sm_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsm_ri    = 0.D0
	ELSE
	  dsm_ri = (sm_2(lri1,lrid0) - sm_2(lri0,lrid0))/ &
             deltarita
	END IF
	sm     = sm_2(lri0,lrid0) + &
            dsm_ri*deltari + dsm_rid*deltarid
!	sh
	IF(lrid1.EQ.lrid0) THEN
	  dsh_rid   = 0.D0
	ELSE
	  dsh_rid = (sh_2(lri0,lrid1) - sh_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dsh_ri    = 0.D0
	ELSE
	  dsh_ri = (sh_2(lri1,lrid0) - sh_2(lri0,lrid0))/ &
              deltarita
	END IF
	sh     = sh_2(lri0,lrid0) + &
            dsh_ri*deltari + dsh_rid*deltarid
!	ss
	IF(lrid1.EQ.lrid0) THEN
	  dss_rid   = 0.D0
	ELSE
	  dss_rid = (ss_2(lri0,lrid1) - ss_2(lri0,lrid0))/ &
              deltaridta
	END IF
	IF(lri1.EQ.lri0) THEN
	  dss_ri    = 0.D0
	ELSE
	  dss_ri = (ss_2(lri1,lrid0) - ss_2(lri0,lrid0))/ &
             deltarita
	END IF
!*****C
	ss     = ss_2(lri0,lrid0) + &
            dss_ri*deltari + dss_rid*deltarid
!*****C


	RETURN
	END SUBROUTINE interp2d_expabs
	
END MODULE vmix_canuto

