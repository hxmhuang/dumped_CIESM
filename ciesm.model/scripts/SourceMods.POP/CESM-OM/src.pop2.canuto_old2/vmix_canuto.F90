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
    use const_for_nasagissmixing, only : hwy_density_interface, &
                        hwy_alpha_beta_interface

	  
	public :: vmix_coeffs_canuto
                	  
    contains
	
	
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


	real(kind=r8), dimension(nx_block,ny_block,km),intent(inout) ::      &
				VVC      ! viscosity for momentum diffusion

	real(kind=r8), dimension(nx_block,ny_block,km,2),intent(inout) :: &
				VDC      ! diffusivity for tracer diffusion 1 for t ;2 for s

	integer:: ii,jj,kk,na,nmax,n,bid,k,nml_error,i,j
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
	integer, parameter :: write_option = 0
	integer :: i_min, j_min, k_min
	integer :: i_max, j_max, k_max
	integer :: k_111, k_222
    character (500) :: i_string, j_string, test_data
    integer, parameter :: column_out = 0 
	  

    real (kind=r8) :: rho_up, rho_down, t_up, t_down, s_up, s_down
    real (kind=r8) :: alpha_up, alpha_down, beta_up, beta_down
    real (kind=r8) :: u_up, u_down, v_up, v_down
    real (kind=r8) :: hwy0,hwy1,hwy2,hwy3
	  
!***********************************************************************  
	write(bid_string,*) my_task 
	bid_fi_string    =   "diag_fi_"//trim(adjustl(bid_string))//".txt"
	bid_fi_r_string  =   "turb_nd2on2_"//trim(adjustl(bid_string))//".txt"
	bid = this_block%local_id

	TMIX_(:,:,:,1)=TMIX(:,:,:,1)
	TMIX_(:,:,:,2)=(TMIX(:,:,:,2)*1000.D0-35.D0)/1000.D0


	!call calculate_s2(VSHEAR,bid,UMIX,VMIX)


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

	!call state(1,1,WORK,TMIX(:,:,1,2),this_block, &
	!          RHOFULL=RHO1, DRHODT=TALPHA, DRHODS=SBETA)
    Bo(:,:) = 0.0d0
    Bosol(:,:)  =   0.0d0
    do i = 1, nx_block
    do j = 1, ny_block
        if (kmt(i,j,bid) .lt. 1) cycle
        call hwy_alpha_beta_interface(talpha(i,j),sbeta(i,j),work(i,j),TMIX(i,j,1,2),zt(1),1,1)
        Bo (i,j)    =   (grav/1.026d0)*talpha(i,j)*STF(i,j,1)
        Bosol(i,j)  =   (grav/1.026d0)*talpha(i,j)*SHF_QSW(i,j)
    end do
    end do

    an2_xy(:,:,:)  =   0.0d0
    rid_xy(:,:,:)  =   0.0d0
    do i = 1, nx_block
    do j = 1, ny_block
        if (kmt(i,j,bid) .lt. 1) cycle
        do k = 1, kmt(i,j,bid) -1
            t_up   = merge(-c2,TMIX(i,j,k,  1),TMIX(i,j,k,  1) < -c2)
            t_down = merge(-c2,TMIX(i,j,k+1,1),TMIX(i,j,k+1,1) < -c2)
            s_up   = TMIX(i,j,k,  2)*1.0d+3
            s_down = TMIX(i,j,k+1,2)*1.0d+3
            call hwy_density_interface(rho_up,  t_up,  TMIX(i,j,k,  2),  zt(k),1,1,1) 
            call hwy_density_interface(rho_down,t_down,TMIX(i,j,k+1,2),  zt(k),1,1,1) 
            call hwy_alpha_beta_interface(alpha_up,beta_up,t_up,TMIX(i,j,k,  2),zt(k),1,1)
            call hwy_alpha_beta_interface(alpha_down,beta_down,t_down,TMIX(i,j,k+1,  2),zt(k+1),1,1)
            rid_xy(i,j,k)  =   grav*alpha_down*(t_up-t_down)/dzw(k) +&
                               grav*beta_down *(s_up-s_down)/dzw(k)
            an2_xy(i,j,k)   =   grav*alpha_down*(t_up-t_down)/dzw(k) -&
                               grav*beta_down *(s_up-s_down)/dzw(k)
        end do
    end do
    end do
!-----------------------------------------------------------------------
!
!  compute turbulent and radiative sfc buoyancy forcing
!
!-----------------------------------------------------------------------
    vshear(:,:,:)  = 0.0d0
    do i = 2, nx_block 
    do j = 2, ny_block
        if (kmt(i,j,bid) .lt. 1) cycle
        do k = 1, kmt(i,j,bid) -1
            hwy0 = 0.0d0
            hwy1 = 0.0d0
            hwy2 = 0.0d0
            hwy3 = 0.0d0
            if (k .le. kmu(i,  j,    bid)) hwy0 = 1.0d0
            if (k .le. kmu(i-1,j-1,  bid)) hwy1 = 1.0d0
            if (k .le. kmu(i-1,j,    bid)) hwy2 = 1.0d0
            if (k .le. kmu(i  ,j-1,  bid)) hwy3 = 1.0d0
            u_up        = umix(i,  j,  k) *  hwy0+&
                          umix(i-1,j-1,k) *  hwy1+&
                          umix(i-1,j  ,k) *  hwy2+&
                          umix(i,  j-1,k) *  hwy3
            v_up        = vmix(i,  j,  k) *  hwy0+& 
                          vmix(i-1,j-1,k) *  hwy1+&
                          vmix(i-1,j  ,k) *  hwy2+&
                          vmix(i,  j-1,k) *  hwy3 
            u_down      = umix(i,  j,  k+1) *  hwy0+& 
                          umix(i-1,j-1,k+1) *  hwy1+&
                          umix(i-1,j  ,k+1) *  hwy2+&
                          umix(i,  j-1,k+1) *  hwy3
            v_down      = vmix(i,  j,  k+1) *  hwy0+& 
                          vmix(i-1,j-1,k+1) *  hwy1+&
                          vmix(i-1,j  ,k+1) *  hwy2+&
                          vmix(i,  j-1,k+1) *  hwy3
            u_up        =   u_up  /(hwy0+hwy1+hwy2+hwy3)
            v_up        =   v_up  /(hwy0+hwy1+hwy2+hwy3)
            u_down      =   u_down/(hwy0+hwy1+hwy2+hwy3)
            v_down      =   v_down/(hwy0+hwy1+hwy2+hwy3)
            vshear(i,j,k)  =   (u_up-u_down)**2+(v_up-v_down)**2
            vshear(i,j,k)  =   vshear(i,j,k)/(dzw(k)**2)

        end do
    end do
    end do

    vvc = 0.0d0
    vdc = 0.0d0
	
	do ii=1,nx_block
	do jj=1,ny_block

		if (kmt(ii,jj,bid) .le. 1) cycle


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
				ri (k)=an2_XY(ii,jj,k)/VSHEAR(ii,jj,k)
			else
				rid(k)=RID_XY(ii,jj,k)/1.0d-25
				ri (k)=an2_XY(ii,jj,k)/1.0d-25
			end if
		end do

	 
		do k=1,na+1
			z(k) = zw(k) 
		end do

        call linear_interpolation_zt2zw(t(1:na),t0(1:na+1),zt(1:na+1),zw(1:na),na)
        call linear_interpolation_zt2zw(s(1:na),s0(1:na+1),zt(1:na+1),zw(1:na),na)
        call linear_interpolation_zt2zw(rh(1:na),rh0(1:na+1),zt(1:na+1),zw(1:na),na)


		call nasa_giss_mixing_column(Coriol, na, z(1:na+1), &
                    t(1:na), s(1:na), an2(1:na), s2(1:na), ri(1:na), rid(1:na), &
                    ustar_,buoytur,buoysol, &
                    amld(ii,jj),akm,akh,aks, my_task)
        if (column_out .eq. 1) then
            write(i_string,"(I4.4)") ii
            write(j_string,"(I4.4)") jj
            test_data = "test_column_"//trim(adjustl(bid_string))//"_" &
                        //trim(adjustl(i_string))//"_"  &
                        //trim(adjustl(j_string))//".txt"
            open(unit=100,position='Append',file=test_data)
               write(100,*) "Coriol "
               write(100,*)  Coriol 
               write(100,*) "na "
               write(100,*)  NA

               write(100,*) "z "
               do k = 1, km 
                   write(100,*) z(k)
               end do

               write(100,*) "ustar_"
               write(100,*) ustar_ 
               write(100,*) "buoytur "
               write(100,*) buoytur 
               write(100,*) "buoysol"
               write(100,*) buoysol 
               write(100,*) "amld"
               write(100,*) amld(ii,jj)

               write(100,*) "temp"
                do k = 1, km 
                    write(100,*) t(k)
                end do

               write(100,*) "salt"
                do k = 1, km 
                    write(100,*) s(k)
                end do

               write(100,*) "an2"
                do k = 1, km
                write(100,*) an2(k)
                end do

               write(100,*) "s2"
                do k = 1, km 
               write(100,*) s2(k)
                end do

               write(100,*) "ri"
                do k = 1, km 
               write(100,*) ri(k)
                end do

               write(100,*) "rid"
                do k = 1, km 
               write(100,*) rid(k)
                end do

               write(100,*) "akm"
                do k = 1, km 
               write(100,*) akm(k)
                end do

               write(100,*) "akt"
                do k = 1, km 
               write(100,*) akh(k)
                end do

               write(100,*) "aks"
                do k = 1, km 
               write(100,*) aks(k)
                end do

            close(100)
        end if
	
		if(write_option.eq.1)then
			open(unit=10,position='Append',file=bid_fi_string)

			write(10,*), "ii, jj: ", ii, jj
			write(10,*), "lat,lon:",TLATD(ii,jj,bid),TLOND(ii,jj,bid)
			write(10,*), "na: ", na
			write(10,*), "n: ", n
			write(10,*), "nmax: ", nmax
			write(10,*), "ustar_:",ustar_
			write(10,*), "buoytur:",buoytur
			write(10,*), "buoysol:",buoysol," alpha:",TALPHA(ii,jj), " SHF_QSW:",SHF_QSW(ii,jj)/hflux_factor
			write(10,*), "Coriol:",Coriol
			write(10,*), "z(1:na+1): ", z(1:na+1)

			write(10,*), "t(1:na): ", t(1:na)
			write(10,*), "s(1:na): ", s(1:na)
			write(10,*), "rh(1:na): ", rh(1:na)

			write(10,*), "ri(1:na): ", ri(1:na)
			write(10,*), "rid(1:na): ", rid(1:na)
			write(10,*), "s2(1:na): ", s2(1:na)
			write(10,*), "an2(1:na): ", an2(1:na)


			write(10,*), "amld(ii,jj)", amld(ii,jj)
			write(10,*), "akm(1:na):", akm(1:na)
			write(10,*), "akh(1:na):", akh(1:na)
			write(10,*), "aks(1:na):", aks(1:na)
			close(10)
		end if

		vvc(ii,jj,:)=akm
		vdc(ii,jj,:,1)=akh
		vdc(ii,jj,:,2)=aks

	end do
	end do
	
	
		
    end subroutine vmix_coeffs_canuto

    function linear_interpolation_2pt(x,y,x0)
        implicit none
        real (r8), intent(in) :: x(2)
        real (r8), intent(in) :: y(2)
        real (r8), intent(in) :: x0
        real (r8) :: linear_interpolation_2pt

        real (r8) :: temp,trend

        trend = y(2)-y(1)
        trend = trend/(x(2)-x(1))
        
        temp    =   y(1)+(x0-x(1))*trend

        linear_interpolation_2pt = temp 
        return
    end function linear_interpolation_2pt

    subroutine linear_interpolation_zt2zw(x_zw,x_zt,zt,zw,na)
        implicit none
        integer, intent(in) :: na
        real (r8), intent(in) :: zt(1:na+1)
        real (r8), intent(in) :: zw(1:na)
        real (r8), intent(in) :: x_zt(1:na+1)
        real (r8), intent(out) :: x_zw(1:na)

        integer :: k


        do k = 1, na
            x_zw(k) =   linear_interpolation_2pt(zt(k:k+1),x_zt(k:k+1),zw(k))
        end do
        
        return
    end subroutine linear_interpolation_zt2zw
	
	
	
	
	
	
	  
	  



end module vmix_canuto
