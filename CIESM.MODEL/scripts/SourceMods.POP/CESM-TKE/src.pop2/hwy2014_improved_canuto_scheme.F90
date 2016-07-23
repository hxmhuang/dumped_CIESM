module hwy2014_improved_canuto_scheme
! Purpose:
! This module is devised for ...
!
!
! Record of revisions:
!    Date             Programmer                Description of change  
!==============    ===================   =============================
! 2016/4/10          WenYu Huang                Original Code
! ...
!
!Any comments please send to huangwenyu@mail.tsinghua.edu.cn                      
    !use module_precision
    !use parameters
    !use utils
    implicit none
    private
    public canuto_2010_main

    integer, parameter :: r8   =   8
    integer, parameter :: charlen  =   500
    real(r8), parameter :: Re = 6.371d6 ! earth radius
    real(r8), parameter :: pi = 4*atan(1.0d0)
    real(r8), parameter :: torad = pi/180.0d0
    real(r8), parameter :: todeg = 180.0d0/pi
    real(r8), parameter :: Cp = 3901.0d0
    real(r8), parameter :: abs_error    =   0.001d0
    real(r8), parameter :: very_small   =   1.d-15
    real(r8), parameter :: d0   =   1.26d+3
    real(r8), parameter :: od0  =   1.0d0/d0
    real(r8), parameter :: g0   =   9.8d0
    real(r8), parameter :: karman = 0.4d0
    real(r8), parameter :: Omega  = 7.29d-5

    contains
    subroutine shengjin_calc(array,nmax,real_x)
    implicit none
    integer,  intent(in)  :: nmax
    real(r8), intent(in)  :: array(nmax) 
    real(r8), intent(out) :: real_x
    !real(r8), intent(out) :: out0
    real(r8) :: sa, sb, sc, sd 
    real(r8) :: A, B, C, T, the
    real(r8) :: delta, x(3), y(2)
    integer :: k, counter 

    sa = array(1)
    sb = array(2)
    sc = array(3)
    sd = array(4)

    if (abs(sa) <= very_small .and. abs(sb) <= very_small) then
            if ((-sd/sc)<0.0d0) then
                    print *,  "No real positive root at position 1"
                    return
            else
                    print *, "Only one real positive root exists 1!"
                    real_x = -sd/sc 
                    return
            end if
    end if

    if (abs(sa) <= very_small .and. abs(sb) >= very_small) then
            delta = sc**2 - 4.0d0*sb*sd
            if (delta < 0.0d0) then
                    print *, "No real roots in this secondary order equation 2!"
                    return
            else if (delta >= 0.0d0) then
                    X(1) = (-sc + sqrt(delta))/(2.0d0*sb)
                    X(2) = (-sc - sqrt(delta))/(2.0d0*sb)
                    if (X(1) < 0.0d0 .and. X(2) < 0.0d0) then
                            print *,  "No real positive root at position 2"
                            return
                    else 
                            print *,  "Two real with one positive at least 2!"
                            real_x = minval(x(1:2),mask = x(1:2) > 0.0d0)
                            return 
                    end if
            end if
    end if

    A   =   sb**2-3.0d0*sa*sc
    B   =   sb*sc-9.0d0*sa*sd
    C   =   sc**2-3.0d0*sb*sd 

    if (abs(A) <= very_small .and. abs(B) <= very_small) then
            if ((-sc/sb) < 0.0d0) then
                    print *,  "No real positive root at position 3"
                    return
            else
                    print *, "Only one real positive root exists 3!"
                    real_x = -sc/sb 
                    return
            end if
    end if

    delta   =   B**2-4.0d0*A*C

    if (delta > 0.0d0) then
            Y(1) = A*sb+3.0d0*sa*(-B+sqrt(delta))/2.0d0
            Y(2) = A*sb+3.0d0*sa*(-B-sqrt(delta))/2.0d0
            real_x = (-sb-sign(abs(y(1))**(1.0d0/3.0d0), y(1))-sign(abs(y(2))**(1.0d0/3.0d0), y(2)))/(3.0d0*sa)
            if ((real_x) < 0.0d0) then
                    print *,  "No real positive root at position 4"
                    return
            else
                    print *, "Only one real positive root exists 4!"
                    print *,  sign(abs(y(1))**(1.0d0/3.0d0), y(1))
                    print *, y(1), y(2),sb, sa,real_x
                    real_x = real_x 
                    return
            end if
    end if


    T   =   (2.d0*A*sb-3.d0*sa*B)/(2.0d0*A**(1.5))
    the =   acos(T)
    x(1) = (-sb-2.0d0*sqrt(A)*cos(the/3.0d0))/3.0d0/sa
    x(2) = (-sb+sqrt(A)*(cos(the/3.0d0)+sqrt(3.0d0)*sin(the/3.0d0)))/3.0d0/sa
    x(3) = (-sb+sqrt(A)*(cos(the/3.0d0)-sqrt(3.0d0)*sin(the/3.0d0)))/3.0d0/sa
    real_x = 500.0d0
    counter = 0
    do k = 1, 3
            if (x(k) >= 0.0d0 .and. x(k) <= 1.d3) then
                    if (x(k) < real_x) then
                            real_x = x(k)
                    end if
                    counter = counter + 1
            end if
    end do

    real_x = minval(x(:),mask = x(:) > 0.0d0)
    return 
    end subroutine shengjin_calc 

    subroutine best_initial_guess(array,nmax,out0)
        implicit none
        integer,  intent(in)  :: nmax
        real(r8), intent(in)  :: array(nmax) 
        real(r8), intent(out) :: out0

        integer,  parameter   :: num = 100000
        real(r8) :: guess_value(-num:num)
        real(r8) :: res_value(-num:num)
        real(r8) :: x, ref_value
        integer :: j, k
        
        ref_value   =   1.d8 
        do k = -num, num, 1
            guess_value(k) = dble(k)*0.1d0
            x   =   guess_value(k)
            res_value(k) = 0.0d0
            do j = 1, nmax,1
                res_value(k) = res_value(k) + array(j)*x**(nmax-j)
            end do
        end do

        do k = -num, num
            if (abs(res_value(k)) <= ref_value .and. guess_value(k) >= 0.0d0) then
                if (k>-num .and. k<num) then
                    if (res_value(k-1)*res_value(k+1)<0.0d0) then  
                        out0 = guess_value(k)
                        ref_value = abs(res_value(k))
                    end if
                end if
            end if
        end do
    end subroutine best_initial_guess

    subroutine newton_iteration(array,nmax,out0)
        implicit none
        integer,  intent(in)  :: nmax
        real(r8), intent(in)  :: array(nmax)
        real(r8), intent(out) :: out0

        integer,  parameter   :: max_iter    =   1000
        real(r8)              :: start0      
        real(r8), parameter   :: eps         =   1.0d-8 

        integer     :: counter, k
        real(r8)    :: x, fx, dfx
        logical     :: if_converge
        
        call best_initial_guess(array(:),nmax,start0)
        call shengjin_calc(array(:),nmax,x)
        out0 = x
        return
        !start0  =   10.0d0
1000    continue
        if_converge = .false.
        counter = 0
        x = start0
        
        do while(.not. if_converge .and. counter < max_iter)
            fx  =  0.0d0
            dfx =  0.0d0
            do k = 1, nmax
                fx  =   fx  +  array(k)*x**(nmax-k)  
            end do
            do k = 1, nmax -1
                dfx =   dfx +  dble(nmax-k)*array(k)*x**(nmax-1-k)
            end do
            x   =      x    -   fx/dfx
            counter = counter + 1
            if_converge = abs(fx) <= eps
        end do

        out0 = x

        if (out0 < 0)  then
            write(*,*) "Strong Error: Gm < 0"
        end if
        if (out0 < 0)  then
            start0 = start0 + 1000.0d0
        !    go to 1000
        end if
        if (.not. if_converge) write(*,*) "after all the iteration times, the program not converged"
    end subroutine

    Subroutine prepare_pi(Ri,Rrho,out_pi)
        implicit none
        real(r8), intent(in)  :: Ri, Rrho
        real(r8), intent(out) :: out_pi(5)
        
        real(r8), parameter   :: a  =   10.0d0    
        real(r8), parameter   :: Ko =   1.66d0 
        real(r8), parameter   :: sig_t  =    0.72d0
        real(r8), parameter   :: zero   =    1.d-10

        real(r8), parameter   :: pi0_1   =    (27.0d0/5.0d0*Ko**3)**(-0.5d0)*(1+sig_t**(-1))**(-1)
        real(r8), parameter   :: pi0_2   =    1.0d0/3.0d0 
        real(r8), parameter   :: pi0_3   =    sig_t 
        real(r8), parameter   :: pi0_4   =    pi0_1 
        real(r8), parameter   :: pi0_5   =    sig_t 

        if (Ri > zero .and. Rrho > zero) then
            out_pi(1) = pi0_1*(1.0d0 + Ri*Rrho/(a + Rrho))**(-1)
            out_pi(2) = pi0_2*((1.0d0 + Ri)**(-1))*(1.0d0+2.0d0*Ri*Rrho*(1.0d0+Rrho**2)**(-1)) 
            out_pi(3) = pi0_3
            out_pi(4) = pi0_4*(1.0d0 + Ri/(1.0d0 + a*Rrho))**(-1)
            out_pi(5) = pi0_5
        end if

        if (Ri > zero .and. Rrho <= zero) then
            out_pi(1) = pi0_1*(1.0d0 + Ri)**(-1)
            out_pi(2) = pi0_2
            out_pi(3) = pi0_3
            out_pi(4) = pi0_4*(1.0d0 + Ri)**(-1)
            out_pi(5) = pi0_5
        end if

        if (Ri <= zero) then
            out_pi(1) = pi0_1
            out_pi(2) = pi0_2
            out_pi(3) = pi0_3
            out_pi(4) = pi0_4
            out_pi(5) = pi0_5
        end if

    End Subroutine prepare_pi

    Subroutine dyn_time_scale_calc(Ri,Rrho,pi_n,Gm)
        implicit none
        real(r8), intent(in)  :: Ri, Rrho, pi_n(5)
        real(r8), intent(out) :: Gm
        real(r8) :: A1, A2, A3, A4, A5, A6
        real(r8) :: B1, B2, B3, B4, B5, B6
        real(r8) :: cubic_coef(4)
        real(r8) :: temp1,  temp2,  temp3,  temp4,  temp5
        real(r8) :: temp6,  temp7,  temp8,  temp9,  temp10
        real(r8) :: temp11, temp12, temp13, temp14, temp15

    !    call prepare_pi(Ri,Rrho,pi_n(:))

    !    if (abs(1.0d0-Rrho) <= very_small) then
    !        Gm = 60.0d0
    !        return
    !    end if

        temp1   =   pi_n(1)*pi_n(4)*(pi_n(4)-pi_n(1)*Rrho)
        temp2   =   pi_n(2)*(15.0d0*pi_n(3)+7.0d0)*(Rrho**2+1.0d0)
        temp3   =   (14.0d0*(pi_n(2)-pi_n(3))-15.0d0*pi_n(3)**2)*Rrho
        temp4   =   150.0d0*(1.0d0-Rrho)**3 !hwy_need_check

        A1  =   temp1*(temp2+temp3)/temp4
        B1  =   pi_n(1)*(pi_n(4)**2)*(temp2+temp3)/150.0d0 ! A1*temp4/(1-Rrho)
        
        temp1   =   pi_n(1)*pi_n(4)
        temp2   =   pi_n(2)*(210.0d0*pi_n(1)-150.0d0*pi_n(3) + 7.0d0)*(Rrho**2+1.0d0)
        temp3   =   14.0d0*(pi_n(2)-pi_n(3))*(1.0d0+15.0d0*pi_n(1)+15.0d0*pi_n(4))
        temp4   =   150.0d0*pi_n(3)**2
        temp5   =   (temp3+temp4)*Rrho
        temp6   =   210.0d0*pi_n(2)*(pi_n(4)-pi_n(1))

        temp7   =   9000.0d0*(1.0d0-Rrho)**2 !hwy_need_check

        A2  =   temp1*(temp2+temp5+temp6)/temp7
        B2  =   temp1*(temp2+temp5+temp6)/9000.0d0

        temp1   =   pi_n(1)
        temp2   =   5.0d0*pi_n(2)*pi_n(4)*(30.0d0*pi_n(3)+17.0d0)  
        temp3   =   pi_n(1)*(15.0d0*pi_n(3)+7.0d0)
        temp4   =   Rrho**2+1.0d0
        temp5   =   temp1*(temp2+temp3)*temp4
        temp6   =   -(15.0d0*pi_n(3)+7.0d0)*(pi_n(1)**2-pi_n(4)**2)
        temp7   =   10.0d0*pi_n(1)*pi_n(3)*pi_n(4)*(15.0d0*pi_n(3)+17.0d0)
        temp8   =   15.0d0*pi_n(2)*(pi_n(1)**2+pi_n(4)**2)
        temp9   =   14.0d0*pi_n(1)*pi_n(4)*(1.0d0-10.0d0*pi_n(2))
        temp10  =   -(temp7+temp8+temp9)*Rrho
        temp11  =   150.0d0*(1.0d0-Rrho)**2 !hwy_need_check
        
        A3  =   (temp5 + temp6 + temp10)/temp11
        B3  =   (temp5 + temp6 + temp10)/150.0d0
        
        temp1   =   150.0d0*(pi_n(1)*pi_n(3)+pi_n(2)*pi_n(4))
        temp2   =   -7.0d0*pi_n(1)*(1.0d0+30.0*pi_n(1))
        temp3   =   (temp1+temp2)*Rrho
        temp4   =   -150.0d0*(pi_n(1)*pi_n(2)+pi_n(3)*pi_n(4))
        temp5   =   7.0d0*pi_n(4)*(1.0d0+30.0d0*pi_n(4))
        temp6   =   9000.0d0*(1.0d0-Rrho) !hwy_need_check
        A4  =   (temp3 + temp4 + temp5)/temp6
        B4  =   (temp1 + temp2)*(-1) 

        temp1   =   -30.0d0*(pi_n(1)*pi_n(3)+pi_n(2)*pi_n(4))
        temp2   =   -17.0d0*pi_n(1)
        temp3   =   (temp1+temp2)*Rrho
        temp4   =   30.0d0*(pi_n(1)*pi_n(2)+pi_n(3)*pi_n(4))
        temp5   =   17.0d0*pi_n(4)
        temp6   =   30.0d0*(1.0d0-Rrho) ! hwy_need_check 
        A5  =   (temp3 + temp4 + temp5)/temp6
        B5  =   (temp3 + temp4 + temp5)

        A6  =   -1.0d0/60.0d0
        B6  =   A6
        
        cubic_coef(1)   =   A1*Ri**3 + A2*Ri**2
        cubic_coef(2)   =   A3*Ri**2 + A4*Ri
        cubic_coef(3)   =   A5*Ri    + A6
        cubic_coef(4)   =   1.0d0
        
        !call newton_iteration(cubic_coef,4,Gm)
        call shengjin_calc(cubic_coef,4,Gm)
        !call check_single_value_real_r8_1d("A",[B1,B2,B3,B4,B5,B6],6)
    
    End Subroutine dyn_time_scale_calc

    Subroutine struct_function(Ri,Rrho,pi_n,Gm,struct_m,struct_h,struct_s,struct_rho,r_out,r_out_new)
        implicit none
        real(r8), intent(in)  :: Ri, Rrho
        real(r8), intent(in)  :: pi_n(5)
        real(r8), intent(in)  :: Gm
        real(r8), intent(out) :: struct_m, struct_h, struct_s,struct_rho, r_out, r_out_new

        real(r8) :: x, p, q, r, ratio, Ah, As, Am, LX
        real(r8) :: Am1, Am2
        real(r8) :: r_new, Am1_new, Am2_new, Am_new, Ah_new, As_new, LX_new, ratio_new
        real(r8) :: ratio_m, ratio_s, ratio_h, ratio_rho

        x   =   Ri*Gm*(1.0d0-Rrho)**(-1) ! hwy_need_check
        p   =   pi_n(4)*pi_n(5) - pi_n(4)*pi_n(2)*(1.0d0+Rrho)
        q   =   pi_n(1)*pi_n(2)*(1.0d0+Rrho) - pi_n(1)*pi_n(3)*Rrho
        
        r   =   (pi_n(4)/pi_n(1))/Rrho*(1.0d0+q*x)/(1.0d0+p*x)     !hwy_need_check
        r_new   =   (pi_n(4)/pi_n(1))*(1.0d0+q*x)/(1.0d0+p*x)     !hwy_need_check
        r_out = r
        r_out_new = r_new

        Ah  =   pi_n(4)/(1.0d0+p*x+pi_n(4)*pi_n(2)*x*(1.0d0-1.0d0/r))
        As  =   Ah/r_new


        Am1  =    4.0d0/5.0d0 - (pi_n(4)-pi_n(1)+(pi_n(1)-1.0d0/150.0d0)*(1.0d0 -1.0d0/r))*x*Ah
        Am2  =   10.0d0 + (pi_n(4)-pi_n(1)*Rrho)*x + 1.0d0/50.0d0*Gm
        Am   =   Am1/Am2
        
        LX  =   (1.0d0 - 1.0d0/r)*x*Ah

        ratio   =  2.0d0/3.0d0/(1.0d0+2.0d0/15.0d0*LX + 1.0d0/10.0d0*Am*Gm)

        struct_m    =   ratio*Am1/Am2
        struct_h    =   ratio*Ah
        struct_s    =   ratio*As
        struct_rho  =   (struct_h-struct_s*Rrho)/(1.0d0-Rrho)

    End Subroutine struct_function
    !---------------------------------------------------------------------
    !                      Subroutine mixing_efficiency
    !---------------------------------------------------------------------
    subroutine mixing_efficiency(mix_eff_var,struct_var,Ri,Gm)
        real(r8), intent(in) :: Ri, Gm
        real(r8), intent(in) :: struct_var
        real(r8), intent(out) :: mix_eff_var 
        mix_eff_var = 0.5d0*Ri*Gm*struct_var
        return
    end subroutine mixing_efficiency
    !---------------------------------------------------------------------
    !                 End of Subroutine mixing_efficiency
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                       Subroutine find_mix_layer_depth
    !---------------------------------------------------------------------
    subroutine find_mix_layer_depth(mld_out,mld_lev_out,den,delta,zlev_in,n)
        integer, intent(in) :: n ! not the whole may be the whole
        integer, intent(out) :: mld_lev_out
        real(r8), intent(in) :: delta
        real(r8), intent(in) :: zlev_in(n), den(n)
        real(r8), intent(out) :: mld_out

        integer :: k
        real(r8) :: den_m, zlev(n)

        zlev(:) = abs(zlev_in(:))

        do k = 1, n
            if (abs(den(k)-den(1)) > delta) then
                den_m = den(1) - sign(delta,den(1)-den(k))
                mld_out = zlev(k) + (zlev(k-1)-zlev(k))*(den_m-den(k))/(den(k-1)-den(k)+very_small)                       
                mld_lev_out = k
                exit
            end if
            mld_out = zlev(n)
            mld_lev_out = n
        end do
    
    end subroutine find_mix_layer_depth
    !---------------------------------------------------------------------
    !                   End of Subroutine find_mix_layer_depth
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                       Subroutine canuto_2010_main
    !---------------------------------------------------------------------
    subroutine canuto_2010_main(Ri_in, Rrho_in, &
                                mld_in, lev_in, N2_in, shear2_in, &
                                Gm,struct_m,struct_h,struct_s,struct_rho,&
                                mix_eff_m,mix_eff_h,mix_eff_s,mix_eff_rho,&
                                Km_out,Kh_out,Ks_out,Kd_out)
        real(r8), intent(in) :: Ri_in, Rrho_in
        real(r8), intent(out) :: struct_m,struct_h,struct_s, struct_rho
        real(r8), intent(out) :: mix_eff_m, mix_eff_h, mix_eff_s, mix_eff_rho
        real(r8), intent(out) :: Gm
        real(r8), intent(in) :: mld_in, lev_in, N2_in, shear2_in
        real(r8), intent(out) :: Km_out,Kh_out,Ks_out,Kd_out 
        real(r8) :: out_pi(5), Gm0, r, rnew, rf, tke_mld
        
        real(r8) :: struct_m_inf,struct_h_inf,struct_s_inf, struct_rho_inf
        real(r8) :: mix_eff_m_inf, mix_eff_h_inf, mix_eff_s_inf, mix_eff_rho_inf
        real(r8) :: Gm_inf
        real(r8) :: out_pi_inf(5), Gm0_inf, r_inf, rnew_inf, rf_inf


        real(r8), parameter :: Rrho_bound = 1.d-3 
        real(r8), parameter :: Ri_low   =  1.0d-2 
        real(r8), parameter :: Ri_high   = 1.0d+3
        real(r8) :: Ri, Rrho

        Ri = Ri_in
        Rrho = Rrho_in

        !Rrho = max(Rrho_in,Rrho_low)
        !Rrho  = Rrho_in
        if (Ri_in > Ri_high) then
            Ri = Ri_high
        else if (Ri_in < Ri_low) then
            Ri = Ri_low
        else
            Ri = Ri_in
        end if

        if (Rrho_in > 0.96d0) then
            Rrho    =   0.96d0
        else if (Rrho_in < 0.0d0) then
            Rrho    =   0.0d0
        else
            Rrho    =   Rrho_in
        end if



        call prepare_pi(1.0d10,Rrho,out_pi_inf(:))
        call dyn_time_scale_calc(1.0d10,Rrho,out_pi_inf(:),Gm_inf)
        call struct_function(1.0d10,Rrho,out_pi_inf(:),Gm_inf,&
                            struct_m_inf,struct_h_inf,struct_s_inf,struct_rho_inf,&
                            R_inf,Rnew_inf)
        call Rf_calc(Rf_inf,1.0d10,Rrho,struct_h_inf,struct_m_inf,Rnew_inf)

        call prepare_pi(Ri,Rrho,out_pi(:))
        call dyn_time_scale_calc(Ri,Rrho,out_pi(:),Gm)
        call struct_function(Ri,Rrho,out_pi(:),Gm,struct_m,struct_h,struct_s,struct_rho,r,rnew)
        call Rf_calc(Rf,Ri,Rrho,struct_h,struct_m,Rnew)
        call mixing_efficiency(mix_eff_m,struct_m,Ri,Gm) 
        call mixing_efficiency(mix_eff_h,struct_h,Ri,Gm) 
        call mixing_efficiency(mix_eff_s,struct_s,Ri,Gm) 
        call mixing_efficiency(mix_eff_rho,struct_rho,Ri,Gm) 

        if ( lev_in <= mld_in) then
            call mixed_layer_TKE_calc(TKE_mld,mld_in,Gm,shear2_in,lev_in,Rf,Rf_inf)
            Km_out = min(mix_eff_m*TKE_mld/(N2_in+very_small),9999.0d0)
            Kh_out = min(mix_eff_h*TKE_mld/(N2_in+very_small),9999.0d0)
            Ks_out = min(mix_eff_s*TKE_mld/(N2_in+very_small),9999.0d0)
            Kd_out = min(mix_eff_rho*TKE_mld/(N2_in+very_small),9999.0d0)
            !print *, "tke_mld: ", tke_mld, "mld: ", mld_in, "Gm: ", gm, "shear2: ",shear2_in, "Km: ",Km_out, "Ks: ",ks_out
        end if
        if (km_out .ge. 9998.0d0 .or. ks_out .ge. 9998.0d0 &
            .or. kh_out .ge. 9998.0d0 .or. kd_out .ge. 9998.0d0) then
            print *, "km: ",km_out, "ks: ",ks_out,&
                     "kh: ",kh_out, "kd: ",kd_out,&
                     "N2: ",N2_in,  "shear2_in: ",shear2_in,&
                     "mld: ",mld_in, "lev: ",lev_in,&
                     "Ri: ",Ri, "Rrho: ",rrho, " tke: ",tke_mld

        end if
        
    end subroutine canuto_2010_main
    !---------------------------------------------------------------------
    !                     End of Subroutine canuto_2010_main
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                          Subroutine mixed_layer_TKE_calc
    !---------------------------------------------------------------------
    subroutine mixed_layer_TKE_calc(TKE_out,mld_in,Gm_in,s2_in,lev_in,Rf_in,Rf_inf_in)
        real(r8), intent(in) :: mld_in, Gm_in, s2_in, lev_in,Rf_in,Rf_inf_in
        real(r8), intent(out) :: TKE_out

        real(r8) :: l0, B1, lB, l
        

        l0 = 0.17d0*mld_in
        B1 = 21.6d0
        lB = karman*abs(lev_in*l0)/(l0+karman*abs(lev_in)) 
        l  = lB*((1.0d0-Rf_in/Rf_inf_in)**4)**(1.0d0/3.0d0)

        TKE_out = (B1**2)*((Gm_in)**(-3.0d0/2.0d0))*(l**2)*(s2_in**(3.0d0/2.0d0))
    
    end subroutine mixed_layer_TKE_calc
    !---------------------------------------------------------------------
    !                      End of Subroutine mixed_layer_TKE_calc
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                Subroutine Thermocline_mixing_coeff_calc
    !---------------------------------------------------------------------
    subroutine Thermocline_mixing_coeff_calc(Km_out, Kh_out, Ks_out, Kd_out,&
                                             mix_eff_m,mix_eff_h,mix_eff_s,mix_eff_d,&
                                             n2_in,&
                                             lat_in)
        real(r8), intent(in) :: mix_eff_m,mix_eff_h,mix_eff_s,mix_eff_d
        real(r8), intent(in) :: lat_in, n2_in
        real(r8), intent(out) :: Km_out, Kh_out, Ks_out, Kd_out

        real(r8) :: f_lat, L_lat
        real(r8), parameter :: N0 = 5.24d-3
        real(r8), parameter :: f30 = 2.0d0*Omega*sin(30.0d0*torad)
        real(r8), parameter :: Tke_n2 = 0.288d-4
        
        real(r8) :: cond1, cond2
        
        f_lat = abs(2.0d0*Omega*sin(lat_in*torad)+very_small)
        cond1 = max(sqrt(N2_in)/f_lat,1.0d0)

        L_lat = (f_lat*acosh(cond1))/&
                (f30*acosh(N0/(f30)))
        
        Km_out = L_lat*Tke_n2*mix_eff_m
        Kh_out = L_lat*Tke_n2*mix_eff_h
        Ks_out = L_lat*Tke_n2*mix_eff_s
        Kd_out = L_lat*Tke_n2*mix_eff_d
        return
    end subroutine Thermocline_mixing_coeff_calc
    !---------------------------------------------------------------------
    !            End of Subroutine Thermocline_mixing_coeff_calc
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                          Subroutine Rf_calc
    !---------------------------------------------------------------------
    subroutine Rf_calc(Rf_out,Ri_in,Rrho_in,struct_h,struct_m,r_in)
        real(r8), intent(in) :: Ri_in, Rrho_in, struct_h, struct_m, r_in
        real(r8), intent(out) :: Rf_out
        Rf_out = Ri_in*struct_h/struct_m*(1.0d0-rrho_in/r_in)!/(1.0d0-Rrho_in)
    end subroutine Rf_calc
    !---------------------------------------------------------------------
    !                      End of Subroutine Rf_calc
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                  Subroutine canuto_2010_interface
    !---------------------------------------------------------------------
    subroutine canuto_2010_interface(Km_out,    Kh_out,     Ks_out,     Kd_out,     mld_out,&
                                     ts_in,     ss_in,      rho_in,     ri_in,      rrho_in,&
                                     n2_in,     s2_in,      lat_in,     lev_in,     num_lev)

        integer, intent(in) :: num_lev
        real(r8), intent(in) :: lev_in(num_lev), s2_in(num_lev-1), n2_in(num_lev-1), lat_in
        real(r8), intent(in) :: ts_in(num_lev), ss_in(num_lev), rho_in(num_lev), ri_in(num_lev-1), rrho_in(num_lev-1)
        real(r8), intent(out) :: Km_out(num_lev), Kh_out(num_lev), Ks_out(num_lev), Kd_out(num_lev), mld_out
        real(r8) :: struct_m,struct_h,struct_s, struct_rho
        real(r8) :: struct_m_inf,struct_h_inf,struct_s_inf, struct_rho_inf
        real(r8) :: mix_eff_m, mix_eff_h, mix_eff_s, mix_eff_rho
        real(r8) :: Rf, Rf_inf, R, Rnew, R_inf, Rnew_inf, TKE_mld
        real(r8), parameter :: Rrho_bound = 1.d-3 
        real(r8), parameter :: Ri_low   = -1.d+10
        real(r8), parameter :: Ri_high   = 1.d+10
        real(r8) :: Ri(num_lev-1), Rrho(num_lev-1), Gm(num_lev-1), out_pi(5), Gm0
        integer :: mld_lev, k

        call find_mix_layer_depth(mld_out,mld_lev,rho_in,0.03d0,lev_in,num_lev)
    

       
        Ri(:) = Ri_in(:)
        Rrho(:) = Rrho_in(:)

        where (Ri_in > Ri_high) 
            Ri = Ri_high
        elsewhere (Ri_in < Ri_low)
            Ri = Ri_low
        endwhere

        do k = 1, num_lev - 1
            if (abs(Rrho_in(k)-1.0d0) < Rrho_bound) then
                if (Rrho_in(k) >= 1.0d0) then
                    Rrho(k)  =  1.001d0
                else if (Rrho_in(k) <= 1.0d0) then
                    Rrho(k)  =  0.999d0 
                end if
            end if
        end do

        do k = 1, num_lev - 1
            call prepare_pi(1.0d10,Rrho(k),out_pi(:))
            call dyn_time_scale_calc(1.0d10,Rrho(k),out_pi(:),Gm(k))
            call struct_function(1.0d10,Rrho(k),out_pi(:),Gm(k),struct_m_inf,struct_h_inf,struct_s_inf,&
                                                                struct_rho_inf,R_inf,Rnew_inf)
            call Rf_calc(Rf_inf,1.0d10,Rrho(k),struct_h_inf,struct_m_inf,Rnew_inf) 
            call prepare_pi(Ri(k),Rrho(k),out_pi(:))
            call dyn_time_scale_calc(Ri(k),Rrho(k),out_pi(:),Gm(k))
            call struct_function(Ri(k),Rrho(k),out_pi(:),Gm(k),struct_m,struct_h,struct_s,struct_rho,R,Rnew)
            call mixing_efficiency(mix_eff_m,struct_m,Ri(k),Gm(k)) 
            call mixing_efficiency(mix_eff_h,struct_h,Ri(k),Gm(k)) 
            call mixing_efficiency(mix_eff_s,struct_s,Ri(k),Gm(k)) 
            call mixing_efficiency(mix_eff_rho,struct_rho,Ri(k),Gm(k)) 
            call Rf_calc(Rf,Ri(k),Rrho(k),struct_h,struct_m,Rnew)
            
            if ( k <= mld_lev) then
                call mixed_layer_TKE_calc(TKE_mld,mld_out,Gm(k),s2_in(k),lev_in(k),Rf,Rf_inf)
                Km_out(k) = mix_eff_m*TKE_mld/(N2_in(k)+very_small)
                Kh_out(k) = mix_eff_h*TKE_mld/(N2_in(k)+very_small)
                Ks_out(k) = mix_eff_s*TKE_mld/(N2_in(k)+very_small)
                Kd_out(k) = mix_eff_rho*TKE_mld/(N2_in(k)+very_small)
            end if


            if (k > mld_lev) then
                call Thermocline_mixing_coeff_calc(Km_out(k), Kh_out(k), Ks_out(k), Kd_out(k),&
                                                  mix_eff_m,mix_eff_h,mix_eff_s,mix_eff_rho,&
                                                  N2_in(k),&
                                                  lat_in)
            end if
        print*,"!======================================================" 
        !call check_single_value_real_r8("GM",gm(k))
        !call check_single_value_real_r8("Rrho",Rrho(k))
        !call check_single_value_real_r8("Ri",Ri(k))
        !call check_single_value_real_r8("struct_m",struct_m)
        !call check_single_value_real_r8("struct_s",struct_s)
        !call check_single_value_real_r8("struct_h",struct_h)
        !call check_single_value_real_r8("struct_rho",struct_rho)
        !call check_single_value_real_r8("mix_eff_m",mix_eff_m)
        !call check_single_value_real_r8("mix_eff_s",mix_eff_s)
        !call check_single_value_real_r8("mix_eff_h",mix_eff_h)
        !call check_single_value_real_r8("mix_eff_rho",mix_eff_rho)
        !call check_single_value_real_r8("R_inf",R_inf)
        !call check_single_value_real_r8("Rf_inf",Rf_inf)
        !call check_single_value_real_r8("Rnew_inf",Rnew_inf)
        !call check_single_value_real_r8("R",R)
        !call check_single_value_real_r8("Rf",Rf)
        !call check_single_value_real_r8("Rnew",Rnew)
        !call check_single_value_real_r8("TKE_mld",TKE_mld)
        !call check_single_value_real_r8("Km",Km_out(k))
        !call check_single_value_real_r8("Kh",Kh_out(k))
        !call check_single_value_real_r8("Ks",Ks_out(k))
        !call check_single_value_real_r8("Kd",Kd_out(k))
        !print*,"!======================================================" 
        end do
        
    end subroutine canuto_2010_interface
    !---------------------------------------------------------------------
    !               End of Subroutine canuto_2010_interface
    !---------------------------------------------------------------------
     

end module hwy2014_improved_canuto_scheme
