#define pop_g16
#undef licom
#undef pop_s05
module const_for_nasagissmixing
    implicit none
    integer, parameter   :: r8       =  selected_real_kind(12)
    real(r8), parameter  :: fricmx   =  1000.0 
    real(r8), parameter  :: wndmix   =  10.0 
#if (defined licom)
    integer, parameter   :: nmax     =  30
#endif
#if (defined pop_g16)
    integer, parameter   :: nmax     =  60
#endif
#if (defined pop_s05)
    integer, parameter   :: nmax     =  40
#endif
    integer, parameter   :: isurfuse    =  1
    integer, parameter   :: ifextermld  =  0
    integer, parameter   :: ifoutput    =  0

    contains

    !---------------------------------------------------------------------
    !                   Subroutine hwy_density_interface 
    !---------------------------------------------------------------------
    subroutine hwy_density_interface(rho_out,theta_in,s_in,p_in,if_pops_in,if_cm_in,ifgcm3_out)
        
        real(r8),   intent(in)  ::  theta_in, s_in, p_in
        real(r8),   intent(out) ::  rho_out
        integer,    intent(in)  ::  if_pops_in, if_cm_in, ifgcm3_out


        real (r8)   ::  theta, s, p

        theta   =   theta_in
        s       =   s_in
        p       =   p_in

        if (if_pops_in .eq. 1)  then
            s   =   s_in*1.0d+3     ! back to psu
        end if
        if (if_cm_in .eq. 1) then
            p   =   p_in*1.0d-2     ! cm to m
        end if

        call density_mcdougall2003(rho_out,theta,s,p)

        if (ifgcm3_out .eq. 1) then
            rho_out =   rho_out*1.0d-3
        end if

        return
    end subroutine hwy_density_interface
    !---------------------------------------------------------------------
    !               End of Subroutine hwy_density_interface 
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !                   Subroutine hwy_alpha_beta_interface 
    !---------------------------------------------------------------------
    subroutine hwy_alpha_beta_interface(alpha,beta,theta_in,s_in,p_in,if_pops_in,if_cm_in)

        real(r8), intent(in)  ::  theta_in, s_in, p_in
        real(r8), intent(out) ::  alpha,  beta
        integer, intent(in)   ::  if_pops_in, if_cm_in

        real (r8)   ::  theta, s, p

        theta   =   theta_in
        s       =   s_in
        p       =   p_in

        if (if_pops_in .eq. 1)  then
            s   =   s_in*1.0d+3     ! back to psu
        end if
        if (if_cm_in .eq. 1) then
            p   =   p_in*1.0d-2     ! cm to m
        end if

        call alpha_beta(theta,s,p,alpha,beta)
      
        return
    end subroutine hwy_alpha_beta_interface
    !---------------------------------------------------------------------
    !             End of Subroutine hwy_alpha_beta_interface 
    !---------------------------------------------------------------------



    !---------------------------------------------------------------------
    !                   Subroutine density_mcdougall2003
    !---------------------------------------------------------------------
    subroutine density_mcdougall2003(rho_out,theta_in,s_in,p_in)
      real(r8), intent(in)  :: theta_in, s_in, p_in
      real(r8), intent(out) :: rho_out
      real(r8) :: p_dbar, a(0:11), b(0:12), P1, P2
      ! the p_in input must be in m or dbar
      p_dbar = abs(p_in)

      a = (/ +9.99843699d+2, +7.35212840d+0, -5.45928211d-2, &
             +3.98476704d-4, +2.96928239d+0, -7.23268813d-3, &
             +2.12382341d-3, +1.04004591d-2, +1.03970529d-7, &
             +5.18761880d-6, -3.24041825d-8, -1.23869360d-11  /)

      b = (/ +1.0d0,          +7.28606739d-3,  -4.60835542d-5,  &
             +3.68390573d-7,  +1.80809186d-10, +2.14691708d-3,  &
             -9.27062484d-6,  -1.78343643d-10, +4.76534122d-6,  &
             +1.63410736d-9,  +5.30848875d-6,  -3.03175128d-16, &
             -1.27934137d-17  /)
      
      P1 = a(0) + a(1)*theta_in**1 & 
                + a(2)*theta_in**2 &
                + a(3)*theta_in**3 &
                + a(4)*s_in &
                + a(5)*s_in*theta_in &
                + a(6)*s_in**2 &
                + a(7)*p_dbar &
                + a(8)*p_dbar*theta_in**2 &
                + a(9)*p_dbar*s_in &
                + a(10)*p_dbar**2 &
                + a(11)*p_dbar**2*theta_in**2
      

      P2 = b(0) + b(1)*theta_in**1 &
                + b(2)*theta_in**2 &
                + b(3)*theta_in**3 &
                + b(4)*theta_in**4 &
                + b(5)*s_in        &
                + b(6)*s_in*theta_in**1 &
                + b(7)*s_in*theta_in**3 &
                + b(8)*s_in**(3.0d0/2.0d0) &
                + b(9)*s_in**(3.0d0/2.0d0)*theta_in**2 &
                + b(10)*p_dbar &
                + b(11)*p_dbar**2*theta_in**3 &
                + b(12)*p_dbar**3*theta_in

      rho_out = P1/P2
    end subroutine density_mcdougall2003
    !---------------------------------------------------------------------
    !               End of Subroutine density_mcdougall2003
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    !                   Subroutine alpha_beta
    !---------------------------------------------------------------------
    subroutine alpha_beta(theta_in,s_in,p_in,alpha, beta)

      real(r8), intent(in)  :: theta_in, s_in, p_in
      real(r8), intent(out) :: alpha, beta
      real(r8) :: rho_t_1, rho_t_2, rho_s_2
      real(r8), parameter :: dt = 0.01
      real(r8), parameter :: ds = 0.001

      

      call density_mcdougall2003(rho_t_1 ,theta_in,s_in,p_in)
      call density_mcdougall2003(rho_t_2 ,theta_in+dt,s_in,p_in)

!       call density_mcdougall2003(rho_s_1 ,theta_in,s_in,p_in)
      call density_mcdougall2003(rho_s_2 ,theta_in,s_in+ds,p_in)

      alpha   = (rho_t_2 - rho_t_1) / (dt * rho_t_1)
      beta    = (rho_s_2 - rho_t_1) / (ds * rho_t_1)

      alpha   = alpha * (-1.0D0)

    end subroutine alpha_beta
    !---------------------------------------------------------------------
    !                   End of Subroutine alpha_beta
    !---------------------------------------------------------------------


end module const_for_nasagissmixing

