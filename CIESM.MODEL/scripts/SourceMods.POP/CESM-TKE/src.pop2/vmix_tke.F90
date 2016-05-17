!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module vmix_tke


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

   implicit none
   private
   save

   public :: init_vmix_tke, &
             vmix_coeffs_tke

   logical (log_kind), parameter :: &
        ln_mxl0   = .FALSE., & ! mixing length scale surface value as function of wind stress or not
        ln_lc     = .FALSE.    ! Langmuir cells (LC) as a source term of TKE or not
   integer (int_kind), parameter :: &
                       nn_mxl    =  2, & ! type of mixing length (=0/1/2/3)  
                       nn_pdl    =  1, & ! Prandtl number or not (ratio avt/avm) (=0/1)
                       nn_etau   =  0, & ! type of depth penetration of surface tke (=0/1/2/3) 
                       nn_htau   =  0    ! type of tke profile of penetration (=0/1)
    
   real (r8), parameter ::  &
                       rau0     =   1026.0_r8,  &
                       vkarmn   =    0.4_r8!,    & grav     =   9.80665_r8

   real (r8), parameter ::  &
           rn_mxl0      =   0.04_r8,    &! surface  min value of mixing length (kappa*z_o=0.4*0.1 m)  [m]
           rn_ediff     =   0.1_r8,     &! coefficient for avt: avt=rn_ediff*mxl*sqrt(e)
           rn_ediss     =   0.7_r8,     &! coefficient of the Kolmogoroff dissipation 
           rn_ebb       =   3.75_r8,    &! coefficient of the surface input of tke
           rn_emin      =   0.7071e-6_r8,&! minimum value of tke [m2/s2]
           rn_emin0     =   1.e-4_r8,   &! surface minimum value of tke [m2/s2]
           rn_bshear    =   1.e-20_r8,  &! background shear (>0) currently a numerical threshold (do not change it) 
           rn_efr       =   1.0_r8,     &! fraction of TKE surface value which penetrates in the ocean   
           rhftau_add   =   1.e-3_r8,   &   ! add offset   applied to HF part of taum  (nn_etau=3)
           rhftau_scl   =   1.0_r8          ! scale factor applied to HF part of taum  (nn_etau=3)

    real (r8) :: &
           ri_cri,    &! critic Richardson number (deduced from rn_ediff and rn_ediss values)
           rmxl_min  ! minimum mixing length value (deduced from rn_ediff and rn_emin values)  [m]




   contains
    subroutine init_vmix_tke(TKE)
        implicit none
        real (r8), dimension(:,:,:,:), intent(inout) :: &
            TKE        !turbulent kinetic energy 
        integer (int_kind) ::   &
            n,                  &! local dummy index for tracer
            k,                  &! local dummy index for vertical lvl
            i, j, iblock,       &! local dummy indexes
            nml_error            ! namelist i/o error flag

        character (16), parameter :: &
            fmt_real = '(a30,2x,1pe12.5)'

        character (11), parameter :: &
            fmt_log  = '(a30,2x,l7)'

        character (11), parameter :: &
            fmt_int  = '(a30,2x,i5)'

        character (char_len) :: &
            string, string2


        ri_cri  =   2.0_r8/(2.0_r8+rn_ediss/rn_ediff)
        rmxl_min    =   1.e-6_r8/(rn_ediff*sqrt(rn_emin))

        TKE =   c0

        if (my_task == master_task) then
            write(stdout,*) ' vmix_tke settings:'
            write(stdout,fmt_real) '  rn_mxl0               =', rn_mxl0 
            write(stdout,fmt_real) '  rn_ediff              =', rn_ediff 
            write(stdout,fmt_real) '  rn_ediss              =', rn_ediss 
            write(stdout,fmt_real) '  rn_ebb                =', rn_ebb 
            write(stdout,fmt_real) '  rn_emin               =', rn_emin 
            write(stdout,fmt_real) '  rn_emin0              =', rn_emin0 
            write(stdout,fmt_real) '  rn_bshear             =', rn_bshear
            write(stdout,fmt_real) '  rn_efr                =', rn_efr
            write(stdout,fmt_real) '  rhftau_add            =', rhftau_add
            write(stdout,fmt_real) '  rhftau_scl            =', rhftau_scl
            write(stdout,fmt_real) '  ri_cri                =', ri_cri 
            write(stdout,fmt_real) '  rmxl_min              =', rmxl_min
            write(stdout,fmt_real) '  time_step_tracer      =', dtt
            write(stdout,fmt_real) '  time_step_momentum    =', dtu
            write(stdout,fmt_real) '  time_step_barotrophic =', dtp
        end if
    end subroutine init_vmix_tke
 
    subroutine vmix_coeffs_tke(VDC, VVC, TRCR, UUU, VVV, UCUR, VCUR, RHOMIX, STF, SHF_QSW, &
                            this_block, convect_diff, convect_visc, &
                            SMF, SMFT)
        real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
                    TRCR                ! tracers at current time

        real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
                    UUU, VVV,           &! velocities at mix time
                    UCUR,VCUR           ! velocities at current time
        
        real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
                    RHOMIX              ! density at mix time

        real (r8), dimension(nx_block,ny_block,nt), intent(in) :: &
                    STF                 ! surface forcing for all tracers

        real (r8), dimension(nx_block,ny_block), intent(in) :: &
                    SHF_QSW             ! short-wave forcing

        real (r8), dimension(nx_block,ny_block,2), intent(in), optional :: &
                    SMF,               &! surface momentum forcing at U points
                    SMFT                ! surface momentum forcing at T points
                                        ! *** either one or the other (not
                                        ! *** both) should be passed

        real (r8), intent(in) :: &
                convect_diff,      &! diffusivity to mimic convection
                convect_visc        ! viscosity   to mimic convection

        type (block), intent(in) :: &
                this_block          ! block information for current block

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block,km), intent(inout) ::      &
                VVC        ! viscosity for momentum diffusion

        real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
                VDC        ! diffusivity for tracer diffusion

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

        character (char_len) ::  &
                error_string

        integer (int_kind) :: &
                k,                 &! vertical level index 
                i,j,               &! horizontal loop indices
                n,                 &! tracer index
                mt2,               &! index for separating temp from other trcrs
                bid                 ! local block address for this block

        integer (int_kind), dimension(nx_block,ny_block) :: &
                KBL                   ! index of first lvl below hbl

        real (r8), dimension(nx_block,ny_block) :: &
                USTAR,      &! surface friction velocity
                taum,       &! wind stress
                WORK,       &
                BFSFC,      &! surface buoyancy forcing
                WORK1,WORK2,&! temporary storage
                FCON,       &! convection temporary
                STABLE,     &! = 1 for stable forcing; = 0 for unstable forcing
                KE_mix,     &! kinetic energy at mix time
                KE_cur,     &! kinetic energy at cur time
                En           ! En for boundary layer kinetic energy
        
        real (r8), dimension(nx_block,ny_block,km) :: &
                DBLOC,      &! buoyancy difference between adjacent levels
                DBSFC,      &! buoyancy difference between level and surface
                GHAT         ! non-local mixing coefficient

        real (r8), dimension(nx_block,ny_block,0:km+1) :: &
                VISC        ! local temp for viscosity

        real (r8) ::  &
                factor

        character (16), parameter :: &
            fmt_real = '(a30,2x,1pe12.5)'

        character (11), parameter :: &
            fmt_log  = '(a30,2x,l7)'

        character (11), parameter :: &
            fmt_int  = '(a30,2x,i5)'

        real (r8) :: zbbrau     !local
        real (r8) :: zraug      !local
        real (r8), dimension(nx_block,ny_block) :: tke_local, zmxlm
        real (r8), dimension(nx_block,ny_block) :: zat, zav 
!-----------------------------------------------------------------------
!
!  initialize  and consistency checks
!
!-----------------------------------------------------------------------

        bid = this_block%local_id

        if (.not. present(SMF) .and. .not. present(SMFT)) then
            error_string = 'ERROR tke: must supply either SMF or SMFT'
            call document ('vmix_coeffs_tke',  trim(error_string))
            call exit_POP(sigAbort, trim(error_string))
        endif

        if (present(SMF)) then
            if (my_task == master_task) then
                write(stdout,*) ' vmix_tke settings: smf present'
            end if
        elseif (present(SMFT)) then
            if (my_task == master_task) then
                write(stdout,*) ' vmix_tke settings: smft present'
            end if
        end if

        if (present(SMFT)) then
            taum = 0.1*sqrt(SMFT(:,:,1)**2 + SMFT(:,:,2)**2)
        else
            WORK = 0.1*sqrt(SMF(:,:,1)**2 + SMF(:,:,2)**2)
            call ugrid_to_tgrid(taum,WORK,bid)
        endif

        !uuu =   0.01_r8*uuu
        !vvv =   0.01_r8*vvv
        !ucur=   0.01_r8*ucur
        !vcur=   0.01_r8*vcur
        !rhomix  =   1000.0_r8*rhomix
        !trcr(:,:,:,2)   =   trcr(:,:,:,:2)/ppt_to_salt

        if (my_task == master_task) then
            write(stdout,*) 'my_task:  ', my_task
            write(stdout,*) 'block_id: ', bid, ' max and min of u', maxval(uuu),minval(uuu)
            write(stdout,*) 'block_id: ', bid, ' max and min of v', maxval(vvv),minval(vvv)
            write(stdout,*) 'block_id: ', bid, ' max and min of ucur', maxval(ucur),minval(ucur)
            write(stdout,*) 'block_id: ', bid, ' max and min of vcur', maxval(vcur),minval(vcur)
            write(stdout,*) 'block_id: ', bid, ' max and min of rhomix', maxval(rhomix),minval(rhomix)
            write(stdout,*) 'block_id: ', bid, ' max and min of trcr1', &
                                        maxval(trcr(:,:,:,1)),minval(trcr(:,:,:,1))
            write(stdout,*) 'block_id: ', bid, ' max and min of trcr2', &
                                        maxval(trcr(:,:,:,2)),minval(trcr(:,:,:,2))
            write(stdout,*) 'block_id: ', bid, ' max and min of trcr3', &
                                        maxval(trcr(:,:,:,3)),minval(trcr(:,:,:,3))
            write(stdout,*) 'block_id: ', bid, ' max and min of stf1', &
                                        maxval(stf(:,:,1)),minval(stf(:,:,1))
            write(stdout,*) 'block_id: ', bid, ' max and min of stf2', &
                                        maxval(stf(:,:,2)),minval(stf(:,:,2))
            write(stdout,*) 'block_id: ', bid, ' max and min of stf3', &
                                        maxval(stf(:,:,3)),minval(stf(:,:,3))
            write(stdout,*) 'block_id: ', bid, ' max and min of SHF_QSW', &
                                        maxval(shf_qsw),minval(shf_qsw)
            if (present(smf)) then
                write(stdout,*) 'block_id: ', bid, ' max and min of smf1', &
                                        maxval(smf(:,:,1)),minval(smf(:,:,1))
                write(stdout,*) 'block_id: ', bid, ' max and min of smf2', &
                                        maxval(smf(:,:,2)),minval(smf(:,:,2))
            else
                write(stdout,*) 'block_id: ', bid, ' max and min of smft1', &
                                        maxval(smft(:,:,1)),minval(smft(:,:,1))
                write(stdout,*) 'block_id: ', bid, ' max and min of smft2', &
                                        maxval(smft(:,:,2)),minval(smft(:,:,2))
            end if
            write(stdout,*) 'max and min of taum', maxval(taum), minval(taum) 
        end if


        zbbrau  =   rn_ebb/rau0
        do j=1,ny_block
        do i=1,nx_block
            if (KMT(i,j,bid) == 0) then
                tke_local(i,j)  =   0.0_r8
                zmxlm(i,j)      =   0.0_r8
                zat(i,j)        =   0.0_r8
            else
                tke_local(i,j)  =   max(rn_emin0, zbbrau * taum(i,j)) 
                zraug   =   vkarmn * 2.e5_r8 / ( rau0 * 9.80665_r8 )
                zmxlm(i,j) = max(rn_mxl0, zraug * taum(i,j))
                zat(i,j)    =   10000.0_r8*rn_ediff * zmxlm(i,j) * sqrt(tke_local(i,j)) 
            end if
        end do
        end do

        call tgrid_to_ugrid(zav,zat,bid)

        if (my_task == master_task) then
            write(stdout,*) 'max and min of zat', maxval(zat), minval(zat) 
            write(stdout,*) 'max and min of zav', maxval(zav), minval(zav) 
            write(stdout,*) 'max and min of surf_vvc', maxval(vvc(:,:,1)), minval(vvc(:,:,1)) 
            write(stdout,*) 'max and min of surf_vdc_1', maxval(vdc(:,:,1,1)), minval(vdc(:,:,1,1)) 
            write(stdout,*) 'max and min of surf_vdc_2', maxval(vdc(:,:,1,2)), minval(vdc(:,:,1,2)) 
        end if
        vvc(:,:,1)  =   max(vvc(:,:,1), zav)
        vdc(:,:,1,1)    =   max(vdc(:,:,1,1),zat)
        vdc(:,:,1,2)    =   max(vdc(:,:,1,2),zat)

 

    end subroutine vmix_coeffs_tke


   end module vmix_tke 
