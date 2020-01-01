module random_number_generator
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to the random number generator
!
! Author: Yong Wang, 9/6/2015
!---------------------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8, i8 => shr_kind_i8
    use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
    use abortutils,   only: endrun
    use cam_logfile,  only : iulog


    implicit none
    public
    save

    ! Public methods
    public :: random_number_init ! set random numbers for stochastic deep convection scheme
    public :: random_number_set  ! set random numbers for stochastic deep convection scheme

    real(r8), allocatable, dimension(:,:,:) :: rannumbers          !(pcols,500,begchunk:endchunk)

    real(r8), allocatable, dimension(:,:,:)   :: rannumbers1       !(pcols,100,begchunk:endchunk)
    real(r8), allocatable, dimension(:,:,:)   :: rannumbers_tot1   !(plon,100,plat)
    real(r8), allocatable, dimension(:,:,:,:) :: rannumbers_tot    !(plon,100,plat,500)

contains
    
   
    subroutine random_number_init()

        use dyn_grid,     only : get_dyn_grid_parm

        implicit none

        integer :: plon
        integer :: plat
        integer :: istat

        plon = get_dyn_grid_parm('plon')
        plat = get_dyn_grid_parm('plat')

        !Initialize random numbers
        !if first call additional random numbers for random cloud lifetime
        !Allocate temporary random field on whole domain
        allocate(rannumbers_tot1(plon, 100, plat), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'random_number_init: failed to allocate rannumbers_tot1; error = ',istat
            call endrun
        end if

        allocate(rannumbers_tot(plon, 100, plat, 5), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'random_number_init: failed to allocate rannumbers_tot; error = ',istat
            call endrun
        end if

        !Allocate temporary random field on sub-domain (chunks)
        allocate(rannumbers1(pcols, 100, begchunk:endchunk), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'random_number_init: failed to allocate rannumbers1; error = ',istat
            call endrun
        end if

        allocate(rannumbers(pcols, 500, begchunk:endchunk), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'random_number_init: failed to allocate rannumbers; error = ',istat
            call endrun
        end if

    end subroutine random_number_init

    subroutine random_number_set()
    
    !----------------------------------------------------------------------------------------
    ! Purpose:  set random numbers
    !----------------------------------------------------------------------------------------

        use spmd_utils,     only: masterproc
        use phys_grid, only: scatter_field_to_chunk
        use dyn_grid,     only : get_dyn_grid_parm

        implicit none


        !=================================================================
        !some switches 
        !logical :: zm_stc = .false.
        !==================================================================

        !================================================================== 
        !namelist variables, which will move to other subroutine latter
        integer(i8) :: RANSEED = 0
        !==================================================================
        
        !==================================================================
        !tunalbe parameters
        integer :: MAXRANNUM = 500                       ! maximum random numbers for each column

        !local variables
        integer(i8) :: IRAN
        
        integer :: i, j, k, n

        integer :: plon
        integer :: plat
        
        plon = get_dyn_grid_parm('plon')
        plat = get_dyn_grid_parm('plat')
 
        if (masterproc) then
            CALL SYSTEM_CLOCK(RANSEED)
            RANSEED=mod(RANSEED,2147483646_i8)+1 !SEED between 1 and 2^31-2
            IRAN = -RANSEED
            do n = 1, 5
                do k = 1, 100
                    do i = 1, plon
                        do j = 1, plat
                            rannumbers_tot(i,k,j,n) =RAN1(IRAN)
                        end do
                    end do
                end do
            end do
        end if
        do n = 1, 5
            rannumbers_tot1 = rannumbers_tot(:,:,:,n)
            call scatter_field_to_chunk(1, 100, 1, plon, rannumbers_tot1, &
                rannumbers1)
            !Combine random field in each chunk
            rannumbers(:,(n-1)*100+1:n*100,:) = rannumbers1 
        end do
            

    end subroutine random_number_set


    REAL(r8) FUNCTION RAN1(IDUM)    

!This is contributed code standardized by Yong Wang

! Random number generator taken from Press et al.
!
! Returns numbers in the range 0-->1
!
! Their description...
! "Minimal" random number generator of Park and Miller with Bays-Durham
! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
! (exclusive of the endpoint values). Call with idum a negative integer to
! initialize; thereafter, do not alter idum between successive calls in a
! sequence. RNMX should approximate the largest floating value that is less
! than 1.
       
        use shr_kind_mod,    only: r8 => shr_kind_r8, i8 => shr_kind_i8
    
        INTEGER(i8), PARAMETER:: NTAB = 32,IQ = 127773,IA = 16807,IR = 2836, &
            IM = 2147483647,NDIV = 1+(IM-1)/NTAB

        REAL (r8), PARAMETER:: AM = 1.0/IM,EPS = 1.2E-7,RNMX = 1.0-EPS

        INTEGER (i8), INTENT(INOUT):: IDUM

        INTEGER (i8):: IY
        INTEGER (i8), DIMENSION(NTAB):: IV
        SAVE IV,IY
        DATA IV /NTAB*0/, IY /0/
        INTEGER (i8):: j,k    
        
        !
        IF (IDUM.LE.0.OR.IY.EQ.0) THEN
        ! Initalize
            IDUM = MAX(-IDUM,1)
            DO j = NTAB+8,1,-1
                k = IDUM/IQ
                IDUM = IA*(IDUM-k*IQ)-IR*k
                IF (IDUM.LT.0) IDUM = IDUM+IM
                IF (j.LE.NTAB) IV(j) = IDUM
            END DO
            IY = IV(1)
        END IF
        ! 
        k = IDUM/IQ
        ! Compute IDUM = MOD(IA*IDUM,IM) without overflows by Schrage's method
        IDUM = IA*(IDUM-k*IQ)-IR*k
        IF (IDUM.LT.0) IDUM = IDUM+IM
        ! j will be in the range 1-->NTAB
            j = 1+IY/NDIV
            ! Output previously stored value and refill the shuffle table
            IY = IV(j)
            IV(j) = IDUM
            RAN1 = MIN(AM*IY,RNMX)
    
        RETURN
    
    END FUNCTION RAN1

end module random_number_generator
