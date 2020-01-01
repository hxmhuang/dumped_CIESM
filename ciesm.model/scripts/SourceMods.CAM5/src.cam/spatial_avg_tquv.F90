module spatial_avg_tquv

    use ppgrid,        only: begchunk, endchunk, pcols, pver
    use physics_types, only: physics_state
    use shr_kind_mod,  only: r8 => shr_kind_r8

    implicit none
    public
    save


    ! Public methods
    public :: spatial_avg_tquv_init
    public :: spatial_avg_tquv_set
    
    real(r8), allocatable, dimension(:,:,:) :: t_spatial_avg
    real(r8), allocatable, dimension(:,:,:) :: q_spatial_avg
    real(r8), allocatable, dimension(:,:,:) :: u_spatial_avg
    real(r8), allocatable, dimension(:,:,:) :: v_spatial_avg

contains

    subroutine spatial_avg_tquv_init()
        
        use abortutils,   only: endrun
        use cam_logfile,  only : iulog

        implicit none
        
        integer :: istat

        allocate(t_spatial_avg(pcols, begchunk:endchunk, pver), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'spatial_avg_tquv_init: failed to allocate t_spatial_avg; error = ',istat
            call endrun
        end if

        allocate(q_spatial_avg(pcols, begchunk:endchunk, pver), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'spatial_avg_tquv_init: failed to allocate q_spatial_avg; error = ',istat
            call endrun
        end if

        allocate(u_spatial_avg(pcols, begchunk:endchunk, pver), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'spatial_avg_tquv_init: failed to allocate u_spatial_avg; error = ',istat
            call endrun
        end if

        allocate(v_spatial_avg(pcols, begchunk:endchunk, pver), stat=istat)
        if( istat /= 0 ) then
            write(iulog,*) 'spatial_avg_tquv_init: failed to allocate v_spatial_avg; error = ',istat
            call endrun
        end if

    end subroutine spatial_avg_tquv_init

    subroutine spatial_avg_tquv_set(state)
        
        use phys_grid,    only: gather_chunk_to_field, scatter_field_to_chunk
        use dyn_grid,     only: get_horiz_grid_dim_d, get_dyn_grid_parm
        use spmd_utils,   only: masterproc

    !
    ! Arguments
    !
        type(physics_state), intent(in) :: state(begchunk:endchunk)

    
    ! Local workspace
    
        real(r8) :: t(pcols,begchunk:endchunk,pver)         ! input array, chunked
        real(r8) :: q(pcols,begchunk:endchunk,pver)         ! input array, chunked
        real(r8) :: u(pcols,begchunk:endchunk,pver)         ! input array, chunked
        real(r8) :: v(pcols,begchunk:endchunk,pver)         ! input array, chunked
    
        integer :: hdim1, hdim2                             ! dimensions of rectangular horizontal 
                                                            ! grid data structure, If 1D data 
                                                            ! structure, then hdim2_d == 1.
        integer :: ngcols                                   ! global column count (all)
        integer :: i, j, k, n                               ! longitude, latitude, level, 
                                                            ! and global column indices
        integer :: c
        integer :: ii, jj
        integer :: plon, plat
        integer :: count
    
        !integer, parameter :: ninp_av = 3
        integer, parameter :: ninp_av = 1
    
        ! rectangular version of t, q, u, v
        real(r8), allocatable :: t_field(:,:,:)
        real(r8), allocatable :: q_field(:,:,:)
        real(r8), allocatable :: u_field(:,:,:)
        real(r8), allocatable :: v_field(:,:,:)
    
        real(r8), allocatable :: t_field_avg(:,:,:)
        real(r8), allocatable :: q_field_avg(:,:,:)
        real(r8), allocatable :: u_field_avg(:,:,:)
        real(r8), allocatable :: v_field_avg(:,:,:)
        
        plon = get_dyn_grid_parm('plon')
        plat = get_dyn_grid_parm('plat')
    
        do c = begchunk, endchunk
            t(:,c,:) = state(c)%t(:,:)
            q(:,c,:) = state(c)%q(:,:,1)
            u(:,c,:) = state(c)%u(:,:)
            v(:,c,:) = state(c)%v(:,:)
        end do
    
        call get_horiz_grid_dim_d(hdim1, hdim2)
    
        allocate(t_field(hdim1,hdim2,pver))
        allocate(q_field(hdim1,hdim2,pver))
        allocate(u_field(hdim1,hdim2,pver))
        allocate(v_field(hdim1,hdim2,pver))
    
        allocate(t_field_avg(hdim1,hdim2,pver))
        allocate(q_field_avg(hdim1,hdim2,pver))
        allocate(u_field_avg(hdim1,hdim2,pver))
        allocate(v_field_avg(hdim1,hdim2,pver))
    
        t_field(:,:,:) = 0.0_r8
        q_field(:,:,:) = 0.0_r8
        u_field(:,:,:) = 0.0_r8
        v_field(:,:,:) = 0.0_r8
    
        call gather_chunk_to_field (1, 1, pver, hdim1, t, t_field)
        call gather_chunk_to_field (1, 1, pver, hdim1, q, q_field)
        call gather_chunk_to_field (1, 1, pver, hdim1, u, u_field)
        call gather_chunk_to_field (1, 1, pver, hdim1, v, v_field)
        
        t_field_avg = t_field
        q_field_avg = q_field
        u_field_avg = u_field
        v_field_avg = v_field
    
        if (masterproc) then
    
            ngcols = hdim1*hdim2
            do k=1,pver
                do j=1,hdim2
                    do i=1,hdim1
                       count = 0
                       t_field_avg(i,j,k) = 0.0_r8
                       q_field_avg(i,j,k) = 0.0_r8
                       u_field_avg(i,j,k) = 0.0_r8
                       v_field_avg(i,j,k) = 0.0_r8
                        do ii = -ninp_av, ninp_av
                            if ( (i+ii) .lt. 1 .or. (i+ii) .gt. hdim1 ) then
                                cycle
                            end if
                            do jj = -ninp_av, ninp_av
                                if ( (j+jj) .lt. 1 .or. (j+jj) .gt. hdim2 ) then
                                    cycle
                                end if
                                t_field_avg(i,j,k) = t_field_avg(i,j,k) + t_field(i+ii,j+jj,k)
                                q_field_avg(i,j,k) = q_field_avg(i,j,k) + q_field(i+ii,j+jj,k)
                                u_field_avg(i,j,k) = u_field_avg(i,j,k) + u_field(i+ii,j+jj,k)
                                v_field_avg(i,j,k) = v_field_avg(i,j,k) + v_field(i+ii,j+jj,k)
                                count = count + 1
                            end do
                        end do
                        t_field_avg(i,j,k) = t_field_avg(i,j,k)/count
                        q_field_avg(i,j,k) = q_field_avg(i,j,k)/count
                        u_field_avg(i,j,k) = u_field_avg(i,j,k)/count
                        v_field_avg(i,j,k) = v_field_avg(i,j,k)/count
                        count = 0
                    end do
                end do
            end do
    
    
         end if
    
         call scatter_field_to_chunk(1, 1, pver, plon, t_field_avg, t_spatial_avg)
         call scatter_field_to_chunk(1, 1, pver, plon, q_field_avg, q_spatial_avg)
         call scatter_field_to_chunk(1, 1, pver, plon, u_field_avg, u_spatial_avg)
         call scatter_field_to_chunk(1, 1, pver, plon, v_field_avg, v_spatial_avg)
    
    end subroutine spatial_avg_tquv_set


end module spatial_avg_tquv
