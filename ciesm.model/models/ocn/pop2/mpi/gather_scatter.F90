!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: gather_scatter
!
! Workaround for performance issue on Cray/Gemini
!
#ifdef _NO_MPI_RSEND
#define MPI_RSEND MPI_SEND
#define mpi_rsend mpi_send
#define MPI_IRSEND MPI_ISEND
#define mpi_irsend mpi_isend
#endif

 module gather_scatter

! !DESCRIPTION:
!  This module contains routines for gathering data to a single
!  processor from a distributed array, scattering data from a
!  single processor to a distributed array and changing distribution
!  of blocks of data (eg from baroclinic to barotropic and back).
!
! !REVISION HISTORY:
!  SVN: $Id: gather_scatter.F90 31545 2011-11-01 17:46:11Z jedwards $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use blocks
   use distribution
   use domain
   use domain_size
   use exit_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: gather_global,      &
             scatter_global,     & 
             redistribute_blocks, &
             gather_new_dbl,scatter_new_int,scatter_new_dbl


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  overload module functions
!
!-----------------------------------------------------------------------

   interface gather_global
     module procedure gather_global_dbl,  &
                      gather_global_real, &
                      gather_global_int
   end interface 

   interface scatter_global
     module procedure scatter_global_dbl,  &
                      scatter_global_real, &
                      scatter_global_int,  &
                      scatter_global_log
   end interface 

   interface redistribute_blocks
     module procedure redistribute_blocks_dbl,  &
                      redistribute_blocks_real, &
                      redistribute_blocks_int
   end interface 

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:

subroutine logpr(a, logp, r)

  integer (int_kind), intent(in) :: a
  integer (int_kind), intent(inout) :: logp, r
  integer (int_kind) :: val, pp
  pp = 0
  val = 1
  do while (val <= a) 
    val = val * 2
    pp = pp + 1
  enddo
  logp = pp - 1
  r = a - (val / 2)
end subroutine

subroutine gather_new_dbl(sendbuf,recvbuf,root,comm,ierr) 
   include "mpif.h"
     real (r8), dimension(:,:,:), intent(in) :: sendbuf
     real (r8), dimension(:,:,:,:), intent(inout):: recvbuf
     integer (int_kind) :: pes, typesize
     integer (int_kind), dimension(:), allocatable :: counts, recv_off, &
     recv_count
     real(r8), dimension(:), allocatable :: bufp
     integer (int_kind) :: i, j, buf_size, logp, r, base, need_malloc , &
     recv_times, send_times, relative_myid, target1, req, tmpcnt,off,&
     flag,k,m, root, comm, ierr, sta(mpi_status_size),n
     

     pes = nblocks_tot
     typesize = 8
     need_malloc = 0
     recv_times = 0
     send_times = 0
     call logpr(pes, logp, r)
     buf_size = 0
     !tmpcnt = nx_block * ny_block
     tmpcnt = 360
     if(r /= 0) logp = logp + 1
     allocate (recv_off(0 : pes - 1), recv_count(0 : pes - 1))
     recv_off(0) = tmpcnt
     buf_size = tmpcnt
     relative_myid = my_task - root
     if(relative_myid < 0) relative_myid = relative_myid + pes
     !print *, "mytask=",my_task, pes, logp, r, tmpcnt
     base = 1
     do i = 1, logp
       if( (mod(relative_myid / base, 2 ) == 0) .and. (mod(relative_myid, base) == 0) &
           .and. (relative_myid + base < pes)) then
           recv_count(recv_times) = 0
           do j=base, base * 2 - 1
             if(relative_myid + j < pes) then
               recv_count(recv_times) = recv_count(recv_times) + tmpcnt
             endif
           enddo
           buf_size = buf_size + recv_count(recv_times)
           recv_off(recv_times + 1) = buf_size
           recv_times = recv_times + 1
       else if((mod(relative_myid / base, 2) == 1) .and. (mod(relative_myid, base) == 0)) then
         send_times = send_times + 1
       endif
       base = base * 2
     enddo

     flag = 0
     if(recv_times /= 0) then
       allocate (bufp(0 : buf_size - 1))
       flag = 1
     endif
     recv_times = 0
     base = 1
     do i=1,logp
       if((mod(relative_myid / base, 2)) == 0 .and. (mod(relative_myid, base)) == 0 &
         .and. (relative_myid + base < pes)) then
         call mpi_irecv(bufp(recv_off(recv_times)), recv_count(recv_times), mpi_double,&
         mod(my_task + base,pes), my_task, comm, req, ierr)
         if(recv_times == 0) then

           do j = 1, 2
             do k = 1, 3
               do n = 1,60
               bufp((j - 1) * 180 + (k - 1) * 60 + n - 1) = sendbuf(n, k, j)
               enddo
             enddo
           enddo
         endif
         call mpi_wait(req, sta, ierr)
         recv_times = recv_times + 1

       else if((mod(relative_myid / base, 2) == 1) .and. (mod(relative_myid, base) == 0)) then
         if(recv_times /=0) then
           call MPI_Send(bufp, buf_size, MPI_DOUBLE, mod(my_task - base + pes, pes), & 
           mod(my_task - base + pes, pes), comm, ierr)
         else 
           call MPI_Send(sendbuf, tmpcnt, MPI_DOUBLE, mod(my_task - base + pes, pes), &
           mod(my_task - base + pes, pes), comm, ierr)
         endif

       endif
       base = base * 2
     enddo
     off = 0
     if(my_task == root) then
       i = root
       do j = 0, pes - 1
         do k = 1, 2
           do m = 1, 3
             do n = 1,60
             recvbuf(n,m, k, i + 1) = bufp(off)
             off = off + 1
           enddo
           enddo
         enddo
         i = mod(i + 1,pes)
       enddo
    
     endif
     
     deallocate(recv_off, recv_count)
     if(flag == 1) deallocate (bufp)

  end subroutine

 subroutine scatter_new_int(ARRAY, ARRAY_G, src_task, dst_dist)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status
   integer (int_kind), dimension(:,:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0
   if (my_task == src_task) then

      allocate (msg_buffer(nx_block,ny_block, nblocks_tot))
      msg_buffer = 0
      do n=1,nblocks_tot

        !*** copy local blocks
!
        if (dst_dist%proc(n) == my_task+1) then
!
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY(i,j) = ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) 
          end do
          end do

        !*** fill land blocks with zeroes
        endif

        if (dst_dist%proc(n) == 0) then

          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            msg_buffer(i, j, n) = undefined_nf_int
          end do
          end do
        else 
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            msg_buffer(i, j, n) =  ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j))
          end do
          end do

        endif

      end do
    endif
    call mpi_scatter(msg_buffer, nx_block*ny_block,MPI_INT, &
    array, nx_block*ny_block,MPI_INT,src_task,MPI_COMM_OCN,ierr)

!
    if (my_task == src_task) then
      deallocate(msg_buffer)
    endif



!-----------------------------------------------------------------------

 end subroutine scatter_new_int

 subroutine scatter_new_dbl(ARRAY, ARRAY_G, src_task, dst_dist)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r8), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real(r8), dimension(:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status
   real(r8), dimension(:,:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = 0
   if (my_task == src_task) then

      allocate (msg_buffer(nx_block,ny_block, nblocks_tot))
      msg_buffer = 0
      do n=1,nblocks_tot

        !*** copy local blocks
!
        if (dst_dist%proc(n) == my_task+1) then
!
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY(i,j) = ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) 
          end do
          end do

        !*** fill land blocks with zeroes
        endif

        if (dst_dist%proc(n) == 0) then

          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            msg_buffer(i, j, n) = undefined_nf
          end do
          end do
        else 
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            msg_buffer(i, j, n) =  ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j))
          end do
          end do

        endif

      end do
    endif
    call mpi_scatter(msg_buffer, nx_block*ny_block,MPI_DOUBLE, &
    array, nx_block*ny_block,MPI_DOUBLE,src_task,MPI_COMM_OCN,ierr)

!
    if (my_task == src_task) then
      deallocate(msg_buffer)
    endif



!-----------------------------------------------------------------------

 end subroutine scatter_new_dbl



 subroutine gather_global_dbl_orig(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays 
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (r8), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   if (my_task == dst_task) then

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(n) == my_task+1) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = &
                  ARRAY(i,j,src_dist%local_block(n))
         end do
         end do

       !*** fill land blocks with zeroes

       else if (src_dist%proc(n) == 0) then

         this_block = get_block(n,n)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = c0
         end do
         end do
       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (src_dist%proc(n) > 0 .and. &
           src_dist%proc(n) /= my_task+1) then

         this_block = get_block(n,n)

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_dbl, src_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

         do j=this_block%jb,this_block%je
         do i=this_block%ib,this_block%ie
           ARRAY_G(this_block%i_glob(i), &
                   this_block%j_glob(j)) = msg_buffer(i,j)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))

     nsends = 0
     do n=1,nblocks_tot
       if (src_dist%proc(n) == my_task+1) then

         nsends = nsends + 1
         src_block = src_dist%local_block(n)
         call MPI_ISEND(ARRAY(1,1,src_block), nx_block*ny_block, &
                     mpi_dbl, dst_task, 3*mpitag_gs+n, &
                     MPI_COMM_OCN, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)

   endif

!-----------------------------------------------------------------------

 end subroutine gather_global_dbl_orig

!***********************************************************************
!BOP
! !IROUTINE: gather_global
! !INTERFACE:


subroutine gather_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

    include 'mpif.h'

! !INPUT PARAMETERS:

    integer (int_kind), intent(in) :: &
      dst_task   ! task to which array should be gathered

    type (distrb), intent(in) :: &
      src_dist   ! distribution of blocks in the source array

    real (r8), dimension(:,:,:), intent(in) :: &
      ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:

    real (r8), dimension(:,:), intent(inout) :: &
      ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer (int_kind) :: &
      i,j,n          ,&! dummy loop counters
      nsends         ,&! number of actual sends
      src_block      ,&! block locator for send
      ierr             ! MPI error flag

    integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
      status

    integer (int_kind), dimension(:), allocatable :: &
      snd_request

    integer (int_kind), dimension(:,:), allocatable :: &
      snd_status

    real (r8), dimension(:,:,:), allocatable :: &
      msg_buffer

    type (block) :: &
      this_block  ! block info for current block
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj


  !  integer (int_kind), save :: &
  !      gather_timer_dbl
  !  logical (log_kind), save :: &
  !      first1 = .true.       ! flag for initializing timers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef _USE_FLOW_CONTROL
    integer (int_kind) :: &
      rcv_request    ,&! request id
      signal           ! MPI handshaking variable

    signal = 1
#endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!    if (first1) then
!      call get_timer(gather_timer_dbl, 'gather_timer_dbl',nblocks_clinic,distrb_clinic%nprocs)
!      first1 = .false.
!    endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (my_task == dst_task) then

      do n=1,nblocks_tot

        !*** copy local blocks
!
        if (src_dist%proc(n) == my_task+1) then
!
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                   ARRAY(i,j,src_dist%local_block(n))
          end do
          end do

        !*** fill land blocks with zeroes

        else if (src_dist%proc(n) == 0) then

          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = undefined_nf
          end do
          end do
        endif

      end do
!
!      !*** receive blocks to fill up the rest
!
      allocate (msg_buffer(nx_block,ny_block, nblocks_tot))
    endif
    !print *, "gather_new_dbl is begin ", my_task
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_start(gather_timer_dbl)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    call mpi_gather(array, nx_block*ny_block,MPI_DOUBLE, &
    msg_buffer, nx_block*ny_block,MPI_DOUBLE,dst_task,MPI_COMM_OCN,ierr)
    !call gather_new_dbl(array, msg_buffer,dst_task, MPI_COMM_OCN, ierr)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_stop(gather_timer_dbl)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !print *, "gather_new_dbl is end ", my_task

!
    if (my_task == dst_task) then
      do n=1,nblocks_tot
        if (src_dist%proc(n) > 0 .and. &
            src_dist%proc(n) /= my_task+1) then

          this_block = get_block(n,n)
          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = msg_buffer(i,j,n)
          end do
          end do
        endif
      enddo
      deallocate(msg_buffer)
    endif

end subroutine gather_global_dbl
!***********************************************************************

 subroutine gather_global_real(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   real (r4), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r4), dimension(:,:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj


   ! integer (int_kind), save :: &
   !     gather_timer_real
   ! logical (log_kind), save :: &
   !     first1 = .true.       ! flag for initializing timers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     rcv_request    ,&! request id
     signal           ! MPI handshaking variable

   signal = 1
#endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   ! if (first1) then
   !   call get_timer(gather_timer_real, 'gather_timer_real',nblocks_clinic,distrb_clinic%nprocs)
   !   first1 = .false.
   ! endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------
    if (my_task == dst_task) then

      do n=1,nblocks_tot

        !*** copy local blocks
!
        if (src_dist%proc(n) == my_task+1) then
!
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                   ARRAY(i,j,src_dist%local_block(n))
          end do
          end do

        !*** fill land blocks with zeroes

        else if (src_dist%proc(n) == 0) then

          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = undefined_nf_r4
          end do
          end do
        endif

      end do
!
!      !*** receive blocks to fill up the rest
!
      allocate (msg_buffer(nx_block,ny_block, nblocks_tot))
    endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_start(gather_timer_real)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    call mpi_gather(array, nx_block*ny_block,MPI_REAL, &
    msg_buffer, nx_block*ny_block,MPI_REAL,dst_task,MPI_COMM_OCN,ierr)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_stop(gather_timer_real)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!
    if (my_task == dst_task) then
      do n=1,nblocks_tot
        if (src_dist%proc(n) > 0 .and. &
            src_dist%proc(n) /= my_task+1) then

          this_block = get_block(n,n)
          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = msg_buffer(i,j,n)
          end do
          end do
        endif
      enddo
      deallocate(msg_buffer)
    endif


!-----------------------------------------------------------------------

 end subroutine gather_global_real

!***********************************************************************

 subroutine gather_global_int(ARRAY_G, ARRAY, dst_task, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   integer (int_kind), dimension(:,:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj


   ! integer (int_kind), save :: &
   !     gather_timer_int
   ! logical (log_kind), save :: &
   !     first1 = .true.       ! flag for initializing timers
   !   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#ifdef _USE_FLOW_CONTROL
   integer (int_kind) :: &
     rcv_request    ,&! request id
     signal           ! MPI handshaking variable

   signal = 1
#endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
   ! if (first1) then
   !   call get_timer(gather_timer_int, 'gather_timer_int',nblocks_clinic,distrb_clinic%nprocs)
   !   first1 = .false.
   ! endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------
    if (my_task == dst_task) then

      do n=1,nblocks_tot

        !*** copy local blocks
!
        if (src_dist%proc(n) == my_task+1) then
!
          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = &
                   ARRAY(i,j,src_dist%local_block(n))
          end do
          end do

        !*** fill land blocks with zeroes

        else if (src_dist%proc(n) == 0) then

          this_block = get_block(n,n)

          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = undefined_nf_int
          end do
          end do
        endif

      end do
!
!      !*** receive blocks to fill up the rest
!
      allocate (msg_buffer(nx_block,ny_block, nblocks_tot))
    endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_start(gather_timer_int)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    call mpi_gather(array, nx_block*ny_block,MPI_INT, &
    msg_buffer, nx_block*ny_block,MPI_INT,dst_task,MPI_COMM_OCN,ierr)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>wwj
!call timer_stop(gather_timer_int)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!
    if (my_task == dst_task) then
      do n=1,nblocks_tot
        if (src_dist%proc(n) > 0 .and. &
            src_dist%proc(n) /= my_task+1) then

          this_block = get_block(n,n)
          do j=this_block%jb,this_block%je
          do i=this_block%ib,this_block%ie
            ARRAY_G(this_block%i_glob(i), &
                    this_block%j_glob(j)) = msg_buffer(i,j,n)
          end do
          end do
        endif
      enddo
      deallocate(msg_buffer)
    endif


 end subroutine gather_global_int

!EOC
!***********************************************************************
!BOP
! !IROUTINE: scatter_global
! !INTERFACE:

 subroutine scatter_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.  Note that this routine only
!  is guaranteed to scatter correct values in the physical domain
!  on each block.  Ghost cells are not filled correctly for tripole
!  boundary conditions due to the complexity of dealing with field
!  locations.  It is wise to call the boundary update routine after
!  each scatter call.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r8), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r8), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_dbl, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         dst_block = dst_dist%local_block(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%local_block(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_dbl, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_dbl

!***********************************************************************

 subroutine scatter_global_real(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (r4), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_real, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         dst_block = dst_dist%local_block(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%local_block(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_real, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_real

!***********************************************************************

 subroutine scatter_global_int(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   integer (int_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = c0

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = c0
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_integer, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         dst_block = dst_dist%local_block(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%local_block(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_integer, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_int

!EOC
!***********************************************************************

 subroutine scatter_global_log(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

!-----------------------------------------------------------------------
!
!  This subroutine scatters a distributed array to a global-sized
!  array on the processor src_task.
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   logical (log_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole bounday conditions
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   logical (log_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to .false. and set up tripole quantities
!
!-----------------------------------------------------------------------

   ARRAY = .false.

   select case (field_loc)
   case (field_loc_center)   ! cell center location
      xoffset = 1
      yoffset = 1
   case (field_loc_NEcorner)   ! cell corner (velocity) location
      xoffset = 0
      yoffset = 0
   case (field_loc_Eface)   ! cell center location
      xoffset = 0
      yoffset = 1
   case (field_loc_Nface)   ! cell corner (velocity) location
      xoffset = 1
      yoffset = 0
   case (field_loc_noupdate) ! ghost cells never used - use cell center
      xoffset = 1
      yoffset = 1
   end select

   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call exit_POP(sigAbort, 'Unknown field type in scatter')
   end select

!-----------------------------------------------------------------------
!
!  if this task is the src_task, copy blocks of global array into 
!  message buffer and send to other processors. also copy local blocks
!
!-----------------------------------------------------------------------

   if (my_task == src_task) then

     !*** send non-local blocks away

     allocate (msg_buffer(nx_block,ny_block))

     do n=1,nblocks_tot
       if (dst_dist%proc(n) > 0 .and. &
           dst_dist%proc(n)-1 /= my_task) then

         msg_buffer = .false.
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                         this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        msg_buffer(i,j) = ARRAY_G(this_block%i_glob(i),&
                                                  this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        msg_buffer(i,j) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif

         call MPI_SEND(msg_buffer, nx_block*ny_block, &
                       mpi_logical, dst_dist%proc(n)-1, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, status, ierr)

       endif
     end do

     deallocate(msg_buffer)

     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         dst_block = dst_dist%local_block(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

            do j=1,ny_block
            do i=1,nx_block
               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                              this_block%j_glob(j))
            end do
            end do

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole

                  jsrc = ny_global + yoffset + &
                         (this_block%j_glob(j) + ny_global)
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        isrc = nx_global + xoffset - this_block%i_glob(i)
                        if (isrc < 1) isrc = isrc + nx_global
                        if (isrc > nx_global) isrc = isrc - nx_global
                        ARRAY(i,j,dst_block) = ARRAY_G(isrc,jsrc)
                     endif
                  end do

               endif
            end do

         endif
       endif
     end do

!-----------------------------------------------------------------------
!
!  otherwise receive data from src_task
!
!-----------------------------------------------------------------------

   else

     allocate (rcv_request(nblocks_tot), &
               rcv_status(MPI_STATUS_SIZE, nblocks_tot))

     rcv_request = 0
     rcv_status  = 0

     nrecvs = 0
     do n=1,nblocks_tot
       if (dst_dist%proc(n) == my_task+1) then
         nrecvs = nrecvs + 1
         dst_block = dst_dist%local_block(n)
         call MPI_IRECV(ARRAY(1,1,dst_block), nx_block*ny_block, &
                       mpi_logical, src_task, 3*mpitag_gs+n, &
                       MPI_COMM_OCN, rcv_request(nrecvs), ierr)
       endif
     end do

     if (nrecvs > 0) &
       call MPI_WAITALL(nrecvs, rcv_request, rcv_status, ierr)

     deallocate(rcv_request, rcv_status)
   endif

!-----------------------------------------------------------------------

 end subroutine scatter_global_log

!EOC
!***********************************************************************
!BOP
! !IROUTINE: redistribute_blocks
! !INTERFACE:

 subroutine redistribute_blocks_dbl(DST_ARRAY, dst_dist, &
                                    SRC_ARRAY, src_dist)

! !DESCRIPTION:
!  This subroutine converts an array distributed in a one decomposition
!  an array distributed in a different decomposition
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays 
!  corresponding to the generic interface scatter_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific interface based
!  on the data type of the input argument).

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
     src_dist    ,&! info on distribution of blocks for source array
     dst_dist      ! info on distribution of blocks for dest   array

   real (r8), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

! !OUTPUT PARAMETERS:

   real (r8), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,nblocks_tot
     src_task = src_dist%proc(n) - 1
     if (dst_dist%proc(n) == my_task+1 .and. src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = dst_dist%local_block(n)

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_dbl, src_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,nblocks_tot

     dst_task = dst_dist%proc(n) - 1

     if (src_dist%proc(n) == my_task+1 .and. dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = src_dist%local_block(n)

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_dbl, dst_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

     if (src_dist%proc(n) == my_task+1 .and. &
         dst_dist%proc(n) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(n)) = &
       SRC_ARRAY(:,:,src_dist%local_block(n))

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_dbl

!***********************************************************************

 subroutine redistribute_blocks_real(DST_ARRAY, dst_dist, &
                                     SRC_ARRAY, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine converts an array distributed in a baroclinic
!  data decomposition to an array in a barotropic decomposition
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(in) :: &
     src_dist      ! info on distribution of blocks for source array

   real (r4), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(inout) :: &
     dst_dist      ! info on dist of blocks for destination array

   real (r4), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,nblocks_tot
     src_task = src_dist%proc(n) - 1
     if (dst_dist%proc(n) == my_task+1 .and. src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = dst_dist%local_block(n)

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_real, src_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,nblocks_tot

     dst_task = dst_dist%proc(n) - 1

     if (src_dist%proc(n) == my_task+1 .and. dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = src_dist%local_block(n)

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_real, dst_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

     if (src_dist%proc(n) == my_task+1 .and. &
         dst_dist%proc(n) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(n)) = &
       SRC_ARRAY(:,:,src_dist%local_block(n))

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_real

!***********************************************************************

 subroutine redistribute_blocks_int(DST_ARRAY, dst_dist, &
                                    SRC_ARRAY, src_dist)

!-----------------------------------------------------------------------
!
!  This subroutine converts an array distributed in a baroclinic
!  data decomposition to an array in a barotropic decomposition
!
!-----------------------------------------------------------------------

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(in) :: &
     src_dist      ! info on distribution of blocks for source array

   integer (int_kind), dimension(:,:,:), intent(in) :: &
     SRC_ARRAY     ! array containing field in source distribution

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   type (distrb), intent(inout) :: &
     dst_dist      ! info on dist of blocks for destination array

   integer (int_kind), dimension(:,:,:), intent(inout) :: &
     DST_ARRAY     ! array containing field in dest distribution

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,n, &
     src_task,       &! task where the block currently resides
     dst_task,       &! task where the block needs to end up
     src_blck,       &! blck where the block currently resides
     dst_blck,       &! blck where the block needs to end up
     num_sends,      &! number of messages sent from this task
     num_recvs,      &! number of messages received by this task
     ierr             ! MPI error flag

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request,    &! request array for receives
     snd_request      ! request array for sends

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status,     &! status array for receives
     snd_status       ! status array for sends

!-----------------------------------------------------------------------
!
!  allocate space for asynchronous send/recv arrays
!
!-----------------------------------------------------------------------

   allocate (rcv_request(nblocks_tot), &
             snd_request(nblocks_tot), &
             rcv_status(MPI_STATUS_SIZE, nblocks_tot), &
             snd_status(MPI_STATUS_SIZE, nblocks_tot))

   rcv_request = 0
   snd_request = 0
   rcv_status = 0
   snd_status = 0

!-----------------------------------------------------------------------
!
!  first determine whether should be receiving messages and post all
!  the receives
!
!-----------------------------------------------------------------------

   num_recvs = 0
   do n=1,nblocks_tot
     src_task = src_dist%proc(n) - 1
     if (dst_dist%proc(n) == my_task+1 .and. src_task /= my_task) then

       num_recvs = num_recvs + 1
       dst_blck = dst_dist%local_block(n)

       call MPI_IRECV(DST_ARRAY(1,1,dst_blck), nx_block*ny_block, &
                      mpi_integer, src_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, rcv_request(num_recvs), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  now determine which sends are required and post the sends
!
!-----------------------------------------------------------------------

   num_sends = 0
   do n=1,nblocks_tot

     dst_task = dst_dist%proc(n) - 1

     if (src_dist%proc(n) == my_task+1 .and. dst_task /= my_task) then

       num_sends = num_sends + 1
       src_blck = src_dist%local_block(n)

       call MPI_ISEND(SRC_ARRAY(1,1,src_blck), nx_block*ny_block, &
                      mpi_integer, dst_task, 3*mpitag_gs+n, &
                      MPI_COMM_OCN, snd_request(num_sends), ierr)
     endif
   end do

!-----------------------------------------------------------------------
!
!  if blocks are local, simply copy the proper buffers
!
!-----------------------------------------------------------------------

   do n=1,nblocks_tot

     if (src_dist%proc(n) == my_task+1 .and. &
         dst_dist%proc(n) == my_task+1) then

       DST_ARRAY(:,:,dst_dist%local_block(n)) = &
       SRC_ARRAY(:,:,src_dist%local_block(n))

     endif
   end do

!-----------------------------------------------------------------------
!
!  finalize all the messages and clean up
!
!-----------------------------------------------------------------------

   if (num_sends /= 0) &
     call MPI_WAITALL(num_sends, snd_request, snd_status, ierr)
   if (num_recvs /= 0) &
     call MPI_WAITALL(num_recvs, rcv_request, rcv_status, ierr)

   deallocate (rcv_request, snd_request, rcv_status, snd_status)

!-----------------------------------------------------------------------

 end subroutine redistribute_blocks_int

!EOC
!***********************************************************************

 end module gather_scatter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
