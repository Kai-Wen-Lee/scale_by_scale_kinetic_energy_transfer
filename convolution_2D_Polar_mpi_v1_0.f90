subroutine convolution2DPolar(matrix, rp, phi, dr, dphi, Delta, comm, filtered_matrix)
  use mpi
  implicit none
!  include 'mpif.h'
  integer :: Nr, Nr2, Nphi, Nphi2
  integer :: ridx, ridx2, ridx3, phiidx, phiidx2, phiidx3
  integer :: phi_staidx, r_staidx, phi_endidx, r_endidx
  integer :: phi_staidx2, r_staidx2, phi_endidx2, r_endidx2
  integer :: ierrmpi, proc_id, num_procs, comm
  integer :: proc_count, block_count, idx
  real(8), intent(in) :: matrix(:,:), rp(:), phi(:), dr(:), dphi, Delta
  real(8), intent(out) :: filtered_matrix(size(matrix,1), size(matrix,2))
  real(8) :: Ubar, window_length, GaussianKernel
  real(8), parameter :: PI =4.d0*datan(1.d0)
  real(8), dimension(:), allocatable :: recv_buff, send_buff
  integer, dimension(:), allocatable :: displs_, recv_count

  ! mpi_comm should be initialised with mpi4py library
  ! -- INTEGER comm is mpi_comm_world
  call mpi_comm_size(comm, num_procs, ierrmpi)
  call mpi_comm_rank(comm, proc_id, ierrmpi)

  Nr = size(rp)
  Nphi = size(phi)

  ! Broadcast information needed across all processes
  call mpi_bcast(matrix, Nr*Nphi, mpi_double_precision, 0, comm, ierrmpi)
  call mpi_bcast(rp, Nr, mpi_double_precision, 0, comm, ierrmpi)
  call mpi_bcast(dr, Nr, mpi_double_precision, 0, comm, ierrmpi)
  call mpi_bcast(dphi, 1, mpi_double_precision, 0, comm, ierrmpi)
  call mpi_bcast(Nr, 1, mpi_integer, 0, comm, ierrmpi)
  call mpi_bcast(Nphi, 1, mpi_integer, 0, comm, ierrmpi)
  call mpi_bcast(Delta, 1, mpi_double_precision, 0, comm, ierrmpi)

  ! Stop if num_proc condition is not met
  ! Number of CPUs specified must be a perfect square
  ! Due to the way the domain is split
  if (mod(num_procs, num_procs/nint(sqrt(real(num_procs)))) /= 0) then
    print *, 'Number of processes MUST be perfect squares!'
    call mpi_abort(comm, 1, ierrmpi)
  endif

  ! Domain is split like so:
  ! _ | dim R                  |
  ! d |-----------|----------- |
  ! i |  CPU 3    |  CPU 3     |
  ! m |           |            |
  !   |-----------|------------|
  ! p |  CPU 1    |  CPU 2     |
  ! h |           |            |
  ! i |---------- | ---------- |

  call process_workload_r(1, Nr, num_procs, proc_id, r_staidx, r_endidx)
  call process_workload_phi(1, Nphi, num_procs, proc_id, phi_staidx, phi_endidx)

  Nr2 = r_endidx-r_staidx+1
  Nphi2 = phi_endidx-phi_staidx+1

!  print *, 'Process #',proc_id,' of ',num_procs,&
!          ' from Phi', phi_staidx,' to ', phi_endidx, 'and R', r_staidx,' to ', r_endidx,&
!          ' Nr', Nr2,' Nphi', Nphi2

  ! ------------------------------------------------------------------------- !
  ! -------------------------- do convolution here -------------------------- !
  ! ------------------------------------------------------------------------- !
  ! Each rank loops over its respective r and phi domain
  ! Where for each r and phi, an integral transform of the form
  !
  !        U_filtered(i,j) = sum_i sum_j G(r-i, phi-j) U(i,j) DA
  !
  ! is applied to each r,phi 
  ! ------------------------------------------------------------------------- !
  allocate(send_buff(Nr2*Nphi2))
  idx=1
  do phiidx = phi_staidx, phi_endidx
    do ridx = r_staidx, r_endidx
      Ubar = 0.d0
      do ridx2 = 1, Nr
        do phiidx2 = 1, Nphi
          window_length = sqrt((rp(ridx)*cos(phi(phiidx)) - rp(ridx2)*cos(phi(phiidx2)))**2.d0 &
                  + (rp(ridx)*sin(phi(phiidx)) - rp(ridx2)*sin(phi(phiidx2)))**2.d0)
          GaussianKernel = sqrt(6.d0/(PI*Delta**2.d0))*exp(-1.d0*(6.d0*window_length**2.d0)/(Delta**2.d0))
          Ubar = Ubar + GaussianKernel * matrix(ridx2,phiidx2) * rp(ridx2) * dr(ridx2) * dphi
        end do
      end do
      send_buff(idx) = Ubar
      idx= idx + 1
    end do
  end do
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !

  ! ---------- ---------- ------- MPI OPERATION ------- ---------- ---------- !
  ! ----------------------- Gather results to process 0 and output ---------- !
  ! ------------------------------------------------------------------------- !
  allocate(recv_buff(Nr*Nphi))
  allocate(recv_count(num_procs))
  allocate(displs_(num_procs))

  ! First gather buffer sizes and buffer position to process zero
  ! send_buff is 1D array of size Nr2*Nphi2
  ! recv_count is a 1D array of size num_procs and contains buffer position of send_buff
  ! from each process
  ! i.e.
  ! Fills displs_ which gives the displacement
  ! 
  ! |      receive buffer       |
  ! | rank 0 data | rank 1 data |
  ! |    num ele  |   num ele   |
  ! ^             ^             ^
  !               displs_= i    displs_= j
  ! --------------------------------------------------------------------------|
  call MPI_GATHER(Nr2*Nphi2, 1, MPI_INTEGER, &
                  recv_count, 1, MPI_INTEGER, &
                  0, comm, ierrmpi)
  if (proc_id == 0) then
    displs_(1) = 0
    do proc_count = 2, num_procs
      displs_(proc_count) = displs_(proc_count-1) + recv_count(proc_count)
    end do
  end if

  ! --------------------------------------------------------------------------|
  ! ----------- Now gather data (send buff) from all processes to  -----------|
  ! ----------- recv_buff at process 0 ---------------------------------------|
  ! --------------------------------------------------------------------------|
  call mpi_gatherv(send_buff, Nr2*Nphi2, mpi_real8,&
    recv_buff, recv_count, displs_, mpi_real8,&
    0, comm, ierrmpi)
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|

  ! ------------Finally, 1D array of recv_buff is organised ------------------|
  ! ------------into output matrix based on displacement array ---------------|
  ! --------------------------------------------------------------------------|
  if (proc_id == 0) then
    block_count = 0

    do proc_count = 0, num_procs-1
      call process_workload_r(1, Nr, num_procs, proc_count, r_staidx2, r_endidx2)
      call process_workload_phi(1, Nphi, num_procs, proc_count, phi_staidx2, phi_endidx2)

      block_count = displs_(proc_count+1)
      idx = 1
      do phiidx3 = phi_staidx2,phi_endidx2
        do ridx3 = r_staidx2, r_endidx2
          filtered_matrix(ridx3, phiidx3) = recv_buff(block_count + idx)
          idx = idx + 1
        end do
      end do
    end do
  endif
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  deallocate(recv_buff)
  deallocate(send_buff)
!  print *, 'Proc ',proc_id,' done'
end subroutine

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine process_workload_phi(n1, n2, nprocs, irank, ista, iend)
  implicit none
  integer, intent(in) :: n1, n2, nprocs, irank
  integer, intent(out):: ista, iend
  integer :: iwork1, iwork2, nprocs_sqrt, irank2
  !split the domain into num_procs areas
  !so then each dimensions should have num_procs**0.5 sections
  !and the num_procs MUST be perfect squares more than 1 
  !i.e. 4,9,16,25,36,49,64,81,100,121,144,169,...
  nprocs_sqrt = nint(sqrt(real(nprocs)))
  irank2 = floor(real(irank/nprocs_sqrt))

  iwork1 = (n2 - n1 + 1) / nprocs_sqrt !divide into sections
  iwork2 = MOD(n2 - n1 + 1, nprocs_sqrt) !get remainders
  ista = irank2 * iwork1 + n1 + MIN(irank2, iwork2)
  iend = ista + iwork1 - 1 
  if (iwork2 > irank2) iend = iend + 1 

end subroutine process_workload_phi
!-------------------------------------------------------------------
!-----------------------------------------------------------------
subroutine process_workload_r(n1, n2, nprocs, irank, ista, iend)
  implicit none
  integer, intent(in)  :: n1, n2, nprocs, irank
  integer, intent(out) :: ista, iend
  integer :: iwork1, iwork2, nprocs_sqrt, irank2

  ! number of processors per dimension
  nprocs_sqrt = nint(sqrt(real(nprocs)))

  ! column index (phi direction)
  irank2 = mod(irank, nprocs_sqrt)

  iwork1 = (n2 - n1 + 1) / nprocs_sqrt
  iwork2 = mod(n2 - n1 + 1, nprocs_sqrt)

  ista = irank2 * iwork1 + n1 + min(irank2, iwork2)
  iend = ista + iwork1 - 1
  if (iwork2 > irank2) iend = iend + 1

end subroutine process_workload_r
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------

