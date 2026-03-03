module domain_decomposition
  use ncinput
  use mpi_context
  implicit none
  integer :: num_procs_z, num_procs_rphi
  integer :: mpi_ridx_sta, mpi_ridx_end
  integer :: mpi_phiidx_sta, mpi_phiidx_end
  integer :: mpi_zidx_sta, mpi_zidx_end
  integer :: mpi_dim_r, mpi_dim_phi, mpi_dim_z

  contains
  subroutine split_domain()
    implicit none
    num_procs_rphi = 4
    if (num_procs >= num_procs_rphi * dim_z) then
      num_procs_z = num_procs - num_procs_rphi * dim_z
      call process_workload_z(1, dim_z, num_procs_z, proc_id, mpi_zidx_sta, mpi_zidx_end)
      call process_workload_phi(1, dim_phi, num_procs_rphi, proc_id, mpi_phiidx_sta, mpi_phiidx_end)
      call process_workload_r(1, dim_r, num_procs_rphi, proc_id, mpi_ridx_sta, mpi_ridx_end)
    else 
      call process_workload_z(1, dim_z, num_procs, proc_id, mpi_zidx_sta, mpi_zidx_end)
      mpi_ridx_sta = 1
      mpi_ridx_end = dim_r
      mpi_phiidx_sta = 1
      mpi_phiidx_end = dim_phi
    end if

    mpi_dim_r = mpi_ridx_end - mpi_ridx_sta + 1
    mpi_dim_phi = mpi_phiidx_end - mpi_phiidx_sta + 1
    mpi_dim_z = mpi_zidx_end - mpi_zidx_sta + 1
  end subroutine split_domain

  subroutine process_workload_z(n1, n2, nprocs, irank, ista, iend)
! Added by Mohammad Emran on June 01, 2016.
    implicit none
    integer, intent(in) :: n1, n2, nprocs, irank
    integer, intent(out):: ista, iend
    integer :: iwork1, iwork2
  
    iwork1 = (n2 - n1 + 1) / nprocs
    iwork2 = MOD(n2 - n1 + 1, nprocs)
    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
    iend = ista + iwork1 - 1 
    if (iwork2 > irank) iend = iend + 1 
  
  end subroutine process_workload_z
!-------------------------------------------------------------------
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

end module domain_decomposition
