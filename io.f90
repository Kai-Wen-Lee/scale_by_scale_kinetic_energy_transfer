module io
  USE mpi_context
  USE domain_decomposition
  implicit none
  contains
  ! add sync_rphi when free
  subroutine sync_z(send_buff, send_buff_size, recv_buff_size, synced_matrix)
  ! ----------------------- Gather results to process 0 and sync across process ---------- !
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
    real(8), allocatable :: recv_buff(:)
    integer, allocatable :: recv_count(:), displs_(:)
    real(8), intent(in) :: send_buff(:)
    real(8), intent(out) :: synced_matrix(:,:,:)
    integer, intent(in) :: send_buff_size, recv_buff_size
    integer :: r_idx_sta, r_idx_end
    integer :: phi_idx_sta, phi_idx_end
    integer :: z_idx_sta, z_idx_end
    integer :: idxi, idxj, idxk, idx
    integer :: proc_count, block_count

    allocate(recv_buff(recv_buff_size))
    allocate(recv_count(num_procs))
    allocate(displs_(num_procs))

    call MPI_GATHER(send_buff_size, 1, MPI_INTEGER, &
                    recv_count, 1, MPI_INTEGER, &
                    0, comm, ierr)

    if (proc_id == 0) then
      displs_(1) = 0
      do proc_count = 2, num_procs
        displs_(proc_count) = displs_(proc_count-1) + recv_count(proc_count)
      end do
    end if


   ! if (proc_id == 0) then
   !   print *, 'recv_count', recv_count
   ! end if
  ! --------------------------------------------------------------------------|
  ! ----------- Now gather data (send buff) from all processes to  -----------|
  ! ----------- recv_buff at process 0 ---------------------------------------|
  ! --------------------------------------------------------------------------|
   ! print *, 'send buff shape', shape(send_buff)
    call mpi_gatherv(send_buff, send_buff_size, mpi_real8,&
      recv_buff, recv_count, displs_, mpi_real8,&
      0, comm, ierr)
   !   print *, 'gather done for ', proc_id
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  ! ------------Finally, 1D array of recv_buff is organised ------------------|
  ! ------------into output matrix based on displacement array ---------------|
  ! --------------------------------------------------------------------------|
    if (proc_id == 0) then
      block_count = 0
!      allocate(synced_matrix(dim_r, dim_phi, dim_z))

      do proc_count = 0, num_procs-1
        call process_workload_z(1, dim_z, num_procs, proc_count, z_idx_sta, z_idx_end)
        block_count = displs_(proc_count+1)

        idx = 1
        do idxi = 1, dim_r
          do idxj = 1, dim_phi
            do idxk = z_idx_sta, z_idx_end
              synced_matrix(idxj, idxi, idxk) = recv_buff(block_count + idx)
              idx = idx + 1
            end do
          end do
        end do
      end do
    endif

  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  ! --------------------------------------------------------------------------|
  end subroutine sync_z
  subroutine flatten(u, uflat)
    real(8), intent(in) :: u(:,:,:)
    real(8), intent(out) :: uflat(:)
    integer :: idxi, idxj, idxk, idx
    integer :: z_idx_sta, z_idx_end

    call process_workload_z(1, dim_z, num_procs, proc_id, z_idx_sta, z_idx_end)
    idx = 1
    do idxi = 1, dim_r
      do idxj = 1, dim_phi
        do idxk = z_idx_sta, z_idx_end
          uflat(idx)=u(idxj, idxi, idxk) 
          idx = idx + 1
        end do
      end do
    end do
  end subroutine flatten
end module io
