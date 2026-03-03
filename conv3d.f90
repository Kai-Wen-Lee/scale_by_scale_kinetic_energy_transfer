module conv3d
  use ncinput
  use domain_decomposition
  use io
  USE OMP_LIB
  implicit none
  contains
  subroutine convolution(u_in, k, u_out)
    implicit none

    real(8), allocatable, intent(out) :: u_out(:,:,:)
    real(8), allocatable :: u_bar(:)
    real(8), intent(in) :: k
    real(8), intent(in) :: u_in(dim_phi,dim_r,dim_z)
    real(8), parameter :: PI = 4.d0*datan(1.d0)
    real(8) :: Delta

    integer :: idxi, idxii, idxj, idxjj, idxk, idxkk
    integer :: idx, local_k
    real(8) :: win_lenx, win_leny, win_lenz
    real(8) :: window_size, GaussianKernel, Ubar

    allocate(u_bar(mpi_dim_r*mpi_dim_phi*mpi_dim_z))
    allocate(u_out(dim_phi, dim_r, dim_z))
    Delta = PI/k

  !$omp parallel do collapse(3) default(shared) private( &
  !$omp idxi,idxj,idxk,idxii,idxjj,idxkk, &
  !$omp idx,local_k,Ubar, &
  !$omp win_lenx,win_leny,win_lenz, &
  !$omp window_size,GaussianKernel )

    do idxi = 1, dim_r
      do idxj = 1, dim_phi
        do idxk = mpi_zidx_sta, mpi_zidx_end

          Ubar = 0.d0

          do idxii = 1, dim_r
            do idxjj = 1, dim_phi
              do idxkk = 1, dim_z

                win_lenx = r_p(idxi)*cos(phi(idxj)) - &
                           r_p(idxii)*cos(phi(idxjj))

                win_leny = r_p(idxi)*sin(phi(idxj)) - &
                           r_p(idxii)*sin(phi(idxjj))

                win_lenz = z_p(idxk) - z_p(idxkk)

                window_size = sqrt(win_lenx**2 + &
                                   win_leny**2 + &
                                   win_lenz**2)

                GaussianKernel = sqrt(6.d0/(PI*Delta**2)) * &
                                 exp(-6.d0*window_size**2/(Delta**2))

                Ubar = Ubar + GaussianKernel * &
                       u_in(idxjj, idxii, idxkk) * &
                       r_p(idxii) * dr_p(idxii) * dphi * dz_p(idxkk)

              end do
            end do
          end do

          local_k = idxk - mpi_zidx_sta + 1

          idx = local_k + mpi_dim_z * &
                ((idxj-1) + mpi_dim_phi*(idxi-1))

          u_bar(idx) = Ubar

        end do
      end do
    end do

  !$omp end parallel do

    ! ---- MPI only happens here (single-threaded) ----

    call sync_z(u_bar,dim_r*dim_phi*mpi_dim_z, &
                dim_r*dim_phi*dim_z, u_out)

  end subroutine convolution
end module conv3d
