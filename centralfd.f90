module centralfd
  USE ncinput
  USE io
  implicit none
  contains
  subroutine gradient(u, r, delta_i, delta_j, delta_k, didu, djdu, dkdu)
    real(8), intent(in) :: u(:,:,:) 
    real(8), intent(in) :: delta_i(:), delta_j, delta_k(:)
    integer :: idxi, idxj, idxk
    real(8), intent(out) :: didu(:,:,:), djdu(:,:,:), dkdu(:,:,:)
    real(8) :: bou_val_L, bou_val_R
    real(8), intent(in) :: r(:)
    do idxj = 1, dim_phi
      do idxk = 1, dim_z
        call get_interior(u(idxj,:,idxk), delta_i, didu(idxj, :,idxk))
        call get_boundary(u(idxj,:,idxk), delta_i, bou_val_L, bou_val_R)
        didu(idxj,1,idxk) = bou_val_L
        didu(idxj, dim_r, idxk) = bou_val_R
      end do
    end do
    do idxi = 1, dim_r
      do idxk = 1, dim_z
          call get_dphi(u(:,idxi,idxk), r(idxi), delta_j, djdu(:,idxi,idxk))
      end do
    end do
    do idxi = 1, dim_r
      do idxj = 1, dim_phi
        call get_interior(u(idxj,idxi,:), delta_k, dkdu(idxj, idxi,:))
        call get_boundary(u(idxj,idxi,:), delta_k, bou_val_L, bou_val_R)
        dkdu(idxj, idxi, 1) = bou_val_L
        dkdu(idxj, idxi, dim_z) = bou_val_R
      end do
    end do
  end subroutine gradient
  subroutine get_interior(u, delta_arr, du)
    ! Compute interior points by 2nd order central differencing
    real(8), intent(in) :: u(:)
    integer :: idx
    real(8), intent(out) :: du(:)
    real(8), intent(in) :: delta_arr(:)
    real(8) :: hs, hd
    do idx = 2, size(u)-1
      hs = delta_arr(idx)
      hd = delta_arr(idx+1)
      du(idx) = ((hs**2.d0)*u(idx+1) + (hd**2.d0 - hs**2.d0)*u(idx) - (hd**2.d0)*u(idx-1))/(hs*hd*(hd+hs))
    end do
  end subroutine get_interior
  subroutine get_boundary(f, dx, fout_L, fout_R)
    ! numerical difference 2nd order edges
    implicit none
    real(8), intent(in) :: f(:)
    real(8), intent(in) :: dx(:)
    real(8), intent(out) :: fout_L, fout_R
    !only output the boundary point
    real(8) :: dx1, dx2
    real(8) :: a, b, c

    ! ---- Left boundary ----
    dx1 = dx(1)
    dx2 = dx(2)

    a = -(2.d0*dx1 + dx2) / (dx1 * (dx1 + dx2))
    b =  (dx1 + dx2) / (dx1 * dx2)
    c = - dx1 / (dx2 * (dx1 + dx2))

    fout_L = a*f(1) + b*f(2) + c*f(3)

    ! ---- Right boundary ----
    dx1 = dx(size(f)-1)
    dx2 = dx(size(f))

    a =  dx2 / (dx1 * (dx1 + dx2))
    b = - (dx1 + dx2) / (dx1 * dx2)
    c = (2.d0*dx2 + dx1) / (dx2 * (dx1 + dx2))

    fout_R = a*f(size(f)-2) + b*f(size(f)-1) + c*f(size(f))
  end subroutine get_boundary
  subroutine get_dphi(u, r, delta, du)
    real(8), intent(in) :: u(:)
    integer :: idx
    real(8), intent(out) :: du(:)
    real(8), intent(in) :: delta, r
    real(8) :: hs, hd
    do idx = 2, size(u)
      if (idx .eq. size(u)) then
        du(idx) = (u(1) - u(idx-1))/(2.d0*delta)
      else
        du(idx) = (u(idx+1) - u(idx-1))/(2.d0*delta)
      end if
    end do
    du(1) = (u(2)-u(size(u)))/(2.d0*delta)
    du = du/r
  end subroutine get_dphi
end module centralfd
