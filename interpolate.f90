module interpolate
  USE io
  USE domain_decomposition
  USE mpi_context
  implicit none
  real(8), allocatable :: ur_p(:,:,:), uz_p(:,:,:)
  contains
  subroutine interp_to_pmesh()
    implicit none
    integer :: i, j, k
    allocate(ur_p(dim_phi, dim_r, dim_z))
    !allocate(uphi_p(mpi_dim_r, mpi_dim_phi, mpi_dim_z))
    allocate(uz_p(dim_phi, dim_r, dim_z))
    ! interp for z
    do i = 1, dim_r
      do j = 1, dim_phi
        call oneD_cubic_interp(z_uz(:),&
        uz(j,i,:),&
        z_p(:),&
        uz_p(j,i,:),dim_z)
        ! oneD_cubic_interp(inputx, inputval, targetx, len(x))
      end do
    end do

    !interp for ur
    do j = 1, dim_phi
      do k = 1, dim_z
        call oneD_cubic_interp(r_ur(:),&
        ur(j,:,k),&
        r_p(:),&
        ur_p(j,:,k),dim_r)
        ! oneD_cubic_interp(inputx, inputval, targetx, len(x))
      end do
    end do
  end subroutine interp_to_pmesh

  subroutine oneD_cubic_interp(k, k_vals, k_target, kout, dim_k)
    implicit none
    integer, intent(in) :: dim_k
    real(8), dimension(dim_k), intent(in) :: k(dim_k), k_vals(dim_k), k_target(dim_k)
    real(8), dimension(dim_k) ::b(dim_k), c(dim_k), d(dim_k)
    real(8), dimension(dim_k), intent(out) :: kout(dim_k)
    integer :: i

    call spline(k, k_vals, b, c, d, dim_k)
    do i= 1,dim_k
      call ispline(kout(i), k_target(i), k, k_vals, b, c, d, dim_k)
    end do
  end subroutine oneD_cubic_interp
  subroutine spline (x, y, b, c, d, n)
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !----------------------------------------------------------------------
  !  input..
  !  x = the arrays of data abscissas (in strictly increasing order)
  !  y = the arrays of data ordinates
  !  n = size of the arrays xi() and yi() (n>=2)
  !  output..
  !  b, c, d  = arrays of spline coefficients
  !  comments ...
  !  spline.f90 program is based on fortran version of program spline.f
  !  the accompanying function fspline can be used for interpolation
  !======================================================================
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: x(n), y(n)
    real(8), dimension(n), intent(out) :: b(n), c(n), d(n)
    integer :: i, j, gap
    real(8) :: h

    gap = n-1
    ! check input
    if ( n < 2 ) return
    if ( n < 3 ) then
      b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
    end if
    !
    ! step 1: preparation
    !
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do
    !
    ! step 2: end conditions 
    !
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0
    c(n) = 0.0
    if(n /= 3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if
    !
    ! step 3: forward elimination 
    !
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do
    !
    ! step 4: back substitution
    !
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    !
    ! step 5: compute spline coefficients
    !
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)
  end subroutine spline
  subroutine ispline(u_out, u, x, y, b, c, d, n)
  !======================================================================
  ! ispline evaluates the cubic spline interpolation at point z
  ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  !----------------------------------------------------------------------
  ! input..
  ! u       = the abscissa at which the spline is to be evaluated
  ! x, y    = the arrays of given data points
  ! b, c, d = arrays of spline coefficients computed by spline
  ! n       = the number of data points
  ! output:
  ! uout = interpolated value at point u
  !=======================================================================
    implicit none
    real(8), intent(out) :: u_out
    integer, intent(in) :: n
    real(8), intent(in) :: u
    real(8), dimension(n), intent(in) :: x(n), y(n), b(n), c(n), d(n)
    integer :: i,j,k
    real(8) :: dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u <= x(1)) then
      u_out = y(1)
      return
    end if
    if(u >= x(n)) then
      u_out = y(n)
      return
    end if
    !*
    !  binary search for for i, such that x(i) <= u <= x(i+1)
    !*
    i = 1
    j = n+1
    do while (j > i+1)
      k = (i+j)/2
      if(u < x(k)) then
        j=k
        else
        i=k
       end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = u - x(i)
    u_out = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  end subroutine ispline
end module interpolate

