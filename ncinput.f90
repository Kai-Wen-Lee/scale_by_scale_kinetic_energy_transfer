module ncinput
  USE mpi
  USE mpi_context
  USE netcdf
  implicit none
  ! Geometry
  integer :: dim_r, dim_phi, dim_z
  real(8), allocatable :: r_p(:), r_ur(:), dr_p(:), dr_ur(:), r_mid(:)
  real(8), allocatable :: z_p(:), z_uz(:), dz_p(:), dz_uz(:)
  real(8), allocatable :: phi(:)
  real(8) :: dphi
  ! Flow
  real(8), allocatable :: t(:,:,:), ur(:,:,:), uphi(:,:,:), uz(:,:,:)

  contains
  ! :::::::::: :::::::::: :::::::::: :::::::::: Read Geom ::::::::::: :::::::::: :::::::::: :::::::::: ::::::::::!
  subroutine readgeom(file_name)
    USE netcdf
    implicit none
    character(len=*), intent(in) :: file_name
    ! Dimensions
    ! ID
    integer :: ncid
    integer :: varid
    integer :: dimid
    ! dummy variables
    real(8) :: phi_
    integer :: i
    ! .......... .......... Open file with mpiposix .......... ...........!
    call check(nf90_open(file_name, ior(nf90_nowrite,nf90_mpiposix), ncid, comm = comm, info = mpi_info_null))

    ! .......... .......... Get dimensions .......... ..........!
    call check(nf90_inq_dimid(ncid, "dim_r", dimid))
    call check(nf90_inquire_dimension(ncid, dimid, len=dim_r))
    call check(nf90_inq_dimid(ncid, "dim_phi", dimid))
    call check(nf90_inquire_dimension(ncid, dimid, len=dim_phi))
    call check(nf90_inq_dimid(ncid, "dim_z", dimid))
    call check(nf90_inquire_dimension(ncid, dimid, len=dim_z))

    ! .......... .......... Allocate arrays .......... .......... !
    allocate(r_p(dim_r))
    allocate(r_ur(dim_r))
    allocate(r_mid(dim_r))
    allocate(dr_p(dim_r))
    allocate(dr_ur(dim_r))

    allocate(z_p(dim_z))
    allocate(z_uz(dim_z))
    allocate(dz_p(dim_z))
    allocate(dz_uz(dim_z))

    allocate(phi(dim_phi-4))

    ! .......... .......... Get geom vars .......... ..........!
    call check(nf90_inq_varid(ncid, 'r_p', varid))
    call check(nf90_get_var(ncid, varid, r_p))

    call check(nf90_inq_varid(ncid, 'r_ur', varid))
    call check(nf90_get_var(ncid, varid, r_ur))

    call check(nf90_inq_varid(ncid, 'r_mid', varid))
    call check(nf90_get_var(ncid, varid, r_mid))

    call check(nf90_inq_varid(ncid, 'delta_r_p', varid))
    call check(nf90_get_var(ncid, varid, dr_p))

    call check(nf90_inq_varid(ncid, 'delta_r_ur', varid))
    call check(nf90_get_var(ncid, varid, dr_ur))

    call check(nf90_inq_varid(ncid, 'z_p', varid))
    call check(nf90_get_var(ncid, varid, z_p))

    call check(nf90_inq_varid(ncid, 'z_uz', varid))
    call check(nf90_get_var(ncid, varid, z_uz))

    call check(nf90_inq_varid(ncid, 'delta_z_uz', varid))
    call check(nf90_get_var(ncid, varid, dz_uz))

    call check(nf90_inq_varid(ncid, 'delta_z_p', varid))
    call check(nf90_get_var(ncid, varid, dz_p))

    call check(nf90_inq_varid(ncid, 'delta_phi', varid))
    call check(nf90_get_var(ncid, varid, dphi))

    call check(nf90_close(ncid))


    ! remove ghost points
    dim_r = dim_r - 2
    dim_z = dim_z - 2
    dim_phi = dim_phi - 4

    r_p = r_p(2:dim_r+1)
    r_ur = r_ur(2:dim_r+1)
    r_mid = r_mid(2:dim_r+1)
    dr_p = dr_p(2:dim_r+1)
    dr_ur = dr_ur(2:dim_r+1)
    
    z_p = z_p(2:dim_z+1)
    z_uz = z_uz(2:dim_z+1)
    dz_p = dz_p(2:dim_z+1)
    dz_uz = dz_uz(2:dim_z+1)

    ! compute derived quantities
    phi_= 0.d0
    phi(1) = 0.d0
    do i = 2, dim_phi
      phi(i) = phi_ + dphi
      phi_ = phi(i)
    end do

  end subroutine readgeom
! :::::::::: :::::::::: :::::::::: :::::::::: ::::::::::: :::::::::: :::::::::: :::::::::: ::::::::::!
  subroutine readflow(file_name)
    USE netcdf
    implicit none
    ! IDs
    integer :: ncid
    integer :: tid, uzid, urid, uphiid
    ! reshape targer array
    real(8), allocatable :: perm(:,:,:), t_(:,:,:), ur_(:,:,:), uphi_(:,:,:), uz_(:,:,:)
    character(len=256), intent(in) :: file_name

    ! .......... .......... Open file with mpiposix .......... ...........!
    call check(nf90_open(file_name, ior(nf90_nowrite,nf90_mpiposix), ncid, comm = comm, info = mpi_info_null))

    ! .......... .......... Allocate arrays .......... .......... !
    allocate(t_(dim_phi,dim_r, dim_z))
    allocate(ur_(dim_phi,dim_r, dim_z))
    allocate(uphi_(dim_phi,dim_r, dim_z))
    allocate(uz_(dim_phi,dim_r, dim_z))

    ! .......... .......... Get arrays .......... .......... !
    call check(nf90_inq_varid(ncid, 'temp', tid))
    call check(nf90_get_var(ncid, tid, t_, start = [3,2,2], count = [dim_phi, dim_r, dim_z]))

    call check(nf90_inq_varid(ncid, 'uz', uzid))
    call check(nf90_get_var(ncid, uzid, uz_, start = [3,2,2], count = [dim_phi, dim_r, dim_z]))

    call check(nf90_inq_varid(ncid, 'ur', urid))
    call check(nf90_get_var(ncid, urid, ur_, start = [3,2,2], count = [dim_phi, dim_r, dim_z]))

    call check(nf90_inq_varid(ncid, 'uphi', uphiid))
    call check(nf90_get_var(ncid, uphiid, uphi_, start = [3,2,2], count = [dim_phi, dim_r, dim_z]))

    call check(nf90_close(ncid))
    ! .......... .......... Reshape .......... .......... !
    allocate(perm(dim_r+2, dim_phi+4, dim_z+2))
    !allocate(t(dim_phi, dim_r, dim_z))
    !allocate(ur(dim_phi, dim_r, dim_z))
    !allocate(uphi(dim_phi, dim_r, dim_z))
    !allocate(uz(dim_phi, dim_r, dim_z))

    !t_ = reshape(t_, shape(perm))
    !ur_ = reshape(ur_, shape(perm))
    !uphi_ = reshape(uphi_, shape(perm))
    !uz_ = reshape(uz_, shape(perm))

    !t = t_(2:dim_r+1, 3:dim_phi+2, 2:dim_z+1)
    !ur = ur_(2:dim_r+1, 3:dim_phi+2, 2:dim_z+1)
    !uphi = uphi_(2:dim_r+1, 3:dim_phi+2, 2:dim_z+1)
    !uz = uz_(2:dim_r+1, 3:dim_phi+2, 2:dim_z+1)
    t = t_
    ur = ur_
    uphi = uphi_
    uz = uz_

  end subroutine readflow
  subroutine check(status)
    USE netcdf
    implicit none
    integer, intent (in):: status
      if(status .ne. nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop 
      end if
  end subroutine check  
end module ncinput
