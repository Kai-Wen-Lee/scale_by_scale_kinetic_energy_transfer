PROGRAM sbs_keflux
  ! .......... .......... Declarations  .......... .......... !
  USE MPI
  USE mpi_context
  USE domain_decomposition
  USE ncinput
  USE ncoutput
  USE interpolate
  USE conv3d
  USE centralfd
  implicit none
  logical :: exist
  character(len=256) :: geom_fname, flow_fname
  real(8), allocatable :: ur_interp(:,:,:), uz_interp(:,:,:)
  real(8), allocatable :: ur_f(:,:,:), uphi_f(:,:,:), uz_f(:,:,:)
  real(8), parameter:: k = 30.d0
  real(8), dimension(:,:,:), allocatable :: u11,u12,u13,u21,u22,u23,u31,u32,u33
  real(8), dimension(:,:,:), allocatable :: u1u1, u1u2, u1u3, u2u1, u2u2, u2u3, u3u1, u3u2, u3u3
  real(8), dimension(:,:,:), allocatable :: P_F
  real(8), dimension(:,:,:), allocatable :: d1u1,d1u2,d1u3,d2u1,d2u2,d2u3,d3u1,d3u2,d3u3
  real(8), dimension(:), allocatable :: P_F_a_av
  real(8), allocatable :: ur_p_flat(:), uz_p_flat(:)
  real(8) :: P_F_v_av, P_F_, A_zi
  integer :: idxiii, idxjjj, idxkkk

  ! .......... .......... MPI INIT .......... .......... !
  call init_mpi()
  if (proc_id == 0) then
    print *, 'mpi init'
  end if
  ! Inner product for resulting tensor to get:
    ! P_F
  ! Outputs: 
    ! One .nc for each n containing:
      ! P_F
      ! Filtered velocity Uz
      ! Filtered velocity Uphi
      ! Filtered velocity Ur
  ! .......... .......... Read geometry data .......... .......... !
  geom_fname = 'geometry.nc'
  inquire(file=geom_fname, exist=exist)
  if (exist) call readgeom(geom_fname)
  ! Read flow data
  flow_fname = 'flow0340.nc'
  inquire(file=flow_fname, exist=exist)
  if (exist) call readflow(flow_fname)

  if (proc_id == 0) then
    print *, 'geom  and dataloaded'
  end if
  ! Domain decomposition for mpi
  call split_domain()
  if (proc_id == 0) then
    print *, 'domain splitted'
  end if
  !verify
  !print *, 'Process', proc_id, ' of ', num_procs&
  !,'z0', mpi_zidx_sta, 'zn', mpi_zidx_end&
  !, 'r0', mpi_ridx_sta, 'rn', mpi_ridx_end&
  !, 'phi0', mpi_phiidx_sta, 'phin', mpi_phiidx_end
  ! .......... .......... Interpolate quantities to pmesh .......... ..........
  call interp_to_pmesh()
  if (proc_id == 0) then
    print *, 'interp done'
  end if
  ! .......... .......... Sync .......... .......... !
  allocate(ur_interp(dim_phi, dim_r, dim_z))
  allocate(uz_interp(dim_phi, dim_r, dim_z))

  allocate(ur_p_flat(dim_r*dim_phi*mpi_dim_z))
  allocate(uz_p_flat(dim_r*dim_phi*mpi_dim_z))
  ! check reshape ordering, might be wrong!
  call flatten(ur, ur_p_flat)
  call flatten(uz, uz_p_flat)
  call sync_z(ur_p_flat, dim_r*dim_phi*mpi_dim_z, dim_r*dim_phi*dim_z, ur_interp)
  call sync_z(uz_p_flat, dim_r*dim_phi*mpi_dim_z, dim_r*dim_phi*dim_z, uz_interp)
  call mpi_bcast(ur_interp, size(ur_interp), mpi_double_precision, 0, comm, ierr)
  call mpi_bcast(uz_interp, size(uz_interp), mpi_double_precision, 0, comm, ierr)
  if (proc_id == 0) then
    print *, 'sync '
  end if

  ! Main loop
  ! .......... .......... Convolution .......... ..........!
  if (proc_id == 0) then
    print *, 'conv ur started'
  end if
  call convolution(ur_interp, k, ur_f)
  call mpi_bcast(ur_f, size(ur_f), mpi_double_precision, 0, comm, ierr)
  if (proc_id == 0) then
    print *, 'conv uz started'
  end if
  call convolution(uz_interp, k, uz_f)
  call mpi_bcast(uz_f, size(uz_f), mpi_double_precision, 0, comm, ierr)
  if (proc_id == 0) then
    print *, 'conv uphi started'
  end if
  call convolution(uphi, k, uphi_f)
  call mpi_bcast(uphi_f, size(uphi_f), mpi_double_precision, 0, comm, ierr)
  ! Get filtered Tres
  allocate(u11(dim_phi, dim_r, dim_z))
  allocate(u12(dim_phi, dim_r, dim_z))
  allocate(u13(dim_phi, dim_r, dim_z))

  allocate(u21(dim_phi, dim_r, dim_z))
  allocate(u22(dim_phi, dim_r, dim_z))
  allocate(u23(dim_phi, dim_r, dim_z))

  allocate(u31(dim_phi, dim_r, dim_z))
  allocate(u32(dim_phi, dim_r, dim_z))
  allocate(u33(dim_phi, dim_r, dim_z))

  allocate(u1u1(dim_phi, dim_r, dim_z))
  allocate(u1u2(dim_phi, dim_r, dim_z))
  allocate(u1u3(dim_phi, dim_r, dim_z))

  allocate(u2u1(dim_phi, dim_r, dim_z))
  allocate(u2u2(dim_phi, dim_r, dim_z))
  allocate(u2u3(dim_phi, dim_r, dim_z))

  allocate(u3u1(dim_phi, dim_r, dim_z))
  allocate(u3u2(dim_phi, dim_r, dim_z))
  allocate(u3u3(dim_phi, dim_r, dim_z))

  if (proc_id == 0) then
    print *, 'Calculating Tres'
  end if
  call convolution(ur_interp*ur_interp, k, u11)
  call mpi_bcast(u11, size(u11), mpi_double_precision, 0, comm, ierr)
  call convolution(ur_interp*uphi, k, u12)
  call mpi_bcast(u12, size(u12), mpi_double_precision, 0, comm, ierr)
  call convolution(ur_interp*uz_interp, k, u13)
  call mpi_bcast(u13, size(u13), mpi_double_precision, 0, comm, ierr)
  !call convolution(uphi*ur_interp, k, u21)
  u21 = u12
  call convolution(uphi*uphi, k, u22)
  call mpi_bcast(u22, size(u22), mpi_double_precision, 0, comm, ierr)
  call convolution(uphi*uz_interp, k, u23)
  call mpi_bcast(u23, size(u23), mpi_double_precision, 0, comm, ierr)
  !call convolution(uz_interp*ur_interp, k, u31)
  !call convolution(uz_interp*uphi, k, u32)
  u31 = u13
  u32 = u23
  call convolution(uz_interp*uz_interp, k, u33)
  call mpi_bcast(u33, size(u33), mpi_double_precision, 0, comm, ierr)
  u1u1 = ur_f*ur_f
  u1u2 = ur_f*uphi_f
  u1u3 = ur_f*uz_f

  u2u1 = uphi_f*ur_f
  u2u2 = uphi_f*uphi_f
  u2u3 = uphi_f*uz_f

  u3u1 = uz_f*ur_f
  u3u2 = uz_f*uphi_f
  u3u3 = uz_f*uz_f

  if (proc_id == 0) then
    print *, 'Calculating Sijk'
  end if

  ! Get filtered Sijk
  allocate(d1u1(dim_phi, dim_r, dim_z))
  allocate(d1u2(dim_phi, dim_r, dim_z))
  allocate(d1u3(dim_phi, dim_r, dim_z))

  allocate(d2u1(dim_phi, dim_r, dim_z))
  allocate(d2u2(dim_phi, dim_r, dim_z))
  allocate(d2u3(dim_phi, dim_r, dim_z))

  allocate(d3u1(dim_phi, dim_r, dim_z))
  allocate(d3u2(dim_phi, dim_r, dim_z))
  allocate(d3u3(dim_phi, dim_r, dim_z))

  call gradient(ur_f, r_p, dr_ur, dphi, dz_uz, d1u1, d2u1, d3u1)
  call gradient(uphi_f, r_p, dr_ur, dphi, dz_uz, d1u2, d2u2, d3u2)
  call gradient(uz_f, r_p, dr_ur, dphi, dz_uz, d1u3, d2u3, d3u3)

  ! Inner product
  if (proc_id == 0) then
    print *, 'Calculating Inner product'
  end if
  allocate(P_F(dim_phi, dim_r, dim_z))
  P_F(:,:,:) =&
  (u11-u1u1)*0.5d0*(d1u1+d1u1)+&
  (u12-u1u2)*0.5d0*(d1u2+d2u1)+&
  (u13-u1u3)*0.5d0*(d1u3+d3u1)+&
  (u21-u2u1)*0.5d0*(d2u1+d1u2)+&
  (u22-u2u2)*0.5d0*(d2u2+d2u2)+&
  (u23-u2u3)*0.5d0*(d2u3+d3u2)+&
  (u31-u3u1)*0.5d0*(d3u1+d1u3)+&
  (u32-u3u2)*0.5d0*(d3u2+d2u3)+&
  (u33-u3u3)*0.5d0*(d3u3+d3u3)

  allocate(P_F_a_av(dim_z))
  ! Area average
  do idxkkk = 1, dim_z
    P_F_=0.d0
    A_zi = 0.d0
    do idxiii = 1, dim_r
      do idxjjj = 1, dim_phi
        P_F_a_av(idxkkk) = P_F_+P_F(idxjjj, idxiii, idxkkk)*r_p(idxiii)*dphi*dr_p(idxiii)
        A_zi = A_zi + r_p(idxiii)*dphi*dr_p(idxiii)
      end do
    end do
    P_F_a_av(idxkkk) = P_F_a_av(idxkkk)/A_zi
  end do

  P_F_v_av=0.d0
  do idxkkk = 1,dim_z
    P_F_v_av=P_F_v_av+P_F_a_av(idxkkk)*dz_p(idxkkk)
  end do
  P_F_v_av=P_F_v_av/maxval(z_p)

  ! Save
  if (proc_id == 0) then
    print *, 'Saving'
  end if
  call writePF('KET_flow0346_k30.nc', P_F, P_F_a_av, P_F_v_av, k, ur_f, uphi_f, uz_f)
  call finalize_mpi()
END PROGRAM

