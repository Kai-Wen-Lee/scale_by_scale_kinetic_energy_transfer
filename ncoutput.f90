module ncoutput
  USE mpi
  USE netcdf
  USE ncinput
  USE mpi_context
  contains
  subroutine writePF(file_name,P_F, PF_A_avg, PF_V_avg, k, ur_f, uphi_f, uz_f)
    implicit none
    character(len=*), intent(in) :: file_name
    real(8), intent(in) :: P_F(:,:,:), PF_A_avg(:), PF_V_avg, k
    real(8), intent(in) :: ur_f(:, :, :), uphi_f(:, :, :), uz_f(:, :, :)
    integer :: ncid
    integer :: urf_id, uphif_id, uzf_id
    integer :: PF_id, PF_A_avg_id, PF_V_avg_id, kid
    integer :: r_dimid, phi_dimid, z_dimid
    integer :: drpid, dphiid, zpid, dzpid, rpid, rurid, rmidid
    integer :: dimids(3)
    integer :: cmode

    cmode = ior(NF90_HDF5, NF90_MPIPOSIX)
    cmode = ior(cmode, NF90_CLOBBER)

    call check(nf90_create(file_name,cmode=cmode, ncid=ncid,comm=comm, info=MPI_INFO_NULL))
                           
  ! Define the dimensions (inluding boundaries). NetCDF will hand back an ID for
  ! each. Metadata operations must take place on all processors.
    call check(nf90_def_dim(ncid, 'dim_r',   dim_r,   r_dimid))
    call check(nf90_def_dim(ncid, 'dim_phi', dim_phi, phi_dimid))
    call check(nf90_def_dim(ncid, 'dim_z',   dim_z,   z_dimid))

    dimids = [phi_dimid, r_dimid, z_dimid]
    
  ! Define the variables.
    call check(nf90_def_var(ncid, 'r_p',   nf90_double, r_dimid, rpid))
    call check(nf90_def_var(ncid, 'r_ur',  nf90_double, r_dimid, rurid))
    call check(nf90_def_var(ncid, 'r_mid', nf90_double, r_dimid, rmidid))
    call check(nf90_def_var(ncid, 'delta_r_p',  nf90_double, r_dimid, drpid))
    call check(nf90_def_var(ncid, 'delta_phi', nf90_double, dphiid))
    call check(nf90_def_var(ncid, 'z_p',  nf90_double, z_dimid, zpid))
    call check(nf90_def_var(ncid, 'delta_z_p',  nf90_double, z_dimid, dzpid))


    call check(nf90_def_var(ncid, 'PF', nf90_double, dimids, PF_id))
    call check(nf90_def_var(ncid, 'ur_f', nf90_double, dimids, urf_id))
    call check(nf90_def_var(ncid, 'uphi_f', nf90_double, dimids, uphif_id))
    call check(nf90_def_var(ncid, 'uz_f', nf90_double, dimids, uzf_id))
    call check(nf90_def_var(ncid, 'PF_A_avg', nf90_double, dimids(3), PF_A_avg_id))
    call check(nf90_def_var(ncid, 'PF_V_avg', nf90_double, PF_V_avg_id))
    call check(nf90_def_var(ncid, 'k', nf90_double, kid))
    call check(nf90_enddef(ncid))

    !only let 1st processor write
    if (proc_id .eq. 0) then

      call check(nf90_put_var(ncid, rpid, r_p))
      call check(nf90_put_var(ncid, rurid, r_ur))
      call check(nf90_put_var(ncid, rmidid, r_mid))
      call check(nf90_put_var(ncid, drpid, dr_p))

      call check(nf90_put_var(ncid, dphiid, dphi))

      call check(nf90_put_var(ncid, zpid, z_p))
      call check(nf90_put_var(ncid, dzpid, dz_p))

      call check(nf90_put_var(ncid, kid, k))

      call check(nf90_put_var(ncid, PF_id, P_F, start=[1,1,1], count=[dim_phi, dim_r, dim_z]))
      call check(nf90_put_var(ncid, urf_id, ur_f, start=[1,1,1], count=[dim_phi, dim_r, dim_z]))
      call check(nf90_put_var(ncid, uphif_id, uphi_f, start=[1,1,1], count=[dim_phi, dim_r, dim_z]))
      call check(nf90_put_var(ncid, uzf_id, uz_f, start=[1,1,1], count=[dim_phi, dim_r, dim_z]))
      call check(nf90_put_var(ncid, PF_A_avg_id, PF_A_avg))
      call check(nf90_put_var(ncid, PF_V_avg_id, PF_V_avg))

    end if
    call check(nf90_close(ncid))
  end subroutine
end module ncoutput
