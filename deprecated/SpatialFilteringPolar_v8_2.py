# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import griddata,RegularGridInterpolator

import netCDF4

from pathlib import Path
import gc
import time

from mpi4py import MPI
import conv2D
# conv2D is the .so resulting from MPIFORT compilation of the
# 2D Convolution .f90 subroutine
# Called as:
# filtered, raw = conv2D.convolution2dpolar(data_xy, rr, phi, dr3, dphi3, delta, fcomm)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# INITIALISED MPI
fcomm = MPI.COMM_WORLD.py2f()
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def convert_to_array(maskedarr):
    computed_array = maskedarr.compute()
    tempndarr = computed_array.data
    return da.from_array(tempndarr)
# NEED TO TREAT R=0 AT SOME POINT !!
def getGeom(geomDir):
    geom = netCDF4.Dataset(geomDir, 'r')
    geom.set_auto_mask(False) 
    zu = np.asarray(geom['z_uz'][1:-1])
    zp = np.asarray(geom['z_p'][1:-1])
    dz = np.asarray(geom['delta_z_p'][1:-1])
    rp = np.asarray(geom['r_p'][1:-1])
    rw = np.asarray(geom['r_ur'][1:-1])
    dr = np.asarray(geom['delta_r_p'][1:-1])
    dphi = np.asarray(geom['delta_phi'][:])
    dimz = len(geom.dimensions['dim_z'])-2
    dimphi = len(geom.dimensions['dim_phi'])-4
    dimr = len(geom.dimensions['dim_r'])-2
    height = np.max(zp)
    phi=np.linspace(-1*np.pi, np.pi, dimphi)
    geom.close()
    return [zu, zp, dz, rp, rw, dr, dphi, dimz, dimphi, dimr, height, phi]

def getFlow(flowDir):
    flow = netCDF4.Dataset(flowDir, 'r')
    flow.set_auto_mask(False) 
    temp = np.asarray(flow['temp'][1:-1,1:-1,2:-2])  # Assuming (z,r,phi)
    uz = np.asarray(flow['uz'][1:-1,1:-1,2:-2])
    uphi = np.asarray(flow['uphi'][1:-1,1:-1,2:-2])
    ur = np.asarray(flow['ur'][1:-1,1:-1,2:-2])
    flow.close()
    return [temp, uz, uphi, ur]

def save_variable_arrays(output_dir: Path, case_id: str, data_dict: dict):
    for key, value in data_dict.items():
        try:
            np.save(output_dir / f"{key}_{case_id}.npy", value)
        except Exception:
            #np.save(output_dir / f"{key}_{case_id}.npy", value)
            np.save(output_dir / f"{key}_{case_id}.npy", value)
            
def load_case_results(input_dir: Path, case_id: str, variable_names: list):
    data = {}
    for var in variable_names:
        file_path = input_dir / f"{var}_{case_id}.npy"
        data[var] = np.load(file_path)
    return data

def areaAv(timeAvgField, geom):
    zu, zp, dz, rp, rw, dr, dphi, dimz, dimphi, dimr, height, phi = geom
    # Compute area weights per (r, phi) cella
    area_i = rp[:] * dr[:] * dphi   # shape (len(rp)-1,)
    area_slice = np.broadcast_to(area_i[:, np.newaxis], (len(rp), len(phi)))
    area_weights = area_slice[:, :]  # (Nr, Nphi-2, 1)
    # Total area
    total_area_per_plane = area_slice.sum()
    
    # Area-weighted average over (Nr, Nphi) for each time/z
    out = (timeAvgField * area_weights).sum(axis=(0,1)) / total_area_per_plane #sum R then sum Phi
    return out

def lineAv_cartesian(PIAreaAvField, z):
    out = np.sum(PIAreaAvField * z) / np.sum(z)
    return out
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def process_case(case_id: str, base_dir: Path, output_dir: Path, k_, fileidx):
    comm=MPI.COMM_WORLD
    if comm.Get_rank() == 0:
        print(f'{case_id} file {fileidx} at k= {k_} started', flush=True)
    start=time.time()
    # ---------- ---------- ---------- Load data ---------- ---------- ----------#
    # NO LOCK AS NOT PARALLELIZED IN CASES ANYMORE
    geom = getGeom(base_dir / 'geometry.nc')
    flow_file = sorted(base_dir.glob("flow*.nc"))[fileidx]
    flowDir = base_dir / flow_file
    flow = getFlow(flowDir)  # This reads from file(s)   

    # 0 , 1  , 2  , 3  , 4  , 5  , 6    , 7    , 8      , 9    , 10     , 11
    #zu_, zp_, dz_, rp_, rw_, dr_, dphi_, dimz_, dimphi_, dimr_, height_, phi_ = geom
    # ---------- ---------- ---------- Interpolate ---------- ---------- ----------#
    Rp, Phi, Zp = np.meshgrid(geom[3], geom[-1], geom[1] , indexing='ij')  # shape: (len(z), 149)
    Uz_interp_rpz = RegularGridInterpolator((geom[3],geom[-1],geom[0]),np.asarray(flow[1].transpose(1,2,0)), method='linear', bounds_error=False, fill_value=None)
    Ur_interp_rpz = RegularGridInterpolator((geom[4],geom[-1],geom[1]),np.asarray(flow[3].transpose(1,2,0)), method='linear', bounds_error=False, fill_value=None)
    Uphi_interp_rpz = RegularGridInterpolator((geom[3],geom[-1],geom[1]),np.asarray(flow[2].transpose(1,2,0)), method='linear', bounds_error=False, fill_value=None)
    def run_interp(interp, pts):
        return interp(pts)
    Uz = Uz_interp_rpz((Rp, Phi, Zp))
    Ur = Ur_interp_rpz((Rp, Phi, Zp))
    Uphi = Uphi_interp_rpz((Rp, Phi, Zp))
    
#    print(f'{case_id} file {fileidx} at k= {k_}: \t Interpolated to P-mesh ::: \t took {end - start:.2f} seconds', flush=True)

    # ---------- ---------- ---------- Filter velocities ---------- ---------- ----------#
    Ur_filtered = np.zeros(Ur.shape)
    Uphi_filtered = np.zeros(Uphi.shape)
    Uz_filtered = np.zeros(Uz.shape)
    
    for z_index in range(len(geom[1])):
        Delta = np.pi/k_
        Ur_filtered[:,:,z_index] = conv2D.convolution2dpolar(Ur[:, :, z_index], geom[3], geom[-1], geom[5], geom[6],Delta, fcomm)
        Uphi_filtered[:,:,z_index] = conv2D.convolution2dpolar(Uphi[:, :, z_index], geom[3],geom[-1], geom[5], geom[6],Delta, fcomm)
        Uz_filtered[:,:,z_index] = conv2D.convolution2dpolar(Uz[:, :, z_index], geom[3],geom[-1], geom[5], geom[6],Delta, fcomm)

    comm.Barrier()
    comm.Bcast(Ur_filtered, root=0)
    comm.Bcast(Uphi_filtered, root=0)
    comm.Bcast(Uz_filtered, root=0)
    gc.collect()
#    print(f'{case_id} file {fileidx} at k= {k_}: \t Ur, Uy, Uz filter done ::: \t took {end - start:.2f} seconds', flush=True)

    # ---------- ---------- ---------- Calculate Rate of Strain ---------- ---------- ----------#
    def GetS(Ur_f,Uphi_f,Uz_f,offset=4.0): 
        # Gradients
        r_ = np.asarray(geom[3])
        phi_ = np.asarray(geom[-1])
        z_ = np.asarray(geom[1])
        d1u1, d2u1, d3u1 = np.gradient(np.asarray(Ur_f), r_+offset, phi_, z_,edge_order=2)
        d1u2, d2u2, d3u2 = np.gradient(np.asarray(Uphi_f),  r_+offset, phi_, z_,edge_order=2)
        d1u3, d2u3, d3u3 = np.gradient(np.asarray(Uz_f),  r_+offset, phi_, z_,edge_order=2)
        R2, Phi2, Z2 = np.meshgrid(r_, phi_, z_, indexing = 'ij')
        ''' Anti-symmetric part of the Jacobian (rotation)
        S = np.array([
            [np.zeros(Rp.shape), 0.5*(d1u2-(1/Rp)*d2u1), 0.5*(d1u3-d3u1)],
            [0.5*(1/Rp*d2u1-d1u2), np.zeros(Rp.shape), 0.5*((1/Rp)*d2u3-d3u2)],
            [0.5*(d3u1-d1u3), 0.5*(d3u2-(1/Rp)*d2u3), np.zeros(Rp.shape)]])
        '''
        # Symmetic part of the Jacobian (strain rate)
        S = np.array([
        [d1u1, 0.5*(d1u2+(1/R2)*d2u1), 0.5*(d1u3+d3u1)],
        [0.5*(1/R2*d2u1+d1u2), (1/R2)*d2u2, 0.5*((1/R2)*d2u3+d3u2)],
        [0.5*(d3u1+d1u3), 0.5*(d3u2+(1/R2)*d2u3), d3u3]])
        del(d1u1,d2u1,d3u1,d1u2,d2u2,d3u2,d1u3,d2u3,d3u3)
        return S

    S = GetS(Ur_filtered,Uphi_filtered,Uz_filtered)
    gc.collect()
    comm.Barrier()
#    print(f'{case_id} file {fileidx} at k= {k_}: \t Sij done ::: \t took {end - start:.2f} seconds', flush=True)

    # ---------- ---------- ---------- Calculate Residual stress tensor ---------- ---------- ----------#
    Tres=np.zeros([3,3,Ur.shape[0],Ur.shape[1], Ur.shape[2]])

    def GetTres(z_index):
        Delta = np.pi/k_

        u11 = conv2D.convolution2dpolar(Ur[:,:,z_index]*Ur[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta, fcomm)
        u12 = conv2D.convolution2dpolar(Ur[:,:,z_index]*Uphi[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta, fcomm)
        u13 = conv2D.convolution2dpolar(Ur[:,:,z_index]*Uz[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta, fcomm)
        
        u21 = u12 #convolution2DPolar(Uphi[:,:,z_index]*Ur[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta)
        u22 = conv2D.convolution2dpolar(Uphi[:,:,z_index]*Uphi[:,:,z_index], geom[3],geom[-1],geom[5], geom[6],  Delta, fcomm)
        u23 = conv2D.convolution2dpolar(Uphi[:,:,z_index]*Uz[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta, fcomm)

        u31 = u13 #convolution2DPolar(Uz[:,:,z_index]*Ur[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta)
        u32 = u23 #convolution2DPolar(Uz[:,:,z_index]*Uphi[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta)
        u33 = conv2D.convolution2dpolar(Uz[:,:,z_index]*Uz[:,:,z_index], geom[3],geom[-1], geom[5], geom[6], Delta, fcomm)

        UUf = np.array([
            [u11,u12,u13],
            [u21,u22,u23],
            [u31,u32,u33]])

        u1u1 = Ur_filtered[:,:,z_index]*Ur_filtered[:,:,z_index]
        u1u2 = Ur_filtered[:,:,z_index]*Uphi_filtered[:,:,z_index]
        u1u3 = Ur_filtered[:,:,z_index]*Uz_filtered[:,:,z_index]

        u2u1 = Uphi_filtered[:,:,z_index]*Ur_filtered[:,:,z_index]
        u2u2 = Uphi_filtered[:,:,z_index]*Uphi_filtered[:,:,z_index]
        u2u3 = Uphi_filtered[:,:,z_index]*Uz_filtered[:,:,z_index]

        u3u1 = Uz_filtered[:,:,z_index]*Ur_filtered[:,:,z_index]
        u3u2 = Uz_filtered[:,:,z_index]*Uphi_filtered[:,:,z_index]
        u3u3 = Uz_filtered[:,:,z_index]*Uz_filtered[:,:,z_index]

        UfUf = np.array([
            [u1u1,u1u2,u1u3],
            [u2u1,u2u2,u2u3],
            [u3u1,u3u2,u3u3]])

        Tres_slice = UUf-UfUf
        del(u11,u12,u13,u21,u22,u23,u31,u32,u33,u1u1,u1u2,u1u3,u2u1,u2u2,u2u3,u3u1,u3u2,u3u3)
        return Tres_slice

    for z_index in range(len(geom[0])):
        Tres[:,:,:,:,z_index] = GetTres(z_index)

    gc.collect()
    comm.Barrier()
#    print(f'{case_id} file {fileidx} at k= {k_}: \t Tres done ::: \t took {end - start:.2f} seconds', flush=True)

    # ---------- ---------- ---------- Calculate KE Flux at k_ ---------- ---------- ----------#
    KEFlux = np.einsum('ijxyz,ijxyz->xyz', -1*Tres, S)  # Pi = < -Tres : S >_A,t

    PI_area_avg = np.zeros(len(geom[0]))
    for z_index in range(len(geom[1])):
        PI_area_avg[z_index] = areaAv(KEFlux[:,:,z_index],geom)

    gc.collect()
    PI = lineAv_cartesian(PI_area_avg,np.asarray(geom[1]))
#    print(f'{case_id} file {fileidx} at k= {k_}: \t Contraction, Avg done ::: \t took {end - start:.2f} seconds', flush=True)

    # ---------- ---------- ---------- Saving ---------- ---------- ----------#
    if comm.Get_rank()==0:
        result_dict = {
        f'TKEFluxZ_k{k_}_n{fileidx}' : PI_area_avg,
        f'TKEFlux_k{k_}_n{fileidx}' : PI,
        }
        save_variable_arrays(output_dir, case_id, result_dict)
    #    print(f'{case_id} file {fileidx} at k= {k_}: \t Saved', flush=True)
        del(result_dict)

    del (geom, Rp, Phi, Zp, 
         Uz_interp_rpz, Uphi_interp_rpz, Ur_interp_rpz,
         Uz, Ur, Uphi,
         Uz_filtered, Ur_filtered, Uphi_filtered,
         Tres, S, KEFlux, PI_area_avg, PI)
    gc.collect()
    comm.Barrier()

    # ---------- ---------- ---------- Cleanup ---------- ---------- ----------#
    end=time.time()

    if comm.Get_rank() == 0:
        print(f'{case_id} file {fileidx} at k= {k_} ended, took ', end-start,' seconds', flush=True)

# ---------- ---------- ---------- MAIN PROGRAM ---------- ---------- ----------#
# ---------- ---------- ---------- Specify cases ---------- ---------- ----------#
base_cases = { 
#    "C13_Ro01_Ra1e9": Path(r'/home/leek198/C1_3_v13_Ro01_Ra1e9/'),
#    "C12_Ro02_Ra1e8": Path(r'/home/leek198/C1_2_v14_Ro02_Ra1e8/'),
    "C11_Ro02_Ra1e8": Path(r'/home/leek198/C1_1_v9_Ro02_Ra1e8/'),
#    "C11_RoInf_Ra1e9": Path(r'/home/leek198/C1_1_v9_RoInf_Ra1e8/'),
}

# ---------- ---------- ---------- Specify k values ---------- ---------- ----------#
karr=np.unique(np.int64(np.logspace(0,4,128)))

# ---------- ---------- ---------- Specify output dir ---------- ---------- ----------#
output_dir = Path("../results")

# ---------- ---------- ---------- Specify file spacing ---------- ---------- ----------#
# 20 files spaced 3 files apart (i.e. 1.5 time-units apart)
fileidx = np.int64(-1 - 3*np.arange(20))

# ---------- ---------- ---------- Loop over K ---------- ---------- ----------#
# Convolution takes the most time, so a FORTRAN module is now compiled and
# imported as a python module using numpy.f2py to allow for multi-node computation
for case_id, case_path in base_cases.items():
    for idx in fileidx:
        for k in karr:
            process_case(case_id, case_path, output_dir, k, idx)
