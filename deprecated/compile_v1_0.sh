#!/bin/bash

rm -rf build

F90="mpif90" FC="mpif90" CC="mpicc" FFLAGS="-g -O3 -Wall -I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi/" LD_LIBRARY_PATH="/usr/lib/:/usr/local/lib/" python3 -m numpy.f2py -c convolution_2D_Polar_mpi_v1_0.f90 -m conv2D -I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15/openmpi/ --backend meson --build-dir ./build/ -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lm
