# Compiler
FC = gfortran
MPIFC = mpif90
#MPIFC = /home/klee/Tau/tau-2.35.1/x86_64/bin/tau_f90.sh

# Flags
FFLAGS = -g -O3 -fopenmp #-fcheck=all -fbacktrace
#FFLAGS = -g -o3 -fopenmp -tau_makefile=/home/klee/Tau/tau-2.35.1/x86_64/lib/Makefile.tau-mpi -tau_options=-optKeepFiles -tau_options=-optCompInst
# enable fcheck for debugging and backtrace
FFTW_INC = -I/usr/include
INCLUDE = -I/usr/local/include

# Libraries
LIBS = -L/usr/local/lib -lnetcdff -lnetcdf -lpthread -llapack -lblas -lfftw3

# Object files
OBJ = mpi_context.o ncinput.o ncoutput.o domain_decomposition.o io.o interpolate.o conv3d.o centralfd.o SBS_KEFlux_v9.o

# Executable name
TARGET = SBS_KEFlux_Polar_3D

# =========================
# Build executable
# =========================
all: $(TARGET)

$(TARGET): $(OBJ)
	$(MPIFC) $(FFLAGS) -o $(TARGET) $(OBJ) $(LIBS)

# =========================
# Compile object file
# =========================
mpi_context.o: mpi_context.f90
	$(MPIFC) $(FFLAGS) -c mpi_context.f90 $(INCLUDE) $(LIBS)

ncinput.o: ncinput.f90
	$(MPIFC) $(FFLAGS) -c ncinput.f90 $(INCLUDE)

ncoutput.o: ncoutput.f90
	$(MPIFC) $(FFLAGS) -c ncoutput.f90 $(INCLUDE)

domain_decomposition.o: domain_decomposition.f90
	$(MPIFC) $(FFLAGS) -c domain_decomposition.f90 $(INCLUDE)

io.o: io.f90
	$(MPIFC) $(FFLAGS) -c io.f90 $(INCLUDE)

interpolate.o: interpolate.f90
	$(MPIFC) $(FFLAGS) -c interpolate.f90 $(INCLUDE)

conv3d.o: conv3d.f90
	$(MPIFC) $(FFLAGS) -c conv3d.f90 $(INCLUDE)

centralfd.o: centralfd.f90
	$(MPIFC) $(FFLAGS) -c centralfd.f90 $(INCLUDE)

SBS_KEFlux_v9.o: SBS_KEFlux_v9.f90
	$(MPIFC) $(FFLAGS) -c SBS_KEFlux_v9.f90 $(FFTW_INC) $(INCLUDE)
# =========================
# Clean
# =========================
clean:
	rm -f profile.* *.o *.mod $(TARGET)
