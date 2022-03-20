# ######## GNU ########
#F90 = gfortran
#FFLAGS = -O3 -Wall
#FFLAGS = -g -Wall
# -ffpe-trap=zero,invalid,underflow -fbacktrace
#LDFLAGS = -O3

# ######## INTEL ########
# F90 = ifort
# FFLAGS = -O3
# #FFLAGS = -g -Wall
# LDFLAGS = -O3

# ######## NVIDIA stdpar multicore ########
#F90 = nvfortran
#FFLAGS = -O3 -stdpar=multicore -acc=multicore  -Minform=warn -Minfo
#LDFLAGS = -O3 -acc=multicore

# ######## NVIDIA stdpar gpu ########
F90 = nvfortran
ARCH_GPU=cc75
CUDA_VERSION=11.6
FFLAGS = -O3 -stdpar=gpu -gpu=$(ARCH_GPU),cuda$(CUDA_VERSION) -acc=gpu -gpu=rdc -Minform=warn -Minfo
LDFLAGS = -O3 -stdpar=gpu -acc=gpu -gpu=rdc

SRCDIR = .

CUDA_SRC = \
	m_HydroPrecision.f90 \
	m_HydroConstants.f90 \
	m_HydroParameters.f90 \
	m_Monitoring.f90 \
	m_HydroUtils.f90 \
	m_HydroRun.f90 \
	main.f90

CUDA_OBJ = $(CUDA_SRC:.f90=.o)

all: euler2d

euler2d: $(CUDA_OBJ)
	$(F90) $(LDFLAGS) $(CUDA_OBJ) -o $@

clean:
	rm -f *.o *.mod euler2d

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS) -c $<
