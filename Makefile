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

# ######## NVIDIA ########
F90 = nvfortran
#FFLAGS = -O3 -stdpar=multicore -Minform=warn -Minfo
FFLAGS = -O3 -stdpar=gpu -Minform=warn -Minfo
LDFLAGS = -O3

ARCH_GPU=cc75
#FFLAGS_CUDA = -v -O3 -Mcuda -Mcuda=6.5 -Mcuda=$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,fastmath,mul24,maxregcount:48,time -Minfo=all -Mpreprocess  -Mcuda=rdc
#FFLAGS_CUDA = -v -O3 -Mcuda -Mcuda=cuda7.5,$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,maxregcount:48,time -Minfo=all -Mpreprocess  -Mcuda=rdc -acc
#-Mcuda=keepgpu
#LDFLAGS_CUDA = -O3 -Mcuda=cuda7.5,$(ARCH_GPU) -Mcuda=rdc -acc

# ARCH_GPU=cc35
# FFLAGS_CUDA = -g -v -Mcuda -Mcuda=5.5 -Mcuda=$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,fastmath,mul24,maxregcount:48,time -Minfo=all -Mpreprocess -Mcuda=keepgpu
# LDFLAGS_CUDA = -Mcuda=5.5 -Mcuda=$(ARCH_GPU)

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
	$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@

clean:
	rm -f *.o *.mod euler2d

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS) -c $<

