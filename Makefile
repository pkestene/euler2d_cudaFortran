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

# ######## PGI ########
F90 = pgfortran
FFLAGS = -O4 -fast -Minform=warn -Minfo
LDFLAGS = -O4

ARCH_GPU=cc75
#FFLAGS_CUDA = -v -O3 -Mcuda -Mcuda=6.5 -Mcuda=$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,fastmath,mul24,maxregcount:48,time -Minfo=all -Mpreprocess  -Mcuda=rdc
FFLAGS_CUDA = -v -O3 -Mcuda -Mcuda=cuda11.4,$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,maxregcount:48,time -Minfo=all -Mpreprocess  -Mcuda=rdc -acc
#-Mcuda=keepgpu
LDFLAGS_CUDA = -O3 -Mcuda=cuda11.4,$(ARCH_GPU) -Mcuda=rdc -acc

# PGI Debug #
# F90 = pgfortran
# FFLAGS = -fast -Minform=warn -Minfo -g
# LDFLAGS =

# ARCH_GPU=cc35
# FFLAGS_CUDA = -g -v -Mcuda -Mcuda=5.5 -Mcuda=$(ARCH_GPU) -Mcuda=ptxinfo -ta=nvidia,fastmath,mul24,maxregcount:48,time -Minfo=all -Mpreprocess -Mcuda=keepgpu
# LDFLAGS_CUDA = -Mcuda=5.5 -Mcuda=$(ARCH_GPU)

SRCDIR = .

CUDA_SRC = \
	m_HydroPrecision.cuf \
	m_HydroConstants.cuf \
	m_HydroParameters.cuf \
	m_Monitoring_gpu.cuf \
	m_HydroParameters_gpu.cuf \
	m_HydroUtils_gpu.cuf \
	m_HydroRun_gpu.cuf \
	main_gpu.cuf

CUDA_OBJ = $(CUDA_SRC:.cuf=.o)

all: euler2d_gpu

euler2d_gpu: $(CUDA_OBJ)
	$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@


clean:
	rm -f *.o *.mod euler2d_cpu euler2d_gpu

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS) -c $<

%.o:    $(SRCDIR)/%.cuf
	$(F90) $(FFLAGS_CUDA) -c $<
