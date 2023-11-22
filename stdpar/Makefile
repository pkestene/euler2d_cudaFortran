# ######## NVHPC ########
F90 = nvfortran
FFLAGS = -O4 -fast -Minform=warn -Minfo
LDFLAGS = -O4

#CUDA_ARCH ?= cc75
#CUDA_VERSION ?= 11.6
FFLAGS_CUDA = -O3 -stdpar=gpu -gpu=ccnative -acc=gpu -gpu=rdc,managed -Minform=warn -Minfo -Mpreprocess
LDFLAGS_CUDA = -O3 -stdpar=gpu -gpu=ccnative -acc=gpu -gpu=rdc,managed

NVTX_LDFLAGS = -L$(NVHPC_ROOT)/cuda/lib64 -lnvToolsExt

SRCDIR = .

SRC = \
	m_nvtx.f90 \
	m_precision.f90 \
	m_constants.f90 \
	m_parameters.f90 \
	m_monitoring.f90 \
	m_utils.f90 \
	m_run.f90 \
	main.f90

CUDA_OBJ = $(SRC:.f90=.o)

all: euler2d_stdpar_gpu

euler2d_stdpar_gpu: $(CUDA_OBJ)
	$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@ $(NVTX_LDFLAGS)

clean:
	rm -f *.o *.mod euler2d_stdpar_cpu euler2d_stdpar_gpu

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS_CUDA) -c $<
