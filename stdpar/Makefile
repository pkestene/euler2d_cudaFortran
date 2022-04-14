# ######## PGI ########
F90 = nvfortran
FFLAGS = -O4 -fast -Minform=warn -Minfo
LDFLAGS = -O4

CUDA_ARCH ?= cc75
CUDA_VERSION ?= 11.6
FFLAGS_CUDA = -O3 -stdpar=gpu -gpu=cuda$(CUDA_VERSION),$(CUDA_ARCH) -acc=gpu -gpu=rdc,managed -Minform=warn -Minfo -Mpreprocess
LDFLAGS_CUDA = -O3 -stdpar=gpu -gpu=cuda$(CUDA_VERSION),$(CUDA_ARCH) -acc=gpu -gpu=rdc,managed

SRCDIR = .

SRC = \
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
	$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@

clean:
	rm -f *.o *.mod euler2d_stdpar_cpu euler2d_stdpar_gpu

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS_CUDA) -c $<
