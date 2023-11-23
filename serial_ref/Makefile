# ######## GNU ########
F90 = gfortran
FFLAGS = -O3 -Wall -Wextra -cpp -ffixed-line-length-none -ffree-line-length-none
LDFLAGS = -O3

#NVTX_LDFLAGS = -L$(NVHPC_ROOT)/cuda/lib64 -lnvToolsExt

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

all: euler2d

euler2d: $(CUDA_OBJ)
	#$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@ $(NVTX_LDFLAGS)
	$(F90) $(LDFLAGS_CUDA) $(CUDA_OBJ) -o $@

clean:
	rm -f *.o *.mod euler2d euler2d

cleanall: clean
	rm -f *.vti *.*.gpu *.*.h crt1.reg.c

%.o:    $(SRCDIR)/%.f90
	$(F90) $(FFLAGS) -c $<
