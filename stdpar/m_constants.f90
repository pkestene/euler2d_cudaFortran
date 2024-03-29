!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_constants

  use m_precision

  implicit none

  ! index to hydro variables (start at 0 to ease porting from C)
  integer(int_kind), parameter :: ID=1
  integer(int_kind), parameter :: IP=2
  integer(int_kind), parameter :: IU=3
  integer(int_kind), parameter :: IV=4

  !! other parameters
  integer(int_kind),parameter  :: NBVAR=4

  ! Riemann solver type for hydro fluxes
  ! used to initialize riemannSolverType in HydroParameters module
  integer(int_kind), parameter :: APPROX=0 !< quasi-exact Riemann solver
  integer(int_kind), parameter :: HLL=1    !< HLL  hydro Riemann solver
  integer(int_kind), parameter :: HLLC=2   !< HLLC hydro Riemann solver

  ! BoundaryConditionType
  integer(int_kind), parameter :: BC_UNDEFINED=0
  integer(int_kind), parameter :: BC_DIRICHLET=1 !< reflecting border
  integer(int_kind), parameter :: BC_NEUMANN=2   !< absorbing border
  integer(int_kind), parameter :: BC_PERIODIC=3  !< periodic border
  integer(int_kind), parameter :: BC_COPY=4      !< only used in MPI parallelized version

  ! Boundary direction
  integer(int_kind), parameter :: XDIR=0
  integer(int_kind), parameter :: YDIR=1

  ! component index
  integer(int_kind), parameter :: IX=1
  integer(int_kind), parameter :: IY=2

  ! variables names (used in file IO)
  character(3), dimension(4), parameter :: varNames = (/ 'rho', '  E', ' mx', ' my' /)

  ! named constants
  real(fp_kind), parameter :: ZERO_F       = 0.0
  real(fp_kind), parameter :: HALF_F       = 0.5
  real(fp_kind), parameter :: ONE_FOURTH_F = 0.25
  real(fp_kind), parameter :: ONE_F        = 1.0
  real(fp_kind), parameter :: TWO_F        = 2.0

  ! face Id's
  integer(int_kind), parameter :: FACE_XMIN = 0
  integer(int_kind), parameter :: FACE_XMAX = 1
  integer(int_kind), parameter :: FACE_YMIN = 2
  integer(int_kind), parameter :: FACE_YMAX = 3
  integer(int_kind), parameter :: FACE_MIN = 0
  integer(int_kind), parameter :: FACE_MAX = 1

end module m_constants
