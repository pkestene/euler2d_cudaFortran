!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_parameters

  use m_precision

  implicit none

  type HydroPar
     !! hydro (numerical scheme) parameters
     real   (fp_kind) :: gamma0        !< specific heat capacity ratio (adiabatic index)
     real   (fp_kind) :: gamma6
     real   (fp_kind) :: cfl           !< Courant-Friedrich-Lewy parameter.
     integer(int_kind):: slope_type    !< type of slope computation (2 for second order scheme).
     integer(int_kind):: iorder        !<
     real   (fp_kind) :: smallr        !< small density cut-off
     real   (fp_kind) :: smallc        !< small speed of sound cut-off
     real   (fp_kind) :: smallp        !< small pressure cut-off
     real   (fp_kind) :: smallpp       !< smallp times smallr
     integer(int_kind) :: niter_riemann=10 !< number of iteration usd in quasi-exact riemann solver
     integer(int_kind) :: isize=0  !< total size (in cell unit) along X direction with ghosts.
     integer(int_kind) :: jsize=0  !< total size (in cell unit) along Y direction with ghosts.

     integer(int_kind) :: nx=0     !< logical size along X (without ghost cells).
     integer(int_kind) :: ny=0     !< logical size along Y (without ghost cells).

     real(fp_kind)     :: dx       !< x resolution
     real(fp_kind)     :: dy       !< y resolution

     integer(int_kind) :: boundary_type_xmin=1
     integer(int_kind) :: boundary_type_xmax=1
     integer(int_kind) :: boundary_type_ymin=1
     integer(int_kind) :: boundary_type_ymax=1

  end type HydroPar

end module m_parameters
