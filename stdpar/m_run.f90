!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_run

  use m_precision
  use m_constants
  use m_parameters
  use m_utils
  use m_nvtx

  implicit none

  real(fp_kind), dimension(:,:,:), allocatable :: u
  real(fp_kind), dimension(:,:,:), allocatable :: u2
  real(fp_kind), dimension(:,:,:), allocatable :: q
  real(fp_kind), dimension(:,:,:), allocatable :: fx,fy

  type(HydroPar)    :: hydroParams
  !$acc declare copyin(hydroParams)

  !! run parameters
  integer(int_kind ) :: nStepmax !< maximun number of time steps.
  real   (fp_kind)   :: tEnd     !< end of simulation time.
  integer(int_kind ) :: nOutput  !< number of time steps between 2 consecutive outputs.

  !! geometry parameters
  integer(int_kind), parameter :: ghostWidth=2
  integer(int_kind) :: imin=0   !< index minimum at X border
  integer(int_kind) :: imax=0   !< index maximum at X border
  integer(int_kind) :: jmin=0   !< index minimum at Y border
  integer(int_kind) :: jmax=0   !< index maximum at Y border

  real(fp_kind)   :: xmin=0.0
  real(fp_kind)   :: xmax=1.0
  real(fp_kind)   :: ymin=0.0
  real(fp_kind)   :: ymax=1.0

  !! IO parameters
  logical :: ioVTK=.true.   !< enable VTK  output file format (using VTI).
  logical :: ioHDF5=.false. !< enable HDF5 output file format.

  integer(int_kind) :: riemannSolverType=0
  character(LEN=20) :: problem='implode'

contains

  subroutine initHydroParameters()

    implicit none

    ! local variables
    integer(int_kind) :: narg
    character(LEN=80) :: inputFilename

    integer(int_kind) :: nx
    integer(int_kind) :: ny
    real   (fp_kind) :: gamma0
    real   (fp_kind) :: cfl
    real   (fp_kind) :: smallr
    real   (fp_kind) :: smallc
    integer(int_kind) :: niter_riemann=10
    integer(int_kind):: iorder
    real   (fp_kind) :: slope_type
    integer(int_kind) :: boundary_type_xmin=1
    integer(int_kind) :: boundary_type_xmax=1
    integer(int_kind) :: boundary_type_ymin=1
    integer(int_kind) :: boundary_type_ymax=1

    ! declare namelist
    namelist/run/tEnd,nStepmax,nOutput
    namelist/mesh/nx,ny,boundary_type_xmin,boundary_type_xmax,boundary_type_ymin,boundary_type_ymax
    namelist/hydro/gamma0,cfl,smallr,smallc,niter_riemann, &
         &         iorder,slope_type

    narg = COMMAND_ARGUMENT_COUNT()
    if (narg .ne. 1) then
       write(*,*)'You should provide an input parameter file (with a namelist)'
       write(*,*)'Usage:'
       write(*,*)'   ./euler2d param.nml'
       stop
    end if

    ! retrieve input parameter file name
    call GET_COMMAND_ARGUMENT(1,inputFilename)

    ! read namelist
    open(1,file=inputFilename)
    read(1,NML=run)
    read(1,NML=mesh)
    read(1,NML=hydro)
    close(1)

    ! set other parameters
    hydroParams%nx = nx
    hydroParams%ny = ny

    hydroParams%boundary_type_xmin = boundary_type_xmin
    hydroParams%boundary_type_xmax = boundary_type_xmax
    hydroParams%boundary_type_ymin = boundary_type_ymin
    hydroParams%boundary_type_ymax = boundary_type_ymax

    imin = 1
    jmin = 1
    imax = hydroParams%nx+2*ghostWidth
    jmax = hydroParams%ny+2*ghostWidth

    hydroParams%isize = imax - imin + 1
    hydroParams%jsize = jmax - jmin + 1

    hydroParams%dx = (xmax - xmin) / nx
    hydroParams%dy = (ymax - ymin) / ny

    hydroParams%gamma0        = gamma0
    hydroParams%cfl           = cfl
    hydroParams%smallr        = smallr
    hydroParams%smallc        = smallc
    hydroParams%niter_riemann = niter_riemann
    hydroParams%iorder        = iorder
    hydroParams%slope_type    = slope_type

    hydroParams%smallp  = hydroParams%smallc * hydroParams%smallc / hydroParams%gamma0
    hydroParams%smallpp = hydroParams%smallr * hydroParams%smallp
    hydroParams%gamma6  = (hydroParams%gamma0 + ONE_F) / (TWO_F * hydroParams%gamma0)

    !$acc update device(hydroParams)

  end subroutine initHydroParameters

  ! dump hydro parameters on stdout
  subroutine printHydroParameters()

    implicit none

    write(*,*) 'Simulation run parameters:'
    write(*,*) 'nx : ', hydroParams%nx
    write(*,*) 'ny : ', hydroParams%ny
    write(*,*) 'dx : ', hydroParams%dx
    write(*,*) 'dy : ', hydroParams%dy
    write(*,*) 'imin : ', imin, 'imax : ', imax
    write(*,*) 'jmin : ', jmin, 'jmax : ', jmax
    write(*,*) 'nStepmax,tEnd,nOutput : ',nStepmax,tEnd,nOutput
    write(*,*) 'gamma0,cfl : ',hydroParams%gamma0,hydroParams%cfl
    write(*,*) 'smallr,smallc,niter_riemann :',hydroParams%smallr,hydroParams%smallc,hydroParams%niter_riemann
    write(*,*) 'iorder,slope_type :',hydroParams%iorder,hydroParams%slope_type

  end subroutine printHydroParameters

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine run_init()

    implicit none

    allocate( u (hydroParams%isize, hydroParams%jsize, NBVAR) )
    allocate( u2(hydroParams%isize, hydroParams%jsize, NBVAR) )
    allocate( q (hydroParams%isize, hydroParams%jsize, NBVAR) )
    allocate( fx(hydroParams%isize, hydroParams%jsize, NBVAR) )
    allocate( fy(hydroParams%isize, hydroParams%jsize, NBVAR) )

    ! initialize u (at t=0)
    select case (problem)
    case('implode') ! discontinuity line along the domain diagonal
       call init_implode(hydroParams,u)
    case default
       write(*,*) 'Unknown problem; default to implode'
       call init_implode(hydroParams,u)
    end select

  end subroutine run_init

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine run_cleanup

    implicit none

    deallocate(u)
    deallocate(u2)
    deallocate(q)
    deallocate(fx)
    deallocate(fy)

  end subroutine run_cleanup

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Wrapper to the actual computation routine
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine godunov_unsplit(nStep, dt)

    implicit none

    ! dummy variables
    integer      , intent(in) :: nStep
    real(fp_kind), intent(in) :: dt

    call nvtxStartRange("godunov_compute",0)

    if ( modulo(nStep,2) == 0 ) then
       call godunov_unsplit_cpu(hydroParams, u , u2, dt)
    else
       call godunov_unsplit_cpu(hydroParams, u2, u , dt)
    end if

    call nvtxEndRange()

  end subroutine godunov_unsplit

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Actual CPU computation of Godunov scheme
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine godunov_unsplit_cpu(params, data_in, data_out, dt)

    use m_monitoring

    implicit none

    ! dummy variables
    type(HydroPar)                                             , intent(in)    :: params
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(inout) :: data_in
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(inout) :: data_out
    real(fp_kind)                                              , intent(in)    :: dt

    ! local variables
    integer :: i, j, ii, jj, iVar

    ! fill ghost cell in data_in
    call timerStart(boundaries_timer)
    call make_boundaries(params, data_in)
    call timerStop(boundaries_timer)

    ! copy data_in into data_out (not necessary) and if it where necessary
    ! it should be done in a parallel region !
    !data_out = data_in

    ! start main computation
    call timerStart(godunov_timer)

    ! convert conservative variable into primitives ones for the entire domain
    call convertToPrimitives(params, data_in, q)

    call computeAndStoreFluxes(params, data_in, dt)

    call updateHydro(params, data_in, data_out)

    call timerStop(godunov_timer)

  end subroutine godunov_unsplit_cpu

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute time step reduction
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_dt(params, dt, useU)
    implicit none

    ! dummy variables
    type(HydroPar)   , intent(in)   :: params
    real   (fp_kind) , intent(out)  :: dt
    integer(int_kind), intent(in)   :: useU

    ! local variables
    real(fp_kind) :: invDt=0.0
    real(fp_kind) :: vx,vy
    integer :: i,j
    real(fp_kind), dimension(nbVar) :: uLoc
    real(fp_kind), dimension(nbVar) :: qLoc
    real(fp_kind)                   :: c

    call nvtxStartRange("compute_dt",1)

    ! for loop over inner region
    if (useU == 0) then

       do concurrent (j=ghostWidth+1:params%jsize-ghostWidth-1, i=ghostWidth+1:params%isize-ghostWidth-1) reduce(max:invDt) &
            & local(uLoc,qLoc,vx,vy,c) shared(params)

          ! retrieve conservative variables in current cell
          uLoc(ID) = u(i,j,ID)
          uLoc(IP) = u(i,j,IP)
          uLoc(IU) = u(i,j,IU)
          uLoc(IV) = u(i,j,IV)

          call computePrimitives(uLoc,qLoc,c,params)
          vx = c + abs(qLoc(IU))
          vy = c + abs(qLoc(IV))
          invDt = max(invDt, vx/params%dx + vy/params%dy)

       end do

    else

       do concurrent (j=ghostWidth+1:params%jsize-ghostWidth-1, i=ghostWidth+1:params%isize-ghostWidth-1) reduce(max:invDt) &
            & local(uLoc,qLoc,vx,vy,c) shared(params)

          ! retrieve conservative variables in current cell
          uLoc(ID) = u2(i,j,ID)
          uLoc(IP) = u2(i,j,IP)
          uLoc(IU) = u2(i,j,IU)
          uLoc(IV) = u2(i,j,IV)

          call computePrimitives(uLoc,qLoc,c,params)
          vx = c + abs(qLoc(IU))
          vy = c + abs(qLoc(IV))
          invDt = max(invDt, vx/params%dx + vy/params%dy)

       end do
    end if


    dt = params%cfl / invDt

    call nvtxEndRange()

  end subroutine compute_dt

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Convert conservative variables array into primitive var array (q)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine convertToPrimitives(params,uData,qData)

    implicit none

    ! dummy variables
    type(HydroPar)                                             , intent(in)    :: params
    real(fp_kind), dimension(params%isize, params%jsize, NBVAR), intent(inout) :: uData
    real(fp_kind), dimension(params%isize, params%jsize, NBVAR), intent(inout) :: qData

    ! local variables

    ! conservative variable state vector
    real(fp_kind), dimension(NBVAR) :: uLoc

    ! primitive variable state vector
    real(fp_kind), dimension(NBVAR) :: qLoc
    real(fp_kind)                   :: c
    integer :: i,j
    real(fp_kind)                   :: eken ! kinetic energy
    real(fp_kind)                   :: e    ! total energy

    call nvtxStartRange("convertToPrimitives",2)

    do concurrent (j=1:params%jsize, i=1:params%isize) local(uLoc,qLoc,c,e,eken) shared(params)

       uLoc(ID) = uData(i,j,ID)
       uLoc(IP) = uData(i,j,IP)
       uLoc(IU) = uData(i,j,IU)
       uLoc(IV) = uData(i,j,IV)

       call computePrimitives(uLoc, qLoc, c, params)

       qData(i,j,ID) = qLoc(ID)
       qData(i,j,IP) = qLoc(IP)
       qData(i,j,IU) = qLoc(IU)
       qData(i,j,IV) = qLoc(IV)

    end do ! end do i,j

    call nvtxEndRange()

  end subroutine convertToPrimitives

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute flux via Riemann solver and update (time integration)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine computeAndStoreFluxes(params, data, dt)

    implicit none

    ! dummy variables
    type(HydroPar)                                             , intent(in)    :: params
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(inout) :: data
    real(fp_kind)                                              , intent(in)    :: dt

    ! local variables
    integer :: i,j,ivar
    real(fp_kind), dimension(nbVar)     :: qleft, qright
    real(fp_kind), dimension(nbVar)     :: flux_x, flux_y
    real(fp_kind), dimension(nbVar)     :: qgdnv
    real(fp_kind), dimension(nbVar)     :: qLoc, qLocN, qN0, qN1, qN2, qN3
    real(fp_kind), dimension(2,nbVar)   :: dq, dqN

    real(fp_kind) :: tmp

    real(fp_kind) :: dtdx
    real(fp_kind) :: dtdy

    dtdx = dt / params%dx
    dtdy = dt / params%dy

    call nvtxStartRange("computeAndStoreFluxes",3)

    do concurrent(j=ghostWidth+1:params%jsize-ghostWidth+1, i=ghostWidth+1:params%isize-ghostWidth+1) &
         & local(ivar,qleft,qright,flux_x,flux_y,qgdnv,qLoc,qLocN,qN0,qN1,qN2,qN3,dq,dqN) &
         & shared(params,dtdx,dtdy)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! deal with left interface along X !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !! get primitive variables state vector
       do ivar=1,nbVar
          qLoc(ivar) = q(i  ,j  , ivar)
          qN0(ivar)  = q(i+1,j  , ivar)
          qN1(ivar)  = q(i-1,j  , ivar)
          qN2(ivar)  = q(i  ,j+1, ivar)
          qN3(ivar)  = q(i  ,j-1, ivar)
       end do

       call slope_unsplit_hydro_2d(qLoc, qN0, qN1, qN2, qN3, dq, params)

       ! slopes at left neighbor along X
       do ivar=1,nbVar
          qLocN(ivar) = q(i-1,j  , ivar)
          qN0(ivar)   = q(i  ,j  , ivar)
          qN1(ivar)   = q(i-2,j  , ivar)
          qN2(ivar)   = q(i-1,j+1, ivar)
          qN3(ivar)   = q(i-1,j-1, ivar)
       end do

       call slope_unsplit_hydro_2d(qLocN, qN0, qN1, qN2, qN3, dqN, params)

       !
       ! compute reconstructed states at left interface along X
       !

       ! left interface : right state
       call trace_unsplit_2d_along_dir(qLoc, dq(1,:), dq(2,:), dtdx, dtdy, FACE_XMIN, qright, params)

       ! left interface : left state
       call trace_unsplit_2d_along_dir(qLocN, dqN(1,:), dqN(2,:), dtdx, dtdy, FACE_XMAX, qleft, params)

       ! compute hydro flux_x
       call riemann_2d(qleft,qright,qgdnv,flux_x, params)

       !
       ! store fluxes X
       !
       do ivar=1,nbVar
          fx(i  ,j  , ivar) = flux_x(ivar) * dtdx
       end do

       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! deal with left interface along Y !
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !! slopes at left neighbor along Y
       do ivar=1,nbVar
          qLocN(ivar) = q(i  ,j-1, ivar)
          qN0(ivar)   = q(i+1,j-1, ivar)
          qN1(ivar)   = q(i-1,j-1, ivar)
          qN2(ivar)   = q(i  ,j  , ivar)
          qN3(ivar)   = q(i  ,j-2, ivar)
       end do

       call slope_unsplit_hydro_2d(qLocN, qN0, qN1, qN2, qN3, dqN, params)

       !!
       !! compute reconstructed states at left interface along Y
       !!

       !! left interface : right state
       call trace_unsplit_2d_along_dir(qLoc, dq(1,:), dq(2,:), dtdx, dtdy, FACE_YMIN, qright, params)

       !! left interface : left state
       call trace_unsplit_2d_along_dir(qLocN, dqN(1,:), dqN(2,:), dtdx, dtdy, FACE_YMAX, qleft, params)

       ! swap IU / IV
       tmp = qleft(IU); qleft(IU) = qleft(IV); qleft(IV) = tmp
       tmp = qright(IU); qright(IU) = qright(IV); qright(IV) = tmp

       ! compute hydro flux_y
       call riemann_2d(qleft,qright,qgdnv,flux_y, params)

       tmp = flux_y(IU); flux_y(IU) = flux_y(IV); flux_y(IV) = tmp

       !
       ! store fluxes Y
       !
       do ivar=1,nbVar
          fy(i  ,j  , ivar) = flux_y(ivar) * dtdy
       end do

    end do

    call nvtxEndRange()

  end subroutine computeAndStoreFluxes

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Update hydro
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine updateHydro(params, dataIn, dataOut)

    implicit none

    ! dummy variables
    type(HydroPar)                                             , intent(in)  :: params
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(in)  :: dataIn
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(out) :: dataOut

    ! local variables
    integer :: i,j,ivar
    real(fp_kind) :: flux_tot

    call nvtxStartRange("updateHydro",4)

    do concurrent(j=ghostWidth+1:params%jsize-ghostWidth, i=ghostWidth+1:params%isize-ghostWidth) &
         & local(ivar,flux_tot)

       do ivar=1,nbVar
          flux_tot  = fx(i  ,j  ,ivar) - fx(i+1,j  ,ivar)
          flux_tot  = fy(i  ,j  ,ivar) - fy(i  ,j+1,ivar) + flux_tot
          dataOut(i, j, ivar) = dataIn(i, j, ivar) + flux_tot
       end do

    end do

    call nvtxEndRange()

  end subroutine updateHydro

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Hydrodynamical Implosion Test :
  !! http://www.astro.princeton.edu/~jstone/tests/implode/Implode.html
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_implode(params, data)

    implicit none

    ! dummy variables
    type(HydroPar)                                             , intent(in)    :: params
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(inout) :: data

    ! local variables
    integer :: i,j
    real(fp_kind) :: tmp

    do concurrent (j=jmin:jmax, i=imin:imax) local(tmp) shared(params)

       tmp = 1.0*(i-ghostWidth-1)/params%nx + 1.0*(j-ghostWidth-1)/params%ny

       if (tmp .gt. 0.5) then
          data(i,j,ID)=1.0
          data(i,j,IP)=1.0/(params%gamma0-1.0)
          data(i,j,IU)=0.0
          data(i,j,IV)=0.0
       else
          data(i,j,ID)=0.125
          data(i,j,IP)=0.14/(params%gamma0-1.0)
          data(i,j,IU)=0.0
          data(i,j,IV)=0.0
       end if
    end do

  end subroutine init_implode

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Fill ghost cells according to border condition :
  !! absorbant, reflexive or periodic
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_boundaries(params, data)
    implicit none
    ! dummy variables
    type(HydroPar)                                             , intent(in)    :: params
    real(fp_kind), dimension(params%isize, params%jsize, nbVar), intent(inout) :: data

    ! local variables
    integer ::i,j,i0,j0,iVar
    real(fp_kind) :: dir

    call nvtxStartRange("make_boundaries",5)

    ! boundary xmin
    do concurrent(j=jmin+ghostWidth:jmax-ghostWidth) local(iVar,i,dir,i0) shared(params)

       do iVar=1,nbVar
          do i=1,ghostWidth
             dir=1.0
             if(params%boundary_type_xmin==1)then
                i0=2*ghostWidth+1-i
                if(iVar==IU)dir=-1.0
             else if(params%boundary_type_xmin==2)then
                i0=ghostWidth+1
             else ! periodic
                i0=params%nx+i
             end if
             data(i,j,iVar)=data(i0,j,iVar)*dir
          end do
       end do
    end do

    ! boundary xmax
    do concurrent(j=jmin+ghostWidth:jmax-ghostWidth) local(iVar,i,dir,i0)
       do iVar=1,nbVar
          do i=params%nx+ghostWidth+1,params%nx+2*ghostWidth
             dir=1.0
             if(params%boundary_type_xmax==1)then
                i0=2*params%nx+2*ghostWidth+1-i
                if(iVar==IU)dir=-1.0
             else if(params%boundary_type_xmax==2)then
                i0=params%nx+ghostWidth
             else ! periodic
                i0=i-params%nx
             end if
             data(i,j,iVar)=data(i0,j,iVar)*dir
          end do
       end do
    end do

    ! boundary ymin
    do concurrent(i=imin+ghostWidth:imax-ghostWidth) local(iVar,j,dir,j0)
       do iVar=1,nbVar
          do j=1,ghostWidth
             dir=1.0
             if(params%boundary_type_ymin==1)then
                j0=2*ghostWidth+1-j
                if(iVar==IV)dir=-1.0
             else if(params%boundary_type_ymin==2)then
                j0=ghostWidth+1
             else ! periodic
                j0=params%ny+j
             end if
             data(i,j,iVar)=data(i,j0,iVar)*dir
          end do
       end do
    end do

    ! boundary ymax
    do concurrent(i=imin+ghostWidth:imax-ghostWidth) local(iVar,j,dir,j0)
       do iVar=1,nbVar
          do j=params%ny+ghostWidth+1,params%ny+2*ghostWidth
             dir=1.0
             if(params%boundary_type_ymax==1)then
                j0=2*params%ny+2*ghostWidth+1-j
                if(iVar==IV)dir=-1.0
             else if(params%boundary_type_ymax==2)then
                j0=params%ny+ghostWidth
             else ! periodic
                j0=j-params%ny
             end if
             data(i,j,iVar)=data(i,j0,iVar)*dir
          end do
       end do
    end do

    call nvtxEndRange()

  end subroutine make_boundaries

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! output routine (VTK file format, ASCII, VtkImageData)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine saveVTK(data,iStep)

    implicit none

    ! dummy variables
    real   (fp_kind), dimension(hydroParams%isize, hydroParams%jsize, NBVAR), intent(inout) :: data
    integer(int_kind) :: iStep

    ! local variables
    integer :: j,iVar
    integer :: error
    character(LEN=80) :: filename
    character(LEN=8)  :: filenum
    character(1), parameter :: endl  = char(10)  ! end of line
    character(500) :: charBuf
    character(7) :: floatType

    if (useDoublePrecision()) then
       write(floatType,'(a)') 'Float64'
    else
       write(floatType,'(a)') 'Float32'
    end if
    write (filenum,'(i8.8)') iStep
    filename='euler2d_' // filenum // '.vti'

    !open(10,file=filename,status='replace',form='unformatted',action='write',iostat=error)
    !open(10,file=filename,status='replace',access='stream',action='write',iostat=error)
    open(10,file=filename,iostat=error)

    ! write header
    write(10,'(a)') '<?xml version="1.0"?>'//endl
    if (isBigEndian()) then
       write(10,'(a)') '<VTKFile type="ImageData" version="0.1" byte_order="BigEndian">'//endl
    else
       write(10,'(a)') '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//endl
    end if

    ! write mesh extent
    write(charBuf,fmt='(6(I7))',iostat=error) 1,hydroParams%nx+1,1,hydroParams%ny+1,1,2
    write(10,'(a)') repeat(' ',2)//'<ImageData WholeExtent="'//trim(charBuf)//'"'
    write(10,'(a)') ' Origin="0 0 0" Spacing="1 1 1">'//endl
    write(10,'(a)') repeat(' ',2)//'<Piece Extent="'//trim(charBuf)//'">'//endl

    write(10,'(a)') repeat(' ',3)//'<PointData>'//endl
    write(10,'(a)') repeat(' ',4)//'</PointData>'//endl
    write(10,'(a)') repeat(' ',4)//'<CellData>'//endl

    ! write data array (ascii), remove ghost cells
    do iVar=1,NBVAR
       write(10,'(a)') repeat(' ',4)//'<DataArray type="'//trim(floatType)// &
            & '" Name="'//varNames(iVar)//'" format="ascii" >'//endl

       do j=jmin+ghostWidth,jmax-ghostWidth
          write(10,*) data(imin+ghostWidth:imax-ghostWidth,j,iVar)
       end do

       write(10,'(a)') repeat(' ',4)//'</DataArray>'//endl
    end do

    ! write footer
    write(10,'(a)') repeat(' ',4)//'</CellData>'//endl
    write(10,'(a)') repeat(' ',2)//'</Piece>'//endl
    write(10,'(a)') repeat(' ',2)//'</ImageData>'//endl
    write(10,'(a)') '</VTKFile>'//endl

    close(10)

  end subroutine saveVTK

end module m_run
