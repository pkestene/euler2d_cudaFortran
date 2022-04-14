!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module m_utils

  use m_constants
  use m_parameters

contains

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! equation of state (ideal gas)
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(eos)
  pure subroutine eos(rho,eint,p,c,params)

    implicit none

    ! dummy variables
    real(fp_kind), intent(in)    :: rho  !< density
    real(fp_kind), intent(in)    :: eint !< internal energy
    real(fp_kind), intent(inout) :: p    !< pressure
    real(fp_kind), intent(out)   :: c    !< speed of sound
    type(HydroPar), intent(in)   :: params !< hydro parameters

    p = max ( (params%gamma0 - 1.0) * rho * eint, rho * params%smallp)
    c = sqrt( params%gamma0 * p / rho)

  end subroutine eos

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute primitive variables
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(computePrimitives)
  pure subroutine computePrimitives(uLoc,qLoc,c,params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(NBVAR) , intent(in)    :: uLoc ! conservative var
    real(fp_kind), dimension(NBVAR) , intent(out)   :: qLoc ! primitive var
    real(fp_kind)                   , intent(out)   :: c
    type(HydroPar)                  , intent(in)    :: params !< hydro parameters

    ! local variables
    real(fp_kind)                   :: eken ! kinetic energy
    real(fp_kind)                   :: e    ! total energy

    ! compute primitive variables
    qLoc(ID) = max(uLoc(ID), params%smallr)
    qLoc(IU) = uLoc(IU) / qLoc(ID)
    qLoc(IV) = uLoc(IV) / qLoc(ID)

    eken = 0.5 * (qLoc(IU) * qloc(IU) + qLoc(IV) * qLoc(IV))
    e = uLoc(IP) / qLoc(ID) - eken
    ! if (e < 0) then
    !    write(*,*) 'FATAL ERROR : hydro eint < 0  : e ', uLoc(IP), 'eken ',eken,' d ',uLoc(ID),' u ',uLoc(IU),' v ',uLoc(IV)
    !    stop
    ! end if

    ! retrieve pressure from equation of state
    call eos(qLoc(ID), e, qLoc(IP), c, params)

  end subroutine computePrimitives

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute trace hydro 2d
  !! Slopes are compute outside, so slope are input to this
  !! routine
  !!
  !! \param[in]  q  primitive variable state vector
  !! \param[in]  dq primitive variable slopes
  !! \param[in]  dtdx dt divided by dx
  !! \param[in]  dtdy dt divided by dy
  !! \param[out] qm
  !! \param[out] qp
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(trace_unsplit_hydro_2d)
  pure subroutine trace_unsplit_hydro_2d(qLoc, dq, dtdx, dtdy, qm, qp, params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar)   , intent(in)  :: qLoc
    real(fp_kind), dimension(2,nbVar) , intent(in)  :: dq
    real(fp_kind)                     , intent(in)  :: dtdx
    real(fp_kind)                     , intent(in)  :: dtdy
    real(fp_kind), dimension(2, nbVar), intent(out) :: qm
    real(fp_kind), dimension(2, nbVar), intent(out) :: qp
    type(HydroPar)                    , intent(in)  :: params

    ! local variables
    real(fp_kind) :: r,p,u,v
    real(fp_kind) :: sr0,sp0,su0,sv0

    ! TVD slopes in all directions
    real(fp_kind) :: drx,dpx,dux,dvx
    real(fp_kind) :: dry,dpy,duy,dvy

    ! initialize cell centered values
    r = qLoc(ID)
    p = qLoc(IP)
    u = qLoc(IU)
    v = qLoc(IV)

    ! Init cell centered TVD slopes in X direction
    drx = dq(IX,ID);  drx = drx * HALF_F
    dpx = dq(IX,IP);  dpx = dpx * HALF_F
    dux = dq(IX,IU);  dux = dux * HALF_F
    dvx = dq(IX,IV);  dvx = dvx * HALF_F

    ! Init cell centered TVD slopes in Y direction
    dry = dq(IY,ID);  dry = dry * HALF_F
    dpy = dq(IY,IP);  dpy = dpy * HALF_F
    duy = dq(IY,IU);  duy = duy * HALF_F
    dvy = dq(IY,IV);  dvy = dvy * HALF_F


    sr0 = (-u*drx-dux*r)       *dtdx + (-v*dry-dvy*r)       *dtdy
    su0 = (-u*dux-dpx/r)       *dtdx + (-v*duy      )       *dtdy
    sv0 = (-u*dvx      )       *dtdx + (-v*dvy-dpy/r)       *dtdy
    sp0 = (-u*dpx-dux*params%gamma0*p)*dtdx + (-v*dpy-dvy*params%gamma0*p)*dtdy

    ! Update in time the  primitive variables
    r = r + sr0
    u = u + su0
    v = v + sv0
    p = p + sp0

    ! Face averaged right state at left interface
    qp(1,ID) = r - drx
    qp(1,IU) = u - dux
    qp(1,IV) = v - dvx
    qp(1,IP) = p - dpx
    qp(1,ID) = max(params%smallR,  qp(1,ID))
    qp(1,IP) = max(params%smallp * qp(1,ID), qp(1,IP))

    ! Face averaged left state at right interface
    qm(1,ID) = r + drx
    qm(1,IU) = u + dux
    qm(1,IV) = v + dvx
    qm(1,IP) = p + dpx
    qm(1,ID) = max(params%smallR,  qm(1,ID))
    qm(1,IP) = max(params%smallp * qm(1,ID), qm(1,IP))

    ! Face averaged top state at bottom interface
    qp(2,ID) = r - dry
    qp(2,IU) = u - duy
    qp(2,IV) = v - dvy
    qp(2,IP) = p - dpy
    qp(2,ID) = max(params%smallR,  qp(2,ID))
    qp(2,IP) = max(params%smallp * qp(2,ID), qp(2,IP))

    ! Face averaged bottom state at top interface
    qm(2,ID) = r + dry
    qm(2,IU) = u + duy
    qm(2,IV) = v + dvy
    qm(2,IP) = p + dpy
    qm(2,ID) = max(params%smallR,  qm(2,ID))
    qm(2,IP) = max(params%smallp * qm(2,ID), qm(2,IP))

  end subroutine trace_unsplit_hydro_2d

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Trace computations for unsplit Godunov scheme.
  !!
  !! \param[in] q          : Primitive variables state.
  !! \param[in] dqX        : slope along X
  !! \param[in] dqY        : slope along Y
  !! \param[in] dtdx       : dt over dx
  !! \param[in] dtdy       : dt over dy
  !! \param[in] faceId     : which face will be reconstructed
  !! \param[out] qface     : q reconstructed state at cell interface
  !!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(trace_unsplit_2d_along_dir)
  pure subroutine trace_unsplit_2d_along_dir(q, dqX, dqY, dtdx, dtdy, faceId, qface, params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar), intent(in)  :: q
    real(fp_kind), dimension(nbVar), intent(in)  :: dqX, dqY
    real(fp_kind)                  , intent(in)  :: dtdx, dtdy
    integer(int_kind)              , intent(in)  :: faceId
    real(fp_kind), dimension(nbVar), intent(out) :: qface
    type(HydroPar)                 , intent(in)  :: params

    ! local variables
    real(fp_kind) :: r,p,u,v
    real(fp_kind) :: drx,dpx,dux,dvx
    real(fp_kind) :: dry,dpy,duy,dvy
    real(fp_kind) :: sr0,sp0,su0,sv0

    ! Cell centered values
    r =  q(ID);
    p =  q(IP);
    u =  q(IU);
    v =  q(IV);

    ! TVD slopes in all directions
    drx = dqX(ID);
    dpx = dqX(IP);
    dux = dqX(IU);
    dvx = dqX(IV);

    dry = dqY(ID);
    dpy = dqY(IP);
    duy = dqY(IU);
    dvy = dqY(IV);

    ! source terms (with transverse derivatives)
    sr0 = -u*drx-v*dry - (dux+dvy)*r;
    sp0 = -u*dpx-v*dpy - (dux+dvy)*params%gamma0*p;
    su0 = -u*dux-v*duy - (dpx    )/r;
    sv0 = -u*dvx-v*dvy - (dpy    )/r;

    if (faceId == FACE_XMIN) then
       ! Right state at left interface
       qface(ID) = r - HALF_F*drx + sr0*dtdx*HALF_F;
       qface(IP) = p - HALF_F*dpx + sp0*dtdx*HALF_F;
       qface(IU) = u - HALF_F*dux + su0*dtdx*HALF_F;
       qface(IV) = v - HALF_F*dvx + sv0*dtdx*HALF_F;
       qface(ID) = max(params%smallr, qface(ID));
    end if

    if (faceId == FACE_XMAX) then
       ! Left state at right interface
       qface(ID) = r + HALF_F*drx + sr0*dtdx*HALF_F;
       qface(IP) = p + HALF_F*dpx + sp0*dtdx*HALF_F;
       qface(IU) = u + HALF_F*dux + su0*dtdx*HALF_F;
       qface(IV) = v + HALF_F*dvx + sv0*dtdx*HALF_F;
       qface(ID) = max(params%smallr, qface(ID));
    end if

    if (faceId == FACE_YMIN) then
       ! Top state at bottom interface
       qface(ID) = r - HALF_F*dry + sr0*dtdy*HALF_F;
       qface(IP) = p - HALF_F*dpy + sp0*dtdy*HALF_F;
       qface(IU) = u - HALF_F*duy + su0*dtdy*HALF_F;
       qface(IV) = v - HALF_F*dvy + sv0*dtdy*HALF_F;
       qface(ID) = max(params%smallr, qface(ID));
    end if

    if (faceId == FACE_YMAX) then
       ! Bottom state at top interface
       qface(ID) = r + HALF_F*dry + sr0*dtdy*HALF_F;
       qface(IP) = p + HALF_F*dpy + sp0*dtdy*HALF_F;
       qface(IU) = u + HALF_F*duy + su0*dtdy*HALF_F;
       qface(IV) = v + HALF_F*dvy + sv0*dtdy*HALF_F;
       qface(ID) = max(params%smallr, qface(ID));
    end if

  end subroutine trace_unsplit_2d_along_dir

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Compute hydro slopes
  !!
  !! \param[in]  q       : current primitive variable state
  !! \param[in]  qPlusX  : state in the next neighbor cell along XDIR
  !! \param[in]  qMinusX : state in the previous neighbor cell along XDIR
  !! \param[in]  qPlusY  : state in the next neighbor cell along YDIR
  !! \param[in]  qMinusY : state in the previous neighbor cell along YDIR
  !! \param[out] dq      : reference to an array returning the X and Y slopes
  !!
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(slope_unsplit_hydro_2d)
  pure subroutine slope_unsplit_hydro_2d(q, qPlusX, qMinusX, qPlusY, qMinusY, dq, params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar)   , intent(in)  :: q
    real(fp_kind), dimension(nbVar)   , intent(in)  :: qPlusX
    real(fp_kind), dimension(nbVar)   , intent(in)  :: qMinusX
    real(fp_kind), dimension(nbVar)   , intent(in)  :: qPlusY
    real(fp_kind), dimension(nbVar)   , intent(in)  :: qMinusY
    real(fp_kind), dimension(2,nbVar) , intent(out) :: dq
    type(HydroPar)                    , intent(in)  :: params

    ! local variables
    integer :: nvar
    real(fp_kind) ::  dlft, drgt, dcen, dsgn, slop, dlim

    ! alias to slopes vector (1 per dimension)
    real(fp_kind), dimension(nbvar) :: dqX
    real(fp_kind), dimension(nbVar) :: dqY

    if (params%slope_type==0) then

       do nVar=1,nbVar
          dqX(nVar) = ZERO_F
          dqY(nVar) = ZERO_F
       end do

    else if (params%slope_type==1 .or. params%slope_type==2) then ! minmod or average

       do nVar=1,nbVar

          ! slopes in first coordinate direction
          dlft = params%slope_type*(q     (nVar) - qMinusX(nVar))
          drgt = params%slope_type*(qPlusX(nVar) - q      (nVar))
          dcen = HALF_F           *(qPlusX(nVar) - qMinusX(nVar))

          dsgn = -ONE_F
          if (dcen >= ZERO_F) dsgn = ONE_F

          slop = min( abs(dlft), abs(drgt) )
          dlim = slop
          if ( (dlft*drgt) <= ZERO_F ) dlim = ZERO_F

          dqX(nVar) = dsgn * min( dlim, abs(dcen) )

          ! slopes in second coordinate direction
          dlft = params%slope_type*(q     (nVar) - qMinusY(nVar))
          drgt = params%slope_type*(qPlusY(nVar) - q      (nVar))
          dcen = HALF_F           *(qPlusY(nVar) - qMinusY(nVar))

          dsgn = -ONE_F
          if (dcen >= ZERO_F)  dsgn = ONE_F

          slop = min( abs(dlft), abs(drgt) )
          dlim = slop
          if ( (dlft*drgt) <= ZERO_F ) dlim = ZERO_F

          dqY(nVar) = dsgn * min( dlim, abs(dcen) )

       end do ! end for nVar

    end if  ! end params%slope_type == 1 or 2

    dq(1,:) = dqX
    dq(2,:) = dqY

  end subroutine slope_unsplit_hydro_2d

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Riemann solver (compute fluxes at cell
  !! this is actually the approx solver
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(riemann_2d)
  pure subroutine riemann_2d(qleft,qright,qgdnv,flux,params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar), intent(in)  :: qleft
    real(fp_kind), dimension(nbVar), intent(in)  :: qright
    real(fp_kind), dimension(nbVar), intent(out) :: qgdnv
    real(fp_kind), dimension(nbVar), intent(out) :: flux
    type(HydroPar)                 , intent(in)  :: params

    ! local variables
    integer :: iter

    ! Left variables
    real(fp_kind) :: rl,pl,ul

    ! Right variables
    real(fp_kind) :: rr,pr,ur

    ! other var
    real(fp_kind) :: cr, cl, wl, wr, pstar, pold, conv
    real(fp_kind) :: wwl, wwr, ql, qr, usl, usr, delp, sgnm
    real(fp_kind) :: ro, uo, po, wo, co
    real(fp_kind) :: ustar, rstar, cstar, spout, spin, ushock, scr, frac

    ! Pressure, density and velocity
    rl = max(qleft (ID), params%smallr)
    ul =     qleft (IU)
    pl = max(qleft (IP), rl*params%smallp)
    rr = max(qright(ID), params%smallr)
    ur =     qright(IU)
    pr = max(qright(IP), rr*params%smallp)

    ! Lagrangian sound speed
    cl = params%gamma0*pl*rl
    cr = params%gamma0*pr*rr

    ! First guess
    wl = sqrt(cl)
    wr = sqrt(cr)
    pstar = max(((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr), ZERO_F)
    pold = pstar
    conv = ONE_F

    !  Newton-Raphson iterations to find pstar at the required accuracy
    iter = 0
    do while (iter < params%niter_riemann .and. conv > 1e-6)
       wwl = sqrt(cl*(ONE_F+params%gamma6*(pold-pl)/pl))
       wwr = sqrt(cr*(ONE_F+params%gamma6*(pold-pr)/pr))
       ql = 2.0*wwl*wwl*wwl/(wwl*wwl+cl)
       qr = 2.0*wwr*wwr*wwr/(wwr*wwr+cr)
       usl = ul-(pold-pl)/wwl
       usr = ur+(pold-pr)/wwr
       delp = max(qr*ql/(qr+ql)*(usl-usr),-pold)

       pold = pold+delp
       conv = abs(delp/(pold+params%smallpp)) ! Convergence indicator
       iter = iter + 1
    end do

    ! Star region pressure
    ! for a two-shock Riemann problem
    pstar = pold
    wl = sqrt(cl*(ONE_F+params%gamma6*(pstar-pl)/pl))
    wr = sqrt(cr*(ONE_F+params%gamma6*(pstar-pr)/pr))

    ! Star region velocity
    ! for a two shock Riemann problem
    ustar = HALF_F * (ul + (pl-pstar)/wl + ur - (pr-pstar)/wr)

    ! Left going or right going contact wave
    if (ustar >= 0) then
       sgnm = 1.0
    else
       sgnm = -1.0
    end if

    ! Left or right unperturbed state
    if (sgnm > ZERO_F) then
       ro = rl
       uo = ul
       po = pl
       wo = wl
    else
       ro = rr
       uo = ur
       po = pr
       wo = wr
    end if
    co = max(params%smallc, sqrt(abs(params%gamma0*po/ro)))

    ! Star region density (Shock, max prevents vacuum formation in star region)
    rstar = max( ro/(ONE_F+ro*(po-pstar)/(wo*wo)), params%smallr)
    ! Star region sound speed
    cstar = max( params%smallc, sqrt(abs(params%gamma0*pstar/rstar)) )

    ! Compute rarefaction head and tail speed
    spout  = co    - sgnm*uo
    spin   = cstar - sgnm*ustar

    ! Compute shock speed
    ushock = wo/ro - sgnm*uo

    if(pstar >= po) then
       spin  = ushock
       spout = ushock
    end if

    ! Sample the solution at x/t=0
    scr = max(spout-spin, params%smallc+abs(spout+spin))
    frac = HALF_F * (ONE_F + (spout + spin)/scr)

    !if (isnan(frac)) then
    !   frac = ZERO_F
    !else
    !   frac = SATURATE(frac)

    qgdnv(ID) = frac*rstar + (ONE_F-frac)*ro
    qgdnv(IU) = frac*ustar + (ONE_F-frac)*uo
    qgdnv(IP) = frac*pstar + (ONE_F-frac)*po

    if(spout < ZERO_F) then
       qgdnv(ID) = ro
       qgdnv(IU) = uo
       qgdnv(IP) = po
    end if

    if(spin > ZERO_F) then
       qgdnv(ID) = rstar
       qgdnv(IU) = ustar
       qgdnv(IP) = pstar
    end if

    ! transverse velocity
    if(sgnm > ZERO_F) then
       qgdnv(IV) = qleft(IV)
    else
       qgdnv(IV) = qright(IV)
    end if

    call cmpflx(qgdnv, flux, params)

  end subroutine riemann_2d

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! convert godunov state into flux
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(cmpflx)
  pure subroutine cmpflx(qgdnv, flux, params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar), intent(in)  :: qgdnv
    real(fp_kind), dimension(nbVar), intent(out) :: flux
    type(HydroPar)                 , intent(in)   :: params

    ! local variables
    real(fp_kind) :: entho, ekin, etot

    ! Compute fluxes
    ! Mass density
    flux(ID) = qgdnv(ID) * qgdnv(IU)

    ! Normal momentum
    flux(IU) = flux(ID) * qgdnv(IU) + qgdnv(IP)

    ! Transverse momentum 1
    flux(IV) = flux(ID) * qgdnv(IV)

    ! Total energy
    entho = 1.0 / (params%gamma0 - 1.0)
    ekin = HALF_F * qgdnv(ID) * (qgdnv(IU)*qgdnv(IU) + qgdnv(IV)*qgdnv(IV))
    etot = qgdnv(IP) * entho + ekin
    flux(IP) = qgdnv(IU) * (etot + qgdnv(IP))

  end subroutine cmpflx

  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Riemann solver (compute fluxes at cell
  !! this is actually the HLLC solver
  !! TODO : add solver approx
  !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !$acc routine(riemann_2d_hllc)
  pure subroutine riemann_2d_hllc(qleft,qright,qgdnv,flux, params)

    implicit none

    ! dummy variables
    real(fp_kind), dimension(nbVar), intent(in)  :: qleft
    real(fp_kind), dimension(nbVar), intent(in)  :: qright
    real(fp_kind), dimension(nbVar), intent(out) :: qgdnv
    real(fp_kind), dimension(nbVar), intent(out) :: flux
    type(HydroPar)                 , intent(in)  :: params

    ! local variables
    real(fp_kind) :: entho
    real(fp_kind) :: ustar, ptotstar, ro, uo, ptoto, etoto

    ! Left variables
    real(fp_kind) :: rl,pl,ul, ecinl, etotl, ptotl
    real(fp_kind) :: cfastl, SL, rcl, rstarl, etotstarl

    ! Right variables
    real(fp_kind) :: rr,pr,ur, ecinr, etotr, ptotr
    real(fp_kind) :: cfastr, SR, rcr, rstarr, etotstarr

    real(fp_kind) :: rcsum

    ! constant
    entho = ONE_F / (params%gamma0 - ONE_F)

    ! compute left and right variables
    rl = max(qleft(ID), params%smallr)
    pl = max(qleft(IP), rl*params%smallp)
    ul =     qleft(IU)

    ecinl =         HALF_F*rl*ul*ul
    ecinl = ecinl + HALF_F*rl*qleft(IV)*qleft(IV)

    etotl = pl*entho+ecinl
    ptotl = pl

    rr = max(qright(ID), params%smallr)
    pr = max(qright(IP), rr*params%smallp)
    ur =     qright(IU)

    ecinr =         HALF_F*rr*ur*ur
    ecinr = ecinr + HALF_F*rr*qright(IV)*qright(IV)

    etotr = pr*entho+ecinr
    ptotr = pr

    ! Find the largest eigenvalues in the normal direction to the interface
    cfastl = sqrt(max(params%gamma0*pl/rl,params%smallc*params%smallc))
    cfastr = sqrt(max(params%gamma0*pr/rr,params%smallc*params%smallc))

    ! Compute HLL wave speed
    SL = min(ul,ur) - max(cfastl,cfastr)
    SR = max(ul,ur) + max(cfastl,cfastr)

    ! Compute lagrangian sound speed
    rcl = rl*(ul-SL)
    rcr = rr*(SR-ur)
    rcsum = max(rcl+rcr,params%smallr)

    ! Compute acoustic star state
    ustar    = (rcr*ur   +rcl*ul   +  (ptotl-ptotr))/(rcr+rcl)
    ptotstar = (rcr*ptotl+rcl*ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

    ! Left star region variables
    rstarl    = rl*(SL-ul)/(SL-ustar)
    etotstarl = ((SL-ul)*etotl-ptotl*ul+ptotstar*ustar)/(SL-ustar)

    ! Right star region variables
    rstarr    = rr*(SR-ur)/(SR-ustar)
    etotstarr = ((SR-ur)*etotr-ptotr*ur+ptotstar*ustar)/(SR-ustar)

    ! Sample the solution at x/t=0

    if (SL > ZERO_F) then
       ro=rl
       uo=ul
       ptoto=ptotl
       etoto=etotl
    else if (ustar > ZERO_F) then
       ro=rstarl
       uo=ustar
       ptoto=ptotstar
       etoto=etotstarl
    else if (SR > ZERO_F) then
       ro=rstarr
       uo=ustar
       ptoto=ptotstar
       etoto=etotstarr
    else
       ro=rr
       uo=ur
       ptoto=ptotr
       etoto=etotr
    end if

    qgdnv=0

    ! Compute the Godunov flux
    flux(ID) = ro*uo
    flux(IU) = ro*uo*uo+ptoto
    flux(IP) = (etoto+ptoto)*uo
    if (flux(ID) > ZERO_F) then
       flux(IV) = flux(ID)*qleft(IV)
    else
       flux(IV) = flux(ID)*qright(IV)
    end if

  end subroutine riemann_2d_hllc

end module m_utils
