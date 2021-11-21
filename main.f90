!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Author: P. Kestener, CEA Saclay
!!
!!
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use,
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info".
!!
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.

!!!! -*- Mode: F90 -*- !!!!
!> \file main.f90
!> \brief 2D Euler solver on CPU/GPU using nvfortran compiler + stdpar

program euler2d

  use HydroParameters  ! get routines initHydroParameters, printHydroParameters
  use HydroRun         ! get computing routines and utilities (init, boundaries, ...)
  use Monitoring       ! get timer routines

  !!use cudafor

  implicit none

  real   (fp_kind)  :: t=0.0
  real   (fp_kind)  :: dt=0.0
  integer           :: istat
  integer           :: nStep=0

  call initHydroParameters()
  call printHydroParameters()

  ! init domain
  call initHydroRun()
  call compute_dt( dt, modulo(nStep,2) )
  write(*,*) 'Initial value for dt',dt

  ! init boundaries
  call make_boundaries(u)

  ! start computation
  write(*,*) 'Start computation....'
  call timerStart(total_timer)

  ! main loop
  do while (t < tEnd .and. nStep < nStepmax)
     ! output
     if ( modulo(nStep,nOutput) == 0) then
        write(*,*) 'Output results at step ',nStep, 'dt ',dt
        call timerStart(io_timer)
        if (modulo(nStep,2) == 0) then
           call saveVTK(u,nStep)
        else
           call saveVTK(u2,nStep)
        end if
        call timerStop(io_timer)
     end if

     ! compute dt
     call compute_dt( dt, modulo(nStep,2) )

     ! perform one step integration
     call godunov_unsplit(nStep,dt)

     nStep = nStep+1

  end do

  ! write finale results
  !write(*,*) 'Output results at step ',nStep, 'dt ',dt
  !h_u = d_u
  !call saveVTK(h_u,nStep)

  ! end of computation
  call timerStop(total_timer)

  call cleanupHydroRun()

  ! print monitoring
  write(*,*) 'total      time : ', total_timer%elapsed,     'secondes'
  write(*,*) 'compute    time : ', godunov_timer%elapsed,   'secondes', 100*godunov_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'io         time : ', io_timer%elapsed,        'secondes', 100*io_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'boundaries time : ', boundaries_timer%elapsed,'secondes', 100*boundaries_timer%elapsed / total_timer%elapsed, '%'

  write(*,*) 'Perf             : ',nStep*isize*jsize/(total_timer%elapsed-io_timer%elapsed), ' number of cell-updates/s'

end program euler2d
