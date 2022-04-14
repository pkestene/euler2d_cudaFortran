!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Author: P. Kestener, CEA
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

!
! nvfortran -O3 -stdpar=gpu -gpu=cc80,cuda11.6 -acc=gpu -gpu=rdc,managed -Minform=warn -Minfo m_precision.f90 m_constants.f90 m_parameters.f90 m_monitoring.f90 m_utils.f90 m_run.f90 main.f90 -o euler2d
!


!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program euler2d

  !use m_constants
  !use m_utils
  use m_run
  use m_monitoring

  implicit none

  real   (fp_kind)  :: t=0.0
  real   (fp_kind)  :: dt=0.0
  integer           :: istat
  integer           :: nStep=0

  ! read namelist and init parameters
  call initHydroParameters()
  call printHydroParameters()

  ! init domain
  call run_init()

  call compute_dt( hydroParams, dt, modulo(nStep,2) )
  write(*,*) 'Initial value for dt',dt

  ! init boundaries
  call make_boundaries(hydroParams, u)

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
     call compute_dt( hydroParams, dt, modulo(nStep,2) )

     ! perform one step integration
     call godunov_unsplit(nStep,dt)

     nStep = nStep+1

  end do


  ! write final results
  write(*,*) 'Output results at step ',nStep, 'dt ',dt
  call saveVTK(u,nStep)

  ! end of computation
  call timerStop(total_timer)

  call run_cleanup()

  ! print monitoring
  write(*,*) 'total      time : ', total_timer%elapsed,     'secondes'
  write(*,*) 'compute    time : ', godunov_timer%elapsed,   'secondes', 100*godunov_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'io         time : ', io_timer%elapsed,        'secondes', 100*io_timer%elapsed / total_timer%elapsed, '%'
  write(*,*) 'boundaries time : ', boundaries_timer%elapsed,'secondes', 100*boundaries_timer%elapsed / total_timer%elapsed, '%'

  write(*,*) 'Perf             : ',nStep*hydroParams%isize*hydroParams%jsize/(total_timer%elapsed-io_timer%elapsed)*1e-6, ' number of Mcell-updates/s'

end program euler2d
