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

module Monitoring

  use HydroPrecision

  ! declare a Timer type
  type Timer
     real(fp_kind) :: elapsed=0.0
     real(fp_kind) :: start=0.0
     real(fp_kind) :: stop=0.0
  end type Timer

  type(Timer) :: total_timer
  type(Timer) :: godunov_timer
  type(Timer) :: boundaries_timer
  type(Timer) :: io_timer

  contains

    ! start timer
    subroutine timerStart(t)
      implicit none
      type(Timer), intent(inout) :: t

      call cpu_time(t%start)

    end subroutine timerStart

    ! stop timer and accumulate timings
    subroutine timerStop(t)
      implicit none
      type(Timer), intent(inout) :: t

      call cpu_time(t%stop)
      t%elapsed = t%elapsed + t%stop - t%start

    end subroutine timerStop

end module Monitoring
