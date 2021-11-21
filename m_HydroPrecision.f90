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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< REAL and INTEGER precision
!! Note that REAL32 and INT32 are not available in gfortran-4.4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module HydroPrecision

  integer, parameter :: sp_kind   = kind(1.0)  !< single precision
  integer, parameter :: dp_kind   = kind(1.d0) !< double precision
  integer, parameter :: int_kind  = kind(1)

  ! default: use single precison
  integer, parameter :: fp_kind = sp_kind

  ! TODO: use iso_fortran_env module (Fortran 2003)

  contains

    function isBigEndian()
      integer, parameter :: short = selected_int_kind(4)
      integer( short ) :: source = 1_short
      logical :: isBigEndian
      isBigEndian = .false.
      if ( iachar( transfer( source, 'a' ) ) == 0 ) isBigEndian = .true.
      return
    end function isBigEndian

    function useDoublePrecision()
      real(fp_kind) :: tmp
      logical       :: useDoublePrecision
      useDoublePrecision = .false.
      if (kind(tmp) == dp_kind ) useDoublePrecision=.true.
    end function useDoublePrecision

end module HydroPrecision
