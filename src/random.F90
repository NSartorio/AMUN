!!******************************************************************************
!!
!! module: random - handles random number generators by Marsaglia & Tsang
!!
!! references:
!! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
!! random variables', J. Statist. Software, v5(8).
!!
!! Copyright (C) 2000 George Marsaglia, Wai Wan Tsang
!! Copyright (C) 2008 John Burkardt
!! Copyright (C) 2010,2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of AMUN.
!!
!!  This program is free software; you can redistribute it and/or
!!  modify it under the terms of the GNU General Public License
!!  as published by the Free Software Foundation; either version 2
!!  of the License, or (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program; if not, write to the Free Software Foundation,
!!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!!
!!******************************************************************************
!!
module random

  implicit none

  private

  integer                                   , save :: nseeds
  integer(kind=4), dimension(:), allocatable, save :: seeds
  integer(kind=4), dimension(128)           , save :: kn
  real   (kind=4), dimension(128)           , save :: fn, wn

  public :: init_generator, nseeds, get_seeds, set_seeds                       &
          , randomu, randomz, randomn

  contains
!
!===============================================================================
!
! init_generator: subroutine initializes the random number generator
!
!===============================================================================
!
  subroutine init_generator()

#ifdef MPI
    use mpitools, only : is_master, mbcasti
#endif /* MPI */

    implicit none

! local variables
!
    integer(kind=4) :: i
    real(kind=8)    :: dn, tn, q

! parameters
!
    real(kind=8), parameter :: m1 = 2.14748364800000d+09
    real(kind=8), parameter :: vn = 9.91256303526217d-03

! OpenMP functions
!
!$  integer :: omp_get_num_threads, omp_get_thread_num
!
!-------------------------------------------------------------------------------
!
! initialize the seeds required by all subroutines
!
    nseeds = 1
!$omp parallel
!$  nseeds = omp_get_num_threads()
!$omp end parallel

! allocate array for seeds
!
    allocate(seeds(0:nseeds-1))

! fill seeds with random numbers
!
#ifdef MPI
    if (is_master()) then
#endif /* MPI */
      do i = 0, nseeds - 1
        call random_number(q)
        seeds(i) = 123456789 * q
      end do
#ifdef MPI
    end if
#endif /* MPI */

#ifdef MPI
! redistribute seeds
!
    call mbcasti(nseeds, seeds)
#endif /* MPI */

! prepare the arrays used by nonunifor distribution generators
!
    dn = 3.442619855899d+00
    tn = 3.442619855899d+00

    q  = vn / exp(-0.5d+00 * dn * dn)

    kn(1)   = int((dn / q) * m1)
    kn(2)   = 0

    wn(1)   = real(q  / m1, kind=4)
    wn(128) = real(dn / m1, kind=4)

    fn(1)   = 1.0e+00
    fn(128) = real(exp(-0.5d+00 * dn * dn), kind=4)

    do i = 127, 2, -1
      dn      = sqrt(-2.0d+00 * log(vn / dn + exp(-0.5d+00 * dn * dn)))
      kn(i+1) = int((dn / tn) * m1)
      tn      = dn
      fn(i)   = real(exp(-0.5d+00*dn*dn), kind=4)
      wn(i)   = real(dn / m1, kind=4)
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_generator
!
!===============================================================================
!
! get_seeds: subroutine returns the generator seeds
!
!===============================================================================
!
  subroutine get_seeds(seed)

    implicit none

! subroutine arguments
!
    integer(kind=4), dimension(0:nseeds-1), intent(out) :: seed
!
!-------------------------------------------------------------------------------
!
    seed(:) = seeds(:)
!
!-------------------------------------------------------------------------------
!
  end subroutine get_seeds
!
!===============================================================================
!
! set_seeds: subroutine sets the generator seeds
!
!===============================================================================
!
  subroutine set_seeds(seed)

    implicit none

! subroutine arguments
!
    integer(kind=4), dimension(0:nseeds-1), intent(in) :: seed
!
!-------------------------------------------------------------------------------
!
    seeds(:) = seed(:)
!
!-------------------------------------------------------------------------------
!
  end subroutine set_seeds
!
!===============================================================================
!
! shr3: subroutine evaluates the SHR3 generator for integers
!
! input arguments:
!
!   np - the process number (for parallel operations)
!
!===============================================================================
!
  function shr3(np) result(jr)

    implicit none

    integer(kind=4) :: np, jp, jt, jr
!
!-------------------------------------------------------------------------------
!
    jp = seeds(np)
    jt = jp

    jt = ieor(jt, ishft(jt,  13))
    jt = ieor(jt, ishft(jt, -17))
    jt = ieor(jt, ishft(jt,   5))

    seeds(np) = jt

    jr = jp + jt

    return
!
!-------------------------------------------------------------------------------
!
  end function shr3
!
!===============================================================================
!
! randomu: function generates uniformly distributed random numbers from the
!          range 0.0...1.0
!
! input arguments:
!
!   np - the process number (for parallel operations)
!
!===============================================================================
!
  function randomu(np) result(val)

    integer(kind=4) :: np
    real            :: val
!
!-------------------------------------------------------------------------------
!
    val = 0.5 + 0.23283064365e-9 * shr3(np)

    return
!
!-------------------------------------------------------------------------------
!
  end function randomu
!
!===============================================================================
!
! randomz: function generates uniformly distributed random numbers from the
!          range -1.0...1.0
!
! input arguments:
!
!   np - the process number (for parallel operations)
!
!===============================================================================
!
  function randomz(np) result(val)

    integer(kind=4) :: np
    real            :: val
!
!-------------------------------------------------------------------------------
!
    val = 0.46566128730e-9 * shr3(np)

    return
!
!-------------------------------------------------------------------------------
!
  end function randomz
!
!===============================================================================
!
! randomn: function generates normally distributed random numbers
!
! input arguments:
!
!   np - the process number (for parallel operations)
!
!===============================================================================
!
  function randomn(np) result(val)

    implicit none

! input/output arguments
!
    integer(kind=4) :: np
    real            :: val

! local variables
!
    integer(kind=4) :: hz, iz
    real(kind=4)    :: x, y, z, t

! parameters
!
    real(kind=4), parameter :: r = 3.442620e+00
!
!-------------------------------------------------------------------------------
!
    hz = shr3(np)
    iz = iand(hz, 127)

    if (abs(hz) .lt. kn(iz+1)) then
      val = real(hz, kind=4) * wn(iz+1)
    else
      do
        if (iz .eq. 0) then
          do
            x = - 0.2904764e+00 * log(randomu(np))
            y = - log(randomu(np))
            if ((x * x) .le. (y + y)) then
              exit
            end if
          end do

          if (hz .le. 0) then
            val = - r - x
          else
            val = + r + x
          end if

          exit
        end if

        x = real(hz, kind=4) * wn(iz+1)
        z = fn(iz+1) + randomu(np) * (fn(iz) - fn(iz+1))
        t = exp(-0.5e+00 * x * x)

        if (z .lt. t) then
          val = x
          exit
        end if

        hz = shr3(np)
        iz = iand(hz, 127)

        if (abs(hz) .lt. kn(iz+1)) then
          val = real(hz, kind=4) * wn(iz+1)
          exit
        end if
      end do
    end if

    return
!
!-------------------------------------------------------------------------------
!
  end function randomn
!
end module random
