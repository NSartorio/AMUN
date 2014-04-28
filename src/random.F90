!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2014 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!! module: RANDOM
!!
!!  This module provides subroutines to random number generators.
!!
!!  References:
!!
!!    [1] Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
!!        random variables', J. Statist. Software, v5(8)
!!
!!******************************************************************************
!
module random

#ifdef PROFILE
! include external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! declare all module variables as implicit
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer            , save :: iri, irc
#endif /* PROFILE */

! random generator type
!
  character(len = 32), save :: gentype = "same"

! variables to store seeds for random generator
!
  integer                                   , save :: kp = 0
  integer                                   , save :: nseeds, lseed
  integer(kind=4), dimension(:), allocatable, save :: seeds

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_random, finalize_random
  public :: nseeds, set_seeds, get_seeds, randomu, randomz, randomn

  contains
!
!===============================================================================
!
! subroutine INITIALIZE_RANDOM:
! ----------------------------
!
!   subroutine initializes random number generator;
!
!===============================================================================
!
  subroutine initialize_random(nprocs, nproc)

! obtain required variables from other modules
!
    use parameters, only : get_parameter_string

! declare all variables as implicit
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: nprocs, nproc

! local variables
!
    integer :: i
    real    :: r
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('random:: initialization'   , iri)
    call set_timer('random:: number generation', irc)

! start accounting time for the random number generator
!
    call start_timer(iri)
#endif /* PROFILE */

! set the processor number
!
    kp = nproc

! obtain the generator type
!
    call get_parameter_string("gentype", gentype)

! calculate the number of seeds
!
    nseeds = 2 * nprocs
    lseed  = nseeds - 1

! allocate seeds for random number generator
!
    allocate(seeds(0:lseed))

! prepare the seeds depending on the type of generator
!
    select case(gentype)
    case('random')
      do i = 0, lseed
        call random_number(r)
        seeds(i) = 123456789 * r
      end do
    case default
      call random_number(r)
      do i = 0, nprocs - 1
        seeds(i) = 123456789 * r
      end do
      call random_number(r)
      do i = nprocs, lseed
        seeds(i) = 123456789 * r
      end do
    end select

#ifdef PROFILE
! stop accounting time for the random number generator
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_random
!
!===============================================================================
!
! subroutine FINALIZE_RANDOM:
! --------------------------
!
!   subroutine releases memory allocated by random number generator variables;
!
!===============================================================================
!
  subroutine finalize_random()

! declare all variables as implicit
!
    implicit none
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generator
!
    call start_timer(iri)
#endif /* PROFILE */

! deallocate seeds if they are allocated
!
    if (allocated(seeds)) deallocate(seeds)

#ifdef PROFILE
! stop accounting time for the random number generator
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_random
!
!===============================================================================
!
! subroutine SET_SEEDS:
! --------------------
!
!   subroutine sets the seeds from the given array;
!
!===============================================================================
!
  subroutine set_seeds(np, seed)

! declare all variables as implicit
!
    implicit none

! input arguments
!
    integer                           , intent(in) :: np
    integer(kind=4), dimension(0:np-1), intent(in) :: seed

! local variables
!
    integer :: i, l
    real    :: r
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generator
!
    call start_timer(iri)
#endif /* PROFILE */

! set the seeds only if the input array and seeds have the same sizes
!
    if (np == nseeds) then

      seeds(0:lseed) = seed(0:lseed)

    else

! if the input array and seeds have different sizes, expand or shrink seeds
!
      select case(gentype)
      case('random')
        l = min(lseed, np - 1)
        seeds(0:l) = seed(0:l)
        if (l < lseed) then
          do i = l + 1, lseed
            call random_number(r)
            seeds(i) = 123456789 * r
          end do
        end if
      case default
        l = nseeds / 2
        do i = 0, l - 1
          seeds(i) = seed(0)
        end do
        do i = l, lseed
          seeds(i) = seed(np-1)
        end do
      end select

    end if

#ifdef PROFILE
! stop accounting time for the random number generator
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine set_seeds
!
!===============================================================================
!
! subroutine GET_SEEDS:
! --------------------
!
!   subroutine returns the seeds through an array;
!
!===============================================================================
!
  subroutine get_seeds(seed)

! declare all variables as implicit
!
    implicit none

! output arguments
!
    integer(kind=4), dimension(0:lseed), intent(out) :: seed
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generator
!
    call start_timer(iri)
#endif /* PROFILE */

    seed(0:lseed) = seeds(0:lseed)

#ifdef PROFILE
! stop accounting time for the random number generator
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine get_seeds
!
!===============================================================================
!
! function RANDOMU:
! ----------------
!
!   function generates uniformly distributed random numbers in range 0.0..1.0;
!
!===============================================================================
!
  function randomu() result(val)

! declare all variables as implicit
!
    implicit none

! output variables
!
    real            :: val

! local variables
!
    integer(kind=4) :: jz, jsr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generation
!
    call start_timer(irc)
#endif /* PROFILE */

    jsr = seeds(kp)
    jz  = jsr

    jsr = ieor( jsr, ishft( jsr,  13 ) )
    jsr = ieor( jsr, ishft( jsr, -17 ) )
    jsr = ieor( jsr, ishft( jsr,   5 ) )

    seeds(kp) = jsr

    val = 0.5 + 0.23283064365e-9 * (jz + jsr)

#ifdef PROFILE
! stop accounting time for the random number generation
!
    call stop_timer(irc)
#endif /* PROFILE */

    return

!-------------------------------------------------------------------------------
!
  end function randomu
!
!===============================================================================
!
! function RANDOMZ:
! ----------------
!
!   function generates uniformly distributed random numbers in range -0.5..0.5;
!
!===============================================================================
!
  function randomz() result(val)

! declare all variables as implicit
!
    implicit none

! output variables
!
    real            :: val

! local variables
!
    integer(kind=4) :: jz, jsr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generation
!
    call start_timer(irc)
#endif /* PROFILE */

    jsr = seeds(kp)
    jz  = jsr

    jsr = ieor( jsr, ishft( jsr,  13 ) )
    jsr = ieor( jsr, ishft( jsr, -17 ) )
    jsr = ieor( jsr, ishft( jsr,   5 ) )

    seeds(kp) = jsr

    val = 0.23283064365e-9 * (jz + jsr)

#ifdef PROFILE
! stop accounting time for the random number generation
!
    call stop_timer(irc)
#endif /* PROFILE */

    return

!-------------------------------------------------------------------------------
!
  end function randomz
!
!===============================================================================
!
! function RANDOMN:
! ----------------
!
!   function generates uniformly distributed random numbers in range -1.0..1.0;
!
!===============================================================================
!
  function randomn() result(val)

! declare all variables as implicit
!
    implicit none

! output variables
!
    real            :: val

! local variables
!
    integer(kind=4) :: jz, jsr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the random number generation
!
    call start_timer(irc)
#endif /* PROFILE */

    jsr = seeds(kp)
    jz  = jsr

    jsr = ieor( jsr, ishft( jsr,  13 ) )
    jsr = ieor( jsr, ishft( jsr, -17 ) )
    jsr = ieor( jsr, ishft( jsr,   5 ) )

    seeds(kp) = jsr

    val = 0.46566128730e-9 * (jz + jsr)

#ifdef PROFILE
! stop accounting time for the random number generation
!
    call stop_timer(irc)
#endif /* PROFILE */

    return

!-------------------------------------------------------------------------------
!
  end function randomn

!===============================================================================
!
end module random
