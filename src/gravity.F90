!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2017 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: GRAVITY
!!
!!  This modules handles the calculation of gravitational acceleration, static
!!  or time dependent.
!!
!!******************************************************************************
!
module gravity

#ifdef PROFILE
! include external procedures
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer, save :: imi, imc
#endif /* PROFILE */

! flag indicating if the gravitational source term is enabled
!
  logical, save :: gravity_enabled = .false.

! pointer to the gravitational acceleration subroutine
!
  procedure(gacc_none), pointer, save :: gravitational_acceleration => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_gravity, finalize_gravity
  public :: gravitational_acceleration
  public :: gravity_enabled

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!!
!!***  PUBLIC SUBROUTINES  *****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine INITIALIZE_GRAVITY:
! -----------------------------
!
!   Subroutine initializes module GRAVITY.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_gravity(verbose, iret)

! include external procedures and variables
!
    use parameters     , only : get_parameter_string
    use user_problem   , only : gravitational_acceleration_user                &
                              , gravity_enabled_user

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64) :: problem_name   = "none"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('gravity:: initialize'  , imi)
    call set_timer('gravity:: acceleration', imc)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the problem name
!
    call get_parameter_string("problem", problem_name)

! select the shape update subroutine depending on the problem
!
    select case(trim(problem_name))

    case("rt", "rayleightaylor", "rayleigh-taylor")
      gravitational_acceleration => gacc_rayleigh_taylor
      gravity_enabled = .true.

    case default

! in case of other problems, gravity is calculated by user
!
      gravitational_acceleration => gravitational_acceleration_user
      gravity_enabled            =  gravity_enabled_user

    end select

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_gravity
!
!===============================================================================
!
! subroutine FINALIZE_GRAVITY:
! ---------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_gravity(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! nullify procedure pointers
!
    nullify(gravitational_acceleration)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_gravity
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine GACC_NONE:
! --------------------
!
!   Subroutine does nothing, but it is required to define the interface for
!   other specific gracitational acceleration subroutines.
!
!   Arguments:
!
!     t, dt   - time and the time increment;
!     x, y, z - rectangular coordinates;
!     gacc    - vector of the gravitational acceleration;
!
!===============================================================================
!
  subroutine gacc_none(t, dt, x, y, z, gacc)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)              , intent(in)  :: t, dt
    real(kind=8)              , intent(in)  :: x, y, z
    real(kind=8), dimension(3), intent(out) :: gacc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the gravitational acceleration calculation
!
    call start_timer(imc)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the gravitational acceleration calculation
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine gacc_none
!
!===============================================================================
!
! subroutine GACC_RAYLEIGH_TAYLOR:
! -------------------------------
!
!   Subroutine returns the gravitational acceleration for the Rayleigh-Taylor
!   instability problem.
!
!   Arguments:
!
!     t, dt   - time and the time increment;
!     x, y, z - rectangular coordinates;
!     gacc    - vector of the gravitational acceleration;
!
!===============================================================================
!
  subroutine gacc_rayleigh_taylor(t, dt, x, y, z, gacc)

! include external procedures and variables
!
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)              , intent(in)  :: t, dt
    real(kind=8)              , intent(in)  :: x, y, z
    real(kind=8), dimension(3), intent(out) :: gacc

! gravitational acceleration constant
!
    logical     , save :: first      = .true.
    real(kind=8), save :: gacc_const = -1.0d-01
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the gravitational acceleration calculation
!
    call start_timer(imc)
#endif /* PROFILE */

! read problem parameters during the first execution
!
    if (first) then

      call get_parameter_real("gacc", gacc_const)

      first = .false.

    end if

! calculate gravitational acceleration components
!
    gacc(1) = 0.0d+00
    gacc(2) = gacc_const
    gacc(3) = 0.0d+00

#ifdef PROFILE
! stop accounting time for the gravitational acceleration calculation
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine gacc_rayleigh_taylor

!===============================================================================
!
end module gravity
