!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2012 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: TIMERS
!!
!!  This module handles the execution time counting for different parts of the
!!  program.
!!
!!******************************************************************************
!
module timers

! module variables are not implicit by default
!
  implicit none

! module variables
!
  integer          , parameter                :: ntimers = 128
  integer                              , save :: ntimer
  logical          , dimension(ntimers), save :: flag
  character(len=32), dimension(ntimers), save :: description
  integer(kind=8)  , dimension(ntimers), save :: times, tstart, tstop
  integer(kind=8)                      , save :: ticks, tbegin
  real   (kind=8)                      , save :: conv

! by default everything is private
!
  private

! declare public subroutines and variables
!
  public :: initialize_timers, set_timer, start_timer, stop_timer, get_timer   &
          , get_timer_total, ntimers, timer_enabled, timer_description

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_TIMERS:
! ----------------------------
!
!   Subroutine initializes module TIMERS and allocates memory to store
!   execution times.
!
!===============================================================================
!
  subroutine initialize_timers()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the current time and the number of ticks per second
!
    call system_clock(count=tbegin, count_rate=ticks)

! initialize the flags signifying that the counter is used
!
    flag(:)   = .false.

! initialize flag  desciptions
!
    description(:) = ''

! initialize the next available timer
!
    ntimer    = 1

! reset timers
!
    times(:)  = 0
    tstart(:) = 0
    tstop(:)  = 0

! prepare the conversion factor
!
    conv = 1.0d0 / ticks

!-------------------------------------------------------------------------------
!
  end subroutine initialize_timers
!
!===============================================================================
!
! subroutine SET_TIMER:
! --------------------
!
!   Subroutine sets the timer by giving its desciption.
!
!   Arguments:
!
!     string - the timer description;
!     timer  - the timer index;
!
!===============================================================================
!
  subroutine set_timer(string, timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: timer
!
!-------------------------------------------------------------------------------
!
! increase the timer count
!
    ntimer = ntimer + 1

! check if the timer didn't exceed the number of avalaible timers
!
    if (ntimer .le. ntimers) then

! set the timer description
!
      description(ntimer) = trim(adjustl(string))

! set the timer flag
!
      flag(ntimer)        = .true.

! return the timer index
!
      timer               = ntimer

    end if

!-------------------------------------------------------------------------------
!
  end subroutine set_timer
!
!===============================================================================
!
! subroutine START_TIMER:
! ----------------------
!
!   Subroutine starts accounting time for a specified timer.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  subroutine start_timer(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer
!
!-------------------------------------------------------------------------------
!
    call system_clock(tstart(timer))

!-------------------------------------------------------------------------------
!
  end subroutine start_timer
!
!===============================================================================
!
! subroutine STOP_TIMER:
! ---------------------
!
!   Subroutine stops accounting time of the specified timer.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  subroutine stop_timer(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer
!
!-------------------------------------------------------------------------------
!
    call system_clock(tstop(timer))

    times(timer) = times(timer) + tstop(timer) - tstart(timer)

!-------------------------------------------------------------------------------
!
  end subroutine stop_timer
!
!===============================================================================
!
! function GET_TIMER:
! ------------------
!
!   Function returns accounted time of the specified timer.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  function get_timer(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer

! return variable
!
    real(kind=8)        :: get_timer
!
!-------------------------------------------------------------------------------
!
! estimate the accounted time for the specified timer
!
    get_timer = conv * max(0, times(timer))

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_timer
!
!===============================================================================
!
! function TIMER_ENABLED:
! ----------------------
!
!   Function returns logical flag determining if the specified timer is enabled.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  function timer_enabled(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer

! return variable
!
    logical             :: timer_enabled
!
!-------------------------------------------------------------------------------
!
! estimate the accounted time for the specified timer
!
    timer_enabled = flag(timer)

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function timer_enabled
!
!===============================================================================
!
! function TIMER_DESCRIPTION:
! --------------------------
!
!   Function returns the the specified timer description.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  function timer_description(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer

! return variable
!
    character(len=32)   :: timer_description
!
!-------------------------------------------------------------------------------
!
! estimate the accounted time for the specified timer
!
    timer_description = description(timer)

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function timer_description
!
!===============================================================================
!
! function GET_TIMER_TOTAL:
! ------------------------
!
!   Function returns the total execution time.
!
!===============================================================================
!
  function get_timer_total()

! local variables are not implicit by default
!
    implicit none

! return value
!
    real(kind=8)        :: get_timer_total

! local variables
!
    integer(kind=8)     :: tend
!
!-------------------------------------------------------------------------------
!
! get the current time
!
    call system_clock(count=tend)

! estimate the total execution time
!
    get_timer_total = conv * max(1, tend - tbegin)

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_timer_total

!===============================================================================
!
end module
