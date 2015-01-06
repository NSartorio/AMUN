!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2015 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!  This module handles the execution time counting.  Its general
!!  implementation allows to insert up to 128 counters which will measure
!!  time spent on the execution of the bounded code block.  Each timer
!!  can be described.
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
  integer                              , save :: ntimer, norder
  logical          , dimension(ntimers), save :: ftimer, forder
  character(len=32), dimension(ntimers), save :: description
  integer(kind=8)  , dimension(ntimers), save :: times, tstart, tstop
  integer(kind=4)  , dimension(ntimers), save :: tcount, torder
  integer(kind=8)                      , save :: ticks, tbegin
  real   (kind=8)                      , save :: conv    = 1.0d+00

! by default everything is private
!
  private

! declare public subroutines and variables
!
  public :: initialize_timers, finalize_timers
  public :: set_timer, start_timer, stop_timer
  public :: get_timer, get_count, get_timer_total
  public :: ntimers, timer_enabled, timer_description

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

! initialize flags for enabled timers and timer order
!
    ftimer(:)      = .false.
    forder(:)      = .false.

! initialize flag  desciptions
!
    description(:) = ''

! initialize the next available timer and the number of occupied positions
! in the order array
!
    ntimer         = 1
    norder         = 0

! reset timers
!
    times(:)       = 0
    tstart(:)      = 0
    tstop(:)       = 0
    tcount(:)      = 0
    torder(:)      = 0

! prepare the conversion factor
!
    if (ticks > 0) conv = 1.0d+00 / ticks

!-------------------------------------------------------------------------------
!
  end subroutine initialize_timers
!
!===============================================================================
!
! subroutine FINALIZE_TIMERS:
! --------------------------
!
!   Subroutine finalizes module.
!
!
!===============================================================================
!
  subroutine finalize_timers()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!

!-------------------------------------------------------------------------------
!
  end subroutine finalize_timers
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
    if (ntimer <= ntimers) then

! set the timer description
!
      description(ntimer) = trim(adjustl(string))

! set the timer flag
!
      ftimer(ntimer)      = .true.

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
  subroutine start_timer(itimer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: itimer
!
!-------------------------------------------------------------------------------
!
! start accounting the time
!
    call system_clock(tstart(itimer))

! return, if the timer is already allocated in the order array
!
    if (forder(itimer)) return

! otherwise, increase the order position
!
    norder = norder + 1

! assign the current timer with the order position and switch the flag
! signifying that the timer is already in the order array
!
    torder(norder) = itimer
    forder(itimer) = .true.

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
! get the system clock
!
    call system_clock(tstop(timer))

! increase the timer count
!
    tcount(timer) = tcount(timer) + 1

! add the time increment
!
    times(timer) = times(timer) + (tstop(timer) - tstart(timer))

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
    get_timer = max(0.0d+00, conv * times(timer))

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_timer
!
!===============================================================================
!
! function GET_COUNT:
! ------------------
!
!   Function returns the call count for the specified timer.
!
!   Arguments:
!
!     timer - the timer index;
!
!===============================================================================
!
  function get_count(timer)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, intent(in) :: timer

! return variable
!
    integer(kind=4)     :: get_count
!
!-------------------------------------------------------------------------------
!
! estimate the accounted time for the specified timer
!
    get_count = tcount(timer)

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_count
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
    timer_enabled = ftimer(timer)

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
    get_timer_total = max(0.0d+00, conv * (tend - tbegin))

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_timer_total

!===============================================================================
!
end module
