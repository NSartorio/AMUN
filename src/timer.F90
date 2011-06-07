!!******************************************************************************
!!
!! module: TIMER - handling the timing of the program execution
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!******************************************************************************
!!
!!  This file is part of the AMUN code.
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
!
module timer

  implicit none

  integer        , parameter                :: ntimers = 32
  integer(kind=8), dimension(ntimers), save :: timers, tstart, tstop
  integer(kind=8)                    , save :: ticks, tbegin

  contains
!
!===============================================================================
!
! init_timers: subroutine initializes timers by creating an array of
!              timers and resetting it
!
!===============================================================================
!
  subroutine init_timers

    implicit none

!
!-------------------------------------------------------------------------------
!
! get current time and the number of ticks per second
!
    call system_clock(count=tbegin, count_rate=ticks)

! reset timers
!
    timers(:) = 0
    tstart(:) = 0
    tstop(:)  = 0

!-------------------------------------------------------------------------------
!
  end subroutine init_timers
!
!===============================================================================
!
! start_timer: subroutine starts timing for a specified timer
!
!===============================================================================
!
  subroutine start_timer(timer)

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
! stop_timer: subroutine stops timing for a specified timer
!
!===============================================================================
!
  subroutine stop_timer(timer)

    implicit none

! input arguments
!
    integer, intent(in) :: timer
!
!-------------------------------------------------------------------------------
!
    call system_clock(tstop(timer))

    timers(timer) = timers(timer) + tstop(timer) - tstart(timer)

!-------------------------------------------------------------------------------
!
  end subroutine stop_timer
!
!===============================================================================
!
! get_timer: function returns the value of specified timer
!
!===============================================================================
!
  function get_timer(timer)

    implicit none

! input arguments
!
    integer, intent(in) :: timer

! output arguments
!
    real                :: get_timer
!
!-------------------------------------------------------------------------------
!
    get_timer = (1.0 * timers(timer)) / ticks

!-------------------------------------------------------------------------------
!
  end function get_timer
!
!===============================================================================
!
! get_timer_total: function returns the value of total execution time
!
!===============================================================================
!
  function get_timer_total()

    implicit none

! output arguments
!
    integer(kind=8)     :: tend
    real                :: get_timer_total
!
!-------------------------------------------------------------------------------
!
! get the current time and the number of ticks per second
!
    call system_clock(count=tend)

    get_timer_total = (1.0 * (tend - tbegin)) / ticks

!-------------------------------------------------------------------------------
!
  end function get_timer_total
!
!===============================================================================
!
! timer_to_time: subroutine converts a real timer to days, hours, minutes and
!                seconds
!
!===============================================================================
!
  subroutine timer_to_time(timer, days, hours, mins, secs, csecs)

    implicit none

! input/output arguments
!
    real           , intent(in)  :: timer
    integer(kind=4), intent(out) :: days, hours, mins, secs, csecs

! local variables
!
    integer(kind=4) :: isecs
!
!-------------------------------------------------------------------------------
!
    isecs = int(timer, kind=4)
    csecs = nint(100 * (timer - real(isecs)))
    secs  = mod(isecs, 60)
    mins  = int(mod(isecs / 60, 60))
    hours = int(isecs / 3600)
    days  = int(hours / 24)
    hours = int(mod(hours, 24))
    days  = min(365, days)

!-------------------------------------------------------------------------------
!
  end subroutine timer_to_time

!===============================================================================
!
end module
