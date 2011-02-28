!!******************************************************************************
!!
!! module: timer - handling the timing of the program execution
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of Godunov-AMR.
!!
!!  Godunov-AMR is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Godunov-AMR is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

!===============================================================================
!
end module
