!!*****************************************************************************
!!
!! program: Godunov-AMR
!!
!! Copyright (C) 2008 Grzegorz Kowal <kowal@astro.wisc.edu>
!!
!!*****************************************************************************
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
!!*****************************************************************************
!!
!
program godunov

! modules
!
  use config, only : read_config
  use io    , only : write_data
  use mesh  , only : init_mesh, clear_mesh
  use timer , only : init_timers, start_timer, stop_timer, get_timer &
                   , get_timer_total
!
!----------------------------------------------------------------------
!
! local variables
!
  character(len=60) :: fmt
  real              :: tall
!
!----------------------------------------------------------------------
!
! print info message
!
  write (*,"(1x,78('-'))")
  write (*,"(1x,18('='),4x,a,4x,19('='))") '      Godunov-AMR algorithm      '
  write (*,"(1x,18('='),4x,a,4x,19('='))") 'Copyright (C) 2008 Grzegorz Kowal'
  write (*,"(1x,78('-'))")
  write (*,*)

! read configuration file
!
  call read_config

! initialize timers
!
  call init_timers

! initialize our adaptive mesh, refine that mesh to the desired level
! according to the initialized problem
!
  call start_timer(1)
  call init_mesh
  call stop_timer(1)

! write down the initial state
!
  call start_timer(2)
  call write_data('r', 0, 0)
  call stop_timer(2)

! TODO: main loop, perform one step evolution of the system, do refinement/derefinement
! TODO: get new time step, dump data, print info about the progress

! deallocate and reset mesh
!
  call clear_mesh

! get total time
!
  tall = get_timer_total()

! print info about execution times
!
  write(fmt,"(a,i2,a)") "(a27,1f", max(1, nint(alog10(tall))) + 6, ".4,' secs = ',f7.3,' %')"

  write (*,*)
  write (*,fmt) "Time for initialization : ", get_timer(1), 100.0*get_timer(1)/tall
  write (*,fmt) "Time for data output    : ", get_timer(2), 100.0*get_timer(2)/tall
  write (*,fmt) "EXECUTION TIME          : ", tall, 100.0

!----------------------------------------------------------------------
!
end program
