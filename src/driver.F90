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
  use config   , only : read_config, nmax, tmax, dtini, dtout, ftype, cfl
  use evolution, only : evolve, n, t, dt, dtn
  use io       , only : write_data
  use mesh     , only : init_mesh, clear_mesh
  use timer    , only : init_timers, start_timer, stop_timer          &
                      , get_timer, get_timer_total
!
!-------------------------------------------------------------------------------
!
! local variables
!
  character(len=60) :: fmt
  integer           :: no, ed, eh, em, es, ec
  real              :: tall, tbeg, tcur
!
!-------------------------------------------------------------------------------
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

! reset number of iterations and time, etc.
!
  n   = 0
  t   = 0.0
  dt  = cfl*dtini
  dtn = dtini
  no  = 0
  tbeg = 0.0

! initialize our adaptive mesh, refine that mesh to the desired level
! according to the initialized problem
!
  call start_timer(1)
  call init_mesh
  call stop_timer(1)

! write down the initial state
!
  call start_timer(3)
  call write_data(ftype, no, 0)
  call stop_timer(3)

! TODO: main loop, perform one step evolution of the system, do refinement/derefinement
! TODO: get new time step, dump data, print info about the progress

! print information
!
  write(*,*)
  write(*,"(1x,a)"   ) "Evolving system:"

! main loop
!
  do while((n .lt. nmax) .and. (t .le. tmax))

! compute new time step
!
    dt = min(cfl*dtn, 2.*dt)

! advance the iteration number and time
!
    t  = t + dt
    n  = n + 1

! performe one step evolution
!
    call start_timer(2)
    call evolve
    call stop_timer(2)

! store data
!
    call start_timer(3)
    if (dtout .gt. 0.0 .and. no .lt. (int(t/dtout))) then
      no = no + 1
      call write_data(ftype, no, 0)
    endif
    call stop_timer(3)

! get current time in seconds
!
    tcur = get_timer_total()
    ec   = int((tmax - t)/(t - tbeg)*tcur, kind=4)
    es   = max(0, int(mod(ec,60)))
    em   = int(mod(ec/60,60))
    eh   = int(ec/3600)
    ed   = int(eh/24)
    eh   = int(mod(eh, 24))
    ed   = min(9999,ed)

! print progress information
!
    if (mod(n, 50) .eq. 1) &
      write(*,'(4x,a4,3(3x,a9,3x),4x,a12)') 'iter', 'time ', 'dt ', 'dtnew ', 'remain. time'
    write(*,'(i8,3(1x,1pe14.6),2x,1i4.1,"d",1i2.2,"h",1i2.2,"m",1i2.2,"s")') n, t, dt, dtn, ed, eh, em, es

  end do

! write down the final state
!
  call start_timer(3)
  no = no + 1
  call write_data(ftype, no, 0)
  call stop_timer(3)

! deallocate and reset mesh
!
  call start_timer(1)
  call clear_mesh
  call stop_timer(1)

! get total time
!
  tall = get_timer_total()

! print info about execution times
!
  write(fmt,"(a,i2,a)") "(a27,1f", max(1, nint(alog10(tall))) + 6, ".4,' secs = ',f7.3,' %')"

  write (*,*)
  write (*,fmt) "Time for initialization : ", get_timer(1), 100.0*get_timer(1)/tall
  write (*,fmt) "Time for evolution      : ", get_timer(2), 100.0*get_timer(2)/tall
  write (*,fmt) "Time for data output    : ", get_timer(3), 100.0*get_timer(3)/tall
  write (*,fmt) "Time for boundary update: ", get_timer(4), 100.0*get_timer(4)/tall
  write (*,fmt) "Time for mesh update    : ", get_timer(5), 100.0*get_timer(5)/tall
  write (*,fmt) "EXECUTION TIME          : ", tall, 100.0

!-------------------------------------------------------------------------------
!
end program
