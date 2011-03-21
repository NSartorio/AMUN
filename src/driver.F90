!!******************************************************************************
!!
!! program: Godunov-AMR
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
program godunov

! modules
!
  use blocks   , only : nleafs
  use config   , only : read_config, nmax, tmax, dtini, dtout, ftype, cfl
#ifdef FORCE
  use config   , only : fdt
#endif /* FORCE */
  use evolution, only : evolve, update_maximum_speed, n, t, dt, dtn
#ifdef FORCE
  use forcing  , only : init_forcing, clear_forcing
#endif /* FORCE */
  use integrals, only : init_integrals, clear_integrals, store_integrals
  use io       , only : write_data
  use mesh     , only : init_mesh, clear_mesh
  use mpitools , only : ncpu, ncpus, init_mpi, clear_mpi, is_master
  use random   , only : init_generator
  use timer    , only : init_timers, start_timer, stop_timer, get_timer        &
                      , get_timer_total
!
!-------------------------------------------------------------------------------
!
! local variables
!
  character(len=60) :: fmt
  integer           :: no, ed, eh, em, es, ec
  real              :: tall, tbeg, tcur, per
!
!-------------------------------------------------------------------------------
!
! initialize MPI
!
  call init_mpi()

! print info message
!
  if (is_master()) then
    write (*,"(1x,78('-'))")
    write (*,"(1x,18('='),4x,a,4x,19('='))") '      Godunov-AMR algorithm      '
    write (*,"(1x,16('='),4x,a,4x,16('='))")                                   &
                                        'Copyright (C) 2008-2011 Grzegorz Kowal'
#ifdef MPI
    write (*,"(1x,18('='),4x,a,i5,a,4x,19('='))") 'MPI enabled with ', ncpus   &
            , ' processors'
#endif /* MPI */
    write (*,"(1x,78('-'))")
    write (*,*)
    write (*,"(1x,a)"         ) "Physics:"
    write (*,"(4x,a,1x,a)"    ) "equations              =",                    &
#ifdef HYDRO
    "HD"
#endif /* HYDRO */
#ifdef MHD
    "MHD"
#endif /* MHD */
    write (*,"(4x,a,1x,a)"    ) "equation of state      =",                    &
#ifdef ADI
    "adiabalic"
#endif /* ADI */
#ifdef ISO
    "isothermal"
#endif /* ISO */

    write (*,*)
  end if

! read configuration file
!
  call read_config()

! initialize timers
!
  call init_timers()

! initialize random number generator
!
  call init_generator()

! initialize our adaptive mesh, refine that mesh to the desired level
! according to the initialized problem
!
  call start_timer(1)
  call init_mesh()
  call stop_timer(1)

! reset number of iterations and time, etc.
!
  n    = 0
  t    = 0.0
  dt   = cfl * dtini
  dtn  = dtini
  no   = 0
  tbeg = 0.0

#ifdef FORCE
! if the forcing time step is larger than the initial time step decrease it
!
  fdt = min(fdt, 0.1d0 * dtini)

! round the initial time step to the integer number of forcing time steps
!
  dt  = fdt * floor(dt / fdt)

! initialize forcing module
!
  call init_forcing()

#endif /* FORCE */
! initialize the integrals module
!
  call init_integrals()

! update the maximum speed in the system
!
  call start_timer(2)
  call update_maximum_speed()
  call stop_timer(2)

! write down the initial state
!
  call start_timer(3)
  call write_data(ftype, no, ncpu)
  call stop_timer(3)

! print information
!
  if (is_master()) then
    write(*,"(1x,a)"   ) "Evolving the system:"
    write(*,'(4x,a4,5x,a4,11x,a2,12x,a6,7x,a3)') 'step', 'time', 'dt'          &
                                                 , 'blocks', 'ETA'
    write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",a1,$)') n, t, dt, nleafs, ed, eh, em, es, char(13)
  end if

! main loop
!
  do while((n .lt. nmax) .and. (t .le. tmax))

! compute new time step
!
    dt = min(cfl * dtn, 2.0d0 * dt)

#ifdef FORCE
! limit the time step to the integer number of forcing time steps
!
    dt = fdt * floor(dt / fdt)
#endif /* FORCE */

! advance the iteration number and time
!
    t  = t + dt
    n  = n + 1

! performe one step evolution
!
    call start_timer(2)
    call evolve()
    call stop_timer(2)

! store integrals
!
    call store_integrals()

! store data
!
    call start_timer(3)
    if (dtout .gt. 0.0 .and. no .lt. (int(t/dtout))) then
      no = no + 1
      call write_data(ftype, no, ncpu)
    end if
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
    if (is_master())                                                           &
      write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //    &
              ',1i2.2,"s",a1,$)') n, t, dt, nleafs, ed, eh, em, es, char(13)
  end do

! add one empty line
!
  if (is_master()) write(*,*)

! write down the final state
!
  call start_timer(3)
  no = no + 1
  call write_data(ftype, no, ncpu)
  call stop_timer(3)

#ifdef FORCE
! finalize forcing module
!
  call clear_forcing()

#endif /* FORCE */
! clear up the integrals module
!
  call clear_integrals()

! deallocate and reset mesh
!
  call start_timer(1)
  call clear_mesh()
  call stop_timer(1)

! get total time
!
  tall = get_timer_total()
  per  = 100.0 / tall

! print info about execution times
!
  if (is_master()) then
    write(fmt,"(a,i2,a)") "(a27,1f", max(1, nint(alog10(tall))) + 5            &
             , ".3,' secs = ',f6.2,' %')"

    write (*,*)
    write(*,"(1x,a)") "Job timings:"
    write (*,fmt) "Initialization        : ", get_timer(1), per * get_timer(1)
    write (*,fmt) "Evolution             : ", get_timer(2), per * get_timer(2)
    write (*,fmt) "Data output           : ", get_timer(3), per * get_timer(3)
    write (*,fmt) "Boundary update       : ", get_timer(4), per * get_timer(4)
    write (*,fmt) "Mesh update           : ", get_timer(5), per * get_timer(5)
#ifdef FORCE
    write (*,fmt) "External forcing      : ", get_timer(10), per * get_timer(10)
    write (*,fmt) " - initialization     : ", get_timer(11), per * get_timer(11)
    write (*,fmt) " - evolution          : ", get_timer(12), per * get_timer(12)
    write (*,fmt) " - real to fourier    : ", get_timer(13), per * get_timer(13)
    write (*,fmt) " - fourier to real    : ", get_timer(14), per * get_timer(14)
#endif /* FORCE */
    write (*,fmt) "EXECUTION TIME        : ", tall        , 100.0
    write (*,*)
  end if

! close access to the MPI
!
  call clear_mpi()

!-------------------------------------------------------------------------------
!
end program
