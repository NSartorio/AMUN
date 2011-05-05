!!******************************************************************************
!!
!! program: AMUN
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
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
!
program amun

! modules
!
  use blocks   , only : init_blocks, get_nleafs
  use config   , only : read_config, nmax, tmax, dtini, dtout, cfl, nres
#ifdef FORCE
  use config   , only : fdt
#endif /* FORCE */
  use evolution, only : evolve, update_maximum_speed, n, t, dt, dtn
#ifdef FORCE
  use forcing  , only : init_forcing, clear_forcing
#endif /* FORCE */
  use integrals, only : init_integrals, clear_integrals, store_integrals
  use io       , only : nfile, write_data, write_restart_data, restart_job
  use mesh     , only : init_mesh, generate_mesh, store_mesh_stats, clear_mesh
  use mpitools , only : ncpu, ncpus, init_mpi, clear_mpi, is_master, mfindmaxi
  use random   , only : init_generator
  use timer    , only : init_timers, start_timer, stop_timer, get_timer        &
                      , get_timer_total
!
!-------------------------------------------------------------------------------
!
! local variables
!
  character(len=60) :: fmt
  integer(kind=4)   :: iterm = 0
  integer           :: ed, eh, em, es, ec
  real              :: tall, tbeg, tcur, per
#ifdef SIGNALS

! commons required to share iterm flag
!
  common /signals/ iterm

! references to functions handling signals
!
  intrinsic signal
  external  terminate

! signal definitions
!
  integer, parameter :: SIGINT = 2, SIGABRT = 6, SIGTERM = 15
#endif /* SIGNALS */
!
!-------------------------------------------------------------------------------
!
#ifdef SIGNALS
! assign function terminate with signals
!
  call signal(SIGINT , terminate)
  call signal(SIGABRT, terminate)
  call signal(SIGTERM, terminate)

#endif /* SIGNALS */
! initialize MPI
!
  call init_mpi()

! print info message
!
  if (is_master()) then
    write (*,"(1x,78('-'))")
    write (*,"(1x,18('='),4x,a,4x,19('='))") '             A M U N             '
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
    "adiabatic"
#endif /* ADI */
#ifdef ISO
    "isothermal"
#endif /* ISO */
    write (*,"(4x,a,1x,a)"    ) "geometry               =",                    &
    "rectangular"
  end if

! read configuration file
!
  call read_config()

! reset number of iterations and time, etc.
!
  n    = 0
  t    = 0.0
  dt   = cfl * dtini
  dtn  = dtini
  tbeg = 0.0
  ed   = 9999
  eh   = 23
  em   = 59
  es   = 59

! initialize timers
!
  call init_timers()

! start the initialization timer
!
  call start_timer(1)

! initialize random number generator
!
  call init_generator()

! initialize block module
!
  call init_blocks()

! check if we initiate new problem or restart previous job
!
  if (nres .lt. 0) then

! initialize the mesh module
!
    call init_mesh(.true.)

! initialize the integrals module
!
    call init_integrals(.true.)

! generate the initial mesh, refine that mesh to the desired level according to
! the initialized problem
!
    call generate_mesh()

! store mesh statistics
!
    call store_mesh_stats(n, t)

! update the maximum speed in the system
!
    call update_maximum_speed()

! store integrals
!
    call store_integrals()

! write down the initial state
!
    call write_data()

  else

! initialize the mesh module
!
    call init_mesh(.false.)

! initialize the integrals module
!
    call init_integrals(.false.)

! reconstruct the meta and data block structures from a given restart file
!
    call restart_job()

! reset the start time for the execution time estimate
!
    tbeg = t

  end if

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

! stop the initialization timer
!
  call stop_timer(1)

! print information
!
  if (is_master()) then
    write(*,*          )
    write(*,"(1x,a)"   ) "Evolving the system:"
    write(*,'(4x,a4,5x,a4,11x,a2,12x,a6,7x,a3)') 'step', 'time', 'dt'          &
                                                 , 'blocks', 'ETA'
    write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",a1,$)') n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
  end if

! main loop
!
  do while((n .lt. nmax) .and. (t .le. tmax) .and. (iterm .eq. 0))

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
    call evolve()

! store mesh statistics
!
    call store_mesh_stats(n, t)

! store integrals
!
    call store_integrals()

! store data
!
    if (dtout .gt. 0.0 .and. nfile .lt. (int(t/dtout))) then
      call write_data()
    end if

! get current time in seconds
!
    tcur = get_timer_total()
    ec   = int((tmax - t) / (t - tbeg) * tcur, kind=4)
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
              ',1i2.2,"s",a1,$)') n, t, dt, get_nleafs(), ed, eh, em, es       &
                                , char(13)

#if defined SIGNALS && defined MPI
! reduce termination flag over all processors
!
    call mfindmaxi(iterm)

#endif /* SIGNALS & MPI */
  end do

! add one empty line
!
  if (is_master()) write(*,*)

! write down the restart dump
!
  call write_restart_data()

! start the initialization timer
!
  call start_timer(1)

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
  call clear_mesh()

! get total time
!
  tall = get_timer_total()
  per  = 100.0 / tall

#ifdef SIGNALS
! print info about termination due to a signal
!
  if (is_master()) then
    if (iterm .eq. 1) then
      write(*,*)
      write(*,"(1x,a)") "Program terminated due to received signal."
      write(*,"(1x,a)") "Restart files have been successfully written."
    end if
  end if
#endif /* SIGNALS */

! stop the initialization timer
!
  call stop_timer(1)

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
    write (*,fmt) "Maximum speed estim.  : ", get_timer(6), per * get_timer(6)
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
#ifdef SIGNALS
!
!===============================================================================
!
! terminate: function handles properly the signals sent to the program
!
!===============================================================================
!
function terminate(snum)

  implicit none

! input arguments
!
  integer(4), intent(in) :: snum

! local and shared variables
!
  integer(4)             :: terminate, iterm

! commons to declate shared variables
!
  common /signals/ iterm
!
!-------------------------------------------------------------------------------
!
  iterm     = 1
  terminate = 1

!-------------------------------------------------------------------------------
!
end
#endif /* SIGNALS */

!===============================================================================
!
