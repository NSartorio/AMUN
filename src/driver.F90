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
!! program: AMUN
!!
!!******************************************************************************
!
program amun

! modules
!
  use blocks   , only : init_blocks, get_nleafs
  use config   , only : read_config
  use config   , only : iterm, nmax, tmax, trun, tsav, dtini, dtout, cfl, nres
#ifdef FORCE
  use config   , only : fdt
#endif /* FORCE */
  use coords   , only : init_coords, clear_coords
  use evolution, only : evolve, find_new_timestep, n, t, dt, dtn
#ifdef FORCE
  use forcing  , only : init_forcing, clear_forcing
#endif /* FORCE */
  use integrals, only : init_integrals, clear_integrals, store_integrals
  use io       , only : nfile, write_data, write_restart_data, restart_job
  use mesh     , only : init_mesh, generate_mesh, store_mesh_stats, clear_mesh
#ifdef MPI
  use mesh     , only : redistribute_blocks
#endif /* MPI */
  use mpitools , only : initialize_mpi, finalize_mpi
#ifdef MPI
  use mpitools      , only : bcast_integer_variable
  use mpitools      , only : reduce_maximum_integer
#endif /* MPI */
  use mpitools , only : master, nprocs
  use parameters    , only : read_parameters, finalize_parameters
#ifdef MPI
  use parameters    , only : redistribute_parameters
#endif /* MPI */
  use random   , only : init_generator
  use timers   , only : initialize_timers, start_timer, stop_timer             &
                      , set_timer, get_timer, get_timer_total                  &
                      , timer_enabled, timer_description, ntimers
!
!-------------------------------------------------------------------------------
!
! local variables
!
  character(len=80) :: fmt, tmp
  integer           :: ed, eh, em, es, ec, ln
  real              :: tall, tcur, tpre, thrs, per

! timer indices
!
  integer           :: nsteps = 1
  integer           :: iin, iev, itm
  real(kind=8)      :: tm_curr, tm_exec, tm_conv

! an array to store execution times
!
  real(kind=8), dimension(ntimers) :: tm

#ifdef SIGNALS
! references to functions handling signals
!
  integer(kind=4) :: iret
#ifdef GNU
  intrinsic signal
#endif /* GNU */
  external  terminate

! signal definitions
!
  integer, parameter :: SIGINT = 2, SIGABRT = 6, SIGTERM = 15
#endif /* SIGNALS */
!
!-------------------------------------------------------------------------------
!
! initialize timers
!
  call initialize_timers()

! set timer descriptions
!
  call set_timer('INITIALIZATION'   , iin)
  call set_timer('EVOLUTION'        , iev)
  call set_timer('TERMINATION'      , itm)

! start time accounting for the initialization
!
  call start_timer(iin)

#ifdef SIGNALS
! assign function terminate with signals
!
#ifdef GNU
  iret = signal(SIGINT , terminate)
  iret = signal(SIGABRT, terminate)
  iret = signal(SIGTERM, terminate)
#else /* GNU */
  iret = signal(SIGINT , terminate, -1)
  iret = signal(SIGABRT, terminate, -1)
  iret = signal(SIGTERM, terminate, -1)
#endif /* GNU */
#endif /* SIGNALS */

! initialize MPI
!
  call initialize_mpi()

! print info message
!
  if (master) then
    write (*,"(1x,78('-'))")
    write (*,"(1x,18('='),4x,a,4x,19('='))") '             A M U N             '
    write (*,"(1x,16('='),4x,a,4x,16('='))")                                   &
                                        'Copyright (C) 2008-2011 Grzegorz Kowal'
#ifdef MPI
    write (*,"(1x,18('='),4x,a,i5,a,4x,19('='))") 'MPI enabled with ', nprocs   &
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

! initialize and read parameters from the parameter file
!
  if (master) call read_parameters(iterm)

#ifdef MPI
! broadcast the termination flag
!
  call bcast_integer_variable(iterm, iret)

! check if the termination flag was broadcaster successfully
!
  if (iterm .gt. 0) go to 20

! reset the termination flag
!
  iterm = 0

! redistribute parameters among all processors
!
  call redistribute_parameters(iterm)

! reduce the termination flag over all processors to check if everything is fine
!
  call reduce_maximum_integer(iterm, iret)
#endif /* MPI */

! read configuration file
!
  call read_config()

! reset number of iterations and time, etc.
!
  n    = 0
  t    = 0.0
  dt   = cfl * dtini
  dtn  = dtini
  ed   = 9999
  eh   = 23
  em   = 59
  es   = 59

! initialize random number generator
!
  call init_generator()

! initialize block module
!
  call init_blocks()

! initialize the coordinate module
!
  call init_coords(master)

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

! find new time step
!
    call find_new_timestep()

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

#ifdef MPI
! redistribute blocks between processors in case the number of processors has
! changed
!
    call redistribute_blocks()
#endif /* MPI */

! find new time step
!
    call find_new_timestep()

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

! stop time accounting for the initialization
!
  call stop_timer(iin)

! start time accounting for the evolution
!
  call start_timer(iev)

! print information
!
  if (master) then
    write(*,*          )
    write(*,"(1x,a)"   ) "Evolving the system:"
    write(*,'(4x,a4,5x,a4,11x,a2,12x,a6,7x,a3)') 'step', 'time', 'dt'          &
                                                 , 'blocks', 'ETA'
#if defined INTEL || defined PATHSCALE
    write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",15x,a1,$)')                                            &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#else /* INTEL | PATHSCALE */
    write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",15x,a1)',advance="no")                                 &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#endif /* INTEL | PATHSCALE */
  end if

! set the previous time needed to estimate the execution time
!
  tpre = get_timer_total()

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
    ec   = int((tmax - t) * (tcur - tpre) / dt, kind=4)
    es   = max(0, int(mod(ec,60)))
    em   = int(mod(ec/60,60))
    eh   = int(ec/3600)
    ed   = int(eh/24)
    eh   = int(mod(eh, 24))
    ed   = min(9999,ed)

! print progress information
!
    if (master)                                                           &
#if defined INTEL || defined PATHSCALE
      write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //    &
              ',1i2.2,"s",15x,a1,$)')                                          &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#else /* INTEL | PATHSCALE */
      write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //    &
              ',1i2.2,"s",15x,a1)',advance="no")                               &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#endif /* INTEL | PATHSCALE */

! obtain the time in hours
!
    thrs = ((2.0 * tcur - tpre) / 60.0 + tsav) / 60.0

! terminate if the (thrs + tsav) > trun
!
    if (thrs .gt. trun) iterm = 1

! update the previous time
!
    tpre = tcur

#if defined SIGNALS && defined MPI
! reduce termination flag over all processors
!
    call reduce_maximum_integer(iterm, iret)
    iterm = iret
#endif /* SIGNALS & MPI */
  end do

! add one empty line
!
  if (master) write(*,*)

! stop time accounting for the evolution
!
  call stop_timer(iev)

! start time accounting for the termination
!
  call start_timer(itm)

! write down the restart dump
!
  call write_restart_data()

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

! finalize the coordinate module
!
  call clear_coords()

! stop time accounting for the termination
!
  call stop_timer(itm)

! get total time
!
  tm(1) = get_timer_total()

! get subtasks timers
!
  do i = 2, ntimers
    tm(i) = get_timer(i)
  end do

! print timings only on master processor
!
  if (master) then

! print one empty line
!
    write (*,'(a)') ''

! calculate the conversion factor
!
    tm_conv = 100.0 / tm(1)

! get the execution time
!
    tm_exec = get_timer_total()

! prepare the string formatting
!
    write (tmp,"(i64)") int(tm(1))
    write (tmp,"(i64)") len_trim(adjustl(tmp)) + 6

! print the execution times
!
    write (fmt,"(a)") "(2x,a32,1x,':',1x,1f" // trim(adjustl(tmp)) //          &
                      ".3,' secs = ',f6.2,' %')"

    write (*,'(1x,a)') 'EXECUTION TIMINGS'
    do i = 2, ntimers
     if (timer_enabled(i)) write (*,fmt) timer_description(i), tm(i)           &
                                                             , tm_conv * tm(i)
    end do

! print the CPU times
!
    write (tmp,"(a)") "(1x,a14,20x,':',1x,1f" // trim(adjustl(tmp)) //         &
                      ".3,' secs = ',f6.2,' %')"
    write (*,tmp) 'TOTAL CPU TIME', tm(1)         , 100.0
    write (*,tmp) 'TIME PER STEP ', tm(1) / nsteps, 100.0 / nsteps
#ifdef MPI
    write (*,tmp) 'TIME PER CPU  ', tm(1) / nprocs, 100.0 / nprocs
#endif /* MPI */

! print the execution time

    write (tmp,"(i64)") int(tm(1))
    write (tmp,"(i64)") len_trim(adjustl(tmp)) + 6
    write (tmp,"(a)") "(1x,a14,20x,':',1x,1f" // trim(adjustl(tmp)) //         &
                      ".3,' secs')"
    write (*,tmp) 'EXECUTION TIME', tm_exec

  end if

! a label to go to if there are any problems
!
20 continue

  if (master) then

! print info about termination due to a signal
!
    if (iterm .ge. 1 .and. iterm .le. 32) then
      write (*,'(a)') ''
      write (*,"(1x,a,i2)") "Terminating program after receiving a" //         &
                                         " termination signal no. ", iterm
      write (*,"(1x,a)") "Restart files have been successfully written."
    end if
    if (iterm .eq. 100) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Terminating program after exceeding the" //          &
                                                            " execution time."
      write (*,"(1x,a)") "Restart files have been successfully written."
    end if
    if (iterm .ge. 101) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "The initial conditions for the selected problem" //  &
                         " could not be set."
      write (*,"(1x,a)") "Program has been terminated."
    end if
    if (iterm .ge. 120 .and. iterm .lt. 125) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Problem with restarting job from restart snapshots."
      write (*,"(1x,a)") "Program has been terminated."
    end if
    if (iterm .ge. 125 .and. iterm .lt. 130) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Problem with storing snapshots."
      write (*,"(1x,a)") "Program has been terminated."
    end if

! print one empty line
!
    write (*,'(a)') ''

  end if

! finalize parameters
!
  call finalize_parameters()

! finalize module mpitools
!
  call finalize_mpi()

!-------------------------------------------------------------------------------
!
end program
#ifdef SIGNALS
!
!===============================================================================
!
! terminate: subroutine sets the iterm variable after receiving a signal
!
!===============================================================================
!
integer(kind=4) function terminate(sig_num)

  use config, only : iterm

  implicit none

! input arguments
!
  integer(kind=4), intent(in) :: sig_num

!-------------------------------------------------------------------------------
!
#ifdef INTEL
  iterm     = sig_num
#else /* INTEL */
  iterm     = 15
#endif /* INTEL */
  terminate = 1

!-------------------------------------------------------------------------------
!
end
#endif /* SIGNALS */

!===============================================================================
!
