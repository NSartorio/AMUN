!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!  AMUN is a code to perform numerical simulations in a fluid approximation on
!!  adaptive mesh for Newtonian and relativistic environments with or without
!!  magnetic field.
!!
!!
!!******************************************************************************
!
program amun

! include external subroutines used in this module
!
  use blocks        , only : initialize_blocks, finalize_blocks, get_nleafs
  use boundaries    , only : initialize_boundaries
  use coordinates   , only : initialize_coordinates, finalize_coordinates
  use equations     , only : initialize_equations, finalize_equations
  use evolution     , only : initialize_evolution, finalize_evolution
  use evolution     , only : advance, n, t, dt, dtn, cfl
#ifdef FORCE
  use forcing       , only : initialize_forcing, clear_forcing
#endif /* FORCE */
  use integrals     , only : init_integrals, clear_integrals, store_integrals
  use interpolations, only : initialize_interpolations
  use io            , only : initialize_io, write_data, write_restart_data     &
                           , restart_job
  use mesh          , only : initialize_mesh, clear_mesh
  use mesh          , only : generate_mesh, store_mesh_stats
#ifdef MPI
  use mesh          , only : redistribute_blocks
#endif /* MPI */
  use mpitools      , only : initialize_mpitools, finalize_mpitools, setup_mpi
#ifdef MPI
  use mpitools      , only : bcast_integer_variable
  use mpitools      , only : reduce_maximum_integer
#endif /* MPI */
  use mpitools      , only : master, nprocs, nproc
  use parameters    , only : read_parameters, finalize_parameters
#ifdef MPI
  use parameters    , only : redistribute_parameters
#endif /* MPI */
  use parameters    , only : get_parameter_integer, get_parameter_real         &
                           , get_parameter_string
  use problems      , only : initialize_problems
  use random        , only : initialize_random, finalize_random
  use refinement    , only : initialize_refinement
  use schemes       , only : initialize_schemes, finalize_schemes
  use timers        , only : initialize_timers, finalize_timers
  use timers        , only : start_timer, stop_timer, set_timer, get_timer
  use timers        , only : get_timer_total, timer_enabled, timer_description
  use timers        , only : ntimers

! module variables are not implicit by default
!
  implicit none

! default parameters
!
  integer, dimension(3) :: div = 1
  logical, dimension(3) :: per = .true.
  integer               :: nmax  = 0, ndat = 1, nres = -1, ires = -1
  real                  :: dtini = 1.0d-8
  real                  :: tmax  = 0.0d0, trun = 9999.0d0, tsav = 20.0d0

! temporary variables
!
  character(len=64)     :: lbnd, ubnd

! the termination and status flags
!
  integer               :: iterm, iret

! timer indices
!
  integer               :: iin, iev, itm
#ifdef PROFILE
  integer               :: ipr, ipi
#endif /* PROFILE */

! local snapshot file counters
!
  integer               :: nrun = 1

! iteration and time variables
!
  integer               :: i, ed, eh, em, es, ec
  integer               :: nsteps = 1
  character(len=80)     :: fmt, tmp

  real                  :: tbeg, thrs
  real(kind=8)          :: tm_curr, tm_exec, tm_conv

#ifdef INTEL
! the type of the function SIGNAL should be defined for Intel compiler
!
  integer(kind=4)       :: signal
#endif /* INTEL */

#ifdef SIGNALS
! references to functions handling signals
!
#ifdef GNU
  intrinsic signal
#endif /* GNU */
  external  terminate

! signal definitions
!
  integer, parameter    :: SIGINT = 2, SIGABRT = 6, SIGTERM = 15
#endif /* SIGNALS */

! an array to store execution times
!
  real(kind=8), dimension(ntimers) :: tm

! common block
!
  common /termination/ iterm
!
!-------------------------------------------------------------------------------
!
! initialize the termination flag
!
  iterm = 0

! initialize module TIMERS
!
  call initialize_timers()

! set timer descriptions
!
  call set_timer('INITIALIZATION', iin)
  call set_timer('EVOLUTION'     , iev)
  call set_timer('TERMINATION'   , itm)

! start time accounting for the initialization
!
  call start_timer(iin)

#ifdef SIGNALS
! assign function terminate() with signals
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

! initialize module MPITOOLS
!
  call initialize_mpitools()

! print the welcome message
!
  if (master) then

    write (*,"(1x,78('-'))")
    write (*,"(1x,18('='),17x,a,17x,19('='))") 'A M U N'
    write (*,"(1x,16('='),4x,a,4x,16('='))")                                   &
                                        'Copyright (C) 2008-2013 Grzegorz Kowal'
    write (*,"(1x,18('='),9x,a,9x,19('='))")                                 &
                                        'under GNU GPLv3 license'
    write (*,"(1x,78('-'))")

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
  if (iterm > 0) go to 20

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

! check if the termination flag was broadcaster successfully
!
  if (iterm > 0) go to 20

! reset the termination flag
!
  iterm = 0

! check if the domain is periodic
!
  lbnd = "periodic"
  ubnd = "periodic"
  call get_parameter_string("xlbndry" , lbnd)
  call get_parameter_string("xubndry" , ubnd)
  per(1) = (lbnd == "periodic") .and. (ubnd == "periodic")
  lbnd = "periodic"
  ubnd = "periodic"
  call get_parameter_string("ylbndry" , lbnd)
  call get_parameter_string("yubndry" , ubnd)
  per(2) = (lbnd == "periodic") .and. (ubnd == "periodic")
#ifdef R3D
  lbnd = "periodic"
  ubnd = "periodic"
  call get_parameter_string("zlbndry" , lbnd)
  call get_parameter_string("zubndry" , ubnd)
  per(3) = (lbnd == "periodic") .and. (ubnd == "periodic")
#endif /* R3D */

! get the execution termination parameters
!
  call get_parameter_integer("nmax" , nmax)
  call get_parameter_real   ("tmax" , tmax)
  call get_parameter_real   ("trun" , trun)
  call get_parameter_real   ("tsav" , tsav)

! get the initial time step
!
  call get_parameter_real   ("dtini", dtini)

! get integral calculation interval
!
  call get_parameter_integer("ndat" , ndat)

! get counter and interval for restart snapshots
!
  call get_parameter_integer("nres" , nres)
  call get_parameter_integer("ires" , ires)

! set up the MPI geometry
!
  call setup_mpi(div(:), per(:), .false.)

! initialize the random number generator
!
  call initialize_random(nprocs, nproc)

! initialize geometry modules and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Geometry:"
  end if

! initialize module COORDINATES
!
  call initialize_coordinates(master, iret)

! jump to the end if the equations could not be initialized
!
  if (iret > 0) go to 60

! initialize physics modules and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Physics:"
  end if

! initialize module EQUATIONS
!
  call initialize_equations(master, iret)

! jump to the end if the equations could not be initialized
!
  if (iret > 0) go to 50

! initialize methods modules and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Methods:"
  end if

! initialize evolution
!
  call initialize_evolution(master, iret)

! jump to the end if the schemes could not be initialized
!
  if (iret > 0) go to 40

! initialize module SCHEMES
!
  call initialize_schemes(master, iret)

! jump to the end if the schemes could not be initialized
!
  if (iret > 0) go to 30

! initialize block module
!
  call initialize_blocks()

! initialize module INTERPOLATIONS
!
  call initialize_interpolations()

! initialize boundaries
!
  call initialize_boundaries()

! initialize module PROBLEMS
!
  call initialize_problems()

! initialize module REFINEMENT
!
  call initialize_refinement()

! initialize module IO
!
  call initialize_io()

! reset number of iterations and time, etc.
!
  n    = 0
  t    = 0.0
  dt   = cfl * dtini
  dtn  = dtini

! check if we initiate new problem or restart previous job
!
  if (nres < 0) then

! initialize the mesh module
!
    call initialize_mesh(.true.)

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

  else

! initialize the mesh module
!
    call initialize_mesh(.false.)

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
  call initialize_forcing()
#endif /* FORCE */

#ifdef MPI
! reduce termination flag over all processors
!
  call reduce_maximum_integer(iterm, iret)
#endif /* MPI */

! check if the problem was successfully initialized or restarted
!
  if (iterm > 0) go to 10

! store integrals and data to a file
!
  if (nres < 0) then

    call store_integrals()
    call write_data()

#ifdef MPI
! reduce termination flag over all processors
!
    call reduce_maximum_integer(iterm, iret)
#endif /* MPI */

  end if

! if the initial data were not successfully writen, exit the program
!
  if (iterm > 0) go to 10

! reset the termination flag
!
  iterm = 0
  iret  = 0

! stop time accounting for the initialization
!
  call stop_timer(iin)

! start time accounting for the evolution
!
  call start_timer(iev)

! get current time in seconds
!
  if (master) &
    tbeg = t

! print progress info on master processor
!
  if (master) then

! initialize estimated remaining time of calculations
!
    ed = 9999
    eh =   23
    em =   59
    es =   59

! print progress info
!
    write(*,*)
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

! main loop
!
  do while((nsteps < nmax) .and. (t <= tmax) .and. (iterm == 0))

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
    nsteps = nsteps + 1

! performe one step evolution
!
    call advance()

! store mesh statistics
!
    call store_mesh_stats(n, t)

! store integrals
!
    call store_integrals()

! store data
!
    call write_data()

! get current time in seconds
!
    tm_curr = get_timer_total()

! compute elapsed time
!
    thrs = (tm_curr / 60.0 + tsav) / 60.0

! check if the time exceeds execution time limit
!
    if (thrs >= trun) iterm = 100

! print progress info to console
!
    if (master) then

! calculate days, hours, seconds
!
      ec   = int(tm_curr * (tmax - t) / max(1.0e-8, t - tbeg), kind = 4)
      es   = max(0, int(mod(ec, 60)))
      em   = int(mod(ec / 60, 60))
      eh   = int(ec / 3600)
      ed   = int(eh / 24)
      eh   = int(mod(eh, 24))
      ed   = min(9999,ed)

#if defined INTEL || defined PATHSCALE
      write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //    &
              ',1i2.2,"s",15x,a1,$)')                                          &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#else /* INTEL | PATHSCALE */
      write(*,'(i8,2(1x,1pe14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //    &
              ',1i2.2,"s",15x,a1)',advance="no")                               &
                              n, t, dt, get_nleafs(), ed, eh, em, es, char(13)
#endif /* INTEL | PATHSCALE */

    end if

! prepare iterm
!
    iterm = max(iterm, iret)

#ifdef MPI
! reduce termination flag over all processors
!
    call reduce_maximum_integer(iterm, iret)
#endif /* MPI */

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

! a label to go to if there are any problems, but since all modules have been
! initialized, we have to finalize them first
!
10 continue

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

! deallocate block structure
!
  call finalize_blocks()

! finalize the random number generator
!
  call finalize_random()

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
    if (iterm >= 1 .and. iterm <= 32) then
      write (*,'(a)') ''
      write (*,"(1x,a,i2)") "Terminating program after receiving a" //         &
                                         " termination signal no. ", iterm
      write (*,"(1x,a)") "Restart files have been successfully written."
    end if
    if (iterm == 100) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Terminating program after exceeding the" //          &
                                                            " execution time."
      write (*,"(1x,a)") "Restart files have been successfully written."
    end if
    if (iterm >= 101) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "The initial conditions for the selected problem" //  &
                         " could not be set."
      write (*,"(1x,a)") "Program has been terminated."
    end if
    if (iterm >= 120 .and. iterm < 125) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Problem with restarting job from restart snapshots."
      write (*,"(1x,a)") "Program has been terminated."
    end if
    if (iterm >= 125 .and. iterm < 130) then
      write (*,'(a)') ''
      write (*,"(1x,a)") "Problem with storing snapshots."
      write (*,"(1x,a)") "Program has been terminated."
    end if

! print one empty line
!
    write (*,'(a)') ''

  end if

! finalize module SCHEMES
!
  call finalize_schemes(iret)

! jump point
!
  30 continue

! finalize module EVOLUTION
!
  call finalize_evolution(iret)

! jump point
!
  40 continue

! finalize module EQUATIONS
!
  call finalize_equations(iret)

! jump point
!
  50 continue

! finalize module COORDINATES
!
  call finalize_coordinates(iret)

! jump point
!
  60 continue

! finalize parameters
!
  call finalize_parameters()

! finalize module MPITOOLS
!
  call finalize_mpitools()

! finalize module TIMERS
!
  call finalize_timers()

!-------------------------------------------------------------------------------
!
end program

#ifdef SIGNALS
!
!===============================================================================
!
! function TERMINATE:
! ------------------
!
!   Function sets variable iterm after receiving a signal.
!
!
!===============================================================================
!
integer(kind=4) function terminate(sig_num)

  implicit none

! input arguments
!
  integer(kind=4), intent(in) :: sig_num
  integer                     :: iterm

! common block
!
  common /termination/ iterm

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
