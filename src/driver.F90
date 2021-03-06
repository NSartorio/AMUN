!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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
  use blocks        , only : build_leaf_list
  use boundaries    , only : initialize_boundaries, finalize_boundaries
  use boundaries    , only : boundary_variables
  use coordinates   , only : initialize_coordinates, finalize_coordinates
  use equations     , only : initialize_equations, finalize_equations
  use evolution     , only : initialize_evolution, finalize_evolution
  use evolution     , only : advance, new_time_step
  use evolution     , only : step, time, dt
  use gravity       , only : initialize_gravity, finalize_gravity
  use integrals     , only : initialize_integrals, finalize_integrals
  use integrals     , only : store_integrals
  use interpolations, only : initialize_interpolations, finalize_interpolations
  use io            , only : initialize_io, finalize_io
  use io            , only : restart_from_snapshot
  use io            , only : read_restart_snapshot, write_restart_snapshot
  use io            , only : write_snapshot, next_tout
  use mesh          , only : initialize_mesh, finalize_mesh
  use mesh          , only : generate_mesh, store_mesh_stats
  use mpitools      , only : initialize_mpitools, finalize_mpitools
#ifdef MPI
  use mpitools      , only : bcast_integer_variable
  use mpitools      , only : reduce_maximum_integer, reduce_sum_real_array
#endif /* MPI */
  use mpitools      , only : master, nprocs, nproc
  use operators     , only : initialize_operators, finalize_operators
  use parameters    , only : read_parameters, finalize_parameters
#ifdef MPI
  use parameters    , only : redistribute_parameters
#endif /* MPI */
  use parameters    , only : get_parameter_integer, get_parameter_real         &
                           , get_parameter_string
  use problems      , only : initialize_problems, finalize_problems
  use random        , only : initialize_random, finalize_random
  use refinement    , only : initialize_refinement, finalize_refinement
  use schemes       , only : initialize_schemes, finalize_schemes
  use shapes        , only : initialize_shapes, finalize_shapes
  use sources       , only : initialize_sources, finalize_sources
  use timers        , only : initialize_timers, finalize_timers
  use timers        , only : start_timer, stop_timer, set_timer, get_timer
  use timers        , only : get_timer_total, timer_enabled, timer_description
  use timers        , only : get_count, ntimers
  use user_problem  , only : initialize_user_problem, finalize_user_problem

! module variables are not implicit by default
!
  implicit none

! default parameters
!
  integer               :: nmax  = huge(1), ndat = 1
  real(kind=8)          :: tmax  = 0.0d+00, trun = 9.999d+03, tsav = 3.0d+01
  real(kind=8)          :: dtnext = 0.0d+00

! flag to adjust time precisely to the snapshots
!
  logical        , save :: precise_snapshots = .false.
  character(len=255)    :: prec_snap         = "off"

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

  real(kind=8)          :: tbeg, thrs
  real(kind=8)          :: tm_curr, tm_exec, tm_conv, tm_last = 0.0d+00

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
  integer, parameter    :: SIGERR = -1
  integer, parameter    :: SIGINT =  2, SIGABRT =  6, SIGTERM = 15
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
! initialize the termination and return flags
!
  iterm = 0
  iret  = 0

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
! assign function terminate() to handle signals
!
#ifdef GNU
  if (signal(SIGINT , terminate)     == SIGERR) iret = 1
  if (signal(SIGABRT, terminate)     == SIGERR) iret = 2
  if (signal(SIGTERM, terminate)     == SIGERR) iret = 3
#endif /* GNU */
#ifdef INTEL
  if (signal(SIGINT , terminate, -1) == SIGERR) iret = 1
  if (signal(SIGABRT, terminate, -1) == SIGERR) iret = 2
  if (signal(SIGTERM, terminate, -1) == SIGERR) iret = 3
#endif /* INTEL */

! in the case of problems with signal handler assignment, quit the program
!
  if (iret > 0) go to 200
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
                                        'Copyright (C) 2008-2018 Grzegorz Kowal'
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
  if (iterm > 0) go to 100

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
  if (iterm > 0) go to 100

! reset the termination flag
!
  iterm = 0

! get the execution termination parameters
!
  call get_parameter_integer("nmax" , nmax)
  call get_parameter_real   ("tmax" , tmax)
  call get_parameter_real   ("trun" , trun)
  call get_parameter_real   ("tsav" , tsav)

! correct the run time by the save time
!
  trun   = trun - tsav / 6.0d+01

! initialize dtnext
!
  dtnext = 2.0d+00 * tmax

! get the precise snapshot flag
!
  call get_parameter_string ("precise_snapshots", prec_snap)

! set the precise snapshot flag
!
  if (prec_snap == "on"  ) precise_snapshots = .true.
  if (prec_snap == "ON"  ) precise_snapshots = .true.
  if (prec_snap == "true") precise_snapshots = .true.
  if (prec_snap == "TRUE") precise_snapshots = .true.
  if (prec_snap == "yes" ) precise_snapshots = .true.
  if (prec_snap == "YES" ) precise_snapshots = .true.

! get integral calculation interval
!
  call get_parameter_integer("ndat" , ndat)

! initialize MPI module and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"        ) "Parallelization:"
    write (*,"(4x,a,1x,i6 )" ) "MPI processes          =", nprocs
  end if

! initialize the random number generator (passes the number of OpenMP threads
! and the current thread number)
!
  call initialize_random(1, 0)

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
  if (iret > 0) go to 100

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
  if (iret > 0) go to 60

! print information about the problem
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Problem:"
  end if

! initialize module USER_PROBLEM
!
  call initialize_user_problem(master, iret)

! initialize refinement module and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Refinement:"
  end if

! jump to the end if the refinement could not be initialized
!
  if (iret > 0) go to 50

! initialize module REFINEMENT
!
  call initialize_refinement(master, iret)

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

! initialize module INTERPOLATIONS
!
  call initialize_interpolations(master, iret)

! initialize module OPERATORS
!
  call initialize_operators(master, iret)

! initialize block module
!
  call initialize_blocks(master, iret)

! initialize boundaries module and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Boundaries:"
  end if

! initialize boundaries
!
  call initialize_boundaries(master, iret)

! initialize module PROBLEMS
!
  call initialize_problems(master, iret)

! initialize module SHAPES
!
  call initialize_shapes(master, iret)

! print header for source term information
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"         ) "Source terms:"
  end if

! initialize module GRAVITY
!
  call initialize_gravity(master, iret)

! initialize module SOURCES
!
  call initialize_sources(master, iret)

! initialize boundaries module and print info
!
  if (master) then
    write (*,*)
    write (*,"(1x,a)"              ) "Snapshots:"
    write (*,"(4x,a22,1x,'=',1x,a)") "precise snapshot times", trim(prec_snap)
  end if

! initialize module IO
!
  call initialize_io(master, nrun, iret)

! check if we initiate new problem or restart previous job
!
  if (restart_from_snapshot()) then

! increase the run number
!
    nrun  = nrun + 1

! initialize the mesh module
!
    call initialize_mesh(nrun, master, iret)

! quit if mesh couldn't be initialized
!
    if (iret > 0) go to 10

! reconstruct the meta and data block structures from a given restart file
!
    call read_restart_snapshot(iterm)

#ifdef MPI
! reduce termination flag over all processors
!
  call reduce_maximum_integer(iterm, iret)
#endif /* MPI */

! quit if there was a problem with reading restart snapshots
!
    if (iterm > 0) go to 10

! update the list of leafs
!
    call build_leaf_list()

  else

! initialize the mesh module
!
    call initialize_mesh(nrun, master, iret)

! quit if mesh couldn't be initialized
!
    if (iret > 0) go to 10

! generate the initial mesh, refine that mesh to the desired level according to
! the initialized problem
!
    call generate_mesh()

! update boundaries
!
    call boundary_variables(0.0d+00, dtnext)

! calculate new timestep
!
    call new_time_step(dtnext)

  end if

! initialize the integrals module
!
  call initialize_integrals(master, nrun, iret)

! check if the module was successfully initialized
!
  if (iret > 0) go to 10

! store mesh statistics
!
  call store_mesh_stats(step, time)

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
  if (.not. restart_from_snapshot()) then

    call store_integrals()
    call write_snapshot()

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
    tbeg = time

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
#ifdef INTEL
    write(*,'(i8,2(1x,1es14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",15x,a1,$)')                                            &
                        step, time, dt, get_nleafs(), ed, eh, em, es, char(13)
#else /* INTEL */
    write(*,'(i8,2(1x,1es14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //      &
            ',1i2.2,"s",15x,a1)',advance="no")                                 &
                        step, time, dt, get_nleafs(), ed, eh, em, es, char(13)
#endif /* INTEL */

  end if

! main loop
!
  do while((nsteps <= nmax) .and. (time < tmax) .and. (iterm == 0))

! get the next snapshot time
!
    if (precise_snapshots) dtnext = next_tout() - time

! performe one step evolution
!
    call advance(dtnext)

! advance the iteration number and time
!
    time   = time   + dt
    step   = step   +  1
    nsteps = nsteps +  1

! get current time in seconds
!
    tm_curr = get_timer_total()

! compute elapsed time
!
    thrs = tm_curr / 3.6d+03

! store mesh statistics
!
    call store_mesh_stats(step, time)

! store integrals
!
    call store_integrals()

! write down the restart snapshot
!
    call write_restart_snapshot(thrs, nrun, iret)

! store data
!
    call write_snapshot()

! check if the time exceeds execution time limit
!
    if (thrs > trun) iterm = 100

! print progress info to console, but not too often
!
    if (master) then
      if (time >= tmax .or. (tm_curr - tm_last) >= 1.0d+00) then

! calculate days, hours, seconds
!
        ec   = int(tm_curr * (tmax - time) / max(1.0d-08, time - tbeg), kind = 4)
        es   = max(0, int(mod(ec, 60)))
        em   = int(mod(ec / 60, 60))
        eh   = int(ec / 3600)
        ed   = int(eh / 24)
        eh   = int(mod(eh, 24))
        ed   = min(9999,ed)

#ifdef INTEL
        write(*,'(i8,2(1x,1es14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //  &
                ',1i2.2,"s",15x,a1,$)')                                        &
                        step, time, dt, get_nleafs(), ed, eh, em, es, char(13)
#else /* INTEL */
        write(*,'(i8,2(1x,1es14.6),2x,i8,2x,1i4.1,"d",1i2.2,"h",1i2.2,"m"' //  &
                ',1i2.2,"s",15x,a1)',advance="no")                             &
                        step, time, dt, get_nleafs(), ed, eh, em, es, char(13)
#endif /* INTEL */

! update the timestamp
!
        tm_last = tm_curr

      end if
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

! write down the restart snapshot
!
  call write_restart_snapshot(1.0d+16, nrun, iret)

! a label to go to if there are any problems, but since all modules have been
! initialized, we have to finalize them first
!
10 continue

! start time accounting for the termination
!
  call start_timer(itm)

! finalize integrals module
!
  call finalize_integrals()

! finalize I/O module
!
  call finalize_io(iret)

! finalize module PROBLEMS
!
  call finalize_problems(iret)

! finalize module SOURCES
!
  call finalize_sources(iret)

! finalize module GRAVITY
!
  call finalize_gravity(iret)

! finalize module SHAPES
!
  call finalize_shapes(iret)

! finalize the mesh module
!
  call finalize_mesh(iret)

! finalize boundary module
!
  call finalize_boundaries(iret)

! deallocate block structure
!
  call finalize_blocks(iret)

! finalize the random number generator
!
  call finalize_random()

! finalize the user problem module
!
  call finalize_user_problem(iret)

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

#ifdef MPI
! sum up timers from all processes
!
  call reduce_sum_real_array(ntimers, tm(:), iret)
#endif /* MPI */

! print timings only on master processor
!
  if (master) then

! calculate the conversion factor
!
    tm_conv = 1.0d+02 / tm(1)

! print one empty line
!
    write (*,'(a)') ''

! print the execution times
!
    write (tmp,"(a)") "(2x,a32,1x,':',3x,1f16.3,' secs = ',f6.2,' %')"

    write (*,'(1x,a)') 'EXECUTION TIMINGS'
    do i = 2, ntimers
     if (timer_enabled(i)) then
       if (get_count(i) > 0) then
         write (*,tmp) timer_description(i), tm(i), tm_conv * tm(i)
       end if ! timer counter > 0
     end if ! enabled
    end do

! print the execution times
!
    write (tmp,"(a)") "(1x,a14,20x,':',3x,1f16.3,' secs = ',f6.2,' %')"
    write (*,tmp) 'TOTAL CPU TIME', tm(1)         , 1.0d+02
    write (*,tmp) 'TIME PER STEP ', tm(1) / nsteps, 1.0d+02 / nsteps
#ifdef MPI
    write (*,tmp) 'TIME PER CPU  ', tm(1) / nprocs, 1.0d+02 / nprocs
#endif /* MPI */

! get the execution time
!
    tm_exec = get_timer_total()

! convert the execution time to days, hours, minutes, and seconds and print it
!
    tm(1) = tm_exec / 8.64d+04
    tm(2) = mod(tm_exec / 3.6d+03, 2.4d+01)
    tm(3) = mod(tm_exec / 6.0d+01, 6.0d+01)
    tm(4) = mod(tm_exec, 6.0d+01)
    tm(5) = nint((tm_exec - floor(tm_exec)) * 1.0d+03)
    write (tmp,"(a)") "(1x,a14,20x,':',3x,i14,'d'" //                          &
                                       ",i3.2,'h',i3.2,'m',i3.2,'.',i3.3,'s')"
    write (*,tmp) 'EXECUTION TIME', int(tm(1:5))

  end if

! finalize module OPERATORS
!
  call finalize_operators(iret)

! finalize module INTERPOLATIONS
!
  call finalize_interpolations(iret)

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

! finalize module REFINEMENT
!
  call finalize_refinement(iret)

! jump point
!
  50 continue

! finalize module EQUATIONS
!
  call finalize_equations(iret)

! jump point
!
  60 continue

! finalize module COORDINATES
!
  call finalize_coordinates(iret)

! a label to go to if there are any problems
!
  100 continue

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
    if (iterm >= 101 .and. iterm < 120) then
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

! finalize modules PARAMETERS
!
  call finalize_parameters()

! finalize module MPITOOLS
!
  call finalize_mpitools()

#ifdef SIGNALS
! a label to go to in the case of signal handler problems
!
  200 continue
#endif /* SIGNALS */

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
