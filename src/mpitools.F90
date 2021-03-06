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
!! module: MPITOOLS
!!
!!  This module provides wrapper subroutines handling the parallel execution
!!  with the Message Passing Interface protocol.
!!
!!
!!******************************************************************************
!
module mpitools

! include external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer

! module variables are not implicit by default
!
  implicit none

! timer indices
!
  integer        , save                 :: imi, imc
#ifdef PROFILE
  integer        , save                 :: imb, imm, ims, imr, ime
#endif /* PROFILE */

! MPI global variables
!
  integer(kind=4), save                 :: comm
  integer(kind=4), save                 :: nproc, nprocs, npmax, npairs
  logical        , save                 :: master = .true.

! allocatable array for processor pairs
!
  integer(kind=4), dimension(:,:), allocatable, save :: pairs

! by default everything is public
!
  public

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!!
!!***  PUBLIC SUBROUTINES  *****************************************************
!!
!===============================================================================
!
! subroutine INITIALIZE_MPITOOLS:
! ------------------------------
!
!   Subroutine initializes the MPITOOLS modules.
!
!
!===============================================================================
!
  subroutine initialize_mpitools()

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
#ifdef MPI
    use mpi            , only : mpi_comm_world, mpi_success
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local variables
!
#ifdef MPI
    integer :: mprocs, i, j, l, n, iret

! allocatable array for processors order
!
    integer(kind=4), dimension(:), allocatable :: procs
#endif /* MPI */

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::initialize_mpitools()'
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
! set timer descriptions
!
    call set_timer('MPI initialization'  , imi)
    call set_timer('MPI communication'   , imc)
#ifdef PROFILE
    call set_timer('mpitools:: broadcast', imb)
    call set_timer('mpitools:: reduce'   , imm)
    call set_timer('mpitools:: send'     , ims)
    call set_timer('mpitools:: receive'  , imr)
    call set_timer('mpitools:: exchange' , ime)
#endif /* PROFILE */

! start time accounting for the MPI initialization
!
    call start_timer(imi)
#endif /* MPI */

! initialize parralel execution parameters
!
    nproc        =  0
    nprocs       =  1
    npmax        =  0

#ifdef MPI
! initialize the MPI interface
!
    call mpi_init(iret)

! check if the MPI interface was initialized successfully
!
    if (iret /= mpi_success) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc)                         &
                               , "The MPI interface could not be initializes!"
    end if

! obtain the total number of processes
!
    call mpi_comm_size(mpi_comm_world, nprocs, iret)

! check if the total number of processes could be obtained
!
    if (iret /= mpi_success) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc)                         &
                               , "The MPI process ID could not be obtained!"
    end if

! obtain the current process identifier
!
    call mpi_comm_rank(mpi_comm_world, nproc , iret)

! check if the process ID was return successfully
!
    if (iret /= mpi_success) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc)                         &
                               , "The MPI process ID could not be obtained!"
    end if

! set the master flag
!
    master = nproc == 0

! calculate the index of the last processor
!
    npmax  = nprocs - 1

! round up the number of processors to even number
!
    mprocs = nprocs + mod(nprocs, 2)

! calculate the number of processor pairs for data exchange
!
    npairs = nprocs * npmax / 2

! allocate space for all processor pairs
!
    allocate(pairs(2 * npairs, 2))

! allocate space for the processor order
!
    allocate(procs(mprocs))

! fill the processor order array
!
    procs(:) = (/(l, l = 0, mprocs - 1)/)

! generate processor pairs
!
    n = 0

! iterate over turns
!
    do l = 1, mprocs - 1

! generate pairs for a given turn
!
      do i = 1, mprocs / 2

! calculate the pair for the current processor
!
        j = mprocs - i + 1

! continue, if the process number is correct (for odd nprocs case)
!
        if (procs(i) < nprocs .and. procs(j) < nprocs) then

! increase the pair number
!
          n = n + 1

! substitute the processor numbers for the current pair
!
          pairs(n,1:2) = (/ procs(i), procs(j) /)

        end if ! max(procs(i), procs(j)) < nprocs

      end do ! i = 1, mprocs / 2

! shift elements in the processor order array
!
      procs(2:mprocs) = cshift(procs(2:mprocs), -1)

    end do ! l = 1, mprocs - 1

! fill out the remaining pairs (swapped)
!
    pairs(npairs+1:2*npairs,1:2) = pairs(1:npairs,2:1:-1)

! allocate space for the processor order
!
    deallocate(procs)

! store the MPI communicator
!
    comm = mpi_comm_world

! stop time accounting for the MPI initialization
!
    call stop_timer(imi)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_mpitools
!
!===============================================================================
!
! subroutine FINALIZE_MPITOOLS:
! ----------------------------
!
!   Subroutine finalizes the MPITOOLS modules.
!
!
!===============================================================================
!
  subroutine finalize_mpitools()

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
#ifdef MPI
    use mpi            , only : mpi_comm_world, mpi_success
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local variables
!
#ifdef MPI
    integer :: iret
#endif /* MPI */

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::finalize_mpitools()'
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
! start time accounting for the MPI initialization
!
    call start_timer(imi)

! initialize the MPI interface
!
    call mpi_finalize(iret)

! check if the MPI interface was finalizes successfully
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! deallocate space used for processor pairs
!
    if (allocated(pairs)) deallocate(pairs)

! stop time accounting for the MPI initialization
!
    call stop_timer(imi)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_mpitools
#ifdef MPI
!
!===============================================================================
!
! subroutine BCAST_INTEGER_VARIABLE:
! ---------------------------------
!
!   Subroutine broadcast an integer variable from the master process to all
!   other processes.
!
!===============================================================================
!
  subroutine bcast_integer_variable(ibuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_integer, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: ibuf
    integer, intent(inout) :: iret

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::bcast_integer_variable()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI broadcast
!
    call start_timer(imb)
#endif /* PROFILE */

    call mpi_bcast(ibuf, 1, mpi_integer, 0, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI broadcast
!
    call stop_timer(imb)
#endif /* PROFILE */

    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine bcast_integer_variable
!
!===============================================================================
!
! subroutine BCAST_REAL_VARIABLE:
! ------------------------------
!
!   Subroutine broadcast a real variable from the master process to all
!   other processes.
!
!===============================================================================
!
  subroutine bcast_real_variable(rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(inout) :: rbuf
    integer     , intent(inout) :: iret

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::bcast_real_variable()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI broadcast
!
    call start_timer(imb)
#endif /* PROFILE */

    call mpi_bcast(rbuf, 1, mpi_real8, 0, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI broadcast
!
    call stop_timer(imb)
#endif /* PROFILE */

    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine bcast_real_variable
!
!===============================================================================
!
! subroutine BCAST_STRING_VARIABLE:
! --------------------------------
!
!   Subroutine broadcast a string variable from the master process to all
!   other processes.
!
!===============================================================================
!
  subroutine bcast_string_variable(sbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_character, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(inout) :: sbuf
    integer         , intent(out)   :: iret

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::bcast_string_variable()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI broadcast
!
    call start_timer(imb)
#endif /* PROFILE */

    call mpi_bcast(sbuf, len(sbuf), mpi_character, 0, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI broadcast
!
    call stop_timer(imb)
#endif /* PROFILE */

    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine bcast_string_variable
!
!===============================================================================
!
! subroutine REDUCE_MINIMUM_INTEGER:
! ---------------------------------
!
!   Subroutine finds the minimum value among the integer values from all
!   processes.
!
!===============================================================================
!
  subroutine reduce_minimum_integer(ibuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_integer, mpi_min, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: ibuf
    integer, intent(out)   :: iret

! local variables
!
    integer                :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_minimum_integer()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_min, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_minimum_integer
!
!===============================================================================
!
! subroutine REDUCE_MINIMUM_REAL:
! ------------------------------
!
!   Subroutine finds the minimum value among the real values from all processes.
!
!===============================================================================
!
  subroutine reduce_minimum_real(rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_min, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(inout) :: rbuf
    integer     , intent(out)   :: iret

! local variables
!
    real(kind=8) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_minimum_real()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_min, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_minimum_real
!
!===============================================================================
!
! subroutine REDUCE_MAXIMUM_INTEGER:
! ---------------------------------
!
!   Subroutine find the maximum value among the integer values from all
!   processes.
!
!===============================================================================
!
  subroutine reduce_maximum_integer(ibuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_integer, mpi_max, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: ibuf
    integer, intent(out)   :: iret

! local variables
!
    integer                :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_maximum_integer()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_max, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_maximum_integer
!
!===============================================================================
!
! subroutine REDUCE_MAXIMUM_REAL:
! ------------------------------
!
!   Subroutine find the maximum value among the values from all processes.
!
!===============================================================================
!
  subroutine reduce_maximum_real(rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_max, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(inout) :: rbuf
    integer     , intent(out)   :: iret

! local variables
!
    real(kind=8) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_maximum_real()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_max, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_maximum_real
!
!===============================================================================
!
! subroutine REDUCE_SUM_INTEGER:
! -----------------------------
!
!   Subroutine finds the sum from all integer values from all processes.
!
!===============================================================================
!
  subroutine reduce_sum_integer(ibuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_integer, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: ibuf
    integer, intent(out)   :: iret

! local variables
!
    integer                :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_sum_integer()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_sum, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_sum_integer
!
!===============================================================================
!
! subroutine REDUCE_SUM_REAL:
! --------------------------
!
!   Subroutine sums the values from all processes.
!
!===============================================================================
!
  subroutine reduce_sum_real(rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(inout) :: rbuf
    integer     , intent(out)   :: iret

! local variables
!
    real(kind=8) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_sum_real()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_sum, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_sum_real
!
!===============================================================================
!
! subroutine REDUCE_MINIMUM_REAL_ARRAY:
! ------------------------------------
!
!   Subroutine find the minimum value for each array element among the
!   corresponding values from all processes.
!
!===============================================================================
!
  subroutine reduce_minimum_real_array(n, rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_min, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(inout) :: rbuf
    integer                   , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_minimum_real_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_min, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_minimum_real_array
!
!===============================================================================
!
! subroutine REDUCE_MAXIMUM_REAL_ARRAY:
! ------------------------------------
!
!   Subroutine find the maximum value for each array element among the
!   corresponding values from all processes.
!
!===============================================================================
!
  subroutine reduce_maximum_real_array(n, rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_max, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(inout) :: rbuf
    integer                   , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_maximum_real_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_max, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_maximum_real_array
!
!===============================================================================
!
! subroutine REDUCE_SUM_INTEGER_ARRAY:
! -----------------------------------
!
!   Subroutine sums the values for each array element from the corresponding
!   values from all processes.
!
!===============================================================================
!
  subroutine reduce_sum_integer_array(n, ibuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_integer, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)    :: n
    integer, dimension(n), intent(inout) :: ibuf
    integer              , intent(out)   :: iret

! local variables
!
    integer, dimension(n)                :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_sum_integer_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(ibuf, tbuf, n, mpi_integer, mpi_sum, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    ibuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_sum_integer_array
!
!===============================================================================
!
! subroutine REDUCE_SUM_REAL_ARRAY:
! --------------------------------
!
!   Subroutine sums the values for each array element from the corresponding
!   values from all processes.
!
!===============================================================================
!
  subroutine reduce_sum_real_array(n, rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(inout) :: rbuf
    integer                   , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n) :: tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_sum_real_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_sum, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_sum_real_array
!
!===============================================================================
!
! subroutine REDUCE_SUM_COMPLEX_ARRAY:
! -----------------------------------
!
!   Subroutine sums the values for each array element from the corresponding
!   complex values from all processes.
!
!===============================================================================
!
  subroutine reduce_sum_complex_array(n, cbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)    :: n
    complex, dimension(n), intent(inout) :: cbuf
    integer              , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n)           :: rbuf, ibuf, tbuf

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::reduce_sum_complex_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI reduce
!
    call start_timer(imm)
#endif /* PROFILE */

    tbuf(:) = real(cbuf(:))
    call mpi_allreduce(tbuf, rbuf, n, mpi_real8, mpi_sum, comm, iret)
    tbuf(:) = aimag(cbuf(:))
    call mpi_allreduce(tbuf, ibuf, n, mpi_real8, mpi_sum, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI reduce
!
    call stop_timer(imm)
#endif /* PROFILE */

! substitute the result
!
    cbuf(1:n) = cmplx(rbuf(1:n), ibuf(1:n))

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', a)") trim(loc), "Operation failed!"
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine reduce_sum_complex_array
!
!===============================================================================
!
! subroutine SEND_REAL_ARRAY:
! --------------------------
!
!   Subroutine sends an arrays of real values to another process.
!
!   Arguments:
!
!     n    - the number of array elements;
!     dst  - the ID of the destination process;
!     tag  - the tag identifying this operation;
!     rbuf - the real array to send;
!     iret - the result flag identifying if the operation was successful;
!
!===============================================================================
!
  subroutine send_real_array(n, dst, tag, rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n, dst, tag
    real(kind=8), dimension(n), intent(in)  :: rbuf
    integer                   , intent(out) :: iret

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::send_real_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI send
!
    call start_timer(ims)
#endif /* PROFILE */

    call mpi_send(rbuf, n, mpi_real8, dst, tag, comm, iret)

#ifdef PROFILE
! stop time accounting for the MPI send
!
    call stop_timer(ims)
#endif /* PROFILE */

! check if the operation was successful
!
    if (iret /= mpi_success .and. master) then
      write(error_unit,"('[', a, ']: ', 2(a, i9))") trim(loc)                  &
                    , "Could not send real array from ", nproc, " to ", dst
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine send_real_array
!
!===============================================================================
!
! subroutine RECEIVE_REAL_ARRAY:
! -----------------------------
!
!   Subroutine receives an arrays of real values from another process.
!
!   Arguments:
!
!     n    - the number of array elements;
!     src  - the7 ID of the source process;
!     tag  - the tag identifying this operation;
!     rbuf - the received real array;
!     iret - the result flag identifying if the operation was successful;
!
!===============================================================================
!
  subroutine receive_real_array(n, src, tag, rbuf, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_success, mpi_status_size

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n, src, tag
    real(kind=8), dimension(n), intent(out) :: rbuf
    integer                   , intent(out) :: iret

! local variables
!
    integer :: status(mpi_status_size)

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::receive_real_array()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI receive
!
    call start_timer(imr)
#endif /* PROFILE */

    call mpi_recv(rbuf, n, mpi_real8, src, tag, comm, status, iret)

#ifdef PROFILE
! stop time accounting for the MPI receive
!
    call stop_timer(imr)
#endif /* PROFILE */

! check if the operation was successful
!
    if (iret /= mpi_success) then
      write(error_unit,"('[', a, ']: ', 2(a, i9))") trim(loc)                  &
                    , "Could not receive real array from ", src, " to ", nproc
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine receive_real_array
!
!===============================================================================
!
! subroutine EXCHANGE_REAL_ARRAYS:
! -------------------------------
!
!   Subroutine exchanges real data buffers between two processes.
!
!   Arguments:
!
!     sproc - the process number to which send the buffer sbuf;
!     stag  - the tag identifying the send operation;
!     ssize - the size of the send buffer sbuf;
!     sbuf  - the real array buffer to send;
!     rproc - the process number from which receive the buffer rbuf;
!     rtag  - the tag identifying the receive operation;
!     rsize - the size of the receive buffer rbuf;
!     rbuf  - the real array buffer to receive;
!     iret  - the result flag identifying if the operation was successful;
!
!===============================================================================
!
  subroutine exchange_real_arrays(sproc, stag, ssize, sbuffer                  &
                                , rproc, rtag, rsize, rbuffer, iret)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit
    use mpi            , only : mpi_real8, mpi_success, mpi_status_size

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                       , intent(in)  :: sproc, rproc
    integer                       , intent(in)  :: stag , rtag
    integer                       , intent(in)  :: ssize, rsize
    real(kind=8), dimension(ssize), intent(in)  :: sbuffer
    real(kind=8), dimension(rsize), intent(in)  :: rbuffer
    integer                       , intent(out) :: iret

! local variables
!
    integer :: status(mpi_status_size)

! local parameters
!
    character(len=*), parameter :: loc = 'MPITOOLS::exchange_real_arrays()'
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

#ifdef PROFILE
! start time accounting for the MPI buffer exchange
!
    call start_timer(ime)
#endif /* PROFILE */

! send sbuf and receive rbuf
!
    call mpi_sendrecv(sbuffer(:), ssize, mpi_real8, sproc, stag                &
                    , rbuffer(:), rsize, mpi_real8, rproc, rtag                &
                                                         , comm, status, iret)

#ifdef PROFILE
! stop time accounting for the MPI buffer exchange
!
    call stop_timer(ime)
#endif /* PROFILE */

! check if the operation was successful
!
    if (iret /= mpi_success) then
      write(error_unit,"('[', a, ']: ', 2(a, i9))") trim(loc)                  &
                             , "Could not exchange real data buffers between " &
                             , sproc, "and", rproc
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine exchange_real_arrays
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
#endif /* MPI */

!===============================================================================
!
end module mpitools
