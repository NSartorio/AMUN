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
!! module: MPITOOLS
!!
!!  This module provides wrapper subroutines handling the parallel execution
!!  with the Message Passing Interface protocol.
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

! MPI global variables
!
  integer(kind=4), save                 :: comm3d
  integer(kind=4), save                 :: nproc, nprocs
  integer(kind=4), save, dimension(3)   :: pdims, pcoords, pparity
  integer(kind=4), save, dimension(3,2) :: pneighs
  logical        , save, dimension(3)   :: periodic
  logical        , save                 :: master = .true.

! by default everything is public
!
  public

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_MPI:
! -------------------------
!
!   Subroutine initializes the MPITOOLS modules.
!
!===============================================================================
!
  subroutine initialize_mpi()

! include external procedures and variables
!
#ifdef MPI
    use mpi, only : mpi_comm_world, mpi_success
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local variables
!
#ifdef MPI
    integer :: iret
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
! set timer descriptions
!
    call set_timer('MPI initialization', imi)
    call set_timer('MPI communication' , imc)

! start time accounting for the MPI initialization
!
    call start_timer(imi)
#endif /* MPI */

! initialize parralel execution parameters
!
    nproc        =  0
    nprocs       =  1
    pdims(:)     =  1
    pcoords(:)   =  0
    pparity(:)   =  0
    pneighs(:,:) = -1
    periodic(:)  = .false.

#ifdef MPI
! initialize the MPI interface
!
    call mpi_init(iret)

! check if the MPI interface was initialized successfully
!
    if (iret .ne. mpi_success) then
      write(*,*) 'The MPI interface could not be initializes! Exiting...'
      write(*,*)
      stop
    end if

! obtain the total number of processes
!
    call mpi_comm_size(mpi_comm_world, nprocs, iret)

! check if the total number of processes could be obtained
!
    if (iret .ne. mpi_success) then
      write(*,*) 'The MPI process ID could not be obtained! Exiting...'
      write(*,*)
      stop
    end if

! obtain the current process identificator
!
    call mpi_comm_rank(mpi_comm_world, nproc , iret)

! check if the process ID was return successfully
!
    if (iret .ne. mpi_success) then
      write(*,*) 'The MPI process ID could not be obtained! Exiting...'
      write(*,*)
      stop
    end if

! set the master flag
!
    master = (nproc .eq. 0)

! store the MPI pool handles
!
    comm3d = mpi_comm_world

! stop time accounting for the MPI initialization
!
    call stop_timer(imi)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_mpi
!
!===============================================================================
!
! subroutine FINALIZE_MPI:
! -----------------------
!
!   Subroutine finalizes the MPITOOLS modules.
!
!===============================================================================
!
  subroutine finalize_mpi()

! include external procedures and variables
!
#ifdef MPI
    use mpi, only : mpi_comm_world, mpi_success
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local variables
!
#ifdef MPI
    integer :: iret
#endif /* MPI */
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
    if (iret .ne. mpi_success) then
      if (master) then
        write(*,*) 'The MPI interface could not be finalized! Exiting...'
        write(*,*)
      end if
      stop
    end if

! stop time accounting for the MPI initialization
!
    call stop_timer(imi)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_mpi
!
!===============================================================================
!
! subroutine SETUP_MPI:
! --------------------
!
!   Subroutine sets the MPI geometry.
!
!===============================================================================
!
  subroutine setup_mpi(div, per)

! include external procedures and variables
!
#ifdef MPI
    use mpi, only : mpi_comm_world, mpi_success
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer, dimension(3), intent(in) :: div
    logical, dimension(3), intent(in) :: per

! local variables
!
    integer :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
! start time accounting for the MPI initialization
!
    call start_timer(imi)

! check if the total number of chunks in division corresponds to the number of
! processes, if not try to find the best division
!
    if (nprocs .ne. product(div(:))) then

      if (master) then
        write(*,*) 'The number of MPI processes does not correspond to'        &
                                              // ' the number of domain chunks!'
        write(*,*) 'Looking for the best division...'
      end if

! try to find the best division
!
      pdims(:) = 1
      iret      = 0

      do while(product(pdims(:)) .lt. nprocs)
#ifdef R3D
        iret = mod(iret, 3) + 1
#else /* R3D */
        iret = mod(iret, 2) + 1
#endif /* R3D */
        pdims(iret) = 2 * pdims(iret)
      end do

! check if the best division found
!
      if (product(pdims(:)) .ne. nprocs) then

        if (master) then
          write(*,*) 'Improssible to find the best domain division! Exiting...'
          write(*,*)
        end if

        call finalize_mpi()
        stop

      end if

      if (master) then
        write(*,*) 'Found the best division:', pdims(:)
        write(*,*)
      end if

    else

! substitute div(:) to pdims(:)
!
      pdims(:) = div(:)

    end if

! set the periodic flag
!
    periodic(:) = per(:)

! set up the Cartesian geometry
!
    call mpi_cart_create(mpi_comm_world, 3, pdims(:), periodic(:)              &
                                                       , .true., comm3d, iret)

    if (iret .ne. mpi_success) then

      if (master) then
        write(*,*) 'The MPI could not create the Cartesian geometry! Exiting...'
        write(*,*)
      end if
      stop

    end if

! assign process coordinate
!
    call mpi_cart_coords(comm3d, nproc, 3, pcoords(:), iret)

    if (iret .ne. mpi_success) then

      if (master) then
        write(*,*) 'The MPI could not assign process coordinates! Exiting...'
        write(*,*)
      end if
      stop

    end if

! set the neighbors
!
    if (pdims(1) .gt. 1) then
      call mpi_cart_shift(comm3d, 0, 1, pneighs(1,1), pneighs(1,2), iret)
    end if
    if (pdims(2) .gt. 1) then
      call mpi_cart_shift(comm3d, 1, 1, pneighs(2,1), pneighs(2,2), iret)
    end if
    if (pdims(3) .gt. 1) then
      call mpi_cart_shift(comm3d, 2, 1, pneighs(3,1), pneighs(3,2), iret)
    end if

! set parity flag
!
    pparity(1) = mod(pcoords(1), 2)
    pparity(2) = mod(pcoords(2), 2)
    pparity(3) = mod(pcoords(3), 2)

! stop time accounting for the MPI initialization
!
    call stop_timer(imi)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine setup_mpi
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
    use mpi, only : mpi_integer, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: ibuf
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_bcast(ibuf, 1, mpi_integer, 0, comm3d, iret)

    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not broadcast an integer variable!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real   , intent(inout) :: rbuf
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_bcast(rbuf, 1, mpi_real8, 0, comm3d, iret)

    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not broadcast an integer variable!'
      write(*,*)
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
    use mpi, only : mpi_character, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(inout) :: sbuf
    integer         , intent(out)   :: iret
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_bcast(sbuf, len(sbuf), mpi_character, 0, comm3d, iret)

    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not broadcast a string variable!'
      write(*,*)
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
    use mpi, only : mpi_integer, mpi_min, mpi_success

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
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_min, comm3d, iret)

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the minimum value!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_min, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real   , intent(inout) :: rbuf
    integer, intent(out)   :: iret

! local variables
!
    real                   :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_min, comm3d, iret)

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the minimum value!'
      write(*,*)
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
    use mpi, only : mpi_integer, mpi_max, mpi_success

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
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_max, comm3d, iret)

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maximum value!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_max, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real   , intent(inout) :: rbuf
    integer, intent(out)   :: iret

! local variables
!
    real                   :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_max, comm3d, iret)

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maximum value!'
      write(*,*)
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
    use mpi, only : mpi_integer, mpi_sum, mpi_success

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
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(ibuf, tbuf, 1, mpi_integer, mpi_sum, comm3d, iret)

! substitute the result
!
    ibuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maximum value!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real   , intent(inout) :: rbuf
    integer, intent(out)   :: iret

! local variables
!
    real                   :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, 1, mpi_real8, mpi_sum, comm3d, iret)

! substitute the result
!
    rbuf = tbuf

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not sum the values from all processes!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_min, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: rbuf
    integer              , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n)           :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_min, comm3d, iret)

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the minima for all array elements!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_max, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: rbuf
    integer              , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n)           :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_max, comm3d, iret)

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maxima for all array elements!'
      write(*,*)
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
    use mpi, only : mpi_integer, mpi_sum, mpi_success

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
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(ibuf, tbuf, n, mpi_integer, mpi_sum, comm3d, iret)

! substitute the result
!
    ibuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maxima for all array elements!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_sum, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)    :: n
    real   , dimension(n), intent(inout) :: rbuf
    integer              , intent(out)   :: iret

! local variables
!
    real(kind=8), dimension(n)           :: tbuf
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_allreduce(rbuf, tbuf, n, mpi_real8, mpi_sum, comm3d, iret)

! substitute the result
!
    rbuf(:) = tbuf(:)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maxima for all array elements!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_sum, mpi_success

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
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    tbuf(:) = real(cbuf(:))
    call mpi_allreduce(tbuf, rbuf, n, mpi_real8, mpi_sum, comm3d, iret)
    tbuf(:) = aimag(cbuf(:))
    call mpi_allreduce(tbuf, ibuf, n, mpi_real8, mpi_sum, comm3d, iret)

! substitute the result
!
    cbuf(1:n) = cmplx(rbuf(1:n), ibuf(1:n))

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not find the maxima for all array elements!'
      write(*,*)
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
    use mpi, only : mpi_real8, mpi_success

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n, dst, tag
    real   , dimension(n), intent(in)  :: rbuf
    integer              , intent(out) :: iret
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_send(rbuf, n, mpi_real8, dst, tag, comm3d, iret)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not send the real array to another process!'
      write(*,*)
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
!     src  - the ID of the source process;
!     tag  - the tag identifying this operation;
!     rbuf - the received real array;
!     iret - the result flag identifying if the operation was successful;
!
!===============================================================================
!
  subroutine receive_real_array(n, src, tag, rbuf, iret)

! include external procedures and variables
!
    use mpi, only : mpi_real8, mpi_success, mpi_status_size

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n, src, tag
    real   , dimension(n), intent(out) :: rbuf
    integer              , intent(out) :: iret

! local variables
!
    integer :: status(mpi_status_size)
!
!-------------------------------------------------------------------------------
!
! start time accounting for the MPI communication
!
    call start_timer(imc)

    call mpi_recv(rbuf, n, mpi_real8, src, tag, comm3d, status, iret)

! check if the operation was successful
!
    if (iret .ne. mpi_success .and. master) then
      write(*,*) 'The MPI could not send the real array to another process!'
      write(*,*)
    end if

! stop time accounting for the MPI communication
!
    call stop_timer(imc)

!-------------------------------------------------------------------------------
!
  end subroutine receive_real_array
#endif /* MPI */

!===============================================================================
!
end module mpitools
