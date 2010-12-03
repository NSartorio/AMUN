!!******************************************************************************
!!
!! module: mpitools - subroutines for MPI communication
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module mpitools

  implicit none

! MPI global variables
!
  integer        , save                 :: comm3d
  integer(kind=4), save                 :: ncpu, ncpus

  contains
!
!===============================================================================
!
! init_mpi: subroutine initializes the MPI variables
!
!===============================================================================
!
  subroutine init_mpi

#ifdef MPI
    use mpi, only : mpi_comm_world
#endif /* MPI */

    implicit none
#ifdef MPI
! local variables
!
    integer :: err
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
    ncpu  = 0
    ncpus = 1

#ifdef MPI
!  initialize the MPI interface
!
    call mpi_init(err)

! get the current process id and the total number of processes
!
    call mpi_comm_rank(mpi_comm_world, ncpu , err)
    call mpi_comm_size(mpi_comm_world, ncpus, err)

    comm3d = mpi_comm_world

#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine init_mpi
!
!===============================================================================
!
! clear_mpi: subroutine clears the MPI variables
!
!===============================================================================
!
  subroutine clear_mpi

    implicit none
#ifdef MPI
! local variables
!
    integer :: err
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
!  finalize the MPI interface
!
    call mpi_finalize(err)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine clear_mpi
!
!===============================================================================
!
! mbarrier: subroutine synchronizes processes
!
!===============================================================================
!
  subroutine mbarrier

    implicit none
#ifdef MPI
! local variables
!
    integer :: err
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
!  finalize the MPI interface
!
    call mpi_barrier(comm3d, err)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mbarrier
!
!===============================================================================
!
! is_master: function returns true if it is the master node, otherwise it
!            returns false
!
!===============================================================================
!
  function is_master()

    implicit none

! return value
!
    logical :: is_master
!
!-------------------------------------------------------------------------------
!
    is_master = ncpu .eq. 0

!-------------------------------------------------------------------------------
!
  end function is_master
!
!===============================================================================
!
! msendi: subroutine sends an array
!
!===============================================================================
!
  subroutine msendi(n, dst, tag, buf)

#ifdef MPI
    use mpi, only : mpi_integer
#endif /* MPI */

    implicit none

! arguments
!
    integer              , intent(in) :: n, dst, tag
    integer, dimension(n), intent(in) :: buf

#ifdef MPI
! local variables
!
    integer :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_send(buf, n, mpi_integer, dst, tag, comm3d, err)
    if (err .ne. 0) print *, 'msendi: error', err
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine msendi
!
!===============================================================================
!
! mrecvi: subroutine receives an array
!
!===============================================================================
!
  subroutine mrecvi(n, src, tag, buf)

#ifdef MPI
    use mpi, only : mpi_status_size, mpi_integer
#endif /* MPI */

! arguments
!
    integer              , intent(in)  :: n, src, tag
    integer, dimension(n), intent(out) :: buf

#ifdef MPI
! local variables
!
    integer :: err, status(mpi_status_size)
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
    buf(:)    = 0
#ifdef MPI
    err       = 0
    status(:) = 0
    call mpi_recv(buf, n, mpi_integer, src, tag, comm3d, status, err)
    if (err .ne. 0) print *, 'mrecvi: error', err
#endif /* MPI */

  end subroutine mrecvi
!
!===============================================================================
!
! msendf: subroutine sends an array
!
!===============================================================================
!
  subroutine msendf(n, dst, tag, buf)

#ifdef MPI
    use mpi, only : mpi_real8
#endif /* MPI */

    implicit none

! arguments
!
    integer                   , intent(in)    :: n, dst, tag
    real(kind=8), dimension(n), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer :: err
!
!-------------------------------------------------------------------------------
!
    call mpi_send(buf, n, mpi_real8, dst, tag, comm3d, err)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine msendf
!
!===============================================================================
!
! mrecvf: subroutine receives an array
!
!===============================================================================
!
  subroutine mrecvf(n, src, tag, buf)

#ifdef MPI
    use mpi, only : mpi_status_size, mpi_real8
#endif /* MPI */

! arguments
!
    integer                   , intent(in)    :: n, src, tag
    real(kind=8), dimension(n), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer :: err, status(mpi_status_size)
!
!-------------------------------------------------------------------------------
!
    call mpi_recv(buf, n, mpi_real8, src, tag, comm3d, status, err)
#endif /* MPI */

  end subroutine mrecvf
!
!===============================================================================
!
! mallreducesuml: subroutine adds values over all proceeses
!
!===============================================================================
!
  subroutine mallreducesuml(n, buf)

#ifdef MPI
    use mpi, only : mpi_integer, mpi_sum
#endif /* MPI */

! arguments
!
    integer              , intent(in)    :: n
    integer, dimension(n), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer, dimension(n) :: tbuf
    integer               :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, n, mpi_integer, mpi_sum, comm3d, err)
    buf(1:n) = tbuf(1:n)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mallreducesuml
!
!===============================================================================
!
! mallreducesuml: subroutine adds values over all proceeses
!
!===============================================================================
!
  subroutine mallreduceprodl(n, buf)

#ifdef MPI
    use mpi, only : mpi_integer, mpi_prod
#endif /* MPI */

! arguments
!
    integer              , intent(in)    :: n
    integer, dimension(n), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer, dimension(n) :: tbuf
    integer               :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, n, mpi_integer, mpi_prod, comm3d, err)
    buf(1:n) = tbuf(1:n)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mallreduceprodl
!
!===============================================================================
!
! mallreducemaxl: subroutine finds maximum values over all proceeses
!
!===============================================================================
!
  subroutine mallreducemaxl(n, buf)

#ifdef MPI
    use mpi, only : mpi_integer, mpi_max
#endif /* MPI */

! arguments
!
    integer              , intent(in)    :: n
    integer, dimension(n), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer, dimension(n) :: tbuf
    integer               :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, n, mpi_integer, mpi_max, comm3d, err)
    buf(1:n) = tbuf(1:n)
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mallreducemaxl
!
!===============================================================================
!
! mallreduceminr: subroutine finds the minimum value over all proceeses
!
!===============================================================================
!
  subroutine mallreduceminr(buf)

#ifdef MPI
    use mpi, only : mpi_real8, mpi_min
#endif /* MPI */

! arguments
!
    real(kind=8), intent(inout) :: buf

#ifdef MPI
! local variables
!
    real(kind=8)        :: tbuf
    integer             :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, 1, mpi_real8, mpi_min, comm3d, err)
    buf = tbuf
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mallreduceminr
!
!===============================================================================
!
! mallreducemaxr: subroutine reduces the maximum value over all processes
!
!===============================================================================
!
  subroutine mallreducemaxr(buf)

#ifdef MPI
    use mpi, only : mpi_real8, mpi_max
#endif /* MPI */

! arguments
!
    real(kind=8), intent(inout) :: buf

#ifdef MPI
! local variables
!
    real(kind=8)        :: tbuf
    integer             :: err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, 1, mpi_real8, mpi_max, comm3d, err)
    buf = tbuf
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine mallreducemaxr
!
!===============================================================================
!
! mfindmaxi: subroutine finds the maximum integer value across all proceeses
!
!===============================================================================
!
  subroutine mfindmaxi(buf)

#ifdef MPI
    use mpi, only : mpi_integer, mpi_max
#endif /* MPI */

! arguments
!
    integer(kind=4), intent(inout) :: buf

#ifdef MPI
! local variables
!
    integer(kind=4)     :: tbuf, err
!
!-------------------------------------------------------------------------------
!
    err = 0
    call mpi_allreduce(buf, tbuf, 1, mpi_integer, mpi_max, comm3d, err)
    buf = tbuf
#endif /* MPI */
!-------------------------------------------------------------------------------
!
  end subroutine mfindmaxi

end module
