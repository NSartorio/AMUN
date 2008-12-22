!!*****************************************************************************
!!
!! module: mpitools - subroutines for MPI communication
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
module mpitools

  implicit none

! MPI global variables
!
  integer(kind=4), save                 :: comm3d
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
!----------------------------------------------------------------------
!
    is_master = ncpu .eq. 0

!-------------------------------------------------------------------------------
!
  end function is_master

end module
