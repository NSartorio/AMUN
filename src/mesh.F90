!!*****************************************************************************
!!
!! module: mesh - handling adaptive mesh structure
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
module mesh

  implicit none

  contains
!
!======================================================================
!
! init_mesh: subroutine initializes mesh by creating blocks according
!            to the geometry, initial problem and refinement criterium
!
!======================================================================
!
  subroutine init_mesh

    use config , only : iblocks, jblocks, kblocks                     &
                      , xmin, xmax, ymin, ymax, zmin, zmax
    use blocks , only : list_allocated, init_blocks, clear_blocks     &
                      , allocate_blocks, block
    use error  , only : print_info
    use problem, only : init_problem

    implicit none

! local variables
!
    type(block), pointer :: pgroup, pblock

!----------------------------------------------------------------------
!
! check if the list is allocated, if yes deallocate it
!
    if (list_allocated()) then
      call print_info("mesh::init_mesh", "Block list is allocated, deallocate it!")

      call clear_blocks
    endif

! allocate initial structure of blocks according the the defined geometry
!
! TODO: by default we initiate 2x2=4 blocks in N configuration
! TODO: in the future allow user to define an arbitrary shape
!
    call init_blocks
    call allocate_blocks('N', pgroup, xmin, xmax, ymin, ymax, zmin, zmax)

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out these blocks with the initial condition
!
    pblock => pgroup
    do while (associated(pblock))

! set level
!
      pblock%level  = 1
      pblock%refine = .false.

! set initial conditions
!
      call init_problem(pblock)

! TODO: for all allocated blocks check criterium
!

! assign pointer to the current chunk
!
      pblock => pblock%next

    end do

! at this point the inital blocks are allocated and set for refinement,
! so iterate over all levels from 1 to maxlevel and create sub-blocks,
! set the initial conditions for each, check criterium and set for
! refinement according to the criterium fullfilment
!
! TODO: iterate over all levels, set neighbors of refined blocks for
!       refinement too, to keep the level jump not larger then 2,
!       refine blocks, set inital conditions at newly created block,
!       and finally check the criterium

! at this point the initial structure of blocks should be ready, and
! the problem should be set and refined

!----------------------------------------------------------------------
!
  end subroutine init_mesh
!
!======================================================================
!
! clears_mesh: subroutine deallocates mesh, removing blocks
!
!======================================================================
!
  subroutine clear_mesh

    use blocks, only : clear_blocks
    use error , only : print_info

    implicit none

!----------------------------------------------------------------------
!
! deallocate block structure
!
    call clear_blocks

!----------------------------------------------------------------------
!
  end subroutine clear_mesh

!======================================================================
!
end module
