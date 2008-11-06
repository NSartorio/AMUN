!!**********************************************************************
!!
!! module: blocks - handling adaptive mesh structure
!!
!! Copyright (C) 2008 Grzegorz Kowal <kowal@astro.wisc.edu>
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
module blocks

  implicit none

! parameters
!
  integer(kind=4), parameter :: nghost =  2
  integer(kind=4), parameter :: ncells = 12
  integer(kind=4), parameter :: ngrids = ncells + 2*nghost
#ifdef R3D
  integer(kind=4), parameter :: ndims  = 3
  integer(kind=4), parameter :: ngridx = ngrids
  integer(kind=4), parameter :: ngridy = ngrids
  integer(kind=4), parameter :: ngridz = ngrids
#else /* R3D */
  integer(kind=4), parameter :: ndims  = 2
  integer(kind=4), parameter :: ngridx = ngrids
  integer(kind=4), parameter :: ngridy = ngrids
  integer(kind=4), parameter :: ngridz = 1
#endif /* R3D */
  integer(kind=4), parameter :: nchild = 2**ndims

! define block type
!
  type block
    type(block), pointer :: next, prev

    integer(kind=4)      :: id, level, parent
    integer(kind=4)      :: neigh(ndims,2), child(nchild)

    real                 :: xmin, xmax, ymin, ymax, zmin, zmax

    real                 :: dn(ngridx,ngridy,ngridz)
    real                 :: mx(ngridx,ngridy,ngridz)
    real                 :: my(ngridx,ngridy,ngridz)
    real                 :: mz(ngridx,ngridy,ngridz)
#ifndef ISO
    real                 :: en(ngridx,ngridy,ngridz)
#endif /* !ISO */
#ifdef MHD
    real                 :: bx(ngridx,ngridy,ngridz)
    real                 :: by(ngridx,ngridy,ngridz)
    real                 :: bz(ngridx,ngridy,ngridz)
#endif /* MHD */
  end type block

! stored pointers
!
  type(block), pointer, save :: pfirst, plast

  contains
!
!===============================================================================
!
! init_blocks: subroutine initializes the structure of blocks
!
!===============================================================================
!
  subroutine init_blocks

    implicit none

! pointers
!
    type(block), pointer :: pcurr
!
!----------------------------------------------------------------------
!
! first check if block list is empty
!
    if (associated(pfirst)) then
      write(*,*) 'ERROR: Block list already associated!'
    endif

! nullify all pointers
!
    nullify(pfirst)

! create the first block
!
    allocate(pcurr)

! fill block structure
!
    pcurr%id         =  1
    pcurr%level      =  1
    pcurr%parent     = -1

    pcurr%neigh(:,:) = -1
    pcurr%child(:,:) = -1

! nullify the prev and next fields
!
    nullify(pcurr%prev)
    nullify(pcurr%next)

! add the block to the list
!
    pfirst => pcurr

  end subroutine init_blocks
!
!===============================================================================
!
! clear_blocks: subroutine clears the structure of blocks
!
!===============================================================================
!
  subroutine clear_blocks

    implicit none

! pointers
!
    type(block), pointer :: pcurr
!
!----------------------------------------------------------------------
!
! until list is free, reiterate over all chunks
!
    do while(associated(pfirst))

! assign temporary pointer to the next chunk
!
      pcurr => pfirst%next

! deallocate the content of current block
!

! deallocate and nullify the current block
!
      deallocate(pfirst)
      nullify(pfirst)

! assign pointer to the current chunk
!
      pfirst => pcurr

    end do

  end subroutine clear_blocks

end module
