!!*****************************************************************************
!!
!! module: blocks - handling adaptive mesh structure
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
module blocks

  implicit none

! parameters
!
#ifdef R3D
  integer(kind=4), parameter :: ndims  = 3
#else /* R3D */
  integer(kind=4), parameter :: ndims  = 2
#endif /* R3D */
  integer(kind=4), parameter :: nchild = 2**ndims

! define block type
!
  type block
    type(block), pointer :: next, prev

    integer(kind=4)      :: id, level, parent
    integer(kind=4)      :: neigh(ndims,2), child(nchild)

    real                 :: xmin, xmax, ymin, ymax, zmin, zmax

    real, dimension(:,:,:), allocatable :: dn, mx, my, mz
#ifndef ISO
    real, dimension(:,:,:), allocatable :: en
#endif /* !ISO */
#ifdef MHD
    real, dimension(:,:,:), allocatable :: bx, by, bz
#endif /* MHD */
  end type block

! stored pointers
!
  type(block), pointer, save :: pfirst, plast
  integer(kind=4)     , save :: nblocks

  contains
!
!======================================================================
!
! allocate_block: subroutine allocates space for one block and returns
!                 pointer to this block
!
!======================================================================
!
  subroutine allocate_block(pblock)

    use config, only : ngrids

    implicit none

! output arguments
!
    type(block), pointer, intent(out) :: pblock
!
!----------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pblock)

! allocate space for variables
!
    if (ndims .eq. 2) then
      allocate(pblock%dn(ngrids,ngrids,1))
      allocate(pblock%mx(ngrids,ngrids,1))
      allocate(pblock%my(ngrids,ngrids,1))
      allocate(pblock%mz(ngrids,ngrids,1))
#ifndef ISO
      allocate(pblock%en(ngrids,ngrids,1))
#endif /* ISO */
#ifdef MHD
      allocate(pblock%bx(ngrids,ngrids,1))
      allocate(pblock%by(ngrids,ngrids,1))
      allocate(pblock%bz(ngrids,ngrids,1))
#endif /* MHD */
    endif
    if (ndims .eq. 3) then
      allocate(pblock%dn(ngrids,ngrids,ngrids))
      allocate(pblock%mx(ngrids,ngrids,ngrids))
      allocate(pblock%my(ngrids,ngrids,ngrids))
      allocate(pblock%mz(ngrids,ngrids,ngrids))
#ifndef ISO
      allocate(pblock%en(ngrids,ngrids,ngrids))
#endif /* !ISO */
#ifdef MHD
      allocate(pblock%bx(ngrids,ngrids,ngrids))
      allocate(pblock%by(ngrids,ngrids,ngrids))
      allocate(pblock%bz(ngrids,ngrids,ngrids))
#endif /* MHD */
    endif

!----------------------------------------------------------------------
!
  end subroutine allocate_block
!
!======================================================================
!
! deallocate_block: subroutine deallocates space ocuppied by a given
!                   block
!
!======================================================================
!
  subroutine deallocate_block(pblock)

    implicit none

! input arguments
!
    type(block), pointer, intent(inout) :: pblock
!
!----------------------------------------------------------------------
!
! deallocate variables
!
    deallocate(pblock%dn)
    deallocate(pblock%mx)
    deallocate(pblock%my)
    deallocate(pblock%mz)
#ifndef ISO
    deallocate(pblock%en)
#endif /* !ISO */
#ifdef MHD
    deallocate(pblock%bx)
    deallocate(pblock%by)
    deallocate(pblock%bz)
#endif /* MHD */

! free and nullify the block
!
    deallocate(pblock)
    nullify(pblock)

!----------------------------------------------------------------------
!
  end subroutine deallocate_block
!
!======================================================================
!
! init_blocks: subroutine initializes the structure of blocks
!
!======================================================================
!
  subroutine init_blocks

    use config, only : iblocks, jblocks, kblocks

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

! reset number of blocks
!
    nblocks = 0

! create the first block
!
    call allocate_block(pcurr)

! fill block structure
!
    pcurr%id         =  1
    pcurr%level      =  1
    pcurr%parent     = -1

    pcurr%neigh(:,:) = -1
    pcurr%child(:)   = -1

! nullify the prev and next fields
!
    nullify(pcurr%prev)
    nullify(pcurr%next)

! add the block to the list
!
    pfirst => pcurr

!----------------------------------------------------------------------
!
  end subroutine init_blocks
!
!======================================================================
!
! clear_blocks: subroutine clears the structure of blocks
!
!======================================================================
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
      call deallocate_block(pfirst)

! assign pointer to the current chunk
!
      pfirst => pcurr

    end do

!----------------------------------------------------------------------
!
  end subroutine clear_blocks

!======================================================================
!
end module
