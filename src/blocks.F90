!!*****************************************************************************
!!
!! module: blocks - handling block storage
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

    logical              :: refine

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

! stored last id (should always increase)
!
  integer(kind=4)     , save :: last_id

  contains
!
!======================================================================
!
! list_allocated: function checks if the block list is empty and
!                 returns true or false
!
!======================================================================
!
  function list_allocated

    implicit none

! output arguments
!
    logical :: list_allocated
!
!----------------------------------------------------------------------
!
    list_allocated = associated(pfirst)

    return

!----------------------------------------------------------------------
!
  end function list_allocated
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

! set unique ID
!
    pblock%id = increase_id()

! reset neighbors
!
    pblock%neigh(:,:) = -1

! allocate space for variables
!
#ifdef R3D
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
#else /* R3D */
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
#endif /* R3D */

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
! append_block: subroutine allocates space for one block and appends
!               it to the list
!
!======================================================================
!
  subroutine append_block(pblock)

    implicit none

! output arguments
!
    type(block), pointer, intent(out) :: pblock
!
!----------------------------------------------------------------------
!
! allocate block
!
    call allocate_block(pblock)

! add to the list
!
    if (list_allocated()) then
      pblock%prev => plast
      nullify(pblock%next)
      plast%next => pblock

      plast => pblock
    else
      pfirst => pblock
      plast => pblock
      nullify(pblock%prev)
      nullify(pblock%next)
    endif

!----------------------------------------------------------------------
!
  end subroutine append_block
!
!======================================================================
!
! allocate_blocks: subroutine allocates a configuration of blocks
!
!======================================================================
!
  subroutine allocate_blocks(block_config, pgroup, xmn, xmx, ymn, ymx &
                           , zmn, zmx)

    use error, only : print_error

    implicit none

! input parameters
!
    character(len=1), intent(in) :: block_config
    real            , intent(in) :: xmn, xmx, ymn, ymx, zmn, zmx

! output arguments
!
    type(block), pointer, intent(out) :: pgroup

! local pointers
!
    type(block), pointer :: pbl, pbr, ptl, ptr

! local variables
!
    real :: xl, xc, xr, yl, yc, yr, zl, zc, zr
!
!----------------------------------------------------------------------
!
    select case(block_config)
    case('n', 'N')

! TODO: create 4 blocks in N configuration; set pointers, neighbors, etc.
!       return pointer to the allocated chain

! create bottom left block
!
      call append_block(pbl)
      call append_block(ptl)
      call append_block(ptr)
      call append_block(pbr)

! set neighbors
!
      pbl%neigh(1,2) = pbr%id
      pbl%neigh(2,2) = ptl%id

      pbr%neigh(1,1) = pbl%id
      pbr%neigh(2,2) = ptr%id

      ptl%neigh(1,2) = ptr%id
      ptl%neigh(2,1) = pbl%id

      ptr%neigh(1,1) = ptl%id
      ptr%neigh(2,1) = pbr%id

! set block bounds
!
      xl = xmn
      xc = xmn + (xmx - xmn) / 2
      xr = xmx
      yl = ymn
      yc = ymn + (ymx - ymn) / 2
      yr = ymx

      pbl%xmin = xl
      pbl%xmax = xc
      pbl%ymin = yl
      pbl%ymax = yc

      ptl%xmin = xl
      ptl%xmax = xc
      ptl%ymin = yc
      ptl%ymax = yr

      ptr%xmin = xc
      ptr%xmax = xr
      ptr%ymin = yc
      ptr%ymax = yr

      pbr%xmin = xc
      pbr%xmax = xr
      pbr%ymin = yl
      pbr%ymax = yc

    case default
      call print_error("blocks::allocate_blocks","Configuration '" // block_config // "' not supported! Terminating!")
    end select

! copy pointer of the first block in chain
!
    pgroup => pbl

!----------------------------------------------------------------------
!
  end subroutine allocate_blocks
!
!======================================================================
!
! init_blocks: subroutine initializes the structure of blocks
!
!======================================================================
!
  subroutine init_blocks()

    use config, only : ngrids, iblocks, jblocks, kblocks, nghost, ncells
!     use problem, only : init_problem

    implicit none

! local variables
!
    integer(kind=4) :: i, j, k

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

! reset ID
!
    last_id = 0

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

      print *, pfirst%id, pfirst%xmin, pfirst%xmax, pfirst%ymin, pfirst%ymax

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
!
!======================================================================
!
! increase_id: function increases last ID by 1 and returns it
!
!======================================================================
!
  function increase_id

    implicit none

! return variable
!
    integer(kind=4) :: increase_id
!
!----------------------------------------------------------------------
!
! increase ID by 1
!
    last_id = last_id + 1

! return ID
!
    increase_id = last_id

    return

!----------------------------------------------------------------------
!
  end function increase_id

!======================================================================
!
end module
