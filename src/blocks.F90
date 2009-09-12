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
  integer(kind=4), parameter :: ndims  = NDIMS
  integer(kind=4), parameter :: nsides = 2
  integer(kind=4), parameter :: nfaces = 2**(ndims-1)
  integer(kind=4), parameter :: nchild = 2**ndims
  integer(kind=4), parameter :: idn = 1, imx = 2, imy = 3, imz = 4 &
                              , ivx = 2, ivy = 3, ivz = 4
#ifdef HYDRO
#ifdef ISO
  integer(kind=4), parameter :: nvars  = 4
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien = 5, ipr = 5
  integer(kind=4), parameter :: nvars  = 5
#endif /* ADI */
#endif /* HYDRO */
#ifdef MHD
#ifdef ISO
  integer(kind=4), parameter :: ibx = 5, iby = 6, ibz = 7
  integer(kind=4), parameter :: nvars  = 7
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien = 5, ibx = 6, iby = 7, ibz = 8
  integer(kind=4), parameter :: nvars  = 8
#endif /* ADI */
#endif /* MHD */
  integer(kind=4), parameter :: maxid = 1000000

! define block type
!
  type blockptr
    type(block), pointer :: ptr
  end type blockptr

  type blockref
    integer(kind=4)      :: cpu, id
  end type blockref

! define pointers to block_meta and block_data structures
!
  type pointer_meta
    type(block_meta), pointer :: ptr
  end type pointer_meta

  type pointer_data
    type(block_data), pointer :: ptr
  end type pointer_data

! define block_meta structure
!
  type block_meta
    type(block_meta)  , pointer :: prev             ! pointer to the previous block
    type(block_meta)  , pointer :: next             ! pointer to the next block
    type(block_meta)  , pointer :: parent           ! pointer to the parent block
    type(pointer_meta)          :: child(nchild)    ! pointers to children
    type(pointer_meta)          :: neigh(ndims,2,2) ! pointers to neighbors

    type(block_data)  , pointer :: data             ! pointer to the data block

    integer(kind=4)             :: id               ! block identificator
    integer(kind=4)             :: cpu              ! the cpu id of the block
    integer(kind=4)             :: level            ! refinement level
    integer(kind=4)             :: config           ! configuration flag
    integer(kind=4)             :: refine           ! refinement flag:
!                                                       -1 - derefine
!                                                        0 - do nothing
!                                                        1 - refine

    logical                     :: leaf             ! leaf flag
  end type block_meta

! define block_data structure
!
  type block_data
    type(block_data), pointer :: prev             ! pointer to the previous block
    type(block_data), pointer :: next             ! pointer to the next block

    type(block_meta), pointer :: meta             ! pointer to the metadata block

    real                      :: xmin, xmax       ! bounds for the x direction
    real                      :: ymin, ymax       ! bounds for the y direction
    real                      :: zmin, zmax       ! bounds for the z direction

    real, dimension(:,:,:,:), allocatable :: u    ! variable array
    real, dimension(:,:,:)  , allocatable :: c    ! criterion array
  end type block_data


  type block
    type(block), pointer :: next, prev
    type(blockref)       :: parent, child(nchild), neigh(ndims,2,2)

    logical              :: leaf
    integer(kind=4)      :: cpu, id, level

    character            :: config
    integer(kind=4)      :: refine
    integer(kind=4)      :: pos(ndims)


    real                 :: xmin, xmax, ymin, ymax, zmin, zmax

    real, dimension(:,:,:,:), allocatable :: u
    real, dimension(:,:,:)  , allocatable :: c
  end type block

! array of ID to pointer conversion
!
  type(blockptr), dimension(:), allocatable, save :: idtoptr

! chains of meta blocks and data blocks
!
  type(block_meta), pointer, save :: list_meta, last_meta
  type(block_data), pointer, save :: list_data, last_data

! stored pointers
!
  type(block), pointer, save :: plist, plast

! stored last id (should always increase)
!
  integer(kind=4)     , save :: last_id

! store number of allocated blocks and leafs
!
  integer(kind=4)     , save :: nblocks, dblocks, nleafs

  contains
!
!===============================================================================
!
! init_blocks: subroutine initializes the block variables
!
!===============================================================================
!
  subroutine init_blocks()

    use error, only : print_warning

    implicit none

! local variables
!
    integer(kind=4) :: p
!
!-------------------------------------------------------------------------------
!
! check if metadata list is empty
!
    if (associated(list_meta)) &
      call print_warning("blocks::init_blocks", "Block metadata list is already associated!")

! check if data list is empty
!
    if (associated(list_data)) &
      call print_warning("blocks::init_blocks", "Block data list is already associated!")

! nullify all pointers
!
    nullify(list_meta)
    nullify(list_data)

! reset number of blocks and leafs
!
    nblocks = 0
    dblocks = 0
    nleafs  = 0

! reset ID
!
    last_id = 0

!-------------------------------------------------------------------------------
!
  end subroutine init_blocks
!
!===============================================================================
!
! clear_blocks: subroutine clears the block variables
!
!===============================================================================
!
  subroutine clear_blocks

    implicit none

! pointers
!
    type(block_meta), pointer :: pblock_meta
    type(block_data), pointer :: pblock_data
!
!-------------------------------------------------------------------------------
!
! clear all data blocks
!
    pblock_data => list_data
    do while(associated(pblock_data))
      call deallocate_datablock(pblock_data)
      pblock_data => list_data
    end do

! clear all meta blocks
!
    pblock_meta => list_meta
    do while(associated(pblock_meta))
      call deallocate_metablock(pblock_meta)
      pblock_meta => list_meta
    end do

!-------------------------------------------------------------------------------
!
  end subroutine clear_blocks
!
!===============================================================================
!
! list_allocated: function returns true if the block list is allocated,
!                 otherwise it returns false
!
!===============================================================================
!
  function list_allocated

    implicit none

! output arguments
!
    logical :: list_allocated
!
!-------------------------------------------------------------------------------
!
    list_allocated = associated(plist)

    return

!-------------------------------------------------------------------------------
!
  end function list_allocated
!
!===============================================================================
!
! increase_id: function increases the last ID by 1 and returns it
!
!===============================================================================
!
  function increase_id

    implicit none

! return variable
!
    integer(kind=4) :: increase_id
!
!-------------------------------------------------------------------------------
!
! increase ID by 1
!
    last_id = last_id + 1

! return ID
!
    increase_id = last_id

    return

!-------------------------------------------------------------------------------
!
  end function increase_id
!
!===============================================================================
!
! get_pointer: function returns a pointer to the block with requested ID
!
!===============================================================================
!
  function get_pointer(id)

    implicit none

! input argument
!
    integer(kind=4) :: id

! return variable
!
    type(block), pointer :: get_pointer
!
!-------------------------------------------------------------------------------
!
    nullify(get_pointer)

    if (id .ge. 1 .and. id .le. maxid) &
      get_pointer => idtoptr(id)%ptr

    return

!-------------------------------------------------------------------------------
!
  end function get_pointer
!
!===============================================================================
!
! allocate_block: subroutine allocates space for one block and returns the
!                 pointer to this block
!
!===============================================================================
!
  subroutine allocate_block(pblock)

    use config  , only : im, jm, km
    use mpitools, only : ncpu

    implicit none

! output arguments
!
    type(block), pointer, intent(out) :: pblock

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pblock)

! set unique ID
!
    pblock%id = increase_id()   ! TODO: replace with get_free_id() which return the first free id

! set configuration and leaf flags
!
    pblock%config = 'N'         ! TODO: replace with an integer number
    pblock%leaf   = .false.

! set the cpu of current block
!
    pblock%cpu    = ncpu

! initialize the refinement flag
!
    pblock%refine = 0

! nullify pointers
!
    nullify(pblock%next)
    nullify(pblock%prev)

! reset parent block
!
    pblock%parent%cpu = -1
    pblock%parent%id  = -1

! reset neighbors
!
    pblock%neigh(:,:,:)%cpu = -1
    pblock%neigh(:,:,:)%id  = -1

! allocate space for variables
!
    allocate(pblock%u(nvars,im,jm,km))
    allocate(pblock%c(im,jm,km))

! set the correspponding pointer in the ID to pointer array to the current block
!
    idtoptr(pblock%id)%ptr => pblock

! increase the number of allocated blocks
!
    nblocks = nblocks + 1

!-------------------------------------------------------------------------------
!
  end subroutine allocate_block
!
!===============================================================================
!
! allocate_metablock: subroutine allocates space for one meta block and returns
!                     the pointer to this block
!
!===============================================================================
!
  subroutine allocate_metablock(pblock)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(out) :: pblock

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pblock)

! nullify pointers
!
    nullify(pblock%prev)
    nullify(pblock%next)
    nullify(pblock%parent)
    nullify(pblock%data)
    do i = 1, nchild
      nullify(pblock%child(i)%ptr)
    end do
    do k = 1, nfaces
      do j = 1, nsides
        do i = 1, ndims
          nullify(pblock%neigh(i,j,k)%ptr)
        end do
      end do
    end do

! set unique ID
!
    pblock%id = increase_id()

! unset the CPU number of current block, level, the configuration, refine and
! leaf flags
!
    pblock%cpu    = -1
    pblock%level  = -1
    pblock%config = -1
    pblock%refine =  0
    pblock%leaf   = .false.

! increase the number of allocated meta blocks
!
    nblocks = nblocks + 1

!-------------------------------------------------------------------------------
!
  end subroutine allocate_metablock
!
!===============================================================================
!
! allocate_datablock: subroutine allocates space for one data block and returns
!                     the pointer to this block
!
!===============================================================================
!
  subroutine allocate_datablock(pblock)

    use config, only : im, jm, km

    implicit none

! output arguments
!
    type(block_data), pointer, intent(out) :: pblock

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pblock)

! nullify pointers
!
    nullify(pblock%prev)
    nullify(pblock%next)
    nullify(pblock%meta)

! allocate space for variables
!
    allocate(pblock%u(nvars,im,jm,km))
    allocate(pblock%c(im,jm,km))

! initialize bounds of the block
!
    pblock%xmin = 0.0
    pblock%xmax = 1.0
    pblock%ymin = 0.0
    pblock%ymax = 1.0
    pblock%zmin = 0.0
    pblock%zmax = 1.0

! increase the number of allocated meta blocks
!
    dblocks = dblocks + 1

!-------------------------------------------------------------------------------
!
  end subroutine allocate_datablock
!
!===============================================================================
!
! deallocate_block: subroutine deallocates space ocuppied by a given block
!
!===============================================================================
!
  subroutine deallocate_block(pblock)

    implicit none

! input arguments
!
    type(block), pointer, intent(inout) :: pblock
!
!-------------------------------------------------------------------------------
!
    if (associated(pblock)) then

! if this is the first block in the list, update the plist pointer
!
      if (pblock%id .eq. plist%id) &
        plist => pblock%next

! update the pointer of previous and next blocks
!
      if (associated(pblock%prev)) &
        pblock%prev%next => pblock%next

      if (associated(pblock%next)) &
        pblock%next%prev => pblock%prev

! deallocate variables
!
      deallocate(pblock%u)
      deallocate(pblock%c)

! nullify pointers
!
      nullify(pblock%next)
      nullify(pblock%prev)

! nullify the corresponding pointer in the ID to pointer array
!
      nullify(idtoptr(pblock%id)%ptr)

! free and nullify the block
!
      deallocate(pblock)
      nullify(pblock)

! decrease the number of allocated blocks
!
      nblocks = nblocks - 1
    endif

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_block
!
!===============================================================================
!
! deallocate_metablock: subroutine deallocates space ocuppied by a given metablock
!
!===============================================================================
!
  subroutine deallocate_metablock(pblock)

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
    if (associated(pblock)) then

! if this is the first block in the list, update the plist pointer
!
      if (pblock%id .eq. list_meta%id) &
        list_meta => pblock%next

! update the pointer of previous and next blocks
!
      if (associated(pblock%prev)) &
        pblock%prev%next => pblock%next

      if (associated(pblock%next)) &
        pblock%next%prev => pblock%prev

! nullify children
!
      do i = 1, nchild
        nullify(pblock%child(i)%ptr)
      end do

! nullify neighbors
!
      do k = 1, nfaces
        do j = 1, nsides
          do i = 1, ndims
            nullify(pblock%neigh(i,j,k)%ptr)
          end do
        end do
      end do

! nullify pointers
!
      nullify(pblock%next)
      nullify(pblock%prev)
      nullify(pblock%data)
      nullify(pblock%parent)

! free and nullify the block
!
      deallocate(pblock)
      nullify(pblock)
    endif

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_metablock
!
!===============================================================================
!
! deallocate_datablock: subroutine deallocates space ocuppied by a given datablock
!
!===============================================================================
!
  subroutine deallocate_datablock(pblock)

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock
!
!-------------------------------------------------------------------------------
!
    if (associated(pblock)) then

! if this is the first block in the list, update the plist pointer
!
      if (pblock%meta%id .eq. list_data%meta%id) &
        list_data => pblock%next

! update the pointer of previous and next blocks
!
      if (associated(pblock%prev)) &
        pblock%prev%next => pblock%next

      if (associated(pblock%next)) &
        pblock%next%prev => pblock%prev

! deallocate variables
!
      deallocate(pblock%u)
      deallocate(pblock%c)

! nullify pointers
!
      nullify(pblock%next)
      nullify(pblock%prev)
      nullify(pblock%meta)

! free and nullify the block
!
      deallocate(pblock)
      nullify(pblock)
    endif

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_datablock
!
!===============================================================================
!
! append_block: subroutine allocates space for one block and appends it to
!               the list
!
!===============================================================================
!
  subroutine append_block(pblock)

    implicit none

! output arguments
!
    type(block), pointer, intent(out) :: pblock
!
!-------------------------------------------------------------------------------
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
      plist => pblock
      plast => pblock
      nullify(pblock%prev)
      nullify(pblock%next)
    endif

!-------------------------------------------------------------------------------
!
  end subroutine append_block
!
!===============================================================================
!
! append_metablock: subroutine allocates space for one meta block and appends it
!                   to the meta block list
!
!===============================================================================
!
  subroutine append_metablock(pblock)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(out) :: pblock
!
!-------------------------------------------------------------------------------
!
! allocate block
!
    call allocate_metablock(pblock)

! add to the list
!
    if (associated(list_meta)) then
      pblock%prev => last_meta
      last_meta%next => pblock
    else
      list_meta => pblock
    endif

! set the pointer to the last block in the list
!
    last_meta => pblock

!-------------------------------------------------------------------------------
!
  end subroutine append_metablock
!
!===============================================================================
!
! append_datablock: subroutine allocates space for one data block and appends it
!                   to the data block list
!
!===============================================================================
!
  subroutine append_datablock(pblock)

    implicit none

! output arguments
!
    type(block_data), pointer, intent(out) :: pblock
!
!-------------------------------------------------------------------------------
!
! allocate block
!
    call allocate_datablock(pblock)

! add to the list
!
    if (associated(list_data)) then
      pblock%prev => last_data
      last_data%next => pblock
    else
      list_data => pblock
    endif

! set the pointer to the last block in the list
!
    last_data => pblock

!-------------------------------------------------------------------------------
!
  end subroutine append_datablock
!
!===============================================================================
!
! associate_blocks: subroutine associates a pair of meta and data blocks
!
!===============================================================================
!
  subroutine associate_blocks(pblock_meta, pblock_data)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(inout) :: pblock_meta
    type(block_data), pointer, intent(inout) :: pblock_data
!
!-------------------------------------------------------------------------------
!
    pblock_meta%data => pblock_data
    pblock_data%meta => pblock_meta
!
!-------------------------------------------------------------------------------
!
  end subroutine associate_blocks
!
!===============================================================================
!
! metablock_setleaf: subroutine sets the leaf flag of data block
!
!===============================================================================
!
  subroutine metablock_setleaf(pblock)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pblock
!
!-------------------------------------------------------------------------------
!
! set the leaf flag
!
    pblock%leaf = .true.

! increase the number of leafs
!
    nleafs = nleafs + 1
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_setleaf
!
!===============================================================================
!
! metablock_unsetleaf: subroutine unsets the leaf flag of data block
!
!===============================================================================
!
  subroutine metablock_unsetleaf(pblock)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pblock
!
!-------------------------------------------------------------------------------
!
! set the leaf flag
!
    pblock%leaf = .false.

! decrease the number of leafs
!
    nleafs = nleafs - 1
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_unsetleaf
!
!===============================================================================
!
! metablock_setconfig: subroutine sets the config flag of data block
!
!===============================================================================
!
  subroutine metablock_setconfig(pblock, config)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pblock
    integer(kind=4)          , intent(in)    :: config
!
!-------------------------------------------------------------------------------
!
! set the config flag
!
    pblock%config = config
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_setconfig
!
!===============================================================================
!
! metablock_setlevel: subroutine sets the level of data block
!
!===============================================================================
!
  subroutine metablock_setlevel(pblock, level)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pblock
    integer(kind=4)          , intent(in)    :: level
!
!-------------------------------------------------------------------------------
!
! set the refinement level
!
    pblock%level = level
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_setlevel
!
!===============================================================================
!
! datablock_setbounds: subroutine sets the bounds of data block
!
!===============================================================================
!
  subroutine datablock_setbounds(pblock, xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

! input/output arguments
!
    type(block_data), pointer, intent(inout) :: pblock
    real                     , intent(in)    :: xmin, xmax, ymin, ymax, zmin, zmax
!
!-------------------------------------------------------------------------------
!
! set bounds of the block
!
    pblock%xmin = xmin
    pblock%xmax = xmax
    pblock%ymin = ymin
    pblock%ymax = ymax
    pblock%zmin = zmin
    pblock%zmax = zmax
!
!-------------------------------------------------------------------------------
!
  end subroutine datablock_setbounds
!
!===============================================================================
!
! refine_block: subroutine refines selected block
!
!===============================================================================
!
  subroutine refine_block(pblock, falloc_data)

    use error   , only : print_error
#ifdef MPI
    use mpitools, only : ncpu
#endif /* MPI */

    implicit none

! input parameters
!
    type(block_meta), pointer, intent(inout) :: pblock
    logical                  , intent(in)    :: falloc_data

! pointers
!
    type(block_meta), pointer :: pneigh, pchild, pfirst, plast
    type(block_data), pointer :: pdata

! local arrays
!
    integer, dimension(nchild) :: config, order

! local variables
!
    integer :: p, i, j, k
    real    :: xln, yln, zln, xmn, xmx, ymn, ymx, zmn, zmx
!
!-------------------------------------------------------------------------------
!
! check if pointer is associated
!
    if (associated(pblock)) then

! unset block leaf flag
!
      call metablock_unsetleaf(pblock)

! reset refinement flag
!
      pblock%refine = 0

! iterate over all child blocks
!
      do p = 1, nchild

! create child meta and data blocks
!
        call allocate_metablock(pblock%child(p)%ptr)

! set it as a leaf
!
        call metablock_setleaf(pblock%child(p)%ptr)

! assign pointer to the parent block
!
        pblock%child(p)%ptr%parent => pblock

! increase the refinement level
!
        pblock%child(p)%ptr%level = pblock%level + 1

! copy the parent cpu number to each child
!
        pblock%child(p)%ptr%cpu = pblock%cpu

      end do

! assign neighbors of the child blocks
!
      do p = 1, nfaces

! X direction (left side)
!
        pblock%child(1)%ptr%neigh(1,1,p)%ptr => pblock%neigh(1,1,1)%ptr
        pblock%child(2)%ptr%neigh(1,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(3)%ptr%neigh(1,1,p)%ptr => pblock%neigh(1,1,2)%ptr
        pblock%child(4)%ptr%neigh(1,1,p)%ptr => pblock%child(3)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(1,1,p)%ptr => pblock%neigh(1,1,3)%ptr
        pblock%child(6)%ptr%neigh(1,1,p)%ptr => pblock%child(5)%ptr
        pblock%child(7)%ptr%neigh(1,1,p)%ptr => pblock%neigh(1,1,4)%ptr
        pblock%child(8)%ptr%neigh(1,1,p)%ptr => pblock%child(7)%ptr
#endif /* NDIMS == 3 */

! X direction (right side)
!
        pblock%child(1)%ptr%neigh(1,2,p)%ptr => pblock%child(2)%ptr
        pblock%child(2)%ptr%neigh(1,2,p)%ptr => pblock%neigh(1,2,1)%ptr
        pblock%child(3)%ptr%neigh(1,2,p)%ptr => pblock%child(4)%ptr
        pblock%child(4)%ptr%neigh(1,2,p)%ptr => pblock%neigh(1,2,2)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(1,2,p)%ptr => pblock%child(6)%ptr
        pblock%child(6)%ptr%neigh(1,2,p)%ptr => pblock%neigh(1,2,3)%ptr
        pblock%child(7)%ptr%neigh(1,2,p)%ptr => pblock%child(8)%ptr
        pblock%child(8)%ptr%neigh(1,2,p)%ptr => pblock%neigh(1,2,4)%ptr
#endif /* NDIMS == 3 */

! Y direction (left side)
!
        pblock%child(1)%ptr%neigh(2,1,p)%ptr => pblock%neigh(2,1,1)%ptr
        pblock%child(2)%ptr%neigh(2,1,p)%ptr => pblock%neigh(2,1,2)%ptr
        pblock%child(3)%ptr%neigh(2,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(4)%ptr%neigh(2,1,p)%ptr => pblock%child(2)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(2,1,p)%ptr => pblock%neigh(2,1,3)%ptr
        pblock%child(6)%ptr%neigh(2,1,p)%ptr => pblock%neigh(2,1,4)%ptr
        pblock%child(7)%ptr%neigh(2,1,p)%ptr => pblock%child(5)%ptr
        pblock%child(8)%ptr%neigh(2,1,p)%ptr => pblock%child(6)%ptr
#endif /* NDIMS == 3 */

! Y direction (right side)
!
        pblock%child(1)%ptr%neigh(2,2,p)%ptr => pblock%child(3)%ptr
        pblock%child(2)%ptr%neigh(2,2,p)%ptr => pblock%child(4)%ptr
        pblock%child(3)%ptr%neigh(2,2,p)%ptr => pblock%neigh(2,2,1)%ptr
        pblock%child(4)%ptr%neigh(2,2,p)%ptr => pblock%neigh(2,2,2)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(2,2,p)%ptr => pblock%child(7)%ptr
        pblock%child(6)%ptr%neigh(2,2,p)%ptr => pblock%child(8)%ptr
        pblock%child(7)%ptr%neigh(2,2,p)%ptr => pblock%neigh(2,2,3)%ptr
        pblock%child(8)%ptr%neigh(2,2,p)%ptr => pblock%neigh(2,2,4)%ptr
#endif /* NDIMS == 3 */

#if NDIMS == 3
! Z direction (left side)
!
        pblock%child(1)%ptr%neigh(3,1,p)%ptr => pblock%neigh(3,1,1)%ptr
        pblock%child(2)%ptr%neigh(3,1,p)%ptr => pblock%neigh(3,1,2)%ptr
        pblock%child(3)%ptr%neigh(3,1,p)%ptr => pblock%neigh(3,1,3)%ptr
        pblock%child(4)%ptr%neigh(3,1,p)%ptr => pblock%neigh(3,1,4)%ptr
        pblock%child(5)%ptr%neigh(3,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(6)%ptr%neigh(3,1,p)%ptr => pblock%child(2)%ptr
        pblock%child(7)%ptr%neigh(3,1,p)%ptr => pblock%child(3)%ptr
        pblock%child(8)%ptr%neigh(3,1,p)%ptr => pblock%child(4)%ptr

! Z direction (right side)
!
        pblock%child(1)%ptr%neigh(3,2,p)%ptr => pblock%child(5)%ptr
        pblock%child(2)%ptr%neigh(3,2,p)%ptr => pblock%child(6)%ptr
        pblock%child(3)%ptr%neigh(3,2,p)%ptr => pblock%child(7)%ptr
        pblock%child(4)%ptr%neigh(3,2,p)%ptr => pblock%child(8)%ptr
        pblock%child(5)%ptr%neigh(3,2,p)%ptr => pblock%neigh(3,2,1)%ptr
        pblock%child(6)%ptr%neigh(3,2,p)%ptr => pblock%neigh(3,2,2)%ptr
        pblock%child(7)%ptr%neigh(3,2,p)%ptr => pblock%neigh(3,2,3)%ptr
        pblock%child(8)%ptr%neigh(3,2,p)%ptr => pblock%neigh(3,2,4)%ptr
#endif /* NDIMS == 3 */

      end do

! set neighbor pointers of the neighbors to the current children
!
! X direction (left side)
!
      pneigh => pblock%neigh(1,1,1)%ptr
      pchild => pblock%child(1)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(1,1,2)%ptr
      pchild => pblock%child(3)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,2,p)%ptr => pchild
          end do
        endif
      endif

#if NDIMS == 3
      pneigh => pblock%neigh(1,1,3)%ptr
      pchild => pblock%child(5)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(1,1,4)%ptr
      pchild => pblock%child(7)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,2,p)%ptr => pchild
          end do
        endif
      endif
#endif /* NDIMS == 3 */

! X direction (right side)
!
      pneigh => pblock%neigh(1,2,1)%ptr
      pchild => pblock%child(2)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(1,2,2)%ptr
      pchild => pblock%child(4)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,1,p)%ptr => pchild
          end do
        endif
      endif

#if NDIMS == 3
      pneigh => pblock%neigh(1,2,3)%ptr
      pchild => pblock%child(6)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(1,2,4)%ptr
      pchild => pblock%child(8)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(1,1,p)%ptr => pchild
          end do
        endif
      endif
#endif /* NDIMS == 3 */

! Y direction (left side)
!
      pneigh => pblock%neigh(2,1,1)%ptr
      pchild => pblock%child(1)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(2,1,2)%ptr
      pchild => pblock%child(2)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,2,p)%ptr => pchild
          end do
        endif
      endif

#if NDIMS == 3
      pneigh => pblock%neigh(2,1,3)%ptr
      pchild => pblock%child(5)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(2,1,4)%ptr
      pchild => pblock%child(6)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,2,p)%ptr => pchild
          end do
        endif
      endif
#endif /* NDIMS == 3 */

! Y direction (right side)
!
      pneigh => pblock%neigh(2,2,1)%ptr
      pchild => pblock%child(2)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(2,2,2)%ptr
      pchild => pblock%child(4)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,1,p)%ptr => pchild
          end do
        endif
      endif

#if NDIMS == 3
      pneigh => pblock%neigh(2,2,3)%ptr
      pchild => pblock%child(6)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(2,2,4)%ptr
      pchild => pblock%child(8)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(2,1,p)%ptr => pchild
          end do
        endif
      endif
#endif /* NDIMS == 3 */

#if NDIMS == 3
! Z direction (left side)
!
      pneigh => pblock%neigh(3,1,1)%ptr
      pchild => pblock%child(1)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,2,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,1,2)%ptr
      pchild => pblock%child(2)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,2,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,1,3)%ptr
      pchild => pblock%child(3)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,2,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,2,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,1,4)%ptr
      pchild => pblock%child(4)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,2,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,2,p)%ptr => pchild
          end do
        endif
      endif

! Z direction (right side)
!
      pneigh => pblock%neigh(3,2,1)%ptr
      pchild => pblock%child(5)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,1,1)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,2,2)%ptr
      pchild => pblock%child(6)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,1,2)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,2,3)%ptr
      pchild => pblock%child(7)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,1,3)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,1,p)%ptr => pchild
          end do
        endif
      endif

      pneigh => pblock%neigh(3,2,4)%ptr
      pchild => pblock%child(8)%ptr
      if (associated(pneigh)) then
        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(3,1,4)%ptr => pchild
        endif
        if (pneigh%level .eq. pchild%level) then
          do p = 1, nfaces
            pneigh%neigh(3,1,p)%ptr => pchild
          end do
        endif
      endif
#endif /* NDIMS == 3 */

! set corresponding configuration of the new blocks
!
      select case(pblock%config)
      case(0) ! 'Z'

#if NDIMS == 2
        config(:) = (/ 0, 0, 0, 0 /)
        order (:) = (/ 1, 2, 3, 4 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
        order (:) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
#endif /* NDIMS == 3 */

      case(1) ! 'N'

#if NDIMS == 2
        config(:) = (/ 2, 1, 1, 3 /)
        order (:) = (/ 1, 3, 4, 2 /)
#endif /* NDIMS == 2 */

      case(2) ! 'D'

#if NDIMS == 2
        config(:) = (/ 1, 2, 2, 4 /)
        order (:) = (/ 1, 2, 4, 3 /)
#endif /* NDIMS == 2 */

      case(3) ! 'C'

#if NDIMS == 2
        config(:) = (/ 4, 3, 3, 1 /)
        order (:) = (/ 2, 1, 3, 4 /)
#endif /* NDIMS == 2 */

      case(4) ! 'U'

#if NDIMS == 2
        config(:) = (/ 3, 4, 4, 2 /)
        order (:) = (/ 3, 1, 2, 4 /)
#endif /* NDIMS == 2 */

      end select

! set blocks configurations
!
      do p = 1, nchild
        pblock%child(p)%ptr%config = config(p)
      end do

! connect blocks in chain
!
      do p = 2, nchild
        pblock%child(order(p  ))%ptr%prev => pblock%child(order(p-1))%ptr
        pblock%child(order(p-1))%ptr%next => pblock%child(order(p  ))%ptr
      end do

! insert this chain after the parent block
!
      pneigh => pblock%next
      pfirst => pblock%child(order(     1))%ptr
      plast  => pblock%child(order(nchild))%ptr
      if (associated(pneigh)) then
        pneigh%prev => plast
        plast%next  => pneigh
      else
        last_meta => plast
        nullify(plast%next)
      endif

      pblock%next => pfirst
      pfirst%prev => pblock

! allocate data blocks if necessary
!
      if (falloc_data) then

! calculate the size of new blocks
!
        xln = 0.5 * (pblock%data%xmax - pblock%data%xmin)
        yln = 0.5 * (pblock%data%ymax - pblock%data%ymin)
#if NDIMS == 3
        zln = 0.5 * (pblock%data%zmax - pblock%data%zmin)
#else /* NDIMS == 3 */
        zln =       (pblock%data%zmax - pblock%data%zmin)
#endif /* NDIMS == 3 */

! iterate over all children and allocate data blocks
!
        do p = 1, nchild

! assign a pointer to the current child
!
          pchild => pblock%child(order(p))%ptr

! allocate data block
!
          call allocate_datablock(pdata)

! calculate block bounds
!
          i   = mod((p - 1)    ,2)
          j   = mod((p - 1) / 2,2)
          k   = mod((p - 1) / 4,2)

          xmn = pblock%data%xmin + xln * i
          ymn = pblock%data%xmin + yln * j
          zmn = pblock%data%xmin + zln * k

          xmx = xmn + xln
          ymx = ymn + yln
          zmx = zmn + zln

! set block bounds
!
          call datablock_setbounds(pdata, xmn, xmx, ymn, ymx, zmn, zmx)

! associate with the meta block
!
          call associate_blocks(pchild, pdata)

        end do

! connect blocks in chain
!
        do p = 2, nchild
          pblock%child(order(p  ))%ptr%data%prev => pblock%child(order(p-1))%ptr%data
          pblock%child(order(p-1))%ptr%data%next => pblock%child(order(p  ))%ptr%data
        end do

! insert this chain after the parent block
!
        pdata => pblock%data%next

        pfirst => pblock%child(order(     1))%ptr
        plast  => pblock%child(order(nchild))%ptr

        if (associated(pdata)) then
          pdata%prev => plast%data
          plast%data%next  => pdata
        else
          last_data => plast%data
          nullify(plast%data%next)
        endif

        pblock%data%next => pfirst%data
        pfirst%data%prev => pblock%data

      end if

! point the current block to the last created one
!
      pblock => plast

! ! depending on the configuration of the parent block
! !
!       select case(pblock%config)
!       case('z', 'Z')
!
! ! set blocks configurations
! !
!         pbl%config = 'Z'
!         pbr%config = 'Z'
!         ptl%config = 'Z'
!         ptr%config = 'Z'
!
! ! connect blocks in a chain
! !
!         pbl%next => pbr
!         pbr%next => ptl
!         ptl%next => ptr
!
!         pbr%prev => pbl
!         ptl%prev => pbr
!         ptr%prev => ptl
!
! ! insert this chain after the parent block
! !
!         pb => pblock%next
!         if (associated(pb)) then
!           pb%prev => ptr
!           ptr%next => pb
!         else
!           plast => ptr
!           nullify(ptr%next)
!         endif
!         pblock%next => pbl
!         pbl%prev => pblock
!
!         pblock => ptr
!
!       case('n', 'N')
!
! ! set blocks configurations
! !
!         pbl%config = 'D'
!         ptl%config = 'N'
!         ptr%config = 'N'
!         pbr%config = 'C'
!
! ! connect blocks in a chain
! !
!         pbl%next => ptl
!         ptl%next => ptr
!         ptr%next => pbr
!
!         ptl%prev => pbl
!         ptr%prev => ptl
!         pbr%prev => ptr
!
! ! insert this chain after the parent the block
! !
!         pb => pblock%next
!         if (associated(pb)) then
!           pb%prev => pbr
!           pbr%next => pb
!         endif
!         pbl%prev => pblock
!         pblock%next => pbl
!
!         pblock => pbr
!
!       case('d', 'D')
!
! ! set blocks configurations
! !
!         pbl%config = 'N'
!         pbr%config = 'D'
!         ptr%config = 'D'
!         ptl%config = 'U'
!
! ! connect blocks in a chain
! !
!         pbl%next => pbr
!         pbr%next => ptr
!         ptr%next => ptl
!
!         pbr%prev => pbl
!         ptr%prev => pbr
!         ptl%prev => ptr
!
! ! insert this chain in the block list
! !
!         pb => pblock%next
!         if (associated(pb)) then
!           pb%prev => ptl
!           ptl%next => pb
!         endif
!         pbl%prev => pblock
!         pblock%next => pbl
!
!         pblock => ptl
!
!       case('c', 'C')
!
! ! set blocks configurations
! !
!         ptr%config = 'U'
!         ptl%config = 'C'
!         pbl%config = 'C'
!         pbr%config = 'N'
!
! ! connect blocks in a chain
! !
!         ptr%next => ptl
!         ptl%next => pbl
!         pbl%next => pbr
!
!         ptl%prev => ptr
!         pbl%prev => ptl
!         pbr%prev => pbl
!
! ! insert this chain in the block list
! !
!         pb => pblock%next
!         if (associated(pb)) then
!           pb%prev => pbr
!           pbr%next => pb
!         endif
!         ptr%prev => pblock
!         pblock%next => ptr
!
!         pblock => pbr
!
!       case('u', 'U')
!
! ! set blocks configurations
! !
!         ptr%config = 'C'
!         pbr%config = 'U'
!         pbl%config = 'U'
!         ptl%config = 'D'
!
! ! connect blocks in a chain
! !
!         ptr%next => pbr
!         pbr%next => pbl
!         pbl%next => ptl
!
!         pbr%prev => ptr
!         pbl%prev => pbr
!         ptl%prev => pbl
!
! ! insert this chain in the block list
! !
!         pb => pblock%next
!         if (associated(pb)) then
!           pb%prev => ptl
!           ptl%next => pb
!         endif
!         ptr%prev => pblock
!         pblock%next => ptr
!
!         pblock => ptl
!
!       end select

    else

! terminate program if the pointer passed by argument is not associated
!
      call print_error("blocks::refine_blocks","Input pointer is not associated! Terminating!")
    endif

!-------------------------------------------------------------------------------
!
  end subroutine refine_block
!
!===============================================================================
!
! derefine_block: subroutine derefines selected block
!
!===============================================================================
!
  subroutine derefine_block(pblock)

#ifdef MPI
    use mpitools, only : ncpu
#endif /* MPI */

    implicit none

! input parameters
!
    type(block), pointer, intent(inout) :: pblock

! local variables
!
    integer :: p

! pointers
!
    type(block), pointer :: pb, pbl, pbr, ptl, ptr, pneigh
!
!-------------------------------------------------------------------------------
!
! prepare pointers to children
!
    pbl => get_pointer(pblock%child(1)%id)
    pbr => get_pointer(pblock%child(2)%id)
    ptl => get_pointer(pblock%child(3)%id)
    ptr => get_pointer(pblock%child(4)%id)

! prepare neighbors
!
#ifdef MPI
    pblock%neigh(1,1,1)%cpu = pbl%neigh(1,1,1)%cpu
    pblock%neigh(1,1,2)%cpu = ptl%neigh(1,1,2)%cpu
    pblock%neigh(1,2,1)%cpu = pbr%neigh(1,2,1)%cpu
    pblock%neigh(1,2,2)%cpu = ptr%neigh(1,2,2)%cpu
    pblock%neigh(2,1,1)%cpu = pbl%neigh(2,1,1)%cpu
    pblock%neigh(2,1,2)%cpu = pbr%neigh(2,1,2)%cpu
    pblock%neigh(2,2,1)%cpu = ptl%neigh(2,2,1)%cpu
    pblock%neigh(2,2,2)%cpu = ptr%neigh(2,2,2)%cpu
#endif /* MPI */
    pblock%neigh(1,1,1)%id = pbl%neigh(1,1,1)%id
    pblock%neigh(1,1,2)%id = ptl%neigh(1,1,2)%id
    pblock%neigh(1,2,1)%id = pbr%neigh(1,2,1)%id
    pblock%neigh(1,2,2)%id = ptr%neigh(1,2,2)%id
    pblock%neigh(2,1,1)%id = pbl%neigh(2,1,1)%id
    pblock%neigh(2,1,2)%id = pbr%neigh(2,1,2)%id
    pblock%neigh(2,2,1)%id = ptl%neigh(2,2,1)%id
    pblock%neigh(2,2,2)%id = ptr%neigh(2,2,2)%id

#ifdef MPI
    if (pblock%neigh(1,1,1)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(1,1,1)%id)
      if (associated(pneigh)) then
        pneigh%neigh(1,2,1)%id = pblock%id
        pneigh%neigh(1,2,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(1,1,2)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(1,1,2)%id)
      if (associated(pneigh)) then
        pneigh%neigh(1,2,1)%id = pblock%id
        pneigh%neigh(1,2,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(1,2,1)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(1,2,1)%id)
      if (associated(pneigh)) then
        pneigh%neigh(1,1,1)%id = pblock%id
        pneigh%neigh(1,1,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(1,2,2)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(1,2,2)%id)
      if (associated(pneigh)) then
        pneigh%neigh(1,1,1)%id = pblock%id
        pneigh%neigh(1,1,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(2,1,1)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(2,1,1)%id)
      if (associated(pneigh)) then
        pneigh%neigh(2,2,1)%id = pblock%id
        pneigh%neigh(2,2,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(2,1,2)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(2,1,2)%id)
      if (associated(pneigh)) then
        pneigh%neigh(2,2,1)%id = pblock%id
        pneigh%neigh(2,2,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(2,2,1)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(2,2,1)%id)
      if (associated(pneigh)) then
        pneigh%neigh(2,1,1)%id = pblock%id
        pneigh%neigh(2,1,2)%id = pblock%id
      endif
    endif

    if (pblock%neigh(2,2,2)%cpu .eq. ncpu) then
      pneigh => get_pointer(pblock%neigh(2,2,2)%id)
      if (associated(pneigh)) then
        pneigh%neigh(2,1,1)%id = pblock%id
        pneigh%neigh(2,1,2)%id = pblock%id
      endif
    endif
#else /* MPI */
    pneigh => get_pointer(pblock%neigh(1,1,1)%id)
    if (associated(pneigh)) then
      pneigh%neigh(1,2,1)%id = pblock%id
      pneigh%neigh(1,2,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(1,1,2)%id)
    if (associated(pneigh)) then
      pneigh%neigh(1,2,1)%id = pblock%id
      pneigh%neigh(1,2,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(1,2,1)%id)
    if (associated(pneigh)) then
      pneigh%neigh(1,1,1)%id = pblock%id
      pneigh%neigh(1,1,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(1,2,2)%id)
    if (associated(pneigh)) then
      pneigh%neigh(1,1,1)%id = pblock%id
      pneigh%neigh(1,1,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(2,1,1)%id)
    if (associated(pneigh)) then
      pneigh%neigh(2,2,1)%id = pblock%id
      pneigh%neigh(2,2,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(2,1,2)%id)
    if (associated(pneigh)) then
      pneigh%neigh(2,2,1)%id = pblock%id
      pneigh%neigh(2,2,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(2,2,1)%id)
    if (associated(pneigh)) then
      pneigh%neigh(2,1,1)%id = pblock%id
      pneigh%neigh(2,1,2)%id = pblock%id
    endif

    pneigh => get_pointer(pblock%neigh(2,2,2)%id)
    if (associated(pneigh)) then
      pneigh%neigh(2,1,1)%id = pblock%id
      pneigh%neigh(2,1,2)%id = pblock%id
    endif
#endif /* MPI */

! set the leaf flag for children
!
    pblock%leaf = .true.
    pbl%leaf    = .false.
    pbr%leaf    = .false.
    ptl%leaf    = .false.
    ptr%leaf    = .false.

! prepare next and prev pointers
!
    select case(pblock%config)
      case('z', 'Z')
        pb => ptr%next
        if (associated(pb)) then
          pb%prev => pblock
          pblock%next => pb
        else
          nullify(pblock%next)
        endif
      case('n', 'N')
        pb => pbr%next
        if (associated(pb)) then
          pb%prev => pblock
          pblock%next => pb
        else
          nullify(pblock%next)
        endif
      case('d', 'D')
        pb => ptl%next
        if (associated(pb)) then
          pb%prev => pblock
          pblock%next => pb
        else
          nullify(pblock%next)
        endif
      case('c', 'C')
        pb => pbr%next
        if (associated(pb)) then
          pb%prev => pblock
          pblock%next => pb
        else
          nullify(pblock%next)
        endif
      case('u', 'U')
        pb => ptl%next
        if (associated(pb)) then
          pb%prev => pblock
          pblock%next => pb
        else
          nullify(pblock%next)
        endif
    end select

    pblock%refine = 0
    pbl%refine = 0
    pbr%refine = 0
    ptl%refine = 0
    ptr%refine = 0

    call deallocate_block(pbl)
    call deallocate_block(pbr)
    call deallocate_block(ptl)
    call deallocate_block(ptr)

!-------------------------------------------------------------------------------
!
  end subroutine derefine_block

!===============================================================================
!
end module
