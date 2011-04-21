!!******************************************************************************
!!
!! module: blocks - handling block storage
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
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
  integer(kind=4), parameter :: ndims  = NDIMS
  integer(kind=4), parameter :: nsides = 2
  integer(kind=4), parameter :: nfaces = 2**(ndims-1)
  integer(kind=4), parameter :: nchild = 2**ndims

!! BLOCK STRUCTURE POINTERS (have to be defined before structures)
!!
! define pointers to block_meta and block_data structures
!
  type pointer_meta
    type(block_meta), pointer :: ptr
  end type pointer_meta

  type pointer_data
    type(block_data), pointer :: ptr
  end type pointer_data

  type pointer_info
    type(block_info), pointer :: ptr
  end type pointer_info

!! BLOCK STRUCTURES
!
! define block_meta structure
!
  type block_meta
    type(block_meta)  , pointer :: prev             ! pointer to the previous block
    type(block_meta)  , pointer :: next             ! pointer to the next block
    type(block_meta)  , pointer :: parent           ! pointer to the parent block
    type(pointer_meta)          :: child(nchild)    ! pointers to children
    type(pointer_meta)          :: neigh(ndims,nsides,nfaces)
                                                    ! pointers to neighbors

    type(block_data)  , pointer :: data             ! pointer to the data block

    integer(kind=4)             :: id               ! block identificator
    integer(kind=4)             :: cpu              ! the cpu id of the block
    integer(kind=4)             :: level            ! refinement level
    integer(kind=4)             :: config           ! configuration flag
    integer(kind=4)             :: refine           ! refinement flag:
!                                                       -1 - derefine
!                                                        0 - do nothing
!                                                        1 - refine

    integer(kind=4)             :: pos(ndims)       ! the position in the parent block
    integer(kind=4)             :: coord(ndims)     ! coordinate of the lower
!                                                     corner in the effective
!                                                     resolution

    logical                     :: leaf             ! leaf flag

    real(kind=8)                :: xmin, xmax       ! bounds for the x direction
    real(kind=8)                :: ymin, ymax       ! bounds for the y direction
    real(kind=8)                :: zmin, zmax       ! bounds for the z direction
  end type block_meta

! define block_data structure
!
  type block_data
    type(block_data), pointer :: prev             ! pointer to the previous block
    type(block_data), pointer :: next             ! pointer to the next block

    type(block_meta), pointer :: meta             ! pointer to the metadata block

    real, dimension(:,:,:,:)  , allocatable :: u  ! the array of the conserved variables
  end type block_data

! define block_info structure for boundary exchange
!
  type block_info
    type(block_info)  , pointer :: prev             ! pointer to the previous block
    type(block_info)  , pointer :: next             ! pointer to the next block
    type(block_meta)  , pointer :: block            ! pointer to the meta block
    type(block_meta)  , pointer :: neigh            ! pointer to the neighbor block
    integer(kind=4)             :: direction        ! direction of the neighbor block
    integer(kind=4)             :: side             ! side of the neighbor block
    integer(kind=4)             :: face             ! face of the neighbor block
    integer(kind=4)             :: level_difference ! the difference of levels
  end type block_info

! chains of meta blocks and data blocks
!
  type(block_meta), pointer, save :: list_meta, last_meta
  type(block_data), pointer, save :: list_data, last_data

! stored last id (should always increase)
!
  integer(kind=4)     , save :: last_id

! store number of allocated blocks and leafs
!
  integer(kind=4)     , save :: mblocks, dblocks, nleafs

! the effective resolution at all levels
!
  integer(kind=4), dimension(:), allocatable, save :: res

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
!
!-------------------------------------------------------------------------------
!
! nullify list pointers
!
    nullify(list_meta)
    nullify(list_data)

! reset the number of meta and data blocks, and leafs
!
    mblocks = 0
    dblocks = 0
    nleafs  = 0

! reset ID
!
    last_id = 0
!
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
  subroutine clear_blocks()

    implicit none

! pointers
!
    type(block_meta), pointer :: pmeta
!
!-------------------------------------------------------------------------------
!
! clear all meta blocks
!
    pmeta => list_meta
    do while(associated(pmeta))
      call deallocate_metablock(pmeta)
      pmeta => list_meta
    end do

! deallocating coordinate variables
!
    if (allocated(res))  deallocate(res)
!
!-------------------------------------------------------------------------------
!
  end subroutine clear_blocks
!
!===============================================================================
!
! increase_id: function increases the last ID by 1 and returns it
!
!===============================================================================
!
  function increase_id()

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
!
!-------------------------------------------------------------------------------
!
  end function increase_id
!
!===============================================================================
!
! allocate_metablock: subroutine allocates space for one meta block and returns
!                     the pointer to this block
!
!===============================================================================
!
  subroutine allocate_metablock(pmeta)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(out) :: pmeta

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pmeta)

! nullify pointers
!
    nullify(pmeta%prev)
    nullify(pmeta%next)
    nullify(pmeta%parent)
    nullify(pmeta%data)
    do i = 1, nchild
      nullify(pmeta%child(i)%ptr)
    end do
    do k = 1, nfaces
      do j = 1, nsides
        do i = 1, ndims
          nullify(pmeta%neigh(i,j,k)%ptr)
        end do
      end do
    end do

! set unique ID
!
    pmeta%id = increase_id()

! unset the CPU number of current block, level, the configuration, refine and
! leaf flags
!
    pmeta%cpu    = -1
    pmeta%level  = -1
    pmeta%config = -1
    pmeta%refine =  0
    pmeta%leaf   = .false.

! initialize the position in the parent block
!
    pmeta%pos(:)   = -1

! initialize the coordinate
!
    pmeta%coord(:) = 0

! initialize bounds of the block
!
    pmeta%xmin = 0.0
    pmeta%xmax = 1.0
    pmeta%ymin = 0.0
    pmeta%ymax = 1.0
    pmeta%zmin = 0.0
    pmeta%zmax = 1.0

! increase the number of allocated meta blocks
!
    mblocks = mblocks + 1
!
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
  subroutine allocate_datablock(pdata)

    use config   , only : im, jm, km
    use variables, only : nqt, nfl

    implicit none

! output arguments
!
    type(block_data), pointer, intent(out) :: pdata
!
!-------------------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pdata)

! nullify pointers
!
    nullify(pdata%prev)
    nullify(pdata%next)
    nullify(pdata%meta)

! allocate space for the conserved variables
!
    allocate(pdata%u(      nqt,im,jm,km))

! increase the number of allocated meta blocks
!
    dblocks = dblocks + 1
!
!-------------------------------------------------------------------------------
!
  end subroutine allocate_datablock
!
!===============================================================================
!
! deallocate_metablock: subroutine deallocates space ocuppied by a given metablock
!
!===============================================================================
!
  subroutine deallocate_metablock(pmeta)

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
    if (associated(pmeta)) then

! if this is the first block in the list, update the list_meta pointer
!
      if (pmeta%id .eq. list_meta%id) &
        list_meta => pmeta%next

! if this is the last block in the list, update the last_meta pointer
!
      if (pmeta%id .eq. last_meta%id) &
        last_meta => pmeta%prev

! update the pointer of previous and next blocks
!
      if (associated(pmeta%prev)) &
        pmeta%prev%next => pmeta%next

      if (associated(pmeta%next)) &
        pmeta%next%prev => pmeta%prev

! nullify children
!
      do i = 1, nchild
        nullify(pmeta%child(i)%ptr)
      end do

! nullify neighbors
!
      do k = 1, nfaces
        do j = 1, nsides
          do i = 1, ndims
            nullify(pmeta%neigh(i,j,k)%ptr)
          end do
        end do
      end do

! if corresponding data block is allocated, deallocate it too
!
      if (associated(pmeta%data)) &
        call deallocate_datablock(pmeta%data)

! nullify pointers
!
      nullify(pmeta%next)
      nullify(pmeta%prev)
      nullify(pmeta%data)
      nullify(pmeta%parent)

! free and nullify the block
!
      deallocate(pmeta)
      nullify(pmeta)

! decrease the number of allocated blocks
!
      mblocks = mblocks - 1

    end if
!
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
  subroutine deallocate_datablock(pdata)

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
    if (associated(pdata)) then

! if this is the first block in the list, update the list_data pointer
!
      if (pdata%meta%id .eq. list_data%meta%id) &
        list_data => pdata%next

! if this is the last block in the list, update the last_data pointer
!
      if (pdata%meta%id .eq. last_data%meta%id) &
        last_data => pdata%prev

! update the pointer of previous and next blocks
!
      if (associated(pdata%prev)) &
        pdata%prev%next => pdata%next

      if (associated(pdata%next)) &
        pdata%next%prev => pdata%prev

! deallocate variables
!
      deallocate(pdata%u)

! nullify pointers
!
      nullify(pdata%next)
      nullify(pdata%prev)
      nullify(pdata%meta)

! free and nullify the block
!
      deallocate(pdata)
      nullify(pdata)

! decrease the number of allocated blocks
!
      dblocks = dblocks - 1

    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine deallocate_datablock
!
!===============================================================================
!
! append_metablock: subroutine allocates space for one meta block and appends it
!                   to the meta block list
!
!===============================================================================
!
  subroutine append_metablock(pmeta)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(out) :: pmeta
!
!-------------------------------------------------------------------------------
!
! allocate block
!
    call allocate_metablock(pmeta)

! add to the list
!
    if (associated(last_meta)) then
      pmeta%prev => last_meta
      last_meta%next => pmeta
    else
      list_meta => pmeta
    end if

! set the pointer to the last block in the list
!
    last_meta => pmeta
!
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
  subroutine append_datablock(pdata)

    implicit none

! output arguments
!
    type(block_data), pointer, intent(out) :: pdata
!
!-------------------------------------------------------------------------------
!
! allocate block
!
    call allocate_datablock(pdata)

! add to the list
!
    if (associated(last_data)) then
      pdata%prev => last_data
      last_data%next => pdata
    else
      list_data => pdata
    end if

! set the pointer to the last block in the list
!
    last_data => pdata
!
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
  subroutine associate_blocks(pmeta, pdata)

    implicit none

! output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
    pmeta%data => pdata
    pdata%meta => pmeta
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
  subroutine metablock_setleaf(pmeta)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! set the leaf flag
!
    pmeta%leaf = .true.

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
  subroutine metablock_unsetleaf(pmeta)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! set the leaf flag
!
    pmeta%leaf = .false.

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
  subroutine metablock_setconfig(pmeta, config)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: config
!
!-------------------------------------------------------------------------------
!
! set the config flag
!
    pmeta%config = config
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
  subroutine metablock_setlevel(pmeta, level)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: level
!
!-------------------------------------------------------------------------------
!
! set the refinement level
!
    pmeta%level = level
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_setlevel
!
!===============================================================================
!
! metablock_set_position: subroutine sets the position of the meta block in the
!                         parent block
!
!===============================================================================
!
  subroutine metablock_set_position(pmeta, px, py, pz)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: px, py, pz
!
!-------------------------------------------------------------------------------
!
! set the position in the parent block
!
    pmeta%pos(1) = px
    pmeta%pos(2) = py
#if NDIMS == 3
    pmeta%pos(3) = pz
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_position
!
!===============================================================================
!
! metablock_set_coord: subroutine sets the coordinates of the meta block
!
!===============================================================================
!
  subroutine metablock_set_coord(pmeta, px, py, pz)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: px, py, pz
!
!-------------------------------------------------------------------------------
!
! set the coordintaes
!
    pmeta%coord(1) = px
    pmeta%coord(2) = py
#if NDIMS == 3
    pmeta%coord(3) = pz
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_coord
!
!===============================================================================
!
! metablock_setbounds: subroutine sets the bounds of data block
!
!===============================================================================
!
  subroutine metablock_setbounds(pmeta, xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    real                     , intent(in)    :: xmin, xmax, ymin, ymax, zmin, zmax
!
!-------------------------------------------------------------------------------
!
! set bounds of the block
!
    pmeta%xmin = xmin
    pmeta%xmax = xmax
    pmeta%ymin = ymin
    pmeta%ymax = ymax
    pmeta%zmin = zmin
    pmeta%zmax = zmax
!
!-------------------------------------------------------------------------------
!
  end subroutine metablock_setbounds
!
!===============================================================================
!
! refine_block: subroutine refines selected block
!
!===============================================================================
!
  subroutine refine_block(pblock, falloc_data)

    use error   , only : print_error

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
    integer, dimension(nchild)              :: config, order
    integer, dimension(ndims,nsides,nfaces) :: set

! local variables
!
    integer :: p, i, j, k, ic, jc, kc
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
! interior of the block
!
      do p = 1, nfaces

! X direction (left side)
!
        pblock%child(2)%ptr%neigh(1,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(4)%ptr%neigh(1,1,p)%ptr => pblock%child(3)%ptr
#if NDIMS == 3
        pblock%child(6)%ptr%neigh(1,1,p)%ptr => pblock%child(5)%ptr
        pblock%child(8)%ptr%neigh(1,1,p)%ptr => pblock%child(7)%ptr
#endif /* NDIMS == 3 */
        pneigh => pblock%neigh(1,1,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(1)%ptr%neigh(1,1,p)%ptr => pblock%child(2)%ptr
            pblock%child(3)%ptr%neigh(1,1,p)%ptr => pblock%child(4)%ptr
#if NDIMS == 3
            pblock%child(5)%ptr%neigh(1,1,p)%ptr => pblock%child(6)%ptr
            pblock%child(7)%ptr%neigh(1,1,p)%ptr => pblock%child(8)%ptr
#endif /* NDIMS == 3 */
          end if
        end if

! X direction (right side)
!
        pblock%child(1)%ptr%neigh(1,2,p)%ptr => pblock%child(2)%ptr
        pblock%child(3)%ptr%neigh(1,2,p)%ptr => pblock%child(4)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(1,2,p)%ptr => pblock%child(6)%ptr
        pblock%child(7)%ptr%neigh(1,2,p)%ptr => pblock%child(8)%ptr
#endif /* NDIMS == 3 */
        pneigh => pblock%neigh(1,2,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(2)%ptr%neigh(1,2,p)%ptr => pblock%child(1)%ptr
            pblock%child(4)%ptr%neigh(1,2,p)%ptr => pblock%child(3)%ptr
#if NDIMS == 3
            pblock%child(6)%ptr%neigh(1,2,p)%ptr => pblock%child(5)%ptr
            pblock%child(8)%ptr%neigh(1,2,p)%ptr => pblock%child(7)%ptr
#endif /* NDIMS == 3 */
          end if
        end if

! Y direction (left side)
!
        pblock%child(3)%ptr%neigh(2,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(4)%ptr%neigh(2,1,p)%ptr => pblock%child(2)%ptr
#if NDIMS == 3
        pblock%child(7)%ptr%neigh(2,1,p)%ptr => pblock%child(5)%ptr
        pblock%child(8)%ptr%neigh(2,1,p)%ptr => pblock%child(6)%ptr
#endif /* NDIMS == 3 */
        pneigh => pblock%neigh(2,1,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(1)%ptr%neigh(2,1,p)%ptr => pblock%child(3)%ptr
            pblock%child(2)%ptr%neigh(2,1,p)%ptr => pblock%child(4)%ptr
#if NDIMS == 3
            pblock%child(5)%ptr%neigh(2,1,p)%ptr => pblock%child(7)%ptr
            pblock%child(6)%ptr%neigh(2,1,p)%ptr => pblock%child(8)%ptr
#endif /* NDIMS == 3 */
          end if
        end if

! Y direction (right side)
!
        pblock%child(1)%ptr%neigh(2,2,p)%ptr => pblock%child(3)%ptr
        pblock%child(2)%ptr%neigh(2,2,p)%ptr => pblock%child(4)%ptr
#if NDIMS == 3
        pblock%child(5)%ptr%neigh(2,2,p)%ptr => pblock%child(7)%ptr
        pblock%child(6)%ptr%neigh(2,2,p)%ptr => pblock%child(8)%ptr
#endif /* NDIMS == 3 */
        pneigh => pblock%neigh(2,2,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(3)%ptr%neigh(2,2,p)%ptr => pblock%child(1)%ptr
            pblock%child(4)%ptr%neigh(2,2,p)%ptr => pblock%child(2)%ptr
#if NDIMS == 3
            pblock%child(7)%ptr%neigh(2,2,p)%ptr => pblock%child(5)%ptr
            pblock%child(8)%ptr%neigh(2,2,p)%ptr => pblock%child(6)%ptr
#endif /* NDIMS == 3 */
          end if
        end if

#if NDIMS == 3
! Z direction (left side)
!
        pblock%child(5)%ptr%neigh(3,1,p)%ptr => pblock%child(1)%ptr
        pblock%child(6)%ptr%neigh(3,1,p)%ptr => pblock%child(2)%ptr
        pblock%child(7)%ptr%neigh(3,1,p)%ptr => pblock%child(3)%ptr
        pblock%child(8)%ptr%neigh(3,1,p)%ptr => pblock%child(4)%ptr
        pneigh => pblock%neigh(3,1,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(1)%ptr%neigh(3,1,p)%ptr => pblock%child(5)%ptr
            pblock%child(2)%ptr%neigh(3,1,p)%ptr => pblock%child(6)%ptr
            pblock%child(3)%ptr%neigh(3,1,p)%ptr => pblock%child(7)%ptr
            pblock%child(4)%ptr%neigh(3,1,p)%ptr => pblock%child(8)%ptr
          end if
        end if

! Z direction (right side)
!
        pblock%child(1)%ptr%neigh(3,2,p)%ptr => pblock%child(5)%ptr
        pblock%child(2)%ptr%neigh(3,2,p)%ptr => pblock%child(6)%ptr
        pblock%child(3)%ptr%neigh(3,2,p)%ptr => pblock%child(7)%ptr
        pblock%child(4)%ptr%neigh(3,2,p)%ptr => pblock%child(8)%ptr
        pneigh => pblock%neigh(3,2,1)%ptr
        if (associated(pneigh)) then
          if (pneigh%id .eq. pblock%id) then
            pblock%child(5)%ptr%neigh(3,2,p)%ptr => pblock%child(1)%ptr
            pblock%child(6)%ptr%neigh(3,2,p)%ptr => pblock%child(2)%ptr
            pblock%child(7)%ptr%neigh(3,2,p)%ptr => pblock%child(3)%ptr
            pblock%child(8)%ptr%neigh(3,2,p)%ptr => pblock%child(4)%ptr
          end if
        end if
#endif /* NDIMS == 3 */
      end do

! prepare set array
!
#if NDIMS == 2
      set(1,1,:) = (/ 1, 3 /)
      set(1,2,:) = (/ 2, 4 /)
      set(2,1,:) = (/ 1, 2 /)
      set(2,2,:) = (/ 3, 4 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      set(1,1,:) = (/ 1, 3, 5, 7 /)
      set(1,2,:) = (/ 2, 4, 6, 8 /)
      set(2,1,:) = (/ 1, 2, 5, 6 /)
      set(2,2,:) = (/ 3, 4, 7, 8 /)
      set(3,1,:) = (/ 1, 2, 3, 4 /)
      set(3,2,:) = (/ 5, 6, 7, 8 /)
#endif /* NDIMS == 3 */

! set pointers to neighbors and update neighbors pointers
!
      do i = 1, ndims
        do j = 1, nsides
          do k = 1, nfaces
            pneigh => pblock%neigh(i,j,k)%ptr

            if (associated(pneigh)) then
              if (pneigh%id .ne. pblock%id) then

! point to the right neighbor
!
                do p = 1, nfaces
                  pblock%child(set(i,j,k))%ptr%neigh(i,j,p)%ptr => pneigh
                end do

! neighbor level is the same as the refined block
!
                if (pneigh%level .eq. pblock%level) then
                  pneigh%neigh(i,3-j,k)%ptr => pblock%child(set(i,j,k))%ptr
                end if

! neighbor level is the same as the child block
!
                if (pneigh%level .gt. pblock%level) then
                  do p = 1, nfaces
                    pneigh%neigh(i,3-j,p)%ptr => pblock%child(set(i,j,k))%ptr
                  end do
                end if

              end if

            end if

          end do
        end do
      end do

! set corresponding configuration of the new blocks
!
      select case(pblock%config)
      case(0)

#if NDIMS == 2
        config(:) = (/ 0, 0, 0, 0 /)
        order (:) = (/ 1, 2, 3, 4 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
        order (:) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
#endif /* NDIMS == 3 */

      case(12)

#if NDIMS == 2
        config(:) = (/ 13, 12, 12, 42 /)
        order (:) = (/  1,  3,  4,  2 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 13, 15, 15, 78, 78, 62, 62, 42 /)
        order (:) = (/  1,  3,  7,  5,  6,  8,  4,  2 /)
#endif /* NDIMS == 3 */

      case(13)

#if NDIMS == 2
        config(:) = (/ 12, 13, 13, 43 /)
        order (:) = (/  1,  2,  4,  3 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 15, 12, 12, 68, 68, 43, 43, 73 /)
        order (:) = (/  1,  5,  6,  2,  4,  8,  7,  3 /)

      case(15)

        config(:) = (/ 12, 13, 13, 48, 48, 75, 75, 65 /)
        order (:) = (/  1,  2,  4,  3,  7,  8,  6,  5 /)
#endif /* NDIMS == 3 */

      case(42)

#if NDIMS == 2
        config(:) = (/ 43, 42, 42, 12 /)
        order (:) = (/  4,  3,  1,  2 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 48, 43, 43, 75, 75, 12, 12, 62 /)
        order (:) = (/  4,  8,  7,  3,  1,  5,  6,  2 /)
#endif /* NDIMS == 3 */

      case(43)

#if NDIMS == 2
        config(:) = (/ 42, 43, 43, 13 /)
        order (:) = (/  4,  2,  1,  3 /)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        config(:) = (/ 42, 48, 48, 65, 65, 73, 73, 13 /)
        order (:) = (/  4,  2,  6,  8,  7,  5,  1,  3 /)
#endif /* NDIMS == 3 */

#if NDIMS == 3
      case(48)

        config(:) = (/ 43, 42, 42, 15, 15, 68, 68, 78 /)
        order (:) = (/  4,  3,  1,  2,  6,  5,  7,  8 /)

      case(62)

        config(:) = (/ 65, 68, 68, 73, 73, 42, 42, 12 /)
        order (:) = (/  6,  5,  7,  8,  4,  3,  1,  2 /)

      case(65)

        config(:) = (/ 68, 62, 62, 43, 43, 15, 15, 75 /)
        order (:) = (/  6,  8,  4,  2,  1,  3,  7,  5 /)

      case(68)

        config(:) = (/ 62, 65, 65, 13, 13, 78, 78, 48 /)
        order (:) = (/  6,  2,  1,  5,  7,  3,  4 , 8 /)

      case(73)

        config(:) = (/ 78, 75, 75, 62, 62, 13, 13, 43 /)
        order (:) = (/  7,  8,  6,  5,  1,  2,  4,  3 /)

      case(75)

        config(:) = (/ 73, 78, 78, 42, 42, 65, 65, 15 /)
        order (:) = (/  7,  3,  4,  8,  6,  2,  1,  5 /)

      case(78)

        config(:) = (/ 75, 73, 73, 12, 12, 48, 48, 68 /)
        order (:) = (/  7,  5,  1,  3,  4,  2,  6,  8 /)
#endif /* NDIMS == 3 */

      end select

! set blocks configurations
!
      do p = 1, nchild
        pblock%child(order(p))%ptr%config = config(p)
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
      end if

      pblock%next => pfirst
      pfirst%prev => pblock

! calculate the size of new blocks
!
      xln = 0.5 * (pblock%xmax - pblock%xmin)
      yln = 0.5 * (pblock%ymax - pblock%ymin)
#if NDIMS == 3
      zln = 0.5 * (pblock%zmax - pblock%zmin)
#else /* NDIMS == 3 */
      zln =       (pblock%zmax - pblock%zmin)
#endif /* NDIMS == 3 */

! iterate over all children and allocate data blocks
!
      do p = 1, nchild

! assign a pointer to the current child
!
        pchild => pblock%child(p)%ptr

! calculate the block position indices
!
        i   = mod((p - 1)    ,2)
        j   = mod((p - 1) / 2,2)
        k   = mod((p - 1) / 4,2)

! set the block position
!
        call metablock_set_position(pchild, i, j, k)

! set the block coordinates
!
        ic  = pblock%coord(1) + i * res(pchild%level)
        jc  = pblock%coord(2) + j * res(pchild%level)
#if NDIMS == 3
        kc  = pblock%coord(3) + k * res(pchild%level)
#endif /* NDIMS == 3 */
        call metablock_set_coord(pchild, ic, jc, kc)

! calculate block bounds
!
        xmn = pblock%xmin + xln * i
        ymn = pblock%ymin + yln * j
        zmn = pblock%zmin + zln * k

        xmx = xmn + xln
        ymx = ymn + yln
        zmx = zmn + zln

! set block bounds
!
        call metablock_setbounds(pchild, xmn, xmx, ymn, ymx, zmn, zmx)

      end do

! allocate data blocks if necessary
!
      if (falloc_data) then

! iterate over all children and allocate data blocks
!
        do p = 1, nchild

! assign a pointer to the current child
!
          pchild => pblock%child(p)%ptr

! allocate data block
!
          call allocate_datablock(pdata)

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
        end if

        pblock%data%next => pfirst%data
        pfirst%data%prev => pblock%data

      end if

! point the current block to the last created one
!
      pblock => plast

    else

! terminate program if the pointer passed by argument is not associated
!
      call print_error("blocks::refine_blocks","Input pointer is not associated! Terminating!")
    end if
!
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

    implicit none

! input parameters
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k, l, p

! local arrays
!
    integer, dimension(ndims, nsides, nfaces) :: arr

! local pointers
!
    type(block_meta), pointer :: pchild, pneigh
!
!-------------------------------------------------------------------------------
!
! prepare reference array
!
#if NDIMS == 3
    arr(1,1,:) = (/ 1, 3, 5, 7 /)
    arr(1,2,:) = (/ 2, 4, 6, 8 /)
    arr(2,1,:) = (/ 1, 2, 5, 6 /)
    arr(2,2,:) = (/ 3, 4, 7, 8 /)
    arr(3,1,:) = (/ 1, 2, 3, 4 /)
    arr(3,2,:) = (/ 5, 6, 7, 8 /)
#else /* NDIMS == 3 */
    arr(1,1,:) = (/ 1, 3 /)
    arr(1,2,:) = (/ 2, 4 /)
    arr(2,1,:) = (/ 1, 2 /)
    arr(2,2,:) = (/ 3, 4 /)
#endif /* NDIMS == 3 */

! iterate over all boundaries of the parent block
!
    do i = 1, ndims
      do j = 1, nsides
        do k = 1, nfaces

! calculate the right child number
!
          p = arr(i,j,k)

! assign the pointer to the current neighbor
!
          pneigh => pblock%child(p)%ptr%neigh(i,j,k)%ptr

! assign the right neighbor to the current neighbor pointer
!
          pblock%neigh(i,j,k)%ptr => pneigh

! update the neighbor fields of neighbors
!
          if (associated(pneigh)) then
            l = 3 - j
            do p = 1, nfaces
              pneigh%neigh(i,l,p)%ptr => pblock
            end do
          end if

        end do
      end do
    end do

! deallocate child blocks
!
    do p = 1, nchild
      call metablock_unsetleaf(pblock%child(p)%ptr)
      call deallocate_metablock(pblock%child(p)%ptr)
    end do

! set the leaf flag of parent block
!
    call metablock_setleaf(pblock)

! reset the refinement flag of the parent block
!
    pblock%refine = 0
!
!-------------------------------------------------------------------------------
!
  end subroutine derefine_block
#ifdef DEBUG
!
!===============================================================================
!
! check_metablock: subroutine checks if the meta block has proper structure
!
!===============================================================================
!
  subroutine check_metablock(pblock, string)

    implicit none

! input parameters
!
    type(block_meta), pointer, intent(in) :: pblock
    character(len=*)         , intent(in) :: string

! local variables
!
    integer :: p, i, j, k

! local pointers
!
    type(block_meta), pointer :: ptemp
!
!-------------------------------------------------------------------------------
!
! check block ID
!
    ptemp => pblock
    if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
      print *, ''
      print *, ''
      print *, trim(string)
      print *, 'wrong meta block id = ', ptemp%id
      stop
    end if

! check prev ID
!
    ptemp => pblock%prev
    if (associated(ptemp)) then
      if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
        print *, ''
        print *, ''
        print *, trim(string)
        print *, 'wrong previous block id = ', ptemp%id, pblock%id
        stop
      end if
    end if

! check next ID
!
    ptemp => pblock%next
    if (associated(ptemp)) then
      if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
        print *, ''
        print *, ''
        print *, trim(string)
        print *, 'wrong next block id = ', ptemp%id, pblock%id
        stop
      end if
    end if

! check parent ID
!
    ptemp => pblock%parent
    if (associated(ptemp)) then
      if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
        print *, ''
        print *, ''
        print *, trim(string)
        print *, 'wrong parent block id = ', ptemp%id, pblock%id
        stop
      end if
    end if

! check children IDs
!
    do p = 1, nchild
      ptemp => pblock%child(p)%ptr
      if (associated(ptemp)) then
        if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
          print *, ''
          print *, ''
          print *, trim(string)
          print *, 'wrong child block id = ', ptemp%id, pblock%id, p
          stop
        end if
      end if
    end do

! check neighbors IDs
!
    do i = 1, ndims
      do j = 1, nsides
        do k = 1, nfaces
          ptemp => pblock%neigh(i,j,k)%ptr
          if (associated(ptemp)) then
            if (ptemp%id .le. 0 .or. ptemp%id .gt. last_id) then
              print *, ''
              print *, ''
              print *, trim(string)
              print *, 'wrong neighbor id = ', ptemp%id, pblock%id, i, j, k
              stop
            end if
          end if
        end do
      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine check_metablock
#endif /* DEBUG */

!===============================================================================
!
end module
