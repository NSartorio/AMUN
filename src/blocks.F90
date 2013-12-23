!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: BLOCKS
!!
!!  This module allocates, deallocates, and handles blocks of the adaptive mesh
!!  structures.
!!
!!******************************************************************************
!
module blocks

! module variables are not implicit by default
!
  implicit none

! module parameters
!
  integer(kind=4), parameter :: ndims  = NDIMS
  integer(kind=4), parameter :: nsides = 2
  integer(kind=4), parameter :: nfaces = 2**(ndims - 1)
  integer(kind=4), parameter :: nchild = 2**ndims

!! BLOCK STRUCTURE POINTERS (they have to be defined before block structures)
!!
! define pointers to meta, data, and info block structures
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
! define the META block structure
!
  type block_meta
                                 ! pointers to the previous and next meta blocks
                                 !
    type(block_meta)  , pointer :: prev, next

                                 ! a pointer to the parent meta block
                                 !
    type(block_meta)  , pointer :: parent

                                 ! pointers to child meta blocks
                                 !
    type(pointer_meta)          :: child(nchild)

                                 ! pointers to neighbor meta blocks
                                 !
    type(pointer_meta)          :: neigh(ndims,nsides,nfaces)

                                 ! a pointer to the associated data block
                                 !
    type(block_data)  , pointer :: data

                                 ! the block identification
                                 !
    integer(kind=4)             :: id

                                 ! the number of associated cpu
                                 !
    integer(kind=4)             :: cpu

                                 ! the level of refinement
                                 !
    integer(kind=4)             :: level

                                 ! the configuration flag for its children order
                                 !
    integer(kind=4)             :: config

                                 ! the refinement flag, -1, 0, and 1 for
                                 ! derefinement, no change, and refinement,
                                 ! respectively
                                 !
    integer(kind=4)             :: refine

                                 ! the position of the block in its parent block
                                 !
    integer(kind=4)             :: pos(ndims)

                                 ! the coordinate of the lower corner of the
                                 ! block in the effective resolution units
                                 !
    integer(kind=4)             :: coord(ndims)

                                 ! the leaf flag
                                 !
    logical                     :: leaf

                                 ! the block coordinates in the physical units
                                 !
    real                        :: xmin, xmax, ymin, ymax, zmin, zmax

  end type block_meta

! define the DATA block structure
!
  type block_data
                                 ! pointers to the previous and next data blocks
                                 !
    type(block_data), pointer :: prev, next

                                 ! a pointer to the associated meta block
                                 !
    type(block_meta), pointer :: meta

                                 ! a pointer to the array conserved variables
                                 !
    real, dimension(:,:,:,:)  , pointer     :: u

                                 ! an allocatable arrays to store all conserved
                                 ! variables
                                 !
    real, dimension(:,:,:,:)  , allocatable :: u0, u1

                                 ! an allocatable array to store all primitive
                                 ! variables
                                 !
    real, dimension(:,:,:,:)  , allocatable :: q

                                 ! an allocatable array to store all fluxes
                                 !
    real, dimension(:,:,:,:,:), allocatable :: f

#ifdef DEBUG
                                 ! an allocatable array to store refinement
                                 ! values
                                 !
    real, dimension(:,:,:)    , allocatable :: c
#endif /* DEBUG */

  end type block_data

! define the INFO block structure
!
  type block_info
                                 ! pointers to the previous and next info blocks
                                 !
    type(block_info)  , pointer :: prev, next

                                 ! a pointer to the associated meta block
                                 !
    type(block_meta)  , pointer :: block

                                 ! a pointer to the associated neighbor block
                                 !
    type(block_meta)  , pointer :: neigh

                                 ! the direction, side and face indices
                                 !
    integer(kind=4)             :: direction, side, face

                                 ! the level difference between the block and
                                 ! its neighbor
                                 !
    integer(kind=4)             :: level_difference

  end type block_info

!! POINTER TO THE FIST AND LAST BLOCKS IN THE LISTS
!!
! chains of meta blocks and data blocks
!
  type(block_meta), pointer, save :: list_meta, last_meta
  type(block_data), pointer, save :: list_data, last_data

!! MODULE VARIABLES
!!
! the identification of the last allocated block (should always increase)
!
  integer(kind=4)     , save :: last_id

! the numbers of allocated meta and data blocks, and leafs
!
  integer(kind=4)     , save :: mblocks, dblocks, nleafs

! the numbers of variables and fluxes stored in data blocks
!
  integer(kind=4)     , save :: nvars, nflux

! the spacial dimensions of data block allocatable arrays
!
  integer(kind=4)     , save :: nx, ny, nz

! all variables and subroutines are private by default
!
  private

! declare public subroutines
!
  public :: pointer_meta, pointer_info
  public :: block_meta, block_data, block_info
  public :: list_meta, list_data
  public :: nchild, ndims, nsides, nfaces
  public :: initialize_blocks, finalize_blocks
  public :: set_last_id, get_last_id, get_mblocks, get_dblocks, get_nleafs
  public :: link_blocks, unlink_blocks
  public :: append_metablock
  public :: allocate_datablock, deallocate_datablock
  public :: append_datablock, remove_datablock
  public :: metablock_set_id, metablock_set_cpu, metablock_set_refine          &
          , metablock_set_config, metablock_set_level, metablock_set_position  &
          , metablock_set_coord, metablock_set_bounds, metablock_set_leaf
  public :: datablock_set_dims
  public :: refine_block, derefine_block
#ifdef DEBUG
  public :: check_metablock
#endif /* DEBUG */

  contains
!
!!==============================================================================
!!
!! INITIALIZATION/FINALIZATION SUBROUTINES
!!
!
!===============================================================================
!
! subroutine INITIALIZE_BLOCKS:
! ----------------------------
!
!   Subroutine initializes the variables related to the elementary block of
!   the adaptive structure.
!
!===============================================================================
!
  subroutine initialize_blocks()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! nullify list pointers
!
    nullify(list_meta)
    nullify(list_data)
    nullify(last_meta)
    nullify(last_data)

! reset the number of meta blocks, data blocks, and leafs
!
    mblocks = 0
    dblocks = 0
    nleafs  = 0

! set the initial number of variables and fluxes
!
    nvars   = 1
    nflux   = 0

! set the initial data block resolution
!
    nx      = 1
    ny      = 1
    nz      = 1

! reset identification counter
!
    last_id = 0

!-------------------------------------------------------------------------------
!
  end subroutine initialize_blocks
!
!===============================================================================
!
! subroutine FINALIZE_BLOCKS:
! --------------------------
!
!   Subroutine iterates over all meta blocks and first deallocates all
!   associated with them data blocks, and then their metadata structure.
!
!===============================================================================
!
  subroutine finalize_blocks()

! local variables are not implicit by default
!
    implicit none

! a pointer to the current meta block
!
    type(block_meta), pointer :: pmeta
!
!-------------------------------------------------------------------------------
!
! assiociate pmeta pointer with the first block in the list
!
    pmeta => list_meta

    do while(associated(pmeta))

! deallocate current meta block
!
      call deallocate_metablock(pmeta)

! associate pmeta pointer with the next meta block in the list
!
      pmeta => list_meta

    end do

!-------------------------------------------------------------------------------
!
  end subroutine finalize_blocks
!
!!==============================================================================
!!
!! META BLOCK SUBROUTINES
!!
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

!-------------------------------------------------------------------------------
!
  end subroutine allocate_metablock
!
!===============================================================================
!
! deallocate_metablock: subroutine deallocates space occupied by a given meta
!                       block
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
        call remove_datablock(pmeta%data)

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

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_metablock
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

!-------------------------------------------------------------------------------
!
  end subroutine append_metablock
!
!!==============================================================================
!!
!! TOOL SUBROUTINES
!!
!
!===============================================================================
!
! increase_id: function increases the last identification by 1 and returns its
!              value
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

!-------------------------------------------------------------------------------
!
  end function increase_id
!
!===============================================================================
!
! set_last_id: subroutine sets the last identification value
!
!===============================================================================
!
  subroutine set_last_id(id)

    use error, only : print_error

    implicit none

! input argument
!
    integer(kind=4), intent(in) :: id
!
!-------------------------------------------------------------------------------
!
    if (last_id .gt. id) then
      call print_error("blocks::set_last_id"                                   &
                                 , "New last_id must be larger than old one!")
    else
      last_id = id
    end if

!-------------------------------------------------------------------------------
!
  end subroutine set_last_id
!
!===============================================================================
!
! get_last_id: function returns the last identification value
!
!===============================================================================
!
  function get_last_id()

    implicit none

! return variable
!
    integer(kind=4) :: get_last_id
!
!-------------------------------------------------------------------------------
!
    get_last_id = last_id

    return

!-------------------------------------------------------------------------------
!
  end function get_last_id
!
!===============================================================================
!
! get_mblocks: function returns the number of meta blocks
!
!===============================================================================
!
  function get_mblocks()

    implicit none

! return variable
!
    integer(kind=4) :: get_mblocks
!
!-------------------------------------------------------------------------------
!
    get_mblocks = mblocks

    return

!-------------------------------------------------------------------------------
!
  end function get_mblocks
!
!===============================================================================
!
! get_dblocks: function returns the number of data blocks
!
!===============================================================================
!
  function get_dblocks()

    implicit none

! return variable
!
    integer(kind=4) :: get_dblocks
!
!-------------------------------------------------------------------------------
!
    get_dblocks = dblocks

    return

!-------------------------------------------------------------------------------
!
  end function get_dblocks
!
!===============================================================================
!
! get_nleafs: function returns the number of leafs
!
!===============================================================================
!
  function get_nleafs()

    implicit none

! return variable
!
    integer(kind=4) :: get_nleafs
!
!-------------------------------------------------------------------------------
!
    get_nleafs = nleafs

    return

!-------------------------------------------------------------------------------
!
  end function get_nleafs
!
!===============================================================================
!
! subroutine LINK_BLOCKS:
! ----------------------
!
!   Subroutine links a data block to meta block.
!
!   Arguments:
!
!     pmeta - the meta block pointer;
!     pdata - the data block pointer;
!
!===============================================================================
!
  subroutine link_blocks(pmeta, pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! set the pointers
!
    pmeta%data => pdata
    pdata%meta => pmeta

!-------------------------------------------------------------------------------
!
  end subroutine link_blocks
!
!===============================================================================
!
! subroutine UNLINK_BLOCKS:
! ------------------------
!
!   Subroutine unlinks meta block with a data block linked to it.
!
!   Arguments:
!
!     pmeta - the meta block pointer;
!     pdata - the data block pointer;
!
!===============================================================================
!
  subroutine unlink_blocks(pmeta, pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! clean up the block pointers
!
    nullify(pmeta%data)
    nullify(pdata%meta)

!-------------------------------------------------------------------------------
!
  end subroutine unlink_blocks
!
!===============================================================================
!
! metablock_set_id: subroutine sets the identification value
!
!===============================================================================
!
  subroutine metablock_set_id(pmeta, id)

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: id
!
!-------------------------------------------------------------------------------
!
! set the id field
!
    pmeta%id = id

! check if the id is larger then last_id, if so reset last_id to id
!
    if (last_id .lt. id) last_id = id

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_id
!
!===============================================================================
!
! metablock_set_cpu: subroutine sets the cpu number
!
!===============================================================================
!
  subroutine metablock_set_cpu(pmeta, cpu)

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: cpu
!
!-------------------------------------------------------------------------------
!
! set the cpu field
!
    pmeta%cpu = cpu

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_cpu
!
!===============================================================================
!
! metablock_set_refine: subroutine sets the refine flag
!
!===============================================================================
!
  subroutine metablock_set_refine(pmeta, refine)

    use error, only : print_error

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: refine
!
!-------------------------------------------------------------------------------
!
! check if the refine value is correct
!
    if (abs(refine) .gt. 1) then

! print error about wrong refine flag
!
      call print_error("blocks::metablock_set_refine"                          &
                                            , "New refine flag is incorrect!")

    else

! set the refine field
!
      pmeta%refine = refine

    end if

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_refine
!
!===============================================================================
!
! metablock_set_leaf: subroutine marks the block as a leaf
!
!===============================================================================
!
  subroutine metablock_set_leaf(pmeta)

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

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_leaf
!
!===============================================================================
!
! metablock_unset_leaf: subroutine unmarks the block as a leaf
!
!===============================================================================
!
  subroutine metablock_unset_leaf(pmeta)

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

!-------------------------------------------------------------------------------
!
  end subroutine metablock_unset_leaf
!
!===============================================================================
!
! metablock_set_config: subroutine sets the configuration flag
!
!===============================================================================
!
  subroutine metablock_set_config(pmeta, config)

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

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_config
!
!===============================================================================
!
! metablock_set_level: subroutine sets the level of data block
!
!===============================================================================
!
  subroutine metablock_set_level(pmeta, level)

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

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_level
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
! set the coordinates
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
! metablock_set_bounds: subroutine sets the bounds of data block
!
!===============================================================================
!
  subroutine metablock_set_bounds(pmeta, xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

! input/output arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    real                     , intent(in)    :: xmin, xmax
    real                     , intent(in)    :: ymin, ymax
    real                     , intent(in)    :: zmin, zmax
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

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_bounds
!
!!==============================================================================
!!
!! DATA BLOCK SUBROUTINES
!!
!
!===============================================================================
!
! allocate_datablock: subroutine allocates space for one data block and returns
!                     the pointer to this block
!
!===============================================================================
!
!===============================================================================
!
! subroutine ALLOCATE_DATABLOCK:
! -----------------------------
!
!   Subroutine allocates space for one data block and returns a pointer
!   associated with it.
!
!   Arguments:
!
!     pdata - the pointer associated with the created data block;
!
!===============================================================================
!
  subroutine allocate_datablock(pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(out) :: pdata
!
!-------------------------------------------------------------------------------
!
! allocate the block structure
!
    allocate(pdata)

! nullify all pointers
!
    nullify(pdata%prev)
    nullify(pdata%next)
    nullify(pdata%meta)

! allocate the space for conserved variables
!
    allocate(pdata%u0(nvars,nx,ny,nz))
    allocate(pdata%u1(nvars,nx,ny,nz))

! allocate the space for primitive variables
!
    allocate(pdata%q(nvars,nx,ny,nz))

! initiate the conserved variable pointer
!
    pdata%u => pdata%u0

! allocate the space for numerical fluxes
!
    if (nflux > 0) allocate(pdata%f(ndims,nflux,nx,ny,nz))

#ifdef DEBUG
! allocate the space for the refinement criterion array
!
    allocate(pdata%c(nx,ny,nz))
#endif /* DEBUG */

! increase the number of allocated meta blocks
!
    dblocks = dblocks + 1

!-------------------------------------------------------------------------------
!
  end subroutine allocate_datablock
!
!===============================================================================
!
! subroutine DEALLOCATE_DATABLOCK:
! -------------------------------
!
!   Subroutine deallocates space of the data block associated with the input
!   pointer.
!
!   Arguments:
!
!     pdata - the pointer pointing to the data block for deallocating;
!
!===============================================================================
!
  subroutine deallocate_datablock(pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! check if the input pointer is associated with a data block
!
    if (associated(pdata)) then

! deallocate conservative variables
!
      if (allocated(pdata%u0)) deallocate(pdata%u0)
      if (allocated(pdata%u1)) deallocate(pdata%u1)

! deallocate primitive variables
!
      if (allocated(pdata%q)) deallocate(pdata%q)

! deallocate numerical fluxes
!
      if (allocated(pdata%f)) deallocate(pdata%f)

#ifdef DEBUG
! deallocate the refinement critarion array
!
      if (allocated(pdata%c)) deallocate(pdata%c)
#endif /* DEBUG */

! nullify pointers
!
      nullify(pdata%u)
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

    end if ! pdata associated with a data block

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_datablock
!
!===============================================================================
!
! subroutine APPEND_DATABLOCK:
! ---------------------------
!
!   Subroutine allocates space for one data block and appends it to the data
!   block list returning a pointer associated with it.
!
!   Arguments:
!
!     pdata - the pointer associated with the created data block;
!
!===============================================================================
!
  subroutine append_datablock(pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(out) :: pdata
!
!-------------------------------------------------------------------------------
!
! allocate the data block
!
    call allocate_datablock(pdata)

! add the allocated block to the data block list
!
    if (associated(last_data)) then
      pdata%prev     => last_data
      last_data%next => pdata
    else
      list_data => pdata
    end if

! set the pointer to the last block in the list
!
    last_data => pdata

!-------------------------------------------------------------------------------
!
  end subroutine append_datablock
!
!===============================================================================
!
! subroutine REMOVE_DATABLOCK:
! ---------------------------
!
!   Subroutine removes a data block associated with the input pointer from
!   the data block list, and deallocates space used by this block.
!
!   Arguments:
!
!     pdata - the pointer pointing to the data block for removing;
!
!===============================================================================
!
  subroutine remove_datablock(pdata)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! check if the input pointer is associated with a data block
!
    if (associated(pdata)) then

! remove from the meta block list if the meta pointer is set
!
      if (associated(pdata%meta)) then

! if this is the first block in the list, update the list_data pointer
!
        if (pdata%meta%id == list_data%meta%id) list_data => pdata%next

! if this is the last block in the list, update the last_data pointer
!
        if (pdata%meta%id == last_data%meta%id) last_data => pdata%prev

! update the pointer of previous and next blocks
!
        if (associated(pdata%prev)) pdata%prev%next => pdata%next

        if (associated(pdata%next)) pdata%next%prev => pdata%prev

      end if ! %meta associated

! deallocate the associated data block
!
      call deallocate_datablock(pdata)

    end if ! pdata associated with a data block

!-------------------------------------------------------------------------------
!
  end subroutine remove_datablock
!
!===============================================================================
!
! datablock_set_dims: subroutine sets the number of variables and dimensions
!                     for arrays allocated in data blocks
!
!===============================================================================
!
  subroutine datablock_set_dims(nv, nf, ni, nj, nk)

    implicit none

! input arguments
!
    integer(kind=4), intent(in) :: nv, nf, ni, nj, nk
!
!-------------------------------------------------------------------------------
!
    nvars = nv
    nflux = nf
    nx    = ni
    ny    = nj
#if NDIMS == 3
    nz    = nk
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine datablock_set_dims
!
!!==============================================================================
!!
!! REFINEMENT/DEREFINEMENT SUBROUTINES
!!
!
!===============================================================================
!
! refine_block: subroutine refines selected block
!
!===============================================================================
!
  subroutine refine_block(pblock, res, falloc_data)

    use error, only : print_error

    implicit none

! input parameters
!
    type(block_meta), pointer    , intent(inout) :: pblock
    integer(kind=4), dimension(3), intent(in)    :: res
    logical                      , intent(in)    :: falloc_data

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
      call metablock_unset_leaf(pblock)

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
        call metablock_set_leaf(pblock%child(p)%ptr)

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
            pchild => pblock%child(set(i,j,k))%ptr

            if (associated(pneigh)) then
              if (pneigh%id .ne. pblock%id) then

! point to the right neighbor
!
                do p = 1, nfaces
                  pchild%neigh(i,j,p)%ptr => pneigh
                end do

! neighbor level is the same as the refined block
!
                if (pneigh%level .eq. pblock%level) then
                  pneigh%neigh(i,3-j,k)%ptr => pchild
                end if

! neighbor level is the same as the child block
!
                if (pneigh%level .gt. pblock%level) then
                  do p = 1, nfaces
                    pneigh%neigh(i,3-j,p)%ptr => pchild
                  end do
                end if

              end if

            end if

          end do
        end do
      end do

! reset neighbor pointers of the parent block
!
      do i = 1, ndims
        do j = 1, nsides
          do k = 1, nfaces
            nullify(pblock%neigh(i,j,k)%ptr)
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
        ic  = pblock%coord(1) + i * res(1)
        jc  = pblock%coord(2) + j * res(2)
#if NDIMS == 3
        kc  = pblock%coord(3) + k * res(3)
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
        call metablock_set_bounds(pchild, xmn, xmx, ymn, ymx, zmn, zmx)

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
          call link_blocks(pchild, pdata)

        end do

! connect blocks in chain
!
        do p = 2, nchild
          pblock%child(order(p  ))%ptr%data%prev =>                            &
                                             pblock%child(order(p-1))%ptr%data
          pblock%child(order(p-1))%ptr%data%next =>                            &
                                             pblock%child(order(p  ))%ptr%data
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
      call print_error("blocks::refine_block"                                  &
                            , "Input pointer is not associated! Terminating!")

    end if

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
      call metablock_unset_leaf(pblock%child(p)%ptr)
      call deallocate_metablock(pblock%child(p)%ptr)
    end do

! set the leaf flag of parent block
!
    call metablock_set_leaf(pblock)

! reset the refinement flag of the parent block
!
    pblock%refine = 0

!-------------------------------------------------------------------------------
!
  end subroutine derefine_block
#ifdef DEBUG
!!==============================================================================
!!
!! DEBUG SUBROUTINES
!!
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

!-------------------------------------------------------------------------------
!
  end subroutine check_metablock
#endif /* DEBUG */

!===============================================================================
!
end module
