!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2014 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!  This module provides data structures, variables and subroutines to
!!  construct and dynamically modify the hierarchy of blocks corresponding
!!  to the simulated mesh geometry.
!!
!!******************************************************************************
!
module blocks

#ifdef PROFILE
! import external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer, save              :: imi, ima, imu, imp, imq, imr, imd
#endif /* PROFILE */

! MODULE PARAMETERS:
! =================
!
!   ndims     - the number of dimensions (2 or 3);
!   nsides    - the number of sides along each direction (2);
!   nfaces    - the number of faces at each side (2 for 2D, 4 for 3D);
!   nchildren - the number of child blocks for each block (4 for 2D, 8 for 3D);
!
  integer(kind=4), parameter :: ndims     = NDIMS
  integer(kind=4), parameter :: nsides    = 2
  integer(kind=4), parameter :: nfaces    = 2**(ndims - 1)
  integer(kind=4), parameter :: nchildren = 2**ndims

! MODULE VARIABLES:
! ================
!
! the identification of the last allocated block (always increases)
!
  integer(kind=4), save      :: last_id

! the number of allocated meta and data blocks, and the number of leafs
!
  integer(kind=4), save      :: mblocks, dblocks, nleafs

! the number of variables and fluxes stored in data blocks
!
  integer(kind=4), save      :: nvars, nflux

! the spacial dimensions of allocatable data block arrays
!
  integer(kind=4), save      :: nx, ny, nz

! BLOCK STRUCTURE POINTERS:
! ========================
!
! define pointers to meta, data, and info block structures defined below;
! they have to be defined before block structures
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

! BLOCK STRUCTURES:
! ================
!
! define the META block structure; each process keeps exactly same meta block
! structure all the time, so processes can know how the block structure changes
! and where to move data blocks;
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
    type(pointer_meta)          :: child(nchildren)

                                 ! pointers to neighbor meta blocks
                                 !
    type(pointer_meta)          :: neigh(ndims,nsides,nfaces)

                                 ! a pointer to the associated data block
                                 !
    type(block_data)  , pointer :: data

                                 ! the identification number (unique for each
                                 ! block)
                                 !
    integer(kind=4)             :: id

                                 ! the process number to which the meta block
                                 ! is bounded
                                 !
    integer(kind=4)             :: process

                                 ! the level of refinement
                                 !
    integer(kind=4)             :: level

                                 ! the number describing the configuration of
                                 ! the child meta blocks
                                 !
    integer(kind=4)             :: conf

                                 ! the refinement flag, -1, 0, and 1 for
                                 ! the block marked to be derefined, not
                                 ! changed, and refined, respectively
                                 !
    integer(kind=4)             :: refine

                                 ! the position of the block in its siblings
                                 ! group
                                 !
    integer(kind=4)             :: pos(ndims)

                                 ! the coordinate of the lower corner of the
                                 ! block in the effective resolution units
                                 !
    integer(kind=4)             :: coords(ndims)

                                 ! the leaf flag, signifying that the block is
                                 ! the highest block in the local block
                                 ! structure
                                 !
    logical                     :: leaf

                                 ! the flag indicates that the corresponding
                                 ! data needs to be updated (e.g. boundaries or
                                 ! primitive variables), therefore it is
                                 ! usually .true.
                                 !
    logical                     :: update

                                 ! the block bounds in the coordinate units
                                 !
    real(kind=8)                :: xmin, xmax, ymin, ymax, zmin, zmax

  end type block_meta

! define the DATA block structure; all data blocks are divided between
! processes, therefore the same data block cannot be associated with two
! different processes, but they can be moved from one process to another;
!
  type block_data
                                 ! pointers to the previous and next data blocks
                                 !
    type(block_data), pointer :: prev, next

                                 ! a pointer to the associated meta block
                                 !
    type(block_meta), pointer :: meta

                                 ! a pointer to the current conserved variable
                                 ! array
                                 !
    real, dimension(:,:,:,:)  , pointer     :: u

                                 ! an allocatable arrays to store all conserved
                                 ! variables (required two for Runge-Kutta
                                 ! temporal integration methods)
                                 !
    real, dimension(:,:,:,:)  , allocatable :: u0, u1

                                 ! an allocatable array to store all primitive
                                 ! variables
                                 !
    real, dimension(:,:,:,:)  , allocatable :: q

                                 ! an allocatable array to store all fluxes
                                 !
    real, dimension(:,:,:,:,:), allocatable :: f

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

                                 ! the direction, side and face numbers
                                 ! indicating the neighbor block orientation
                                 ! with respect to the block
                                 !
    integer(kind=4)             :: direction, side, face

                                 ! the level difference between the block and
                                 ! its neighbor
                                 !
    integer(kind=4)             :: level_difference

  end type block_info

! POINTERS TO THE FIST AND LAST BLOCKS IN THE LISTS:
! =================================================
!
! these pointers construct the lists of meta and data blocks;
!
  type(block_meta), pointer, save :: list_meta, last_meta
  type(block_data), pointer, save :: list_data, last_data

! all variables and subroutines are private by default
!
  private

! declare public pointers, structures, and variables
!
  public :: pointer_meta, pointer_info
  public :: block_meta, block_data, block_info
  public :: list_meta, list_data
  public :: ndims, nsides, nfaces, nchildren

! declare public subroutines
!
  public :: initialize_blocks, finalize_blocks
  public :: set_block_dimensions
  public :: append_metablock, remove_metablock
  public :: append_datablock, remove_datablock
  public :: allocate_metablock, deallocate_metablock
  public :: allocate_datablock, deallocate_datablock
  public :: link_blocks, unlink_blocks
  public :: refine_block, derefine_block
  public :: set_last_id, get_last_id, get_mblocks, get_dblocks, get_nleafs
  public :: set_blocks_update
  public :: metablock_set_id, metablock_set_process, metablock_set_level
  public :: metablock_set_configuration, metablock_set_refinement
  public :: metablock_set_position, metablock_set_coordinates
  public :: metablock_set_bounds, metablock_set_leaf, metablock_unset_leaf

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!!
!!***  PUBLIC SUBROUTINES  *****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine INITIALIZE_BLOCKS:
! ----------------------------
!
!   Subroutine initializes the module structures, pointers and variables.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_blocks(verbose, iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('blocks:: initialization'         , imi)
    call set_timer('blocks:: meta block allocation'  , ima)
    call set_timer('blocks:: meta block deallocation', imu)
    call set_timer('blocks:: data block allocation'  , imp)
    call set_timer('blocks:: data block deallocation', imq)
    call set_timer('blocks:: refine'                 , imr)
    call set_timer('blocks:: derefine'               , imd)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! reset identification counter
!
    last_id = 0

! reset the number of meta blocks, data blocks, and leafs
!
    mblocks = 0
    dblocks = 0
    nleafs  = 0

! set the initial number of variables and fluxes
!
    nvars   = 1
    nflux   = 1

! set the initial data block resolution
!
    nx      = 1
    ny      = 1
    nz      = 1

! nullify pointers defining the meta and data lists
!
    nullify(list_meta)
    nullify(list_data)
    nullify(last_meta)
    nullify(last_data)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

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
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_blocks(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout)    :: iret

! local variables
!
    type(block_meta), pointer :: pmeta
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! associate the pointer with the last block on the meta block list
!
    pmeta => last_meta

! iterate until the first block on the list is reached
!
    do while(associated(pmeta))

! deallocate the last meta block
!
      call remove_metablock(pmeta)

! assign the pointer to the last block on the meta block list
!
      pmeta => last_meta

    end do ! meta blocks

! nullify pointers defining the meta and data lists
!
    nullify(list_meta)
    nullify(list_data)
    nullify(last_meta)
    nullify(last_data)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_blocks
!
!===============================================================================
!
! subroutine SET_BLOCK_DIMENSIONS:
! -------------------------------
!
!   Subroutine sets the number of variables, fluxes and block dimensions
!   (without ghost cells) for arrays allocated in data blocks.
!
!   Arguments:
!
!     nv - the number of variables stored in %u and %q;
!     nf - the number of fluxes stored in %f;
!     ni - the block dimension along X;
!     nj - the block dimension along Y;
!     nk - the block dimension along Z;
!
!===============================================================================
!
  subroutine set_block_dimensions(nv, nf, ni, nj, nk)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(kind=4), intent(in) :: nv, nf, ni, nj, nk
!
!-------------------------------------------------------------------------------
!
! set the number of variables and fluxes
!
    nvars = nv
    nflux = nf

! set the block dimensions
!
    nx    = ni
    ny    = nj
#if NDIMS == 3
    nz    = nk
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine set_block_dimensions
!
!===============================================================================
!
! subroutine APPEND_METABLOCK:
! ---------------------------
!
!   Subroutine allocates memory for one meta block, appends it to the meta
!   block list and returns a pointer associated with it.
!
!   Arguments:
!
!     pmeta - the pointer associated with the newly appended meta block;
!
!===============================================================================
!
  subroutine append_metablock(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(out) :: pmeta
!
!-------------------------------------------------------------------------------
!
! allocate memory for the new meta block
!
    call allocate_metablock(pmeta)

! check if there are any blocks in the meta block list
!
    if (associated(last_meta)) then

! add the new block to the end of the list
!
      pmeta%prev     => last_meta
      last_meta%next => pmeta

    else

! there are no blocks in the list, so add this one as the first block
!
      list_meta      => pmeta

    end if

! update the pointer to the last block on the list
!
    last_meta => pmeta

!-------------------------------------------------------------------------------
!
  end subroutine append_metablock
!
!===============================================================================
!
! subroutine REMOVE_METABLOCK:
! ---------------------------
!
!   Subroutine removes a meta block associated with the input pointer from
!   the meta block list, and deallocates space used by it.
!
!   Arguments:
!
!     pmeta - the pointer pointing to the meta block which will be removed;
!
!===============================================================================
!
  subroutine remove_metablock(pmeta)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! check if the pointer is actually associated with any block
!
    if (associated(pmeta)) then

! if this is the first block in the list, update the list_meta pointer
!
      if (pmeta%id == list_meta%id) list_meta => pmeta%next

! if this is the last block in the list, update the last_meta pointer
!
      if (pmeta%id == last_meta%id) last_meta => pmeta%prev

! update the %next and %prev pointers of the previous and next blocks,
! respectively
!
      if (associated(pmeta%prev)) pmeta%prev%next => pmeta%next
      if (associated(pmeta%next)) pmeta%next%prev => pmeta%prev

! deallocate memory used by the meta block
!
      call deallocate_metablock(pmeta)

    else

! the argument contains a null pointer, so print an error
!
      call print_error("blocks::remove_metablock"                              &
                                     , "Null pointer argument to meta block!")
    end if

!-------------------------------------------------------------------------------
!
  end subroutine remove_metablock
!
!===============================================================================
!
! subroutine APPEND_DATABLOCK:
! ---------------------------
!
!   Subroutine allocates memory for one data block, appends it to the data
!   block list and returns a pointer associated with it.
!
!   Arguments:
!
!     pdata - the pointer associated with the newly appended data block;
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
! allocate memory for the new data block
!
    call allocate_datablock(pdata)

! check if there are any blocks in the data block list
!
    if (associated(last_data)) then

! add the new block to the end of the list
!
      pdata%prev     => last_data
      last_data%next => pdata

    else

! there are no blocks in the list, so add this one as the first block
!
      list_data      => pdata

    end if

! update the pointer to the last block on the list
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
!   the data block list, and deallocates space used by it.
!
!   Arguments:
!
!     pdata - the pointer pointing to the data block which will be removed;
!
!===============================================================================
!
  subroutine remove_datablock(pdata)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! check if the pointer is actually associated with any block
!
    if (associated(pdata)) then

! check if the data block has associated meta block
!
      if (associated(pdata%meta)) then

! if this is the first block in the list, update the list_data pointer
!
        if (pdata%meta%id == list_data%meta%id) list_data => pdata%next

! if this is the last block in the list, update the last_data pointer
!
        if (pdata%meta%id == last_data%meta%id) last_data => pdata%prev

! update the %next and %prev pointers of the previous and next blocks,
! respectively
!
        if (associated(pdata%prev)) pdata%prev%next => pdata%next
        if (associated(pdata%next)) pdata%next%prev => pdata%prev

      else ! %meta associated

! there is no meta block associated, so print an error
!
        call print_error("blocks::remove_datablock"                            &
                            , "No meta block associated with the data block!")

      end if ! %meta associated

! deallocate the associated data block
!
      call deallocate_datablock(pdata)

    else

! the argument contains a null pointer, so print an error
!
      call print_error("blocks::remove_datablock"                              &
                                     , "Null pointer argument to data block!")
    end if

!-------------------------------------------------------------------------------
!
  end subroutine remove_datablock
!
!===============================================================================
!
! subroutine ALLOCATE_METABLOCK:
! -----------------------------
!
!   Subroutine allocates memory for one meta block, initializes its fields
!   and returns a pointer associated with it.
!
!   Arguments:
!
!     pmeta - the pointer associated with the newly allocated meta block;
!
!===============================================================================
!
  subroutine allocate_metablock(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(out) :: pmeta

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the meta block allocation
!
    call start_timer(ima)
#endif /* PROFILE */

! allocate the meta block structure for one object
!
    allocate(pmeta)

! nullify fields pointing to previous and next block on the meta block list
!
    nullify(pmeta%prev)
    nullify(pmeta%next)

! nullify the field pointing to the parent
!
    nullify(pmeta%parent)

! nullify fields pointing to children
!
    do i = 1, nchildren
      nullify(pmeta%child(i)%ptr)
    end do

! nullify fields pointing to neighbors
!
    do i = 1, ndims
      do j = 1, nsides
        do k = 1, nfaces
          nullify(pmeta%neigh(i,j,k)%ptr)
        end do
      end do
    end do

! nullify the field pointing to the associated data block
!
    nullify(pmeta%data)

! set unique ID
!
    pmeta%id        = increase_id()

! unset the process number, level, the children configuration, refine, leaf,
! and update flags
!
    pmeta%process   = -1
    pmeta%level     = -1
    pmeta%conf      = -1
    pmeta%refine    =  0
    pmeta%leaf      = .false.
    pmeta%update    = .true.

! initialize the position in the parent block
!
    pmeta%pos(:)    = -1

! initialize the effective coordinates
!
    pmeta%coords(:) = 0

! initialize coordinate bounds of the block
!
    pmeta%xmin      = 0.0d+00
    pmeta%xmax      = 1.0d+00
    pmeta%ymin      = 0.0d+00
    pmeta%ymax      = 1.0d+00
    pmeta%zmin      = 0.0d+00
    pmeta%zmax      = 1.0d+00

! increase the number of allocated meta blocks
!
    mblocks         = mblocks + 1

#ifdef PROFILE
! stop accounting time for the meta block allocation
!
    call stop_timer(ima)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine allocate_metablock
!
!===============================================================================
!
! subroutine DEALLOCATE_METABLOCK:
! -------------------------------
!
!   Subroutine releases memory used by the meta block associated with
!   the pointer argument.
!
!   Arguments:
!
!     pmeta - the pointer associated with the meta block which will be
!             deallocated;
!
!===============================================================================
!
  subroutine deallocate_metablock(pmeta)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the meta block deallocation
!
    call start_timer(imu)
#endif /* PROFILE */

! check if the pointer is actually associated with any block
!
    if (associated(pmeta)) then

! decrease the number of leafs
!
      if (pmeta%leaf) nleafs = nleafs - 1

! decrease the number of allocated meta blocks
!
      mblocks = mblocks - 1

! nullify fields pointing to previous and next block on the meta block list
!
      nullify(pmeta%prev)
      nullify(pmeta%next)

! nullify the field pointing to the parent
!
      nullify(pmeta%parent)

! nullify fields pointing to children
!
      do i = 1, nchildren
        nullify(pmeta%child(i)%ptr)
      end do

! nullify fields pointing to neighbors
!
      do i = 1, ndims
        do j = 1, nsides
          do k = 1, nfaces
            nullify(pmeta%neigh(i,j,k)%ptr)
          end do
        end do
      end do

! if there is a data block is associated, remove it
!
      if (associated(pmeta%data)) call remove_datablock(pmeta%data)

! nullify the field pointing to the associated data block
!
      nullify(pmeta%data)

! release the memory occupied by the block
!
      deallocate(pmeta)

! nullify the pointer to the deallocated meta block
!
      nullify(pmeta)

    else

! the argument contains a null pointer, so print an error
!
      call print_error("blocks::deallocate_metablock"                          &
                                     , "Null pointer argument to meta block!")
    end if

#ifdef PROFILE
! stop accounting time for the meta block deallocation
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_metablock
!
!===============================================================================
!
! subroutine ALLOCATE_DATABLOCK:
! -----------------------------
!
!   Subroutine allocates memory for one data block, initializes its fields
!   and returns a pointer associated with it.
!
!   Arguments:
!
!     pdata - the pointer associated with the newly allocated data block;
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
#ifdef PROFILE
! start accounting time for the data block allocation
!
    call start_timer(imp)
#endif /* PROFILE */

! allocate the block structure
!
    allocate(pdata)

! nullify field pointing to the previous and next blocks on the data block list
!
    nullify(pdata%prev)
    nullify(pdata%next)

! nullify the field pointing to the associate meta block list
!
    nullify(pdata%meta)

! allocate space for conserved variables
!
    allocate(pdata%u0(nvars,nx,ny,nz), pdata%u1(nvars,nx,ny,nz))

! allocate space for primitive variables
!
    allocate(pdata%q(nvars,nx,ny,nz))

! allocate space for numerical fluxes
!
    allocate(pdata%f(ndims,nflux,nx,ny,nz))

! initiate the conserved variable pointer
!
    pdata%u => pdata%u0

! increase the number of allocated data blocks
!
    dblocks = dblocks + 1

#ifdef PROFILE
! stop accounting time for the data block allocation
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine allocate_datablock
!
!===============================================================================
!
! subroutine DEALLOCATE_DATABLOCK:
! -------------------------------
!
!   Subroutine releases memory used by the data block associated with
!   the pointer argument.
!
!   Arguments:
!
!     pdata - the pointer associated with the data block which will be
!             deallocated;
!
!===============================================================================
!
  subroutine deallocate_datablock(pdata)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the data block deallocation
!
    call start_timer(imq)
#endif /* PROFILE */

! check if the pointer is actually associated with any block
!
    if (associated(pdata)) then

! decrease the number of allocated data blocks
!
      dblocks = dblocks - 1

! nullify field pointing to the previous and next blocks on the data block list
!
      nullify(pdata%prev)
      nullify(pdata%next)

! nullify the field pointing to the associate meta block list
!
      nullify(pdata%meta)

! nullify pointer to the current conserved variable array
!
      nullify(pdata%u)

! deallocate conserved variables
!
      if (allocated(pdata%u0)) deallocate(pdata%u0)
      if (allocated(pdata%u1)) deallocate(pdata%u1)

! deallocate primitive variables
!
      if (allocated(pdata%q )) deallocate(pdata%q )

! deallocate numerical fluxes
!
      if (allocated(pdata%f )) deallocate(pdata%f )

! release the memory occupied by the block
!
      deallocate(pdata)

! nullify the pointer to the deallocated meta block
!
      nullify(pdata)

    else

! the argument contains a null pointer, so print an error
!
      call print_error("blocks::deallocate_datablock"                          &
                                     , "Null pointer argument to data block!")

    end if ! pdata associated with a data block

#ifdef PROFILE
! stop accounting time for the data block deallocation
!
    call stop_timer(imq)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_datablock
!
!===============================================================================
!
! subroutine LINK_BLOCKS:
! ----------------------
!
!   Subroutine links meta and data blocks.
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
! associate the corresponging pointers
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
!   Subroutine unlinks meta and data blocks.
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
! nullify the corresponging pointers
!
    nullify(pmeta%data)
    nullify(pdata%meta)

!-------------------------------------------------------------------------------
!
  end subroutine unlink_blocks
!
!===============================================================================
!
! subroutine REFINE_BLOCK:
! -----------------------
!
!   Subroutine creates children of the current block and initializes their
!   configuration, pointers and fields.
!
!   Arguments:
!
!     pmeta - a pointer to meta block for which children will be created;
!     res   - the resolution of the block;
!     fdata - a flag indicating if data blocks for children should be allocated;
!
!===============================================================================
!
  subroutine refine_block(pmeta, res, fdata)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer    , intent(inout) :: pmeta
    integer(kind=4), dimension(3), intent(in)    :: res
    logical                      , intent(in)    :: fdata

! pointers
!
    type(block_meta), pointer :: pnext, pneigh, pchild
    type(block_data), pointer :: pdata

! local variables
!
    integer :: p, q, i, j, k, ic, jc, kc
    real    :: xln, yln, zln, xmn, xmx, ymn, ymx, zmn, zmx

! local arrays
!
    integer, dimension(nchildren)           :: config, order
    integer, dimension(ndims,nsides,nfaces) :: set
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the block refinement
!
    call start_timer(imr)
#endif /* PROFILE */

! check if pointer is associated
!
    if (associated(pmeta)) then

! store the pointer to the next block on the list
!
      pnext => pmeta%next

!! PREPARE CHILD CONFIGURATION PARAMETERS
!!
! set corresponding configuration of the new blocks
!
      select case(pmeta%conf)
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

! calculate sizes of the child blocks
!
      xln = 0.5d+00 * (pmeta%xmax - pmeta%xmin)
      yln = 0.5d+00 * (pmeta%ymax - pmeta%ymin)
#if NDIMS == 3
      zln = 0.5d+00 * (pmeta%zmax - pmeta%zmin)
#else /* NDIMS == 3 */
      zln =           (pmeta%zmax - pmeta%zmin)
#endif /* NDIMS == 3 */

!! ALLOCATE CHILDREN AND APPEND THEM TO THE META LIST
!!
! iterate over the number of children in the reverse configuration order, i.e.
! the allocated blocks are inserted after the parent block following
! the reversed Hilbert curve
!
      do p = nchildren, 1, -1

! insert a new meta block after pmeta and associate it with pchild
!
        call insert_metablock_after(pmeta, pchild)

! set the child configuration number
!
        call metablock_set_configuration(pchild, config(p))

! associate the parent's children array element with the freshly created
! meta block
!
        pmeta%child(order(p))%ptr => pchild

      end do ! nchildren

! iterate over all children
!
      do p = 1, nchildren

! associate a pointer with the current child
!
        pchild        => pmeta%child(p)%ptr

! associate the parent field with pmeta
!
        pchild%parent => pmeta

! mark the child as the leaf
!
        call metablock_set_leaf(pchild)

! set the child refinement level
!
        call metablock_set_level(pchild, pmeta%level + 1)

! set the child process number
!
        call metablock_set_process(pchild, pmeta%process)

! calculate the block position indices
!
        q   = p - 1
        i   = mod(q    ,2)
        j   = mod(q / 2,2)
        k   = mod(q / 4,2)

! calculate the block coordinates in effective resolution units
!
        ic  = pmeta%coords(1) + i * res(1)
        jc  = pmeta%coords(2) + j * res(2)
#if NDIMS == 3
        kc  = pmeta%coords(3) + k * res(3)
#endif /* NDIMS == 3 */

! calculate block bounds
!
        xmn = pmeta%xmin + xln * i
        ymn = pmeta%ymin + yln * j
        zmn = pmeta%zmin + zln * k

        xmx = xmn        + xln
        ymx = ymn        + yln
        zmx = zmn        + zln

! set the block position
!
        call metablock_set_position(pchild, i, j, k)

! set the effective resolution coordinates
!
        call metablock_set_coordinates(pchild, ic, jc, kc)

! set the child block bounds
!
        call metablock_set_bounds(pchild, xmn, xmx, ymn, ymx, zmn, zmx)

      end do ! nchildren

!! ASSIGN PROPER NEIGHBORS FOR THE CHILDREN IN THE INTERIOR OF THE PARENT BLOCK
!!
! iterate over faces and update the interior of the block
!
      do p = 1, nfaces

! X direction (left side)
!
        pmeta%child(2)%ptr%neigh(1,1,p)%ptr => pmeta%child(1)%ptr
        pmeta%child(4)%ptr%neigh(1,1,p)%ptr => pmeta%child(3)%ptr
#if NDIMS == 3
        pmeta%child(6)%ptr%neigh(1,1,p)%ptr => pmeta%child(5)%ptr
        pmeta%child(8)%ptr%neigh(1,1,p)%ptr => pmeta%child(7)%ptr
#endif /* NDIMS == 3 */

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(1,1,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(1)%ptr%neigh(1,1,p)%ptr => pmeta%child(2)%ptr
            pmeta%child(3)%ptr%neigh(1,1,p)%ptr => pmeta%child(4)%ptr
#if NDIMS == 3
            pmeta%child(5)%ptr%neigh(1,1,p)%ptr => pmeta%child(6)%ptr
            pmeta%child(7)%ptr%neigh(1,1,p)%ptr => pmeta%child(8)%ptr
#endif /* NDIMS == 3 */

          end if

        end if

! X direction (right side)
!
        pmeta%child(1)%ptr%neigh(1,2,p)%ptr => pmeta%child(2)%ptr
        pmeta%child(3)%ptr%neigh(1,2,p)%ptr => pmeta%child(4)%ptr
#if NDIMS == 3
        pmeta%child(5)%ptr%neigh(1,2,p)%ptr => pmeta%child(6)%ptr
        pmeta%child(7)%ptr%neigh(1,2,p)%ptr => pmeta%child(8)%ptr
#endif /* NDIMS == 3 */

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(1,2,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(2)%ptr%neigh(1,2,p)%ptr => pmeta%child(1)%ptr
            pmeta%child(4)%ptr%neigh(1,2,p)%ptr => pmeta%child(3)%ptr
#if NDIMS == 3
            pmeta%child(6)%ptr%neigh(1,2,p)%ptr => pmeta%child(5)%ptr
            pmeta%child(8)%ptr%neigh(1,2,p)%ptr => pmeta%child(7)%ptr
#endif /* NDIMS == 3 */

          end if

        end if

! Y direction (left side)
!
        pmeta%child(3)%ptr%neigh(2,1,p)%ptr => pmeta%child(1)%ptr
        pmeta%child(4)%ptr%neigh(2,1,p)%ptr => pmeta%child(2)%ptr
#if NDIMS == 3
        pmeta%child(7)%ptr%neigh(2,1,p)%ptr => pmeta%child(5)%ptr
        pmeta%child(8)%ptr%neigh(2,1,p)%ptr => pmeta%child(6)%ptr
#endif /* NDIMS == 3 */

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(2,1,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(1)%ptr%neigh(2,1,p)%ptr => pmeta%child(3)%ptr
            pmeta%child(2)%ptr%neigh(2,1,p)%ptr => pmeta%child(4)%ptr
#if NDIMS == 3
            pmeta%child(5)%ptr%neigh(2,1,p)%ptr => pmeta%child(7)%ptr
            pmeta%child(6)%ptr%neigh(2,1,p)%ptr => pmeta%child(8)%ptr
#endif /* NDIMS == 3 */

          end if

        end if

! Y direction (right side)
!
        pmeta%child(1)%ptr%neigh(2,2,p)%ptr => pmeta%child(3)%ptr
        pmeta%child(2)%ptr%neigh(2,2,p)%ptr => pmeta%child(4)%ptr
#if NDIMS == 3
        pmeta%child(5)%ptr%neigh(2,2,p)%ptr => pmeta%child(7)%ptr
        pmeta%child(6)%ptr%neigh(2,2,p)%ptr => pmeta%child(8)%ptr
#endif /* NDIMS == 3 */

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(2,2,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(3)%ptr%neigh(2,2,p)%ptr => pmeta%child(1)%ptr
            pmeta%child(4)%ptr%neigh(2,2,p)%ptr => pmeta%child(2)%ptr
#if NDIMS == 3
            pmeta%child(7)%ptr%neigh(2,2,p)%ptr => pmeta%child(5)%ptr
            pmeta%child(8)%ptr%neigh(2,2,p)%ptr => pmeta%child(6)%ptr
#endif /* NDIMS == 3 */

          end if

        end if

#if NDIMS == 3
! Z direction (left side)
!
        pmeta%child(5)%ptr%neigh(3,1,p)%ptr => pmeta%child(1)%ptr
        pmeta%child(6)%ptr%neigh(3,1,p)%ptr => pmeta%child(2)%ptr
        pmeta%child(7)%ptr%neigh(3,1,p)%ptr => pmeta%child(3)%ptr
        pmeta%child(8)%ptr%neigh(3,1,p)%ptr => pmeta%child(4)%ptr

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(3,1,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(1)%ptr%neigh(3,1,p)%ptr => pmeta%child(5)%ptr
            pmeta%child(2)%ptr%neigh(3,1,p)%ptr => pmeta%child(6)%ptr
            pmeta%child(3)%ptr%neigh(3,1,p)%ptr => pmeta%child(7)%ptr
            pmeta%child(4)%ptr%neigh(3,1,p)%ptr => pmeta%child(8)%ptr

          end if

        end if

! Z direction (right side)
!
        pmeta%child(1)%ptr%neigh(3,2,p)%ptr => pmeta%child(5)%ptr
        pmeta%child(2)%ptr%neigh(3,2,p)%ptr => pmeta%child(6)%ptr
        pmeta%child(3)%ptr%neigh(3,2,p)%ptr => pmeta%child(7)%ptr
        pmeta%child(4)%ptr%neigh(3,2,p)%ptr => pmeta%child(8)%ptr

! associate pneigh with a neighbor
!
        pneigh => pmeta%neigh(3,2,1)%ptr

! if neighbor and associated and points to parent block, it corresponds to
! periodic boundaries at the lowest level
!
        if (associated(pneigh)) then

          if (pneigh%id == pmeta%id) then

            pmeta%child(5)%ptr%neigh(3,2,p)%ptr => pmeta%child(1)%ptr
            pmeta%child(6)%ptr%neigh(3,2,p)%ptr => pmeta%child(2)%ptr
            pmeta%child(7)%ptr%neigh(3,2,p)%ptr => pmeta%child(3)%ptr
            pmeta%child(8)%ptr%neigh(3,2,p)%ptr => pmeta%child(4)%ptr

          end if

        end if
#endif /* NDIMS == 3 */

      end do ! nfaces

!! UPDATE NEIGHBORS AND EXTERNAL NEIGHBORS OF CHILDREN
!!
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

! set pointers to neighbors and update neighbors' pointers
!
      do i = 1, ndims
        do j = 1, nsides

! prepare reverse side index
!
          q = 3 - j

! iterate over all faces
!
          do k = 1, nfaces

! associate pointers with the neighbor and child
!
            pneigh => pmeta%neigh(i,j,k)%ptr
            pchild => pmeta%child(set(i,j,k))%ptr

! check if neighbor is associated
!
            if (associated(pneigh)) then

! check if the parent block does not point to itself (periodic boundaries)
!
              if (pneigh%id /= pmeta%id) then

! point the child neigh field to the right neighbor
!
                do p = 1, nfaces
                  pchild%neigh(i,j,p)%ptr => pneigh
                end do

! update neighbor pointer if it is at the same level
!
                if (pneigh%level == pmeta%level) then
                  pneigh%neigh(i,q,k)%ptr => pchild
                end if

! update neighbor pointer if it is at higher level
!
                if (pneigh%level > pmeta%level) then
                  do p = 1, nfaces
                    pneigh%neigh(i,q,p)%ptr => pchild
                  end do
                end if

! if neighbor has lower level than parent, something is wrong, since lower
! levels should be already refined
!
                if (pneigh%level < pmeta%level) then
                  call print_error("blocks::refine_block"                      &
                                           , "Neighbor found at lower level!")
                end if

              end if ! pmeta and pneigh point to different blocks

            end if ! pneigh is associated

          end do ! nfaces
        end do ! nsides
      end do ! ndims

!! ASSOCIATE DATA BLOCKS IF NECESSARY
!!
! allocate data blocks if requested
!
      if (fdata) then

! iterate over all children
!
        do p = 1, nchildren

! assign a pointer to the current child
!
          pchild => pmeta%child(p)%ptr

! allocate new data block and append it to the data block list
!
          call append_datablock(pdata)

! associate the new data block with the current child
!
          call link_blocks(pchild, pdata)

        end do ! nchildren

      end if ! allocate data blocks for children

!! SET SURROUNDING NEIGHBORS TO BE UPDATED TOO
!!

!! RESET PARENT'S FIELDS
!!
! unset the block leaf flag
!
      call metablock_unset_leaf(pmeta)

! reset the refinement flag
!
      call metablock_set_refinement(pmeta, 0)

! nullify the parent's neighbor pointers
!
      do i = 1, ndims
        do j = 1, nsides
          do k = 1, nfaces
            nullify(pmeta%neigh(i,j,k)%ptr)
          end do
        end do
      end do

! restore the pointer to the current block
!
      if (associated(pnext)) then
        pmeta => pnext%prev
      else
        pmeta => last_meta
      end if

    else ! pmeta is not associated

! it's impossible to refine since there is not block associated with
! the argument pointer
!
      call print_error("blocks::refine_block"                                  &
                           , "No block associated with the argument pointer!")

    end if

#ifdef PROFILE
! stop accounting time for the block refinement
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine refine_block
!
!===============================================================================
!
! subroutine DEREFINE_BLOCK:
! -------------------------
!
!   Subroutine derefines the current block by distrying all its children and
!   restoring the block configuration, pointers and fields.
!
!   Arguments:
!
!     pmeta - a pointer to derefined meta block;
!
!===============================================================================
!
  subroutine derefine_block(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta

! local pointers
!
    type(block_meta), pointer :: pchild, pneigh

! local variables
!
    integer       :: i, j, k, l, p

! local saved variables
!
    logical, save :: first = .true.

! local arrays
!
    integer, dimension(ndims, nsides, nfaces), save :: arr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the block derefinement
!
    call start_timer(imd)
#endif /* PROFILE */

! prepare saved variables at the first execution
!
    if (first) then

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

! reset the first execution flag
!
      first = .false.

    end if

! iterate over dimensions, sides, and faces
!
    do i = 1, ndims
      do j = 1, nsides
        do k = 1, nfaces

! get the current child index
!
          p = arr(i,j,k)

! associate a pointer with the neighbor
!
          pneigh => pmeta%child(p)%ptr%neigh(i,j,k)%ptr

! update the parent neighbor field
!
          pmeta%neigh(i,j,k)%ptr => pneigh

! update the neigh field of the neighbor
!
          if (associated(pneigh)) then

            l = 3 - j
            do p = 1, nfaces
              pneigh%neigh(i,l,p)%ptr => pmeta
            end do

          end if ! pneigh is associated

        end do ! nfaces
      end do ! nsides
    end do ! ndims

! iterate over children
!
    do p = 1, nchildren

! remove the child from the meta block list
!
      call remove_metablock(pmeta%child(p)%ptr)

    end do ! nchild

! update the parent leaf flag
!
    call metablock_set_leaf(pmeta)

! reset the refinement flag of the parent block
!
    call metablock_set_refinement(pmeta, 0)

#ifdef PROFILE
! stop accounting time for the block derefinement
!
    call stop_timer(imd)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine derefine_block
!
!===============================================================================
!
! subroutine SET_LAST_ID:
! ----------------------
!
!   Subroutine sets the last identification number.  This subroutine should
!   be only used when the job is resumed.
!
!   Arguments:
!
!     id - the identification number to set;
!
!===============================================================================
!
  subroutine set_last_id(id)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(kind=4), intent(in) :: id
!
!-------------------------------------------------------------------------------
!
! check if the new last_id is larger than the already existing
!
    if (last_id > id) then

      call print_error("blocks::set_last_id"                                   &
                             , "New last_id must be larger than the old one!")
    else

! set the last identification number
!
      last_id = id

    end if

!-------------------------------------------------------------------------------
!
  end subroutine set_last_id
!
!===============================================================================
!
! function GET_LAST_ID:
! --------------------
!
!   Function returns the last identification number.
!
!
!===============================================================================
!
  function get_last_id() result(id)

! local variables are not implicit by default
!
    implicit none

! return variable
!
    integer(kind=4) :: id
!
!-------------------------------------------------------------------------------
!
! set the return value
!
    id = last_id

!-------------------------------------------------------------------------------
!
  end function get_last_id
!
!===============================================================================
!
! function GET_MBLOCKS:
! --------------------
!
!   Function returns the number of meta blocks.
!
!
!===============================================================================
!
  function get_mblocks() result(nr)

! local variables are not implicit by default
!
    implicit none

! return variable
!
    integer(kind=4) :: nr
!
!-------------------------------------------------------------------------------
!
! set the return value
!
    nr = mblocks

!-------------------------------------------------------------------------------
!
  end function get_mblocks
!
!===============================================================================
!
! function GET_DBLOCKS:
! --------------------
!
!   Function returns the number of data blocks.
!
!
!===============================================================================
!
  function get_dblocks() result(nr)

! local variables are not implicit by default
!
    implicit none

! return variable
!
    integer(kind=4) :: nr
!
!-------------------------------------------------------------------------------
!
! set the return value
!
    nr = dblocks

!-------------------------------------------------------------------------------
!
  end function get_dblocks
!
!===============================================================================
!
! function GET_NLEAFS:
! -------------------
!
!   Function returns the number of leafs.
!
!
!===============================================================================
!
  function get_nleafs() result(nr)

! local variables are not implicit by default
!
    implicit none

! return variable
!
    integer(kind=4) :: nr
!
!-------------------------------------------------------------------------------
!
! set the return value
!
    nr = nleafs

!-------------------------------------------------------------------------------
!
  end function get_nleafs
!
!===============================================================================
!
! subroutine SET_BLOCKS_UPDATE:
! ----------------------------
!
!   Subroutine sets the update flag of all meta block in the list.
!
!   Arguments:
!
!     flag - the flag to be set;
!
!===============================================================================
!
  subroutine set_blocks_update(flag)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)       :: flag

! local pointers
!
    type(block_meta), pointer :: pmeta
!
!-------------------------------------------------------------------------------
!
! associate the pointer with the first block on the meta block list
!
    pmeta => list_meta

! iterate over all blocks in the list
!
    do while(associated(pmeta))

! mark the block for update
!
      pmeta%update = flag

! associate the pointer with the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

!-------------------------------------------------------------------------------
!
  end subroutine set_blocks_update
!
!===============================================================================
!
! subroutine METABLOCK_SET_ID:
! ---------------------------
!
!   Subroutine sets the identification number of the meta block pointed by
!   the input argument. This subroutine should be used only when resuming jobs.
!
!   Arguments:
!
!     pmeta - a pointer to the updated meta block;
!     id    - the identification number to set;
!
!===============================================================================
!
  subroutine metablock_set_id(pmeta, id)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: id
!
!-------------------------------------------------------------------------------
!
! set the meta block %id field
!
    pmeta%id = id

! check if the last identification number is smaller than id, if so set
! the value of last_id to id
!
    if (last_id < id) last_id = id

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_id
!
!===============================================================================
!
! subroutine METABLOCK_SET_PROCESS:
! --------------------------------
!
!   Subroutine sets the process number of the meta block pointed by
!   the input argument.
!
!   Arguments:
!
!     pmeta - a pointer to the updated meta block;
!     np    - the process number;
!
!===============================================================================
!
  subroutine metablock_set_process(pmeta, np)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: np
!
!-------------------------------------------------------------------------------
!
! set the block's %process field
!
    pmeta%process = np

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_process
!
!===============================================================================
!
! subroutine METABLOCK_SET_LEVEL:
! ------------------------------
!
!   Subroutine sets the refinement level number of the meta block pointed
!   by the input argument.
!
!   Arguments:
!
!     pmeta - a pointer to the updated meta block;
!     lv    - the refinement level number;
!
!===============================================================================
!
  subroutine metablock_set_level(pmeta, lv)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: lv
!
!-------------------------------------------------------------------------------
!
! set the block's refinement level
!
    pmeta%level = lv

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_level
!
!===============================================================================
!
! subroutine METABLOCK_SET_CONFIGURATION:
! --------------------------------------
!
!   Subroutine sets the children block configuration number of the meta block
!   pointed by the input argument.
!
!   Arguments:
!
!     pmeta - a pointer to the updated meta block;
!     cf    - the configuration number;
!
!===============================================================================
!
  subroutine metablock_set_configuration(pmeta, cf)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: cf
!
!-------------------------------------------------------------------------------
!
! set the block's children configuration number
!
    pmeta%conf = cf

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_configuration
!
!===============================================================================
!
! subroutine METABLOCK_SET_REFINEMENT:
! -----------------------------------
!
!   Subroutine sets the refinement flag of the meta block pointed by
!   the input argument.
!
!   Arguments:
!
!     pmeta - a pointer to the updated meta block;
!     rf    - the refinement flag;
!
!===============================================================================
!
  subroutine metablock_set_refinement(pmeta, rf)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: rf
!
!-------------------------------------------------------------------------------
!
! check if the refinement value is correct
!
    if (abs(rf) > 1) then

! print error about wrong refine flag
!
      call print_error("blocks::metablock_set_refinement"                      &
                                           , "The refinement value is wrong!")

    else

! set the block's refinement field
!
      pmeta%refine = rf

    end if

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_refinement
!
!===============================================================================
!
! subroutine METABLOCK_SET_POSITION:
! ---------------------------------
!
!   Subroutine sets the position coordinates in the parent block of
!   the meta block pointed by the input argument.
!
!   Arguments:
!
!     pmeta      - a pointer to the updated meta block;
!     px, py, pz - the block position coordinates;
!
!===============================================================================
!
  subroutine metablock_set_position(pmeta, px, py, pz)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: px, py, pz
!
!-------------------------------------------------------------------------------
!
! set the block's position in the parent block
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
! subroutine METABLOCK_SET_COORDINATES:
! ------------------------------------
!
!   Subroutine sets the effective resolution coordinates of the meta block
!   pointed by the input argument.
!
!   Arguments:
!
!     pmeta      - a pointer to the updated meta block;
!     px, py, pz - the effective resolution coordinates;
!
!===============================================================================
!
  subroutine metablock_set_coordinates(pmeta, px, py, pz)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    integer(kind=4)          , intent(in)    :: px, py, pz
!
!-------------------------------------------------------------------------------
!
! set the block's effective resolution coordinates
!
    pmeta%coords(1) = px
    pmeta%coords(2) = py
#if NDIMS == 3
    pmeta%coords(3) = pz
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_coordinates
!
!===============================================================================
!
! subroutine METABLOCK_SET_BOUNDS:
! -------------------------------
!
!   Subroutine sets the physical coordinate bounds of the meta block pointed
!   by the input argument.
!
!   Arguments:
!
!     pmeta    - a pointer to the updated meta block;
!     xmn, xmx - the coordinate bounds along X;
!     ymn, ymx - the coordinate bounds along Y;
!     zmn, zmx - the coordinate bounds along Z;
!
!===============================================================================
!
  subroutine metablock_set_bounds(pmeta, xmn, xmx, ymn, ymx, zmn, zmx)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
    real(kind=8)             , intent(in)    :: xmn, xmx
    real(kind=8)             , intent(in)    :: ymn, ymx
    real(kind=8)             , intent(in)    :: zmn, zmx
!
!-------------------------------------------------------------------------------
!
! set the block's coordinate bounds
!
    pmeta%xmin = xmn
    pmeta%xmax = xmx
    pmeta%ymin = ymn
    pmeta%ymax = ymx
    pmeta%zmin = zmn
    pmeta%zmax = zmx

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_bounds
!
!===============================================================================
!
! subroutine METABLOCK_SET_LEAF:
! -----------------------------
!
!   Subroutine marks the meta block pointed by the input argument as the leaf.
!
!   Arguments:
!
!     pmeta    - a pointer to the updated meta block;
!
!===============================================================================
!
  subroutine metablock_set_leaf(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! set the block's leaf flag
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
! subroutine METABLOCK_UNSET_LEAF:
! -------------------------------
!
!   Subroutine marks the meta block pointed by the input argument as non-leaf.
!
!   Arguments:
!
!     pmeta    - a pointer to the updated meta block;
!
!===============================================================================
!
  subroutine metablock_unset_leaf(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! unset the block's leaf flag
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
! subroutine METABLOCK_SET_UPDATE:
! -------------------------------
!
!   Subroutine marks the meta block pointed by the input argument to be updated.
!
!   Arguments:
!
!     pmeta    - a pointer to the updated meta block;
!
!===============================================================================
!
  subroutine metablock_set_update(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! set the block's update flag
!
    pmeta%update = .true.

!-------------------------------------------------------------------------------
!
  end subroutine metablock_set_update
!
!===============================================================================
!
! subroutine METABLOCK_UNSET_UPDATE:
! ---------------------------------
!
!   Subroutine marks the meta block pointed by the input argument to not
!   be updated.
!
!   Arguments:
!
!     pmeta    - a pointer to the updated meta block;
!
!===============================================================================
!
  subroutine metablock_unset_update(pmeta)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta
!
!-------------------------------------------------------------------------------
!
! unset the block's update flag
!
    pmeta%update = .false.

!-------------------------------------------------------------------------------
!
  end subroutine metablock_unset_update
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine INSERT_METABLOCK_AFTER:
! ---------------------------------
!
!   Subroutine allocates memory for one meta block, inserts it to the meta
!   block list after the provided pointer and returns a pointer associated
!   with it.
!
!   Arguments:
!
!     pmeta - the pointer associated with the newly appended meta block;
!     pprev - the pointer after which the new block has to be inserted;
!
!===============================================================================
!
  subroutine insert_metablock_after(pprev, pmeta)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(in)  :: pprev
    type(block_meta), pointer, intent(out) :: pmeta
!
!-------------------------------------------------------------------------------
!
! allocate memory for the new meta block
!
    call allocate_metablock(pmeta)

! if pprev is associated, insert the new block after it
!
    if (associated(pprev)) then

! associate the %prev and %next pointers
!
      pmeta%prev => pprev
      pmeta%next => pprev%next

! update the pointer of the next and previous blocks
!
      if (associated(pprev%next)) pprev%next%prev => pmeta
      pprev%next => pmeta

! check if last_meta is associated
!
      if (associated(last_meta)) then

! update the last_meta pointer if necessary
!
        if (pprev%id == last_meta%id) last_meta => pmeta

      else

! strange situation, pprev is associated, but last_meta not
!
        call print_error("blocks::intert_metablock_after"                      &
                       , "Argument pprev is associated but last_meta is not!")

      end if

    else

! if pprev is null and list_meta is associated, there is something wrong
!
      if (associated(list_meta)) then

! strange situation, pprev is null but list_meta is associated
!
        call print_error("blocks::intert_metablock_after"                      &
                      , "Argument pprev is null but list_meta is associated!")

      else

! pprev and list_meta are nulls, so add the first block to the list by
! associating list_meta and last_meta
!
        list_meta => pmeta
        last_meta => pmeta

      end if

    end if

!-------------------------------------------------------------------------------
!
  end subroutine insert_metablock_after
!
!===============================================================================
!
! subroutine INSERT_METABLOCK_BEFORE:
! ----------------------------------
!
!   Subroutine allocates memory for one meta block, inserts it to the meta
!   block list before the provided pointer and returns a pointer associated
!   with it.
!
!   Arguments:
!
!     pmeta - the pointer associated with the newly appended meta block;
!     pnext - the pointer before which the new block has to be inserted;
!
!===============================================================================
!
  subroutine insert_metablock_before(pnext, pmeta)

! import external procedures
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(in)  :: pnext
    type(block_meta), pointer, intent(out) :: pmeta
!
!-------------------------------------------------------------------------------
!
! allocate memory for the new meta block
!
    call allocate_metablock(pmeta)

! if pnext is associated, insert the new block before it
!
    if (associated(pnext)) then

! associate the %prev and %next pointers
!
      pmeta%prev => pnext%prev
      pmeta%next => pnext

! update the pointer of the next and previous blocks
!
      if (associated(pnext%prev)) pnext%prev%next => pmeta
      pnext%prev => pmeta

! check if list_meta is associated
!
      if (associated(list_meta)) then

! update the list_meta pointer if necessary
!
        if (pnext%id == list_meta%id) list_meta => pmeta

      else

! strange situation, pnext is associated, but list_meta not
!
        call print_error("blocks::intert_metablock_before"                     &
                       , "Argument pnext is associated but list_meta is not!")

      end if

    else

! if pnext is null and last_meta is associated, there is something wrong
!
      if (associated(last_meta)) then

! strange situation, pnext is null but last_meta is associated
!
        call print_error("blocks::intert_metablock_before"                     &
                      , "Argument pnext is null but last_meta is associated!")

      else

! pnext and last_meta are nulls, so add the first block to the list by
! associating list_meta and last_meta
!
        list_meta => pmeta
        last_meta => pmeta

      end if

    end if

!-------------------------------------------------------------------------------
!
  end subroutine insert_metablock_before
!
!===============================================================================
!
! function INCREASE_ID:
! --------------------
!
!   Function increases the last identification number by 1 and returns its
!   value.
!
!
!===============================================================================
!
  function increase_id() result(id)

! local variables are not implicit by default
!
    implicit none

! return variable
!
    integer(kind=4) :: id
!
!-------------------------------------------------------------------------------
!
! increase the last identification number by 1
!
    last_id = last_id + 1

! return its value
!
    id      = last_id

!-------------------------------------------------------------------------------
!
  end function increase_id

!===============================================================================
!
end module
