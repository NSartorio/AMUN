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

  type block
    type(block), pointer :: next, prev
    type(blockref)       :: parent, child(nchild), neigh(ndims,2,2)

    logical              :: leaf
    integer(kind=4)      :: cpu, id, level

    character            :: config
    integer(kind=4)      :: refine


    real                 :: xmin, xmax, ymin, ymax, zmin, zmax

    real, dimension(:,:,:,:), allocatable :: u
    real, dimension(:,:,:)  , allocatable :: c
  end type block

! array of ID to pointer conversion
!
  type(blockptr), dimension(:), allocatable, save :: idtoptr

! stored pointers
!
  type(block), pointer, save :: plist, plast

! stored last id (should always increase)
!
  integer(kind=4)     , save :: last_id

! store number of allocated blocks
!
  integer(kind=4)     , save :: nblocks

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
! first check if block list is empty
!
    if (associated(plist)) &
      call print_warning("blocks::init_blocks", "Block list already associated!")

! nullify all pointers
!
    nullify(plist)

! reset number of blocks
!
    nblocks = 0

! reset ID
!
    last_id = 0

! allocate space for ID to pointer array
!
    allocate(idtoptr(maxid))

! nullify all pointers in ID to pointer array
!
    do p = 1, maxid
      nullify(idtoptr(p)%ptr)
    end do

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
    type(block), pointer :: pnext
!
!-------------------------------------------------------------------------------
!
! untill the list is free, iterate over all chunks and deallocate blocks
!
    do while(associated(plist))

! assign temporary pointer to the next chunk
!
      pnext => plist%next

! deallocate and nullify the current block
!
      call deallocate_block(plist)

! assign pointer to the current chunk
!
      plist => pnext

    end do

! deallocate ID to pointer conversion array
!
    deallocate(idtoptr)

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

    use config, only : im, jm, km

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

!-------------------------------------------------------------------------------
!
  end subroutine deallocate_block
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
! refine_block: subroutine refines selected block
!
!===============================================================================
!
  subroutine refine_block(pblock)

    use error, only : print_error

    implicit none

! input parameters
!
    type(block), pointer, intent(inout) :: pblock

! pointers
!
    type(block), pointer :: pb, pbl, pbr, ptl, ptr, pneigh
!
!-------------------------------------------------------------------------------
!
! check if pointer is associated
!
    if (associated(pblock)) then

! unset the refinement and leaf flags for the parent block
!
      pblock%refine = 0
      pblock%leaf   = .false.

! create 4 blocks
!
      call allocate_block(pbl)
      call allocate_block(pbr)
      call allocate_block(ptl)
      call allocate_block(ptr)

! set parent
!
      pbl%parent%id = pblock%id
      pbr%parent%id = pblock%id
      ptl%parent%id = pblock%id
      ptr%parent%id = pblock%id

! set level
!
      pbl%level = pblock%level + 1
      pbr%level = pbl%level
      ptl%level = pbl%level
      ptr%level = pbl%level

! set leaf flags
!
      pbl%leaf   = .true.
      pbr%leaf   = .true.
      ptl%leaf   = .true.
      ptr%leaf   = .true.

! set bounds
!
      pbl%xmin = pblock%xmin
      pbl%xmax = 0.5 * (pblock%xmax + pblock%xmin)
      pbl%ymin = pblock%ymin
      pbl%ymax = 0.5 * (pblock%ymax + pblock%ymin)

      pbr%xmin = pbl%xmax
      pbr%xmax = pblock%xmax
      pbr%ymin = pblock%ymin
      pbr%ymax = pbl%ymax

      ptl%xmin = pblock%xmin
      ptl%xmax = pbl%xmax
      ptl%ymin = pbl%ymax
      ptl%ymax = pblock%ymax

      ptr%xmin = ptl%xmax
      ptr%xmax = pblock%xmax
      ptr%ymin = ptl%ymin
      ptr%ymax = pblock%ymax

! set neighbors to the refined blocks
!
      pbl%neigh(1,2,1)%id = pbr%id  ! BL right  -> BR
      pbl%neigh(1,2,2)%id = pbr%id
      pbl%neigh(2,2,1)%id = ptl%id  ! BL top    -> TL
      pbl%neigh(2,2,2)%id = ptl%id

      pbr%neigh(1,1,1)%id = pbl%id  ! BR left   -> BL
      pbr%neigh(1,1,2)%id = pbl%id
      pbr%neigh(2,2,1)%id = ptr%id  ! BR top    -> TR
      pbr%neigh(2,2,2)%id = ptr%id

      ptl%neigh(1,2,1)%id = ptr%id  ! TL right  -> TR
      ptl%neigh(1,2,2)%id = ptr%id
      ptl%neigh(2,1,1)%id = pbl%id  ! TL bottom -> BL
      ptl%neigh(2,1,2)%id = pbl%id

      ptr%neigh(1,1,1)%id = ptl%id  ! TR left   -> TL
      ptr%neigh(1,1,2)%id = ptl%id
      ptr%neigh(2,1,1)%id = pbr%id  ! TR bottom -> BR
      ptr%neigh(2,1,2)%id = pbr%id

! set pointer to the neighbors of the parent block
!
      pneigh => get_pointer(pblock%neigh(1,1,1)%id) ! left lower neighbor
      if (associated(pneigh)) then
        pbl%neigh(1,1,1)%id = pneigh%id
        pbl%neigh(1,1,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,1)%id = pbl%id
          pneigh%neigh(1,2,2)%id = ptl%id
        endif
        if (pneigh%level .eq. pbl%level) then
          pneigh%neigh(1,2,1)%id = pbl%id
          pneigh%neigh(1,2,2)%id = pbl%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(1,1,2)%id) ! left upper neighbor
      if (associated(pneigh)) then
        ptl%neigh(1,1,1)%id = pneigh%id
        ptl%neigh(1,1,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,2,1)%id = pbl%id
          pneigh%neigh(1,2,2)%id = ptl%id
        endif
        if (pneigh%level .eq. ptl%level) then
          pneigh%neigh(1,2,1)%id = ptl%id
          pneigh%neigh(1,2,2)%id = ptl%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(1,2,1)%id) ! right lower neighbor
      if (associated(pneigh)) then
        pbr%neigh(1,2,1)%id = pneigh%id
        pbr%neigh(1,2,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,1)%id = pbr%id
          pneigh%neigh(1,1,2)%id = ptr%id
        endif
        if (pneigh%level .eq. pbr%level) then
          pneigh%neigh(1,1,1)%id = pbr%id
          pneigh%neigh(1,1,2)%id = pbr%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(1,2,2)%id) ! right upper neighbor
      if (associated(pneigh)) then
        ptr%neigh(1,2,1)%id = pneigh%id
        ptr%neigh(1,2,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(1,1,1)%id = pbr%id
          pneigh%neigh(1,1,2)%id = ptr%id
        endif
        if (pneigh%level .eq. ptr%level) then
          pneigh%neigh(1,1,1)%id = ptr%id
          pneigh%neigh(1,1,2)%id = ptr%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(2,1,1)%id)  ! bottom left neighbor
      if (associated(pneigh)) then
        pbl%neigh(2,1,1)%id = pneigh%id
        pbl%neigh(2,1,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,1)%id = pbl%id
          pneigh%neigh(2,2,2)%id = pbr%id
        endif
        if (pneigh%level .eq. pbl%level) then
          pneigh%neigh(2,2,1)%id = pbl%id
          pneigh%neigh(2,2,2)%id = pbl%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(2,1,2)%id)  ! bottom right neighbor
      if (associated(pneigh)) then
        pbr%neigh(2,1,1)%id = pneigh%id
        pbr%neigh(2,1,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,2,1)%id = pbl%id
          pneigh%neigh(2,2,2)%id = pbr%id
        endif
        if (pneigh%level .eq. pbr%level) then
          pneigh%neigh(2,2,1)%id = pbr%id
          pneigh%neigh(2,2,2)%id = pbr%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(2,2,1)%id)  ! top left neighbor
      if (associated(pneigh)) then
        ptl%neigh(2,2,1)%id = pneigh%id
        ptl%neigh(2,2,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,1)%id = ptl%id
          pneigh%neigh(2,1,2)%id = ptr%id
        endif
        if (pneigh%level .eq. ptl%level) then
          pneigh%neigh(2,1,1)%id = ptl%id
          pneigh%neigh(2,1,2)%id = ptl%id
        endif
      endif

      pneigh => get_pointer(pblock%neigh(2,2,2)%id)  ! top right neighbor
      if (associated(pneigh)) then
        ptr%neigh(2,2,1)%id = pneigh%id
        ptr%neigh(2,2,2)%id = pneigh%id

        if (pneigh%level .eq. pblock%level) then
          pneigh%neigh(2,1,1)%id = ptl%id
          pneigh%neigh(2,1,2)%id = ptr%id
        endif
        if (pneigh%level .eq. ptr%level) then
          pneigh%neigh(2,1,1)%id = ptr%id
          pneigh%neigh(2,1,2)%id = ptr%id
        endif
      endif

! set children
!
      pblock%child(1)%id = pbl%id
      pblock%child(2)%id = pbr%id
      pblock%child(3)%id = ptl%id
      pblock%child(4)%id = ptr%id

! depending on the configuration of the parent block
!
      select case(pblock%config)
      case('z', 'Z')

! set blocks configurations
!
        pbl%config = 'Z'
        pbr%config = 'Z'
        ptl%config = 'Z'
        ptr%config = 'Z'

! connect blocks in a chain
!
        pbl%next => pbr
        pbr%next => ptl
        ptl%next => ptr

        pbr%prev => pbl
        ptl%prev => pbr
        ptr%prev => ptl

! insert this chain after the parent block
!
        pb => pblock%next
        if (associated(pb)) then
          pb%prev => ptr
          ptr%next => pb
        else
          plast => ptr
          nullify(ptr%next)
        endif
        pblock%next => pbl
        pbl%prev => pblock

        pblock => ptr

      case('n', 'N')

! set blocks configurations
!
        pbl%config = 'D'
        ptl%config = 'N'
        ptr%config = 'N'
        pbr%config = 'C'

! connect blocks in a chain
!
        pbl%next => ptl
        ptl%next => ptr
        ptr%next => pbr

        ptl%prev => pbl
        ptr%prev => ptl
        pbr%prev => ptr

! insert this chain after the parent the block
!
        pb => pblock%next
        if (associated(pb)) then
          pb%prev => pbr
          pbr%next => pb
        endif
        pbl%prev => pblock
        pblock%next => pbl

        pblock => pbr

      case('d', 'D')

! set blocks configurations
!
        pbl%config = 'N'
        pbr%config = 'D'
        ptr%config = 'D'
        ptl%config = 'U'

! connect blocks in a chain
!
        pbl%next => pbr
        pbr%next => ptr
        ptr%next => ptl

        pbr%prev => pbl
        ptr%prev => pbr
        ptl%prev => ptr

! insert this chain in the block list
!
        pb => pblock%next
        if (associated(pb)) then
          pb%prev => ptl
          ptl%next => pb
        endif
        pbl%prev => pblock
        pblock%next => pbl

        pblock => ptl

      case('c', 'C')

! set blocks configurations
!
        ptr%config = 'U'
        ptl%config = 'C'
        pbl%config = 'C'
        pbr%config = 'N'

! connect blocks in a chain
!
        ptr%next => ptl
        ptl%next => pbl
        pbl%next => pbr

        ptl%prev => ptr
        pbl%prev => ptl
        pbr%prev => pbl

! insert this chain in the block list
!
        pb => pblock%next
        if (associated(pb)) then
          pb%prev => pbr
          pbr%next => pb
        endif
        ptr%prev => pblock
        pblock%next => ptr

        pblock => pbr

      case('u', 'U')

! set blocks configurations
!
        ptr%config = 'C'
        pbr%config = 'U'
        pbl%config = 'U'
        ptl%config = 'D'

! connect blocks in a chain
!
        ptr%next => pbr
        pbr%next => pbl
        pbl%next => ptl

        pbr%prev => ptr
        pbl%prev => pbr
        ptl%prev => pbl

! insert this chain in the block list
!
        pb => pblock%next
        if (associated(pb)) then
          pb%prev => ptl
          ptl%next => pb
        endif
        ptr%prev => pblock
        pblock%next => ptr

        pblock => ptl

      end select

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
    pblock%neigh(1,1,1)%id = pbl%neigh(1,1,1)%id
    pblock%neigh(1,1,2)%id = ptl%neigh(1,1,2)%id
    pblock%neigh(1,2,1)%id = pbr%neigh(1,2,1)%id
    pblock%neigh(1,2,2)%id = ptr%neigh(1,2,2)%id
    pblock%neigh(2,1,1)%id = pbl%neigh(2,1,1)%id
    pblock%neigh(2,1,2)%id = pbr%neigh(2,1,2)%id
    pblock%neigh(2,2,1)%id = ptl%neigh(2,2,1)%id
    pblock%neigh(2,2,2)%id = ptr%neigh(2,2,2)%id

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
!
!===============================================================================
!
! allocate_blocks: subroutine allocates a configuration of blocks
!
! TODO: move to/replace with a subroutine in problem module; this will allow
!       a user to define his/her own initial mesh
!
!===============================================================================
!
  subroutine allocate_blocks(block_config, xmn, xmx, ymn, ymx &
                           , zmn, zmx)

    use config, only : xlbndry, xubndry, ylbndry, yubndry
    use error , only : print_error

    implicit none

! input parameters
!
    character(len=1), intent(in) :: block_config
    real            , intent(in) :: xmn, xmx, ymn, ymx, zmn, zmx

! local pointers
!
    type(block), pointer :: pbl, pbr, ptl, ptr

! local variables
!
    real :: xl, xc, xr, yl, yc, yr, zl, zc, zr
!
!-------------------------------------------------------------------------------
!
    select case(block_config)
    case('z', 'Z')

! create root blocks
!
      call append_block(pbl)
      call append_block(pbr)
      call append_block(ptl)
      call append_block(ptr)

! set configurations
!
      pbl%config = 'Z'
      pbr%config = 'Z'
      ptl%config = 'Z'
      ptr%config = 'Z'

! copy pointer of the first block in chain
!
      plist => pbl

    case('n', 'N')

! create root blocks
!
      call append_block(pbl)
      call append_block(ptl)
      call append_block(ptr)
      call append_block(pbr)

! set configurations
!
      pbl%config = 'D'
      ptl%config = 'N'
      ptr%config = 'N'
      pbr%config = 'C'

! copy pointer of the first block in chain
!
      plist => pbl

    case default
      call print_error("blocks::allocate_blocks","Configuration '" // block_config // "' not supported! Terminating!")
    end select

! set leaf flags
!
    pbl%leaf   = .true.
    pbr%leaf   = .true.
    ptl%leaf   = .true.
    ptr%leaf   = .true.

! set neighbors
!
    if (xlbndry .eq. 'periodic') &
      pbl%neigh(1,1,:)%id = pbr%id
    pbl%neigh(1,2,:)%id = pbr%id
    if (ylbndry .eq. 'periodic') &
      pbl%neigh(2,1,:)%id = ptl%id
    pbl%neigh(2,2,:)%id = ptl%id

    pbr%neigh(1,1,:)%id = pbl%id
    if (xubndry .eq. 'periodic') &
      pbr%neigh(1,2,:)%id = pbl%id
    if (ylbndry .eq. 'periodic') &
      pbr%neigh(2,1,:)%id = ptr%id
    pbr%neigh(2,2,:)%id = ptr%id

    if (xlbndry .eq. 'periodic') &
      ptl%neigh(1,1,:)%id = ptr%id
    ptl%neigh(1,2,:)%id = ptr%id
    ptl%neigh(2,1,:)%id = pbl%id
    if (yubndry .eq. 'periodic') &
      ptl%neigh(2,2,:)%id = pbl%id

    ptr%neigh(1,1,:)%id = ptl%id
    if (xubndry .eq. 'periodic') &
      ptr%neigh(1,2,:)%id = ptl%id
    ptr%neigh(2,1,:)%id = pbr%id
    if (yubndry .eq. 'periodic') &
      ptr%neigh(2,2,:)%id = pbr%id

! set block bounds
!
    xl = xmn
    xc = 0.5 * (xmx + xmn)
    xr = xmx
    yl = ymn
    yc = 0.5 * (ymx + ymn)
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

!-------------------------------------------------------------------------------
!
  end subroutine allocate_blocks

!===============================================================================
!
end module
