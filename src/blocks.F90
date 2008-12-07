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
  type bpointer
    type(block), pointer :: p
  end type bpointer

  type block
    type(block), pointer :: next, prev, parent
    type(bpointer)       :: child(nchild), pneigh(ndims,2,2)

    character            :: config, leaf
    integer(kind=4)      :: refine

    integer(kind=4)      :: id, level
    integer(kind=4)      :: neigh(ndims,2,2)

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

! local variables
!
    integer :: i, j, k
!
!----------------------------------------------------------------------
!
! allocate block structure
!
    allocate(pblock)

! set unique ID
!
    pblock%id = increase_id()

! set configuration and leaf flags
!
    pblock%config = 'N'
    pblock%leaf   = 'F'

! initialize refinement flag
!
    pblock%refine = 0

! nullify pointers
!
    nullify(pblock%next)
    nullify(pblock%prev)
    nullify(pblock%parent)

! reset neighbors
!
    pblock%neigh(:,:,:) = -1
    do i = 1, ndims
      do j = 1, 2
        do k = 1, 2
          nullify(pblock%pneigh(i,j,k)%p)
        end do
      end do
    end do

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
      pgroup => pbl

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
      pgroup => pbl

    case default
      call print_error("blocks::allocate_blocks","Configuration '" // block_config // "' not supported! Terminating!")
    end select

! set leaf flags
!
    pbl%leaf   = 'T'
    pbr%leaf   = 'T'
    ptl%leaf   = 'T'
    ptr%leaf   = 'T'

! set neighbors
!
    pbl%neigh(1,2,:) = pbr%id
    pbl%neigh(2,2,:) = ptl%id

    pbr%neigh(1,1,:) = pbl%id
    pbr%neigh(2,2,:) = ptr%id

    ptl%neigh(1,2,:) = ptr%id
    ptl%neigh(2,1,:) = pbl%id

    ptr%neigh(1,1,:) = ptl%id
    ptr%neigh(2,1,:) = pbr%id

! set neighbor pointers
!
    pbl%pneigh(1,2,1)%p => pbr  ! BL right  -> BR
    pbl%pneigh(1,2,2)%p => pbr
    pbl%pneigh(2,2,1)%p => ptl  ! BL top    -> TL
    pbl%pneigh(2,2,2)%p => ptl

    pbr%pneigh(1,1,1)%p => pbl  ! BR left   -> BL
    pbr%pneigh(1,1,2)%p => pbl
    pbr%pneigh(2,2,1)%p => ptr  ! BR top    -> TR
    pbr%pneigh(2,2,2)%p => ptr

    ptl%pneigh(1,2,1)%p => ptr  ! TL right  -> TR
    ptl%pneigh(1,2,2)%p => ptr
    ptl%pneigh(2,1,1)%p => pbl  ! TL bottom -> BL
    ptl%pneigh(2,1,2)%p => pbl

    ptr%pneigh(1,1,1)%p => ptl  ! TR left   -> TL
    ptr%pneigh(1,1,2)%p => ptl
    ptr%pneigh(2,1,1)%p => pbr  ! TR bottom -> BR
    ptr%pneigh(2,1,2)%p => pbr

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

! assign temporary pointer to the next chunk
!
      pcurr => pfirst%next

! deallocate the content of current block
!
!       write (*,"(i9.9,2x,i2,1x,4(1x,f6.3))") pfirst%id, pfirst%level, pfirst%xmin, pfirst%xmax, pfirst%ymin, pfirst%ymax

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
!
!======================================================================
!
! refine_block: subroutine refines selected block
!
!======================================================================
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
!----------------------------------------------------------------------
!
! check if pointer is associated
!
    if (associated(pblock)) then

! unset the refinement and leaf flags for the parent block
!
      pblock%refine = 0
      pblock%leaf   = 'F'

! create 4 blocks
!
      call allocate_block(pbl)
      call allocate_block(pbr)
      call allocate_block(ptl)
      call allocate_block(ptr)

! set parent
!
      pbl%parent => pblock
      pbr%parent => pblock
      ptl%parent => pblock
      ptr%parent => pblock

! set level
!
      pbl%level = pblock%level + 1
      pbr%level = pbl%level
      ptl%level = pbl%level
      ptr%level = pbl%level

! set leaf flags
!
      pbl%leaf   = 'T'
      pbr%leaf   = 'T'
      ptl%leaf   = 'T'
      ptr%leaf   = 'T'

! set bounds
!
      pbl%xmin = pblock%xmin
      pbl%xmax = 0.5 * (pblock%xmin + pblock%xmax)
      pbl%ymin = pblock%ymin
      pbl%ymax = 0.5 * (pblock%ymin + pblock%ymax)

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

! set neighbor pointers to the refined blocks
!
      pbl%pneigh(1,2,1)%p => pbr  ! BL right  -> BR
      pbl%pneigh(1,2,2)%p => pbr
      pbl%pneigh(2,2,1)%p => ptl  ! BL top    -> TL
      pbl%pneigh(2,2,2)%p => ptl

      pbr%pneigh(1,1,1)%p => pbl  ! BR left   -> BL
      pbr%pneigh(1,1,2)%p => pbl
      pbr%pneigh(2,2,1)%p => ptr  ! BR top    -> TR
      pbr%pneigh(2,2,2)%p => ptr

      ptl%pneigh(1,2,1)%p => ptr  ! TL right  -> TR
      ptl%pneigh(1,2,2)%p => ptr
      ptl%pneigh(2,1,1)%p => pbl  ! TL bottom -> BL
      ptl%pneigh(2,1,2)%p => pbl

      ptr%pneigh(1,1,1)%p => ptl  ! TR left   -> TL
      ptr%pneigh(1,1,2)%p => ptl
      ptr%pneigh(2,1,1)%p => pbr  ! TR bottom -> BR
      ptr%pneigh(2,1,2)%p => pbr

! set pointer to the neighbors of the parent block
!
      pneigh => pblock%pneigh(1,1,1)%p ! left lower neighbor
      if (associated(pneigh)) then
        pbl%pneigh(1,1,1)%p => pneigh
        pbl%pneigh(1,1,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(1,2,1)%p => pbl
          pneigh%pneigh(1,2,2)%p => ptl
        endif
        if (pneigh%level .eq. pbl%level) then
          pneigh%pneigh(1,2,1)%p => pbl
          pneigh%pneigh(1,2,2)%p => pbl
        endif
      endif

      pneigh => pblock%pneigh(1,1,2)%p ! left upper neighbor
      if (associated(pneigh)) then
        ptl%pneigh(1,1,1)%p => pneigh
        ptl%pneigh(1,1,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(1,2,1)%p => pbl
          pneigh%pneigh(1,2,2)%p => ptl
        endif
        if (pneigh%level .eq. ptl%level) then
          pneigh%pneigh(1,2,1)%p => ptl
          pneigh%pneigh(1,2,2)%p => ptl
        endif
      endif

      pneigh => pblock%pneigh(1,2,1)%p ! right lower neighbor
      if (associated(pneigh)) then
        pbr%pneigh(1,2,1)%p => pneigh
        pbr%pneigh(1,2,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(1,1,1)%p => pbr
          pneigh%pneigh(1,1,2)%p => ptr
        endif
        if (pneigh%level .eq. pbr%level) then
          pneigh%pneigh(1,1,1)%p => pbr
          pneigh%pneigh(1,1,2)%p => pbr
        endif
      endif

      pneigh => pblock%pneigh(1,2,2)%p ! right upper neighbor
      if (associated(pneigh)) then
        ptr%pneigh(1,2,1)%p => pneigh
        ptr%pneigh(1,2,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(1,1,1)%p => pbr
          pneigh%pneigh(1,1,2)%p => ptr
        endif
        if (pneigh%level .eq. ptr%level) then
          pneigh%pneigh(1,1,1)%p => ptr
          pneigh%pneigh(1,1,2)%p => ptr
        endif
      endif

      pneigh => pblock%pneigh(2,1,1)%p  ! bottom left neighbor
      if (associated(pneigh)) then
        pbl%pneigh(2,1,1)%p => pneigh
        pbl%pneigh(2,1,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(2,2,1)%p => pbl
          pneigh%pneigh(2,2,2)%p => pbr
        endif
        if (pneigh%level .eq. pbl%level) then
          pneigh%pneigh(2,2,1)%p => pbl
          pneigh%pneigh(2,2,2)%p => pbl
        endif
      endif

      pneigh => pblock%pneigh(2,1,2)%p  ! bottom right neighbor
      if (associated(pneigh)) then
        pbr%pneigh(2,1,1)%p => pneigh
        pbr%pneigh(2,1,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(2,2,1)%p => pbl
          pneigh%pneigh(2,2,2)%p => pbr
        endif
        if (pneigh%level .eq. pbr%level) then
          pneigh%pneigh(2,2,1)%p => pbr
          pneigh%pneigh(2,2,2)%p => pbr
        endif
      endif

      pneigh => pblock%pneigh(2,2,1)%p  ! top left neighbor
      if (associated(pneigh)) then
        ptl%pneigh(2,2,1)%p => pneigh
        ptl%pneigh(2,2,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(2,1,1)%p => ptl
          pneigh%pneigh(2,1,2)%p => ptr
        endif
        if (pneigh%level .eq. ptl%level) then
          pneigh%pneigh(2,1,1)%p => ptl
          pneigh%pneigh(2,1,2)%p => ptl
        endif
      endif

      pneigh => pblock%pneigh(2,2,2)%p  ! top right neighbor
      if (associated(pneigh)) then
        ptr%pneigh(2,2,1)%p => pneigh
        ptr%pneigh(2,2,2)%p => pneigh

        if (pneigh%level .eq. pblock%level) then
          pneigh%pneigh(2,1,1)%p => ptl
          pneigh%pneigh(2,1,2)%p => ptr
        endif
        if (pneigh%level .eq. ptr%level) then
          pneigh%pneigh(2,1,1)%p => ptr
          pneigh%pneigh(2,1,2)%p => ptr
        endif
      endif

! set children
!
      pblock%child(1)%p => pbl
      pblock%child(2)%p => pbr
      pblock%child(3)%p => ptl
      pblock%child(4)%p => ptr

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

!----------------------------------------------------------------------
!
  end subroutine refine_block

!======================================================================
!
end module
