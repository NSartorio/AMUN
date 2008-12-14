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

! minimum grid step
!
  real, save :: dx_min = 1.0

! space steps for all levels of refinements
!
  real, dimension(:,:), allocatable, save :: ax  , ay  , az
  real, dimension(:  ), allocatable, save :: adx , ady , adz
  real, dimension(:  ), allocatable, save :: adxi, adyi, adzi

  contains
!
!===============================================================================
!
! init_mesh: subroutine initializes mesh by creating blocks according
!            to the geometry, initial problem and refinement criterium
!
!===============================================================================
!
  subroutine init_mesh

    use config , only : iblocks, jblocks, kblocks, ncells             &
                      , xmin, xmax, ymin, ymax, zmin, zmax, maxlev, ngrids
    use blocks , only : list_allocated, init_blocks, clear_blocks     &
                      , allocate_blocks, refine_block, block, nchild, ndims, plist, last_id
    use error  , only : print_info
    use problem, only : init_problem, check_ref

    implicit none

! local variables
!
    type(block), pointer :: pblock, pparent, pchild, pneigh
    integer(kind=4)      :: l, p, i, j, k, n
    character(len=32)    :: bstr, tstr

!----------------------------------------------------------------------
!
! check if the list is allocated, if yes deallocate it
!
    if (list_allocated()) then
      call print_info("mesh::init_mesh", "Block list is allocated, deallocate it!")

      call clear_blocks
    endif

! print information
!
    write(*,"(1x,a)"   ) "Generating initial mesh:"
    write(*,"(4x,a,1x,i6)") "refining to max. level =", maxlev
    write(*,"(4x,a,1x,i6)") "effective resolution   =", ncells*2**maxlev

! allocate initial structure of blocks according the the defined geometry
!
! TODO: by default we initiate 2x2=4 blocks in N configuration
! TODO: in the future allow user to define an arbitrary shape
!
    call init_blocks
    call allocate_blocks('N', xmin, xmax, ymin, ymax, zmin, zmax)

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out these blocks with the initial condition
!
    pblock => plist
    do while (associated(pblock))

! set level
!
      pblock%level    = 1
      pblock%refine   = 0

! set initial conditions
!
      call init_problem(pblock)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! at this point the inital blocks are allocated and set for refinement,
! so iterate over all levels from 1 to maxlevel and create sub-blocks,
! set the initial conditions for each, check criterium and set for
! refinement according to the criterium fullfilment
!
! TODO: iterate over all levels, set neighbors of refined blocks for
!       refinement too, to keep the level jump not larger then 1,
!       refine blocks, set inital conditions at newly created block,
!       and finally check the criterium
!
    write(*,"(4x,a,$)") "refining level         =    "
    do l = 1, maxlev-1

! print information
!
      write(*,"(1x,i2,$)") l

! iterate over all blocks and check refinement criterion
!
      pblock => plist
      do while (associated(pblock))

! check refinements criterion for the current block
!
        if (pblock%leaf .eq. 'T' .and. pblock%level .lt. maxlev) then
          pblock%refine = check_ref(pblock)

! do not allow for derefinement initially
!
          if (pblock%refine .eq. -1) &
            pblock%refine = 0
        else
          pblock%refine = 0
        endif

! assign pointer to the next block
!
        pblock => pblock%next

      end do

! iterate over all blocks and select the neighbors of refined blocks for refinement
!
      do n = l, 2, -1
        pblock => plist
        do while (associated(pblock))

! check if block is a leaf and selected for refinement
!
          if (pblock%refine .eq. 1 .and. pblock%level .eq. n) then

! iterate over all neighbors
!
            do i = 1, ndims
              do j = 1, 2
                do k = 1, 2

                  pneigh => pblock%pneigh(i,j,k)%p

! check if neighbor is associated
!
                  if (associated(pneigh)) then
                    if (pneigh%leaf .eq. 'T') then
                      if (pneigh%level .lt. pblock%level) &
                        pneigh%refine = 1
                    else
                      print *, pneigh%id, 'is not a leaf'
                    endif
                  endif

                end do
              end do
            end do

          endif

! assign pointer to the next block
!
          pblock => pblock%next

        end do
      end do

! iterate over all blocks and refine selected ones
!
      do n = 1, l
        pblock => plist
        do while (associated(pblock))

! check if block needs to be refined and if it is a leaf
!
          if (pblock%refine .eq. 1 .and. pblock%level .eq. n) then

! refine this block
!
            call refine_block(pblock)

! iterate over all children
!
            pparent => pblock%parent
            do p = 1, nchild

              pchild => pparent%child(p)%p

! initialize problem for that children and check the refinement criterion
!
              if (associated(pchild)) &
                call init_problem(pchild)

            end do

          endif

! assign pointer to the next block
!
          pblock => pblock%next

        end do
      end do

    end do

! print information
!
    write(bstr,"(i)") last_id
    write(tstr,"(i)") (2**maxlev)**ndims
    write(*,*)
    write(*,"(4x,a,1x,a6,' / ',a,' = ',f8.4,' %')") "allocated/total blocks =", trim(adjustl(bstr)),trim(adjustl(tstr)), (100.0*last_id)/(2**maxlev)**ndims

! allocating space for coordinate variables
!
    allocate(ax  (maxlev, ngrids))
    allocate(ay  (maxlev, ngrids))
    allocate(az  (maxlev, ngrids))
    allocate(adx (maxlev))
    allocate(ady (maxlev))
    allocate(adz (maxlev))
    allocate(adxi(maxlev))
    allocate(adyi(maxlev))
    allocate(adzi(maxlev))

! generating coordinates for all levels
!
    do l = 1, maxlev
      adx (l) = (xmax - xmin) / (ncells*2**l)
      adxi(l) = 1.0 / adx(l)
      ady (l) = (ymax - ymin) / (ncells*2**l)
      adyi(l) = 1.0 / ady(l)
#if NDIMS == 3
      adz (l) = (zmax - zmin) / (ncells*2**l)
#else
      adz (l) = 1.0
#endif
      adzi(l) = 1.0 / adz(l)
    end do

! get the minimum grid step
!
    dx_min = 0.5*min(adx(maxlev), ady(maxlev), adz(maxlev))

!-------------------------------------------------------------------------------
!
  end subroutine init_mesh
!
!===============================================================================
!
! update_mesh: subroutine check the refinement criterion for each block,
!              refines or derefines it if necessary, and restricts or
!              prolongates all data to the newly created blocks
!
!===============================================================================
!
  subroutine update_mesh(ref)

    use config , only : maxlev
    use blocks , only : block, plist, ndims, nchild, refine_block
    use error  , only : print_info
    use problem, only : check_ref

    implicit none

    integer, intent(in) :: ref

! local variables
!
    type(block), pointer :: pblock, pneigh, pchild, pparent
    integer(kind=4)      :: l, i, j, k, n, p

!----------------------------------------------------------------------
!
!     if (ref .eq. 1) then

    do l = 1, maxlev

! check refinement criterion
!
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l) then
          pblock%refine = 0

          if (pblock%leaf .eq. 'T') then
            pblock%refine = check_ref(pblock)

!             if (pblock%refine .eq. -1) &
!               pblock%refine = 0
            if (pblock%refine .eq. -1 .and. pblock%level .eq. 1) &
              pblock%refine = 0

            if (pblock%refine .eq. 1 .and. pblock%level .eq. maxlev) &
              pblock%refine = 0
          endif
        endif

        pblock => pblock%next
      end do

! refinement conditions for blocks:
!
! - all neighbors must be at the same or lower level
!
      do n = l, 2, -1
        pblock => plist
        do while (associated(pblock))

          if (pblock%level .eq. n) then
            if (pblock%refine .eq. 1 .and. pblock%leaf .eq. 'T') then
              do i = 1, ndims
                do j = 1, 2
                  do k = 1, 2
                    pneigh => pblock%pneigh(i,j,k)%p
                    if (associated(pneigh)) then
                      if (pneigh%level .lt. pblock%level) &
                        pneigh%refine = 1
                    endif

                  end do
                end do
              end do
            endif
          endif

          pblock => pblock%next
        end do
      end do

      do n = 1, l
        pblock => plist
        do while (associated(pblock))

          if (pblock%level .eq. n) then
            if (pblock%refine .eq. 1 .and. pblock%leaf .eq. 'T') then

              pparent => pblock

              print *, 'refine block ', pparent%id
              call refine_block(pblock)
              call prolong_block(pparent)
            endif
          endif

          pblock => pblock%next

        end do
      end do

    end do

! else

    do l = maxlev, 1, -1

! check refinement criterion
!
!       pblock => plist
!       do while (associated(pblock))
!
!         if (pblock%level .eq. l) then
!           pblock%refine = 0
!
!           if (pblock%leaf .eq. 'T') then
!             pblock%refine = check_ref(pblock)
!
!             if (pblock%refine .eq. -1 .and. pblock%level .eq. 1) &
!               pblock%refine = 0
!
!             if (pblock%refine .eq. 1) &
!               pblock%refine = 0
! !             if (pblock%refine .eq. 1 .and. pblock%level .eq. maxlev) &
! !               pblock%refine = 0
!           endif
!         endif
!
!         pblock => pblock%next
!       end do

! derefinement conditions for blocks:
!
! - all neighbors must be at the same or lower level
! - all neighbors at the same level must be not selected for refinement
! - all sibling must be at the same level and marked for derefinement
!
! if all above conditions are fulfilled, select the refinement of
! the parent of the current block to -1
!
! if at least one is not fulfilled, set all sibling to not do the refinement
!
      do n = l, maxlev
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. n) then
          if (pblock%refine .eq. -1 .and. pblock%leaf .eq. 'T') then
            do i = 1, ndims
              do j = 1, 2
                do k = 1, 2
                  pneigh => pblock%pneigh(i,j,k)%p
                  if (associated(pneigh)) then
                    if (pneigh%level .gt. pblock%level) &
                      pblock%refine = 0
                    if (pneigh%level .eq. pblock%level .and. pneigh%refine .eq. 1) &
                      pblock%refine = 0
                  endif

                end do
              end do
            end do

            pparent => pblock%parent
            if (associated(pparent)) then
              do p = 1, nchild
                pchild => pparent%child(p)%p

                if (associated(pchild)) then
                  if (pchild%refine .ne. -1 .or. pchild%level .ne. pblock%level) &
                    pblock%refine = 0
                endif
              end do

              if (pblock%refine .ne. -1) then
                do p = 1, nchild
                  pchild => pparent%child(p)%p

                  if (associated(pchild)) then
                    if (pchild%refine .eq. -1 .and. pchild%level .eq. pblock%level) &
                      pchild%refine = 0
                  endif
                end do
              endif
            endif

          endif
        endif

        pparent => pblock%parent
        if (associated(pparent) .and. pblock%refine .eq. -1) &
          pparent%refine = -1

! point to the next block
!
        pblock => pblock%next
      end do
      end do

! perform derefinements for all selected blocks
!
      do n = maxlev, l, -1
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. n) then
          if (pblock%refine .eq. -1 .and. pblock%leaf .eq. 'T') then

            pparent => pblock%parent

            if (associated(pparent) .and. pparent%refine .eq. -1) then
              print *, 'derefine block ', pparent%id
              call restrict_block(pparent)
              pblock => pparent
            endif
          endif
        endif

        pblock => pblock%next
      end do
      end do

    end do

! endif

!-------------------------------------------------------------------------------
!
  end subroutine update_mesh
!
!===============================================================================
!
! prolong_block: subroutine expands the block data and copy them to children
!
!===============================================================================
!
  subroutine prolong_block(pblock)

    use blocks       , only : block, nv => nvars
    use config       , only : in, jn, kn, im, jm, km, ng
    use interpolation, only : expand

    implicit none

! input arguments
!
    type(block), pointer :: pblock, pchild

! local variables
!
    integer :: q, dm(3), fm(3), il, iu, jl, ju, k

! local arrays
!
    real, dimension(nv,2*im,2*jm,2*km) :: u

!-------------------------------------------------------------------------------
!
    dm(1) = im
    dm(2) = jm
    dm(3) = km
    fm(:) = 2 * (dm(:) - ng)
    if (km .eq. 1) &
      fm(3) = 1

! first expand variables
!
    do q = 1, nv
      call expand(dm,fm,ng,pblock%u(q,:,:,:),u(q,:,:,:),'t','t','t')
    end do

! fill values of children
!
    il = 1
    iu = il + im - 1
    jl = 1
    ju = jl + jm - 1
    do k = 1, km
      pblock%child(1)%p%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = in + 1
    iu = il + im - 1
    jl = 1
    ju = jl + jm - 1
    do k = 1, km
      pblock%child(2)%p%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = 1
    iu = il + im - 1
    jl = in + 1
    ju = jl + jm - 1
    do k = 1, km
      pblock%child(3)%p%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = in + 1
    iu = il + im - 1
    jl = in + 1
    ju = jl + jm - 1
    do k = 1, km
      pblock%child(4)%p%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

!-------------------------------------------------------------------------------
!
  end subroutine prolong_block
!
!===============================================================================
!
! restrict_block: subroutine shrinks the block data and copy them from children
!
!===============================================================================
!
  subroutine restrict_block(pblock)

    use blocks       , only : block, nv => nvars, derefine_block
    use config       , only : in, jn, kn, im, jm, km, ng

    implicit none

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

! local variables
!
    integer :: q, i, j, k, il, iu, jl, ju, i1, i2, j1, j2

! local pointers
!
    type(block), pointer :: pbl, pbr, ptl, ptr, pb

!-------------------------------------------------------------------------------
!
    pbl => pblock%child(1)%p
    pbr => pblock%child(2)%p
    ptl => pblock%child(3)%p
    ptr => pblock%child(4)%p

! BL
!
    il = ng / 2 + 1
    iu = im / 2
    jl = ng / 2 + 1
    ju = jm / 2
    do k = 1, km
      do j = jl, ju
        j2 = 2 * j - ng
        j1 = j2 - 1
        do i = il, iu
          i2 = 2 * i - ng
          i1 = i2 - 1

          pblock%u(1:nv,i,j,k) = 0.25 * (pbl%u(1:nv,i1,j1,k) + pbl%u(1:nv,i1,j2,k) &
                                       + pbl%u(1:nv,i2,j1,k) + pbl%u(1:nv,i2,j2,k))
        end do
      end do
    end do

! BR
!
    il = im / 2 + 1
    iu = im - ng / 2
    jl = ng / 2 + 1
    ju = jm / 2
    do k = 1, km
      do j = jl, ju
        j2 = 2 * j - ng
        j1 = j2 - 1
        do i = il, iu
          i2 = 2 * i - im + ng
          i1 = i2 - 1

          pblock%u(1:nv,i,j,k) = 0.25 * (pbr%u(1:nv,i1,j1,k) + pbr%u(1:nv,i1,j2,k) &
                                       + pbr%u(1:nv,i2,j1,k) + pbr%u(1:nv,i2,j2,k))
        end do
      end do
    end do

! TL
!
    il = ng / 2 + 1
    iu = im / 2
    jl = jm / 2 + 1
    ju = jm - ng / 2
    do k = 1, km
      do j = jl, ju
        j2 = 2 * j - jm + ng
        j1 = j2 - 1
        do i = il, iu
          i2 = 2 * i - ng
          i1 = i2 - 1

          pblock%u(1:nv,i,j,k) = 0.25 * (ptl%u(1:nv,i1,j1,k) + ptl%u(1:nv,i1,j2,k) &
                                       + ptl%u(1:nv,i2,j1,k) + ptl%u(1:nv,i2,j2,k))
        end do
      end do
    end do

! TR
!
    il = im / 2 + 1
    iu = im - ng / 2
    jl = jm / 2 + 1
    ju = jm - ng / 2
    do k = 1, km
      do j = jl, ju
        j2 = 2 * j - jm + ng
        j1 = j2 - 1
        do i = il, iu
          i2 = 2 * i - im + ng
          i1 = i2 - 1

          pblock%u(1:nv,i,j,k) = 0.25 * (ptr%u(1:nv,i1,j1,k) + ptr%u(1:nv,i1,j2,k) &
                                       + ptr%u(1:nv,i2,j1,k) + ptr%u(1:nv,i2,j2,k))
        end do
      end do
    end do

! derefine block
!
    pb => pblock
    call derefine_block(pb)

    nullify(pb)
    nullify(pbl)
    nullify(pbr)
    nullify(ptl)
    nullify(ptr)

!-------------------------------------------------------------------------------
!
  end subroutine restrict_block
!
!===============================================================================
!
! clears_mesh: subroutine deallocates mesh, removing blocks
!
!===============================================================================
!
  subroutine clear_mesh

    use blocks, only : clear_blocks
    use error , only : print_info

    implicit none

!-------------------------------------------------------------------------------
!
! deallocate block structure
!
    call clear_blocks

! deallocating coordinate variables
!
    deallocate(ax)
    deallocate(ay)
    deallocate(az)
    deallocate(adx)
    deallocate(ady)
    deallocate(adz)
    deallocate(adxi)
    deallocate(adyi)
    deallocate(adzi)

!-------------------------------------------------------------------------------
!
  end subroutine clear_mesh

!===============================================================================
!
end module
