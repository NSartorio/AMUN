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

! spatial coordinates for all levels of refinements
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

    use config  , only : im, jm, km, xmin, xmax, ymin, ymax, zmin, zmax        &
                       , ncells, maxlev
    use blocks  , only : list_allocated, init_blocks, clear_blocks             &
                       , deallocate_block, refine_block, get_pointer           &
                       , block, nchild, ndims, plist, last_id, nblocks, nleafs
    use error   , only : print_info, print_error
    use mpitools, only : is_master, ncpu, ncpus
    use problem , only : init_domain, init_problem, check_ref

    implicit none

! local variables
!
    type(block), pointer :: pblock, pparent, pchild, pneigh, pnext
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
    if (is_master()) then
      write(*,"(1x,a)"      ) "Generating initial mesh:"
      write(*,"(4x,a,1x,i6)") "refining to max. level =", maxlev
      write(*,"(4x,a,1x,i6)") "effective resolution   =", ncells*2**maxlev
    endif

! initialize blocks
!
    call init_blocks

! allocate the initial structure of blocks according to the problem
!
    call init_domain

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out the coarse blocks with the initial condition
!
!     pblock => plist
!     do while (associated(pblock))
!
! ! set level
! !
!       pblock%level    = 1
!       pblock%refine   = 0
!
! ! set initial conditions
! !
!       call init_problem(pblock)
!
! ! assign pointer to the next block
! !
!       pblock => pblock%next
!
!     end do

! TODO: refine blocks on master untill the total number of blocks exceeds
!       the number of MPI processes, then interrupt refining and autobalance
!       the allocated blocks, after that continue refining on all processes
!       autobalancing after each level refinement
!

! at this point the inital blocks are allocated and set for refinement,
! so iterate over all levels from 1 to maxlevel and create sub-blocks,
! set the initial conditions for each, check criterium and set for
! refinement according to the criterium fullfilment
!
! TODO: refine blocks on master untill the total number of blocks exceeds
!       the number of MPI processes, then interrupt refining and autobalance
!       the allocated blocks, after that continue refining on all processes
!       autobalancing after each level of refinement
!
! refine the blocks until the number of blocks is smaller than the number
! of processes
!
    if (is_master()) &
      write(*,"(4x,a,$)") "refining level         =    "

    l = 1
    do while (l .le. maxlev)

! print the level processed
!
      if (is_master()) &
        write(*,"(1x,i2,$)") l

! iterate over all blocks at the current level and initialize problem
!
      pblock => plist
      do while (associated(pblock))

! set initial conditions if the block at current level
!
        if (pblock%level .le. l) &
          call init_problem(pblock)

! assign pointer to the next block
!
        pblock => pblock%next

      end do

! for the maximum level we only initialize problem (with not check of the
! refinements criterion or the refinement)
!
      if (l .lt. maxlev) then

! iterate over all blocks at the current level and check the refinement
! criterion; do not allow for derefinement
!
        pblock => plist
        do while (associated(pblock))

          if (pblock%level .eq. l) &
            pblock%refine = max(0, check_ref(pblock))

! assign pointer to the next block
!
          pblock => pblock%next

        end do

! walk through all levels down and check if the neighbors have to be refined
! too; there is no need for checking the blocks at the lowest level;
!
        do n = l, 2, -1

! iterate over all blocks of the current level and if the current block is
! selected for the refinement and its neighbors are at lower levels select them
! for refinement too;
!
          pblock => plist
          do while (associated(pblock))

! check if the current block is at the current level, a leaf, and selected for
! refimenet
!
            if (pblock%level .eq. n) then
              if (pblock%leaf) then
                if (pblock%refine .eq. 1) then

! iterate over all neighbors
!
                  do i = 1, ndims
                    do j = 1, 2
                      do k = 1, 2

                        pneigh => get_pointer(pblock%neigh(i,j,k)%id)

! check if the neighbor is associated
!
                        if (associated(pneigh)) then
                          if (pneigh%leaf) then


! if th eneighbor has lower level, select it to be refined too
!
                            if (pneigh%level .lt. pblock%level) &
                              pneigh%refine = 1
                          else
                            call print_error("mesh::init_mesh", "Neighbor is not a leaf!")
                          endif
                        endif

                      end do
                    end do
                  end do

                endif
              endif
            endif

! assign pointer to the next block
!
            pblock => pblock%next

          end do
        end do

! walk through the levels starting from the lowest to the current
!
        do n = 1, l

! iterate over all blocks and refine selected
!
          pblock => plist
          do while (associated(pblock))

! check if the current block is at the current level and selected for refinement
! and perform the refinement
!
            if (pblock%level .eq. n .and. pblock%refine .eq. 1) &
              call refine_block(pblock)

! assign pointer to the next block
!
            pblock => pblock%next

          end do

        end do

      endif

      l = l + 1
    end do

#ifdef MPI
! divide all blocks between all processes
!
    l = 0
    pblock => plist
    do while (associated(pblock))

! assign the cpu to the current block
!
      pblock%cpu    = l * ncpus / nblocks

! assign pointer to the next block
!
      pblock => pblock%next

      l = l + 1
    end do

! update the cpu field of the neighbors, parent and children
!
    pblock => plist
    do while (associated(pblock))

! update neighbors
!
      do i = 1, ndims
        do j = 1, 2
          do k = 1, 2

            pneigh => get_pointer(pblock%neigh(i,j,k)%id)

            if (associated(pneigh)) &
              pblock%neigh(i,j,k)%cpu = pneigh%cpu

          end do
        end do
      end do

! update parent
!
      pparent => get_pointer(pblock%parent%id)
      if (associated(pparent)) &
        pblock%parent%cpu = pparent%cpu

! update children
!
      do p = 1, nchild
        pchild => get_pointer(pblock%child(p)%id)

        if (associated(pchild)) &
          pblock%child(p)%cpu = pchild%cpu
      end do

! assign pointer to the next block
!
      pblock => pblock%next
    end do

! remove all blocks which don't belong to the current process
!
    pblock => plist
    do while (associated(pblock))
      pnext => pblock%next

      if (pblock%cpu .ne. ncpu) &
        call deallocate_block(pblock)

      pblock => pnext
    end do
#endif /* MPI */

! print information
!
    if (is_master()) then
      write(bstr,"(i)") nblocks
      write(tstr,"(i)") (2**maxlev)**ndims
      write(*,*)
      write(*,"(4x,a,1x,a6,' / ',a,' = ',f8.4,' %')") "allocated/total blocks =", trim(adjustl(bstr)),trim(adjustl(tstr)), (100.0*last_id)/(2**maxlev)**ndims

      write(bstr,"(i)") nleafs
      write(tstr,"(i)") (2**maxlev)**ndims
      write(*,"(4x,a,1x,a6,' / ',a,' = ',f8.4,' %')") "leafs/cover blocks     =", trim(adjustl(bstr)),trim(adjustl(tstr)), (100.0*nleafs)/(2**maxlev)**ndims
    endif

! allocating space for coordinate variables
!
    allocate(ax  (maxlev, im))
    allocate(ay  (maxlev, jm))
    allocate(az  (maxlev, km))
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

    use config  , only : maxlev
    use blocks  , only : block, plist, ndims, nchild, nblocks, refine_block    &
                       , get_pointer, maxid, last_id
    use error   , only : print_info
#ifdef MPI
    use mpitools, only : ncpus, ncpu, mfindmaxi, mallreducesuml, mallreducemaxl, msendi, mrecvi
#endif /* MPI */
    use problem , only : check_ref

    implicit none

    integer, intent(in) :: ref

! local variables
!
    type(block), pointer :: pblock, pneigh, pchild, pparent
    integer(kind=4)      :: l, i, j, k, m, n, p, q

#ifdef MPI
    logical              :: flag
    integer(kind=4)      :: mblk

! local arrays for MPI
!
    integer(kind=4), dimension(:,:,:)    , allocatable :: iblk, idif, ichl
#endif /* MPI */

!-------------------------------------------------------------------------------
!
#ifdef MPI
! find maximum number of blocks
!
    mblk = last_id
    call mfindmaxi(mblk)

! allocate and initialize temporary arrays
!
    allocate(iblk(0:ncpus-1, mblk, 3))
    allocate(idif(0:ncpus-1, mblk, 2))
    allocate(ichl(0:ncpus-1, mblk, 4))
    iblk(:,:,:) = 0
    idif(:,:,:) = 0
    ichl(:,:,:) = 0

! check the refinement criterion for all leafes
!
    pblock => plist
    do while (associated(pblock))

      if (pblock%leaf) then
        pblock%refine = check_ref(pblock)

        if (pblock%level .eq. maxlev) &
          pblock%refine = min(0, pblock%refine)
      endif

      pblock => pblock%next
    end do

! walk down from the highest level and select neighbors for refinement if they
! are at lower levels
!
    do l = maxlev, 1, -1

! iterate over all blocks at the current level and if the current block is
! selected for refinement set the neighbor at lower level to be refined too
!
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l .and. pblock%leaf) then
          if (pblock%refine .eq. 1) then

            do i = 1, ndims
              do j = 1, 2
                do k = 1, 2
                  if (pblock%neigh(i,j,k)%cpu .eq. ncpu) then
                    pneigh => get_pointer(pblock%neigh(i,j,k)%id)

                    if (pneigh%level .lt. pblock%level) &
                      pneigh%refine = 1

                    if (pneigh%level .eq. pblock%level .and. pneigh%refine .eq. -1) &
                      pneigh%refine = 0
                  endif
                end do
              end do
            end do

          endif
        endif

        pblock => pblock%next
      end do

! check if all children are selected for derefinement; if this condition is not
! fullfiled deselect all children from derefinement
!
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l-1 .and. .not. pblock%leaf) then
          flag = .true.
          do p = 1, nchild
            if (pblock%child(p)%cpu .eq. ncpu) then
              pchild => get_pointer(pblock%child(p)%id)

              flag = flag .and. (pchild%refine .eq. -1)
            endif
          end do
          if (.not. flag) then
            do p = 1, nchild
              if (pblock%child(p)%cpu .eq. ncpu) then
                pchild => get_pointer(pblock%child(p)%id)

                if (pchild%leaf) &
                  pchild%refine = max(0, pchild%refine)
              endif
            end do
          endif
        endif

        pblock => pblock%next
      end do

! deselect neighbors from derefinement if the current block is set for
! refinement or it is at higher level and not selected for derefinement
!
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l .and. pblock%leaf) then

          if (pblock%refine .ne. -1) then

            do i = 1, ndims
              do j = 1, 2
                do k = 1, 2
                  if (pblock%neigh(i,j,k)%cpu .eq. ncpu) then
                    pneigh => get_pointer(pblock%neigh(i,j,k)%id)

                    if (pneigh%refine .eq. -1) then

                      if (pneigh%level .lt. pblock%level) &
                        pneigh%refine = 0

                    endif

                  endif
                end do
              end do
            end do

          endif

        endif

        pblock => pblock%next
      end do

    end do

! ! check blocks at the maximum level for derefinement
! !
!     pblock => plist
!     do while (associated(pblock))
!       if (pblock%leaf) then
!         pblock%refine = check_ref(pblock)
!
!         if (pblock%level .eq. maxlev) &
!           pblock%refine = min(0, pblock%refine)
!
!         if (pblock%level .eq. 1) &
!           pblock%refine = max(0, pblock%refine)
!
!         iblk(ncpu,pblock%id,1) = pblock%refine
!         iblk(ncpu,pblock%id,2) = pblock%level
!         iblk(ncpu,pblock%id,3) = 1
!       endif
!       pblock => pblock%next
!     end do
!
! ! exchange information about the block refinements
! !
!     call mallreducesuml(size(iblk), iblk)
!
! ! go down to lower levels
! !
!     do l = maxlev - 1, 1, -1
!
!       idif(:,:,:) = 0
!
!       pblock => plist
!       do while (associated(pblock))
!         if (pblock%leaf) then
!           if (pblock%level .eq. l) then
!             if (pblock%refine .eq. -1) then
!               do i = 1, ndims
!                 do j = 1, 2
!                   do k = 1, 2
!                     if (pblock%neigh(i,j,k)%id .gt. 0) then
!                       if (iblk(pblock%neigh(i,j,k)%cpu,pblock%neigh(i,j,k)%id,2) .gt. pblock%level) then
!                         idif(ncpu,pblock%id,1) = 1
!                       endif
!                     endif
!                   end do
!                 end do
!               end do
!             endif
!
!             if (pblock%refine .eq. 1) then
!               do i = 1, ndims
!                 do j = 1, 2
!                   do k = 1, 2
!                     if (pblock%neigh(i,j,k)%id .gt. 0) then
!                       if (iblk(pblock%neigh(i,j,k)%cpu,pblock%neigh(i,j,k)%id,2) .lt. pblock%level) then
!                         idif(pblock%neigh(i,j,k)%cpu,pblock%neigh(i,j,k)%id,2) = 1
!                       endif
!                     endif
!                   end do
!                 end do
!               end do
!             endif
!           endif
!         endif
!         pblock => pblock%next
!       end do
!
!       pblock => plist
!       do while (associated(pblock))
!         if (.not. pblock%leaf) then
!           flag = .true.
!           do p = 1, nchild
!             if (pblock%child(p)%id .gt. 0) &
!               flag = flag .and. (iblk(pblock%child(p)%cpu,pblock%child(p)%id,1) .eq. -1) .and. (iblk(pblock%child(p)%cpu,pblock%child(p)%id,3) .eq. 1)
!           end do
!           if (.not. flag) then
!             do p = 1, nchild
!               if (pblock%child(p)%id .gt. 0) &
!                 idif(pblock%child(p)%cpu,pblock%child(p)%id,1) = 1
!             end do
!           endif
!         endif
!         pblock => pblock%next
!       end do
!
!       call mallreducemaxl(size(idif), idif)
!
!       do i = 0, ncpus - 1
!         do j = 1, mblk
!           if (idif(i,j,1) .eq. 1) &
!             iblk(i,j,1) = max(0, iblk(i,j,1))
!           if (idif(i,j,2) .eq. 1) &
!             iblk(i,j,1) = 1
!         end do
!       end do
!
!       do i = 1, mblk
!         if (iblk(ncpu,i,2) .gt. 0) then
!           pblock => get_pointer(i)
!           pblock%refine = iblk(ncpu,i,1)
!         endif
!       end do
!
!     end do

! perform the actual derefinement
!
    do l = maxlev, 2, -1

      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l) then
          if (pblock%leaf) then
            if (pblock%refine .eq. -1) then
              pparent => get_pointer(pblock%parent%id)

              if (associated(pparent)) then
                call restrict_block(pparent)
                pblock => pparent
              endif
            endif
          endif
        endif

        pblock => pblock%next
      end do

    end do

! perform the actual refinement
!
    do l = 1, maxlev - 1

      pblock => plist
      do while (associated(pblock))

        if (pblock%leaf) then
          if (pblock%level .eq. l) then
            if (pblock%refine .eq. 1) then
              pparent => pblock

              call refine_block(pblock)
              call prolong_block(pparent)
            endif
          endif
        endif

        pblock => pblock%next

      end do

    end do


! ! perform the actual refinement
! !
!     do l = 1, 3!maxlev - 1
!       pblock => plist
!       do while (associated(pblock))
!         if (pblock%leaf) then
!           if (pblock%level .eq. l) then
!             if (pblock%refine .eq. 1) then
!               pparent => pblock
!
!               call refine_block(pblock)
!               call prolong_block(pparent)
!
!               do p = 1, nchild
!                 ichl(ncpu,pparent%id,p) = pparent%child(p)%id
!               end do
!             endif
!           endif
!         endif
!         pblock => pblock%next
!       end do
!     end do

! deallocate temporary arrays
!
    deallocate(iblk)
    deallocate(idif)
    deallocate(ichl)
#else /* MPI */
#endif /* MPI */

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

    use blocks       , only : block, nv => nvars, get_pointer
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
    pchild => get_pointer(pblock%child(1)%id)
    do k = 1, km
      pchild%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = in + 1
    iu = il + im - 1
    jl = 1
    ju = jl + jm - 1
    pchild => get_pointer(pblock%child(2)%id)
    do k = 1, km
      pchild%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = 1
    iu = il + im - 1
    jl = in + 1
    ju = jl + jm - 1
    pchild => get_pointer(pblock%child(3)%id)
    do k = 1, km
      pchild%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
    end do

    il = in + 1
    iu = il + im - 1
    jl = in + 1
    ju = jl + jm - 1
    pchild => get_pointer(pblock%child(4)%id)
    do k = 1, km
      pchild%u(:,:,:,k) = u(:,il:iu,jl:ju,k)
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

    use blocks       , only : block, nv => nvars, derefine_block, get_pointer
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
    pbl => get_pointer(pblock%child(1)%id)
    pbr => get_pointer(pblock%child(2)%id)
    ptl => get_pointer(pblock%child(3)%id)
    ptr => get_pointer(pblock%child(4)%id)

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
