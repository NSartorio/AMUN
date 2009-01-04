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

    use config  , only : im, jm, km, xmin, xmax, ymin, ymax, zmin, zmax        &
                       , ncells, maxlev
    use blocks  , only : list_allocated, init_blocks, clear_blocks             &
                       , deallocate_block, refine_block, get_pointer           &
                       , block, nchild, ndims, plist, last_id, nblocks
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
    do while (l .lt. maxlev)

! print the level processed
!
      if (is_master()) &
        write(*,"(1x,i2,$)") l

! iterate over all blocks at the current level and check the refinement
! criterion
!
      pblock => plist
      do while (associated(pblock))

! check refinements criterion for the current block;
! do not allow for derefinement initially;
!
        if (pblock%leaf) then
          if (pblock%level .le. l) then
            pblock%refine = max(0, check_ref(pblock))
          else
            pblock%refine = 0
          endif
        endif

! assign pointer to the next block
!
        pblock => pblock%next

      end do

! iterate over all blocks of the current level and lower levels and select
! the neighbors of refined blocks to be refine too
!
      do n = l, 2, -1
        pblock => plist
        do while (associated(pblock))

! check if the current block is a leaf, at the current level and selected
! to be refined
!
          if (pblock%leaf) then
            if (pblock%level .eq. n) then
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
                          call print_error("mesh::init_mesh", "Block is not a leaf!")
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

! iterate over all blocks and refine selected ones
!
      do n = 1, l
        pblock => plist
        do while (associated(pblock))

! check if the current block is a leaf, is at the current level and needs
! to be refined
!
          if (pblock%leaf) then
            if (pblock%level .eq. n) then
              if (pblock%refine .eq. 1) then

! refine this block
!
                call refine_block(pblock)

! iterate over all children
!
                pparent => get_pointer(pblock%parent%id)
                do p = 1, nchild
                  pchild => get_pointer(pparent%child(p)%id)

! initialize problem for that child
!
                  if (associated(pchild)) &
                    call init_problem(pchild)
                end do
              endif
            endif
          endif

! assign pointer to the next block
!
          pblock => pblock%next

        end do
      end do

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

! update the cpu field of the neighbors
!
    pblock => plist
    do while (associated(pblock))

      do i = 1, ndims
        do j = 1, 2
          do k = 1, 2

            pneigh => get_pointer(pblock%neigh(i,j,k)%id)

            if (associated(pneigh)) &
              pblock%neigh(i,j,k)%cpu = pneigh%cpu

          end do
        end do
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
      write(bstr,"(i)") last_id
      write(tstr,"(i)") (2**maxlev)**ndims
      write(*,*)
      write(*,"(4x,a,1x,a6,' / ',a,' = ',f8.4,' %')") "allocated/total blocks =", trim(adjustl(bstr)),trim(adjustl(tstr)), (100.0*last_id)/(2**maxlev)**ndims
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
                       , get_pointer
    use error   , only : print_info
#ifdef MPI
    use mpitools, only : ncpus, ncpu, mallreducemaxl, msendi, mrecvi
#endif /* MPI */
    use problem , only : check_ref

    implicit none

    integer, intent(in) :: ref

! local variables
!
    type(block), pointer :: pblock, pneigh, pchild, pparent
    integer(kind=4)      :: l, i, j, k, m, n, p, q

#ifdef MPI
    integer(kind=4)      :: itag

! local arrays for MPI
!
    integer(kind=4), dimension(:,:,:)    , allocatable :: iblk, ichl
    integer(kind=4), dimension(:,:)      , allocatable :: cn, cc, ibuf
#endif /* MPI */

!----------------------------------------------------------------------
!
#ifdef MPI
! allocate temporary arrays; since we have two blocks per boundary and 4 sides
! of each block we need to increase the second dimension
!
    allocate(cn  (0:ncpus-1,0:ncpus-1))
    allocate(cc  (0:ncpus-1,0:ncpus-1))
    allocate(iblk(0:ncpus-1,2**(NDIMS-1)*NDIMS*nblocks,3))
    allocate(ichl(0:ncpus-1,2**(NDIMS-1)*NDIMS*nblocks,7))
#endif /* MPI */

    do l = 1, maxlev

#ifdef MPI
! reset the local arrays storing blocks to exchange
!
      cn(:,:)     = 0
      iblk(:,:,:) = 0
#endif /* MPI */

! check refinement criterion
!
      pblock => plist
      do while (associated(pblock))

        if (pblock%level .eq. l) then
          pblock%refine = 0

          if (pblock%leaf) then
            pblock%refine = check_ref(pblock)

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
! - all neighbors must be at the same or one level higher
!
      do n = l, 2, -1
        pblock => plist
        do while (associated(pblock))

          if (pblock%level .eq. n) then
            if (pblock%leaf .and. pblock%refine .eq. 1) then
              do i = 1, ndims
                do j = 1, 2
                  do k = 1, 2
#ifdef MPI
                    if (pblock%neigh(i,j,k)%cpu .eq. ncpu) then
#endif /* MPI */
                      pneigh => get_pointer(pblock%neigh(i,j,k)%id)
                      if (associated(pneigh)) then
                        if (pneigh%level .lt. pblock%level) &
                          pneigh%refine = 1
                      endif
#ifdef MPI
                    else

! TODO: fill temporary array with neighbors of the current block laying on other processes
!       this array should store current block ID and level, and neighbor ID
!
! get the processor number of neighbor
!
                      p = pblock%neigh(i,j,k)%cpu

! increase the number of blocks to retrieve from that CPU
!
                      cn(ncpu,p) = cn(ncpu,p) + 1

! fill out the info array
!
                      iblk(p,cn(ncpu,p),1) = pblock%id
                      iblk(p,cn(ncpu,p),2) = pblock%level
                      iblk(p,cn(ncpu,p),3) = pblock%neigh(i,j,k)%id

                    endif
#endif /* MPI */

                  end do
                end do
              end do
            endif
          endif

          pblock => pblock%next
        end do

#ifdef MPI
! TODO: after sweeping the current level, send the temporary array to neighbors
!       and select them for refinement if their level is lower
!
! update number of blocks across all processes
!
        call mallreducemaxl(size(cn),cn)

        if (maxval(cn) .gt. 0) then

          allocate(ibuf(maxval(cn),3))

          do p = 0, ncpus-1
            do q = 0, ncpus-1

              itag = 10*(p * ncpus + q) + 2111

              if (cn(p,q) .gt. 0) then
                if (ncpu .eq. p) then  ! sender
                  call msendi(size(iblk(q,1:cn(p,q),:)), q, itag, iblk(q,1:cn(p,q),:))
                endif
                if (ncpu .eq. q) then  ! receiver
                  call mrecvi(size(ibuf(1:cn(p,q),:)), p, itag, ibuf(1:cn(p,q),:))

                  do i = 1, cn(p,q)
                    pblock => get_pointer(ibuf(i,3))

                   if (pblock%level .lt. ibuf(i,2)) &
                     pblock%refine = 1
                  end do
                endif
              endif
            end do
          end do

          deallocate(ibuf)
        endif
#endif /* MPI */
      end do

! after selecting blocks for refinement, do refine them gradually starting from
! the lowest level in order to avoid situation when two neighbors have level
! difference larger than 2
!
      do n = 1, l
#ifdef MPI
        cc  (:,:)   = 0
        ichl(:,:,:) = 0
#endif /* MPI */

        pblock => plist
        do while (associated(pblock))

          if (pblock%level .eq. n) then
            if (pblock%leaf .and. pblock%refine .eq. 1) then

              pparent => pblock

              call refine_block(pblock)
              call prolong_block(pparent)

#ifdef MPI
! TODO: collect information about refined block in order to update neighbors
!       laying on other processors

              do i = 1, ndims
                do j = 1, 2
                  do k = 1, 2
                    if (pparent%neigh(i,j,k)%cpu .ne. ncpu) then
                      p = pparent%neigh(i,j,k)%cpu
                      cc(ncpu,p) = cc(ncpu,p) + 1
                      ichl(p,cc(ncpu,p),1) = pparent%neigh(i,j,k)%id
                      ichl(p,cc(ncpu,p),2) = pparent%id
                      ichl(p,cc(ncpu,p),3) = pparent%level
                      ichl(p,cc(ncpu,p),4) = pparent%child(1)%id
                      ichl(p,cc(ncpu,p),5) = pparent%child(2)%id
                      ichl(p,cc(ncpu,p),6) = pparent%child(3)%id
                      ichl(p,cc(ncpu,p),7) = pparent%child(4)%id
                    endif
                  end do
                end do
              end do
#endif /* MPI */

! TODO: remove parent after refinement
!
            endif
          endif

          pblock => pblock%next

        end do

#ifdef MPI
! TODO: exchange information about refined blocks in order to update neighbors
!
        call mallreducemaxl(size(cc),cc)

        if (maxval(cc) .gt. 0) then

          allocate(ibuf(maxval(cc),7))

          do p = 0, ncpus-1
            do q = 0, ncpus-1

              itag = 10*(p * ncpus + q) + 3111

              if (cc(p,q) .gt. 0) then
                if (ncpu .eq. p) then  ! sender
                  call msendi(size(ichl(q,1:cc(p,q),:)), q, itag, ichl(q,1:cc(p,q),:))
                endif
                if (ncpu .eq. q) then  ! receiver
                  call mrecvi(size(ibuf(1:cc(p,q),:)), p, itag, ibuf(1:cc(p,q),:))

                  do m = 1, cn(p,q)
                    pblock => get_pointer(ibuf(m,1))

                    if (pblock%leaf) then
! TODO: if leaf, update its neighbors according to the information sent
!
                      do i = 1, ndims
                        do j = 1, 2
                          do k = 1, 2
                            if (pblock%neigh(i,j,k)%cpu .eq. p) then
                              if (pblock%neigh(i,j,k)%id .eq. ibuf(m,2)) then

! TODO: depending in the level difference and side update neighbors
!
!                                 if (pblock%level .eq. ibuf(m,3) then
!                                 endif
!                                 pblock%neigh(i,j,k)%id =
                              endif
                            endif
                          end do
                        end do
                      end do

                    else
! TODO: if not leaf, it means that the block has been refined in the meantime, so
!       update its children
!
                    endif
                  end do
                endif
              endif
            end do
          end do

          deallocate(ibuf)
        endif
#endif /* MPI */
      end do

    end do

! ! iterate starting from the maximum level and select blocks for derefinement
! ! if they fullfil below conditions
! !
!     do l = maxlev, 1, -1
!
! ! derefinement conditions for blocks:
! !
! ! - all neighbors must be at the same or lower level
! ! - all neighbors at the same level must be not selected for refinement
! ! - all sibling must be at the same level and marked for derefinement
! !
! ! if all above conditions are fulfilled, select the refinement of
! ! the parent of the current block to -1
! !
! ! if at least one is not fulfilled, set all sibling to not do the refinement
! !
!       do n = l, maxlev
!       pblock => plist
!       do while (associated(pblock))
!
!         if (pblock%level .eq. n) then
!           if (pblock%leaf .and. pblock%refine .eq. -1) then
!             do i = 1, ndims
!               do j = 1, 2
!                 do k = 1, 2
!                   pneigh => get_pointer(pblock%neigh(i,j,k)%id)
!                   if (associated(pneigh)) then
!                     if (pneigh%level .gt. pblock%level) &
!                       pblock%refine = 0
!                     if (pneigh%level .eq. pblock%level .and. pneigh%refine .eq. 1) &
!                       pblock%refine = 0
!                   endif
!
!                 end do
!               end do
!             end do
!
!             pparent => get_pointer(pblock%parent%id)
!             if (associated(pparent)) then
!               do p = 1, nchild
!                 pchild => get_pointer(pparent%child(p)%id)
!
!                 if (associated(pchild)) then
!                   if (pchild%refine .ne. -1 .or. pchild%level .ne. pblock%level) &
!                     pblock%refine = 0
!                 endif
!               end do
!
!               if (pblock%refine .ne. -1) then
!                 do p = 1, nchild
!                   pchild => get_pointer(pparent%child(p)%id)
!
!                   if (associated(pchild)) then
!                     if (pchild%refine .eq. -1 .and. pchild%level .eq. pblock%level) &
!                       pchild%refine = 0
!                   endif
!                 end do
!               endif
!             endif
!
!           endif
!         endif
!
!         pparent => get_pointer(pblock%parent%id)
!         if (associated(pparent) .and. pblock%refine .eq. -1) &
!           pparent%refine = -1
!
! ! point to the next block
! !
!         pblock => pblock%next
!       end do
!       end do
!
! ! perform derefinement for all selected blocks
! !
!       do n = maxlev, l, -1
!       pblock => plist
!       do while (associated(pblock))
!
!         if (pblock%level .eq. n) then
!           if (pblock%leaf .and. pblock%refine .eq. -1) then
!
!             pparent => get_pointer(pblock%parent%id)
!
!             if (associated(pparent) .and. pparent%refine .eq. -1) then
!               call restrict_block(pparent)
!               pblock => pparent
!             endif
!           endif
!         endif
!
!         pblock => pblock%next
!       end do
!       end do
!
!     end do

#ifdef MPI
! deallocate temporary arrays
!
    deallocate(cc)
    deallocate(ichl)
    deallocate(cn)
    deallocate(iblk)
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
