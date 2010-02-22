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
                       , ncells, maxlev, rdims
    use blocks  , only : block_meta, block_data, list_meta, list_data          &
                       , init_blocks, clear_blocks, refine_block               &
                       , deallocate_datablock, nblocks, nleafs, dblocks        &
                       , nchild, ndims, nsides, nfaces
    use error   , only : print_info, print_error
    use mpitools, only : is_master, ncpu, ncpus
    use problem , only : init_domain, init_problem, check_ref

    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta_block, pneigh, pnext
    type(block_data), pointer :: pdata_block

! local variables
!
    integer(kind=4)      :: l, p, i, j, k, n
    character(len=64)    :: fmt
    character(len=32)    :: bstr, tstr

!----------------------------------------------------------------------
!
! check if the list is allocated, if yes deallocate it
!
    if (associated(list_meta)) then
      call print_info("mesh::init_mesh", "Block list is allocated, deallocate it!")

      call clear_blocks
    end if

! initialize blocks
!
    call init_blocks

! allocate the initial structure of blocks according to the problem
!
    call init_domain

! print general information about resolutions
!
    if (is_master()) then
      write(*,"(1x,a)"         ) "Generating initial mesh:"
      write(*,"(4x,a,  1x,i6)" ) "refining to max. level  =", maxlev
      write(*,"(4x,a,3(1x,i6))") "lowest level resolution =", rdims(1:ndims) * ncells
      write(*,"(4x,a,3(1x,i6))") "effective resolution    =", rdims(1:ndims) * ncells * 2**(maxlev - 1)
    end if

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out the coarse blocks with the initial condition
!
    if (is_master()) &
      write(*,"(4x,a,$)") "refining level          =    "

    l = 1
    do while (l .le. maxlev)

! print the level currently processed
!
      if (is_master()) &
        write(*,"(1x,i2,$)") l

! iterate over all data blocks at the current level and initialize the problem
!
      pdata_block => list_data
      do while (associated(pdata_block))

! set the initial conditions at the current block
!
        if (pdata_block%meta%level .le. l) &
          call init_problem(pdata_block)

! assign pointer to the next block
!
        pdata_block => pdata_block%next
      end do

! at the maximum level we only initialize the problem (without checking the
! refinement criterion)
!
      if (l .lt. maxlev) then

! iterate over all data blocks at the current level and check the refinement
! criterion; do not allow for derefinement
!
        pdata_block => list_data
        do while (associated(pdata_block))

          if (pdata_block%meta%level .eq. l) then
            pdata_block%meta%refine = max(0, check_ref(pdata_block))

! if there is only one block, and it is set not to be refined, refine it anyway
! because the resolution for the problem initiation may be too small
!
            if (nblocks .eq. 1 .and. l .eq. 1) &
              pdata_block%meta%refine = 1
          end if

! assign pointer to the next block
!
          pdata_block => pdata_block%next
        end do

! walk through all levels down from the current level and check if select all
! neighbors for the refinement if they are at lower level; there is no need for
! checking the blocks at the lowest level;
!
        do n = l, 2, -1

! iterate over all meta blocks at the level n and if the current block is
! selected for the refinement and its neighbors are at lower levels select them
! for refinement too;
!
          pmeta_block => list_meta
          do while (associated(pmeta_block))

! check if the current block is at the level n, is a leaf, and selected for
! refinement
!
            if (pmeta_block%level .eq. n) then
              if (pmeta_block%leaf) then
                if (pmeta_block%refine .eq. 1) then

! iterate over all neighbors
!
                  do i = 1, ndims
                    do j = 1, nsides
                      do k = 1, nfaces

! assign pointer to the neighbor
!
                        pneigh => pmeta_block%neigh(i,j,k)%ptr

! check if the neighbor is associated
!
                        if (associated(pneigh)) then

! check if the neighbor is a leaf, if not something wrong is going on
!
                          if (pneigh%leaf) then

! if the neighbor has lower level, select it to be refined too
!
                            if (pneigh%level .lt. pmeta_block%level) &
                              pneigh%refine = 1

                          else
                            call print_error("mesh::init_mesh", "Neighbor is not a leaf!")
                          end if
                        end if

                      end do
                    end do
                  end do

                end if
              end if
            end if

! assign pointer to the next block
!
            pmeta_block => pmeta_block%next

          end do
        end do

!! refine all selected blocks starting from the lowest level
!!
! walk through the levels starting from the lowest to the current level
!
        do n = 1, l

! iterate over all meta blocks
!
          pmeta_block => list_meta
          do while (associated(pmeta_block))

! check if the current block is at the level n and, if it is selected for
! refinement and if so, perform the refinement on this block
!
            if (pmeta_block%level .eq. n .and. pmeta_block%refine .eq. 1) then

! perform the refinement
!
              call refine_block(pmeta_block, .true.)

            end if

! assign pointer to the next block
!
            pmeta_block => pmeta_block%next
          end do
        end do
      end if

      l = l + 1
    end do

! deallocate data blocks of non leafs
!
    pmeta_block => list_meta
    do while (associated(pmeta_block))

      if (.not. pmeta_block%leaf) &
        call deallocate_datablock(pmeta_block%data)

! assign pointer to the next block
!
      pmeta_block => pmeta_block%next
    end do

#ifdef MPI
! divide blocks between all processes, use the number of data blocks to do this,
! but keep blocks from the top level which have the same parent packed together
!
    l = nleafs / ncpus
    n = 0
    p = 0

    pmeta_block => list_meta
    do while (associated(pmeta_block))

! assign the cpu to the current block
!
      pmeta_block%cpu = n

! increase the number of blocks on the current process; if it exceeds the
! allowed number reset the counter and increase the processor number
!
      if (pmeta_block%leaf) then
        p = p + 1
        if (p .ge. l) then
          n = min(ncpus - 1, n + 1)
          p = 0
        end if
      end if

! assign pointer to the next block
!
      pmeta_block => pmeta_block%next
    end do

! remove all data blocks which do not belong to the current process
!
    pmeta_block => list_meta
    do while (associated(pmeta_block))
      pnext => pmeta_block%next

! if the current block belongs to another process and its data field is
! associated, deallocate its data field
!
      if (pmeta_block%cpu .ne. ncpu .and. associated(pmeta_block%data)) &
        call deallocate_datablock(pmeta_block%data)

! assign pointer to the next block
!
      pmeta_block => pnext
    end do
#endif /* MPI */

! print information about the generated geometry
!
    if (is_master()) then
      p = 0
      do l = 0, maxlev - 1
        k = 2**(ndims * l)
        p = p + k
      end do
      p = p * rdims(1) * rdims(2) * rdims(3)
      k = k * rdims(1) * rdims(2) * rdims(3)

      i = nint(alog10(1.0*nblocks + 1)) + 1
      j = nint(alog10(1.0*p + 1)) + 1

      write(fmt, "(a,i1,a,i1,a)") "(4x,a,1x,i", i, ",' / ',i", j, ",' = ',f8.4,' %')"

      write(*,*)
      write(*,fmt) "leafs    /cover blocks  =", nleafs , k, (100.0 * nleafs ) / k
      write(*,fmt) "allocated/total blocks  =", nblocks, p, (100.0 * nblocks) / p
    end if

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
      p = ncells * 2**(l - 1)

      adx (l) = (xmax - xmin) / (rdims(1) * p)
      adxi(l) = 1.0 / adx(l)
      ady (l) = (ymax - ymin) / (rdims(2) * p)
      adyi(l) = 1.0 / ady(l)
#if NDIMS == 3
      adz (l) = (zmax - zmin) / (rdims(3) * p)
#else
      adz (l) = 1.0
#endif /* NDIMS == 3 */
      adzi(l) = 1.0 / adz(l)
    end do

! get the minimum grid step
!
    dx_min = 0.5 * min(adx(maxlev), ady(maxlev), adz(maxlev))

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
  subroutine update_mesh

    use config  , only : maxlev, im, jm, km
    use blocks  , only : block_meta, block_data, list_meta, list_data          &
                       , nleafs, dblocks, nchild, ndims, nsides, nfaces        &
                       , refine_block, derefine_block, append_datablock        &
                       , associate_blocks, deallocate_datablock, nv => nvars
    use error   , only : print_info
#ifdef MPI
    use mpitools, only : ncpus, ncpu, is_master, mallreducesuml, msendf, mrecvf
#endif /* MPI */
    use problem , only : check_ref

    implicit none

! local variables
!
    logical         :: flag
    integer(kind=4) :: i, j, k, l, n, p

#ifdef MPI
! tag for data exchange
!
    integer(kind=4)                         :: itag

! array for the update of the refinement flag on all processors
!
    integer(kind=4), dimension(nleafs)      :: ibuf

! local buffer for data block exchange
!
    real(kind=8)   , dimension(nv,im,jm,km) :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh, pchild, pparent
    type(block_data), pointer :: pdata

!-------------------------------------------------------------------------------
!
! check the refinement criterion for all leafs
!
    pdata => list_data
    do while (associated(pdata))
      pmeta => pdata%meta

      if (pmeta%leaf) then
        pmeta%refine = check_ref(pdata)

        if (pmeta%level .eq. maxlev) &
          pmeta%refine = min(0, pmeta%refine)
      end if

      pdata => pdata%next
    end do

#ifdef MPI
! store refinement flags for all blocks for exchange between processors
!
    l = 1
    pmeta => list_meta
    do while (associated(pmeta))
      if (pmeta%leaf) then
        ibuf(l) = pmeta%refine
        l = l + 1
      end if
      pmeta => pmeta%next
    end do

! update information across all processors
!
    call mallreducesuml(nleafs, ibuf)

! update non-local block refinement flags
!
    l = 1
    pmeta => list_meta
    do while (associated(pmeta))
      if (pmeta%leaf) then
        if (pmeta%cpu .ne. ncpu) &
          pmeta%refine = ibuf(l)
        l = l + 1
      end if
      pmeta => pmeta%next
    end do
#endif /* MPI */

! walk down from the highest level and select neighbors for refinement if they
! are at lower levels
!
    do l = maxlev, 1, -1

! iterate over all blocks at the current level and if the current block is
! selected for refinement set the neighbor at lower level to be refined too
!
      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%level .eq. l .and. pmeta%leaf) then
          if (pmeta%refine .eq. 1) then

            do i = 1, ndims
              do j = 1, nsides
                do k = 1, nfaces
                  pneigh => pmeta%neigh(i,j,k)%ptr
                  if (associated(pneigh)) then
                    if (pneigh%level .lt. pmeta%level) &
                      pneigh%refine = 1

                    if (pneigh%level .eq. pmeta%level .and. pneigh%refine .eq. -1) &
                      pneigh%refine = 0
                  end if
                end do
              end do
            end do

          end if
        end if

        pmeta => pmeta%next
      end do

! check if all children are selected for derefinement; if this condition is not
! fullfiled deselect all children from derefinement
!
      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%level .eq. (l - 1) .and. .not. pmeta%leaf) then
          flag = .true.
          do p = 1, nchild
            pchild => pmeta%child(p)%ptr

            flag = flag .and. (pchild%refine .eq. -1)
          end do
          if (.not. flag) then
            do p = 1, nchild
              pchild => pmeta%child(p)%ptr
              if (pchild%leaf) &
                pchild%refine = max(0, pchild%refine)
            end do
          end if
        end if

        pmeta => pmeta%next
      end do

! deselect neighbors from derefinement if the current block is set for
! refinement or it is at higher level and not selected for derefinement
!
      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%level .eq. l .and. pmeta%leaf) then

          if (pmeta%refine .ne. -1) then

            do i = 1, ndims
              do j = 1, nsides
                do k = 1, nfaces
                  pneigh => pmeta%neigh(i,j,k)%ptr
                  if (associated(pneigh)) then
                    if (pneigh%refine .eq. -1) then
                      if (pneigh%level .lt. pmeta%level) &
                        pneigh%refine = 0
                    end if
                  end if
                end do
              end do
            end do

          end if

        end if

        pmeta => pmeta%next
      end do

    end do

#ifdef MPI
! find all sibling blocks which are spread over different processors
!
    pmeta => list_meta
    do while (associated(pmeta))
      if (.not. pmeta%leaf) then
        if (pmeta%child(1)%ptr%refine .eq. -1) then

! check if the parent blocks is on the same processor as the next block, if not
! move it to the same processor
!
          if (pmeta%cpu .ne. pmeta%next%cpu) &
            pmeta%cpu = pmeta%next%cpu

! find the case when child blocks are spread across at least 2 processors
!
          flag = .false.
          do p = 1, nchild
            flag = flag .or. (pmeta%child(p)%ptr%cpu .ne. pmeta%cpu)
          end do

          if (flag) then

! iterate over all children
!
            do p = 1, nchild

! generate the tag for communication
!
              itag = pmeta%child(p)%ptr%cpu * ncpus + pmeta%cpu + ncpus + p + 1

! if the current children is not on the same processor, then ...
!
              if (pmeta%child(p)%ptr%cpu .ne. pmeta%cpu) then

! allocate data blocks for children on the processor which will receive data
!
                if (pmeta%cpu .eq. ncpu) then
                  call append_datablock(pdata)
                  call associate_blocks(pmeta%child(p)%ptr, pdata)

! receive the data
!
                  call mrecvf(size(rbuf), pmeta%child(p)%ptr%cpu, itag, rbuf)

! coppy buffer to data
!
                  pmeta%child(p)%ptr%data%u(:,:,:,:) = rbuf(:,:,:,:)
                end if

! send data to the right processor and deallocate data block
!
                if (pmeta%child(p)%ptr%cpu .eq. ncpu) then

! copy data to buffer
!
                  rbuf(:,:,:,:) = pmeta%child(p)%ptr%data%u(:,:,:,:)

! send data
!
                  call msendf(size(rbuf), pmeta%cpu, itag, rbuf)

! deallocate data block
!
                  call deallocate_datablock(pmeta%child(p)%ptr%data)
                end if

! set the current processor of the block
!
                pmeta%child(p)%ptr%cpu = pmeta%cpu
              end if
            end do

          end if
        end if
      end if

      pmeta => pmeta%next
    end do
#endif /* MPI */

! perform the actual derefinement
!
    do l = maxlev, 2, -1

      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%leaf) then
          if (pmeta%level .eq. l) then
            if (pmeta%refine .eq. -1) then
              pparent => pmeta%parent

              if (associated(pparent)) then
#ifdef MPI
                if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */
                  if (.not. associated(pparent%data)) then
                    call append_datablock(pdata)
                    call associate_blocks(pparent, pdata)
                  end if
                  call restrict_block(pparent)
#ifdef MPI
                end if
#endif /* MPI */

                call derefine_block(pparent)
                pmeta => pparent
              end if
            end if
          end if
        end if

        pmeta => pmeta%next
      end do

    end do

! perform the actual refinement
!
    do l = 1, maxlev - 1

      pmeta => list_meta
      do while (associated(pmeta))
        if (pmeta%leaf) then
          if (pmeta%level .eq. l) then
            if (pmeta%refine .eq. 1) then
              pparent => pmeta
#ifdef MPI
              if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */
                call refine_block(pmeta, .true.)
                call prolong_block(pparent)
                call deallocate_datablock(pparent%data)
#ifdef MPI
              else
                call refine_block(pmeta, .false.)
              end if
#endif /* MPI */
            end if
          end if
        end if
        pmeta => pmeta%next
      end do

    end do

#ifdef MPI
!! AUTO BALANCING
!!
! calculate the new division
!
    l = nleafs / ncpus

! iterate over all metablocks and reassign the processor numbers
!
    n = 0
    p = 0

    pmeta => list_meta
    do while (associated(pmeta))

! assign the cpu to the current block
!
      if (pmeta%cpu .ne. n) then

        if (pmeta%leaf) then

! generate the tag for communication
!
          itag = pmeta%cpu * ncpus + n + ncpus + 1

          if (ncpu .eq. pmeta%cpu) then

! copy data to buffer
!
            rbuf(:,:,:,:) = pmeta%data%u(:,:,:,:)

! send data
!
            call msendf(size(rbuf), n, itag, rbuf)

! deallocate data block
!
             call deallocate_datablock(pmeta%data)

! send data block
!
          end if

          if (ncpu .eq. n) then

! allocate data block for the current block
!
            call append_datablock(pdata)
            call associate_blocks(pmeta, pdata)

! receive the data
!
            call mrecvf(size(rbuf), pmeta%cpu, itag, rbuf)

! coppy buffer to data
!
            pmeta%data%u(:,:,:,:) = rbuf(:,:,:,:)

! receive data block
!
          end if
        end if

! set new processor number
!
        pmeta%cpu = n

      end if

! increase the number of blocks on the current process; if it exceeds the
! allowed number reset the counter and increase the processor number
!
      if (pmeta%leaf) then
        p = p + 1
        if (p .ge. l) then
          n = min(ncpus - 1, n + 1)
          p = 0
        end if
      end if

! assign pointer to the next block
!
      pmeta => pmeta%next
    end do
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

    use blocks       , only : block_meta, nvars, nchild, ifl, iqt
#ifdef MHD
    use blocks       , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */
    use config       , only : ng, in, jn, kn, im, jm, km
    use interpolation, only : expand
#ifdef MHD
    use interpolation, only : magtocen
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: q, p
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, js, ks

! local arrays
!
    integer, dimension(3) :: dm, fm

! local allocatable arrays
!
    real, dimension(:,:,:,:), allocatable :: u

! local pointers
!
    type(block_meta), pointer :: pchild

!-------------------------------------------------------------------------------
!
! prepare dimensions
!
    dm(:) = (/ im, jm, km /)
    fm(:) = 2 * (dm(:) - ng)
#if NDIMS == 2
    fm(3) = 1
    ks    = 1
    kl    = 1
    ku    = 1
#endif /* NDIMS == 2 */

! allocate array to the product of expansion
!
    allocate(u(nvars, fm(1), fm(2), fm(3)))

! expand all variables and place them in the array u
!
    do q = 1, nvars
      call expand(dm, fm, ng, pblock%data%u(q,:,:,:), u(q,:,:,:), 't', 't', 't')
    end do
#ifdef MHD
! prolong face centered Bx
!
    call expand(dm, fm, ng, pblock%data%u(ibx,:,:,:), u(ibx,:,:,:), 'c', 'l', 'l')

! prolong face centered By
!
    call expand(dm, fm, ng, pblock%data%u(iby,:,:,:), u(iby,:,:,:), 'l', 'c', 'l')

#if NDIMS == 3
! prolong face centered Bz
!
    call expand(dm, fm, ng, pblock%data%u(ibz,:,:,:), u(ibz,:,:,:), 'l', 'l', 'c')
#endif /* NDIMS == 3 */

! calculate magnetic field at the cell centers
!
    call magtocen(fm(1), fm(2), fm(3), u(ibx:ibz,:,:,:), u(icx:icz,:,:,:))
#endif /* MHD */

! iterate over all children
!
    do p = 1, nchild

! assign pointer to the current child
!
      pchild => pblock%child(p)%ptr

! calculate the position of child in the parent block
!
      is = mod((p - 1)    ,2)
      js = mod((p - 1) / 2,2)
#if NDIMS == 3
      ks = mod((p - 1) / 4,2)
#endif /* NDIMS == 3 */

! calculate indices of the current child subdomain
!
      il = 1 + is * in
      jl = 1 + js * jn
#if NDIMS == 3
      kl = 1 + ks * kn
#endif /* NDIMS == 3 */

      iu = il + im - 1
      ju = jl + jm - 1
#if NDIMS == 3
      ku = kl + km - 1
#endif /* NDIMS == 3 */

! copy data to the current child
!
      pchild%data%u(1:nvars,1:im,1:jm,1:km) = u(1:nvars,il:iu,jl:ju,kl:ku)

    end do

! deallocate local arrays
!
    if (allocated(u)) deallocate(u)

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

    use blocks       , only : block_meta, nvars, nchild, ifl
#ifdef MHD
    use blocks       , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */
    use config       , only : ng, im, jm, km, ib, jb, kb
    use interpolation, only : shrink
#ifdef MHD
    use interpolation, only : magtocen
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k, q, p
    integer :: il, iu, jl, ju, kl, ku, i1, i2, j1, j2, k1, k2
    integer :: is, js, ks

! local arrays
!
    integer, dimension(3)     :: dm, fm, pm
    integer, dimension(3,0:1) :: cm

! local allocatable arrays
!
    real, dimension(:,:,:), allocatable :: u

! local pointers
!
    type(block_meta), pointer :: pchild

!-------------------------------------------------------------------------------
!
! prepare dimensions
!
    dm(:)   = (/ im, jm, km /)
    pm(:)   = dm(:) / 2
    fm(:)   = (dm(:) - ng) / 2
    cm(:,0) = ng
    cm(:,1) = dm(:) - ng
#if NDIMS == 2
    pm(3)   = 1
    fm(3)   = 1
    ks      = 1
    kl      = 1
    ku      = 1
    k1      = 1
    k2      = 1
    k       = 1
#endif /* NDIMS == 2 */

! allocate temporary array
!
    allocate(u(pm(1), pm(2), pm(3)))

! iterate over all children
!
    do p = 1, nchild

! assign pointer to the current child
!
      pchild => pblock%child(p)%ptr

! calculate the position of child in the parent block
!
      is = mod((p - 1)    ,2)
      js = mod((p - 1) / 2,2)
#if NDIMS == 3
      ks = mod((p - 1) / 4,2)
#endif /* NDIMS == 3 */

! calculate the bounds of the input array indices
!
      i1 =  1 + ng / 2 * is
      j1 =  1 + ng / 2 * js
#if NDIMS == 3
      k1 =  1 + ng / 2 * ks
#endif /* NDIMS == 3 */

      i2 = i1 + fm(1) - 1
      j2 = j1 + fm(2) - 1
      k2 = k1 + fm(3) - 1

! calculate the bounds of the destination array indices
!
      il = ib - ng / 2 + is * fm(1)
      jl = jb - ng / 2 + js * fm(2)
#if NDIMS == 3
      kl = kb - ng / 2 + ks * fm(3)
#endif /* NDIMS == 3 */

      iu = il + fm(1) - 1
      ju = jl + fm(2) - 1
      ku = kl + fm(3) - 1

! iterate over all quantities
!
      do q = 1, ifl

! shrink the current child
!
        call shrink(dm, pm, 0, pchild%data%u(q,:,:,:), u(:,:,:), 'm', 'm', 'm')

! fill the parent block
!
        pblock%data%u(q,il:iu,jl:ju,kl:ku) = u(i1:i2,j1:j2,k1:k2)

      end do
#ifdef MHD
! shrink face centered Bx
!
      call shrink(dm, pm, 0, pchild%data%u(ibx,:,:,:), u(:,:,:), 'c', 'l', 'l')
      pblock%data%u(ibx,il:iu,jl:ju,kl:ku) = u(i1:i2,j1:j2,k1:k2)

! shrink face centered By
!
      call shrink(dm, pm, 0, pchild%data%u(iby,:,:,:), u(:,:,:), 'l', 'c', 'l')
      pblock%data%u(iby,il:iu,jl:ju,kl:ku) = u(i1:i2,j1:j2,k1:k2)

#if NDIMS == 3
! shrink face centered Bz
!
      call shrink(dm, pm, 0, pchild%data%u(ibz,:,:,:), u(:,:,:), 'l', 'l', 'c')
      pblock%data%u(ibz,il:iu,jl:ju,kl:ku) = u(i1:i2,j1:j2,k1:k2)
#endif /* NDIMS == 3 */
#endif /* MHD */

    end do

#ifdef MHD
! calculate magnetic field at the cell centers
!
    call magtocen(im, jm, km, pblock%data%u(ibx:ibz,:,:,:)                     &
                            , pblock%data%u(icx:icz,:,:,:))
#endif /* MHD */

! deallocate temporary array
!
    deallocate(u)

! ! prepare dimensions
! !
!     dm(:)   = (/ im, jm, km /)
!     fm(:)   = (dm(:) - ng) / 2
!     cm(:,0) = ng
!     cm(:,1) = dm(:) - ng
! #if NDIMS == 2
!     fm(3)   = 1
!     ks      = 1
!     kl      = 1
!     k1      = 1
!     k2      = 1
!     k       = 1
! #endif /* NDIMS == 2 */
!
! ! iterate over all children
! !
!     do p = 1, nchild
!
! ! assign pointer to the current child
! !
!       pchild => pblock%child(p)%ptr
!
! ! calculate the position of child in the parent block
! !
!       is = mod((p - 1)    ,2)
!       js = mod((p - 1) / 2,2)
! #if NDIMS == 3
!       ks = mod((p - 1) / 4,2)
! #endif /* NDIMS == 3 */
!
! ! calculate the bounds of the destination array indices
! !
!       il = ib - ng / 2 + is * fm(1)
!       jl = jb - ng / 2 + js * fm(2)
! #if NDIMS == 3
!       kl = kb - ng / 2 + ks * fm(3)
! #endif /* NDIMS == 3 */
!
!       iu = il + fm(1) - 1
!       ju = jl + fm(2) - 1
!       ku = kl + fm(3) - 1
!
! ! perform the current block restriction
! !
! #if NDIMS == 3
!       do k = kl, ku
!         k2 = 2 * k - cm(3,ks)
!         k1 = k2 - 1
! #endif /* NDIMS == 3 */
!         do j = jl, ju
!           j2 = 2 * j - cm(2,js)
!           j1 = j2 - 1
!           do i = il, iu
!             i2 = 2 * i - cm(1,is)
!             i1 = i2 - 1
!
!             do q = 1, nvars
!               pblock%data%u(q,i,j,k) = sum(pchild%data%u(q,i1:i2,j1:j2,k1:k2)) / nchild
!             end do
!           end do
!         end do
! #if NDIMS == 3
!       end do
! #endif /* NDIMS == 3 */
!
!     end do
!
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
    if (allocated(ax))   deallocate(ax)
    if (allocated(ay))   deallocate(ay)
    if (allocated(az))   deallocate(az)
    if (allocated(adx))  deallocate(adx)
    if (allocated(ady))  deallocate(ady)
    if (allocated(adz))  deallocate(adz)
    if (allocated(adxi)) deallocate(adxi)
    if (allocated(adyi)) deallocate(adyi)
    if (allocated(adzi)) deallocate(adzi)

!-------------------------------------------------------------------------------
!
  end subroutine clear_mesh

!===============================================================================
!
end module
