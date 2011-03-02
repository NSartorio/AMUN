!!******************************************************************************
!!
!! module: mesh - handling adaptive mesh structure
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
  real, dimension(:  ), allocatable, save :: advol

  contains
!
!===============================================================================
!
! init_mesh: subroutine initializes mesh by creating blocks according
!            to the geometry, initial problem and refinement criterium
!
!===============================================================================
!
  subroutine init_mesh()

    use config  , only : im, jm, km, xmin, xmax, ymin, ymax, zmin, zmax        &
                       , ncells, maxlev, rdims, ng
    use blocks  , only : block_meta, block_data, list_meta, list_data          &
                       , init_blocks, refine_block               &
                       , deallocate_datablock, mblocks, nleafs, dblocks        &
                       , nchild, ndims, nsides, nfaces, res
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
    integer(kind=4)      :: i, j, k, l, n
    character(len=64)    :: fmt
    character(len=32)    :: bstr, tstr
#ifdef MPI
    integer(kind=4), dimension(0:ncpus-1) :: lb
#endif /* MPI */

!-------------------------------------------------------------------------------
!
! initialize blocks
!
    call init_blocks()

! allocate the effective resolution array
!
    allocate(res(maxlev))

! calculate the effective resolution at each level
!
    do l = 1, maxlev
      res(l) = ncells * 2**(maxlev - l)
    end do

! allocate the initial structure of blocks according to the problem
!
    call init_domain

! print general information about resolutions
!
    if (is_master()) then
      write(*,"(1x,a)"         ) "Generating the initial mesh:"
      write(*,"(4x,a,3(1x,i6))") "base configuration     =", rdims(1:ndims)
      write(*,"(4x,a,3(1x,i6))") "base resolution        =", rdims(1:ndims) * ncells
      write(*,"(4x,a,3(1x,i6))") "effective resolution   =", rdims(1:ndims) * res(1)
      write(*,"(4x,a,  1x,i6)" ) "refinement to level    =", maxlev
    end if

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out the coarse blocks with the initial condition
!
    if (is_master()) &
      write(*,"(4x,a,$)") "generating level       =    "

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
            if (mblocks .eq. 1 .and. l .eq. 1) &
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
    l       = mod(nleafs, ncpus) - 1
    lb( : ) = nleafs / ncpus
    lb(0:l) = lb(0:l) + 1

! reset the processor and block numbers
    n = 0
    l = 0

    pmeta_block => list_meta
    do while (associated(pmeta_block))

! assign the cpu to the current block
!
      pmeta_block%cpu = n

! increase the number of blocks on the current process; if it exceeds the
! allowed number reset the counter and increase the processor number
!
      if (pmeta_block%leaf) then
        l = l + 1
        if (l .ge. lb(n)) then
          n = min(ncpus - 1, n + 1)
          l = 0
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
      n = 0
      do l = 0, maxlev - 1
        k = 2**(ndims * l)
        n = n + k
      end do
      n = n * rdims(1) * rdims(2) * rdims(3)
      k = k * rdims(1) * rdims(2) * rdims(3)

      i = nint(alog10(1.0*mblocks + 1)) + 1
      j = nint(alog10(1.0*n + 1)) + 1

      write(fmt, "(a,i1,a,i1,a)") "(4x,a,1x,i", i, ",' / ',i", j, ",' = ',f8.4,' %')"

      write(*,*)
      write(*,fmt) "leafs    /cover blocks =", nleafs , k, (100.0 * nleafs ) / k
      write(*,fmt) "allocated/total blocks =", mblocks, n, (100.0 * mblocks) / n
    end if

! allocating space for coordinate variables
!
    allocate(ax   (maxlev, im))
    allocate(ay   (maxlev, jm))
    allocate(az   (maxlev, km))
    allocate(adx  (maxlev))
    allocate(ady  (maxlev))
    allocate(adz  (maxlev))
    allocate(adxi (maxlev))
    allocate(adyi (maxlev))
    allocate(adzi (maxlev))
    allocate(advol(maxlev))

! reset the coordinate variables
!
    ax(:,:)  = 0.0d0
    ay(:,:)  = 0.0d0
    az(:,:)  = 0.0d0
    adx(:)   = 1.0d0
    ady(:)   = 1.0d0
    adz(:)   = 1.0d0
    adxi(:)  = 1.0d0
    adyi(:)  = 1.0d0
    adzi(:)  = 1.0d0
    advol(:) = 1.0d0

! generating coordinates for all levels
!
    do l = 1, maxlev
      n = ncells * 2**(l - 1)

      adx (l) = (xmax - xmin) / (rdims(1) * n)
      ady (l) = (ymax - ymin) / (rdims(2) * n)
#if NDIMS == 3
      adz (l) = (zmax - zmin) / (rdims(3) * n)
#endif /* NDIMS == 3 */

      ax(l,:) = ((/(i, i = 1, im)/) - ng - 0.5d0) * adx(l)
      ay(l,:) = ((/(j, j = 1, jm)/) - ng - 0.5d0) * ady(l)
#if NDIMS == 3
      az(l,:) = ((/(k, k = 1, km)/) - ng - 0.5d0) * adz(l)
#endif /* NDIMS == 3 */

      adxi(l) = 1.0d0 / adx(l)
      adyi(l) = 1.0d0 / ady(l)
#if NDIMS == 3
      adzi(l) = 1.0d0 / adz(l)
#endif /* NDIMS == 3 */

      advol(l) = adx(l) * ady(l) * adz(l)
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
  subroutine update_mesh()

    use config   , only : maxlev, im, jm, km
    use blocks   , only : block_meta, block_data, list_meta, list_data         &
                        , nleafs, dblocks, nchild, ndims, nsides, nfaces       &
                        , refine_block, derefine_block, append_datablock       &
                        , associate_blocks, deallocate_datablock
    use error    , only : print_info
#ifdef MPI
    use mpitools , only : ncpus, ncpu, is_master, mallreducesuml, msendf, mrecvf
#endif /* MPI */
    use problem  , only : check_ref
    use variables, only : nqt

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

! array for number of data block for autobalancing
!
    integer(kind=4), dimension(0:ncpus-1) :: lb

! local buffer for data block exchange
!
    real(kind=8)   , dimension(nqt,im,jm,km) :: rbuf
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

        if (pmeta%level .eq. 1) &
          pmeta%refine = max(0, pmeta%refine)
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
    l       = mod(nleafs, ncpus) - 1
    lb( : ) = nleafs / ncpus
    lb(0:l) = lb(0:l) + 1

! iterate over all metablocks and reassign the processor numbers
!
    n = 0
    l = 0

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
        l = l + 1
        if (l .ge. lb(n)) then
          n = min(ncpus - 1, n + 1)
          l = 0
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

    use blocks       , only : block_meta, block_data, nchild
    use config       , only : ng, nh, in, jn, kn, im, jm, km
    use interpolation, only : expand
    use variables    , only : nfl, nqt
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
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
    type(block_data), pointer :: pdata

!-------------------------------------------------------------------------------
!
! assign the pdata pointer
!
    pdata => pblock%data

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
    allocate(u(nqt, fm(1), fm(2), fm(3)))

! expand all variables and place them in the array u
!
    do q = 1, nfl
      call expand(dm(:), fm(:), nh, pdata%u(q,:,:,:), u(q,:,:,:))
    end do

#ifdef MHD
! expand the cell centered magnetic field components
!
    do q = ibx, ibz
      call expand(dm(:), fm(:), nh, pdata%u(q,:,:,:), u(q,:,:,:))
    end do
#ifdef GLM

! expand the scalar potential Psi
!
    call expand(dm(:), fm(:), nh, pdata%u(iph,:,:,:), u(iph,:,:,:))
#endif /* GLM */
#endif /* MHD */
! iterate over all children
!
    do p = 1, nchild

! assign pointer to the current child
!
      pchild => pblock%child(p)%ptr

! obtain the position of child in the parent block
!
      is = pchild%pos(1)
      js = pchild%pos(2)
#if NDIMS == 3
      ks = pchild%pos(3)
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
      pchild%data%u(1:nqt,1:im,1:jm,1:km) = u(1:nqt,il:iu,jl:ju,kl:ku)

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

    use blocks       , only : block_meta, block_data, nchild
    use config       , only : ng, in, ih, im, ib, ie, nh, jn, jh, jm, jb, je   &
                                , kn, kh, km, kb, ke
    use variables    , only : nfl
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: p
    integer :: if, jf, kf
    integer :: il, jl, kl, iu, ju, ku
    integer :: ip, jp, kp
    integer :: is, js, ks, it, jt, kt

! local pointers
!
    type(block_data), pointer :: pparent, pchild

!-------------------------------------------------------------------------------
!
! assign pointers
!
    pparent => pblock%data

! iterate over all children
!
    do p = 1, nchild

! assign pointer to the current child
!
      pchild  => pblock%child(p)%ptr%data

! obtain the position of the current child in the parent block
!
      if = pchild%meta%pos(1)
      jf = pchild%meta%pos(2)
#if NDIMS == 3
      kf = pchild%meta%pos(3)
#endif /* NDIMS == 3 */

! calculate the bound indices of the source nad destination arrays
!
      if (if .eq. 0) then
        il = 1
        iu = ie
        is = ib - nh
        it = ih
      else
        il = ib
        iu = im
        is = ih + 1
        it = ie + nh
      end if
      ip = il + 1
      if (jf .eq. 0) then
        jl = 1
        ju = je
        js = jb - nh
        jt = jh
      else
        jl = jb
        ju = jm
        js = jh + 1
        jt = je + nh
      end if
      jp = jl + 1
#if NDIMS == 3
      if (kf .eq. 0) then
        kl = 1
        ku = ke
        ks = kb - nh
        kt = kh
      else
        kl = kb
        ku = km
        ks = kh + 1
        kt = ke + nh
      end if
      kp = kl + 1
#endif /* NDIMS == 3 */

! copy the variables from the current child to the proper location of
! the parent block
!
#if NDIMS == 2
      pparent%u(1:nfl,is:it,js:jt,1) =                                         &
                                   0.25 * (pchild%u(1:nfl,il:iu:2,jl:ju:2,1)   &
                                         + pchild%u(1:nfl,ip:iu:2,jl:ju:2,1)   &
                                         + pchild%u(1:nfl,il:iu:2,jp:ju:2,1)   &
                                         + pchild%u(1:nfl,ip:iu:2,jp:ju:2,1))

#ifdef MHD
      pparent%u(ibx:ibz,is:it,js:jt,1) =                                       &
                                 0.25 * (pchild%u(ibx:ibz,il:iu:2,jl:ju:2,1)   &
                                       + pchild%u(ibx:ibz,ip:iu:2,jl:ju:2,1)   &
                                       + pchild%u(ibx:ibz,il:iu:2,jp:ju:2,1)   &
                                       + pchild%u(ibx:ibz,ip:iu:2,jp:ju:2,1))
#ifdef GLM
      pparent%u(iph    ,is:it,js:jt,1) =                                       &
                                 0.25 * (pchild%u(iph    ,il:iu:2,jl:ju:2,1)   &
                                       + pchild%u(iph    ,ip:iu:2,jl:ju:2,1)   &
                                       + pchild%u(iph    ,il:iu:2,jp:ju:2,1)   &
                                       + pchild%u(iph    ,ip:iu:2,jp:ju:2,1))
#endif /* GLM */
#endif /* MHD */
#endif /* NDIMS == 2 */
#if NDIMS == 3
      pparent%u(1:nfl,is:it,js:jt,ks:kt) =                                     &
                             0.125 * (pchild%u(1:nfl,il:iu:2,jl:ju:2,kl:ku:2)  &
                                    + pchild%u(1:nfl,ip:iu:2,jl:ju:2,kl:ku:2)  &
                                    + pchild%u(1:nfl,il:iu:2,jp:ju:2,kl:ku:2)  &
                                    + pchild%u(1:nfl,ip:iu:2,jp:ju:2,kl:ku:2)  &
                                    + pchild%u(1:nfl,il:iu:2,jl:ju:2,kp:ku:2)  &
                                    + pchild%u(1:nfl,ip:iu:2,jl:ju:2,kp:ku:2)  &
                                    + pchild%u(1:nfl,il:iu:2,jp:ju:2,kp:ku:2)  &
                                    + pchild%u(1:nfl,ip:iu:2,jp:ju:2,kp:ku:2))
#ifdef MHD
      pparent%u(ibx:ibz,is:it,js:jt,ks:kt) =                                   &
                           0.125 * (pchild%u(ibx:ibz,il:iu:2,jl:ju:2,kl:ku:2)  &
                                  + pchild%u(ibx:ibz,ip:iu:2,jl:ju:2,kl:ku:2)  &
                                  + pchild%u(ibx:ibz,il:iu:2,jp:ju:2,kl:ku:2)  &
                                  + pchild%u(ibx:ibz,ip:iu:2,jp:ju:2,kl:ku:2)  &
                                  + pchild%u(ibx:ibz,il:iu:2,jl:ju:2,kp:ku:2)  &
                                  + pchild%u(ibx:ibz,ip:iu:2,jl:ju:2,kp:ku:2)  &
                                  + pchild%u(ibx:ibz,il:iu:2,jp:ju:2,kp:ku:2)  &
                                  + pchild%u(ibx:ibz,ip:iu:2,jp:ju:2,kp:ku:2))
#ifdef GLM
      pparent%u(iph    ,is:it,js:jt,ks:kt) =                                   &
                           0.125 * (pchild%u(iph    ,il:iu:2,jl:ju:2,kl:ku:2)  &
                                  + pchild%u(iph    ,ip:iu:2,jl:ju:2,kl:ku:2)  &
                                  + pchild%u(iph    ,il:iu:2,jp:ju:2,kl:ku:2)  &
                                  + pchild%u(iph    ,ip:iu:2,jp:ju:2,kl:ku:2)  &
                                  + pchild%u(iph    ,il:iu:2,jl:ju:2,kp:ku:2)  &
                                  + pchild%u(iph    ,ip:iu:2,jl:ju:2,kp:ku:2)  &
                                  + pchild%u(iph    ,il:iu:2,jp:ju:2,kp:ku:2)  &
                                  + pchild%u(iph    ,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* GLM */
#endif /* MHD */
#endif /* NDIMS == 3 */
    end do
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
  subroutine clear_mesh()

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
    if (allocated(ax)   ) deallocate(ax)
    if (allocated(ay)   ) deallocate(ay)
    if (allocated(az)   ) deallocate(az)
    if (allocated(adx)  ) deallocate(adx)
    if (allocated(ady)  ) deallocate(ady)
    if (allocated(adz)  ) deallocate(adz)
    if (allocated(adxi) ) deallocate(adxi)
    if (allocated(adyi) ) deallocate(adyi)
    if (allocated(adzi) ) deallocate(adzi)
    if (allocated(advol)) deallocate(advol)

!-------------------------------------------------------------------------------
!
  end subroutine clear_mesh

!===============================================================================
!
end module
