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
!! module: MESH - handling adaptive mesh structure
!!
!!******************************************************************************
!
module mesh

  implicit none

! maximum number of covering blocks
!
  integer                 , save :: tblocks = 1

! log file for the mesh statistics
!
  character(len=32), save :: fname = "mesh.log"
  integer          , save :: funit = 11

  contains
!
!===============================================================================
!
! subroutine INITIALIZE_MESH:
! --------------------------
!
!   Subroutine initializes mesh module and its variables.
!
!===============================================================================
!
  subroutine initialize_mesh(flag)

    use blocks     , only : datablock_set_dims
    use coordinates, only : xmin, xmax, ymin, ymax, zmin, zmax
    use coordinates, only : toplev, im, jm, km
    use mpitools   , only : master, nprocs
    use equations  , only : nv

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    logical, intent(in) :: flag

! local variables
!
    integer(kind=4) :: i, j, k, l, n
    logical         :: info

! local arrays
!
    integer(kind=4), dimension(3) :: bm, rm, dm
!
!-------------------------------------------------------------------------------
!
! set data block dimensions
!
    call datablock_set_dims(nv, nv, im, jm, km)

! print general information about resolutions
!
    if (master) then

! prepare the file for logging mesh statistics; only the master process handles
! this part
!
! check if the mesh statistics file exists
!
      inquire(file = fname, exist = info)

! if flag is set to .true. or the original mesh statistics file does not exist,
! create a new one and write the header in it, otherwise open the original file
! and move the writing position to the end to allow for appending
!
      if (flag .or. .not. info) then

! create a new mesh statistics file
!
        open(unit=funit, file=fname, form='formatted', status='replace')

! write the mesh statistics file header
!
        write(funit,"('#')")
        write(funit,"('#',4x,'step',5x,'time',11x,'leaf',4x,'meta'," //        &
                           "6x,'coverage',8x,'AMR',11x,'block distribution')")
        write(funit,"('#',27x,'blocks',2x,'blocks',4x,'efficiency'," //        &
                                              "6x,'efficiency')",advance='no')
        write(funit,"(1x,a12)",advance="no") 'level = 1'
        do l = 2, toplev
          write(funit,"(2x,i6)",advance="no") l
        end do
#ifdef MPI
        write(funit,"(1x,a10)",advance="no") 'cpu = 1'
        do n = 2, nprocs
          write(funit,"(2x,i6)",advance="no") n
        end do
#endif /* MPI */
        write(funit,"('' )")
        write(funit,"('#')")

      else

! open the mesh statistics file and set the position at the end of file
!
        open(unit=funit, file=fname, form='formatted', position='append')

! write a marker that the job has been restarted from here
!
        write(funit,"('#',1x,a)") "job restarted from this point"

      end if

    end if ! master

!-------------------------------------------------------------------------------
!
  end subroutine initialize_mesh
!
!===============================================================================
!
! generate_mesh: subroutine generate the initial block structure by creating
!                blocks according to the geometry, initial problem and
!                refinement criterium
!
!===============================================================================
!
  subroutine generate_mesh()

    use blocks  , only : block_meta, block_data, list_meta, list_data
    use blocks  , only : refine_block, remove_datablock
    use blocks  , only : nchild, nsides, nfaces
    use blocks  , only : get_mblocks, get_nleafs
    use coordinates, only : minlev, maxlev, res
    use domains    , only : setup_domain
    use error   , only : print_info, print_error
    use mpitools, only : master, nproc, nprocs
    use problems   , only : setup_problem
    use refinement , only : check_refinement_criterion

    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh, pnext
    type(block_data), pointer :: pdata

! local variables
!
    integer(kind=4)                       :: i, j, k, l, n
#ifdef MPI
    integer(kind=4), dimension(0:nprocs-1) :: lb
#endif /* MPI */

!-------------------------------------------------------------------------------
!
! allocate the initial structure of blocks according to the problem
!
    call setup_domain()

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out the coarse blocks with the initial condition
!
    if (master) then
      write(*,*)
      write(*,"(1x,a)"             ) "Generating the initial mesh:"
      write(*,"(4x,a)",advance="no") "generating level       =    "
    end if

    l = 1
    do while (l .le. maxlev)

! print the level currently processed
!
      if (master) &
        write(*,"(1x,i2)",advance="no") l

! iterate over all data blocks at the current level and initialize the problem
!
      pdata => list_data
      do while (associated(pdata))

! set the initial conditions at the current block
!
        if (pdata%meta%level .le. l) &
          call setup_problem(pdata)

! assign pointer to the next block
!
        pdata => pdata%next
      end do

! at the maximum level we only initialize the problem (without checking the
! refinement criterion)
!
      if (l .lt. maxlev) then

! iterate over all data blocks at the current level and check the refinement
! criterion; do not allow for derefinement
!
        pdata => list_data
        do while (associated(pdata))

          if (pdata%meta%level .eq. l) then
            pdata%meta%refine = max(0, check_refinement_criterion(pdata))

! if there is only one block, and it is set not to be refined, refine it anyway
! because the resolution for the problem initiation may be too small
!
            if (get_mblocks() .eq. 1 .and. l .eq. 1) &
              pdata%meta%refine = 1

! if the level is lower then the minimum level set the block to be refined
! anyway
!
            if (l .lt. minlev) pdata%meta%refine = 1

          end if

! assign pointer to the next block
!
          pdata => pdata%next
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
          pmeta => list_meta
          do while (associated(pmeta))

! check if the current block is at the level n, is a leaf, and selected for
! refinement
!
            if (pmeta%level .eq. n) then
              if (pmeta%leaf) then
                if (pmeta%refine .eq. 1) then

! iterate over all neighbors
!
                  do i = 1, NDIMS
                    do j = 1, nsides
                      do k = 1, nfaces

! assign pointer to the neighbor
!
                        pneigh => pmeta%neigh(i,j,k)%ptr

! check if the neighbor is associated
!
                        if (associated(pneigh)) then

! check if the neighbor is a leaf, if not something wrong is going on
!
                          if (pneigh%leaf) then

! if the neighbor has lower level, select it to be refined too
!
                            if (pneigh%level .lt. pmeta%level) &
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
            pmeta => pmeta%next

          end do
        end do

!! refine all selected blocks starting from the lowest level
!!
! walk through the levels starting from the lowest to the current level
!
        do n = 1, l

! iterate over all meta blocks
!
          pmeta => list_meta
          do while (associated(pmeta))

! check if the current block is at the level n and, if it is selected for
! refinement and if so, perform the refinement on this block
!
            if (pmeta%level .eq. n .and. pmeta%refine .eq. 1) then

! perform the refinement
!
              call refine_block(pmeta, res(pmeta%level + 1,:), .true.)

            end if

! assign pointer to the next block
!
            pmeta => pmeta%next
          end do
        end do
      end if

      l = l + 1
    end do

! deallocate data blocks of non leafs
!
    pmeta => list_meta
    do while (associated(pmeta))

      if (.not. pmeta%leaf) &
        call remove_datablock(pmeta%data)

! assign pointer to the next block
!
      pmeta => pmeta%next
    end do

#ifdef MPI
! divide blocks between all processes, use the number of data blocks to do this,
! but keep blocks from the top level which have the same parent packed together
!
    l       = mod(get_nleafs(), nprocs) - 1
    lb( : ) = get_nleafs() / nprocs
    lb(0:l) = lb(0:l) + 1

! reset the processor and block numbers
    n = 0
    l = 0

    pmeta => list_meta
    do while (associated(pmeta))

! assign the cpu to the current block
!
      pmeta%cpu = n

! increase the number of blocks on the current process; if it exceeds the
! allowed number reset the counter and increase the processor number
!
      if (pmeta%leaf) then
        l = l + 1
        if (l .ge. lb(n)) then
          n = min(nprocs - 1, n + 1)
          l = 0
        end if
      end if

! assign pointer to the next block
!
      pmeta => pmeta%next
    end do

! remove all data blocks which do not belong to the current process
!
    pmeta => list_meta
    do while (associated(pmeta))
      pnext => pmeta%next

! if the current block belongs to another process and its data field is
! associated, deallocate its data field
!
      if (pmeta%cpu .ne. nproc .and. associated(pmeta%data)) &
        call remove_datablock(pmeta%data)

! assign pointer to the next block
!
      pmeta => pnext
    end do
#endif /* MPI */

! go to a new line after generating levels
!
    if (master) then
      write(*,*)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine generate_mesh
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

    use blocks   , only : block_meta, block_data, list_meta, list_data         &
                        , nchild, ndims, nsides, nfaces                        &
                        , refine_block, derefine_block, append_datablock       &
                        , associate_blocks, remove_datablock
    use blocks   , only : get_nleafs
    use coordinates, only : minlev, maxlev, toplev, im, jm, km, res
    use error    , only : print_info, print_error
#ifdef MPI
    use mpitools , only : reduce_sum_integer_array
    use mpitools , only : send_real_array, receive_real_array
    use mpitools , only : master, nprocs, nproc
#endif /* MPI */
    use refinement , only : check_refinement_criterion
    use equations  , only : nv

    implicit none

! local variables
!
    logical         :: flag
    integer(kind=4) :: i, j, k, l, n, p
    integer         :: iret

#ifdef MPI
! tag for the MPI data exchange
!
    integer(kind=4)                            :: itag

! array for storing the refinement flags
!
    integer(kind=4), dimension(:), allocatable :: ibuf

! array for number of data block for autobalancing
!
    integer(kind=4), dimension(0:nprocs-1)      :: lb

! local buffer for data block exchange
!
    real(kind=8)   , dimension(nv,im,jm,km)     :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh, pchild, pparent
    type(block_data), pointer :: pdata

!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! check mesh
!
    call check_mesh('before update_mesh')

#endif /* DEBUG */
! iterate over elements of the data block list
!
    pdata => list_data
    do while (associated(pdata))

! assign a pointer to the meta block associated with the current data block
!
      pmeta => pdata%meta

! if the current data block has a meta block associated
!
      if (associated(pmeta)) then

! if the associated meta block is a leaf
!
        if (pmeta%leaf) then

! check the refinement criterion for the current data block
!
          pmeta%refine = check_refinement_criterion(pdata)

! correct the refinement of the block for the base and top levels
!
          if (pmeta%level .lt. minlev) pmeta%refine =  1
          if (pmeta%level .eq. minlev) pmeta%refine = max(0, pmeta%refine)
          if (pmeta%level .eq. maxlev) pmeta%refine = min(0, pmeta%refine)
          if (pmeta%level .gt. maxlev) pmeta%refine = -1

        end if ! pmeta is a leaf
      end if ! pmeta associated

! assign a pointer to the next data block
!
      pdata => pdata%next

   end do

#ifdef MPI
! allocate buffer for the refinement field values
!
    allocate(ibuf(get_nleafs()))

! reset the buffer
!
    ibuf(:) = 0

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

! update refinement flags across all processors
!
    call reduce_sum_integer_array(get_nleafs(), ibuf(1:get_nleafs()), iret)

! update non-local block refinement flags
!
    l = 1
    pmeta => list_meta
    do while (associated(pmeta))
      if (pmeta%leaf) then
        pmeta%refine = ibuf(l)

        l = l + 1
      end if
      pmeta => pmeta%next
    end do

! deallocate the buffer
!
    if (allocated(ibuf)) deallocate(ibuf)

#endif /* MPI */
! iterate over all levels starting from top and correct the refinement of blocks
!
    do l = toplev, 1, -1

! iterate over all meta blocks
!
      pmeta => list_meta
      do while (associated(pmeta))

! check only leafs at the current level
!
        if (pmeta%leaf .and. pmeta%level .eq. l) then

! iterte over all neighbors of the current leaf
!
          do i = 1, ndims
            do j = 1, nsides
              do k = 1, nfaces

! assign a pointer to the current neighbor
!
                pneigh => pmeta%neigh(i,j,k)%ptr

! check if the pointer is associated with any block
!
                if (associated(pneigh)) then

!= conditions for blocks which are selected to be refined
!
                  if (pmeta%refine .eq. 1) then

! if the neighbor is set to be derefined, reset its flags (this applies to
! blocks at the current and lower levels)
!
                    pneigh%refine = max(0, pneigh%refine)

! if the neighbor is at lower level, always set it to be refined
!
                    if (pneigh%level .lt. pmeta%level) pneigh%refine = 1

                  end if ! refine = 1

!= conditions for blocks which stay at the same level
!
                  if (pmeta%refine .eq. 0) then

! if the neighbor lays at lower level and is set to be derefined, cancel its
! derefinement
!
                    if (pneigh%level .lt. pmeta%level)                         &
                                         pneigh%refine = max(0, pneigh%refine)

                  end if ! refine = 0

!= conditions for blocks which are selected to be derefined
!
                  if (pmeta%refine .eq. -1) then

! if the neighbor is at lower level and is set to be derefined, cancel its
! derefinement
!
                    if (pneigh%level .lt. pmeta%level)                         &
                                         pneigh%refine = max(0, pneigh%refine)

! if the neighbor is set to be refined, cancel derefinement of the current block
!
                    if (pneigh%refine .eq. 1) pmeta%refine = 0

                  end if ! refine = -1

                end if ! associated(pneigh)

               end do
             end do
          end do

        end if ! leafs at level l

! assign a pointer to the next block
!
        pmeta => pmeta%next

      end do ! meta blocks

!= now check all derefined block if their siblings are set for derefinement too
!  and are at the same level; check ony levels >= 2
!
      if (l .ge. 2) then

! iterate over all blocks
!
        pmeta => list_meta
        do while (associated(pmeta))

! check only leafs at the current level
!
          if (pmeta%leaf .and. pmeta%level .eq. l) then

! check blocks which are selected to be derefined
!
            if (pmeta%refine .eq. -1) then

! assign a pointer to the parent of the current block
!
              pparent => pmeta%parent

! check if parent is associated with any block
!
              if (associated(pparent)) then

! reset derefinement flag
!
                flag = .true.

! iterate over all children
!
                do p = 1, nchild

! assign a pointer to the current child
!
                  pchild => pparent%child(p)%ptr

! check if the current child is a leaf
!
                  flag = flag .and. (pchild%leaf)

! check if the current child is set to be derefined
!
                  flag = flag .and. (pchild%refine .eq. -1)

                end do ! over all children

! if not all children are proper for derefinement, cancel derefinement of all
! children
!
                if (.not. flag) then

! iterate over all children
!
                  do p = 1, nchild

! assign a pointer to the current child
!
                    pchild => pparent%child(p)%ptr

! reset derefinement of the current child
!
                    pchild%refine = max(0, pchild%refine)

                  end do ! children

                end if ! ~flag

              end if ! pparent is associated

            end if ! refine = -1

          end if ! leafs at level l

! assign a pointer to the next block
!
          pmeta => pmeta%next

        end do ! meta blocks

      end if ! l >= 2

    end do ! levels

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
              itag = pmeta%child(p)%ptr%cpu * nprocs + pmeta%cpu + nprocs + p + 1

! if the current children is not on the same processor, then ...
!
              if (pmeta%child(p)%ptr%cpu .ne. pmeta%cpu) then

! allocate data blocks for children on the processor which will receive data
!
                if (pmeta%cpu .eq. nproc) then
                  call append_datablock(pdata)
                  call associate_blocks(pmeta%child(p)%ptr, pdata)

! receive the data
!
                  call receive_real_array(size(rbuf), pmeta%child(p)%ptr%cpu, itag, rbuf, iret)

! coppy buffer to data
!
                  pmeta%child(p)%ptr%data%u(:,:,:,:) = rbuf(:,:,:,:)
                end if

! send data to the right processor and deallocate data block
!
                if (pmeta%child(p)%ptr%cpu .eq. nproc) then

! copy data to buffer
!
                  rbuf(:,:,:,:) = pmeta%child(p)%ptr%data%u(:,:,:,:)

! send data
!
                  call send_real_array(size(rbuf), pmeta%cpu, itag, rbuf, iret)

! deallocate data block
!
                  call remove_datablock(pmeta%child(p)%ptr%data)
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
    do l = toplev, 2, -1

      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%leaf) then
          if (pmeta%level .eq. l) then
            if (pmeta%refine .eq. -1) then
              pparent => pmeta%parent

              if (associated(pparent)) then
#ifdef MPI
                if (pmeta%cpu .eq. nproc) then
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
              else
                call print_error("mesh::update_mesh"                           &
                               , "Parent of derefined block is not associated!")
              end if
            end if
          end if
        end if

        pmeta => pmeta%next
      end do

    end do

! perform the actual refinement
!
    do l = 1, toplev - 1

      pmeta => list_meta
      do while (associated(pmeta))
        if (pmeta%leaf) then
          if (pmeta%level .eq. l) then
            if (pmeta%refine .eq. 1) then
              pparent => pmeta
#ifdef MPI
              if (pmeta%cpu .eq. nproc) then
#endif /* MPI */
                call refine_block(pmeta, res(pmeta%level + 1,:), .true.)
                call prolong_block(pparent)
                call remove_datablock(pparent%data)
#ifdef MPI
              else
                call refine_block(pmeta, res(pmeta%level + 1,:), .false.)
              end if
#endif /* MPI */
            end if
          end if
        end if
        pmeta => pmeta%next
      end do

    end do

#ifdef MPI
! redistribute blocks equally among all processors
!
    call redistribute_blocks()
#endif /* MPI */

#ifdef DEBUG
! check mesh
!
    call check_mesh('after update_mesh')
#endif /* DEBUG */

!-------------------------------------------------------------------------------
!
  end subroutine update_mesh
#ifdef MPI
!
!===============================================================================
!
! redistribute_blocks: subroutine redistributes blocks equally between
!                      processors
!
!===============================================================================
!
  subroutine redistribute_blocks()

    use blocks   , only : block_meta, block_data, list_meta, list_data
    use blocks   , only : get_nleafs, append_datablock, remove_datablock   &
                        , associate_blocks
    use coordinates, only : im, jm, km
    use mpitools , only : send_real_array, receive_real_array
    use mpitools , only : nprocs, nproc
    use equations, only : nv

    implicit none

! local variables
!
    integer         :: iret
    integer(kind=4) :: l, n
!
! tag for the MPI data exchange
!
    integer(kind=4)                            :: itag

! array for number of data block for autobalancing
!
    integer(kind=4), dimension(0:nprocs-1)      :: lb

! local buffer for data block exchange
!
    real(kind=8)   , dimension(nv,im,jm,km)   :: rbuf

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata

!-------------------------------------------------------------------------------
!
!! AUTO BALANCING
!!
! calculate the new division
!
    l       = mod(get_nleafs(), nprocs) - 1
    lb( : ) = get_nleafs() / nprocs
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
          itag = pmeta%cpu * nprocs + n + nprocs + 1

          if (nproc .eq. pmeta%cpu) then

! copy data to buffer
!
            rbuf(:,:,:,:) = pmeta%data%u(:,:,:,:)

! send data
!
            call send_real_array(size(rbuf), n, itag, rbuf, iret)

! deallocate data block
!
             call remove_datablock(pmeta%data)

! send data block
!
          end if

          if (nproc .eq. n) then

! allocate data block for the current block
!
            call append_datablock(pdata)
            call associate_blocks(pmeta, pdata)

! receive the data
!
            call receive_real_array(size(rbuf), pmeta%cpu, itag, rbuf, iret)

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
          n = min(nprocs - 1, n + 1)
          l = 0
        end if
      end if

! assign pointer to the next block
!
      pmeta => pmeta%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine redistribute_blocks
#endif /* MPI */
!
!===============================================================================
!
! prolong_block: subroutine expands the block data and copy them to children
!
!===============================================================================
!
  subroutine prolong_block(pblock)

    use blocks        , only : block_meta, block_data, nchild
    use coordinates   , only : ng, nh, in, jn, kn, im, jm, km
    use coordinates   , only : ib, ie, jb, je, kb, ke
    use equations     , only : nv
    use interpolations, only : limiter

    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k, q, p
    integer :: il, iu, jl, ju, kl, ku
    integer :: ic, jc, kc, ip, jp, kp
    real    :: dul, dur, dux, duy, duz

! local arrays
!
    integer, dimension(3) :: dm

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
    dm(:) = 2 * ((/ in, jn, kn /) + ng)
#if NDIMS == 2
    dm(3) = 1
#endif /* NDIMS == 2 */

! allocate array for the expanded arrays
!
    allocate(u(nv, dm(1), dm(2), dm(3)))

! prepare indices
!
    il = ib - nh
    iu = ie + nh
    jl = jb - nh
    ju = je + nh
#if NDIMS == 3
    kl = kb - nh
    ku = ke + nh
#endif /* NDIMS == 3 */

! expand the block variables
!
#if NDIMS == 2
    do k = 1, km
      kc = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = kl, ku
      kc = 2 * (k - kl) + 1
      kp = kc + 1
#endif /* NDIMS == 3 */
      do j = jl, ju
        jc = 2 * (j - jl) + 1
        jp = jc + 1
        do i = il, iu
          ic = 2 * (i - il) + 1
          ip = ic + 1

          do p = 1, nv

            dul = pdata%u(p,i  ,j,k) - pdata%u(p,i-1,j,k)
            dur = pdata%u(p,i+1,j,k) - pdata%u(p,i  ,j,k)
            dux = limiter(0.25d+00, dul, dur)

            dul = pdata%u(p,i,j  ,k) - pdata%u(p,i,j-1,k)
            dur = pdata%u(p,i,j+1,k) - pdata%u(p,i,j  ,k)
            duy = limiter(0.25d+00, dul, dur)

#if NDIMS == 3
            dul = pdata%u(p,i,j,k  ) - pdata%u(p,i,j,k-1)
            dur = pdata%u(p,i,j,k+1) - pdata%u(p,i,j,k  )
            duz = limiter(0.25d+00, dul, dur)
#endif /* NDIMS == 3 */

#if NDIMS == 2
            u(p,ic,jc,kc) = pdata%u(p,i,j,k) - (dux + duy)
            u(p,ip,jc,kc) = pdata%u(p,i,j,k) + (dux - duy)
            u(p,ic,jp,kc) = pdata%u(p,i,j,k) + (duy - dux)
            u(p,ip,jp,kc) = pdata%u(p,i,j,k) + (dux + duy)
#endif /* NDIMS == 2 */

#if NDIMS == 3
            u(p,ic,jc,kc) = pdata%u(p,i,j,k) - dux - duy - duz
            u(p,ip,jc,kc) = pdata%u(p,i,j,k) + dux - duy - duz
            u(p,ic,jp,kc) = pdata%u(p,i,j,k) - dux + duy - duz
            u(p,ip,jp,kc) = pdata%u(p,i,j,k) + dux + duy - duz
            u(p,ic,jc,kp) = pdata%u(p,i,j,k) - dux - duy + duz
            u(p,ip,jc,kp) = pdata%u(p,i,j,k) + dux - duy + duz
            u(p,ic,jp,kp) = pdata%u(p,i,j,k) - dux + duy + duz
            u(p,ip,jp,kp) = pdata%u(p,i,j,k) + dux + duy + duz
#endif /* NDIMS == 3 */
          end do
        end do
      end do
    end do

! iterate over all children
!
    do p = 1, nchild

! assign pointer to the current child
!
      pchild => pblock%child(p)%ptr

! obtain the position of child in the parent block
!
      ic = pchild%pos(1)
      jc = pchild%pos(2)
#if NDIMS == 3
      kc = pchild%pos(3)
#endif /* NDIMS == 3 */

! calculate indices of the current child subdomain
!
      il = 1 + ic * in
      jl = 1 + jc * jn
#if NDIMS == 3
      kl = 1 + kc * kn
#endif /* NDIMS == 3 */

      iu = il + im - 1
      ju = jl + jm - 1
#if NDIMS == 3
      ku = kl + km - 1
#endif /* NDIMS == 3 */

! copy data to the current child
!
#if NDIMS == 2
      pchild%data%u(1:nv,1:im,1:jm,1:km) = u(1:nv,il:iu,jl:ju, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      pchild%data%u(1:nv,1:im,1:jm,1:km) = u(1:nv,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

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
! subroutine RESTRICT_BLOCK:
! -------------------------
!
!   Subroutine restricts variables from children data blocks linked to
!   the input meta block and copy the resulting array of variables to
!   the associated data block.  The process of data restriction conserves
!   stored variables.
!
!   Arguments:
!
!     pblock - the input meta block;
!
!===============================================================================
!
  subroutine restrict_block(pblock)

! variables and subroutines imported from other modules
!
    use blocks     , only : ndims
    use blocks     , only : block_meta, block_data, nchild
    use coordinates, only : ng, nh, in, jn, kn, im, jm, km
    use coordinates, only : ih, jh, kh, ib, jb, kb, ie, je, ke
    use equations  , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: p
    integer :: il, jl, kl, iu, ju, ku, ip, jp, kp
    integer :: is, js, ks, it, jt, kt

! local arrays
!
    integer, dimension(ndims) :: pos

! local pointers
!
    type(block_data), pointer :: pparent, pchild

!-------------------------------------------------------------------------------
!
! assign the parent data pointer
!
    pparent => pblock%data

! iterate over all children
!
    do p = 1, nchild

! assign a pointer to the current child
!
      pchild  => pblock%child(p)%ptr%data

! obtain the child position in the parent block
!
      pos(1:ndims) = pchild%meta%pos(1:ndims)

! calculate the bound indices of the source nad destination arrays
!
      if (pos(1) == 0) then
        il =  1
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
      if (pos(2) == 0) then
        jl =  1
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
      if (pos(3) == 0) then
        kl =  1
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

! restrict conserved variables from the current child and copy the resulting
! array to the proper location of the parent data block
!
#if NDIMS == 2
      pparent%u(1:nv,is:it,js:jt, 1   ) =                                      &
                        2.50d-01 *  ((pchild%u(1:nv,il:iu:2,jl:ju:2, 1     )   &
                                 +    pchild%u(1:nv,ip:iu:2,jp:ju:2, 1     ))  &
                                 +   (pchild%u(1:nv,il:iu:2,jp:ju:2, 1     )   &
                                 +    pchild%u(1:nv,ip:iu:2,jl:ju:2, 1     )))
#endif /* NDIMS == 2 */
#if NDIMS == 3
      pparent%u(1:nv,is:it,js:jt,ks:kt) =                                      &
                        1.25d-01 * (((pchild%u(1:nv,il:iu:2,jl:ju:2,kl:ku:2)   &
                                 +    pchild%u(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))  &
                                 +   (pchild%u(1:nv,il:iu:2,jl:ju:2,kp:ku:2)   &
                                 +    pchild%u(1:nv,ip:iu:2,jp:ju:2,kl:ku:2))) &
                                 +  ((pchild%u(1:nv,il:iu:2,jp:ju:2,kp:ku:2)   &
                                 +    pchild%u(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))  &
                                 +   (pchild%u(1:nv,il:iu:2,jp:ju:2,kl:ku:2)   &
                                 +    pchild%u(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

    end do ! p = 1, nchild

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

    use error   , only : print_info
    use mpitools, only : master

    implicit none

!-------------------------------------------------------------------------------
!
! close the handler of the mesh statistics file
!
    if (master) close(funit)

!-------------------------------------------------------------------------------
!
  end subroutine clear_mesh
!
!===============================================================================
!
! store_mesh_stats: subroutine stores mesh statistics
!
!===============================================================================
!
  subroutine store_mesh_stats(n, t)

    use blocks  , only : block_meta, list_meta
    use blocks  , only : get_mblocks, get_nleafs
    use coordinates, only : ng, im, jm, km, toplev, effres
    use mpitools, only : master, nprocs

    implicit none

! arguments
!
    integer     , intent(in) :: n
    real(kind=8), intent(in) :: t

! local variables
!
    integer(kind=4) :: l
    real(kind=4)    :: cov, eff

! local arrays
!
    integer(kind=4), dimension(toplev) :: ldist
#ifdef MPI
    integer(kind=4), dimension(nprocs)  :: cdist
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta

! local variables
!
    integer(kind=4), save :: nm = 0, nl = 0

!-------------------------------------------------------------------------------
!
! store the statistics about mesh
!
    if (master) then

      if (nm .ne. get_mblocks() .or. nl .ne. get_nleafs()) then

! set new numbers of meta blocks and leafs
!
      nm = get_mblocks()
      nl = get_nleafs()

! calculate the coverage
!
      cov = (1.0 * nl) / tblocks
      eff = 1.0 * nl * (im * jm * km) / product(effres(1:NDIMS) + 2 * ng)

! get the block level distribution
!
        ldist(:) = 0
#ifdef MPI
        cdist(:) = 0
#endif /* MPI */
        pmeta => list_meta
        do while(associated(pmeta))
          if (pmeta%leaf) then
            ldist(pmeta%level) = ldist(pmeta%level) + 1
#ifdef MPI
            cdist(pmeta%cpu+1) = cdist(pmeta%cpu+1) + 1
#endif /* MPI */
          end if
          pmeta => pmeta%next
        end do

! write down the statistics
!
        write(funit,"(2x,i8,2x,1pe14.8,2(2x,i6),2(2x,1pe14.8))",advance="no")  &
                                                        n, t, nl, nm, cov, eff
        write(funit,"('   ')",advance="no")
        do l = 1, toplev
          write(funit,"(2x,i6)",advance="no") ldist(l)
        end do
#ifdef MPI
        write(funit,"('   ')",advance="no")
        do l = 1, nprocs
          write(funit,"(2x,i6)",advance="no") cdist(l)
        end do
#endif /* MPI */
        write(funit,"('')")

      end if ! number of blocks or leafs changed

    end if ! is master
!
!-------------------------------------------------------------------------------
!
  end subroutine store_mesh_stats
#ifdef DEBUG
!
!===============================================================================
!
! check_mesh: subroutine checks if the block structure is correct
!
!===============================================================================
!
  subroutine check_mesh(string)

    use blocks, only : block_meta, list_meta
    use blocks, only : check_metablock

    implicit none

! input arguments
!
    character(len=*), intent(in) :: string

! local pointers
!
    type(block_meta), pointer :: pmeta

!-------------------------------------------------------------------------------
!
! check meta blocks
!
    pmeta => list_meta
    do while(associated(pmeta))

! check the current block
!
      call check_metablock(pmeta, string)

      pmeta => pmeta%next
    end do

!-------------------------------------------------------------------------------
!
  end subroutine check_mesh
#endif /* DEBUG */

!===============================================================================
!
end module
