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
!! module: MESH - handling adaptive mesh structure
!!
!!******************************************************************************
!
module mesh

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
  integer            , save :: imi, ims, img, imu, ima, imp, imr
#endif /* PROFILE */

! file handler for the mesh statistics
!
  integer, save :: funit   = 11

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_mesh, finalize_mesh
  public :: generate_mesh, update_mesh
  public :: redistribute_blocks, store_mesh_stats

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
! subroutine INITIALIZE_MESH:
! --------------------------
!
!   Subroutine initializes mesh module and its variables.
!
!   Arguments:
!
!     nrun    - the restarted execution count;
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_mesh(nrun, verbose, iret)

! import external procedures and variables
!
    use blocks         , only : set_block_dimensions
    use coordinates    , only : xmin, xmax, ymin, ymax, zmin, zmax
    use coordinates    , only : toplev, im, jm, km
    use equations      , only : nv
    use mpitools       , only : master, nprocs

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in)    :: nrun
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64) :: fn
    integer(kind=4)   :: i, j, k, l, n

! local arrays
!
    integer(kind=4), dimension(3) :: bm, rm, dm
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('mesh:: initialization'    , imi)
    call set_timer('mesh:: statistics'        , ims)
    call set_timer('mesh:: initial generation', img)
    call set_timer('mesh:: adaptive update'   , imu)
    call set_timer('mesh:: autobalancing'     , ima)
    call set_timer('mesh:: block restriction' , imr)
    call set_timer('mesh:: block prolongation', imp)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! set data block dimensions
!
    call set_block_dimensions(nv, nv, im, jm, km)

! only master prepares the mesh statistics file
!
    if (master) then

! generate the mesh statistics file name
!
      write(fn, "('mesh_',i2.2,'.dat')") nrun

! create a new mesh statistics file
!
#ifdef INTEL
      open(unit = funit, file = fn, form = 'formatted', status = 'replace'     &
                                                           , buffered = 'yes')
#else /* INTEL */
      open(unit = funit, file = fn, form = 'formatted', status = 'replace')
#endif /* INTEL */

! write the mesh statistics header
!
      write(funit, "('#',2(4x,a4),10x,a8,6x,a10,5x,a6,6x,a5,4x,a20)")          &
                                     'step', 'time', 'coverage', 'efficiency'  &
                                   , 'blocks', 'leafs', 'block distributions:'

! write the mesh distributions header
!
      write(funit, "('#',76x,a8)", advance="no") 'level = '
      do l = 1, toplev
        write(funit, "(1x,i9)", advance="no") l
      end do
#ifdef MPI
      write(funit, "(4x,a6)", advance="no") 'cpu = '
      do n = 1, nprocs
        write(funit, "(1x,i9)", advance="no") n
      end do
#endif /* MPI */
      write(funit, "('' )")
      write(funit, "('#')")

    end if ! master

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_mesh
!
!===============================================================================
!
! subroutine FINALIZE_MESH:
! ------------------------
!
!   Subroutine releases memory used by the module variables.
!
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_mesh(iret)

! import external procedures and variables
!
    use mpitools       , only : master

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! only master posses a handler to the mesh statistics file
!
    if (master) then

! close the statistics file
!
      close(funit)

    end if ! master

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_mesh
!
!===============================================================================
!
! subroutine STORE_MESH_STATS:
! ---------------------------
!
!   Subroutine prepares and stores various mesh statistics.
!
!   Arguments:
!
!     step - the integration step;
!     time - the physical time;
!
!===============================================================================
!
  subroutine store_mesh_stats(step, time)

! import external procedures and variables
!
    use blocks         , only : ndims, block_meta, list_meta
    use blocks         , only : get_mblocks, get_nleafs
    use coordinates    , only : ng, nd, in, jn, kn, im, jm, km, ir, jr, kr, toplev
    use mpitools       , only : master, nprocs

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer     , intent(in) :: step
    real(kind=8), intent(in) :: time

! local variables
!
    integer(kind=4) :: l, n
    real(kind=8)    :: cv, ef, ff

! local pointers
!
    type(block_meta), pointer :: pmeta

! local saved variables
!
    logical        , save :: first = .true.
    integer(kind=4), save :: nm = 0, nl = 0
    real(kind=8)   , save :: fcv = 1.0d+00, fef = 1.0d+00

! local arrays
!
    integer(kind=4), dimension(toplev) :: ldist
#ifdef MPI
    integer(kind=4), dimension(nprocs) :: cdist
#endif /* MPI */

!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the mesh statistics
!
    call start_timer(ims)
#endif /* PROFILE */

! store the mesh statistics only on master
!
    if (master) then

! store the mesh statistics only if they changed
!
      if (nm /= get_mblocks() .or. nl /= get_nleafs()) then

! get the mximum block number and maximum level
!
        if (first) then

! reset the maximum number of base blocks
!
          n = 0

! set the pointer to the first block on the meta block list
!
          pmeta => list_meta

! determine the number of base blocks
!
          do while(associated(pmeta))

! if the block is at the lowest level, count it
!
            if (pmeta%level == 1) n = n + 1

! associate the pointer with the next block
!
            pmeta => pmeta%next

          end do ! pmeta

! calculate the maximum block number
!
          n = n * 2**(ndims*(toplev - 1))

! prepare coverage and efficiency factors
!
          ff  = 2**(toplev - 1)
          fcv = 1.0d+00 / n
          fef = 1.0d+00 * im / (ir * in * ff + nd)
          fef =     fef * jm / (jr * jn * ff + nd)
#if NDIMS == 3
          fef =     fef * km / (kr * kn * ff + nd)
#endif /* NDIMS == 3 */

! reset the first execution flag
!
          first = .false.

        end if ! first

! get the new numbers of meta blocks and leafs
!
        nm = get_mblocks()
        nl = get_nleafs()

! calculate the coverage (the number of leafs divided by the maximum
! block number) and the efficiency (the cells count for adaptive mesh
! divided by the cell count for corresponding uniform mesh)
!
        cv = fcv * nl
        ef = fef * nl

! initialize the level and process block counter
!
        ldist(:) = 0
#ifdef MPI
        cdist(:) = 0
#endif /* MPI */

! set the pointer to the first block on the meta block list
!
        pmeta => list_meta

! scan all meta blocks and prepare get the block level and process
! distributions
!
        do while(associated(pmeta))

! process only leafs
!
          if (pmeta%leaf) then

! increase the block level and process counts
!
            ldist(pmeta%level)     = ldist(pmeta%level)     + 1
#ifdef MPI
            cdist(pmeta%process+1) = cdist(pmeta%process+1) + 1
#endif /* MPI */

          end if ! the leaf

! associate the pointer with the next block
!
          pmeta => pmeta%next

        end do ! pmeta

! write down the block statistics
!
        write(funit, "(i9,3e14.6,2(2x,i9))", advance="no")                     &
                                                    step, time, cv, ef, nm, nl

! write down the block level distribution
!
        write(funit,"(12x)", advance="no")
        do l = 1, toplev
          write(funit,"(1x,i9)", advance="no") ldist(l)
        end do ! l = 1, toplev

#ifdef MPI
! write down the process level distribution
!
        write(funit,"(10x)", advance="no")
        do l = 1, nprocs
          write(funit,"(1x,i9)", advance="no") cdist(l)
        end do ! l = 1, nprocs
#endif /* MPI */

! write the new line symbol
!
        write(funit,"('')")

      end if ! number of blocks or leafs changed

    end if ! master

#ifdef PROFILE
! stop accounting time for the mesh statistics
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine store_mesh_stats
!
!===============================================================================
!
! subroutine GENERATE_MESH:
! ------------------------
!
!   Subroutine generates the iinitial block structure by creating blocks
!   according to the geometry, initial problems and refinement criterion.
!
!
!===============================================================================
!
  subroutine generate_mesh()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta, list_data
    use blocks         , only : ndims, nchildren, nsides, nfaces
    use blocks         , only : allocate_datablock, deallocate_datablock
    use blocks         , only : append_datablock
    use blocks         , only : link_blocks, unlink_blocks, refine_block
    use blocks         , only : get_mblocks, get_nleafs
    use blocks         , only : set_neighbors_refine
#ifdef DEBUG
    use blocks         , only : check_neighbors
#endif /* DEBUG */
    use coordinates    , only : minlev, maxlev
    use domains        , only : setup_domain
    use error          , only : print_error
    use mpitools       , only : master, nproc, nprocs
    use problems       , only : setup_problem
    use refinement     , only : check_refinement_criterion

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh, pnext
    type(block_data), pointer :: pdata

! local variables
!
    integer(kind=4)                        :: level, lev, idir, iside, iface
    integer(kind=4)                        :: np, nl
    integer(kind=4), dimension(0:nprocs-1) :: lb

!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the initial mesh generation
!
    call start_timer(img)
#endif /* PROFILE */

! allocate the initial lowest level structure of blocks
!
    call setup_domain()

! at this point we assume, that the initial structure of blocks
! according to the defined geometry is already created; no refinement
! is done yet; we fill out the coarse blocks with the initial condition
!
    if (master) then

! print the levels while generating them
!
      write(*,*)
      write(*,"(1x,a)"             ) "Generating the initial mesh:"
      write(*,"(4x,a)",advance="no") "generating level       =    "

    end if

! allocate temporary data block to determine the initial block structure
!
    call allocate_datablock(pdata)

! reset the currently processed level
!
    level = 1

! iterate over all level up to the maximum one
!
    do while (level < maxlev)

! print the currently processed level
!
      if (master) write(*, "(1x,i2)", advance = "no") level

!! DETERMINE THE REFINEMENT OF ALL BLOCKS AT CURRENT LEVEL
!!
! set the pointer to the first block on the meta block list
!
      pmeta => list_meta

! iterate over all meta blocks at the current level and initialize the problem
! for them
!
      do while (associated(pmeta))

! set the initial conditions at the current block
!
        if (pmeta%level == level) then

! associate the temporary data block with the current meta block
!
          call link_blocks(pmeta, pdata)

! set up the initial conditions for the block
!
          call setup_problem(pdata)

! check the refinement criterion for the current block, and do not allow for
! the derefinement
!
          pmeta%refine = max(0, check_refinement_criterion(pdata))

! if the level is lower then the minimum level set the block to be refined
!
          if (level < minlev) pmeta%refine = 1

! if there is only one block, refine it anyway because the resolution for
! the problem initiation may be too small
!
          if (get_mblocks() == 1 .and. level == 1) pmeta%refine = 1

! unlink the data block from the current meta block
!
          call unlink_blocks(pmeta, pdata)

        end if ! pmeta%level == 1

! assign pointer to the next block
!
        pmeta => pmeta%next

      end do ! pmeta

!! STEP DOWN AND SELECT BLOCKS WHICH NEED TO BE REFINED
!!
! walk through all levels down from the current level and check if neighbors
! of the refined blocks need to be refined as well; there is no need for
! checking the blocks at the lowest level;
!
      do lev = level, 2, -1

! set the pointer to the first block on the meta block list
!
        pmeta => list_meta

! iterate over all meta blocks at the level n and if the current block is
! selected for the refinement and its neighbors are at lower levels select them
! for refinement too;
!
        do while (associated(pmeta))

! check if the current block is a leaf at the level n and selected for
! refinement
!
          if (pmeta%leaf .and. pmeta%level == lev .and. pmeta%refine == 1) then

! select all neighbors which lay on lower levels to be refined as well
!
            call set_neighbors_refine(pmeta)

          end if ! leaf at level n and marked for refinement

! assign pointer to the next block
!
          pmeta => pmeta%next

        end do ! pmeta

      end do ! lev = level, 2, -1

!! REFINE ALL BLOCKS FROM THE LOWEST LEVEL UP
!!
! walk through the levels starting from the lowest to the current level
!
      do lev = 1, level

! set the pointer to the first block on the meta block list
!
        pmeta => list_meta

! iterate over all meta blocks
!
        do while (associated(pmeta))

! check if the current block is at the level lev and refine if if
! it is selected for refinement
!
          if (pmeta%level == lev .and. pmeta%refine == 1) then

! perform the refinement without creating new data blocks
!
            call refine_block(pmeta, .false.)

          end if ! selected for refinement and at the current level

! assign pointer to the next block
!
          pmeta => pmeta%next

        end do ! pmeta

      end do ! lev = 1, level

! increase the level number
!
      level = level + 1

    end do ! level = 1, maxlev

! deallocate temporary data block
!
    call deallocate_datablock(pdata)

! print the currently processed level
!
    if (master) write(*, "(1x,i2)") level

!! ASSOCIATE THE META BLOCKS WITH PROCESSES AND CREATE INITIAL DATA BLOCKS
!!
! divide blocks between all processes, use the number of data blocks to do this,
! but keep blocks from the top level which have the same parent packed together
!
    nl       = mod(get_nleafs(), nprocs) - 1
    lb( :  ) = get_nleafs() / nprocs
    lb(0:nl) = lb(0:nl) + 1

! reset the processor and block numbers
!
    np = 0
    nl = 0

! set the pointer to the first block on the meta block list
!
    pmeta => list_meta

! iterate over all meta blocks
!
    do while (associated(pmeta))

! assign the process number to the current block
!
      pmeta%process = np

! check if the current block is the leaf
!
      if (pmeta%leaf) then

! increase the number of leafs for the current process
!
        nl = nl + 1

! if the block belongs to the current process, append a new data block, link it
! to the current meta block and initialize the problem
!
        if (pmeta%process == nproc) then

! append new data block
!
          call append_datablock(pdata)

! link the data block with the current meta block
!
          call link_blocks(pmeta, pdata)

! set up the initial conditions for the block
!
          call setup_problem(pdata)

        end if ! leaf and belongs to the current process

! if the number of leafs for the current process exceeds the number of assigned
! blocks, reset the counter and increase the process number
!
        if (nl >= lb(np)) then

! reset the leaf counter for the current process
!
          nl = 0

! increase the process number
!
          np = min(nprocs - 1, np + 1)

        end if ! nl >= lb(np)

      end if ! leaf

! assign pointer to the next block
!
      pmeta => pmeta%next

    end do ! pmeta

#ifdef DEBUG
! check if neighbors are consistent after mesh generation
!
    call check_neighbors()
#endif /* DEBUG */

#ifdef PROFILE
! stop accounting time for the initial mesh generation
!
    call stop_timer(img)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine generate_mesh
!
!===============================================================================
!
! subroutine UPDATE_MESH:
! ----------------------
!
!   Subroutine checks the refinement criterion for each block, and refines
!   or restricts it if necessary by prolongating or restricting its data.
!   In the MPI version the data blocks are redistributed among all processes
!   after the mesh update.
!
!
!===============================================================================
!
  subroutine update_mesh()

! import external procedures and variables
!
#ifdef DEBUG
    use blocks         , only : check_neighbors
#endif /* DEBUG */

! local variables are not implicit by default
!
    implicit none

!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the adaptive mesh refinement update
!
    call start_timer(imu)
#endif /* PROFILE */

! check the refinement criterion of all data blocks at the current process
!
    call check_data_block_refinement()

! update neighbor refinement flags, if they need to be refined as well
!
    call update_neighbor_refinement()

! prepare siblings of blocks marked for restriction
!
    call prepare_sibling_derefinement()

! restrict selected blocks
!
    call derefine_selected_blocks()

! prolong selected blocks
!
    call refine_selected_blocks()

#ifdef MPI
! redistribute blocks equally among all processors
!
    call redistribute_blocks()
#endif /* MPI */

#ifdef DEBUG
! check if neighbors are consistent after mesh refinement
!
    call check_neighbors()
#endif /* DEBUG */

#ifdef PROFILE
! stop accounting time for the adaptive mesh refinement update
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_mesh
!
!===============================================================================
!
! subroutine REDISTRIBUTE_BLOCKS:
! ------------------------------
!
!   Subroutine redistributes data blocks between processes.
!
!
!===============================================================================
!
  subroutine redistribute_blocks()

#ifdef MPI
! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta, list_data
    use blocks         , only : get_nleafs
    use blocks         , only : append_datablock, remove_datablock, link_blocks
    use coordinates    , only : im, jm, km
    use equations      , only : nv
    use mpitools       , only : nprocs, nproc
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

#ifdef MPI
! local variables
!
    integer                                   :: iret
    integer(kind=4)                           :: np, nl

! local pointers
!
    type(block_meta), pointer                 :: pmeta
    type(block_data), pointer                 :: pdata

! tag for the MPI data exchange
!
    integer(kind=4)                           :: itag

! array for number of data block for autobalancing
!
    integer(kind=4), dimension(0:nprocs-1)    :: lb

! local buffer for data block exchange
!
    real(kind=8)   , dimension(2,nv,im,jm,km) :: rbuf
#endif /* MPI */

!-------------------------------------------------------------------------------
!
#ifdef MPI
#ifdef PROFILE
! start accounting time for the block redistribution
!
    call start_timer(ima)
#endif /* PROFILE */

! calculate the new data block division between processes
!
    nl       = mod(get_nleafs(), nprocs) - 1
    lb( :  ) = get_nleafs() / nprocs
    lb(0:nl) = lb(0:nl) + 1

! reset the processor and leaf numbers
!
    np = 0
    nl = 0

! set the pointer to the first block on the meta block list
!
    pmeta => list_meta

! iterate over all meta blocks and reassign their process numbers
!
    do while (associated(pmeta))

! consider only meta blocks which belong to active processes
!
      if (pmeta%process < nprocs) then

! check if the block belongs to another process
!
        if (pmeta%process /= np) then

! check if the block is the leaf
!
          if (pmeta%leaf) then

! generate a tag for communication
!
            itag = pmeta%process * nprocs + np + nprocs + 1

! sends the block to the right process
!
            if (nproc == pmeta%process) then

! copy data to buffer
!
              rbuf(1,:,:,:,:) = pmeta%data%u(:,:,:,:)
              rbuf(2,:,:,:,:) = pmeta%data%q(:,:,:,:)

! send data
!
              call send_real_array(size(rbuf), np, itag, rbuf, iret)

! remove data block from the current process
!
              call remove_datablock(pmeta%data)

! send data block
!
            end if ! nproc == pmeta%process

! receive the block from another process
!
            if (nproc == np) then

! allocate a new data block and link it with the current meta block
!
              call append_datablock(pdata)
              call link_blocks(pmeta, pdata)

! receive the data
!
              call receive_real_array(size(rbuf), pmeta%process, itag, rbuf, iret)

! coppy the buffer to data block
!
              pmeta%data%u(:,:,:,:) = rbuf(1,:,:,:,:)
              pmeta%data%q(:,:,:,:) = rbuf(2,:,:,:,:)

            end if ! nproc == n

          end if ! leaf

! set new processor number
!
          pmeta%process = np

        end if ! pmeta%process /= np

! increase the number of blocks on the current process; if it exceeds the
! allowed number reset the counter and increase the processor number
!
        if (pmeta%leaf) then

! increase the number of leafs for the current process
!
          nl = nl + 1

! if the number of leafs for the current process exceeds the number of assigned
! blocks, reset the counter and increase the process number
!
          if (nl >= lb(np)) then

! reset the leaf counter for the current process
!
            nl = 0

! increase the process number
!
            np = min(nprocs - 1, np + 1)

          end if ! l >= lb(n)

        end if ! leaf

      end if ! pmeta%process < nprocs

! assign the pointer to the next meta block
!
      pmeta => pmeta%next

    end do ! pmeta

#ifdef PROFILE
! stop accounting time for the block redistribution
!
    call stop_timer(ima)
#endif /* PROFILE */
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine redistribute_blocks
!
!===============================================================================
!
! subroutine PROLONG_BLOCK:
! ------------------------
!
!   Subroutine prolongs variables from a data blocks linked to the input
!   meta block and copy the resulting array of variables to
!   the newly created children data block.  The process of data restriction
!   conserves stored variables.
!
!   Arguments:
!
!     pblock - the input meta block;
!
!===============================================================================
!
  subroutine prolong_block(pblock)

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, nchildren
    use coordinates    , only : ng, nh, in, jn, kn, im, jm, km
    use coordinates    , only : ib, ie, jb, je, kb, ke
    use equations      , only : nv
    use interpolations , only : limiter

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_meta), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k, q, p
    integer :: il, iu, jl, ju, kl, ku
    integer :: ic, jc, kc, ip, jp, kp
    real    :: dul, dur, dux, duy, duz, du1, du2, du3, du4

! local pointers
!
    type(block_meta), pointer :: pchild
    type(block_data), pointer :: pdata

! local arrays
!
    integer, dimension(3) :: dm

! local allocatable arrays
!
    real, dimension(:,:,:,:), allocatable :: u

!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the block prolongation
!
    call start_timer(imp)
#endif /* PROFILE */

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
            du1 = dux + duy
            du2 = dux - duy
            u(p,ic,jc,kc) = pdata%u(p,i,j,k) - du1
            u(p,ip,jc,kc) = pdata%u(p,i,j,k) + du2
            u(p,ic,jp,kc) = pdata%u(p,i,j,k) - du2
            u(p,ip,jp,kc) = pdata%u(p,i,j,k) + du1
#endif /* NDIMS == 2 */

#if NDIMS == 3
            du1 = dux + duy + duz
            du2 = dux - duy - duz
            du3 = dux - duy + duz
            du4 = dux + duy - duz
            u(p,ic,jc,kc) = pdata%u(p,i,j,k) - du1
            u(p,ip,jc,kc) = pdata%u(p,i,j,k) + du2
            u(p,ic,jp,kc) = pdata%u(p,i,j,k) - du3
            u(p,ip,jp,kc) = pdata%u(p,i,j,k) + du4
            u(p,ic,jc,kp) = pdata%u(p,i,j,k) - du4
            u(p,ip,jc,kp) = pdata%u(p,i,j,k) + du3
            u(p,ic,jp,kp) = pdata%u(p,i,j,k) - du2
            u(p,ip,jp,kp) = pdata%u(p,i,j,k) + du1
#endif /* NDIMS == 3 */
          end do
        end do
      end do
    end do

! iterate over all children
!
    do p = 1, nchildren

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

    end do ! nchildren

! deallocate local arrays
!
    if (allocated(u)) deallocate(u)

#ifdef PROFILE
! stop accounting time for the block prolongation
!
    call stop_timer(imp)
#endif /* PROFILE */

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

! import external procedures and variables
!
    use blocks         , only : ndims
    use blocks         , only : block_meta, block_data, nchildren
    use coordinates    , only : ng, nh, in, jn, kn, im, jm, km
    use coordinates    , only : ih, jh, kh, ib, jb, kb, ie, je, ke
    use equations      , only : nv

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
#ifdef PROFILE
! start accounting time for the block restriction
!
    call start_timer(imr)
#endif /* PROFILE */

! assign the parent data pointer
!
    pparent => pblock%data

! iterate over all children
!
    do p = 1, nchildren

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

    end do ! p = 1, nchildren

#ifdef PROFILE
! stop accounting time for the block restriction
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine restrict_block
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine CHECK_DATA_BLOCK_REFINEMENT:
! --------------------------------------
!
!   Subroutine scans over all data blocks, gets and corrects their refinement
!   flags. If the MPI is used, the refinement flags are syncronized among all
!   processes.
!
!
!===============================================================================
!
  subroutine check_data_block_refinement()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta, list_data
#ifdef MPI
    use blocks         , only : get_nleafs
#endif /* MPI */
    use coordinates    , only : minlev, maxlev
    use error          , only : print_error
#ifdef MPI
#ifdef DEBUG
    use mpitools       , only : nproc
#endif /* DEBUG */
    use mpitools       , only : reduce_sum_integer_array
#endif /* MPI */
    use refinement     , only : check_refinement_criterion

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata

#ifdef MPI
! local variables
!
    integer(kind=4) :: nl, l
    integer         :: iret

! array for storing the refinement flags
!
    integer(kind=4), dimension(:), allocatable :: ibuf
#endif /* MPI */

!-------------------------------------------------------------------------------
!
! 1) reset the refinement flag for all meta blocks
!
! assign pmeta to the first meta block on the list
!
    pmeta => list_meta

! iterate over all meta blocks
!
    do while (associated(pmeta))

! reset the refinement flag of pmeta
!
      pmeta%refine = 0

! assign pmeta to the next meta block
!
      pmeta => pmeta%next

    end do ! iterate over meta blocks

! 2) determine the refinement of data block from the current process
!
! assign pdata to the first data block on the list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! assign pmeta to the meta block associated with pdata
!
      pmeta => pdata%meta

#ifdef DEBUG
! check if pmeta is associated
!
      if (associated(pmeta)) then

! check if pmeta is a leaf
!
        if (pmeta%leaf) then
#endif /* DEBUG */

! check the refinement criterion for the current data block
!
          pmeta%refine = check_refinement_criterion(pdata)

! correct the refinement flag for the minimum and maximum levels
!
          if (pmeta%level <  minlev) pmeta%refine =  1
          if (pmeta%level == minlev) pmeta%refine = max(0, pmeta%refine)
          if (pmeta%level == maxlev) pmeta%refine = min(0, pmeta%refine)
          if (pmeta%level >  maxlev) pmeta%refine = -1

#ifdef DEBUG
        else ! pmeta is a leaf
          call print_error("mesh::check_data_block_refinement"                 &
                                     , "Associated meta block is not a leaf!")
        end if ! pmeta is a leaf

      else ! pmeta associated
        call print_error("mesh::check_data_block_refinement"                   &
                    , "No meta block associated with the current data block!")
      end if ! pmeta associated
#endif /* DEBUG */

! assign pdata to the next data block
!
      pdata => pdata%next

   end do ! iterate over data blocks

#ifdef MPI
! 3) synchronize the refinement flags between all processes
!
! get the number of leafs
!
    nl = get_nleafs()

! allocate a buffer for the refinement flags
!
    allocate(ibuf(nl))

! check if the buffer was allocated successfully
!
    if (allocated(ibuf)) then

! reset the buffer
!
      ibuf(1:nl) = 0

! reset the leaf block counter
!
      l = 0

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! process only leafs
!
        if (pmeta%leaf) then

! increase the leaf block counter
!
          l = l + 1

! store the refinement flag in the buffer
!
          ibuf(l) = pmeta%refine

        end if ! pmeta is a leaf

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

! update refinement flags across all processes
!
      call reduce_sum_integer_array(nl, ibuf(1:nl), iret)

! reset the leaf block counter
!
      l = 0

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! process only leafs
!
        if (pmeta%leaf) then

! increase the leaf block counter
!
          l = l + 1

#ifdef DEBUG
! check if the MPI update process does not change the local refinement flags
!
          if (pmeta%process == nproc .and. pmeta%refine /= ibuf(l)) then
            call print_error("mesh::check_data_block_refinement"               &
                      , "Refinement flag does not match after MPI reduction!")
          end if
#endif /* DEBUG */

! restore the refinement flags
!
          pmeta%refine = ibuf(l)

        end if ! pmeta is a leaf

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

! deallocate the refinement flag buffer
!
      deallocate(ibuf)

    else ! buffer couldn't be allocated
      call print_error("mesh::check_data_block_refinement"                     &
                           , "Refinement flag buffer could not be allocated!")
    end if ! buffer couldn't be allocated
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine check_data_block_refinement
!
!===============================================================================
!
! subroutine UPDATE_NEIGHBOR_REFINEMENT:
! -------------------------------------
!
!   Subroutine scans over all neighbors of blocks selected for refinement or
!   derefinement, and if necessary selects them to be refined as well, or
!   cancels their derefinement.
!
!
!===============================================================================
!
  subroutine update_neighbor_refinement()

! import external procedures and variables
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : set_neighbors_refine
    use coordinates    , only : toplev

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta

! local variables
!
    integer(kind=4) :: l

!-------------------------------------------------------------------------------
!
! iterate down over all levels and correct the refinement of neighbor blocks
!
    do l = toplev, 1, -1

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! check only leafs at the current level
!
        if (pmeta%leaf .and. pmeta%level == l) then

! correct neighbor refinement flags
!
          call set_neighbors_refine(pmeta)

        end if ! the leaf at level l

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

    end do ! levels

!-------------------------------------------------------------------------------
!
  end subroutine update_neighbor_refinement
!
!===============================================================================
!
! subroutine PREPARE_SIBLING_DEREFINEMENT:
! ---------------------------------------
!
!   Subroutine scans over all blocks selected for derefinement and checks if
!   their siblings can be derefined as well. If any of them cannot be
!   derefined, the derefinement of all siblings is canceled. Then, if MPI is
!   used, the subroutine brings back all siblings together to lay on
!   the same process.
!
!   Note: This subroutine sets %refine flag of the parent to -1 to let the next
!         executed subroutine derefine_selected_blocks() which parent block
!         has to be derefined. That subroutine resets %refine flag of the
!         parent after performing full restriction.
!
!===============================================================================
!
  subroutine prepare_sibling_derefinement()

! import external procedures and variables
!
    use blocks         , only : block_meta, list_meta
#ifdef MPI
    use blocks         , only : block_data
#endif /* MPI */
    use blocks         , only : nchildren
    use blocks         , only : set_neighbors_refine
#ifdef MPI
    use blocks         , only : append_datablock, remove_datablock, link_blocks
#endif /* MPI */
    use coordinates    , only : toplev
#ifdef MPI
    use coordinates    , only : im, jm, km
    use equations      , only : nv
#endif /* MPI */
    use error          , only : print_error
#ifdef MPI
    use mpitools       , only : nprocs, nproc
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pparent, pchild
#ifdef MPI
    type(block_data), pointer :: pdata
#endif /* MPI */

! local variables
!
    logical         :: flag
    integer(kind=4) :: l, p

#ifdef MPI
! tag for the MPI data exchange
!
    integer(kind=4) :: itag
    integer         :: iret

! local buffer for data block exchange
!
    real(kind=8), dimension(nv,im,jm,km) :: rbuf
#endif /* MPI */

!-------------------------------------------------------------------------------
!
! 1) check if the siblings of the block selected for derefinement, can be
!    derefined as well, if not cancel the derefinement of all siblings
!
! iterate over levels and check sibling derefinement
!
    do l = 2, toplev

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! check only leafs at the current level
!
        if (pmeta%leaf .and. pmeta%level == l) then

! check if block is selected for derefinement
!
          if (pmeta%refine == -1) then

! assign pparent to the parent block of pmeta
!
            pparent => pmeta%parent

#ifdef DEBUG
! check if pparent is associated
!
            if (associated(pparent)) then
#endif /* DEBUG */

! reset derefinement flag
!
              flag = .true.

! iterate over all children
!
              do p = 1, nchildren

! assign pchild to the current child
!
                pchild => pparent%child(p)%ptr

#ifdef DEBUG
! check if pchild is associated
!
                if (associated(pchild)) then
#endif /* DEBUG */

! check if the current child is a leaf and selected for derefinement as well
!
                  flag = flag .and. (pchild%leaf .and. pchild%refine == -1)

#ifdef DEBUG
                else ! pchild is associated
                  call print_error("mesh::check_children_derefinement"         &
                                                 , "Children does not exist!")
                end if ! pparent is associated
#endif /* DEBUG */

              end do ! over all children

! if children can be derefined, set the refine flag of the parent to -1,
! otherwise, cancel the derefinement of all siblings
!
              if (flag) then
                pparent%refine = -1
              else

! iterate over all children
!
                do p = 1, nchildren

! assign pchild to the current child
!
                  pchild => pparent%child(p)%ptr

! reset its derefinement
!
                  pchild%refine = max(0, pchild%refine)

                end do ! children

              end if ! ~flag

#ifdef DEBUG
            else ! pparent is associated
              call print_error("mesh::check_children_derefinement"             &
                                        , "Current meta block has no parent!")
            end if ! pparent is associated
#endif /* DEBUG */

          end if ! %refine = -1

        end if ! only leafs at level l

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

    end do ! levels

#ifdef MPI
! 2) bring all siblings together to the same process
!
! assign pmeta to the first meta block on the list
!
    pmeta => list_meta

! iterate over all meta blocks
!
    do while (associated(pmeta))

! process only parent blocks (not leafs)
!
      if (.not. pmeta%leaf) then

! check if the first child is selected for derefinement
!
        if (pmeta%refine == -1) then

! assign pchild with the first child
!
          pchild => pmeta%child(1)%ptr

! set the parent process to be the same as the first child
!
          pmeta%process = pchild%process

! iterate over remaining children and if they are not on the same process,
! bring them to the parent's one
!
          do p = 2, nchildren

! assign pchild to the current child
!
            pchild => pmeta%child(p)%ptr

! if pchild belongs to a different process move its data block to the process
! of its parent
!
            if (pchild%process /= pmeta%process) then

! generate the tag for communication
!
              itag = pchild%process * nprocs + pmeta%process + nprocs + p + 1

! send data block from the current child to the parent process and deallocate it
!
              if (pchild%process == nproc) then

! assign pdata to the daba block of the current child
!
                pdata => pchild%data

#ifdef DEBUG
! check if pdata is associated
!
                if (associated(pdata)) then
#endif /* DEBUG */

! copy data to the local buffer
!
                  rbuf(:,:,:,:) = pdata%u(:,:,:,:)

#ifdef DEBUG
                else ! pdata associated
                  call print_error("mesh::check_children_derefinement"         &
                              , "Current child has no data block associated!")
                end if ! pdata associated
#endif /* DEBUG */

! send data
!
                call send_real_array(size(rbuf), pmeta%process                 &
                                                  , itag, rbuf(:,:,:,:), iret)

! deallocate the associated data block (it has to be pchild%data, and not pdata,
! otherwise, pchild%data won't be nullified)
!
                call remove_datablock(pchild%data)

              end if ! pchild%process == nproc

! allocate data block at the curent child, and receive its data from
! a different process
!
              if (pmeta%process == nproc) then

! allocate data block for the current child
!
                call append_datablock(pdata)
                call link_blocks(pchild, pdata)

! receive the data
!
                call receive_real_array(size(rbuf)                             &
                                  , pchild%process, itag, rbuf(:,:,:,:), iret)

! copy buffer to data block
!
                pdata%u(:,:,:,:) = rbuf(:,:,:,:)

              end if ! pmeta%process == nproc

! set the current processor of the block
!
              pchild%process = pmeta%process

            end if ! pchild belongs to a different process

          end do ! children

        end if ! pmeta children are selected for derefinement

      end if ! the block is parent

! assign pmeta to the next meta block
!
      pmeta => pmeta%next

    end do ! iterate over meta blocks
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine prepare_sibling_derefinement
!
!===============================================================================
!
! subroutine DEREFINE_SELECTED_BLOCKS:
! -----------------------------------
!
!   Subroutine scans over all blocks and restrict those selected.
!
!   Note: This subroutine resets the flag %refine set in subroutine
!         prepare_sibling_derefinement().
!
!===============================================================================
!
  subroutine derefine_selected_blocks()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta
    use blocks         , only : append_datablock, link_blocks, derefine_block
    use coordinates    , only : toplev
#ifdef MPI
    use mpitools       , only : nproc
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata

! local variables
!
    integer(kind=4) :: l

!-------------------------------------------------------------------------------
!
! iterate over levels and restrict the blocks selected for restriction
!
    do l = toplev - 1, 1, -1

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! process non-leafs at the current level selected for restriction
!
        if (.not. pmeta%leaf .and. pmeta%level == l                            &
                                                .and. pmeta%refine == -1) then

#ifdef MPI
! check if pmeta belongs to the current process
!
          if (pmeta%process == nproc) then
#endif /* MPI */

! check if a data block is associated with pmeta, if not create one
!
            if (.not. associated(pmeta%data)) then

! append new data block
!
              call append_datablock(pdata)

! link it with the current pmeta
!
              call link_blocks(pmeta, pdata)

            end if ! no data block associated

! perform the block restriction
!
            call restrict_block(pmeta)

#ifdef MPI
          end if ! pmeta belongs to the current process
#endif /* MPI */

! perform the mesh derefinement
!
          call derefine_block(pmeta)

! reset the refinement flag of the current block
!
          pmeta%refine = 0

        end if ! non-leaf at current level selected for derefinement

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

    end do ! levels

!-------------------------------------------------------------------------------
!
  end subroutine derefine_selected_blocks
!
!===============================================================================
!
! subroutine REFINE_SELECTED_BLOCKS:
! ---------------------------------
!
!   Subroutine scans over all blocks and prolongates those selected.
!
!
!===============================================================================
!
  subroutine refine_selected_blocks()

! import external procedures and variables
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : refine_block, remove_datablock
    use coordinates    , only : toplev
#ifdef MPI
    use mpitools       , only : nproc
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pparent

! local variables
!
    integer(kind=4) :: l

!-------------------------------------------------------------------------------
!
! iterate over all levels and prolong those selected for prolongation
!
    do l = 1, toplev - 1

! assign pmeta to the first meta block on the list
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while (associated(pmeta))

! process leafs at the current level selected for prolongation
!
        if (pmeta%leaf .and. pmeta%level == l .and. pmeta%refine == 1) then

! assign pparent with the new parent block
!
          pparent => pmeta

#ifdef MPI
! check if pmeta belongs to the current process
!
          if (pmeta%process == nproc) then
#endif /* MPI */

! prepare child blocks with allocating the data blocks
!
            call refine_block(pmeta, .true.)

! perform the data prolongation
!
            call prolong_block(pparent)

! remove the data block associated with the new parent
!
            call remove_datablock(pparent%data)

#ifdef MPI
          else ! pmeta belongs to the current process

! prepare child blocks without actually allocating the data blocks
!
            call refine_block(pmeta, .false.)

          end if ! pmeta belongs to the current process
#endif /* MPI */

        end if ! leaf at current level selected for refinement

! assign pmeta to the next meta block
!
        pmeta => pmeta%next

      end do ! iterate over meta blocks

    end do ! levels

!-------------------------------------------------------------------------------
!
  end subroutine refine_selected_blocks

!===============================================================================
!
end module
