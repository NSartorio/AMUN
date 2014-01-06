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
!! module: BOUNDARIES
!!
!!  This module handles the boundary synchronization.
!!
!!******************************************************************************
!
module boundaries

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
  integer            , save :: imi, imv, imf, ims, imc, imr, imp
#endif /* PROFILE */

! module parameters for the boundary update order and boundary type
!
  character(len = 32), save     :: xlbndry = "periodic"
  character(len = 32), save     :: xubndry = "periodic"
  character(len = 32), save     :: ylbndry = "periodic"
  character(len = 32), save     :: yubndry = "periodic"
  character(len = 32), save     :: zlbndry = "periodic"
  character(len = 32), save     :: zubndry = "periodic"

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_boundaries, finalize_boundaries
  public :: boundary_variables, boundary_fluxes
  public :: xlbndry, ylbndry, zlbndry, xubndry, yubndry, zubndry

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
! subroutine INITIALIZE_BOUNDARIES:
! --------------------------------
!
!   Subroutine initializes the module BOUNDARIES by setting its parameters.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_boundaries(verbose, iret)

! import external procedures and variables
!
#ifdef MPI
    use mpitools       , only : pdims, pcoords, periodic
#endif /* MPI */
    use parameters     , only : get_parameter_string

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
    call set_timer('boundaries:: initialization', imi)
    call set_timer('boundaries:: variables'     , imv)
    call set_timer('boundaries:: fluxes'        , imf)
    call set_timer('boundaries:: specific'      , ims)
    call set_timer('boundaries:: copy'          , imc)
    call set_timer('boundaries:: restrict'      , imr)
    call set_timer('boundaries:: prolong'       , imp)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get runtime values for the boundary types
!
    call get_parameter_string ("xlbndry" , xlbndry)
    call get_parameter_string ("xubndry" , xubndry)
    call get_parameter_string ("ylbndry" , ylbndry)
    call get_parameter_string ("yubndry" , yubndry)
    call get_parameter_string ("zlbndry" , zlbndry)
    call get_parameter_string ("zubndry" , zubndry)

! print information about the boundary conditions
!
    if (verbose) then

      write (*,"(4x,a10,13x,'=',2(1x,a))") "x-boundary"                        &
                                                , trim(xlbndry), trim(xubndry)
      write (*,"(4x,a10,13x,'=',2(1x,a))") "y-boundary"                        &
                                                , trim(ylbndry), trim(yubndry)
#if NDIMS == 3
      write (*,"(4x,a10,13x,'=',2(1x,a))") "z-boundary"                        &
                                                , trim(zlbndry), trim(zubndry)
#endif /* NDIMS == 3 */

    end if

#ifdef MPI
! change the internal boundaries to "exchange" type for the MPI update
!
    if (pdims(1) > 1) then
      if (periodic(1)) then
        xlbndry       = "exchange"
        xubndry       = "exchange"
      else
        if (pcoords(1) > 0         ) then
          xlbndry       = "exchange"
        end if
        if (pcoords(1) < pdims(1)-1) then
          xubndry       = "exchange"
        end if
      end if
    end if

    if (pdims(2) > 1) then
      if (periodic(2)) then
        ylbndry       = "exchange"
        yubndry       = "exchange"
      else
        if (pcoords(2) > 0         ) then
          ylbndry       = "exchange"
        end if
        if (pcoords(2) < pdims(2)-1) then
          yubndry       = "exchange"
        end if
      end if
    end if

    if (pdims(3) > 1) then
      if (periodic(3)) then
        zlbndry       = "exchange"
        zubndry       = "exchange"
      else
        if (pcoords(3) > 0         ) then
          zlbndry       = "exchange"
        end if
        if (pcoords(3) < pdims(3)-1) then
          zubndry       = "exchange"
        end if
      end if
    end if
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_boundaries
!
!===============================================================================
!
! subroutine FINALIZE_BOUNDARIES:
! ------------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_boundaries(iret)

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

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_boundaries
!
!===============================================================================
!
! subroutine BOUNDARY_VARIABLES:
! -----------------------------
!
!   Subroutine updates the ghost zones of the data blocks from their neighbours
!   or applies the specific boundary conditions.
!
!
!===============================================================================
!
  subroutine boundary_variables()

! import external procedures and variables
!
    use blocks         , only : ndims
    use coordinates    , only : toplev
    use mpitools       , only : periodic

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: idir, ilev
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable boundary update
!
    call start_timer(imv)
#endif /* PROFILE */

! step down from the top level
!
    do ilev = toplev, 1, -1

! iterate over all directions
!
      do idir = 1, ndims

! update boundaries which don't have neighbors and which are not periodic
!
        if (.not. periodic(idir)) call specific_boundaries(ilev, idir)

! copy boundaries between blocks at the same levels
!
        call copy_boundaries(ilev, idir)

      end do ! directions

! restrict blocks from higher level neighbours
!
      do idir = 1, ndims

        call restrict_boundaries(ilev - 1, idir)

      end do ! directions

    end do ! levels

! step up from the first level
!
    do ilev = 1, toplev

! prolong boundaries from lower level neighbours
!
      do idir = 1, ndims

        call prolong_boundaries(ilev, idir)

      end do ! boundaries

    end do ! levels

#ifdef PROFILE
! stop accounting time for variable boundary update
!
    call stop_timer(imv)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_variables
!
!===============================================================================
!
! subroutine BOUNDARY_FLUXES:
! --------------------------
!
!   Subroutine updates the numerical fluxes of neighors at different levels.
!
!
!===============================================================================
!
  subroutine boundary_fluxes()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta
#ifdef MPI
    use blocks         , only : block_info, pointer_info
#endif /* MPI */
    use blocks         , only : ndims, nsides, nfaces
    use coordinates    , only : toplev
#ifdef MPI
    use coordinates    , only : im, jm, km
#endif /* MPI */
    use coordinates    , only : ibl, ie, jbl, je, kbl, ke
#ifdef MPI
    use equations      , only : nv
    use mpitools       , only : nprocs, nproc
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: idir, iside, iface
    integer :: is, js, ks
#ifdef MPI
    integer :: irecv, isend, nblocks, itag, l, iret
#endif /* MPI */

#ifdef MPI
! local pointer arrays
!
    type(pointer_info), dimension(ndims,0:nprocs-1,0:nprocs-1) :: block_array
#endif /* MPI */

#ifdef MPI
! local arrays
!
    integer     , dimension(ndims,0:nprocs-1,0:nprocs-1) :: block_counter
    real(kind=8), dimension(:,:,:,:), allocatable        :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! do not correct fluxes if we do not use adaptive mesh
!
    if (toplev == 1) return

#ifdef PROFILE
! start accounting time for flux boundary update
!
    call start_timer(imf)
#endif /* PROFILE */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the block counter
!
    block_counter(:,:,:) = 0

! nullify pointers to blocks which need to be exchanged between processes
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1
        do idir = 1, ndims
          nullify(block_array(idir,irecv,isend)%ptr)
        end do ! idir
      end do ! isend
    end do ! irecv
#endif /* MPI */

!! 2. UPDATE THE FLUX BOUNDARIES BETWEEN THE LOCAL BLOCKS
!!
! assign the pointer to the first block on the meta list
!
    pmeta => list_meta

! scan all meta blocks in the list
!
    do while(associated(pmeta))

! check if the meta block is the leaf
!
      if (pmeta%leaf) then

! iterate over all block neighbors
!
        do idir = 1, ndims
          do iside = 1, nsides
            do iface = 1, nfaces

! associate a pointer to the current neighbor
!
              pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
              if (associated(pneigh)) then

! check if the neighbor has high level than the current block
!
                if (pmeta%level < pneigh%level) then

#ifdef MPI
! check if the block and neighbor belong to the same process, if so, update
! fluxes directly
!
                  if (pmeta%cpu == nproc .and. pneigh%cpu == nproc) then
#endif /* MPI */

! update directional flux from the neighbor
!
                    select case(idir)
                    case(1)

! prepare the boundary layer index depending on the side
!
                      if (iside == 1) then
                        is = ie
                      else
                        is = ibl
                      end if

! correct the flux from the neighor at higher level
!
                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,is,:,:), idir, iside, iface)

                    case(2)

! prepare the boundary layer index depending on the side
!
                      if (iside == 1) then
                        js = je
                      else
                        js = jbl
                      end if

! correct the flux from the neighor at higher level
!
                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,:,js,:), idir, iside, iface)

#if NDIMS == 3
                    case(3)

! prepare the boundary layer index depending on the side
!
                      if (iside == 1) then
                        ks = ke
                      else
                        ks = kbl
                      end if

! correct the flux from the neighor at higher level
!
                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,:,:,ks), idir, iside, iface)
#endif /* NDIMS == 3 */

                    end select

#ifdef MPI
! block belong to different processes, therefore prepare the block exchange
! arrays
!
                  else

! increase the counter for the number of blocks to exchange
!
                    block_counter(idir,pmeta%cpu,pneigh%cpu) =                 &
                                  block_counter(idir,pmeta%cpu,pneigh%cpu) + 1

! allocate a new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pmeta
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pmeta%level - pneigh%level

! nullify pointer fields
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! check if the list is empty
!
                    if (associated(block_array(idir,pmeta%cpu,pneigh%cpu)%ptr))&
                                                                          then
! if it is, associate the newly created block with it
!
                      pinfo%prev => block_array(idir,pmeta%cpu,pneigh%cpu)%ptr

                    end if ! %ptr associated

! point the list to the newly created block
!
                    block_array(idir,pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! pmeta and pneigh on local process
#endif /* MPI */

                end if ! pmeta level < pneigh level

              end if ! pneigh associated

            end do ! iface
          end do ! iside
        end do ! idir

      end if ! leaf

! associate the pointer with the next block
!
      pmeta => pmeta%next ! assign pointer to the next meta block in the list

    end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processes
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1
        do idir = 1, ndims

! process only pairs which have anything to exchange
!
          if (block_counter(idir,irecv,isend) > 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(idir,irecv,isend)

! prepare the tag for communication
!
            itag = (irecv * nprocs + isend) * nprocs + idir

! allocate the buffer for variables depending on the direction
!
            select case(idir)
            case(1)
              allocate(rbuf(nblocks,nv,jm,km))
            case(2)
              allocate(rbuf(nblocks,nv,im,km))
#if NDIMS == 3
            case(3)
              allocate(rbuf(nblocks,nv,im,jm))
#endif /* NDIMS == 3 */
            end select

! if isend == nproc we are sending data
!
            if (isend == nproc) then

! reset the block counter
!
              l = 0

! fill out the buffer with the data from all blocks depepnding on the direction
!
              select case(idir)
              case(1)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! prepare the ghost layer index depending on the side
!
                  if (pinfo%side == 1) then
                    is = ie
                  else
                    is = ibl
                  end if

! fill the buffer with data from the current block
!
                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,is,:,:)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks

              case(2)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! prepare the ghost layer index depending on the side
!
                  if (pinfo%side == 1) then
                    js = je
                  else
                    js = jbl
                  end if

! fill the buffer with data from the current block
!
                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,:,js,:)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks

#if NDIMS == 3
              case(3)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! prepare the ghost layer index depending on the side
!
                  if (pinfo%side == 1) then
                    ks = ke
                  else
                    ks = kbl
                  end if

! fill the buffer with data from the current block
!
                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,:,:,ks)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks
#endif /* NDIMS == 3 */

              end select

! send the data buffer to another process
!
              call send_real_array(size(rbuf(:,:,:,:)), irecv, itag            &
                                                        , rbuf(:,:,:,:), iret)

            end if ! isend = nproc

! if irecv == nproc we are receiving data
!
            if (irecv == nproc) then

! receive the data buffer
!
              call receive_real_array(size(rbuf(:,:,:,:)), isend, itag         &
                                                        , rbuf(:,:,:,:), iret)

! reset the block counter
!
              l = 0

! iterate over all received blocks and update fluxes depending on the direction
!
              select case(idir)
              case(1)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! set side and face indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! associate pointers to the meta block and neighbor
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks

              case(2)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! set side and face indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! associate pointers to the meta block and neighbor
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks

#if NDIMS == 3
              case(3)

! associate the pointer with the first block in the exchange list
!
                pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the list
!
                do while(associated(pinfo))

! increase the block count
!
                  l = l + 1

! set side and face indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! associate pointers to the meta block and neighbor
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

! associate the pointer with the next block
!
                  pinfo => pinfo%prev

                end do ! %ptr blocks
#endif /* NDIMS == 3 */

              end select

            end if ! irecv = nproc

! deallocate data buffer
!
            deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
            pinfo => block_array(idir,irecv,isend)%ptr

! scan all blocks on the exchange list
!
            do while(associated(pinfo))

! associate the exchange list pointer
!
              block_array(idir,irecv,isend)%ptr => pinfo%prev

! nullify pointer fields
!
              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

! deallocate info block
!
              deallocate(pinfo)

! associate the pointer with the next block
!
              pinfo => block_array(idir,irecv,isend)%ptr

            end do ! %ptr blocks

          end if ! if block_count > 0
        end do ! idir
      end do ! isend
    end do ! irecv
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for flux boundary update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_fluxes
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine UPDATE_CORNERS:
! -------------------------
!
!   Subroutine scans over all data blocks and updates their edges and corners.
!   This is required since the boundary update by restriction leaves the corners
!   untouched in some cases, which may result in unphysical values, like
!   negative density or pressure.  The edge/corner update should not influence
!   the solution, but just assures, that the variables are physical in all
!   cells.
!
!
!===============================================================================
!
  subroutine update_corners()

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use coordinates   , only : ib, jb, kb, ie, je, ke
    use coordinates   , only : ibl, jbl, kbl, iel, jel, kel
    use coordinates   , only : ibu, jbu, kbu, ieu, jeu, keu
    use equations     , only : nv

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: n, i, j, k

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! assign the pointer to the first block on the list
!
    pdata => list_data

! scan all data blocks until the last is reached
!
    do while(associated(pdata))

! iterate over all variables
!
      do n = 1, nv

! edges
!
#if NDIMS == 3
        do i = 1, im

          pdata%u(n,i,  1:jbl,  1:kbl) = pdata%u(n,i,jb,kb)
          pdata%u(n,i,jeu:jm ,  1:kbl) = pdata%u(n,i,je,kb)
          pdata%u(n,i,  1:jbl,keu:km ) = pdata%u(n,i,jb,ke)
          pdata%u(n,i,jeu:jm ,keu:km ) = pdata%u(n,i,je,ke)

        end do

        do j = 1, jm

          pdata%u(n,  1:ibl,j,  1:kbl) = pdata%u(n,ib,j,kb)
          pdata%u(n,ieu:im ,j,  1:kbl) = pdata%u(n,ie,j,kb)
          pdata%u(n,  1:ibl,j,keu:km ) = pdata%u(n,ib,j,ke)
          pdata%u(n,ieu:im ,j,keu:km ) = pdata%u(n,ie,j,ke)

        end do
#endif /* == 3 */

        do k = 1, km

          pdata%u(n,  1:ibl,  1:jbl,k) = pdata%u(n,ib,jb,k)
          pdata%u(n,ieu:im ,  1:jbl,k) = pdata%u(n,ie,jb,k)
          pdata%u(n,  1:ibl,jeu:jm ,k) = pdata%u(n,ib,je,k)
          pdata%u(n,ieu:im ,jeu:jm ,k) = pdata%u(n,ie,je,k)

        end do

! corners
!
#if NDIMS == 3
        pdata%u(n,  1:ibl,  1:jbl,  1:kbl) = pdata%u(n,ib,jb,kb)
        pdata%u(n,ieu:im ,  1:jbl,  1:kbl) = pdata%u(n,ie,jb,kb)
        pdata%u(n,  1:ibl,jeu:jm ,  1:kbl) = pdata%u(n,ib,je,kb)
        pdata%u(n,ieu:im ,jeu:jm ,  1:kbl) = pdata%u(n,ie,je,kb)
        pdata%u(n,  1:ibl,  1:jbl,keu:km ) = pdata%u(n,ib,jb,ke)
        pdata%u(n,ieu:im ,  1:jbl,keu:km ) = pdata%u(n,ie,jb,ke)
        pdata%u(n,  1:ibl,jeu:jm ,keu:km ) = pdata%u(n,ib,je,ke)
        pdata%u(n,ieu:im ,jeu:jm ,keu:km ) = pdata%u(n,ie,je,ke)
#endif /* == 3 */

      end do

! assign the pointer to the next block on the list
!
      pdata => pdata%next

    end do ! data blocks

!-------------------------------------------------------------------------------
!
  end subroutine update_corners
!
!===============================================================================
!
! subroutine SPECIFIC_BOUNDARIES:
! ------------------------------
!
!   Subroutine scans over all leaf blocks in order to find blocks without
!   neighbours, then updates the boundaries for selected type.
!
!   Arguments:
!
!     ilev - the level to be processed;
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine specific_boundaries(ilev, idir)

! import external procedures and variables
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : nsides
#ifdef MPI
    use mpitools       , only : nproc
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in)       :: ilev, idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh

! local variables
!
    integer                   :: iside
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for specific boundary update
!
    call start_timer(ims)
#endif /* PROFILE */

! assign the pointer with the first block on the meta list
!
    pmeta => list_meta

! scan all data blocks until the last is reached
!
    do while(associated(pmeta))

! check if the current meta block is a leaf
!
      if (pmeta%leaf .and. pmeta%level == ilev) then

#ifdef MPI
! check if the current block belongs to the local process
!
        if (pmeta%cpu == nproc) then
#endif /* MPI */

! iterate over all neighbors
!
          do iside = 1, nsides

! assign a neighbour pointer to the current neighbour
!
            pneigh => pmeta%neigh(idir,iside,1)%ptr

! make sure that the neighbour is not associated, then apply specific boundaries
!
            if (.not. associated(pneigh))                                      &
                               call boundary_specific(pmeta%data, idir, iside)

          end do ! sides

#ifdef MPI
        end if ! block belong to the local process
#endif /* MPI */
      end if ! leaf

! assign the pointer to the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef PROFILE
! stop accounting time for specific boundary update
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine specific_boundaries
!
!===============================================================================
!
! subroutine COPY_BOUNDARIES:
! --------------------------
!
!   Subroutine scans over all leaf blocks in order to find neighbours at
!   the same levels, then updates the boundaries between neighbours.
!
!
!===============================================================================
!
  subroutine copy_boundaries(ilev, idir)

! include external procedures
!
#ifdef MPI
    use mpitools      , only : send_real_array, receive_real_array
#endif /* MPI */

! include external variables
!
    use blocks        , only : ndims, nsides, nfaces
    use blocks        , only : block_meta, block_data, list_meta
    use blocks        , only : block_info, pointer_info
    use coordinates   , only : toplev
    use coordinates   , only : ng, nd, nh, im, jm, km
    use coordinates   , only : ib, jb, kb, ie, je, ke
    use coordinates   , only : ibu, jbu, kbu, iel, jel, kel
    use mpitools      , only : periodic
#ifdef MPI
    use mpitools      , only : nproc, nprocs, npmax
    use equations     , only : nv
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: ilev, idir

! local variables
!
    integer :: iside, iface, nside, nface
    integer :: iret
    integer :: il, jl, kl, iu, ju, ku
#ifdef MPI
    integer :: isend, irecv, nblocks, itag, l

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)        :: block_counter
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_data), pointer :: pdata
#ifdef MPI
    type(block_info), pointer :: pinfo
    type(pointer_info), dimension(0:npmax,0:npmax)  :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for copy boundary update
!
    call start_timer(imc)
#endif /* PROFILE */

#ifdef MPI
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(irecv,isend)%ptr)
      end do
    end do
#endif /* MPI */

! assign the pointer with the first block on in the list
!
    pmeta => list_meta

! scan all meta blocks and process blocks at the current level
!
    do while(associated(pmeta))

! check if the block is a leaf at the current level
!
      if (pmeta%leaf .and. pmeta%level .eq. ilev) then

! scan over sides and faces
!
        do iside = 1, nsides
          do iface = 1, nfaces

! assign a pointer to the neighbor
!
            pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at the same level
!
              if (pneigh%level .eq. pmeta%level) then

! copy blocks only for the first face
!
                if (iface .eq. 1) then

#ifdef MPI
! check if the current meta block and its neighbor lay on the same processor
!
                  if (pmeta%cpu .eq. pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                    if (pmeta%cpu .eq. nproc) then
#endif /* MPI */

! assign a pointer to the data structure of the current block
!
                      pdata  => pmeta%data

! update the boundaries of the current block
!
                      select case(idir)
                      case(1)
                        if (iside .eq. 1) then
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,iel:ie,:,:), idir, iside)
                        else
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,ib:ibu,:,:), idir, iside)
                        end if
                      case(2)
                        if (iside .eq. 1) then
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,:,jel:je,:), idir, iside)
                        else
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,:,jb:jbu,:), idir, iside)
                        end if
#if NDIMS == 3
                      case(3)
                        if (iside .eq. 1) then
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,:,:,kel:ke), idir, iside)
                        else
                          call boundary_copy(pdata                         &
                               , pneigh%data%u(:,:,:,kb:kbu), idir, iside)
                        end if
#endif /* NDIMS == 3 */
                      end select

#ifdef MPI
                    end if ! pmeta on the current cpu

                  else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                    block_counter(pmeta%cpu,pneigh%cpu) =                      &
                                       block_counter(pmeta%cpu,pneigh%cpu) + 1

! allocate a new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pmeta
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pmeta%level - pneigh%level

! nullify pointers
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! if the list is not empty append the created block
!
                    if (associated(block_array(pmeta%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(pmeta%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! block and neighbor on different processors
#endif /* MPI */

                end if ! iface = 1

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! faces
        end do ! sides

      end if ! leaf

! associate the pointer with the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
    do irecv = 0, npmax
      do isend = 0, npmax

! process only pairs which have boundaries to exchange
!
        if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
          nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
          itag = 10 * (irecv * nprocs + isend + 1) + 4

! allocate space for variables
!
          select case(idir)
          case(1)
            allocate(rbuf(nblocks,nv,ng,jm,km))
          case(2)
            allocate(rbuf(nblocks,nv,im,ng,km))
#if NDIMS == 3
          case(3)
            allocate(rbuf(nblocks,nv,im,jm,ng))
#endif /* NDIMS == 3 */
          end select

! if isend == nproc we are sending data
!
          if (isend .eq. nproc) then

! iterate over exchange blocks along the current direction and fill out
! the buffer with the block data
!
            select case(idir)
            case(1)
              l = 1
              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                if (pinfo%side .eq. 1) then
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,iel:ie,:,:)
                else
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,ib:ibu,:,:)
                end if

                pinfo => pinfo%prev
                l = l + 1
              end do

            case(2)
              l = 1
              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                if (pinfo%side .eq. 1) then
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,jel:je,:)
                else
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,jb:jbu,:)
                end if

                pinfo => pinfo%prev
                l = l + 1
              end do

#if NDIMS == 3
            case(3)
              l = 1
              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                if (pinfo%side .eq. 1) then
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,kel:ke)
                else
                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,kb:kbu)
                end if

                pinfo => pinfo%prev
                l = l + 1
              end do
#endif /* NDIMS == 3 */
            end select

! send the data buffer
!
            call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

          end if ! isend = nproc

! if irecv == nproc we are receiving data
!
          if (irecv .eq. nproc) then

! receive data
!
            call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:), iret)

! iterate over all received blocks and update boundaries
!
            l = 1
            pinfo => block_array(irecv,isend)%ptr
            do while(associated(pinfo))

! set indices
!
              iside =  pinfo%side

! assign a pointer to the data structure of the current block
!
              pdata => pinfo%block%data

! update the boundaries of the current block
!
              call boundary_copy(pdata, rbuf(l,:,:,:,:), idir, iside)

              pinfo => pinfo%prev
              l = l + 1
            end do

          end if ! irecv = nproc

! deallocate buffers
!
          if (allocated(rbuf)) deallocate(rbuf)

! deallocate info blocks
!
          pinfo => block_array(irecv,isend)%ptr
          do while(associated(pinfo))
            block_array(irecv,isend)%ptr => pinfo%prev

            nullify(pinfo%prev)
            nullify(pinfo%next)
            nullify(pinfo%block)
            nullify(pinfo%neigh)

            deallocate(pinfo)

            pinfo => block_array(irecv,isend)%ptr
          end do

        end if ! if block_count > 0

      end do ! isend
    end do ! irecv
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for copy boundary update
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine copy_boundaries
!
!===============================================================================
!
! subroutine RESTRICT_BOUNDARIES:
! ------------------------------
!
!   Subroutine scans over all leaf blocks in order to find neighbours at
!   different levels, then updates the boundaries of blocks at lower levels by
!   restricting variables from higher level blocks.
!
!
!===============================================================================
!
  subroutine restrict_boundaries(ilev, idir)

! include external procedures
!
#ifdef MPI
    use mpitools      , only : send_real_array, receive_real_array
#endif /* MPI */

! include external variables
!
    use blocks        , only : ndims, nsides, nfaces
    use blocks        , only : block_meta, block_data, list_meta
    use blocks        , only : block_info, pointer_info
    use coordinates   , only : toplev
    use coordinates   , only : ng, nd, nh, im, jm, km
    use coordinates   , only : ib, jb, kb, ie, je, ke
    use coordinates   , only : ibu, jbu, kbu, iel, jel, kel
    use mpitools      , only : periodic
#ifdef MPI
    use mpitools      , only : nproc, nprocs, npmax
    use equations     , only : nv
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: ilev, idir

! local variables
!
    integer :: iside, iface, nside, nface, level
    integer :: iret
    integer :: il, jl, kl, iu, ju, ku
#ifdef MPI
    integer :: isend, irecv, nblocks, itag, l

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)        :: block_counter
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_data), pointer :: pdata
#ifdef MPI
    type(block_info), pointer :: pinfo
    type(pointer_info), dimension(0:npmax,0:npmax)  :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for restrict boundary update
!
    call start_timer(imr)
#endif /* PROFILE */

#ifdef MPI
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1
        nullify(block_array(irecv,isend)%ptr)
      end do
    end do
#endif /* MPI */

! assign the pointer to the first meta block on in the list
!
    pmeta => list_meta

! scan all meta blocks until the last is reached
!
    do while(associated(pmeta))

! process only leafs from the current level
!
      if (pmeta%leaf .and. pmeta%level == ilev) then

! scan over all directions, sides and faces
!
!         do idir = 1, ndims
          do iside = 1, nsides
            do iface = 1, nfaces

! assign a neighbour pointer to the neighbour
!
              pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbour is associated
!
              if (associated(pneigh)) then

! if the neighbour is at the higher level
!
                if (pmeta%level < pneigh%level) then

#ifdef MPI
! check if the current meta block and its neighbour belong to the same processor
!
                  if (pmeta%cpu == pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                    if (pmeta%cpu == nproc) then
#endif /* MPI */

! prepare indices of the neighbour array
!
                      select case(idir)
                        case(1)
                          if (iside .eq. 1) then
                            il = ie - nd + 1
                            iu = ie
                          else
                            il = ib
                            iu = ib + nd - 1
                          end if
                          jl = 1
                          ju = jm
                          kl = 1
                          ku = km

                        case(2)
                          if (iside .eq. 1) then
                            jl = je - nd + 1
                            ju = je
                          else
                            jl = jb
                            ju = jb + nd - 1
                          end if
                          il = 1
                          iu = im
                          kl = 1
                          ku = km

#if NDIMS == 3
                        case(3)
                          if (iside .eq. 1) then
                            kl = ke - nd + 1
                            ku = ke
                          else
                            kl = kb
                            ku = kb + nd - 1
                          end if
                          il = 1
                          iu = im
                          jl = 1
                          ju = jm
#endif /* NDIMS == 3 */
                      end select

! assign a pointer to the data structure of the current block
!
                      pdata  => pmeta%data

! update the boundaries of the current block
!
                      call boundary_restrict(pdata                             &
                                         , pneigh%data%u(:,il:iu,jl:ju,kl:ku)  &
                                                        , idir, iside, iface)

#ifdef MPI
                    end if ! block on the current processor

                  else ! block and neighbour on different processors

! increase the counter for number of blocks to exchange
!
                      block_counter(pmeta%cpu,pneigh%cpu) =                    &
                                       block_counter(pmeta%cpu,pneigh%cpu) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out its fields
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%direction        =  idir
                      pinfo%side             =  iside
                      pinfo%face             =  iface
                      pinfo%level_difference =  pmeta%level - pneigh%level

! nullify pointers
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the created block
!
                      if (associated(block_array(pmeta%cpu,pneigh%cpu)%ptr)) then
                        pinfo%prev => block_array(pmeta%cpu,pneigh%cpu)%ptr
                        nullify(pinfo%next)
                      end if

! point the list to the last created block
!
                      block_array(pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! block and neighbour on different processors
#endif /* MPI */

                end if ! block at lower level than neighbour

              end if ! neighbour associated

            end do ! faces
          end do ! sides
!         end do ! directions

      end if ! leaf

! assign the pointer to the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
        do irecv = 0, npmax
          do isend = 0, npmax

!! process blocks with the neighbors at higher levels
!!
! process only pairs which have boundaries to exchange
!
            if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
              nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
              itag = 10 * (irecv * nprocs + isend + 1) + 2

! allocate space for variables
!
              select case(idir)
              case(1)
                allocate(rbuf(nblocks,nv,nd,jm,km))
              case(2)
                allocate(rbuf(nblocks,nv,im,nd,km))
              case(3)
                allocate(rbuf(nblocks,nv,im,jm,nd))
              end select

! if isend == nproc we are sending data
!
              if (isend .eq. nproc) then

! fill out the buffer with block data
!
                l = 1

                pinfo => block_array(irecv,isend)%ptr
                do while(associated(pinfo))

! prepare indices of the neighbor array
!
                  select case(idir)
                  case(1)
                    if (pinfo%side .eq. 1) then
                      il = ie - nd + 1
                      iu = ie
                    else
                      il = ib
                      iu = ib + nd - 1
                    end if

                    rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,il:iu,:,:)

                  case(2)
                    if (pinfo%side .eq. 1) then
                      jl = je - nd + 1
                      ju = je
                    else
                      jl = jb
                      ju = jb + nd - 1
                    end if

                    rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,jl:ju,:)

                  case(3)
                    if (pinfo%side .eq. 1) then
                      kl = ke - nd + 1
                      ku = ke
                    else
                      kl = kb
                      ku = kb + nd - 1
                    end if

                    rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,kl:ku)
                  end select

                  pinfo => pinfo%prev
                  l = l + 1
                end do

! send data buffer
!
                call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

              end if

! if irecv == nproc we are receiving data
!
              if (irecv .eq. nproc) then

! receive data
!
                call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:), iret)

! iterate over all received blocks and update boundaries
!
                l = 1

                pinfo => block_array(irecv,isend)%ptr
                do while(associated(pinfo))

! set indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! assign a pointer to the data structure of the current block
!
                  pdata => pinfo%block%data

! update the boundaries of the current block
!
                  call boundary_restrict(pdata, rbuf(l,:,:,:,:)                &
                                                         , idir, iside, iface)

                  pinfo => pinfo%prev
                  l = l + 1
                end do

              end if

! deallocate buffers
!
              if (allocated(rbuf)) deallocate(rbuf)

! deallocate info blocks
!
              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))
                block_array(irecv,isend)%ptr => pinfo%prev

                nullify(pinfo%prev)
                nullify(pinfo%next)
                nullify(pinfo%block)
                nullify(pinfo%neigh)

                deallocate(pinfo)

                pinfo => block_array(irecv,isend)%ptr
              end do

            end if ! if block_count > 0

          end do ! isend
        end do ! irecv
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for restrict boundary update
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine restrict_boundaries
!
!===============================================================================
!
! subroutine PROLONG_BOUNDARIES:
! -----------------------------
!
!   Subroutine scans over all leaf blocks in order to find neighbours at
!   different levels, then updates the boundaries of blocks at higher level by
!   prolongating variables from lower level blocks.
!
!
!===============================================================================
!
  subroutine prolong_boundaries(ilev, idir)

! include external procedures
!
#ifdef MPI
    use mpitools      , only : send_real_array, receive_real_array
#endif /* MPI */

! include external variables
!
    use blocks        , only : ndims, nsides, nfaces
    use blocks        , only : block_meta, block_data, list_meta
    use blocks        , only : block_info, pointer_info
    use coordinates   , only : toplev
    use coordinates   , only : ng, nd, nh, im, jm, km
    use coordinates   , only : ib, jb, kb, ie, je, ke
    use coordinates   , only : ibu, jbu, kbu, iel, jel, kel
    use mpitools      , only : periodic
#ifdef MPI
    use mpitools      , only : nproc, nprocs, npmax
    use equations     , only : nv
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: ilev, idir

! local variables
!
    integer :: iside, iface, nside, nface
    integer :: iret
    integer :: il, jl, kl, iu, ju, ku
#ifdef MPI
    integer :: isend, irecv, nblocks, itag, l

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)        :: block_counter
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_data), pointer :: pdata
#ifdef MPI
    type(block_info), pointer :: pinfo
    type(pointer_info), dimension(0:npmax,0:npmax)  :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for prolong boundary update
!
    call start_timer(imp)
#endif /* PROFILE */

#ifdef MPI
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(irecv,isend)%ptr)
      end do
    end do
#endif /* MPI */

! assign the pointer with the first block on in the list
!
    pmeta => list_meta

! scan all meta blocks and process blocks at the current level
!
    do while(associated(pmeta))

! check if the block is a leaf at the current level
!
      if (pmeta%leaf .and. pmeta%level .eq. ilev) then

! scan over sides and faces
!
        do iside = 1, nsides
          do iface = 1, nfaces

! assign a pointer to the neighbor
!
            pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at lower level
!
              if (pneigh%level .lt. pmeta%level) then

! perform update only for the first face, since all faces point the same block
!
                if (iface .eq. 1) then

#ifdef MPI
! check if the current meta block and its neighbor lay on the same processor
!
                  if (pmeta%cpu .eq. pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                    if (pmeta%cpu .eq. nproc) then
#endif /* MPI */

! find the face of the current block which the neighbor belongs to
!
                      nside = 3 - iside
                      nface = 1
                      do while(pmeta%id .ne.                               &
                                    pneigh%neigh(idir,nside,nface)%ptr%id)
                        nface = nface + 1
                      end do

! prepare indices of the neighbor array
!
                      il = 1
                      iu = im
                      jl = 1
                      ju = jm
                      kl = 1
                      ku = km

                      select case(idir)
                      case(1)
                        if (iside .eq. 1) then
                          il = ie - nh
                          iu = ie + 1
                        else
                          il = ib - 1
                          iu = ib + nh
                        end if
                      case(2)
                        if (iside .eq. 1) then
                          jl = je - nh
                          ju = je + 1
                        else
                          jl = jb - 1
                          ju = jb + nh
                        end if
                      case(3)
                        if (iside .eq. 1) then
                          kl = ke - nh
                          ku = ke + 1
                        else
                          kl = kb - 1
                          ku = kb + nh
                        end if
                      end select

! assign a pointer to the data structure of the current block
!
                      pdata  => pmeta%data

! update the boundaries of the current block
!
                      call boundary_prolong(pdata                              &
                                         , pneigh%data%u(:,il:iu,jl:ju,kl:ku)  &
                                                         , idir, iside, nface)

#ifdef MPI
                    end if ! pmeta on the current cpu

                  else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                    block_counter(pmeta%cpu,pneigh%cpu) =                      &
                                       block_counter(pmeta%cpu,pneigh%cpu) + 1

! allocate a new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pmeta
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pmeta%level - pneigh%level

! nullify pointers
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! if the list is not empty append the created block
!
                    if (associated(block_array(pmeta%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(pmeta%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! block and neighbor on different processors
#endif /* MPI */

                end if ! iface = 1

              end if ! neighbor at lower level

            end if ! neighbor associated

          end do ! faces
        end do ! sides

      end if ! leaf

! associate the pointer with the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1

! process only pairs which have boundaries to exchange
!
        if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
          nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
          itag = 10 * (irecv * nprocs + isend + 1) + 3

! allocate space for variables
!
          select case(idir)
          case(1)
            allocate(rbuf(nblocks,nv,nh+2,jm,km))
          case(2)
            allocate(rbuf(nblocks,nv,im,nh+2,km))
          case(3)
            allocate(rbuf(nblocks,nv,im,jm,nh+2))
          end select

! if isend == nproc we are sending data
!
          if (isend .eq. nproc) then

! fill out the buffer with block data
!
            l = 1

            pinfo => block_array(irecv,isend)%ptr
            do while(associated(pinfo))

! prepare indices of the neighbor array
!
              select case(idir)
              case(1)
                if (pinfo%side .eq. 1) then
                  il = ie - nh
                  iu = ie + 1
                else
                  il = ib - 1
                  iu = ib + nh
                end if

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,il:iu,:,:)

              case(2)
                if (pinfo%side .eq. 1) then
                  jl = je - nh
                  ju = je + 1
                else
                  jl = jb - 1
                  ju = jb + nh
                end if

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,jl:ju,:)

              case(3)
                if (pinfo%side .eq. 1) then
                  kl = ke - nh
                  ku = ke + 1
                else
                  kl = kb - 1
                  ku = kb + nh
                end if

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,kl:ku)
              end select

              pinfo => pinfo%prev
              l = l + 1
            end do

! send data buffer
!
            call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

          end if

! if irecv == nproc we are receiving data
!
          if (irecv .eq. nproc) then

! receive data
!
            call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:), iret)

! iterate over all received blocks and update boundaries
!
            l = 1

            pinfo => block_array(irecv,isend)%ptr
            do while(associated(pinfo))

! set indices
!
              iside = pinfo%side
              iface = pinfo%face

! assign pointers to the meta, data and neighbor blocks
!
              pmeta  => pinfo%block
              pdata  => pinfo%block%data
              pneigh => pmeta%neigh(idir,iside,iface)%ptr

! find the face of the current block which the neighbor belongs to
!
              nside = 3 - iside
              nface = 1
              do while(pmeta%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                nface = nface + 1
              end do

! update the boundaries of the current block
!
              call boundary_prolong(pdata, rbuf(l,:,:,:,:)                 &
                                                     , idir, iside, nface)

              pinfo => pinfo%prev
              l = l + 1
            end do

          end if

! deallocate buffers
!
          if (allocated(rbuf)) deallocate(rbuf)

! deallocate info blocks
!
          pinfo => block_array(irecv,isend)%ptr
          do while(associated(pinfo))
            block_array(irecv,isend)%ptr => pinfo%prev

            nullify(pinfo%prev)
            nullify(pinfo%next)
            nullify(pinfo%block)
            nullify(pinfo%neigh)

            deallocate(pinfo)

            pinfo => block_array(irecv,isend)%ptr
          end do

        end if ! if block_count > 0

      end do ! isend
    end do ! irecv
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for prolong boundary update
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prolong_boundaries
!
!===============================================================================
!
! correct_flux: subroutine copies the boundary flux from the neighbor at higher
!               level and updates its own
!
!===============================================================================
!
  subroutine correct_flux(pdata, f, idir, iside, iface)

    use blocks   , only : block_data
    use coordinates, only : ng, in, jn, kn, ih, jh, kh                         &
                          , ib, ie, ibl, jb, je, jbl, kb, ke, kbl

    implicit none

! local variables
!
    integer :: i, ic, it, il, iu, i1, i2
    integer :: j, jc, jt, jl, ju, j1, j2
#if NDIMS == 3
    integer :: k, kc, kt, kl, ku, k1, k2
#endif /* NDIMS == 3 */

! arguments
!
    type(block_data), pointer         , intent(inout) :: pdata
    real            , dimension(:,:,:), intent(in)    :: f
    integer                           , intent(in)    :: idir, iside, iface
!
!-------------------------------------------------------------------------------
!
    select case(idir)

! X direction
!
    case(1)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        it = ibl
      else                   ! right side
        it = ie
      end if

! convert face number to index
!
      jc = mod(iface - 1, 2)
#if NDIMS == 3
      kc =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! bounds for the perpendicular update
!
      jl = jb + (jh - ng) * jc
      ju = jh + (jh - ng) * jc
#if NDIMS == 3
      kl = kb + (kh - ng) * kc
      ku = kh + (kh - ng) * kc
#endif /* NDIMS == 3 */

! iterate over perpendicular direction
!
      do j = jl, ju
        j1 = 2 * (j - jl) + jb
        j2 = j1 + 1

#if NDIMS == 2
        pdata%f(idir,:,it,j,:) = 0.5d0 * (f(:,j1,:) + f(:,j2,:))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,it,j,k) = 0.25d0 * (f(:,j1,k1) + f(:,j2,k1)           &
                                           + f(:,j1,k2) + f(:,j2,k2))
        end do
#endif /* NDIMS == 3 */
      end do

! Y direction
!
    case(2)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        jt = jbl
      else                   ! right side
        jt = je
      end if

! convert face number to index
!
      ic = mod(iface - 1, 2)
#if NDIMS == 3
      kc =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! bounds for the perpendicular update
!
      il = ib + (ih - ng) * ic
      iu = ih + (ih - ng) * ic
#if NDIMS == 3
      kl = kb + (kh - ng) * kc
      ku = kh + (kh - ng) * kc
#endif /* NDIMS == 3 */

! iterate over perpendicular direction
!
      do i = il, iu
        i1 = 2 * (i - il) + ib
        i2 = i1 + 1

#if NDIMS == 2
        pdata%f(idir,:,i,jt,:) = 0.5d0 * (f(:,i1,:) + f(:,i2,:))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,i,jt,k) = 0.25d0 * (f(:,i1,k1) + f(:,i2,k1)           &
                                           + f(:,i1,k2) + f(:,i2,k2))
        end do
#endif /* NDIMS == 3 */
      end do

#if NDIMS == 3
! Z direction
!
    case(3)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        kt = kbl
      else                   ! right side
        kt = ke
      end if

! convert face number to index
!
      ic = mod(iface - 1, 2)
      jc =    (iface - 1) / 2

! bounds for the perpendicular update
!
      il = ib + (ih - ng) * ic
      iu = ih + (ih - ng) * ic
      jl = jb + (jh - ng) * jc
      ju = jh + (jh - ng) * jc

! iterate over perpendicular direction
!
      do i = il, iu
        i1 = 2 * (i - il) + ib
        i2 = i1 + 1

        do j = jl, ju
          j1 = 2 * (j - jl) + jb
          j2 = j1 + 1

          pdata%f(idir,:,i,j,kt) = 0.25d0 * (f(:,i1,j1) + f(:,i2,j1)           &
                                           + f(:,i1,j2) + f(:,i2,j2))
        end do
      end do
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine correct_flux
!
!===============================================================================
!
! boundary_copy: subroutine copies the interior of neighbor to update
!                a boundary of the current block
!
!===============================================================================
!
  subroutine boundary_copy(pdata, u, idir, iside)

    use blocks   , only : block_data
    use coordinates, only : ng, im, jm, km, ibl, ieu, jbl, jeu, kbl, keu
    use equations, only : nv

    implicit none

! input arguments
!
    type(block_data), pointer  , intent(inout) :: pdata
    real   , dimension(:,:,:,:), intent(in)    :: u
    integer                    , intent(in)    :: idir, iside
!
!-------------------------------------------------------------------------------
!
    select case(idir)

    case(1)

      if (iside .eq. 1) then
        pdata%u(1:nv,  1:ibl,1:jm,1:km) = u(1:nv,1:ng,1:jm,1:km)
      else
        pdata%u(1:nv,ieu:im ,1:jm,1:km) = u(1:nv,1:ng,1:jm,1:km)
      end if

    case(2)

      if (iside .eq. 1) then
        pdata%u(1:nv,1:im,  1:jbl,1:km) = u(1:nv,1:im,1:ng,1:km)
      else
        pdata%u(1:nv,1:im,jeu:jm ,1:km) = u(1:nv,1:im,1:ng,1:km)
      end if

#if NDIMS == 3
    case(3)

      if (iside .eq. 1) then
        pdata%u(1:nv,1:im,1:jm,  1:kbl) = u(1:nv,1:im,1:jm,1:ng)
      else
        pdata%u(1:nv,1:im,1:jm,keu:km ) = u(1:nv,1:im,1:jm,1:ng)
      end if
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine boundary_copy
!
!===============================================================================
!
! subroutine BOUNDARY_RESTRICT:
! ----------------------------
!
!   Subroutine restricts variables from the interior of the neighbor data
!   block copies the resulting array of variables to the fills out
!   the proper range of ghost zones.  The process of data restriction
!   conserves stored variables.
!
!   Arguments:
!
!     pdata              - the input data block;
!     u                  - the conserved array obtained from the neighbor;
!     idir, iside, iface - the positions of the neighbor block;
!
!===============================================================================
!
  subroutine boundary_restrict(pdata, u, idir, iside, iface)

! variables and subroutines imported from other modules
!
    use blocks       , only : block_data
    use coordinates  , only : ng, im, ih, ib, ie, ieu           &
                            , nd, jm, jh, jb, je, jeu           &
                            , nh, km, kh, kb, ke, keu
    use equations    , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data)           , pointer, intent(inout) :: pdata
    real(kind=8)    , dimension(:,:,:,:), intent(in)    :: u
    integer                             , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: ic, jc, kc, ip, jp, kp
    integer :: il, jl, kl, iu, ju, ku
    integer :: is, js, ks, it, jt, kt
!
!-------------------------------------------------------------------------------
!
! prepare indices
!
    select case(idir)

    case(1)

! X indices
!
      if (iside == 1) then
        is =  1
        it = ng
      else
        is = ieu
        it = im
      end if

      il =  1
      iu = nd
      ip = il + 1

! Y indices
!
      jc = mod(iface - 1, 2)

      js = jb - nh + (jh - nh) * jc
      jt = jh      + (jh - nh) * jc

      jl =  1 + ng * jc
      ju = je + ng * jc
      jp = jl + 1

#if NDIMS == 3
! Z indices
!
      kc = (iface - 1) / 2

      ks = kb - nh + (kh - nh) * kc
      kt = kh      + (kh - nh) * kc

      kl =  1 + ng * kc
      ku = ke + ng * kc
      kp = kl + 1
#endif /* NDIMS == 3 */

    case(2)

! X indices
!
      ic = mod(iface - 1, 2)

      is = ib - nh + (ih - nh) * ic
      it = ih      + (ih - nh) * ic

      il =  1 + ng * ic
      iu = ie + ng * ic
      ip = il + 1

! Y indices
!
      if (iside == 1) then
        js =  1
        jt = ng
      else
        js = jeu
        jt = jm
      end if

      jl =  1
      ju = nd
      jp = jl + 1

#if NDIMS == 3
! Z indices
!
      kc = (iface - 1) / 2

      ks = kb - nh + (kh - nh) * kc
      kt = kh      + (kh - nh) * kc

      kl =  1 + ng * kc
      ku = ke + ng * kc
      kp = kl + 1
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

! X indices
!
      ic = mod(iface - 1, 2)

      is = ib - nh + (ih - nh) * ic
      it = ih      + (ih - nh) * ic

      il =  1 + ng * ic
      iu = ie + ng * ic
      ip = il + 1

! Y indices
!
      jc = (iface - 1) / 2

      js = jb - nh + (jh - nh) * jc
      jt = jh      + (jh - nh) * jc

      jl =  1 + ng * jc
      ju = je + ng * jc
      jp = jl + 1

! Z indices
!
      if (iside == 1) then
        ks =  1
        kt = ng
      else
        ks = keu
        kt = km
      end if

      kl =  1
      ku = nd
      kp = kl + 1
#endif /* NDIMS == 3 */

    end select

! update boundaries of the conserved variables
!
#if NDIMS == 2
    pdata%u(:,is:it,js:jt, 1   ) =                                             &
                               2.50d-01 *  ((u(1:nv,il:iu:2,jl:ju:2, 1     )   &
                                        +    u(1:nv,ip:iu:2,jp:ju:2, 1     ))  &
                                        +   (u(1:nv,il:iu:2,jp:ju:2, 1     )   &
                                        +    u(1:nv,ip:iu:2,jl:ju:2, 1     )))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(:,is:it,js:jt,ks:kt) =                                             &
                               1.25d-01 * (((u(1:nv,il:iu:2,jl:ju:2,kl:ku:2)   &
                                        +    u(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))  &
                                        +   (u(1:nv,il:iu:2,jl:ju:2,kp:ku:2)   &
                                        +    u(1:nv,ip:iu:2,jp:ju:2,kl:ku:2))) &
                                        +  ((u(1:nv,il:iu:2,jp:ju:2,kp:ku:2)   &
                                        +    u(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))  &
                                        +   (u(1:nv,il:iu:2,jp:ju:2,kl:ku:2)   &
                                        +    u(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_restrict
!
!===============================================================================
!
! boundary_prolong: subroutine copies the prolongated interior of a neighbor in
!                   order to update a boundary of the current block
!
!===============================================================================
!
  subroutine boundary_prolong(pdata, u, idir, iside, iface)

    use blocks        , only : block_data
    use coordinates   , only : ng, im, ih, ib, ie, ieu                         &
                             , nd, jm, jh, jb, je, jeu                         &
                             , nh, km, kh, kb, ke, keu
    use equations     , only : nv
    use interpolations, only : limiter

    implicit none

! arguments
!
    type(block_data), pointer  , intent(inout) :: pdata
    real   , dimension(:,:,:,:), intent(in)    :: u
    integer                    , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: i, j, k, q
    integer :: ic, jc, kc, ip, jp, kp
    integer :: il, jl, kl, iu, ju, ku
    integer :: is, js, ks, it, jt, kt
    real    :: dul, dur, dux, duy, duz
!
!-------------------------------------------------------------------------------
!
! prepare indices
!
    select case(idir)

    case(1)

! X indices
!
      if (iside .eq. 1) then
        is = 1
      else
        is = ieu
      end if

      il = 2
      iu = 1 + nh

! Y indices
!
      jc = mod(iface - 1, 2)
      js = 1
      jl = jb - nh + (jh - ng) * jc
      ju = jh + nh + (jh - ng) * jc

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2
      ks = 1
      kl = kb - nh + (kh - ng) * kc
      ku = kh + nh + (kh - ng) * kc
#endif /* NDIMS == 3 */

! update variable boundaries
!
#if NDIMS == 2
      do k = 1, km
        kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do k = kl, ku
        kt = 2 * (k - kl) + ks
        kp = kt + 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          jt = 2 * (j - jl) + js
          jp = jt + 1
          do i = il, iu
            it = 2 * (i - il) + is
            ip = it + 1

            pdata%u(:,it,jt,kt) = u(:,i,j,k)
            pdata%u(:,ip,jt,kt) = u(:,i,j,k)
            pdata%u(:,it,jp,kt) = u(:,i,j,k)
            pdata%u(:,ip,jp,kt) = u(:,i,j,k)
#if NDIMS == 3
            pdata%u(:,it,jt,kp) = u(:,i,j,k)
            pdata%u(:,ip,jt,kp) = u(:,i,j,k)
            pdata%u(:,it,jp,kp) = u(:,i,j,k)
            pdata%u(:,ip,jp,kp) = u(:,i,j,k)
#endif /* NDIMS == 3 */
          end do
        end do
      end do

! Y indices
!
      jc = mod(iface - 1, 2)
      js = jb
      jl = jb + (jh - ng) * jc
      ju = jh + (jh - ng) * jc

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2
      ks = kb
      kl = kb + (kh - ng) * kc
      ku = kh + (kh - ng) * kc
#endif /* NDIMS == 3 */

    case(2)

! X indices
!
      ic = mod(iface - 1, 2)
      is = 1
      il = ib - nh + (ih - ng) * ic
      iu = ih + nh + (ih - ng) * ic

! Y indices
!
      if (iside .eq. 1) then
        js = 1
      else
        js = jeu
      end if

      jl = 2
      ju = 1 + nh

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2
      ks = 1
      kl = kb - nh + (kh - ng) * kc
      ku = kh + nh + (kh - ng) * kc
#endif /* NDIMS == 3 */

! update variable boundaries
!
#if NDIMS == 2
      do k = 1, km
        kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do k = kl, ku
        kt = 2 * (k - kl) + ks
        kp = kt + 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          jt = 2 * (j - jl) + js
          jp = jt + 1
          do i = il, iu
            it = 2 * (i - il) + is
            ip = it + 1

            pdata%u(:,it,jt,kt) = u(:,i,j,k)
            pdata%u(:,ip,jt,kt) = u(:,i,j,k)
            pdata%u(:,it,jp,kt) = u(:,i,j,k)
            pdata%u(:,ip,jp,kt) = u(:,i,j,k)
#if NDIMS == 3
            pdata%u(:,it,jt,kp) = u(:,i,j,k)
            pdata%u(:,ip,jt,kp) = u(:,i,j,k)
            pdata%u(:,it,jp,kp) = u(:,i,j,k)
            pdata%u(:,ip,jp,kp) = u(:,i,j,k)
#endif /* NDIMS == 3 */
          end do
        end do
      end do

! X indices
!
      ic = mod(iface - 1, 2)
      is = ib
      il = ib + (ih - ng) * ic
      iu = ih + (ih - ng) * ic

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2
      ks = kb
      kl = kb + (kh - ng) * kc
      ku = kh + (kh - ng) * kc
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

! X indices
!
      ic = mod(iface - 1, 2)
      is = 1
      il = ib - nh + (ih - ng) * ic
      iu = ih + nh + (ih - ng) * ic

! Y indices
!
      jc =    (iface - 1) / 2
      js = 1
      jl = jb - nh + (jh - ng) * jc
      ju = jh + nh + (jh - ng) * jc

! Z indices
!
      if (iside .eq. 1) then
        ks = 1
      else
        ks = keu
      end if

      kl = 2
      ku = 1 + nh

! update variable boundaries
!
#if NDIMS == 2
      do k = 1, km
        kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do k = kl, ku
        kt = 2 * (k - kl) + ks
        kp = kt + 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          jt = 2 * (j - jl) + js
          jp = jt + 1
          do i = il, iu
            it = 2 * (i - il) + is
            ip = it + 1

            pdata%u(:,it,jt,kt) = u(:,i,j,k)
            pdata%u(:,ip,jt,kt) = u(:,i,j,k)
            pdata%u(:,it,jp,kt) = u(:,i,j,k)
            pdata%u(:,ip,jp,kt) = u(:,i,j,k)
#if NDIMS == 3
            pdata%u(:,it,jt,kp) = u(:,i,j,k)
            pdata%u(:,ip,jt,kp) = u(:,i,j,k)
            pdata%u(:,it,jp,kp) = u(:,i,j,k)
            pdata%u(:,ip,jp,kp) = u(:,i,j,k)
#endif /* NDIMS == 3 */
          end do
        end do
      end do

! X indices
!
      ic = mod(iface - 1, 2)
      is = ib
      il = ib + (ih - ng) * ic
      iu = ih + (ih - ng) * ic

! Y indices
!
      jc =    (iface - 1) / 2
      js = jb
      jl = jb + (jh - ng) * jc
      ju = jh + (jh - ng) * jc
#endif /* NDIMS == 3 */

    end select

! update variable boundaries with the linear interpolation
!
#if NDIMS == 2
    do k = 1, km
      kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = kl, ku
      kt = 2 * (k - kl) + ks
      kp = kt + 1
#endif /* NDIMS == 3 */
      do j = jl, ju
        jt = 2 * (j - jl) + js
        jp = jt + 1
        do i = il, iu
          it = 2 * (i - il) + is
          ip = it + 1

          do q = 1, nv

            dul = u(q,i  ,j,k) - u(q,i-1,j,k)
            dur = u(q,i+1,j,k) - u(q,i  ,j,k)
            dux = limiter(0.25d+00, dul, dur)

            dul = u(q,i,j  ,k) - u(q,i,j-1,k)
            dur = u(q,i,j+1,k) - u(q,i,j  ,k)
            duy = limiter(0.25d+00, dul, dur)

#if NDIMS == 3
            dul = u(q,i,j,k  ) - u(q,i,j,k-1)
            dur = u(q,i,j,k+1) - u(q,i,j,k  )
            duz = limiter(0.25d+00, dul, dur)
#endif /* NDIMS == 3 */

#if NDIMS == 2
            pdata%u(q,it,jt,kt) = u(q,i,j,k) - (dux + duy)
            pdata%u(q,ip,jt,kt) = u(q,i,j,k) + (dux - duy)
            pdata%u(q,it,jp,kt) = u(q,i,j,k) + (duy - dux)
            pdata%u(q,ip,jp,kt) = u(q,i,j,k) + (dux + duy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pdata%u(q,it,jt,kt) = u(q,i,j,k) - dux - duy - duz
            pdata%u(q,ip,jt,kt) = u(q,i,j,k) + dux - duy - duz
            pdata%u(q,it,jp,kt) = u(q,i,j,k) - dux + duy - duz
            pdata%u(q,ip,jp,kt) = u(q,i,j,k) + dux + duy - duz
            pdata%u(q,it,jt,kp) = u(q,i,j,k) - dux - duy + duz
            pdata%u(q,ip,jt,kp) = u(q,i,j,k) + dux - duy + duz
            pdata%u(q,it,jp,kp) = u(q,i,j,k) - dux + duy + duz
            pdata%u(q,ip,jp,kp) = u(q,i,j,k) + dux + duy + duz
#endif /* NDIMS == 3 */
          end do
        end do
      end do
    end do

!-------------------------------------------------------------------------------
!
  end subroutine boundary_prolong
!
!===============================================================================
!
! boundary_specific: subroutine applies specific boundary conditions to the
!                    current block
!
!===============================================================================
!
  subroutine boundary_specific(pdata, idir, iside)

    use blocks       , only : block_data
    use coordinates  , only : ng, im, jm, km, ib, ibl, ie, ieu, jb    &
                            , jbl, je, jeu, kb, kbl, ke, keu
    use equations    , only : idn, imx, imy, imz, ibx, iby, ibz, ibp
    use error        , only : print_warning
    use equations    , only : nv

    implicit none

! arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    integer                  , intent(in)    :: idir, iside

! local variables
!
    integer :: ii, i, j, k, it, jt, kt, is, js, ks
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 10 * idir + iside

! perform update according to the flag
!
    select case(ii)

! left side along the X direction
!
    case(11)

      select case(xlbndry)

      case("reflecting", "reflect")

        do i = 1, ng
          it = ib  - i
          is = ibl + i

          pdata%u(  :,it,:,:) =   pdata%u(  :,is,:,:)
          pdata%u(imx,it,:,:) = - pdata%u(imx,is,:,:)
        end do

      case default ! "open" as default boundary conditions

        do i = 1, ng
          pdata%u(  :,i,:,:) = pdata%u(:,ib,:,:)
        end do

      end select

! right side along the X direction
!
    case(12)

      select case(xubndry)

      case("reflecting", "reflect")

        do i = 1, ng
          it = ie  + i
          is = ieu - i

          pdata%u(  :,it,:,:) =   pdata%u(  :,is,:,:)
          pdata%u(imx,it,:,:) = - pdata%u(imx,is,:,:)
        end do

      case default ! "open" as default boundary conditions

        do i = ieu, im
          pdata%u(  :,i,:,:) = pdata%u(:,ie,:,:)
        end do

      end select

! left side along the Y direction
!
    case(21)

      select case(ylbndry)

      case("reflecting", "reflect")

        do j = 1, ng
          jt = jb  - j
          js = jbl + j

          pdata%u(  :,:,jt,:) =   pdata%u(  :,:,js,:)
          pdata%u(imy,:,jt,:) = - pdata%u(imy,:,js,:)
        end do

      case default ! "open" as default boundary conditions

        do j = 1, ng
          pdata%u(  :,:,j,:) = pdata%u(:,:,jb,:)
        end do

      end select

! right side along the Y direction
!
    case(22)

      select case(yubndry)

      case("reflecting", "reflect")

        do j = 1, ng
          jt = je  + j
          js = jeu - j

          pdata%u(  :,:,jt,:) =   pdata%u(  :,:,js,:)
          pdata%u(imy,:,jt,:) = - pdata%u(imy,:,js,:)
        end do

      case default ! "open" as default boundary conditions

        do j = jeu, jm
          pdata%u(  :,:,j,:) = pdata%u(:,:,je,:)
        end do

      end select

#if NDIMS == 3
! left side along the Z direction
!
    case(31)

      select case(zlbndry)

      case("reflecting", "reflect")

        do k = 1, ng
          kt = kb  - k
          ks = kbl + k

          pdata%u(  :,:,:,kt) =   pdata%u(  :,:,:,ks)
          pdata%u(imz,:,:,kt) = - pdata%u(imz,:,:,ks)
        end do

      case default ! "open" as default boundary conditions

        do k = 1, ng
          pdata%u(  :,:,:,k) = pdata%u(:,:,:,kb)
        end do

      end select

! right side along the Z direction
!
    case(32)

      select case(zubndry)

      case("reflecting", "reflect")

        do k = 1, ng
          kt = ke  + k
          ks = keu - k

          pdata%u(  :,:,:,kt) =   pdata%u(  :,:,:,ks)
          pdata%u(imz,:,:,kt) = - pdata%u(imz,:,:,ks)
        end do

      case default ! "open" as default boundary conditions

        do k = keu, km
          pdata%u(  :,:,:,k) = pdata%u(:,:,:,ke)
        end do

      end select
#endif /* NDIMS == 3 */

    case default
      call print_warning("boundaries::boundary_specific"                       &
                       , "Wrong direction or side of the boundary condition!")

    end select

!-------------------------------------------------------------------------------
!
  end subroutine boundary_specific

!===============================================================================
!
end module
