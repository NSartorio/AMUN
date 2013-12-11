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
!! module: BOUNDARIES
!!
!!  This module handles the boundary synchronization.
!!
!!******************************************************************************
!
module boundaries

! module variables are not implicit by default
!
  implicit none

! module parameters for the boundary update order and boundary type
!
  character(len = 32), save     :: xlbndry = "periodic"
  character(len = 32), save     :: xubndry = "periodic"
  character(len = 32), save     :: ylbndry = "periodic"
  character(len = 32), save     :: yubndry = "periodic"
  character(len = 32), save     :: zlbndry = "periodic"
  character(len = 32), save     :: zubndry = "periodic"

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_BOUNDARIES:
! --------------------------------
!
!   Subroutine initializes the module BOUNDARIES by setting its parameters.
!
!
!===============================================================================
!
  subroutine initialize_boundaries()

! include external procedures
!
    use parameters    , only : get_parameter_string

! include external variables
!
#ifdef MPI
    use mpitools      , only : pdims, pcoords, periodic
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get runtime values for the boundary types
!
    call get_parameter_string ("xlbndry" , xlbndry)
    call get_parameter_string ("xubndry" , xubndry)
    call get_parameter_string ("ylbndry" , ylbndry)
    call get_parameter_string ("yubndry" , yubndry)
    call get_parameter_string ("zlbndry" , zlbndry)
    call get_parameter_string ("zubndry" , zubndry)

#ifdef MPI
! change the internal boundaries to "exchange" type for the MPI update
!
    if (pdims(1) .gt. 1) then
      if (periodic(1)) then
        xlbndry       = "exchange"
        xubndry       = "exchange"
      else
        if (pcoords(1) .gt. 0         ) then
          xlbndry       = "exchange"
        end if
        if (pcoords(1) .lt. pdims(1)-1) then
          xubndry       = "exchange"
        end if
      end if
    end if

    if (pdims(2) .gt. 1) then
      if (periodic(2)) then
        ylbndry       = "exchange"
        yubndry       = "exchange"
      else
        if (pcoords(2) .gt. 0         ) then
          ylbndry       = "exchange"
        end if
        if (pcoords(2) .lt. pdims(2)-1) then
          yubndry       = "exchange"
        end if
      end if
    end if

    if (pdims(3) .gt. 1) then
      if (periodic(3)) then
        zlbndry       = "exchange"
        zubndry       = "exchange"
      else
        if (pcoords(3) .gt. 0         ) then
          zlbndry       = "exchange"
        end if
        if (pcoords(3) .lt. pdims(3)-1) then
          zubndry       = "exchange"
        end if
      end if
    end if
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_boundaries
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

! include external variables
!
    use coordinates   , only : toplev
    use mpitools      , only : periodic

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: idir, ilev
!
!-------------------------------------------------------------------------------
!
! 1. FIRST FILL OUT THE CORNERS WITH THE EXTRAPOLATED VARIABLES
!
    call update_corners()

! step down from the top level
!
    do ilev = toplev, 1, -1

! iterate over all directions
!
      do idir = 1, NDIMS

! update boundaries which don't have neighbors and which are not periodic
!
        if (.not. periodic(idir)) call specific_boundaries(ilev, idir)

! copy boundaries between blocks at the same levels
!
        call copy_boundaries(ilev, idir)

      end do ! directions

! restrict blocks from higher level neighbours
!
      do idir = 1, NDIMS

        call restrict_boundaries(ilev - 1, idir)

      end do

    end do ! levels

! step up from the first level
!
    do ilev = 1, toplev

! prolong boundaries from lower level neighbours
!
      do idir = 1, NDIMS

        call prolong_boundaries(ilev, idir)

      end do

      do idir = 1, NDIMS

! update boundaries which don't have neighbors and which are not periodic
!
        if (.not. periodic(idir)) call specific_boundaries(ilev, idir)

! copy boundaries between blocks at the same levels
!
        call copy_boundaries(ilev, idir)

      end do

    end do ! levels

!-------------------------------------------------------------------------------
!
  end subroutine boundary_variables
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
!
!===============================================================================
!
  subroutine specific_boundaries(ilev, idir)

! include external variables
!
    use blocks        , only : block_meta, list_meta
    use blocks        , only : nsides
#ifdef MPI
    use mpitools      , only : nproc
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
! assign the pointer to the first block on the list
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
        end if ! local process
#endif /* MPI */
      end if ! leaf

! assign the pointer to the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

!-------------------------------------------------------------------------------
!
  end subroutine specific_boundaries
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

!-------------------------------------------------------------------------------
!
  end subroutine prolong_boundaries
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

!-------------------------------------------------------------------------------
!
  end subroutine copy_boundaries
!
!===============================================================================
!
! boundary_fluxes: subroutine sweeps over all leaf blocks and if it
!                          finds that two neighbors lay at different levels, it
!                          corrects the numerical fluxes of block at lower level
!                          copying the flux from higher level neighbor
!
!===============================================================================
!
  subroutine boundary_fluxes()

    use blocks   , only : block_meta, block_data, list_meta
    use blocks   , only : nsides, nfaces
    use coordinates, only : toplev
    use coordinates, only : ibl, ie, jbl, je, kbl, ke
#ifdef MPI
    use blocks   , only : block_info, pointer_info
    use coordinates, only : im, jm, km
    use mpitools , only : send_real_array, receive_real_array
    use mpitools , only : nprocs, nproc
    use equations, only : nv
#endif /* MPI */

    implicit none

! local variables
!
    integer :: idir, iside, iface
    integer :: is, js, ks
#ifdef MPI
    integer :: irecv, isend, nblocks, itag, l, iret

! local arrays
!
    integer     , dimension(NDIMS,0:nprocs-1,0:nprocs-1) :: block_counter
    real(kind=8), dimension(:,:,:,:), allocatable      :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo

! local pointer arrays
!
    type(pointer_info), dimension(NDIMS,0:nprocs-1,0:nprocs-1) :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! do not correct fluxes if we do not use adaptive mesh
!
    if (toplev .eq. 1) return

#ifdef MPI
! reset the block counter
!
    block_counter(:,:,:) = 0

! nullify info pointers
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1
        do idir = 1, NDIMS
          nullify(block_array(idir,irecv,isend)%ptr)
        end do
      end do
    end do
#endif /* MPI */

! scan all meta blocks in the list
!
    pmeta => list_meta
    do while(associated(pmeta))

! check if the meta block is a leaf
!
      if (pmeta%leaf) then

! iterate over all neighbors
!
        do idir = 1, NDIMS
          do iside = 1, nsides
            do iface = 1, nfaces

! associate pneigh with the current neighbor
!
              pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
              if (associated(pneigh)) then

! check if the current block lays at lower level than its neighbor
!
                if (pmeta%level .lt. pneigh%level) then

#ifdef MPI
! check if the block and neighbor are on the local processors
!
                  if (pmeta%cpu .eq. nproc .and. pneigh%cpu .eq. nproc) then
#endif /* MPI */

! update directional flux from the neighbor
!
                    select case(idir)
                    case(1)

                      if (iside .eq. 1) then
                        is = ie
                      else
                        is = ibl
                      end if

                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,is,:,:), idir, iside, iface)

                    case(2)

                      if (iside .eq. 1) then
                        js = je
                      else
                        js = jbl
                      end if

                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,:,js,:), idir, iside, iface)

#if NDIMS == 3
                    case(3)

                      if (iside .eq. 1) then
                        ks = ke
                      else
                        ks = kbl
                      end if

                      call correct_flux(pmeta%data                             &
                           , pneigh%data%f(idir,:,:,:,ks), idir, iside, iface)
#endif /* NDIMS == 3 */
                    end select

#ifdef MPI
                  else

! increase the counter for number of blocks to exchange
!
                    block_counter(idir,pmeta%cpu,pneigh%cpu) =                 &
                                  block_counter(idir,pmeta%cpu,pneigh%cpu) + 1

! allocate new info object
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

! if the list is not emply append the created block
!
                    if (associated(block_array(idir,pmeta%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(idir,pmeta%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(idir,pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! pmeta and pneigh on local cpu
#endif /* MPI */

                end if ! pmeta level < pneigh level

              end if ! pneigh associated

            end do ! iface
          end do ! iside
        end do ! idir

      end if ! leaf

      pmeta => pmeta%next ! assign pointer to the next meta block in the list
    end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
    do irecv = 0, nprocs - 1
      do isend = 0, nprocs - 1
        do idir = 1, NDIMS

! process only pairs which have boundaries to exchange
!
          if (block_counter(idir,irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(idir,irecv,isend)

! prepare the tag for communication
!
            itag = (irecv * nprocs + isend) * nprocs + idir

! allocate space for variables
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
            if (isend .eq. nproc) then

! fill out the buffer with block data
!
              l = 1

              select case(idir)
              case(1)

                pinfo => block_array(idir,irecv,isend)%ptr
                do while(associated(pinfo))

                  if (pinfo%side .eq. 1) then
                    is = ie
                  else
                    is = ibl
                  end if

                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,is,:,:)

                  pinfo => pinfo%prev
                  l = l + 1
                end do

              case(2)

                pinfo => block_array(idir,irecv,isend)%ptr
                do while(associated(pinfo))

                  if (pinfo%side .eq. 1) then
                    js = je
                  else
                    js = jbl
                  end if

                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,:,js,:)

                  pinfo => pinfo%prev
                  l = l + 1
                end do

#if NDIMS == 3
              case(3)

                pinfo => block_array(idir,irecv,isend)%ptr
                do while(associated(pinfo))

                  if (pinfo%side .eq. 1) then
                    ks = ke
                  else
                    ks = kbl
                  end if

                  rbuf(l,:,:,:) = pinfo%neigh%data%f(idir,:,:,:,ks)

                  pinfo => pinfo%prev
                  l = l + 1
                end do
#endif /* NDIMS == 3 */
              end select

! send data buffer
!
              call send_real_array(size(rbuf(:,:,:,:)), irecv, itag, rbuf(:,:,:,:), iret)

            end if ! isend = nproc

! if irecv == nproc we are receiving data
!
            if (irecv .eq. nproc) then

! receive data
!
              call receive_real_array(size(rbuf(:,:,:,:)), isend, itag, rbuf(:,:,:,:), iret)

! iterate over all received blocks and update fluxes
!
              l = 1

              select case(idir)
              case(1)

                pinfo => block_array(idir,irecv,isend)%ptr

                do while(associated(pinfo))

! set indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! set pointers
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

                  pinfo => pinfo%prev

                  l = l + 1

                end do

              case(2)

                pinfo => block_array(idir,irecv,isend)%ptr

                do while(associated(pinfo))

! set indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! set pointers
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

                  pinfo => pinfo%prev

                  l = l + 1

                end do

#if NDIMS == 3
              case(3)

                pinfo => block_array(idir,irecv,isend)%ptr

                do while(associated(pinfo))

! set indices
!
                  iside = pinfo%side
                  iface = pinfo%face

! set pointers
!
                  pmeta  => pinfo%block
                  pneigh => pmeta%neigh(idir,iside,iface)%ptr

! update fluxes
!
                  call correct_flux(pmeta%data, rbuf(l,:,:,:)                  &
                                                         , idir, iside, iface)

                  pinfo => pinfo%prev

                  l = l + 1

                end do
#endif /* NDIMS == 3 */
              end select

            end if ! irecv = nproc

! deallocate buffers
!
            deallocate(rbuf)

! deallocate info blocks
!
            pinfo => block_array(idir,irecv,isend)%ptr
            do while(associated(pinfo))
              block_array(idir,irecv,isend)%ptr => pinfo%prev

              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

              deallocate(pinfo)

              pinfo => block_array(idir,irecv,isend)%ptr
            end do

          end if ! if block_count > 0
        end do ! idir
      end do ! isend
    end do ! irecv
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_fluxes
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
! boundary_restrict: subroutine copies the restricted interior of the neighbor
!                    in order to update a boundary of the current block
!
!===============================================================================
!
  subroutine boundary_restrict(pdata, u, idir, iside, iface)

    use blocks   , only : block_data
    use coordinates, only : ng, im, ih, ib, ie, ieu           &
                        , nd, jm, jh, jb, je, jeu           &
                        , nh, km, kh, kb, ke, keu
    use equations, only : nv

    implicit none

! arguments
!
    type(block_data), pointer  , intent(inout) :: pdata
    real   , dimension(:,:,:,:), intent(in)    :: u
    integer                    , intent(in)    :: idir, iside, iface

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
      if (iside .eq. 1) then
        is = 1
        it = ng
      else
        is = ieu
        it = im
      end if

      il = 1
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
      kc =    (iface - 1) / 2

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
      if (iside .eq. 1) then
        js = 1
        jt = ng
      else
        js = jeu
        jt = jm
      end if

      jl = 1
      ju = nd
      jp = jl + 1

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2

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
      jc =    (iface - 1) / 2

      js = jb - nh + (jh - nh) * jc
      jt = jh      + (jh - nh) * jc

      jl =  1 + ng * jc
      ju = je + ng * jc
      jp = jl + 1

! Z indices
!
      if (iside .eq. 1) then
        ks = 1
        kt = ng
      else
        ks = keu
        kt = km
      end if

      kl = 1
      ku = nd
      kp = kl + 1
#endif /* NDIMS == 3 */

    end select

! update variable boundaries
!
#if NDIMS == 2
    pdata%u(:,is:it,js:jt,:) = 0.25d0 * (u(:,il:iu:2,jl:ju:2,:)                &
                                       + u(:,ip:iu:2,jl:ju:2,:)                &
                                       + u(:,il:iu:2,jp:ju:2,:)                &
                                       + u(:,ip:iu:2,jp:ju:2,:))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(:,is:it,js:jt,ks:kt) = 0.125d0 * (u(:,il:iu:2,jl:ju:2,kl:ku:2)     &
                                            + u(:,ip:iu:2,jl:ju:2,kl:ku:2)     &
                                            + u(:,il:iu:2,jp:ju:2,kl:ku:2)     &
                                            + u(:,ip:iu:2,jp:ju:2,kl:ku:2)     &
                                            + u(:,il:iu:2,jl:ju:2,kp:ku:2)     &
                                            + u(:,ip:iu:2,jl:ju:2,kp:ku:2)     &
                                            + u(:,il:iu:2,jp:ju:2,kp:ku:2)     &
                                            + u(:,ip:iu:2,jp:ju:2,kp:ku:2))
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
    use interpolations, only : minmod3
    use equations     , only : nv

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
            dux = 0.25d0 * minmod3(dul, dur)

            dul = u(q,i,j  ,k) - u(q,i,j-1,k)
            dur = u(q,i,j+1,k) - u(q,i,j  ,k)
            duy = 0.25d0 * minmod3(dul, dur)

#if NDIMS == 3
            dul = u(q,i,j,k  ) - u(q,i,j,k-1)
            dur = u(q,i,j,k+1) - u(q,i,j,k  )
            duz = 0.25d0 * minmod3(dul, dur)
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
