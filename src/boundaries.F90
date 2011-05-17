!!******************************************************************************
!!
!! module: BOUNDARIES - subroutines for handling the boundary conditions
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!******************************************************************************
!!
!!  This file is part of the AMUN code.
!!
!!  This program is free software; you can redistribute it and/or
!!  modify it under the terms of the GNU General Public License
!!  as published by the Free Software Foundation; either version 2
!!  of the License, or (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program; if not, write to the Free Software Foundation,
!!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!!
!!******************************************************************************
!!
!
module boundaries

  implicit none

  contains
!
!===============================================================================
!
! boundary_variables: subroutine sweeps over all leaf blocks and performs
!                     their conserved variables boundary update
!
!===============================================================================
!
  subroutine boundary_variables()

    use blocks   , only : ndims, nsides, nfaces
    use blocks   , only : block_meta, block_data, block_info, pointer_info     &
                        , list_meta
    use config   , only : periodic, ng, nd, nh
    use config   , only : im, jm, km
    use config   , only : ib, ibu, iel, ie, jb, jbu, jel, je, kb, kbu, kel, ke
    use timer    , only : start_timer, stop_timer
#ifdef MPI
    use mpitools , only : ncpu, ncpus, is_master, msendf, mrecvf
    use variables, only : nqt
#endif /* MPI */

! declare variables
!
    implicit none

! local variables
!
    integer :: idir, iside, iface, nside, nface
    integer :: il, jl, kl, iu, ju, ku
#ifdef MPI
    integer :: isend, irecv, nblocks, itag, l

! local arrays
!
    integer     , dimension(3,0:ncpus-1,0:ncpus-1)    :: block_counter
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_data), pointer :: pdata
#ifdef MPI
    type(block_info), pointer :: pinfo
    type(pointer_info), dimension(3,0:ncpus-1,0:ncpus-1) :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! start the boundary timer
!
    call start_timer(4)

! update boundaries direction by direction
!
    do idir = 1, ndims

! first update boundaries which don't have neighbors and which are not periodic
!
      if (.not. periodic(idir)) then

! associate a pointer with the first meta block
!
        pmeta => list_meta

! iterate over all meta blocks
!
        do while(associated(pmeta))

! check if the current meta block is a leaf
!
          if (pmeta%leaf) then

#ifdef MPI
! check if the current block belongs to the current processor
!
            if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */

! iterate over all neighbors
!
              do iside = 1, nsides

! associate a pointer with the current neighbor
!
                pneigh => pmeta%neigh(idir,iside,1)%ptr

! check if the neighbor is not associated; it means that this is a specific
! boundary
!
                if (.not. associated(pneigh))                                  &
                               call boundary_specific(pmeta%data, idir, iside)

              end do ! sides
#ifdef MPI
            end if ! current cpu
#endif /* MPI */
          end if ! leaf

! associate the pointer with the next meta block
!
          pmeta => pmeta%next

        end do ! meta blocks

      end if ! not periodic

#ifdef MPI
! reset the exchange block counters
!
      block_counter(:,:,:) = 0

! nullify the info pointers
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1
          nullify(block_array(1,irecv,isend)%ptr)
          nullify(block_array(2,irecv,isend)%ptr)
          nullify(block_array(3,irecv,isend)%ptr)
        end do
      end do
#endif /* MPI */

!! now, update boundaries from the neighbors at the same or higher levels
!!
! associate a pointer with the first meta block
!
      pmeta => list_meta

! iterate over all meta blocks
!
      do while(associated(pmeta))

! check if the current meta block is a leaf
!
        if (pmeta%leaf) then

! iterate over all sides and faces along the current direction
!
          do iside = 1, nsides
            do iface = 1, nfaces

! associate a pointer with the current neighbor
!
              pneigh => pmeta%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
              if (associated(pneigh)) then

! if the neighbor is at the same level
!
                if (pmeta%level .eq. pneigh%level) then

! perform update only for the first face, since all faces point the same block
!
                  if (iface .eq. 1) then

#ifdef MPI
! check if the current meta block and its neighbor lay on the same processor
!
                    if (pmeta%cpu .eq. pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                      if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */

! assign a pointer to the data structure of the current block
!
                        pdata  => pmeta%data

! update the boundaries of the current block
!
                        select case(idir)
                        case(1)
                          if (iside .eq. 1) then
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,iel:ie,:,:), idir, iside)
                          else
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,ib:ibu,:,:), idir, iside)
                          end if
                        case(2)
                          if (iside .eq. 1) then
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,:,jel:je,:), idir, iside)
                          else
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,:,jb:jbu,:), idir, iside)
                          end if
#if NDIMS == 3
                        case(3)
                          if (iside .eq. 1) then
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,:,:,kel:ke), idir, iside)
                          else
                            call boundary_copy(pdata                           &
                                   , pneigh%data%u(:,:,:,kb:kbu), idir, iside)
                          end if
#endif /* NDIMS == 3 */
                        end select

#ifdef MPI
                      end if ! pmeta on the current cpu

                    else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                      block_counter(1,pmeta%cpu,pneigh%cpu) =                  &
                                     block_counter(1,pmeta%cpu,pneigh%cpu) + 1

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
                      if (associated(block_array(1,pmeta%cpu,pneigh%cpu)%ptr)) then
                        pinfo%prev => block_array(1,pmeta%cpu,pneigh%cpu)%ptr
                        nullify(pinfo%next)
                      end if

! point the list to the last created block
!
                      block_array(1,pmeta%cpu,pneigh%cpu)%ptr => pinfo

                    end if ! block and neighbor on different processors
#endif /* MPI */

                  end if ! iface = 1

                end if ! block and neighbor at the same level

! if the neighbor is at the higher level
!
                if (pmeta%level .lt. pneigh%level) then

#ifdef MPI
! check if the current meta block and its neighbor lay on the same processor
!
                  if (pmeta%cpu .eq. pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                    if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */

! prepare indices of the neighbor array
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
                    end if ! pmeta on the current cpu

                  else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                    block_counter(2,pmeta%cpu,pneigh%cpu) =                    &
                                     block_counter(2,pmeta%cpu,pneigh%cpu) + 1

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
                    if (associated(block_array(2,pmeta%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(2,pmeta%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(2,pmeta%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! block and neighbor on different processors
#endif /* MPI */

                end if ! block at lower level than neighbor

! if the neighbor is at lower level
!
                if (pmeta%level .gt. pneigh%level) then

! perform update only for the first face, since all faces point the same block
!
                  if (iface .eq. 1) then

#ifdef MPI
! check if the current meta block and its neighbor lay on the same processor
!
                    if (pmeta%cpu .eq. pneigh%cpu) then

! check if the current meta block lays on the current processors
!
                      if (pmeta%cpu .eq. ncpu) then
#endif /* MPI */

! find the face of the current block which the neighbor belongs to
!
                        nside = 3 - iside
                        nface = 1
                        do while(pmeta%id .ne.                                 &
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
                            il = ie - nh + 1
                            iu = ie
                          else
                            il = ib
                            iu = ib + nh - 1
                          end if
                        case(2)
                          if (iside .eq. 1) then
                            jl = je - nh + 1
                            ju = je
                          else
                            jl = jb
                            ju = jb + nh - 1
                          end if
                        case(3)
                          if (iside .eq. 1) then
                            kl = ke - nh + 1
                            ku = ke
                          else
                            kl = kb
                            ku = kb + nh - 1
                          end if
                        end select

! assign a pointer to the data structure of the current block
!
                        pdata  => pmeta%data

! update the boundaries of the current block
!
                        call boundary_prolong(pdata                            &
                                         , pneigh%data%u(:,il:iu,jl:ju,kl:ku)  &
                                                         , idir, iside, nface)

#ifdef MPI
                      end if ! pmeta on the current cpu

                    else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                      block_counter(3,pmeta%cpu,pneigh%cpu) =                  &
                                     block_counter(3,pmeta%cpu,pneigh%cpu) + 1

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
                      if (associated(block_array(3,pmeta%cpu,pneigh%cpu)%ptr)) then
                        pinfo%prev => block_array(3,pmeta%cpu,pneigh%cpu)%ptr
                        nullify(pinfo%next)
                      end if

! point the list to the last created block
!
                      block_array(3,pmeta%cpu,pneigh%cpu)%ptr => pinfo

                    end if ! block and neighbor on different processors
#endif /* MPI */

                  end if ! iface = 1

                end if ! neighbor at lower level

              end if ! if neighbor is associated
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
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1

!! process blocks at the same level
!!
! process only pairs which have boundaries to exchange
!
          if (block_counter(1,irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(1,irecv,isend)

! prepare the tag for communication
!
            itag = 10 * (irecv * ncpus + isend + 1) + 1

! allocate space for variables
!
            select case(idir)
            case(1)
              allocate(rbuf(nblocks,nqt,ng,jm,km))
            case(2)
              allocate(rbuf(nblocks,nqt,im,ng,km))
#if NDIMS == 3
            case(3)
              allocate(rbuf(nblocks,nqt,im,jm,ng))
#endif /* NDIMS == 3 */
            end select

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! iterate over exchange blocks along the current direction and fill out
! the buffer with the block data
!
              select case(idir)
              case(1)
                l = 1
                pinfo => block_array(1,irecv,isend)%ptr
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
                pinfo => block_array(1,irecv,isend)%ptr
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
                pinfo => block_array(1,irecv,isend)%ptr
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
              call msendf(size(rbuf), irecv, itag, rbuf(:,:,:,:,:))

            end if ! isend = ncpu

! if irecv == ncpu we are receiving data
!
            if (irecv .eq. ncpu) then

! receive data
!
              call mrecvf(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:))

! iterate over all received blocks and update boundaries
!
              l = 1
              pinfo => block_array(1,irecv,isend)%ptr
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

            end if ! irecv = ncpu

! deallocate buffers
!
            if (allocated(rbuf)) deallocate(rbuf)

! deallocate info blocks
!
            pinfo => block_array(1,irecv,isend)%ptr
            do while(associated(pinfo))
              block_array(1,irecv,isend)%ptr => pinfo%prev

              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

              deallocate(pinfo)

              pinfo => block_array(1,irecv,isend)%ptr
            end do

          end if ! if block_count > 0

!! process blocks with the neighbors at higher levels
!!
! process only pairs which have boundaries to exchange
!
          if (block_counter(2,irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(2,irecv,isend)

! prepare the tag for communication
!
            itag = 10 * (irecv * ncpus + isend + 1) + 2

! allocate space for variables
!
            select case(idir)
            case(1)
              allocate(rbuf(nblocks,nqt,nd,jm,km))
            case(2)
              allocate(rbuf(nblocks,nqt,im,nd,km))
            case(3)
              allocate(rbuf(nblocks,nqt,im,jm,nd))
            end select

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! fill out the buffer with block data
!
              l = 1

              pinfo => block_array(2,irecv,isend)%ptr
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
              call msendf(size(rbuf), irecv, itag, rbuf(:,:,:,:,:))

            end if

! if irecv == ncpu we are receiving data
!
            if (irecv .eq. ncpu) then

! receive data
!
              call mrecvf(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:))

! iterate over all received blocks and update boundaries
!
              l = 1

              pinfo => block_array(2,irecv,isend)%ptr
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
                call boundary_restrict(pdata, rbuf(l,:,:,:,:)                  &
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
            pinfo => block_array(2,irecv,isend)%ptr
            do while(associated(pinfo))
              block_array(2,irecv,isend)%ptr => pinfo%prev

              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

              deallocate(pinfo)

              pinfo => block_array(2,irecv,isend)%ptr
            end do

          end if ! if block_count > 0

!! process blocks with the neighbors at lower levels
!!
! process only pairs which have boundaries to exchange
!
          if (block_counter(3,irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(3,irecv,isend)

! prepare the tag for communication
!
            itag = 10 * (irecv * ncpus + isend + 1) + 3

! allocate space for variables
!
            select case(idir)
            case(1)
              allocate(rbuf(nblocks,nqt,nh,jm,km))
            case(2)
              allocate(rbuf(nblocks,nqt,im,nh,km))
            case(3)
              allocate(rbuf(nblocks,nqt,im,jm,nh))
            end select

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! fill out the buffer with block data
!
              l = 1

              pinfo => block_array(3,irecv,isend)%ptr
              do while(associated(pinfo))

! prepare indices of the neighbor array
!
                select case(idir)
                case(1)
                  if (pinfo%side .eq. 1) then
                    il = ie - nh + 1
                    iu = ie
                  else
                    il = ib
                    iu = ib + nh - 1
                  end if

                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,il:iu,:,:)

                case(2)
                  if (pinfo%side .eq. 1) then
                    jl = je - nh + 1
                    ju = je
                  else
                    jl = jb
                    ju = jb + nh - 1
                  end if

                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,jl:ju,:)

                case(3)
                  if (pinfo%side .eq. 1) then
                    kl = ke - nh + 1
                    ku = ke
                  else
                    kl = kb
                    ku = kb + nh - 1
                  end if

                  rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,kl:ku)
                end select

                pinfo => pinfo%prev
                l = l + 1
              end do

! send data buffer
!
              call msendf(size(rbuf), irecv, itag, rbuf(:,:,:,:,:))

            end if

! if irecv == ncpu we are receiving data
!
            if (irecv .eq. ncpu) then

! receive data
!
              call mrecvf(size(rbuf(:,:,:,:,:)), isend, itag, rbuf(:,:,:,:,:))

! iterate over all received blocks and update boundaries
!
              l = 1

              pinfo => block_array(3,irecv,isend)%ptr
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
                call boundary_prolong(pdata, rbuf(l,:,:,:,:)                   &
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
            pinfo => block_array(3,irecv,isend)%ptr
            do while(associated(pinfo))
              block_array(3,irecv,isend)%ptr => pinfo%prev

              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

              deallocate(pinfo)

              pinfo => block_array(3,irecv,isend)%ptr
            end do

          end if ! if block_count > 0

        end do ! isend
      end do ! irecv
#endif /* MPI */

    end do ! directions

! stop the boundary timer
!
    call stop_timer(4)

!-------------------------------------------------------------------------------
!
  end subroutine boundary_variables
#ifdef CONSERVATIVE
!
!===============================================================================
!
! boundary_correct_fluxes: subroutine sweeps over all leaf blocks and if it
!                          finds that two neighbors lay at different levels, it
!                          corrects the numerical fluxes of block at lower level
!                          copying the flux from higher level neighbor
!
!===============================================================================
!
  subroutine boundary_correct_fluxes()

    use blocks   , only : block_meta, block_data, list_meta
    use blocks   , only : nsides, nfaces
    use config   , only : maxlev
    use config   , only : ibl, ie, jbl, je, kbl, ke
    use timer    , only : start_timer, stop_timer
#ifdef MPI
    use blocks   , only : block_info, pointer_info
    use config   , only : im, jm, km
    use mpitools , only : ncpus, ncpu, msendf, mrecvf
    use variables, only : nqt
#endif /* MPI */

    implicit none

! local variables
!
    integer :: idir, iside, iface
    integer :: is, js, ks
#ifdef MPI
    integer :: irecv, isend, nblocks, itag, l

! local arrays
!
    integer     , dimension(NDIMS,0:ncpus-1,0:ncpus-1) :: block_counter
    real(kind=8), dimension(:,:,:,:), allocatable      :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo

! local pointer arrays
!
    type(pointer_info), dimension(NDIMS,0:ncpus-1,0:ncpus-1) :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! do not correct fluxes if we do not use adaptive mesh
!
    if (maxlev .eq. 1) return

! start the boundary timer
!
    call start_timer(4)

#ifdef MPI
! reset the block counter
!
    block_counter(:,:,:) = 0

! nullify info pointers
!
    do irecv = 0, ncpus - 1
      do isend = 0, ncpus - 1
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
                  if (pmeta%cpu .eq. ncpu .and. pneigh%cpu .eq. ncpu) then
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
    do irecv = 0, ncpus - 1
      do isend = 0, ncpus - 1
        do idir = 1, NDIMS

! process only pairs which have boundaries to exchange
!
          if (block_counter(idir,irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(idir,irecv,isend)

! prepare the tag for communication
!
            itag = (irecv * ncpus + isend) * ncpus + idir

! allocate space for variables
!
            select case(idir)
            case(1)
              allocate(rbuf(nblocks,nqt,jm,km))
            case(2)
              allocate(rbuf(nblocks,nqt,im,km))
#if NDIMS == 3
            case(3)
              allocate(rbuf(nblocks,nqt,im,jm))
#endif /* NDIMS == 3 */
            end select

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

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
              call msendf(size(rbuf(:,:,:,:)), irecv, itag, rbuf(:,:,:,:))

            end if ! isend = ncpu

! if irecv == ncpu we are receiving data
!
            if (irecv .eq. ncpu) then

! receive data
!
              call mrecvf(size(rbuf(:,:,:,:)), isend, itag, rbuf(:,:,:,:))

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

            end if ! irecv = ncpu

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

! stop the boundary timer
!
    call stop_timer(4)

!-------------------------------------------------------------------------------
!
  end subroutine boundary_correct_fluxes
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
    use config   , only : ng, in, jn, kn, ib, ie, ibl, jb, je, jbl, kb, ke, kbl

    implicit none

! local variables
!
    integer :: i, ip, it, ih, il, iu, i1, i2
    integer :: j, jp, jt, jh, jl, ju, j1, j2
#if NDIMS == 3
    integer :: k, kp, kt, kh, kl, ku, k1, k2
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
      else ! right side
        it = ie
      end if

! convert face number to index
!
      jp = mod(iface - 1, 2)
#if NDIMS == 3
      kp = (iface - 1) / 2
#endif /* NDIMS == 3 */

! bounds for the perpendicular update
!
      jh = jn / 2
      jl = jb + jh * jp
      ju = jl + jh - 1
#if NDIMS == 3
      kh = kn / 2
      kl = kb + kh * kp
      ku = kl + kh - 1
#endif /* NDIMS == 3 */

! iterate over perpendicular direction
!
#if NDIMS == 2
      do j = jl, ju
        j1 = 2 * (j - jl) + jb
        j2 = j1 + 1

        pdata%f(idir,:,it,j,:) = 0.5d0 * (f(:,j1,:) + f(:,j2,:))
      end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do j = jl, ju
        j1 = 2 * (j - jl) + jb
        j2 = j1 + 1
        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,it,j,k) = 0.25d0 * (f(:,j1,k1) + f(:,j2,k1)           &
                                           + f(:,j1,k2) + f(:,j2,k2))
        end do
      end do
#endif /* NDIMS == 3 */

! Y direction
!
    case(2)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        jt = jbl
      else ! right side
        jt = je
      end if

! convert face number to index
!
      ip = mod(iface - 1, 2)
#if NDIMS == 3
      kp = (iface - 1) / 2
#endif /* NDIMS == 3 */

! bounds for the perpendicular update
!
      ih = in / 2
      il = ib + ih * ip
      iu = il + ih - 1
#if NDIMS == 3
      kh = kn / 2
      kl = kb + kh * kp
      ku = kl + kh - 1
#endif /* NDIMS == 3 */

! iterate over perpendicular direction
!
#if NDIMS == 2
      do i = il, iu
        i1 = 2 * (i - il) + ib
        i2 = i1 + 1

        pdata%f(idir,:,i,jt,:) = 0.5d0 * (f(:,i1,:) + f(:,i2,:))
      end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do i = il, iu
        i1 = 2 * (i - il) + ib
        i2 = i1 + 1

        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,i,jt,k) = 0.25d0 * (f(:,i1,k1) + f(:,i2,k1)           &
                                           + f(:,i1,k2) + f(:,i2,k2))
        end do
      end do
#endif /* NDIMS == 3 */

#if NDIMS == 3
! Z direction
!
    case(3)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        kt = kbl
      else ! right side
        kt = ke
      end if

! convert face number to index
!
      ip = mod(iface - 1, 2)
      jp = (iface - 1) / 2

! bounds for the perpendicular update
!
      ih = in / 2
      il = ib + ih * ip
      iu = il + ih - 1
      jh = jn / 2
      jl = jb + jh * jp
      ju = jl + jh - 1

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
#endif /* CONSERVTIVE */
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
    use config   , only : ng, im, jm, km, ibl, ieu, jbl, jeu, kbl, keu
    use variables, only : nqt

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
        pdata%u(1:nqt,  1:ibl,1:jm,1:km) = u(1:nqt,1:ng,1:jm,1:km)
      else
        pdata%u(1:nqt,ieu:im ,1:jm,1:km) = u(1:nqt,1:ng,1:jm,1:km)
      end if

    case(2)

      if (iside .eq. 1) then
        pdata%u(1:nqt,1:im,  1:jbl,1:km) = u(1:nqt,1:im,1:ng,1:km)
      else
        pdata%u(1:nqt,1:im,jeu:jm ,1:km) = u(1:nqt,1:im,1:ng,1:km)
      end if

#if NDIMS == 3
    case(3)

      if (iside .eq. 1) then
        pdata%u(1:nqt,1:im,1:jm,  1:kbl) = u(1:nqt,1:im,1:jm,1:ng)
      else
        pdata%u(1:nqt,1:im,1:jm,keu:km ) = u(1:nqt,1:im,1:jm,1:ng)
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
    use config   , only : ng, im, ih, ib, ie, ieu           &
                        , nd, jm, jh, jb, je, jeu           &
                        , nh, km, kh, kb, ke, keu
    use variables, only : nqt, nfl

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

    use blocks   , only : block_data
    use config   , only : ng, im, ih, ib, ie, ieu           &
                        , nd, jm, jh, jb, je, jeu           &
                        , nh, km, kh, kb, ke, keu

    implicit none

! arguments
!
    type(block_data), pointer  , intent(inout) :: pdata
    real   , dimension(:,:,:,:), intent(in)    :: u
    integer                    , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: i, j, k
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
      else
        is = ieu
      end if

      il = 1
      iu = nh

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

      jl = 1
      ju = nh

#if NDIMS == 3
! Z indices
!
      kc =    (iface - 1) / 2
      ks = 1
      kl = kb - nh + (kh - ng) * kc
      ku = kh + nh + (kh - ng) * kc
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

      kl = 1
      ku = nh
#endif /* NDIMS == 3 */

    end select

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
    use config       , only : xlbndry, xubndry, ylbndry, yubndry, zlbndry      &
                            , zubndry, ng, im, jm, km, ib, ibl, ie, ieu, jb    &
                            , jbl, je, jeu, kb, kbl, ke, keu
    use error        , only : print_warning
    use interpolation, only : limiter
    use variables    , only : nvr, nfl, imx, imy, imz
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    integer                  , intent(in)    :: idir, iside

! local variables
!
    integer :: ii, i, j, k, it, jt, kt, is, js, ks
#ifdef MHD
    integer :: im1, ip1
    integer :: jm1, jp1
#if NDIMS == 3
    integer :: km1, kp1
#endif /* NDIMS == 3 */
    real    :: dbm, dbp, dbx, dby, dbz
#endif /* MHD */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do i = 1, ng
          pdata%u(1:nfl,i,:,:) = pdata%u(1:nfl,ib,:,:)
          pdata%u(  iby,i,:,:) = pdata%u(  iby,ib,:,:)
          pdata%u(  ibz,i,:,:) = pdata%u(  ibz,ib,:,:)
#ifdef GLM
          pdata%u(  iph,i,:,:) = pdata%u(  iph,ib,:,:)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do j = 1, jm
            jm1 = max( 1, j-1)
            jp1 = min(jm, j+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(iby,ib,jp1,k) - pdata%u(iby,ib,j  ,k)
            dbm = pdata%u(iby,ib,j  ,k) - pdata%u(iby,ib,jm1,k)
            dby = limiter(dbp, dbm)
#if NDIMS == 3
            dbp = pdata%u(ibz,ib,j,kp1) - pdata%u(ibz,ib,j,k  )
            dbm = pdata%u(ibz,ib,j,k  ) - pdata%u(ibz,ib,j,km1)
            dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

! update normal component of magnetic field from the divergence free condition
!
            do i = ng, 1, -1
#if NDIMS == 2
              pdata%u(ibx,i,j,k) = pdata%u(ibx,i+1,j,k) + dby
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pdata%u(ibx,i,j,k) = pdata%u(ibx,i+1,j,k) + dby + dbz
#endif /* NDIMS == 3 */
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do i = 1, ng
          pdata%u(  :,i,:,:) = pdata%u(:,ib,:,:)
#if defined MHD && defined GLM
          pdata%u(iph,i,:,:) = 0.0d0
#endif /* MHD & GLM */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do i = ieu, im
          pdata%u(1:nfl,i,:,:) = pdata%u(1:nfl,ie,:,:)
          pdata%u(  iby,i,:,:) = pdata%u(  iby,ie,:,:)
          pdata%u(  ibz,i,:,:) = pdata%u(  ibz,ie,:,:)
#ifdef GLM
          pdata%u(  iph,i,:,:) = pdata%u(  iph,ie,:,:)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do j = 1, jm
            jm1 = max( 1, j-1)
            jp1 = min(jm, j+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(iby,ie,jp1,k) - pdata%u(iby,ie,j  ,k)
            dbm = pdata%u(iby,ie,j  ,k) - pdata%u(iby,ie,jm1,k)
            dby = limiter(dbp, dbm)
#if NDIMS == 3
            dbp = pdata%u(ibz,ie,j,kp1) - pdata%u(ibz,ie,j,k  )
            dbm = pdata%u(ibz,ie,j,k  ) - pdata%u(ibz,ie,j,km1)
            dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

! update normal component of magnetic field from the divergence free condition
!
            do i = ieu, im
#if NDIMS == 2
              pdata%u(ibx,i,j,k) = pdata%u(ibx,i-1,j,k) - dby
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pdata%u(ibx,i,j,k) = pdata%u(ibx,i-1,j,k) - dby - dbz
#endif /* NDIMS == 3 */
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do i = ieu, im
          pdata%u(  :,i,:,:) = pdata%u(:,ie,:,:)
#if defined MHD && defined GLM
          pdata%u(iph,i,:,:) = 0.0d0
#endif /* MHD & GLM */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do j = 1, ng
          pdata%u(1:nfl,:,j,:) = pdata%u(1:nfl,:,jb,:)
          pdata%u(  ibx,:,j,:) = pdata%u(  ibx,:,jb,:)
          pdata%u(  ibz,:,j,:) = pdata%u(  ibz,:,jb,:)
#ifdef GLM
          pdata%u(  iph,:,j,:) = pdata%u(  iph,:,jb,:)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(ibx,ip1,jb,k) - pdata%u(ibx,i  ,jb,k)
            dbm = pdata%u(ibx,i  ,jb,k) - pdata%u(ibx,im1,jb,k)
            dbx = limiter(dbp, dbm)
#if NDIMS == 3
            dbp = pdata%u(ibz,i,jb,kp1) - pdata%u(ibz,i,jb,k  )
            dbm = pdata%u(ibz,i,jb,k  ) - pdata%u(ibz,i,jb,km1)
            dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

! update normal component of magnetic field from the divergence free condition
!
            do j = ng, 1, -1
#if NDIMS == 2
              pdata%u(iby,i,j,k) = pdata%u(iby,i,j+1,k) + dbx
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pdata%u(iby,i,j,k) = pdata%u(iby,i,j+1,k) + dbx + dbz
#endif /* NDIMS == 3 */
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do j = 1, ng
          pdata%u(  :,:,j,:) = pdata%u(:,:,jb,:)
#if defined MHD && defined GLM
          pdata%u(iph,:,j,:) = 0.0d0
#endif /* MHD & GLM */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do j = jeu, jm
          pdata%u(1:nfl,:,j,:) = pdata%u(1:nfl,:,je,:)
          pdata%u(  ibx,:,j,:) = pdata%u(  ibx,:,je,:)
          pdata%u(  ibz,:,j,:) = pdata%u(  ibz,:,je,:)
#ifdef GLM
          pdata%u(  iph,:,j,:) = pdata%u(  iph,:,je,:)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(ibx,ip1,je,k) - pdata%u(ibx,i  ,je,k)
            dbm = pdata%u(ibx,i  ,je,k) - pdata%u(ibx,im1,je,k)
            dbx = limiter(dbp, dbm)
#if NDIMS == 3
            dbp = pdata%u(ibz,i,je,kp1) - pdata%u(ibz,i,je,k  )
            dbm = pdata%u(ibz,i,je,k  ) - pdata%u(ibz,i,je,km1)
            dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

! update normal component of magnetic field from the divergence free condition
!
            do j = jeu, jm
#if NDIMS == 2
              pdata%u(iby,i,j,k) = pdata%u(iby,i,j-1,k) - dbx
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pdata%u(iby,i,j,k) = pdata%u(iby,i,j-1,k) - dbx - dbz
#endif /* NDIMS == 3 */
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do j = jeu, jm
          pdata%u(  :,:,j,:) = pdata%u(:,:,je,:)
#if defined MHD && defined GLM
          pdata%u(iph,:,j,:) = 0.0d0
#endif /* MHD & GLM */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do k = 1, ng
          pdata%u(1:nfl,:,:,k) = pdata%u(1:nfl,:,:,kb)
          pdata%u(  ibx,:,:,k) = pdata%u(  ibx,:,:,kb)
          pdata%u(  iby,:,:,k) = pdata%u(  iby,:,:,kb)
#ifdef GLM
          pdata%u(  iph,:,:,k) = pdata%u(  iph,:,:,kb)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do j = 1, jm
          jm1 = max( 1, j-1)
          jp1 = min(jm, j+1)
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(ibx,ip1,j,kb) - pdata%u(ibx,i  ,j,kb)
            dbm = pdata%u(ibx,i  ,j,kb) - pdata%u(ibx,im1,j,kb)
            dbx = limiter(dbp, dbm)

            dbp = pdata%u(iby,i,jp1,kb) - pdata%u(iby,i,j  ,kb)
            dbm = pdata%u(iby,i,j  ,kb) - pdata%u(iby,i,jm1,kb)
            dbz = limiter(dbp, dbm)

! update normal component of magnetic field from the divergence free condition
!
            do k = ng, 1, -1
              pdata%u(ibz,i,j,k) = pdata%u(ibz,i,j,k+1) + dbx + dbz
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do k = 1, ng
          pdata%u(  :,:,:,k) = pdata%u(:,:,:,kb)
#if defined MHD && defined GLM
          pdata%u(iph,:,:,k) = 0.0d0
#endif /* MHD & GLM */
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

#ifdef MHD
      case("open", "divb")

! update fluid variables and transverse components of magnetic field
!
        do k = keu, km
          pdata%u(1:nfl,:,:,k) = pdata%u(1:nfl,:,:,ke)
          pdata%u(  ibx,:,:,k) = pdata%u(  ibx,:,:,ke)
          pdata%u(  iby,:,:,k) = pdata%u(  iby,:,:,ke)
#ifdef GLM
          pdata%u(  iph,:,:,k) = pdata%u(  iph,:,:,ke)
#endif /* GLM */
        end do

! update normal component of magnetic field
!
        do j = 1, jm
          jm1 = max( 1, j-1)
          jp1 = min(jm, j+1)
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)

! calculate derivatives of transverse components of magnetic field
!
            dbp = pdata%u(ibx,ip1,j,ke) - pdata%u(ibx,i  ,j,ke)
            dbm = pdata%u(ibx,i  ,j,ke) - pdata%u(ibx,im1,j,ke)
            dbx = limiter(dbp, dbm)

            dbp = pdata%u(iby,i,jp1,ke) - pdata%u(iby,i,j  ,ke)
            dbm = pdata%u(iby,i,j  ,ke) - pdata%u(iby,i,jm1,ke)
            dbz = limiter(dbp, dbm)

! update normal component of magnetic field from the divergence free condition
!
            do k = keu, km
              pdata%u(ibz,i,j,k) = pdata%u(ibz,i,j,k-1) - dbx - dbz
            end do
          end do
        end do
#endif /* MHD */

      case default ! "open" as default boundary conditions

        do k = keu, km
          pdata%u(  :,:,:,k) = pdata%u(:,:,:,ke)
#if defined MHD && defined GLM
          pdata%u(iph,:,:,k) = 0.0d0
#endif /* MHD & GLM */
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
