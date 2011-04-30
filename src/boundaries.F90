!!******************************************************************************
!!
!! module: boundaries - routines for handling the boundary conditions
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of AMUN.
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
!                    their boundary update
!
!===============================================================================
!
  subroutine boundary_variables()

! references to other modules
!
    use blocks   , only : block_meta, block_data, block_info, pointer_info     &
                        , list_meta
    use blocks   , only : ndims, nsides, nfaces
    use config   , only : periodic
#ifdef MPI
    use config   , only : im, jm, km
    use mpitools , only : ncpu, ncpus, is_master, msendf, mrecvf
    use variables, only : nqt
#endif /* MPI */

! declare variables
!
    implicit none

! local variables
!
    integer :: idir, iside, iface, nside, nface
#ifdef MPI
    integer :: isend, irecv, nblocks, itag, l

! local arrays
!
    integer     , dimension(0:ncpus-1,0:ncpus-1)    :: block_counter
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pblock, pneigh
    type(block_data), pointer :: pbdata, pndata
#ifdef MPI
    type(block_info), pointer :: pinfo
    type(pointer_info), dimension(0:ncpus-1,0:ncpus-1) :: block_array
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! update boundaries direction by direction
!
    do idir = 1, ndims

! first update boundaries which don't have neighbors and which are not periodic
!
      if (.not. periodic(idir)) then

        pblock => list_meta
        do while(associated(pblock))

#ifdef MPI
          if (pblock%leaf .and. pblock%cpu .eq. ncpu) then
#else /* MPI */
          if (pblock%leaf) then
#endif /* MPI */
            do iside = 1, nsides
              do iface = 1, nfaces
                pneigh => pblock%neigh(idir,iside,iface)%ptr
                if (.not. associated(pneigh) .and. iface .eq. 1)               &
                  call bnd_spec(pblock%data, idir, iside)
              end do ! faces
            end do ! sides
          end if ! leaf

          pblock => pblock%next ! assign pointer to the next block
        end do ! meta blocks

      end if

#ifdef MPI
! reset the block counter
!
      block_counter(:,:) = 0

! nullify info pointers
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1
          nullify(block_array(irecv,isend)%ptr)
        end do
      end do

#endif /* MPI */
! now, update boundaries from neighbors
!
      pblock => list_meta
      do while(associated(pblock))

        if (pblock%leaf) then

! iterate over all sides and faces
!
          do iside = 1, nsides
            do iface = 1, nfaces

! assign pointer to the neighbor
!
              pneigh => pblock%neigh(idir,iside,iface)%ptr

! process only associated neighbors
!
              if (associated(pneigh)) then

#ifdef MPI
! process the block and its neighbor belong to the same processor
!
                if (pblock%cpu .eq. pneigh%cpu) then

! process only blocks belonding to the current processor
!
                  if (pblock%cpu .eq. ncpu) then

! assign pointers to data structures of the current block and neighbor
!
                    pbdata => pblock%data
                    pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                    if (pblock%level .lt. pneigh%level)                        &
                      call bnd_rest(pblock%data, pneigh%data%u, idir, iside    &
                                                                    , iface)

                    if (pblock%level .eq. pneigh%level .and. iface .eq. 1)     &
                      call bnd_copy(pbdata, pndata%u, idir, iside)

                  end if ! if neighbors on the current cpu
                else ! if block and neighbor are on the same cpu

! the neighbor and current block are on different processors, so we need to prepare
! information about the all blocks belonging to different processors
!
                  if (pblock%level .le. pneigh%level) then

! increase the counter for number of blocks to exchange
!
                    block_counter(pblock%cpu,pneigh%cpu) =                     &
                                     block_counter(pblock%cpu,pneigh%cpu) + 1

! allocate new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pblock
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pblock%level - pneigh%level

! nullify pointers
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! if the list is not emply append the created block
!
                    if (associated(block_array(pblock%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(pblock%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(pblock%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! process only copying and restriction
                end if ! if block and neighbor are on the same cpu
#else /* MPI */
! assign pointers to data structures of the current block and neighbor
!
                pbdata => pblock%data
                pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                if (pblock%level .lt. pneigh%level)                            &
                  call bnd_rest(pblock%data, pneigh%data%u, idir, iside, iface)

                if (pblock%level .eq. pneigh%level .and. iface .eq. 1)         &
                  call bnd_copy(pbdata, pndata%u, idir, iside)
#endif /* MPI */

              end if ! if neighbor is associated
            end do ! faces
          end do ! sides
        end if ! leaf

! assign pointer to the next block
!
        pblock => pblock%next

      end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1

! process only pairs which have boundaries to exchange
!
          if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
            itag = irecv * ncpus + isend + ncpus + 1

! allocate space for variables
!
            allocate(rbuf(nblocks,nqt,im,jm,km))

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! fill out the buffer with block data
!
              l = 1

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,:)

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

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

! set indices
!
                iside = pinfo%side
                iface = pinfo%face

! update boundaries
!
                if (pinfo%level_difference .eq. -1)                            &
                  call bnd_rest(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside &
                                                                       , iface)
                if (pinfo%level_difference .eq. 0 .and. iface .eq. 1)          &
                  call bnd_copy(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside)

                pinfo => pinfo%prev
                l = l + 1
              end do

            end if

! deallocate buffers
!
            deallocate(rbuf)

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

#ifdef MPI
! reset the block counter
!
      block_counter(:,:) = 0

! nullify info pointers
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1
          nullify(block_array(irecv,isend)%ptr)
        end do
      end do

#endif /* MPI */
! perform the boundary prolongation at the end, since its interpolation
! requires values from the boundaries
!
      pblock => list_meta
      do while(associated(pblock))

        if (pblock%leaf) then

! iterate over all sides and faces
!
          do iside = 1, nsides
            do iface = 1, nfaces

! assign pointer to the neighbor
!
              pneigh => pblock%neigh(idir,iside,iface)%ptr

! process only associated neighbors
!
              if (associated(pneigh)) then

#ifdef MPI
! process the block and its neighbor belong to the same processor
!
                if (pblock%cpu .eq. pneigh%cpu) then

! process only blocks belonding to the current processor
!
                  if (pblock%cpu .eq. ncpu) then

! assign pointers to data structures of the current block and neighbor
!
                    pbdata => pblock%data
                    pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                    if (pblock%level .gt. pneigh%level .and. iface .eq. 1) then
                      nside = 3 - iside
                      nface = 1
                      do while(pblock%id .ne.                                  &
                                          pneigh%neigh(idir,nside,nface)%ptr%id)
                        nface = nface + 1
                      end do

                      pbdata => pblock%data
                      pndata => pneigh%data
                      call bnd_prol(pbdata, pndata%u, idir, iside, nface)
                    end if

                  end if ! if neighbors on the current cpu
                else ! if block and neighbor are on the same cpu

! the neighbor and current block are on different processors, so we need to
! prepare information about the all blocks belonging to different processors
!
                  if (pblock%level .gt. pneigh%level) then
! increase the counter for number of blocks to exchange
!
                    block_counter(pblock%cpu,pneigh%cpu) =                     &
                                     block_counter(pblock%cpu,pneigh%cpu) + 1

! allocate new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pblock
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pblock%level - pneigh%level

! nullify pointers
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! if the list is not emply append the created block
!
                    if (associated(block_array(pblock%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(pblock%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(pblock%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! process only prolongation
                end if ! if block and neighbor are on the same cpu
#else /* MPI */
! assign pointers to data structures of the current block and neighbor
!
                pbdata => pblock%data
                pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                if (pblock%level .gt. pneigh%level .and. iface .eq. 1) then
                  nside = 3 - iside
                  nface = 1
                  do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                    nface = nface + 1
                  end do

                  pbdata => pblock%data
                  pndata => pneigh%data
                  call bnd_prol(pbdata, pndata%u, idir, iside, nface)
                end if
#endif /* MPI */

              end if ! if neighbor is associated
            end do ! faces
          end do ! sides
        end if ! leaf

! assign pointer to the next block
!
        pblock => pblock%next

      end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1

! process only pairs which have boundaries to exchange
!
          if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
            itag = irecv * ncpus + isend + ncpus + 1

! allocate space for variables
!
            allocate(rbuf(nblocks,nqt,im,jm,km))

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! fill out the buffer with block data
!
              l = 1

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,:)

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

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

! set indices
!
                iside = pinfo%side
                iface = pinfo%face

! update boundaries
!
                if (pinfo%level_difference .eq. 1 .and. iface .eq. 1) then
                  pblock => pinfo%block
                  pneigh => pblock%neigh(idir,iside,iface)%ptr

                  nside = 3 - iside
                  nface = 1
                  do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                    nface = nface + 1
                  end do

                  call bnd_prol(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside &
                                                                       , nface)
                end if

                pinfo => pinfo%prev
                l = l + 1
              end do

            end if

! deallocate buffers
!
            deallocate(rbuf)

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
    end do ! directions

! repeat the boundary prolongation once again since the interpolation requires
! also perpendicular directions
!
    do idir = 1, ndims

#ifdef MPI
! reset the block counter
!
      block_counter(:,:) = 0

! nullify info pointers
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1
          nullify(block_array(irecv,isend)%ptr)
        end do
      end do

#endif /* MPI */
      pblock => list_meta
      do while(associated(pblock))

        if (pblock%leaf) then

! iterate over all sides and faces
!
          do iside = 1, nsides
            do iface = 1, nfaces

! assign pointer to the neighbor
!
              pneigh => pblock%neigh(idir,iside,iface)%ptr

! process only associated neighbors
!
              if (associated(pneigh)) then

#ifdef MPI
! process the block and its neighbor belong to the same processor
!
                if (pblock%cpu .eq. pneigh%cpu) then

! process only blocks belonding to the current processor
!
                  if (pblock%cpu .eq. ncpu) then

! assign pointers to data structures of the current block and neighbor
!
                    pbdata => pblock%data
                    pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                    if (pblock%level .gt. pneigh%level .and. iface .eq. 1) then
                      nside = 3 - iside
                      nface = 1
                      do while(pblock%id .ne.                                  &
                                          pneigh%neigh(idir,nside,nface)%ptr%id)
                        nface = nface + 1
                      end do

                      pbdata => pblock%data
                      pndata => pneigh%data
                      call bnd_prol(pbdata, pndata%u, idir, iside, nface)
                    end if

                  end if ! if neighbors on the current cpu
                else ! if block and neighbor are on the same cpu

! the neighbor and current block are on different processors, so we need to
! prepare information about the all blocks belonging to different processors
!
                  if (pblock%level .gt. pneigh%level) then

! increase the counter for number of blocks to exchange
!
                    block_counter(pblock%cpu,pneigh%cpu) =                     &
                                     block_counter(pblock%cpu,pneigh%cpu) + 1

! allocate new info object
!
                    allocate(pinfo)

! fill out its fields
!
                    pinfo%block            => pblock
                    pinfo%neigh            => pneigh
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
                    pinfo%level_difference =  pblock%level - pneigh%level

! nullify pointers
!
                    nullify(pinfo%prev)
                    nullify(pinfo%next)

! if the list is not emply append the created block
!
                    if (associated(block_array(pblock%cpu,pneigh%cpu)%ptr)) then
                      pinfo%prev => block_array(pblock%cpu,pneigh%cpu)%ptr
                      nullify(pinfo%next)
                    end if

! point the list to the last created block
!
                    block_array(pblock%cpu,pneigh%cpu)%ptr => pinfo

                  end if ! process only prolongation
                end if ! if block and neighbor are on the same cpu
#else /* MPI */
! assign pointers to data structures of the current block and neighbor
!
                pbdata => pblock%data
                pndata => pneigh%data

! depending on the level difference perform the proper boundary update
!
                if (pblock%level .gt. pneigh%level .and. iface .eq. 1) then
                  nside = 3 - iside
                  nface = 1
                  do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                    nface = nface + 1
                  end do

                  pbdata => pblock%data
                  pndata => pneigh%data
                  call bnd_prol(pbdata, pndata%u, idir, iside, nface)
                end if
#endif /* MPI */

              end if ! if neighbor is associated
            end do ! faces
          end do ! sides
        end if ! leaf

! assign pointer to the next block
!
        pblock => pblock%next

      end do ! meta blocks

#ifdef MPI
! iterate over sending and receiving processors
!
      do irecv = 0, ncpus - 1
        do isend = 0, ncpus - 1

! process only pairs which have boundaries to exchange
!
          if (block_counter(irecv,isend) .gt. 0) then

! obtain the number of blocks to exchange
!
            nblocks = block_counter(irecv,isend)

! prepare the tag for communication
!
            itag = irecv * ncpus + isend + ncpus + 1

! allocate space for variables
!
            allocate(rbuf(nblocks,nqt,im,jm,km))

! if isend == ncpu we are sending data
!
            if (isend .eq. ncpu) then

! fill out the buffer with block data
!
              l = 1

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

                rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,:)

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

              pinfo => block_array(irecv,isend)%ptr
              do while(associated(pinfo))

! set indices
!
                iside = pinfo%side
                iface = pinfo%face

! update boundaries
!
                if (pinfo%level_difference .eq. 1 .and. iface .eq. 1) then

                  pblock => pinfo%block
                  pneigh => pblock%neigh(idir,iside,iface)%ptr

                  nside = 3 - iside
                  nface = 1
                  do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                    nface = nface + 1
                  end do

                  call bnd_prol(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside &
                                                                       , nface)

                end if

                pinfo => pinfo%prev
                l = l + 1
              end do

            end if

! deallocate buffers
!
            deallocate(rbuf)

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
    end do ! directions
!
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

    implicit none

! local variables
!
    integer :: idir, iside, iface

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
!
!-------------------------------------------------------------------------------
!
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

! update directional flux from the neighbor
!
                  call correct_flux(pmeta%data, pneigh%data%f                  &
                                                         , idir, iside, iface)

                end if ! pmeta level < pneigh level

              end if ! pneigh associated

            end do ! iface
          end do ! iside
        end do ! idir

      end if ! leaf

      pmeta => pmeta%next ! assign pointer to the next meta block in the list
    end do ! meta blocks

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
    integer :: i, ip, is, it, ih, il, iu, i1, i2
    integer :: j, jp, js, jt, jh, jl, ju, j1, j2
#if NDIMS == 3
    integer :: k, kp, ks, kt, kh, kl, ku, k1, k2
#endif /* NDIMS == 3 */

! arguments
!
    type(block_data), pointer             , intent(inout) :: pdata
    real            , dimension(:,:,:,:,:), intent(in)    :: f
    integer                               , intent(in)    :: idir, iside, iface
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
        is = ie
        it = ibl
      else ! right side
        is = ibl
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

        pdata%f(idir,:,it,j,:) = 0.5d0 * (f(idir,:,is,j1,:) + f(idir,:,is,j2,:))
      end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do j = jl, ju
        j1 = 2 * (j - jl) + jb
        j2 = j1 + 1
        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,it,j,k) = (f(idir,:,is,j1,k1) + f(idir,:,is,j2,k1)    &
                           + f(idir,:,is,j1,k2) + f(idir,:,is,j2,k2)) / 0.25d0
        end do
      end do
#endif /* NDIMS == 3 */

! Y direction
!
    case(2)

! index of the slice which will be updated
!
      if (iside .eq. 1) then ! left side
        js = je
        jt = jbl
      else ! right side
        js = jbl
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

        pdata%f(idir,:,i,jt,:) = 0.5d0 * (f(idir,:,i1,js,:) + f(idir,:,i2,js,:))
      end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
      do i = il, iu
        i1 = 2 * (i - il) + ib
        i2 = i1 + 1

        do k = kl, ku
          k1 = 2 * (k - kl) + kb
          k2 = k1 + 1

          pdata%f(idir,:,i,jt,k) = (f(idir,:,i1,js,k1) + f(idir,:,i2,js,k1)    &
                           + f(idir,:,i1,js,k2) + f(idir,:,i2,js,k2)) * 0.25d0
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
        ks = ke
        kt = kbl
      else ! right side
        ks = kbl
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

          pdata%f(idir,:,i,j,kt) = (f(idir,:,i1,j1,ks) + f(idir,:,i2,j1,ks)    &
                           + f(idir,:,i1,j2,ks) + f(idir,:,i2,j2,ks)) * 0.25d0
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
! bnd_copy: subroutine copies the interior of neighbor to update the
!           boundaries of the current block
!
!===============================================================================
!
  subroutine bnd_copy(pdata, u, idir, iside)

    use blocks   , only : block_data
    use config   , only : im, ib, ibl, ibu, ie, iel, ieu                       &
                        , jm, jb, jbl, jbu, je, jel, jeu                       &
                        , km, kb, kbl, kbu, ke, kel, keu
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! arguments
!
    type(block_data), pointer       , intent(inout) :: pdata
    real   , dimension(nqt,im,jm,km), intent(in)    :: u
    integer                         , intent(in)    :: idir, iside
!
!-------------------------------------------------------------------------------
!
    select case(idir)

    case(1)

      if (iside .eq. 1) then
        pdata%u(  1:nfl,1:ibl  ,1:jm,1:km) = u(  1:nfl,iel:ie  ,1:jm,1:km)
#ifdef MHD
        pdata%u(ibx:ibz,1:ibl  ,1:jm,1:km) = u(ibx:ibz,iel:ie  ,1:jm,1:km)
#ifdef GLM
        pdata%u(iph    ,1:ibl  ,1:jm,1:km) = u(iph    ,iel:ie  ,1:jm,1:km)
#endif /* GLM */
#endif /* MHD */
      else
        pdata%u(  1:nfl,ieu:im,1:jm,1:km) = u(  1:nfl,ib:ibu,1:jm,1:km)
#ifdef MHD
        pdata%u(ibx:ibz,ieu:im,1:jm,1:km) = u(ibx:ibz,ib:ibu,1:jm,1:km)
#ifdef GLM
        pdata%u(iph    ,ieu:im,1:jm,1:km) = u(iph    ,ib:ibu,1:jm,1:km)
#endif /* GLM */
#endif /* MHD */
      end if

    case(2)

      if (iside .eq. 1) then
        pdata%u(  1:nfl,1:im,1:jbl  ,1:km) = u(  1:nfl,1:im,jel:je  ,1:km)
#ifdef MHD
        pdata%u(ibx:ibz,1:im,1:jbl  ,1:km) = u(ibx:ibz,1:im,jel:je  ,1:km)
#ifdef GLM
        pdata%u(iph    ,1:im,1:jbl  ,1:km) = u(iph    ,1:im,jel:je  ,1:km)
#endif /* GLM */
#endif /* MHD */
      else
        pdata%u(  1:nfl,1:im,jeu:jm,1:km) = u(  1:nfl,1:im,jb:jbu,1:km)
#ifdef MHD
        pdata%u(ibx:ibz,1:im,jeu:jm,1:km) = u(ibx:ibz,1:im,jb:jbu,1:km)
#ifdef GLM
        pdata%u(iph    ,1:im,jeu:jm,1:km) = u(iph    ,1:im,jb:jbu,1:km)
#endif /* GLM */
#endif /* MHD */
      end if

#if NDIMS == 3
    case(3)

      if (iside .eq. 1) then
        pdata%u(  1:nfl,1:im,1:jm,1:kbl  ) = u(  1:nfl,1:im,1:jm,kel:ke  )
#ifdef MHD
        pdata%u(ibx:ibz,1:im,1:jm,1:kbl  ) = u(ibx:ibz,1:im,1:jm,kel:ke  )
#ifdef GLM
        pdata%u(iph    ,1:im,1:jm,1:kbl  ) = u(iph    ,1:im,1:jm,kel:ke  )
#endif /* GLM */
#endif /* MHD */
      else
        pdata%u(  1:nfl,1:im,1:jm,keu:km ) = u(  1:nfl,1:im,1:jm,kb:kbu)
#ifdef MHD
        pdata%u(ibx:ibz,1:im,1:jm,keu:km) = u(ibx:ibz,1:im,1:jm,kb:kbu)
#ifdef GLM
        pdata%u(iph    ,1:im,1:jm,keu:km) = u(iph    ,1:im,1:jm,kb:kbu)
#endif /* GLM */
#endif /* MHD */
      end if
#endif /* NDIMS == 3 */

    end select
!
!-------------------------------------------------------------------------------
!
  end subroutine bnd_copy
!
!===============================================================================
!
! bnd_rest: subroutine copies the restricted interior of the neighbor to update
!           the boundaries of the current block
!
!===============================================================================
!
  subroutine bnd_rest(pdata, u, idir, iside, iface)

    use blocks   , only : block_data
    use config   , only : ng, in, im, ih, ib, ibl, ibu, ie, iel, ieu           &
                        , nd, jn, jm, jh, jb, jbl, jbu, je, jel, jeu           &
                        , nh, kn, km, kh, kb, kbl, kbu, ke, kel, keu
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! arguments
!
    type(block_data), pointer       , intent(inout) :: pdata
    real   , dimension(nqt,im,jm,km), intent(in)    :: u
    integer                         , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: if, jf, kf, ip, jp, kp
    integer :: il, jl, kl, iu, ju, ku
    integer :: is, js, ks, it, jt, kt
!
!-------------------------------------------------------------------------------
!
#if NDIMS == 2
! prepare temporary indices
!
    kl   = 1
    ku   = 1
    ks   = 1
    kt   = 1
#endif /* NDIMS == 2 */

! prepare indices
!
    select case(idir)

    case(1)

      if =     iside - 1
      jf = mod(iface - 1, 2)
#if NDIMS == 3
      kf =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! indices of the source and destination arrays
!
      if (if .eq. 0) then
        il = ie - nd + 1
        iu = ie
        is = 1
        it = ibl
      else
        il = ib
        iu = ib + nd - 1
        is = ie + 1
        it = im
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

    case(2)

      if = mod(iface - 1, 2)
      jf =     iside - 1
#if NDIMS == 3
      kf =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! indices of the source and destination arrays
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
        jl = je - nd + 1
        ju = je
        js = 1
        jt = jbl
      else
        jl = jb
        ju = jb + nd - 1
        js = jeu
        jt = jm
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

#if NDIMS == 3
    case(3)

      if = mod(iface - 1, 2)
      jf =    (iface - 1) / 2
      kf =     iside - 1

! indices of the source and destination arrays
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
      if (kf .eq. 0) then
        kl = ke - nd + 1
        ku = ke
        ks = 1
        kt = kbl
      else
        kl = kb
        ku = kb + nd - 1
        ks = keu
        kt = km
      end if
      kp = kl + 1
#endif /* NDIMS == 3 */

    end select

! update field variable boundaries
!
#if NDIMS == 2
    pdata%u(1:nfl,is:it,js:jt,1) = 0.25 * (u(1:nfl,il:iu:2,jl:ju:2,1)          &
                                         + u(1:nfl,ip:iu:2,jl:ju:2,1)          &
                                         + u(1:nfl,il:iu:2,jp:ju:2,1)          &
                                         + u(1:nfl,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(1:nfl,is:it,js:jt,ks:kt) =                                         &
                                  0.125 * (u(1:nfl,il:iu:2,jl:ju:2,kl:ku:2)    &
                                         + u(1:nfl,ip:iu:2,jl:ju:2,kl:ku:2)    &
                                         + u(1:nfl,il:iu:2,jp:ju:2,kl:ku:2)    &
                                         + u(1:nfl,ip:iu:2,jp:ju:2,kl:ku:2)    &
                                         + u(1:nfl,il:iu:2,jl:ju:2,kp:ku:2)    &
                                         + u(1:nfl,ip:iu:2,jl:ju:2,kp:ku:2)    &
                                         + u(1:nfl,il:iu:2,jp:ju:2,kp:ku:2)    &
                                         + u(1:nfl,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#ifdef MHD

! update magnetic field components
!
#if NDIMS == 2
    pdata%u(ibx:ibz,is:it,js:jt,1) = 0.25d0 * (u(ibx:ibz,il:iu:2,jl:ju:2,1)    &
                                            +  u(ibx:ibz,ip:iu:2,jl:ju:2,1)    &
                                            +  u(ibx:ibz,il:iu:2,jp:ju:2,1)    &
                                            +  u(ibx:ibz,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(ibx:ibz,is:it,js:jt,ks:kt) =                                       &
                                0.125d0 * (u(ibx:ibz,il:iu:2,jl:ju:2,kl:ku:2)  &
                                        +  u(ibx:ibz,ip:iu:2,jl:ju:2,kl:ku:2)  &
                                        +  u(ibx:ibz,il:iu:2,jp:ju:2,kl:ku:2)  &
                                        +  u(ibx:ibz,ip:iu:2,jp:ju:2,kl:ku:2)  &
                                        +  u(ibx:ibz,il:iu:2,jl:ju:2,kp:ku:2)  &
                                        +  u(ibx:ibz,ip:iu:2,jl:ju:2,kp:ku:2)  &
                                        +  u(ibx:ibz,il:iu:2,jp:ju:2,kp:ku:2)  &
                                        +  u(ibx:ibz,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#ifdef GLM

! update the scalar potential
!
#if NDIMS == 2
    pdata%u(iph,is:it,js:jt,1) = 0.25d0 * (u(iph,il:iu:2,jl:ju:2,1)            &
                                        +  u(iph,ip:iu:2,jl:ju:2,1)            &
                                        +  u(iph,il:iu:2,jp:ju:2,1)            &
                                        +  u(iph,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(iph,is:it,js:jt,ks:kt) = 0.125d0 * (u(iph,il:iu:2,jl:ju:2,kl:ku:2) &
                                             +  u(iph,ip:iu:2,jl:ju:2,kl:ku:2) &
                                             +  u(iph,il:iu:2,jp:ju:2,kl:ku:2) &
                                             +  u(iph,ip:iu:2,jp:ju:2,kl:ku:2) &
                                             +  u(iph,il:iu:2,jl:ju:2,kp:ku:2) &
                                             +  u(iph,ip:iu:2,jl:ju:2,kp:ku:2) &
                                             +  u(iph,il:iu:2,jp:ju:2,kp:ku:2) &
                                             +  u(iph,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#endif /* GLM */
#endif /* MHD */
!
!-------------------------------------------------------------------------------
!
  end subroutine bnd_rest
!
!===============================================================================
!
! bnd_prol: subroutine copies the prolongated interior of the neighbor to update
!           the boundaries the of current block
!
!===============================================================================
!
  subroutine bnd_prol(pdata, ub, idir, iside, iface)

    use blocks       , only : block_data
    use config       , only : ng, in, im, ib, ibl, ibu, ie, iel, ieu           &
                            , jn, jm, jb, jbl, jbu, je, jel, jeu               &
                            , kn, km, kb, kbl, kbu, ke, kel, keu
    use error        , only : print_warning
    use interpolation, only : expand
    use variables    , only : nvr, nqt, nfl
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real       , dimension(nqt,im,jm,km), intent(in)    :: ub
    integer                             , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: q, i, j, k
    integer :: ih, jh, kh, nh
    integer :: if, jf, kf, id, jd, kd
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt

! local arrays
!
    integer, dimension(3)                    :: dm, cm, pm
    real   , dimension(:,:,:)  , allocatable :: u

! parameters
!
    integer :: dl = 2
!
!-------------------------------------------------------------------------------
!
! prepare useful variables
!
    nh = max(1, ng / 2)
    ih = max(1, in / 2)
    jh = max(1, jn / 2)
    kh = max(1, kn / 2)

! prepare dimensions of the input arrays
!
    dm(:)    = (/ im, jm, km /)
    dm(idir) = ng
    pm(:)    = (/ ih, jh, kh /) + ng
    pm(idir) = nh
    cm(:)    = pm(:) + 2 * dl
#if NDIMS == 2
    dm(3)    = 1
    pm(3)    = 1
    cm(3)    = 1

    kl = 1
    ku = 1
#endif /* NDIMS == 2 */

! prepare indices of for the destination array
!
    is = 1
    it = im
    js = 1
    jt = jm
    ks = 1
    kt = km

! prepare indices
!
    select case(idir)

    case(1)

      if =     iside - 1
      jf = mod(iface - 1, 2)
      kf =    (iface - 1) / 2

! indices of the input array
!
      if (if .eq. 0) then
        il = ie - nh - dl + 1
        iu = ie      + dl
      else
        il = ib      - dl
        iu = ib + nh + dl - 1
      end if
      jl = jb - nh - dl + jh * jf
      ju = jl + cm(2) - 1
#if NDIMS == 3
      kl = kb - nh - dl + kh * kf
      ku = kl + cm(3) - 1
#endif /* NDIMS == 3 */

! indices of the destination array
!
      if (if .eq. 0) then
        is = 1
        it = ibl
      else
        is = ieu
        it = im
      end if

    case(2)

      if = mod(iface - 1, 2)
      jf =     iside - 1
      kf =    (iface - 1) / 2

! indices of the input array
!
      il = ib - nh - dl + ih * if
      iu = il + cm(1) - 1
      if (jf .eq. 0) then
        jl = je - nh - dl + 1
        ju = je      + dl
      else
        jl = jb      - dl
        ju = jb + nh + dl - 1
      end if
#if NDIMS == 3
      kl = kb - nh - dl + kh * kf
      ku = kl + cm(3) - 1
#endif /* NDIMS == 3 */

! indices of the destination array
!
      if (jf .eq. 0) then
        js = 1
        jt = jbl
      else
        js = jeu
        jt = jm
      end if

#if NDIMS == 3
    case(3)

      if = mod(iface - 1, 2)
      jf =    (iface - 1) / 2
      kf =     iside - 1

! indices of the input array
!
      il = ib - nh - dl + ih * if
      iu = il + cm(1) - 1
      jl = jb - nh - dl + jh * jf
      ju = jl + cm(2) - 1
      if (kf .eq. 0) then
        kl = ke - nh - dl + 1
        ku = ke      + dl
      else
        kl = kb      - dl
        ku = kb + nh + dl - 1
      end if

! indices of the destination array
!
      if (kf .eq. 0) then
        ks = 1
        kt = kbl
      else
        ks = keu
        kt = km
      end if
#endif /* NDIMS == 3 */

    end select

! allocate temporary array
!
    allocate(u(dm(1),dm(2),dm(3)))

! iterate over all variables
!
    do q = 1, nfl

! expand the boundary
!
      call expand(cm(:), dm(:), dl, ub(q,il:iu,jl:ju,kl:ku), u(:,:,:))

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(:,:,:)

    end do
#ifdef MHD
! iterate over magnetic field components
!
    do q = ibx, ibz

! expand the boundary
!
      call expand(cm(:), dm(:), dl, ub(q,il:iu,jl:ju,kl:ku), u(:,:,:))

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(:,:,:)

    end do

#ifdef GLM
! expand and update the scalar potential
!
    call expand(cm(:), dm(:), dl, ub(iph,il:iu,jl:ju,kl:ku), u(:,:,:))

! copy expanded boundary in the proper place of the block
!
    pdata%u(iph,is:it,js:jt,ks:kt) = u(:,:,:)
#endif /* GLM */
#endif /* MHD */

! deallocate temporary array
!
    deallocate(u)
!
!-------------------------------------------------------------------------------
!
  end subroutine bnd_prol
!
!===============================================================================
!
! bnd_spec: subroutine to apply specific boundary conditions
!
!===============================================================================
!
  subroutine bnd_spec(pb, id, il)

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
    type(block_data), pointer, intent(inout) :: pb
    integer                  , intent(in)    :: id, il

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
    ii = 10 * id + il

! perform update according to the flag
!
    select case(ii)
    case(11)
      select case(xlbndry)
      case("reflecting")
        do i = 1, ng
          it = ib  - i
          is = ibl + i
          do k = 1, km
            do j = 1, jm
              pb%u(  :,it,j,k) =  pb%u(   :,is,j,k)
              pb%u(imx,it,j,k) = -pb%u(imx,is,j,k)
          end do
          end do
        end do
#ifdef MHD
      case("divb")
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do j = 1, jm
            jm1 = max( 1, j-1)
            jp1 = min(jm, j+1)
            do i = ng, 1, -1
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,ib,j,k)

              dbp = pb%u(iby,ib,jp1,k) - pb%u(iby,ib,j  ,k)
              dbm = pb%u(iby,ib,j  ,k) - pb%u(iby,ib,jm1,k)
              dby = limiter(dbp, dbm)
#if NDIMS == 3
              dbp = pb%u(ibz,ib,j,kp1) - pb%u(ibz,ib,j,k  )
              dbm = pb%u(ibz,ib,j,k  ) - pb%u(ibz,ib,j,km1)
              dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */
#if NDIMS == 2
              pb%u(ibx,i,j,k) = pb%u(ibx,i+1,j,k) + dby
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pb%u(ibx,i,j,k) = pb%u(ibx,i+1,j,k) + dby + dbz
#endif /* NDIMS == 3 */
              pb%u(iby,i,j,k) = pb%u(iby,ib,j,k)
              pb%u(ibz,i,j,k) = pb%u(ibz,ib,j,k)
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,ib,j,k)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = 1, km
          do j = 1, jm
            do i = 1, ng
              pb%u(:,i,j,k) = pb%u(:,ib,j,k)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
    case(12)
      select case(xubndry)
      case("reflecting")
        do i = 1, ng
          it = ie  + i
          is = ieu - i
          do k = 1, km
            do j = 1, jm
              pb%u(  :,it,j,k) =  pb%u(  :,is,j,k)
              pb%u(imx,it,j,k) = -pb%u(imx,is,j,k)
            end do
          end do
        end do
#ifdef MHD
      case("divb")
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do j = 1, jm
            jm1 = max( 1, j-1)
            jp1 = min(jm, j+1)
            do i = ieu, im
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,ie,j,k)

              dbp = pb%u(iby,ie,jp1,k) - pb%u(iby,ie,j  ,k)
              dbm = pb%u(iby,ie,j  ,k) - pb%u(iby,ie,jm1,k)
              dby = limiter(dbp, dbm)
#if NDIMS == 3
              dbp = pb%u(ibz,ie,j,kp1) - pb%u(ibz,ie,j,k  )
              dbm = pb%u(ibz,ie,j,k  ) - pb%u(ibz,ie,j,km1)
              dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

#if NDIMS == 2
              pb%u(ibx,i,j,k) = pb%u(ibx,i-1,j,k) - dby
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pb%u(ibx,i,j,k) = pb%u(ibx,i-1,j,k) - dby - dbz
#endif /* NDIMS == 3 */
              pb%u(iby,i,j,k) = pb%u(iby,ie,j,k)
              pb%u(ibz,i,j,k) = pb%u(ibz,ie,j,k)
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,ie,j,k)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = 1, km
          do j = 1, jm
            do i = ieu, im
              pb%u(:,i,j,k) = pb%u(:,ie,j,k)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
    case(21)
      select case(ylbndry)
      case("reflecting")
        do j = 1, ng
          jt = jb  - j
          js = jbl + j
          do k = 1, km
            do i = 1, im
              pb%u(  :,i,jt,k) =  pb%u(  :,i,js,k)
              pb%u(imy,i,jt,k) = -pb%u(imy,i,js,k)
            end do
          end do
        end do
#ifdef MHD
      case("divb")
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)
            do j = ng, 1, -1
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,i,jb,k)

              dbp = pb%u(ibx,ip1,jb,k) - pb%u(ibx,i  ,jb,k)
              dbm = pb%u(ibx,i  ,jb,k) - pb%u(ibx,im1,jb,k)
              dbx = limiter(dbp, dbm)
#if NDIMS == 3
              dbp = pb%u(ibz,i,jb,kp1) - pb%u(ibz,i,jb,k  )
              dbm = pb%u(ibz,i,jb,k  ) - pb%u(ibz,i,jb,km1)
              dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

              pb%u(ibx,i,j,k) = pb%u(ibx,i,jb,k)
#if NDIMS == 2
              pb%u(iby,i,j,k) = pb%u(iby,i,j+1,k) + dbx
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pb%u(iby,i,j,k) = pb%u(iby,i,j+1,k) + dbx + dbz
#endif /* NDIMS == 3 */
              pb%u(ibz,i,j,k) = pb%u(ibz,i,jb,k)
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,i,jb,k)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = 1, km
          do j = 1, ng
            do i = 1, im
              pb%u(:,i,j,k) = pb%u(:,i,jb,k)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
    case(22)
      select case(yubndry)
      case("reflecting")
        do j = 1, ng
          jt = je  + j
          js = jeu - j
          do k = 1, km
            do i = 1, im
              pb%u(  :,i,jt,k) =  pb%u(  :,i,js,k)
              pb%u(imy,i,jt,k) = -pb%u(imy,i,js,k)
            end do
          end do
        end do
#ifdef MHD
      case("divb")
        do k = 1, km
#if NDIMS == 3
          km1 = max( 1, k-1)
          kp1 = min(km, k+1)
#endif /* NDIMS == 3 */
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)
            do j = jeu, jm
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,i,je,k)

              dbp = pb%u(ibx,ip1,je,k) - pb%u(ibx,i  ,je,k)
              dbm = pb%u(ibx,i  ,je,k) - pb%u(ibx,im1,je,k)
              dbx = limiter(dbp, dbm)
#if NDIMS == 3
              dbp = pb%u(ibz,i,je,kp1) - pb%u(ibz,i,je,k  )
              dbm = pb%u(ibz,i,je,k  ) - pb%u(ibz,i,je,km1)
              dbz = limiter(dbp, dbm)
#endif /* NDIMS == 3 */

              pb%u(ibx,i,j,k) = pb%u(ibx,i,je,k)
#if NDIMS == 2
              pb%u(iby,i,j,k) = pb%u(iby,i,j-1,k) - dbx
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pb%u(iby,i,j,k) = pb%u(iby,i,j-1,k) - dbx - dbz
#endif /* NDIMS == 3 */
              pb%u(ibz,i,j,k) = pb%u(ibz,i,je,k)
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,i,je,k)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = 1, km
          do j = jeu, jm
            do i = 1, im
              pb%u(:,i,j,k) = pb%u(:,i,je,k)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
#if NDIMS == 3
    case(31)
      select case(zlbndry)
      case("reflecting")
        do k = 1, ng
          kt = kb  - k
          ks = kbl + k
          do j = 1, jm
            do i = 1, im
              pb%u(  :,i,j,kt) =  pb%u(  :,i,j,ks)
              pb%u(imz,i,j,kt) = -pb%u(imz,i,j,ks)
            end do
          end do
        end do
#ifdef MHD
      case("divb")
        do j = 1, jm
          jm1 = max( 1, j-1)
          jp1 = min(jm, j+1)
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)
            do k = ng, 1, -1
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,i,j,kb)

              dbp = pb%u(ibx,ip1,j,kb) - pb%u(ibx,i  ,j,kb)
              dbm = pb%u(ibx,i  ,j,kb) - pb%u(ibx,im1,j,kb)
              dbx = limiter(dbp, dbm)

              dbp = pb%u(iby,i,jp1,kb) - pb%u(iby,i,j  ,kb)
              dbm = pb%u(iby,i,j  ,kb) - pb%u(iby,i,jm1,kb)
              dbz = limiter(dbp, dbm)

              pb%u(ibx,i,j,k) = pb%u(ibx,i,j,kb)
              pb%u(iby,i,j,k) = pb%u(iby,i,j,kb)
              pb%u(ibz,i,j,k) = pb%u(ibz,i,j,k+1) + dbx + dbz
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,i,j,kb)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = 1, ng
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,k) = pb%u(:,i,j,kb)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
    case(32)
      select case(zubndry)
      case("reflecting")
        do k = 1, ng
          kt = ke  + k
          ks = keu - k
          do j = 1, jm
            do i = 1, im
              pb%u(  :,i,j,kt) =  pb%u(  :,i,j,ks)
              pb%u(imz,i,j,kt) = -pb%u(imz,i,j,ks)
            end do
          end do
        end do
#ifdef MHD
      case("divb")
        do j = 1, jm
          jm1 = max( 1, j-1)
          jp1 = min(jm, j+1)
          do i = 1, im
            im1 = max( 1, i-1)
            ip1 = min(im, i+1)
            do k = keu, km
              pb%u(1:nfl,i,j,k) = pb%u(1:nfl,i,j,ke)

              dbp = pb%u(ibx,ip1,j,ke) - pb%u(ibx,i  ,j,ke)
              dbm = pb%u(ibx,i  ,j,ke) - pb%u(ibx,im1,j,ke)
              dbx = limiter(dbp, dbm)

              dbp = pb%u(iby,i,jp1,ke) - pb%u(iby,i,j  ,ke)
              dbm = pb%u(iby,i,j  ,ke) - pb%u(iby,i,jm1,ke)
              dbz = limiter(dbp, dbm)

              pb%u(ibx,i,j,k) = pb%u(ibx,i,j,ke)
              pb%u(iby,i,j,k) = pb%u(iby,i,j,ke)
              pb%u(ibz,i,j,k) = pb%u(ibz,i,j,k-1) - dbx - dbz
#ifdef GLM
              pb%u(iph,i,j,k) = pb%u(iph,i,j,ke)
#endif /* GLM */
            end do
          end do
        end do
#endif /* MHD */
      case default ! set "open" for default
        do k = keu, km
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,k) = pb%u(:,i,j,ke)
#if defined MHD && defined GLM
              pb%u(iph,i,j,k) = 0.0d0
#endif /* MHD & GLM */
            end do
          end do
        end do
      end select
#endif /* NDIMS == 3 */
    case default
      call print_warning("boundaries::bnd_spec", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_spec

!===============================================================================
!
end module
