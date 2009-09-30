!!*****************************************************************************
!!
!! module: boundaries - routines for handling the boundary conditions
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
module boundaries

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! boundary: subroutine sweeps over all leaf blocks and performs the boundary
!           update
!
!===============================================================================
!
  subroutine boundary

    use config  , only : im, jm, km
    use blocks  , only : nv => nvars, ndims, nsides, nfaces, nchild            &
                       , block_meta, block_info, pointer_info, list_meta
    use error   , only : print_error
    use mpitools, only : ncpus, ncpu, msendf, mrecvf

    implicit none

! local variables
!
    integer :: i, j, k, l, m, p, q, dl

#ifdef MPI
! MPI tag
!
    integer(kind=4)                                       :: itag

! pointer to the list of pairs of blocks belonging to different processes for
! which we perform the boundary exchange
!
    type(pointer_info), dimension(0:ncpus-1,0:ncpus-1)    :: block_array

! declare pointer to the info block
!
    type(block_info), pointer                             :: pinfo

! create counter for exchanged blocks
!
    integer(kind=4)   , dimension(0:ncpus-1,0:ncpus-1)    :: block_counter

! local arrays for boundary exchange
!
    real(kind=8)      , dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */

! local pointers
!
    type(block_meta), pointer :: pblock, pneigh
!
!-------------------------------------------------------------------------------
!
#ifdef MPI
! nullify array pointers and reset the counter
!
    do p = 0, ncpus-1
      do q = 0, ncpus-1
        nullify(block_array(p,q)%ptr)
        block_counter(p,q) = 0
      end do
    end do
#endif /* MPI */

! iterate over all leaf blocks and perform boundary update
!
    pblock => list_meta
    do while (associated(pblock))

! if the current block is a leaf...
!
      if (pblock%leaf) then

! iterate over all neighbors
!
        do i = 1, ndims
          do j = 1, nsides
            do k = 1, nfaces

! associate pointer to the neighbor
!
              pneigh => pblock%neigh(i,j,k)%ptr

! check if neighbor is associated, if yes exchange boundaries, if not call
! specific boundary conditions
!
              if (associated(pneigh)) then

#ifdef MPI
! get the processor number of the current block and its neighbor
!
                p = pblock%cpu
                q = pneigh%cpu

! check if the neighbor is on the same cpu
!
                if (p .eq. q) then
                  if (p .eq. ncpu) then
#endif /* MPI */

! calculate the difference of current and neighbor level
!
                    dl = pblock%level - pneigh%level

! depending on the level difference
!
                    select case(dl)
                    case(-1)  ! restriction
                      call bnd_rest(pblock%data, pneigh%data%u, i, j, k)
                    case(0)   ! the same level, copying
                      if (k .eq. 1) &
                        call bnd_copy(pblock%data, pneigh%data%u, i, j, k)
                    case(1)   ! prolongation
                      m = 1
                      do while(pblock%id .ne. pneigh%neigh(i,3-j,m)%ptr%id .and. m .lt. nfaces)
                        m = m + 1
                      end do
                      if (m .le. nfaces) then
                        if (k .eq. 1) &
                          call bnd_prol(pblock%data, pneigh%data%u, i, j, m)
                      else
                        call print_error("boundaries::boundary", "Index m out of the limit!")
                      end if
                    case default
                      call print_error("boundaries::boundary", "Level difference unsupported!")
                   end select

#ifdef MPI
                  end if
                else

! neighbor is on another cpu, which means we have to update it later, so add it
! to the list of blocks to exchange
!
! allocate info block
!
                  allocate(pinfo)

! fill out its fields
!
                  pinfo%block            => pblock
                  pinfo%neigh            => pneigh
                  pinfo%direction        =  i
                  pinfo%side             =  j
                  pinfo%face             =  k
                  pinfo%level_difference = pblock%level - pneigh%level

! nullify pointers
!
                  nullify(pinfo%prev)
                  nullify(pinfo%next)

! if the list is not emply append the created block
!
                  if (associated(block_array(p,q)%ptr)) then
                    pinfo%prev => block_array(p,q)%ptr
                    nullify(pinfo%next)
                  end if

! point the list to the last created block
!
                  block_array(p,q)%ptr => pinfo

! increase the counter
!
                  block_counter(p,q) = block_counter(p,q) + 1
                end if
#endif /* MPI */
              else

! neighbor is not associated, it means that we have non periodic boundary here
!
#ifdef MPI
                if (k .eq. 1 .and. pblock%cpu .eq. ncpu) &
                  call bnd_spec(pblock%data, i, j, k)
#else /* MPI */
                if (k .eq. 1) &
                  call bnd_spec(pblock%data, i, j, k)
#endif /* MPI */

              endif
            end do
          end do
        end do

      endif

! assign pointer to the next block
!
      pblock => pblock%next

    end do

#ifdef MPI
!! TO DO:
!!
!! - send/receive only boundaries, not whole blocks
!!
! iterate over receiving and sending processors
!
    do p = 0, ncpus - 1
      do q = 0, ncpus - 1

        if (block_counter(p,q) .gt. 0) then

! get the number of nodes to exchange
!
           l = block_counter(p,q)

! get the tag for communication
!
          itag = p * ncpus + q + ncpus + 1

! allocate space for variables
!
          allocate(rbuf(l,nv,im,jm,km))

! if p == ncpu we are receiving data
!
          if (p .eq. ncpu) then

! receive data
!
            call mrecvf(size(rbuf(:,:,:,:,:)), q, itag, rbuf(:,:,:,:,:))

! iterate over all received blocks and update boundaries
!
            l = 1

            pinfo => block_array(p,q)%ptr
            do while(associated(pinfo))

! set indices
!
              i  = pinfo%direction
              j  = pinfo%side
              k  = pinfo%face
              dl = pinfo%level_difference

! update boundaries
!
              select case(dl)
              case(-1)  ! restriction
                call bnd_rest(pinfo%block%data,rbuf(l,:,:,:,:),i,j,k)
              case(0)   ! the same level, copying
                if (k .eq. 1) &
                  call bnd_copy(pinfo%block%data,rbuf(l,:,:,:,:),i,j,k)
              case(1)   ! prolongation
                pblock => pinfo%block
                pneigh => pblock%neigh(i,j,k)%ptr
                m = 1
                do while(pblock%id .ne. pneigh%neigh(i,3-j,m)%ptr%id .and. m .lt. nfaces)
                  m = m + 1
                end do
                if (m .le. nfaces) then
                  if (k .eq. 1) &
                    call bnd_prol(pinfo%block%data,rbuf(l,:,:,:,:),i,j,m)
                  else
                    call print_error("boundaries::boundary", "Index m out of the limit!")
                end if
              case default
                call print_error("boundaries::boundary", "Level difference unsupported!")
              end select

              pinfo => pinfo%prev
              l = l + 1
            end do

          end if

! if q == ncpu we are sending data
!
          if (q .eq. ncpu) then

! fill out the buffer with block data
!
            l = 1

            pinfo => block_array(p,q)%ptr
            do while(associated(pinfo))

              rbuf(l,:,:,:,:) = pinfo%neigh%data%u(:,:,:,:)

              pinfo => pinfo%prev
              l = l + 1
            end do

! send data buffer
!
            call msendf(size(rbuf), p, itag, rbuf)

          end if

! deallocate buffers
!
          deallocate(rbuf)

! deallocate info blocks
!
            pinfo => block_array(p,q)%ptr
            do while(associated(pinfo))
              block_array(p,q)%ptr => pinfo%prev

              nullify(pinfo%prev)
              nullify(pinfo%next)
              nullify(pinfo%block)
              nullify(pinfo%neigh)

              deallocate(pinfo)

              pinfo => block_array(p,q)%ptr
            end do

        end if

      end do
    end do
#endif /* MPI */

!-------------------------------------------------------------------------------
!
  end subroutine boundary
!
!===============================================================================
!
! bnd_copy: subroutine copies the interior of neighbor to update the
!           boundaries of current block
!
!===============================================================================
!
  subroutine bnd_copy(pdata, un, id, is, ip)

    use blocks, only : block_data, nv => nvars
    use config, only : im, ib, ibl, ibu, ie, iel, ieu                          &
                     , jm, jb, jbl, jbu, je, jel, jeu                          &
                     , km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: un
    integer                             , intent(in)    :: id, is, ip

! local variables
!
    integer :: ii, i, j, k
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * is

! perform update according to the flag
!
    select case(ii)

    case(110)
      do k = 1, km
        do j = 1, jm
          pdata%u(1:nv,1:ibl,j,k) = un(1:nv,iel:ie,j,k)
        end do
      end do

    case(120)
      do k = 1, km
        do j = 1, jm
          pdata%u(1:nv,ieu:im,j,k) = un(1:nv,ib:ibu,j,k)
        end do
      end do

    case(210)
      do k = 1, km
        do i = 1, im
          pdata%u(1:nv,i,1:jbl,k) = un(1:nv,i,jel:je,k)
        end do
      end do

    case(220)
      do k = 1, km
        do i = 1, im
          pdata%u(1:nv,i,jeu:jm,k) = un(1:nv,i,jb:jbu,k)
        end do
      end do

    case(310)
      do j = 1, jm
        do i = 1, im
          pdata%u(1:nv,i,j,1:kbl) = un(1:nv,i,j,kel:ke)
        end do
      end do

    case(320)
      do j = 1, jm
        do i = 1, im
          pdata%u(1:nv,i,j,keu:km) = un(1:nv,i,j,kb:kbu)
        end do
      end do

    case default
      call print_warning("boundaries::bnd_copy", "Boundary flag unsupported!")

    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_copy
!
!===============================================================================
!
! bnd_rest: subroutine copies the interior of neighbor to update the boundaries
!           of current block
!
!===============================================================================
!
  subroutine bnd_rest(pdata, un, id, is, ip)

    use blocks, only : block_data, nv => nvars, nchild
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                 &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                 &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: un
    integer                             , intent(in)    :: id, is, ip

! local variables
!
    integer :: ii, q, i, j, k
    integer :: il, iu, jl, ju, kl, ku
    integer :: i1, i2, j1, j2, k1, k2
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * is + ip

! perform update according to the flag
!
    select case(ii)

!! X-direction
!!
    case(111)

      il = 1
      iu = ng
      jl = ng / 2 + 1
      ju = jm / 2
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

    case(112)

      il = 1
      iu = ng
      jl = jm / 2 + 1
      ju = jm - ng / 2
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(113)

      il = 1
      iu = ng
      jl = ng / 2 + 1
      ju = jm / 2
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(114)

      il = 1
      iu = ng
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do
#endif /* NDIMS == 3 */


    case(121)

      il = ieu
      iu = im
      jl = ng / 2 + 1
      ju = jm / 2
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

    case(122)

      il = ieu
      iu = im
      jl = jm / 2 + 1
      ju = jm - ng / 2
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(123)

      il = ieu
      iu = im
      jl = ng / 2 + 1
      ju = jm / 2
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(124)

      il = ieu
      iu = im
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do
#endif /* NDIMS == 3 */

!! Y-direction
!!
    case(211)

      il = ng / 2 + 1
      iu = im / 2
      jl = 1
      ju = ng
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = je + 2 * (j - ng)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

    case(212)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = 1
      ju = ng
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = je + 2 * (j - ng)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(213)

      il = ng / 2 + 1
      iu = im / 2
      jl = 1
      ju = ng
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = je + 2 * (j - ng)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(214)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = 1
      ju = ng
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = je + 2 * (j - ng)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do
#endif /* NDIMS == 3 */

    case(221)

      il = ng / 2 + 1
      iu = im / 2
      jl = jeu
      ju = jm
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = ng + 2 * (j - je)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

    case(222)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = jeu
      ju = jm
#if NDIMS == 3
      kl = ng / 2 + 1
      ku = km / 2
#endif /* NDIMS == 3 */
#if NDIMS == 2
      k  = 1
      k1 = 1
      k2 = 1
#endif /* NDIMS == 2 */

#if NDIMS == 3
      do k = kl, ku
        k2 = 2 * k - ng
        k1 = k2 - 1
#endif /* NDIMS == 3 */
        do j = jl, ju
          j2 = ng + 2 * (j - je)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
#if NDIMS == 3
      end do
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(223)

      il = ng / 2 + 1
      iu = im / 2
      jl = jeu
      ju = jm
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = ng + 2 * (j - je)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(224)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = jeu
      ju = jm
      kl = km / 2 + 1
      ku = km - ng / 2

      do k = kl, ku
        k2 = 2 * k - km + ng
        k1 = k2 - 1
        do j = jl, ju
          j2 = ng + 2 * (j - je)
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do


!! Z-direction
!!
    case(311)

      il = ng / 2 + 1
      iu = im / 2
      jl = ng / 2 + 1
      ju = jm / 2
      kl = 1
      ku = ng

      do k = kl, ku
        k2 = ke + 2 * (k - ng)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(312)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = ng / 2 + 1
      ju = jm / 2
      kl = 1
      ku = ng

      do k = kl, ku
        k2 = ke + 2 * (k - ng)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(313)

      il = ng / 2 + 1
      iu = im / 2
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = 1
      ku = ng

      do k = kl, ku
        k2 = ke + 2 * (k - ng)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(314)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = 1
      ku = ng

      do k = kl, ku
        k2 = ke + 2 * (k - ng)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(321)

      il = ng / 2 + 1
      iu = im / 2
      jl = ng / 2 + 1
      ju = jm / 2
      kl = keu
      ku = km

      do k = kl, ku
        k2 = ng + 2 * (k - ke)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(322)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = ng / 2 + 1
      ju = jm / 2
      kl = keu
      ku = km

      do k = kl, ku
        k2 = ng + 2 * (k - ke)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(323)

      il = ng / 2 + 1
      iu = im / 2
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = keu
      ku = km

      do k = kl, ku
        k2 = ng + 2 * (k - ke)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do

    case(324)

      il = im / 2 + 1
      iu = im - ng / 2
      jl = jm / 2 + 1
      ju = jm - ng / 2
      kl = keu
      ku = km

      do k = kl, ku
        k2 = ng + 2 * (k - ke)
        k1 = k2 - 1
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = il, iu
            i2 = 2 * i - im + ng
            i1 = i2 - 1

            do q = 1, nv
              pdata%u(q,i,j,k) = sum(un(q,i1:i2,j1:j2,k1:k2)) / nchild
            end do
          end do
        end do
      end do
#endif /* NDIMS == 3 */

    case default
      call print_warning("boundaries::bnd_rest", "Boundary flag unsupported!")

    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_rest
!
!===============================================================================
!
! bnd_prol: subroutine copies the interior of neighbor to update the boundaries
!           of current block
!
!===============================================================================
!
  subroutine bnd_prol(pdata, un, idir, iside, iface)

    use blocks, only : block_data, nv => nvars
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                  &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                  &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : expand

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: un
    integer                             , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: ii, i, j, k, q
    integer :: i1, i2, j1, j2, k1, k2
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt

! local arrays
!
    integer, dimension(3)              :: dm, cm, fm
    real   , dimension(2*im,2*jm,2*km) :: ux

! parameters
!
    integer :: del = 1
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * idir + 10 * iside + iface

! prepare dimensions
!
    dm(:)  = (/ im, jm, km /)
    dm(idir) = ng
    cm(:)  = dm(:) / 2 + 2 * del
    cm(idir) = ng / 2 + 2 * del
    fm(:)  = 2 * cm(:)
#if NDIMS == 2
    dm(3)  = 1
    cm(3)  = 1
    fm(3)  = 1
#endif /* NDIMS == 2 */

! prepare lower indices of the source array
!
    i1 = (im - ng) / 2 - del + 1
    j1 = (jm - ng) / 2 - del + 1
#if NDIMS == 3
    k1 = (km - ng) / 2 - del + 1
#else /* NDIMS == 3 */
    k1 = 1
#endif /* NDIMS == 3 */

    il = 2 * del + 1
    jl = 2 * del + 1
#if NDIMS == 3
    kl = 2 * del + 1
#else /* NDIMS == 3 */
    kl = 1
#endif /* NDIMS == 3 */

    is = 1
    js = 1
    ks = 1

! update lower indices accoring to the flag
!
    select case(ii)

    case(111)

      i1 =  ie - ng  / 2 - del + 1
      j1 =  jb - ng  / 2 - del
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

    case(112)

      i1 =  ie - ng  / 2 - del + 1
      j1 = (jm - ng) / 2 - del + 1
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(113)

      i1 =  ie - ng  / 2 - del + 1
      j1 =  jb - ng  / 2 - del
      k1 = (km - ng) / 2 - del + 1

    case(114)

      i1 =  ie - ng  / 2 - del + 1
      j1 = (jm - ng) / 2 - del + 1
      k1 = (km - ng) / 2 - del + 1
#endif /* NDIMS == 3 */

    case(121)

      i1 =  ib - del
      j1 =  jb - ng  / 2 - del
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

      is = ieu

    case(122)

      i1 =  ib - del
      j1 = (jm - ng) / 2 - del + 1
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

      is = ieu

#if NDIMS == 3
    case(123)

      i1 =  ib - del
      j1 =  jb - ng  / 2 - del
      k1 = (km - ng) / 2 - del + 1

      is = ieu

    case(124)

      i1 =  ib - del
      j1 = (jm - ng) / 2 - del + 1
      k1 = (km - ng) / 2 - del + 1

      is = ieu
#endif /* NDIMS == 3 */

    case(211)

      i1 =  ib - ng  / 2 - del
      j1 =  je - ng  / 2 - del + 1
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

    case(212)

      i1 = (im - ng) / 2 - del + 1
      j1 =  je - ng  / 2 - del + 1
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(213)

      i1 =  ib - ng  / 2 - del
      j1 =  je - ng  / 2 - del + 1
      k1 = (km - ng) / 2 - del + 1

    case(214)

      i1 = (im - ng) / 2 - del + 1
      j1 =  je - ng  / 2 - del + 1
      k1 = (km - ng) / 2 - del + 1
#endif /* NDIMS == 3 */

    case(221)

      i1 =  ib - ng  / 2 - del
      j1 =  jb - del
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

      js = jeu

    case(222)

      i1 = (im - ng) / 2 - del + 1
      j1 =  jb - del
#if NDIMS == 3
      k1 =  kb - ng  / 2 - del
#else /* NDIMS == 3 */
      k1 = 1
#endif /* NDIMS == 3 */

      js = jeu

#if NDIMS == 3
    case(223)

      i1 =  ib - ng  / 2 - del
      j1 =  jb - del
      k1 = (km - ng) / 2 - del + 1

      js = jeu

    case(224)

      i1 = (im - ng) / 2 - del + 1
      j1 =  jb - del
      k1 = (km - ng) / 2 - del + 1

      js = jeu

    case(311)

      i1 =  ib - ng  / 2 - del
      j1 =  jb - ng  / 2 - del
      k1 =  ke - ng  / 2 - del + 1

    case(312)

      i1 = (im - ng) / 2 - del + 1
      j1 =  jb - ng  / 2 - del
      k1 =  ke - ng  / 2 - del + 1

    case(313)

      i1 =  ib - ng  / 2 - del
      j1 = (jm - ng) / 2 - del + 1
      k1 =  ke - ng  / 2 - del + 1

    case(314)

      i1 = (im - ng) / 2 - del + 1
      j1 = (jm - ng) / 2 - del + 1
      k1 =  ke - ng  / 2 - del + 1

    case(321)

      i1 =  ib - ng  / 2 - del
      j1 =  jb - ng  / 2 - del
      k1 =  kb - del

      ks = keu

    case(322)

      i1 = (im - ng) / 2 - del + 1
      j1 =  jb - ng  / 2 - del
      k1 =  kb - del

      ks = keu

    case(323)

      i1 =  ib - ng  / 2 - del
      j1 = (jm - ng) / 2 - del + 1
      k1 =  kb - del

      ks = keu

    case(324)

      i1 = (im - ng) / 2 - del + 1
      j1 = (jm - ng) / 2 - del + 1
      k1 =  kb - del

      ks = keu

#endif /* NDIMS == 3 */

    case default
      call print_warning("boundaries::bnd_prol", "Boundary flag unsupported!")

    end select

! update upper indices
!
    i2 = i1 + cm(1) - 1
    j2 = j1 + cm(2) - 1
    k2 = k1 + cm(3) - 1

    iu = il + dm(1) - 1
    ju = jl + dm(2) - 1
    ku = kl + dm(3) - 1

    it = is + dm(1) - 1
    jt = js + dm(2) - 1
    kt = ks + dm(3) - 1

! iterate over all variables
!
    do q = 1, nv

! expand the boundary
!
      call expand(cm, fm, 0, un(q,i1:i2,j1:j2,k1:k2), ux(:,:,:), 't', 't', 't')

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = ux(il:iu,jl:ju,kl:ku)

    end do

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
  subroutine bnd_spec(pb, id, il, ip)

    use blocks, only : block_data, nv => nvars
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , ng, im, jm, km, ib, ibl, ie, ieu, jb, jbl, je, jeu      &
                     , kb, kbl, ke, keu
    use error , only : print_warning
    use interpolation, only : expand

    implicit none

! arguments
!
    type(block_data), pointer, intent(inout) :: pb
    integer                  , intent(in)    :: id, il, ip

! local variables
!
    integer :: ii, i, j, k, it, jt, kt, is, js, ks, itm1, itp1, jtm1, jtp1, ktm1, ktp1
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * il

! perform update according to the flag
!
    select case(ii)
    case(110)
      select case(xlbndry)
      case("reflecting")
        do i = 1, ng
          it = ib  - i
          is = ibl + i
          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) =  pb%u(:,is,j,k)
              pb%u(2,it,j,k) = -pb%u(2,is,j,k)
          end do
          end do
        end do
      case default ! set "open" for default
        do i = 1, ng
          it   = ib - i
          itp1 = it + 1

          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) = pb%u(:,itp1,j,k)
            end do
          end do
        end do
      end select
    case(120)
      select case(xubndry)
      case("reflecting")
        do i = 1, ng
          it = ie  + i
          is = ieu - i
          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) =  pb%u(:,is,j,k)
              pb%u(2,it,j,k) = -pb%u(2,is,j,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do i = 1, ng
          it   = ie + i
          itm1 = it - 1

          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) = pb%u(:,itm1,j,k)
            end do
          end do
        end do
      end select
    case(210)
      select case(ylbndry)
      case("reflecting")
        do j = 1, ng
          jt = jb  - j
          js = jbl + j
          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) =  pb%u(:,i,js,k)
              pb%u(3,i,jt,k) = -pb%u(3,i,js,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do j = 1, ng
          jt   = jb - j
          jtp1 = jt + 1

          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) = pb%u(:,i,jtp1,k)
            end do
          end do
        end do
      end select
    case(220)
      select case(yubndry)
      case("reflecting")
        do j = 1, ng
          jt = je  + j
          js = jeu - j
          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) =  pb%u(:,i,js,k)
              pb%u(3,i,jt,k) = -pb%u(3,i,js,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do j = 1, ng
          jt   = je + j
          jtm1 = jt - 1

          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) = pb%u(:,i,jtm1,k)
            end do
          end do
        end do
      end select
    case(310)
      select case(zlbndry)
      case("reflecting")
        do k = 1, ng
          kt = kb  - k
          ks = kbl + k
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) =  pb%u(:,i,j,ks)
              pb%u(4,i,j,kt) = -pb%u(4,i,j,ks)
            end do
          end do
        end do
      case default ! set "open" for default
        do k = 1, ng
          kt   = kb - k
          ktp1 = kt + 1

          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) = pb%u(:,i,j,ktp1)
            end do
          end do
        end do
      end select
    case(320)
      select case(zubndry)
      case("reflecting")
        do k = 1, ng
          kt = ke  + k
          ks = keu - k
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) =  pb%u(:,i,j,ks)
              pb%u(4,i,j,kt) = -pb%u(4,i,j,ks)
            end do
          end do
        end do
      case default ! set "open" for default
        do k = 1, ng
          kt   = ke + k
          ktm1 = kt - 1

          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) = pb%u(:,i,j,ktm1)
            end do
          end do
        end do
      end select
    case default
      call print_warning("boundaries::bnd_spec", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_spec

!===============================================================================
!
end module
