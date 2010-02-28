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
    use blocks  , only : nvr, nqt, ndims, nsides, nfaces, nchild               &
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
! iterate over all directions
!
    do i = 1, ndims

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

              end if
            end do
          end do

        end if

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
            allocate(rbuf(l,nqt,im,jm,km))

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
    end do

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
  subroutine bnd_copy(pdata, u, idir, is, ip)

    use blocks, only : block_data, nvr, nqt
    use config, only : im, ib, ibl, ibu, ie, iel, ieu                          &
                     , jm, jb, jbl, jbu, je, jel, jeu                          &
                     , km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning

    implicit none

! arguments
!
    type(block_data), pointer            , intent(inout) :: pdata
    real        , dimension(nqt,im,jm,km), intent(in)    :: u
    integer                              , intent(in)    :: idir, is, ip

! local variables
!
    integer :: ii
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * idir + 10 * is

! perform update according to the flag
!
    select case(ii)

    case(110)
      pdata%u(1:nqt,  1:ibl,1:jm,1:km) = u(1:nqt,iel:ie ,1:jm,1:km)

    case(120)
      pdata%u(1:nqt,ieu:im ,1:jm,1:km) = u(1:nqt, ib:ibu,1:jm,1:km)

    case(210)
      pdata%u(1:nqt,1:im,  1:jbl,1:km) = u(1:nqt,1:im,jel:je ,1:km)

    case(220)
      pdata%u(1:nqt,1:im,jeu:jm ,1:km) = u(1:nqt,1:im, jb:jbu,1:km)

#if NDIMS == 3
    case(310)
      pdata%u(1:nqt,1:im,1:jm,  1:kbl) = u(1:nqt,1:im,1:jm,kel:ke )

    case(320)
      pdata%u(1:nqt,1:im,1:jm,keu:km ) = u(1:nqt,1:im,1:jm, kb:kbu)
#endif /* NDIMS == 3 */
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
  subroutine bnd_rest(pdata, un, idir, iside, iface)

    use blocks, only : block_data, nvr, nchild, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                 &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                 &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : shrink

    implicit none

! arguments
!
    type(block_data), pointer            , intent(inout) :: pdata
    real(kind=8), dimension(nqt,im,jm,km), intent(in)    :: un
    integer                              , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: ii, q, i, j, k
    integer :: il, iu, jl, ju, kl, ku
    integer :: i1, i2, j1, j2, k1, k2
    integer :: is, it, js, jt, ks, kt

! local arrays
!
    integer, dimension(3)                  :: dm, cm, fm, pm
    real   , dimension(:,:,:), allocatable :: ux

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
    dm(:)    = (/ im, jm, km /)
    dm(idir) = ng
    cm(:)    = dm(:) / 2
    cm(idir) = ng + 2 * del
    fm(:)    = 2 * cm(:)
    pm(:)    = dm(:) / 2 - ng / 2
    pm(idir) = ng
#if NDIMS == 2
    dm(3)  = 1
    cm(3)  = 1
    fm(3)  = 1
    pm(3)  = 1
    ks     = 1
#endif /* NDIMS == 2 */

! initiate lower indices
!
    i1 = 1
    j1 = 1
    k1 = 1

    il = 1
    jl = 1
    kl = 1

    is = ng / 2 + 1
    js = ng / 2 + 1
#if NDIMS == 3
    ks = ng / 2 + 1
#endif /* NDIMS == 3 */

! prepare lower indices depending on the direction and side of the boundary
!
    select case(ii)

!! X-direction
!!
    case(111)

      i1 = ie - 2 * (ng + del) + 1
      il =  1 + del
      is =  1

    case(112)

      i1 = ie - 2 * (ng + del) + 1
      il =  1 + del
      jl = ng / 2 + 1
      is =  1
      js = jm / 2 + 1

#if NDIMS == 3
    case(113)

      i1 = ie - 2 * (ng + del) + 1
      il =  1 + del
      kl = ng / 2 + 1
      is =  1
      ks = km / 2 + 1

    case(114)

      i1 = ie - 2 * (ng + del) + 1
      il =  1 + del
      jl = ng / 2 + 1
      kl = ng / 2 + 1
      is =  1
      js = jm / 2 + 1
      ks = km / 2 + 1
#endif /* NDIMS == 3 */

    case(121)

      i1 = ib - 2 * del
      il =  1 + del
      is = ieu

    case(122)

      i1 = ib - 2 * del
      il =  1 + del
      jl = ng / 2 + 1
      is = ieu
      js = jm / 2 + 1

#if NDIMS == 3
    case(123)

      i1 = ib - 2 * del
      il =  1 + del
      kl = ng / 2 + 1
      is = ieu
      ks = km / 2 + 1

    case(124)

      i1 = ib - 2 * del
      il =  1 + del
      jl = ng / 2 + 1
      kl = ng / 2 + 1
      is = ieu
      js = jm / 2 + 1
      ks = km / 2 + 1
#endif /* NDIMS == 3 */


!! Y-direction
!!
    case(211)

      j1 = je - 2 * (ng + del) + 1
      jl =  1 + del
      js =  1

    case(212)

      j1 = je - 2 * (ng + del) + 1
      il = ng / 2 + 1
      jl =  1 + del
      is = im / 2 + 1
      js =  1

#if NDIMS == 3
    case(213)

      j1 = je - 2 * (ng + del) + 1
      kl = ng / 2 + 1
      jl =  1 + del
      js =  1
      ks = km / 2 + 1

    case(214)

      j1 = je - 2 * (ng + del) + 1
      il = ng / 2 + 1
      jl =  1 + del
      kl = ng / 2 + 1
      is = im / 2 + 1
      js =  1
      ks = km / 2 + 1
#endif /* NDIMS == 3 */

    case(221)

      j1 = jb - 2 * del
      jl =  1 + del
      js = jeu

    case(222)

      j1 = jb - 2 * del
      il = ng / 2 + 1
      jl =  1 + del
      is = im / 2 + 1
      js = jeu

#if NDIMS == 3
    case(223)

      j1 = jb - 2 * del
      jl =  1 + del
      kl = ng / 2 + 1
      js = jeu
      ks = km / 2 + 1

    case(224)

      j1 = jb - 2 * del
      il = ng / 2 + 1
      jl =  1 + del
      kl = ng / 2 + 1
      is = im / 2 + 1
      js = jeu
      ks = km / 2 + 1


!! Z-direction
!!
    case(311)

      k1 = ke - 2 * (ng + del) + 1
      kl =  1 + del
      ks =  1

    case(312)

      k1 = ke - 2 * (ng + del) + 1
      il = ng / 2 + 1
      kl =  1 + del
      is = im / 2 + 1
      ks =  1

    case(313)

      k1 = ke - 2 * (ng + del) + 1
      jl = ng / 2 + 1
      kl =  1 + del
      js = jm / 2 + 1
      ks =  1

    case(314)

      k1 = ke - 2 * (ng + del) + 1
      il = ng / 2 + 1
      jl = ng / 2 + 1
      kl =  1 + del
      is = im / 2 + 1
      js = jm / 2 + 1
      ks =  1

    case(321)

      k1 = kb - 2 * del
      kl =  1 + del
      ks = keu

    case(322)

      k1 = kb - 2 * del
      il = ng / 2 + 1
      kl =  1 + del
      is = im / 2 + 1
      ks = keu

    case(323)

      k1 = kb - 2 * del
      jl = ng / 2 + 1
      kl =  1 + del
      js = jm / 2 + 1
      ks = keu

    case(324)

      k1 = kb - 2 * del
      il = ng / 2 + 1
      jl = ng / 2 + 1
      kl =  1 + del
      is = im / 2 + 1
      js = jm / 2 + 1
      ks = keu
#endif /* NDIMS == 3 */

    case default
      call print_warning("boundaries::bnd_rest", "Boundary flag unsupported!")

    end select

! prepare upper indices
!
      i2 = i1 + fm(1) - 1
      j2 = j1 + fm(2) - 1
      k2 = k1 + fm(3) - 1

      iu = il + pm(1) - 1
      ju = jl + pm(2) - 1
      ku = kl + pm(3) - 1

      it = is + pm(1) - 1
      jt = js + pm(2) - 1
#if NDIMS == 3
      kt = ks + pm(3) - 1
#endif /* NDIMS == 3 */
#if NDIMS == 2
      kt = 1
#endif /* NDIMS == 2 */

! allocate temporary array
!
    allocate(ux(cm(1),cm(2),cm(3)))

! iterate over all variables
!
    do q = 1, nfl

! shrink the boundary
!
      call shrink(fm, cm, 0, un(q,i1:i2,j1:j2,k1:k2), ux(:,:,:), 'm', 'm', 'm')

! copy shrinked boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = ux(il:iu,jl:ju,kl:ku)

    end do

#ifdef MHD
#ifdef FIELDCD
! iterate over magnetic field components
!
    do q = ibx, ibz

! shrink the boundary
!
      call shrink(fm, cm, 0, un(q,i1:i2,j1:j2,k1:k2), ux(:,:,:), 'm', 'm', 'm')

! copy shrinked boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = ux(il:iu,jl:ju,kl:ku)

    end do
#endif /* FIELDCD */
#endif /* MHD */

! deallocate temporary array
!
    deallocate(ux)
!
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

    use blocks, only : block_data, nvr, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                  &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                  &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : expand

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nqt,im,jm,km), intent(in)    :: un
    integer                             , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: ii, i, j, k, q
    integer :: i1, i2, j1, j2, k1, k2
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt

! local arrays
!
    integer, dimension(3)                  :: dm, cm, pm
    real   , dimension(:,:,:), allocatable :: ux

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
    dm(:)    = (/ im, jm, km /)
    dm(idir) = ng
    cm(:)    = dm(:) / 2 + 2 * del
    cm(idir) = ng / 2 + 2 * del
    pm(:)    = 2 * dm(:)
#if NDIMS == 2
    dm(3)    = 1
    cm(3)    = 1
    pm(3)    = 1
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

! allocate temporary array
!
    allocate(ux(pm(1),pm(2),pm(3)))

! iterate over all variables
!
    do q = 1, nfl

! expand the boundary
!
      call expand(cm, pm, 0, un(q,i1:i2,j1:j2,k1:k2), ux(:,:,:), 'a', 'a', 'a')

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = ux(il:iu,jl:ju,kl:ku)

    end do
#ifdef MHD
#ifdef FIELDCD
! iterate over magnetic field components
!
    do q = ibx, ibz

! expand the boundary
!
      call expand(cm, pm, 0, un(q,i1:i2,j1:j2,k1:k2), ux(:,:,:), 'a', 'a', 'a')

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = ux(il:iu,jl:ju,kl:ku)

    end do
#endif /* FIELDCD */
#endif /* MHD */

! deallocate temporary array
!
    deallocate(ux)
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
  subroutine bnd_spec(pb, id, il, ip)

    use blocks, only : block_data, nvr, imx, imy, imz
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , ng, im, jm, km, ib, ibl, ie, ieu, jb, jbl, je, jeu      &
                     , kb, kbl, ke, keu
    use error , only : print_warning

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
              pb%u(  :,it,j,k) =  pb%u(   :,is,j,k)
              pb%u(imx,it,j,k) = -pb%u(imx,is,j,k)
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
              pb%u(  :,it,j,k) =  pb%u(  :,is,j,k)
              pb%u(imx,it,j,k) = -pb%u(imx,is,j,k)
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
              pb%u(  :,i,jt,k) =  pb%u(  :,i,js,k)
              pb%u(imy,i,jt,k) = -pb%u(imy,i,js,k)
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
              pb%u(  :,i,jt,k) =  pb%u(  :,i,js,k)
              pb%u(imy,i,jt,k) = -pb%u(imy,i,js,k)
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
              pb%u(  :,i,j,kt) =  pb%u(  :,i,j,ks)
              pb%u(imz,i,j,kt) = -pb%u(imz,i,j,ks)
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
              pb%u(  :,i,j,kt) =  pb%u(  :,i,j,ks)
              pb%u(imz,i,j,kt) = -pb%u(imz,i,j,ks)
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
