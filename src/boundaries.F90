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

    use config  , only : im, jm, km, maxlev
    use blocks  , only : nvr, nqt, ndims, nsides, nfaces, nchild               &
                       , block_meta, block_info, pointer_info, list_meta
    use error   , only : print_error
    use mpitools, only : ncpus, ncpu, msendf, mrecvf

    implicit none

! local variables
!
    integer :: idir, iside, iface, nside, nface, dlev, l, p, q
    integer :: level

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
! iterate over all levels starting from the top one
!
    do level = maxlev, 1, -1

! iterate over all directions
!
      do idir = 1, ndims

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
          if (pblock%leaf .and. pblock%level .eq. level) then

! iterate over all neighbors
!
            do iside = 1, nsides
              do iface = 1, nfaces

! associate pointer to the neighbor
!
                pneigh => pblock%neigh(idir,iside,iface)%ptr

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
                      dlev = pblock%level - pneigh%level

! depending on the level difference
!
                      select case(dlev)
                      case(-1)  ! neighbor is on higer level -> restriction
                        call bnd_rest(pblock%data, pneigh%data%u, idir, iside, iface)
                      case(0)   ! block are on the same level -> copying
                        if (iface .eq. 1) &
                          call bnd_copy(pblock%data, pneigh%data%u, idir, iside, iface)
                      case(1)   ! neighbor is on lower level -> prolongation
                        if (iface .eq. 1) then
                          nside = 3 - iside
                          nface = 1
                          do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                            nface = nface + 1
                          end do

                          call bnd_prol(pblock%data, pneigh%data%u, idir, iside, nface)
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
                    pinfo%direction        =  idir
                    pinfo%side             =  iside
                    pinfo%face             =  iface
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
                  if (iface .eq. 1 .and. pblock%cpu .eq. ncpu) &
                    call bnd_spec(pblock%data, idir, iside, iface)
#else /* MPI */
                  if (iface .eq. 1) &
                    call bnd_spec(pblock%data, idir, iside, iface)
#endif /* MPI */
                end if

              end do ! iteration over faces
            end do ! iteration over sides

          end if ! if leaf and at the current level

! assign pointer to the next block
!
          pblock => pblock%next

        end do ! iteration over all blocks

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
                  iside = pinfo%side
                  iface = pinfo%face
                  dlev  = pinfo%level_difference

! update boundaries
!
                  select case(dlev)
                  case(-1)  ! neighbor is on higer level -> restriction
                    call bnd_rest(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside, iface)
                  case(0)   ! block are on the same level -> copying
                    if (iface .eq. 1) &
                      call bnd_copy(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside, iface)
                  case(1)   ! neighbor is on lower level -> prolongation
                    if (iface .eq. 1) then
                      pblock => pinfo%block
                      pneigh => pblock%neigh(idir,iside,iface)%ptr

                      nside = 3 - iside
                      nface = 1
                      do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                        nface = nface + 1
                      end do

                      call bnd_prol(pinfo%block%data, rbuf(l,:,:,:,:), idir, iside, nface)
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

            end if ! if block counter > 0

          end do ! iteration over q
        end do ! iteration over p
#endif /* MPI */
      end do ! iteration over directions
    end do ! iteration over levels
!
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
  subroutine bnd_copy(pdata, u, idir, iside, iface)

    use blocks, only : block_data, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : im, ib, ibl, ibu, ie, iel, ieu                          &
                     , jm, jb, jbl, jbu, je, jel, jeu                          &
                     , km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_error

    implicit none

! arguments
!
    type(block_data), pointer       , intent(inout) :: pdata

    real   , dimension(nqt,im,jm,km), intent(in)    :: u
    integer                         , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt
!
!-------------------------------------------------------------------------------
!
! prepare common indices
!
    il = 1
    iu = im
    jl = 1
    ju = jm
    kl = 1
    ku = km

    is = 1
    it = im
    js = 1
    jt = jm
    ks = 1
    kt = km

! prepare source and destination boundary indices
!
    select case(idir)

    case(1)

      if (iside .eq. 1) then
        il = iel
        iu = ie

        is = 1
        it = ibl
      else
        il = ib
        iu = ibu

        is = ieu
        it = im
      end if

    case(2)

      if (iside .eq. 1) then
        jl = jel
        ju = je

        js = 1
        jt = jbl
      else
        jl = jb
        ju = jbu

        js = jeu
        jt = jm
      end if

#if NDIMS == 3
    case(3)

      if (iside .eq. 1) then
        kl = kel
        ku = ke

        ks = 1
        kt = kbl
      else
        kl = kb
        ku = kbu

        ks = keu
        kt = km
      end if
#endif /* NDIMS == 3 */

    case default
      call print_error("boundaries::bnd_copy", "Direction unsupported!")

    end select

! perform update of the fluid variables
!
    pdata%u(1:nfl,is:it,js:jt,ks:kt) = u(1:nfl,il:iu,jl:ju,kl:ku)

#ifdef MHD
#ifdef FIELDCD
! perform update of the cell-centered magnetic field components
!
    pdata%u(ibx:ibz,is:it,js:jt,ks:kt) = u(ibx:ibz,il:iu,jl:ju,kl:ku)
#endif /* FIELDCD */
#ifdef FLUXCT
! perform update of the staggered magnetic field components
!
    if (it .eq. ibl) then
      pdata%u(ibx,is:it-1,js:jt,ks:kt) = u(ibx,il:iu-1,jl:ju,kl:ku)
    else
      pdata%u(ibx,is:it  ,js:jt,ks:kt) = u(ibx,il:iu  ,jl:ju,kl:ku)
    end if

    if (jt .eq. jbl) then
      pdata%u(iby,is:it,js:jt-1,ks:kt) = u(iby,il:iu,jl:ju-1,kl:ku)
    else
      pdata%u(iby,is:it,js:jt  ,ks:kt) = u(iby,il:iu,jl:ju  ,kl:ku)
    end if

#if NDIMS == 3
    if (kt .eq. kbl) then
      pdata%u(ibz,is:it,js:jt,ks:kt-1) = u(ibz,il:iu,jl:ju,kl:ku-1)
    else
      pdata%u(ibz,is:it,js:jt,ks:kt  ) = u(ibz,il:iu,jl:ju,kl:ku  )
    end if
#else /* NDIMS == 3 */
    pdata%u(ibz,is:it,js:jt,ks:kt) = u(ibz,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */
#endif /* FLUXCT */
#endif /* MHD */
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
  subroutine bnd_rest(pdata, ub, idir, iside, iface)

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
    real        , dimension(nqt,im,jm,km), intent(in)    :: ub
    integer                              , intent(in)    :: idir, iside, iface

! local variables
!
    integer ::  q,  i,  j,  k
    integer :: il, iu, jl, ju, kl, ku
    integer :: i1, i2, j1, j2, k1, k2
    integer :: is, it, js, jt, ks, kt
    integer :: if, jf, kf

! local arrays
!
    integer, dimension(3)                  :: dm, cm, fm
    real   , dimension(:,:,:), allocatable :: u

! parameters
!
    integer :: del = 2
!
!-------------------------------------------------------------------------------
!
! prepare dimensions
!
    dm(:)    = (/ im, jm, km /)
    dm(idir) = ng
    cm(:)    = dm(:) / 2
    cm(idir) = ng + del
    fm(:)    = 2 * cm(:)
#if NDIMS == 2
    dm(3)  = 1
    cm(3)  = 1
    fm(3)  = 1
#endif /* NDIMS == 2 */

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
        i1 = ie - 2 * ng - del + 1
        i2 = ie + del
      else
        i1 = ib - del
        i2 = ib + 2 * ng + del - 1
      end if
      j1 =  1
      j2 = jm
      k1 =  1
      k2 = km

! indices for the output array
!
      il = del / 2 + 1
      iu = il + ng - 1
      jl =   1           + ng / 2 * jf
      ju = (jm - ng) / 2 + ng / 2 * jf
#if NDIMS == 3
      kl =   1           + ng / 2 * kf
      ku = (km - ng) / 2 + ng / 2 * kf
#else /* NDIMS == 3 */
      kl = 1
      ku = 1
#endif /* NDIMS == 3 */

! indices of the destination array
!
      if (if .eq. 0) then
        is = 1
        it = ng
      else
        is = ie + 1
        it = im
      end if
      if (jf .eq. 0) then
        js = ng / 2 + 1
        jt = jm / 2
      else
        js = jm / 2 + 1
        jt = jm - ng / 2
      end if
#if NDIMS == 3
      if (kf .eq. 0) then
        ks = ng / 2 + 1
        kt = km / 2
      else
        ks = km / 2 + 1
        kt = km - ng / 2
      end if
#else /* NDIMS == 3 */
      ks = 1
      kt = 1
#endif /* NDIMS == 3 */

    case(2)

      if = mod(iface - 1, 2)
      jf =     iside - 1
      kf =    (iface - 1) /   2

! indices of the input array
!
      i1 = 1
      i2 = im
      if (jf .eq. 0) then
        j1 = je - 2 * ng - del + 1
        j2 = je + del
      else
        j1 = jb - del
        j2 = jb + 2 * ng + del - 1
      end if
      k1 = 1
      k2 = km

! indices for the output array
!
      il =   1           + ng / 2 * if
      iu = (im - ng) / 2 + ng / 2 * if
      jl = del / 2 + 1
      ju = jl + ng - 1
#if NDIMS == 3
      kl =   1           + ng / 2 * kf
      ku = (km - ng) / 2 + ng / 2 * kf
#else /* NDIMS == 3 */
      kl = 1
      ku = 1
#endif /* NDIMS == 3 */

! indices of the destination array
!
      if (if .eq. 0) then
        is = ng / 2 + 1
        it = im / 2
      else
        is = im / 2 + 1
        it = im - ng / 2
      end if
      if (jf .eq. 0) then
        js = 1
        jt = ng
      else
        js = je + 1
        jt = jm
      end if
#if NDIMS == 3
      if (kf .eq. 0) then
        ks = ng / 2 + 1
        kt = km / 2
      else
        ks = km / 2 + 1
        kt = km - ng / 2
      end if
#else /* NDIMS == 3 */
      ks = 1
      kt = 1
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

      if = mod(iface - 1, 2)
      jf =    (iface - 1) / 2
      kf =     iside - 1

! indices of the input array
!
      i1 = 1
      i2 = im
      j1 = 1
      j2 = jm
      if (kf .eq. 0) then
        k1 = ke - 2 * ng - del + 1
        k2 = ke + del
      else
        k1 = kb - del
        k2 = kb + 2 * ng + del - 1
      end if

! indices for the output array
!
      il =   1           + ng / 2 * if
      iu = (im - ng) / 2 + ng / 2 * if
      jl =   1           + ng / 2 * jf
      ju = (jm - ng) / 2 + ng / 2 * jf
      kl = del / 2 + 1
      ku = kl + ng - 1

! indices of the destination array
!
      if (if .eq. 0) then
        is = ng / 2 + 1
        it = im / 2
      else
        is = im / 2 + 1
        it = im - ng / 2
      end if
      if (jf .eq. 0) then
        js = ng / 2 + 1
        jt = jm / 2
      else
        js = jm / 2 + 1
        jt = jm - ng / 2
      end if
      if (kf .eq. 0) then
        ks = 1
        kt = ng
      else
        ks = ke + 1
        kt = km
      end if
#endif /* NDIMS == 3 */

    end select

! allocate temporary array
!
    allocate(u(cm(1),cm(2),cm(3)))

! iterate over all variables
!
    do q = 1, nfl

! shrink the boundary
!
      call shrink(fm, cm, ub(q,i1:i2,j1:j2,k1:k2), u(:,:,:), 'm', 'm', 'm')

! copy shrinked boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(il:iu,jl:ju,kl:ku)

    end do

#ifdef MHD
#ifdef FIELDCD
! iterate over magnetic field components
!
    do q = ibx, ibz

! shrink the boundary
!
      call shrink(fm, cm, ub(q,i1:i2,j1:j2,k1:k2), u(:,:,:), 'm', 'm', 'm')

! copy shrinked boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(il:iu,jl:ju,kl:ku)

    end do
#endif /* FIELDCD */
#ifdef FLUXCT
! X-component of magnetic field
!
    call shrink(fm, cm, ub(ibx,i1:i2,j1:j2,k1:k2), u(:,:,:), 'c', 'm', 'm')

! copy shrinked boundary in the proper place of the block
!
    if (it .eq. ng) then
      pdata%u(ibx,is:it-1,js:jt,ks:kt) = u(il:iu-1,jl:ju,kl:ku)
    else
      pdata%u(ibx,is:it  ,js:jt,ks:kt) = u(il:iu  ,jl:ju,kl:ku)
    end if

! Y-component of magnetic field
!
    call shrink(fm, cm, ub(iby,i1:i2,j1:j2,k1:k2), u(:,:,:), 'm', 'c', 'm')

! copy shrinked boundary in the proper place of the block
!
    if (jt .eq. ng) then
      pdata%u(iby,is:it,js:jt-1,ks:kt) = u(il:iu,jl:ju-1,kl:ku)
    else
      pdata%u(iby,is:it,js:jt  ,ks:kt) = u(il:iu,jl:ju  ,kl:ku)
    end if

! Z-component of magnetic field
!
    call shrink(fm, cm, ub(ibz,i1:i2,j1:j2,k1:k2), u(:,:,:), 'm', 'm', 'c')

! copy shrinked boundary in the proper place of the block
!
#if NDIMS == 3
    if (kt .eq. ng) then
      pdata%u(ibz,is:it,js:jt,ks:kt-1) = u(il:iu,jl:ju,kl:ku-1)
    else
      pdata%u(ibz,is:it,js:jt,ks:kt  ) = u(il:iu,jl:ju,kl:ku  )
    end if
#else /* NDIMS == 3 */
    pdata%u(ibz,is:it,js:jt,ks:kt) = u(il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */
#endif /* FLUXCT */
#endif /* MHD */

! deallocate temporary array
!
    deallocate(u)
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

    use blocks, only : block_data, nvr, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                  &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                  &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : expand
#if defined MHD && defined FLUXCT
    use interpolation, only : expand_mag
#endif /* MHD & FLUXCT */

    implicit none

! arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real       , dimension(nqt,im,jm,km), intent(in)    :: ub
    integer                             , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: q, i, j, k
    integer :: i1, i2, j1, j2, k1, k2
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt
    integer :: if, jf, kf

! local arrays
!
    integer, dimension(3)                    :: dm, cm, pm
    real   , dimension(:,:,:)  , allocatable :: u
    real   , dimension(:,:,:,:), allocatable :: b

! parameters
!
    integer :: del = 1
!
!-------------------------------------------------------------------------------
!
! prepare dimensions
!
    dm(:)    = (/ im, jm, km /)
    dm(idir) = ng
    cm(:)    = dm(:) / 2 + ng
    cm(idir) = 3 * ng / 2 + 2 * del
    pm(:)    = dm(:)
    pm(idir) = ng + 4 * del
#if NDIMS == 2
    dm(3)    = 1
    cm(3)    = 1
    pm(3)    = 1
#endif /* NDIMS == 2 */

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
        i1 = ie - ng + 1 - del
        i2 = ie + ng / 2 + del
      else
        i1 = ib - ng / 2 - del
        i2 = ib + ng - 1 + del
      end if
      if (jf .eq. 0) then
        j1 =  1
        j2 = jm / 2 + ng
      else
        j1 = jm / 2 - ng + 1
        j2 = jm
      end if
#if NDIMS == 3
      if (kf .eq. 0) then
        k1 =  1
        k2 = km / 2 + ng
      else
        k1 = km / 2 - ng + 1
        k2 = km
      end if
#else /* NDIMS == 3 */
      k1 = 1
      k2 = 1
#endif /* NDIMS == 3 */

! indices for the output array
!
      il =  1 + 2 * del
      iu = ng + 2 * del
      jl =  1
      ju = jm
      kl =  1
      ku = km

! indices of the destination array
!
      if (if .eq. 0) then
        is =  1
        it = ng
      else
        is = ie + 1
        it = im
      end if
      js = 1
      jt = jm
      ks = 1
      kt = km

    case(2)

      if = mod(iface - 1, 2)
      jf =     iside - 1
      kf =    (iface - 1) / 2

! indices of the input array
!
      if (if .eq. 0) then
        i1 =  1
        i2 = im / 2 + ng
      else
        i1 = im / 2 - ng + 1
        i2 = im
      end if
      if (jf .eq. 0) then
        j1 = je - ng + 1 - del
        j2 = je + ng / 2 + del
      else
        j1 = jb - ng / 2 - del
        j2 = jb + ng - 1 + del
      end if
#if NDIMS == 3
      if (kf .eq. 0) then
        k1 =  1
        k2 = km / 2 + ng
      else
        k1 = km / 2 - ng + 1
        k2 = km
      end if
#else /* NDIMS == 3 */
      k1 = 1
      k2 = 1
#endif /* NDIMS == 3 */

! indices for the output array
!
      il =  1
      iu = im
      jl =  1 + 2 * del
      ju = ng + 2 * del
      kl =  1
      ku = km

! indices of the destination array
!
      is =  1
      it = im
      if (jf .eq. 0) then
        js =  1
        jt = ng
      else
        js = je + 1
        jt = jm
      end if
      ks =  1
      kt = km

#if NDIMS == 3
    case(3)

      if = mod(iface - 1, 2)
      jf =    (iface - 1) / 2
      kf =     iside - 1

! indices of the input array
!
      if (if .eq. 0) then
        i1 =  1
        i2 = im / 2 + ng
      else
        i1 = im / 2 - ng + 1
        i2 = im
      end if
      if (jf .eq. 0) then
        j1 =  1
        j2 = jm / 2 + ng
      else
        j1 = jm / 2 - ng + 1
        j2 = jm
      end if
      if (kf .eq. 0) then
        k1 = ke - ng + 1 - del
        k2 = ke + ng / 2 + del
      else
        k1 = kb - ng / 2 - del
        k2 = kb + ng - 1 + del
      end if

! indices for the output array
!
      il =  1
      iu = im
      jl =  1
      ju = jm
      kl =  1 + 2 * del
      ku = ng + 2 * del

! indices of the destination array
!
      is =  1
      it = im
      js =  1
      jt = jm
      if (kf .eq. 0) then
        ks =  1
        kt = ng
      else
        ks = ke + 1
        kt = km
      end if
#endif /* NDIMS == 3 */

    end select

! allocate temporary array
!
    allocate(u(pm(1),pm(2),pm(3)))

! iterate over all variables
!
    do q = 1, nfl

! expand the boundary
!
      call expand(cm, pm, ng, ub(q,i1:i2,j1:j2,k1:k2), u(:,:,:), 't', 't', 't')

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(il:iu,jl:ju,kl:ku)

    end do
#ifdef MHD
#ifdef FIELDCD
! iterate over magnetic field components
!
    do q = ibx, ibz

! expand the boundary
!
      call expand(cm, pm, ng, ub(q,i1:i2,j1:j2,k1:k2), u(:,:,:), 't', 't', 't')

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(il:iu,jl:ju,kl:ku)

    end do
#endif /* FIELDCD */
#ifdef FLUXCT
! allocate space for the prolongated magnetic field components
!
    allocate(b(3,pm(1),pm(2),pm(3)))

! prolongate magnetic field preserving the divergence-free condition
!
    call expand_mag(cm, pm, ng, ub(ibx,i1:i2,j1:j2,k1:k2), ub(iby,i1:i2,j1:j2,k1:k2), ub(ibz,i1:i2,j1:j2,k1:k2), b(1,:,:,:), b(2,:,:,:), b(3,:,:,:))

! update the X component of magnetic field
!
    if (it .eq. ng) then
      pdata%u(ibx,is:it-1,js:jt,ks:kt) = b(1,il:iu-1,jl:ju,kl:ku)
    else
      pdata%u(ibx,is:it  ,js:jt,ks:kt) = b(1,il:iu  ,jl:ju,kl:ku)
    end if

! update the Y component of magnetic field
!
    if (jt .eq. ng) then
      pdata%u(iby,is:it,js:jt-1,ks:kt) = b(2,il:iu,jl:ju-1,kl:ku)
    else
      pdata%u(iby,is:it,js:jt  ,ks:kt) = b(2,il:iu,jl:ju  ,kl:ku)
    end if

! update the Z component of magnetic field
!
#if NDIMS == 3
    if (kt .eq. ng) then
      pdata%u(ibz,is:it,js:jt,ks:kt-1) = b(3,il:iu,jl:ju,kl:ku-1)
    else
      pdata%u(ibz,is:it,js:jt,ks:kt  ) = b(3,il:iu,jl:ju,kl:ku  )
    end if
#else /* NDIMS == 3 */
    pdata%u(ibz,is:it,js:jt,ks:kt) = b(3,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

! deallocate space for the prolongated magnetic field components
!
    deallocate(b)
#endif /* FLUXCT */
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
