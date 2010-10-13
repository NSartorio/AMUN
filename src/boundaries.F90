!!******************************************************************************
!!
!! module: boundaries - routines for handling the boundary conditions
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module boundaries

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! boundary_flux: subroutine sweeps over all leaf blocks and updates the flux
!                boundaries in order to keep consistency between different
!                levels
!
! notes: - only update the fluxes in blocks from neighbors at higher level
!
!! TODO:
!!
!! - MPI version
!!
!
!===============================================================================
!
  subroutine boundary_fluxes

    use blocks, only : block_meta, list_meta
    use blocks, only : nsides, nfaces

    implicit none

! local variables
!
    integer :: idir, iside, iface

! local pointers
!
    type(block_meta), pointer :: pblock, pneigh
!
!-------------------------------------------------------------------------------
!
! iterate over all directions
!
    do idir = 1, NDIMS

! iterate over all leaf blocks and perform the boundary update
!
      pblock => list_meta
      do while (associated(pblock))

! check if the current block is a leaf
!
        if (pblock%leaf) then

! iterate over all neighbors (directions, sides, and faces)
!
          do iside = 1, nsides
            do iface = 1, nfaces

! associate the pointer pneigh to the current neighbor
!
              pneigh => pblock%neigh(idir,iside,iface)%ptr

! check if the neighbor is associated
!
              if (associated(pneigh)) then

! if the neighbor is at the higher level update the fluxes of the current block
!
                if (pneigh%level .gt. pblock%level) &
                  call restrict_flux(pneigh%data%f, pneigh%data%e, pblock%data, idir, iside, iface)

              end if

            end do ! iteration over faces
          end do ! iteration over sides

        end if ! condition: if the leaf

! assign the pointer to the next block
!
        pblock => pblock%next

      end do ! iteration over all blocks

    end do ! iteration over directions
!
!-------------------------------------------------------------------------------
!
  end subroutine boundary_fluxes
!
!===============================================================================
!
! boundary: subroutine sweeps over all leaf blocks and performs their boundary
!           update
!
!===============================================================================
!
  subroutine boundary_variables

! references to other modules
!
    use blocks, only : block_meta, block_data, list_meta
    use blocks, only : ndims, nsides, nfaces

! declare variables
!
    implicit none

! local variables
!
    integer :: idir, iside, iface, nside, nface

! local pointers
!
    type(block_meta), pointer :: pblock, pneigh
    type(block_data), pointer :: pbdata, pndata
!
!-------------------------------------------------------------------------------
!
! update boundaries which don't have neighbors
!
    do idir = 1, ndims
      pblock => list_meta
      do while(associated(pblock))
        if (pblock%leaf) then
          do iside = 1, nsides
            do iface = 1, nfaces
              pneigh => pblock%neigh(idir,iside,iface)%ptr
              if (.not. associated(pneigh)) then
                if (iface .eq. 1) &
                  call bnd_spec(pblock%data, idir, iside, iface)
              end if
            end do ! faces
          end do ! sides
        end if ! leaf and current level
        pblock => pblock%next ! assign pointer to the next block
      end do ! metablocks
    end do ! directions

! update boundaries from the blocks at the same level
!
    do idir = 1, ndims
      pblock => list_meta
      do while(associated(pblock))
        if (pblock%leaf) then
          do iside = 1, nsides
            do iface = 1, nfaces
              pneigh => pblock%neigh(idir,iside,iface)%ptr
              if (associated(pneigh)) then
                pbdata => pblock%data
                pndata => pneigh%data

                if (pblock%level .eq. pneigh%level .and. iface .eq. 1) &
                  call bnd_copy(pbdata, pndata%u, idir, iside)
              end if
            end do ! faces
          end do ! sides
        end if ! leaf and current level
        pblock => pblock%next ! assign pointer to the next block
      end do ! metablocks
    end do ! directions

! TODO: update corners from the neighbors at higher level
!

! update boundaries from all neighbors
!
    do idir = 1, ndims
      pblock => list_meta
      do while(associated(pblock))
        if (pblock%leaf) then
          do iside = 1, nsides
            do iface = 1, nfaces
              pneigh => pblock%neigh(idir,iside,iface)%ptr
              if (associated(pneigh)) then
                pbdata => pblock%data
                pndata => pneigh%data

                if (pblock%level .lt. pneigh%level) &
                  call bnd_rest(pblock%data, pneigh%data%u, idir, iside, iface)
                if (pblock%level .eq. pneigh%level .and. iface .eq. 1) &
                  call bnd_copy(pbdata, pndata%u, idir, iside)
                if (pblock%level .gt. pneigh%level) then
                  if (iface .eq. 1) then
                    nside = 3 - iside
                    nface = 1
                    do while(pblock%id .ne. pneigh%neigh(idir,nside,nface)%ptr%id)
                      nface = nface + 1
                    end do

                    pbdata => pblock%data
                    pndata => pneigh%data
                    call bnd_prol(pbdata, pndata%u, idir, iside, nface)
                  end if
                end if
              end if
            end do ! faces
          end do ! sides
        end if ! leaf and current level
        pblock => pblock%next ! assign pointer to the next block
      end do ! metablocks
    end do ! directions
!
!-------------------------------------------------------------------------------
!
  end subroutine boundary_variables
!
!===============================================================================
!
! restrict_flux: subroutine restricts the flux from the neighbor which is at
!                higher level and updates the boundary flux of the current block
!
!===============================================================================
!
  subroutine restrict_flux(f, e, pdata, idir, iside, iface)

    use config, only : ng, nd, nh, im, jm, km, ih, jh, kh, ib, jb, kb  &
                     , ie, je, ke, ibl, jbl, kbl, ieu, jeu, keu
    use blocks, only : block_data, nfl

    implicit none

! arguments
!
    type(block_data), pointer             , intent(inout) :: pdata
    real   , dimension(NDIMS,nfl,im,jm,km), intent(in)    :: f
    real   , dimension(        3,im,jm,km), intent(in)    :: e
    integer                               , intent(in)    :: idir, iside, iface

! local variables
!
    integer :: if, jf, kf
    integer :: il, jl, kl, ip, jp, kp, iu, ju, ku
    integer :: is, js, ks, it, jt, kt
!
!-------------------------------------------------------------------------------
!
#if NDIMS == 2
    kl = 1
    ku = 1
    ks = 1
#endif /* NDIMS == 2 */

    select case(idir)

    case(1)

! calculate the indices determining the position of the neighbor
!
      if =     iside - 1
      jf = mod(iface - 1, 2)
#if NDIMS == 3
      kf =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! indices of the input and destination arrays
!
      if (if .eq. 0) then
        il = ieu - nd
        iu = ie
        is = 1
        it = ibl
      else
        il = ib
        iu = ibl + nd
        is = ieu
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

! update the boundaries of the fluid fluxes
!
#if NDIMS == 2
!       if (it .eq. im) then
!         pdata%f(1,1:nfl,is-1:it,js:jt,1) = 0.5 * (f(1,1:nfl,ip-2:iu:2,jl:ju:2,1)   &
!                                                 + f(1,1:nfl,ip-2:iu:2,jp:ju:2,1))
!       else
!         pdata%f(1,1:nfl,is  :it,js:jt,1) = 0.5 * (f(1,1:nfl,ip  :iu:2,jl:ju:2,1)   &
!                                                 + f(1,1:nfl,ip  :iu:2,jp:ju:2,1))
!       end if
!       pdata%f(2,1:nfl,is:it,js:jt,1) = 0.5 * (f(2,1:nfl,il:iu:2,jp:ju:2,1)     &
!                                             + f(2,1:nfl,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
!       pdata%f(1,1:nfl,is:it,js:jt,ks:kt) =                                     &
!                                 0.25 * (f(1,1:nfl,ip:iu:2,jl:ju:2,kl:ku:2)     &
!                                       + f(1,1:nfl,ip:iu:2,jp:ju:2,kl:ku:2)     &
!                                       + f(1,1:nfl,ip:iu:2,jl:ju:2,kp:ku:2)     &
!                                       + f(1,1:nfl,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#ifdef MHD

! update the boundaries of the electromotive force
!
#if NDIMS == 2
!       pdata%e(1,is:it,js:jt,1) = 0.5 * (e(1,il:iu:2,jp:ju:2,1)                 &
!                                       + e(1,ip:iu:2,jp:ju:2,1))
!       if (it .eq. im) then
!         pdata%e(2,is-1:it,js:jt,1) = 0.5 * (e(2,ip-2:iu:2,jl:ju:2,1)           &
!                                           + e(2,ip-2:iu:2,jp:ju:2,1))
!         pdata%e(3,is-1:it,js:jt,1) =        e(3,ip-2:iu:2,jp:ju:2,1)
!       else
!         pdata%e(2,is  :it,js:jt,1) = 0.5 * (e(2,ip  :iu:2,jl:ju:2,1)           &
!                                           + e(2,ip  :iu:2,jp:ju:2,1))
!         pdata%e(3,is  :it,js:jt,1) =        e(3,ip  :iu:2,jp:ju:2,1)
!       end if
#endif /* NDIMS == 2 */
#if NDIMS == 3
!       pdata%e(1,is:it,js:jt,ks:kt) = 0.5 * (e(1,il:iu:2,jp:ju:2,kp:ku:2)       &
!                                           + e(1,ip:iu:2,jp:ju:2,kp:ku:2))
!       pdata%e(2,is:it,js:jt,ks:kt) = 0.5 * (e(2,ip:iu:2,jl:ju:2,kp:ku:2)       &
!                                           + e(2,ip:iu:2,jp:ju:2,kp:ku:2))
!       pdata%e(3,is:it,js:jt,ks:kt) = 0.5 * (e(3,ip:iu:2,jp:ju:2,kl:ku:2)       &
!                                           + e(3,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#endif /* MHD */

    case(2)

! calculate the indices determining the position of the neighbor
!
      if = mod(iface - 1, 2)
      jf =     iside - 1
#if NDIMS == 3
      kf =    (iface - 1) / 2
#endif /* NDIMS == 3 */

! indices of the input and destination arrays
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
        jl = jeu - nd
        ju = je
        js = 1
        jt = jbl
      else
        jl = jb
        ju = jbl + nd
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

! update the boundaries of the fluid fluxes
!
#if NDIMS == 2
!       pdata%f(1,1:nfl,is:it,js:jt,1) = 0.5 * (f(1,1:nfl,ip:iu:2,jl:ju:2,1)     &
!                                             + f(1,1:nfl,ip:iu:2,jp:ju:2,1))
!       if (jt .eq. jm) then
!         pdata%f(2,1:nfl,is:it,js-1:jt,1) = 0.5 * (f(2,1:nfl,il:iu:2,jp-2:ju:2,1)     &
!                                                 + f(2,1:nfl,ip:iu:2,jp-2:ju:2,1))
!       else
!         pdata%f(2,1:nfl,is:it,js  :jt,1) = 0.5 * (f(2,1:nfl,il:iu:2,jp  :ju:2,1)     &
!                                                 + f(2,1:nfl,ip:iu:2,jp  :ju:2,1))
!       end if
#endif /* NDIMS == 2 */
#if NDIMS == 3
!       pdata%f(2,1:nfl,is:it,js:jt,ks:kt) =                                     &
!                                 0.25 * (f(2,1:nfl,il:iu:2,jp:ju:2,kl:ku:2)     &
!                                       + f(2,1:nfl,ip:iu:2,jp:ju:2,kl:ku:2)     &
!                                       + f(2,1:nfl,il:iu:2,jp:ju:2,kp:ku:2)     &
!                                       + f(2,1:nfl,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#ifdef MHD

! update the boundaries of the electromotive force
!
#if NDIMS == 2
!       if (jt .eq. jm) then
!         pdata%e(1,is:it,js-1:jt,1) = 0.5 * (e(1,il:iu:2,jp-2:ju:2,1)           &
!                                           + e(1,ip:iu:2,jp-2:ju:2,1))
!         pdata%e(2,is:it,js  :jt,1) = 0.5 * (e(2,ip:iu:2,jl  :ju:2,1)           &
!                                           + e(2,ip:iu:2,jp  :ju:2,1))
!         pdata%e(3,is:it,js-1:jt,1) =        e(3,ip:iu:2,jp-2:ju:2,1)
!       else
!         pdata%e(1,is:it,js  :jt,1) = 0.5 * (e(1,il:iu:2,jp  :ju:2,1)           &
!                                           + e(1,ip:iu:2,jp  :ju:2,1))
!         pdata%e(3,is:it,js  :jt,1) =        e(3,ip:iu:2,jp  :ju:2,1)
!       end if
#endif /* NDIMS == 2 */
#if NDIMS == 3
!       pdata%e(1,is:it,js:jt,ks:kt) = 0.5 * (e(1,il:iu:2,jp:ju:2,kp:ku:2)       &
!                                           + e(1,ip:iu:2,jp:ju:2,kp:ku:2))
!       pdata%e(2,is:it,js:jt,ks:kt) = 0.5 * (e(2,ip:iu:2,jl:ju:2,kp:ku:2)       &
!                                           + e(2,ip:iu:2,jp:ju:2,kp:ku:2))
!       pdata%e(3,is:it,js:jt,ks:kt) = 0.5 * (e(3,ip:iu:2,jp:ju:2,kl:ku:2)       &
!                                           + e(3,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#endif /* MHD */

#if NDIMS == 3
    case(3)

! calculate the indices determining the position of the neighbor
!
      if = mod(iface - 1, 2)
      jf =    (iface - 1) / 2
      kf =     iside - 1

! indices of the input and destination arrays
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
        kl = keu - nd
        ku = ke
        ks = 1
        kt = kbl
      else
        kl = kb - 2
        ku = kbl + nd
        ks = ke
        kt = km
      end if
      kp = kl + 1

! update the boundaries of the fluid fluxes
!
      pdata%f(3,1:nfl,is:it,js:jt,ks:kt) =                                     &
                                0.25 * (f(3,1:nfl,il:iu:2,jl:ju:2,kp:ku:2)     &
                                      + f(3,1:nfl,ip:iu:2,jl:ju:2,kp:ku:2)     &
                                      + f(3,1:nfl,il:iu:2,jp:ju:2,kp:ku:2)     &
                                      + f(3,1:nfl,ip:iu:2,jp:ju:2,kp:ku:2))
#ifdef MHD

! update the boundaries of the electromotive force
!
      pdata%e(1,is:it,js:jt,ks:kt) = 0.5 * (e(1,il:iu:2,jp:ju:2,kp:ku:2)       &
                                          + e(1,ip:iu:2,jp:ju:2,kp:ku:2))
      pdata%e(2,is:it,js:jt,ks:kt) = 0.5 * (e(2,ip:iu:2,jl:ju:2,kp:ku:2)       &
                                          + e(2,ip:iu:2,jp:ju:2,kp:ku:2))
      pdata%e(3,is:it,js:jt,ks:kt) = 0.5 * (e(3,ip:iu:2,jp:ju:2,kl:ku:2)       &
                                          + e(3,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* MHD */
#endif /* NDIMS == 3 */

    end select
!
!-------------------------------------------------------------------------------
!
  end subroutine restrict_flux
!
!===============================================================================
!
! bnd_copy: subroutine copies the interior of neighbor to update the
!           boundaries of current block
!
!===============================================================================
!
  subroutine bnd_copy(pdata, u, idir, iside)

    use blocks, only : block_data, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : im, ib, ibl, ibu, ie, iel, ieu                          &
                     , jm, jb, jbl, jbu, je, jel, jeu                          &
                     , km, kb, kbl, kbu, ke, kel, keu

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
        pdata%u(1:nfl,1:ibl  ,1:jm,1:km) = u(1:nfl,iel:ie  ,1:jm,1:km)
#ifdef MHD
#ifdef FLUXCT
        pdata%u(  ibx,1:ibl-1,1:jm,1:km) = u(  ibx,iel:ie-1,1:jm,1:km)
        pdata%u(  iby,1:ibl  ,1:jm,1:km) = u(  iby,iel:ie  ,1:jm,1:km)
        pdata%u(  ibz,1:ibl  ,1:jm,1:km) = u(  ibz,iel:ie  ,1:jm,1:km)
#endif /* FLUXCT */
#endif /* MHD */
      else
        pdata%u(1:nfl,ieu:im ,1:jm,1:km) = u(1:nfl,ib:ibu  ,1:jm,1:km)
#ifdef MHD
#ifdef FLUXCT
        pdata%u(  ibx,ieu:im ,1:jm,1:km) = u(  ibx,ib:ibu  ,1:jm,1:km)
        pdata%u(  iby,ieu:im ,1:jm,1:km) = u(  iby,ib:ibu  ,1:jm,1:km)
        pdata%u(  ibz,ieu:im ,1:jm,1:km) = u(  ibz,ib:ibu  ,1:jm,1:km)
#endif /* FLUXCT */
#endif /* MHD */
      end if

    case(2)

      if (iside .eq. 1) then
        pdata%u(1:nfl,1:im,1:jbl  ,1:km) = u(1:nfl,1:im,jel:je  ,1:km)
#ifdef MHD
#ifdef FLUXCT
        pdata%u(  ibx,1:im,1:jbl  ,1:km) = u(  ibx,1:im,jel:je  ,1:km)
        pdata%u(  iby,1:im,1:jbl-1,1:km) = u(  iby,1:im,jel:je-1,1:km)
        pdata%u(  ibz,1:im,1:jbl  ,1:km) = u(  ibz,1:im,jel:je  ,1:km)
#endif /* FLUXCT */
#endif /* MHD */
      else
        pdata%u(1:nfl,1:im,jeu:jm ,1:km) = u(1:nfl,1:im,jb:jbu  ,1:km)
#ifdef MHD
#ifdef FLUXCT
        pdata%u(  ibx,1:im,jeu:jm ,1:km) = u(  ibx,1:im,jb:jbu  ,1:km)
        pdata%u(  iby,1:im,jeu:jm ,1:km) = u(  iby,1:im,jb:jbu  ,1:km)
        pdata%u(  ibz,1:im,jeu:jm ,1:km) = u(  ibz,1:im,jb:jbu  ,1:km)
#endif /* FLUXCT */
#endif /* MHD */
      end if

#if NDIMS == 3
    case(3)

      if (iside .eq. 1) then
        kl = kel
        ku = ke

        ks = 1
        kt = kbl

        pdata%u(1:nfl,1:im,1:jm,1:kbl  ) = u(1:nfl,1:im,jel:je  ,1:km)
#ifdef MHD
#ifdef FLUXCT
        pdata%u(  ibx,1:im,1:jm,1:kbl  ) = u(  ibx,1:im,jel:je  ,1:km)
        pdata%u(  iby,1:im,1:jm,1:kbl  ) = u(  iby,1:im,jel:je-1,1:km)
        pdata%u(  ibz,1:im,1:jm,1:kbl-1) = u(  ibz,1:im,jel:je  ,1:km)
#endif /* FLUXCT */
#endif /* MHD */
      else
        kl = kb
        ku = kbu

        ks = keu
        kt = km
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

    use blocks, only : block_data, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : ng, in, im, ih, ib, ibl, ibu, ie, iel, ieu              &
                     , nd, jn, jm, jh, jb, jbl, jbu, je, jel, jeu              &
                     , nh, kn, km, kh, kb, kbl, kbu, ke, kel, keu

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
#ifdef FIELDCD
! update magnetic field component boundaries
!
#if NDIMS == 2
    pdata%u(ibx:ibz,is:it,js:jt,1) = 0.25 * (u(ibx:ibz,il:iu:2,jl:ju:2,1)      &
                                           + u(ibx:ibz,ip:iu:2,jl:ju:2,1)      &
                                           + u(ibx:ibz,il:iu:2,jp:ju:2,1)      &
                                           + u(ibx:ibz,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    pdata%u(ibx:ibz,is:it,js:jt,ks:kt) =                                       &
                                  0.125 * (u(ibx:ibz,il:iu:2,jl:ju:2,kl:ku:2)  &
                                         + u(ibx:ibz,ip:iu:2,jl:ju:2,kl:ku:2)  &
                                         + u(ibx:ibz,il:iu:2,jp:ju:2,kl:ku:2)  &
                                         + u(ibx:ibz,ip:iu:2,jp:ju:2,kl:ku:2)  &
                                         + u(ibx:ibz,il:iu:2,jl:ju:2,kp:ku:2)  &
                                         + u(ibx:ibz,ip:iu:2,jl:ju:2,kp:ku:2)  &
                                         + u(ibx:ibz,il:iu:2,jp:ju:2,kp:ku:2)  &
                                         + u(ibx:ibz,ip:iu:2,jp:ju:2,kp:ku:2))
#endif /* NDIMS == 3 */
#endif /* FIELDCD */
#ifdef FLUXCT
! update staggered magnetic field component boundaries
!
#if NDIMS == 2
    if (it .eq. ibl) then
      pdata%u(ibx,is:it-1,js:jt,1) = 0.5 * (u(ibx,ip:iu-2:2,jl:ju:2,1)         &
                                          + u(ibx,ip:iu-2:2,jp:ju:2,1))
    else
      pdata%u(ibx,is:it  ,js:jt,1) = 0.5 * (u(ibx,ip:iu  :2,jl:ju:2,1)         &
                                          + u(ibx,ip:iu  :2,jp:ju:2,1))
    end if

    if (jt .eq. jbl) then
      pdata%u(iby,is:it,js:jt-1,1) = 0.5 * (u(iby,il:iu:2,jp:ju-2:2,1)         &
                                          + u(iby,ip:iu:2,jp:ju-2:2,1))
    else
      pdata%u(iby,is:it,js:jt  ,1) = 0.5 * (u(iby,il:iu:2,jp:ju  :2,1)         &
                                          + u(iby,ip:iu:2,jp:ju  :2,1))
    end if

    pdata%u(ibz,is:it,js:jt,1) = 0.25 * (u(ibz,il:iu:2,jl:ju:2,1)              &
                                       + u(ibz,ip:iu:2,jl:ju:2,1)              &
                                       + u(ibz,il:iu:2,jp:ju:2,1)              &
                                       + u(ibz,ip:iu:2,jp:ju:2,1))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    if (it .eq. ibl) then
      pdata%u(ibx,is:it-1,js:jt,ks:kt) =                                       &
                                     0.25 * (u(ibx,ip:iu-2:2,jl:ju:2,kl:ku:2)  &
                                           + u(ibx,ip:iu-2:2,jp:ju:2,kl:ku:2)  &
                                           + u(ibx,ip:iu-2:2,jl:ju:2,kp:ku:2)  &
                                           + u(ibx,ip:iu-2:2,jp:ju:2,kp:ku:2))
    else
      pdata%u(ibx,is:it  ,js:jt,ks:kt) =                                       &
                                     0.25 * (u(ibx,ip:iu  :2,jl:ju:2,kl:ku:2)  &
                                           + u(ibx,ip:iu  :2,jp:ju:2,kl:ku:2)  &
                                           + u(ibx,ip:iu  :2,jl:ju:2,kp:ku:2)  &
                                           + u(ibx,ip:iu  :2,jp:ju:2,kp:ku:2))
    end if

    if (jt .eq. jbl) then
      pdata%u(iby,is:it,js:jt-1,ks:kt) =                                       &
                                     0.25 * (u(iby,il:iu:2,jp:ju-2:2,kl:ku:2)  &
                                           + u(iby,ip:iu:2,jp:ju-2:2,kl:ku:2)  &
                                           + u(iby,il:iu:2,jp:ju-2:2,kp:ku:2)  &
                                           + u(iby,ip:iu:2,jp:ju-2:2,kp:ku:2))
    else
      pdata%u(iby,is:it,js:jt  ,ks:kt) =                                       &
                                     0.25 * (u(iby,il:iu:2,jp:ju  :2,kl:ku:2)  &
                                           + u(iby,ip:iu:2,jp:ju  :2,kl:ku:2)  &
                                           + u(iby,il:iu:2,jp:ju  :2,kp:ku:2)  &
                                           + u(iby,ip:iu:2,jp:ju  :2,kp:ku:2))
    end if

    if (kt .eq. kbl) then
      pdata%u(ibz,is:it,js:jt,ks:kt-1) =                                       &
                                     0.25 * (u(ibz,il:iu:2,jl:ju:2,kl:ku-2:2)  &
                                           + u(ibz,ip:iu:2,jl:ju:2,kl:ku-2:2)  &
                                           + u(ibz,il:iu:2,jp:ju:2,kl:ku-2:2)  &
                                           + u(ibz,ip:iu:2,jp:ju:2,kl:ku-2:2))
    else
      pdata%u(ibz,is:it,js:jt,ks:kt  ) =                                       &
                                     0.25 * (u(ibz,il:iu:2,jl:ju:2,kl:ku  :2)  &
                                           + u(ibz,ip:iu:2,jl:ju:2,kl:ku  :2)  &
                                           + u(ibz,il:iu:2,jp:ju:2,kl:ku  :2)  &
                                           + u(ibz,ip:iu:2,jp:ju:2,kl:ku  :2))
    end if
#endif /* NDIMS == 3 */
#endif /* FLUXCT */
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

    use blocks, only : block_data, nvr, nfl, nqt
#ifdef MHD
    use blocks, only : ibx, iby, ibz
#endif /* MHD */
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                  &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                  &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : expand_tvd
#if defined MHD && defined FLUXCT
    use interpolation, only : expand_mag_bnd
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
    integer :: ih, jh, kh, nh
    integer :: if, jf, kf, id, jd, kd
    integer :: il, iu, jl, ju, kl, ku
    integer :: is, it, js, jt, ks, kt

! local arrays
!
    integer, dimension(3)                    :: dm, cm, pm
    real   , dimension(:,:,:)  , allocatable :: u
#if defined MHD && defined FLUXCT
    real   , dimension(:,:,:)  , allocatable :: bx, by, bz
    real   , dimension(:,:)    , allocatable :: bn
#endif /* MHD & FLUXCT */

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
      call expand_tvd(cm, dm, dl, ub(q,il:iu,jl:ju,kl:ku), u(:,:,:))

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(:,:,:)

    end do
#ifdef MHD
#ifdef FIELDCD
! iterate over magnetic field components
!
    do q = ibx, ibz

! expand the boundary
!
      call expand_tvd(cm, dm, dl, ub(q,il:iu,jl:ju,kl:ku), u(:,:,:))

! copy expanded boundary in the proper place of the block
!
      pdata%u(q,is:it,js:jt,ks:kt) = u(:,:,:)

    end do
#endif /* FIELDCD */
#ifdef FLUXCT
! allocate space for the prolongated magnetic field components
!
    allocate(bx(dm(1)+1,dm(2)  ,dm(3)  ))
    allocate(by(dm(1)  ,dm(2)+1,dm(3)  ))
    allocate(bz(dm(1)  ,dm(2)  ,dm(3)+1))

! depending on the direction perform the right prolongation
!
    select case(idir)

!! boundary update along the X direction
!!
    case(1)

! prepare the index along the longitudinal direction from which the boundary
! should be copied
!
      id = ibl + if * in

! prepare the indices in ourder to determine the input array
!
      if (if .eq. 0) then
        is = ie - nh
        it = ie
      else
        is = ib      - 1
        it = ib + nh - 1
      end if
      il = is - dl + 1
      iu = it + dl

      jl = jb - nh - dl + jf * jh
      js = jl - 1
      ju = jl + cm(2) - 1

#if NDIMS == 2
      kl = 1
      ks = 1
      ku = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      kl = kb - nh - dl + kf * kh
      ks = kl - 1
      ku = kl + cm(3) - 1
#endif /* NDIMS == 3 */

! allocate the temporary boundary array
!
      allocate(bn(jm,km))

! copy the boundary from the higher refinement level neighbor
!
      bn(1:jm,1:km) = pdata%u(ibx,id,1:jm,1:km)

! prolongate magnetic field preserving the divergence-free condition
!
      call expand_mag_bnd(idir, if, dl, pm, bn(:,:)                            &
                 , ub(ibx,is:it,jl:ju,kl:ku), ub(iby,il:iu,js:ju,kl:ku)        &
                 , ub(ibz,il:iu,jl:ju,ks:ku), bx(:,:,:), by(:,:,:), bz(:,:,:))

! deallocate the boundary array
!
      deallocate(bn)

! update the magnetic field components
!
      if (if .eq. 0) then
        pdata%u(ibx,  1:ibl-1,1:jm,1:km) = bx(2:ng  ,1:jm  ,1:km  )
        pdata%u(iby,  1:ibl  ,1:jm,1:km) = by(1:ng  ,2:jm+1,1:km  )
#if NDIMS == 3
        pdata%u(ibz,  1:ibl  ,1:jm,1:km) = bz(1:ng  ,1:jm  ,2:km+1)
#endif /* NDIMS == 3 */
      else
        pdata%u(ibx,ieu:im   ,1:jm,1:km) = bx(2:ng+1,1:jm  ,1:km  )
        pdata%u(iby,ieu:im   ,1:jm,1:km) = by(1:ng  ,2:jm+1,1:km  )
#if NDIMS == 3
        pdata%u(ibz,ieu:im   ,1:jm,1:km) = bz(1:ng  ,1:jm  ,2:km+1)
#endif /* NDIMS == 3 */
      end if

!! boundary update along the Y direction
!!
    case(2)

! prepare the index along the longitudinal direction from which the boundary
! should be copied
!
      jd = jbl + jf * jn

! prepare the indices in ourder to determine the input array
!
      il = ib - nh - dl + if * ih
      is = il - 1
      iu = il + cm(1) - 1

      if (jf .eq. 0) then
        js = je - nh
        jt = je
      else
        js = jb      - 1
        jt = jb + nh - 1
      end if
      jl = js - dl + 1
      ju = jt + dl

#if NDIMS == 2
      kl = 1
      ks = 1
      ku = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      kl = kb - nh - dl + kf * kh
      ks = kl - 1
      ku = kl + cm(3) - 1
#endif /* NDIMS == 3 */

! allocate the temporary boundary array
!
      allocate(bn(im,km))

! copy the boundary from the higher refinement level neighbor
!
      bn(1:im,1:km) = pdata%u(iby,1:im,jd,1:km)

! prolongate magnetic field preserving the divergence-free condition
!
      call expand_mag_bnd(idir, jf, dl, pm, bn(:,:)                            &
                 , ub(ibx,is:iu,jl:ju,kl:ku), ub(iby,il:iu,js:jt,kl:ku)        &
                 , ub(ibz,il:iu,jl:ju,ks:ku), bx(:,:,:), by(:,:,:), bz(:,:,:))

! deallocate the boundary array
!
      deallocate(bn)

! update the magnetic field components
!
      if (jf .eq. 0) then
        pdata%u(ibx,1:im,  1:jbl  ,1:km) = bx(2:im+1,1:ng  ,1:km  )
        pdata%u(iby,1:im,  1:jbl-1,1:km) = by(1:im  ,2:ng  ,1:km  )
#if NDIMS == 3
        pdata%u(ibz,1:im,  1:jbl  ,1:km) = bz(1:im  ,1:ng  ,2:km+1)
#endif /* NDIMS == 3 */
      else
        pdata%u(ibx,1:im,jeu:jm   ,1:km) = bx(2:im+1,1:ng  ,1:km  )
        pdata%u(iby,1:im,jeu:jm   ,1:km) = by(1:im  ,2:ng+1,1:km  )
#if NDIMS == 3
        pdata%u(ibz,1:im,jeu:jm   ,1:km) = bz(1:im  ,1:ng  ,2:km+1)
#endif /* NDIMS == 3 */
      end if

#if NDIMS == 3
!! boundary update along the Y direction
!!
    case(3)

! prepare the index along the longitudinal direction from which the boundary
! should be copied
!
      kd = kbl + kf * kn

! prepare the indices in ourder to determine the input array
!
      il = ib - nh - dl + if * ih
      is = il - 1
      iu = il + cm(1) - 1

      jl = jb - nh - dl + jf * jh
      js = jl - 1
      ju = jl + cm(2) - 1

      if (kf .eq. 0) then
        ks = ke - nh
        kt = ke
      else
        ks = kb      - 1
        kt = kb + nh - 1
      end if
      kl = ks - dl + 1
      ku = kt + dl

! allocate the temporary boundary array
!
      allocate(bn(im,jm))

! copy the boundary from the higher refinement level neighbor
!
      bn(1:im,1:jm) = pdata%u(ibz,1:im,1:jm,kd)

! prolongate magnetic field preserving the divergence-free condition
!
      call expand_mag_bnd(idir, kf, dl, pm, bn(:,:)                            &
                 , ub(ibx,is:iu,jl:ju,kl:ku), ub(iby,il:iu,js:ju,kl:ku)        &
                 , ub(ibz,il:iu,jl:ju,ks:kt), bx(:,:,:), by(:,:,:), bz(:,:,:))

! deallocate the boundary array
!
      deallocate(bn)

! update the magnetic field components
!
      if (kf .eq. 0) then
        pdata%u(ibx,1:im,1:jm,  1:kbl  ) = bx(2:im+1,1:jm  ,1:ng  )
        pdata%u(iby,1:im,1:jm,  1:kbl  ) = by(1:im  ,2:jm+1,1:ng  )
        pdata%u(ibz,1:im,1:jm,  1:kbl-1) = bz(1:im  ,1:jm  ,2:ng  )
      else
        pdata%u(ibx,1:im,1:jm,keu:km   ) = bx(2:im+1,1:jm  ,1:ng  )
        pdata%u(iby,1:im,1:jm,keu:km   ) = by(1:im  ,2:jm+1,1:ng  )
        pdata%u(ibz,1:im,1:jm,keu:km   ) = bz(1:im  ,1:jm  ,2:ng+1)
      end if
#endif /* NDIMS == 3 */

    case default
    end select

! deallocate space for the prolongated magnetic field components
!
    deallocate(bx)
    deallocate(by)
    deallocate(bz)
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
