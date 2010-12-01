!!******************************************************************************
!!
!! module: scheme - handling the actual solver of the set of equations
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
module scheme

  implicit none

  contains
!
!===============================================================================
!
! numerical_flux: subroutine updates the numerical flux of the current data
!                 block
!
! notes: - the fluid fluxes must be located at cell interface centers
!             [ fx(i+1/2,j,k), fy(i,j+1/2,k), fz(i,j,k+1/2) ]
!        - the magnetic fluxes must be located at the cell edge centers
!             [ ex(i,j+1/2,j+1/2), ey(i+1/2,j,j+1/2), ez(i+1/2,j+1/2,k) ]
!
!===============================================================================
!
  subroutine numerical_flux(u, f, e)

    use blocks       , only : block_data
    use config       , only : im, jm, km
#ifdef FLUXCT
    use interpolation, only : magtocen
#endif /* FLUXCT */
    use variables    , only : nvr, nqt
    use variables    , only : idn, imx, imy, imz
#ifdef ADI
    use variables    , only : ien
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input/output arguments
!
    real, dimension(      nqt,im,jm,km), intent(in)  :: u
    real, dimension(NDIMS,nqt,im,jm,km), intent(out) :: f
    real, dimension(        3,im,jm,km), intent(out) :: e

! local variables
!
    integer :: i, j, k
#ifdef MHD
    integer :: im1, jm1, km1
#endif /* MHD */

! local temporary arrays
!
    real, dimension(nvr,im)       :: ux
    real, dimension(nqt,im)       :: fx
    real, dimension(nvr,jm)       :: uy
    real, dimension(nqt,jm)       :: fy
#if NDIMS == 3
    real, dimension(nvr,km)       :: uz
    real, dimension(nqt,km)       :: fz
#endif /* NDIMS == 3 */
    real, dimension(nvr,im,jm,km) :: w
!
!-------------------------------------------------------------------------------
!
! 0. prepare the numerical flux calculation
!
! - reset the flux array to zero
!
    f(:,:,:,:,:) = 0.0
    e(:,:,:,:)   = 0.0

! - copy the conserved variables to the local array
!
    w(1:nqt,1:im,1:jm,1:km) = u(1:nqt,1:im,1:jm,1:km)
#ifdef FLUXCT
! - interpolate the cell centered components of the magnetic field
!
    call magtocen(w)
#endif /* FLUXCT */

! 1. update along the X direction
!
    do k = 1, km
#if NDIMS == 3 && defined MHD
      km1 = max(k - 1, 1)
#endif /* NDIMS == 3 & MHD */
      do j = 1, jm
#ifdef MHD
        jm1 = max(j - 1, 1)
#endif /* MHD */

! - copy the directional vectors of variables for the one dimensional solver
!
        ux(idn,1:im) = w(idn,1:im,j,k)
        ux(imx,1:im) = w(imx,1:im,j,k)
        ux(imy,1:im) = w(imy,1:im,j,k)
        ux(imz,1:im) = w(imz,1:im,j,k)
#ifdef ADI
        ux(ien,1:im) = w(ien,1:im,j,k)
#endif /* ADI */
#ifdef MHD
        ux(ibx,1:im) = w(ibx,1:im,j,k)
        ux(iby,1:im) = w(iby,1:im,j,k)
        ux(ibz,1:im) = w(ibz,1:im,j,k)
#ifdef FLUXCT
        ux(icx,1:im) = w(icx,1:im,j,k)
        ux(icy,1:im) = w(icy,1:im,j,k)
        ux(icz,1:im) = w(icz,1:im,j,k)
#endif /* FLUXCT */
#endif /* MHD */

! - execute the Riemann solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (im, ux, fx, .true.)
#endif /* HLL */
#ifdef HLLC
        call hllc(im, ux, fx, .true.)
#endif /* HLLC */

! - update the fluxes along the current direction
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
#ifdef ADI
        f(1,ien,1:im,j,k) = fx(ien,1:im)
#endif /* ADI */
#ifdef MHD
#if NDIMS == 2
        e(2,1:im,j  ,k  ) = e(2,1:im,j  ,k  ) +        fx(ibz,1:im)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        e(2,1:im,j  ,k  ) = e(2,1:im,j  ,k  ) + 0.25 * fx(ibz,1:im)
        e(2,1:im,j  ,km1) = e(2,1:im,j  ,km1) + 0.25 * fx(ibz,1:im)
#endif /* NDIMS == 3 */
        e(3,1:im,j  ,k  ) = e(3,1:im,j  ,k  ) - 0.25 * fx(iby,1:im)
        e(3,1:im,jm1,k  ) = e(3,1:im,jm1,k  ) - 0.25 * fx(iby,1:im)
#endif /* MHD */
      end do
    end do

! 2. update along the Y direction
!
    do k = 1, km
#if NDIMS == 3 && defined MHD
      km1 = max(k - 1, 1)
#endif /* NDIMS == 3 & MHD */
      do i = 1, im
#ifdef MHD
        im1 = max(i - 1, 1)
#endif /* MHD */

! - copy the directional vectors of variables for the one dimensional solver
!
        uy(idn,1:jm) = w(idn,i,1:jm,k)
        uy(imx,1:jm) = w(imy,i,1:jm,k)
        uy(imy,1:jm) = w(imz,i,1:jm,k)
        uy(imz,1:jm) = w(imx,i,1:jm,k)
#ifdef ADI
        uy(ien,1:jm) = w(ien,i,1:jm,k)
#endif /* ADI */
#ifdef MHD
        uy(ibx,1:jm) = w(iby,i,1:jm,k)
        uy(iby,1:jm) = w(ibz,i,1:jm,k)
        uy(ibz,1:jm) = w(ibx,i,1:jm,k)
#ifdef FLUXCT
        uy(icx,1:jm) = w(icy,i,1:jm,k)
        uy(icy,1:jm) = w(icz,i,1:jm,k)
        uy(icz,1:jm) = w(icx,i,1:jm,k)
#endif /* FLUXCT */
#endif /* MHD */

! - execute the Riemann solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (jm, uy, fy, .true.)
#endif /* HLL */
#ifdef HLLC
        call hllc(jm, uy, fy, .true.)
#endif /* HLLC */

! - update the fluxes along the current direction
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
#ifdef ADI
        f(2,ien,i,1:jm,k) = fy(ien,1:jm)
#endif /* ADI */
#ifdef MHD
#if NDIMS == 2
        e(1,i  ,1:jm,k  ) = e(1,i  ,1:jm,k  ) -        fy(iby,1:jm)
#endif /* NDIMS == 2 */
#if NDIMS == 3
        e(1,i  ,1:jm,k  ) = e(1,i  ,1:jm,k  ) - 0.25 * fy(iby,1:jm)
        e(1,i  ,1:jm,km1) = e(1,i  ,1:jm,km1) - 0.25 * fy(iby,1:jm)
#endif /* NDIMS == 3 */
        e(3,i  ,1:jm,k  ) = e(3,i  ,1:jm,k  ) + 0.25 * fy(ibz,1:jm)
        e(3,im1,1:jm,k  ) = e(3,im1,1:jm,k  ) + 0.25 * fy(ibz,1:jm)
#endif /* MHD */
      end do
    end do

#if NDIMS == 3
! 3. update along the Z direction
!
    do j = 1, jm
#ifdef MHD
      jm1 = max(j - 1, 1)
#endif /* MHD */
      do i = 1, im
#ifdef MHD
        im1 = max(i - 1, 1)
#endif /* MHD */

! - copy the directional vectors of variables for the one dimensional solver
!
        uz(idn,1:km) = w(idn,i,j,1:km)
        uz(imx,1:km) = w(imz,i,j,1:km)
        uz(imy,1:km) = w(imx,i,j,1:km)
        uz(imz,1:km) = w(imy,i,j,1:km)
#ifdef ADI
        uz(ien,1:km) = w(ien,i,j,1:km)
#endif /* ADI */
#ifdef MHD
        uz(ibx,1:km) = w(ibz,i,j,1:km)
        uz(iby,1:km) = w(ibx,i,j,1:km)
        uz(ibz,1:km) = w(iby,i,j,1:km)
#ifdef FLUXCT
        uz(icx,1:km) = w(icz,i,j,1:km)
        uz(icy,1:km) = w(icx,i,j,1:km)
        uz(icz,1:km) = w(icy,i,j,1:km)
#endif /* FLUXCT */
#endif /* MHD */

! - execute the Riemann solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (km, uz, fz, .true.)
#endif /* HLL */
#ifdef HLLC
        call hllc(km, uz, fz, .true.)
#endif /* HLLC */

! - update the fluxes along the current direction
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
#ifdef ADI
        f(3,ien,i,j,1:km) = fz(ien,1:km)
#endif /* ADI */
#ifdef MHD
        e(1,i  ,j  ,1:km) = e(1,i  ,j  ,1:km) + 0.25 * fz(ibz,1:km)
        e(1,i  ,jm1,1:km) = e(1,i  ,jm1,1:km) + 0.25 * fz(ibz,1:km)
        e(2,i  ,j  ,1:km) = e(2,i  ,j  ,1:km) - 0.25 * fz(iby,1:km)
        e(2,im1,j  ,1:km) = e(2,im1,j  ,1:km) - 0.25 * fz(iby,1:km)
#endif /* MHD */
      end do
    end do
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
  end subroutine numerical_flux
!
!===============================================================================
!
! update: subroutine sweeps over all directions and integrates the directional
!         derivatives of the flux in order to get the increment of solution
!
!===============================================================================
!
  subroutine update(u, du, dxi, dyi, dzi)

    use config       , only : im, jm, km
#ifdef FLUXCT
    use interpolation, only : magtocen
#endif /* FLUXCT */
    use variables    , only : nvr, nqt, nfl
    use variables    , only : idn, imx, imy, imz
#ifdef ADI
    use variables    , only : ien
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input arguments
!
    real, dimension(nqt,im,jm,km)      , intent(in)  :: u
    real, dimension(nqt,im,jm,km)      , intent(out) :: du
    real                               , intent(in)  :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k, im1, jm1, km1, ip1, jp1, kp1

! local temporary arrays
!
    real, dimension(nvr,im)       :: ux
    real, dimension(nqt,im)       :: fx
    real, dimension(nvr,jm)       :: uy
    real, dimension(nqt,jm)       :: fy
#if NDIMS == 3
    real, dimension(nvr,km)       :: uz
    real, dimension(nqt,km)       :: fz
#endif /* NDIMS == 3 */
    real, dimension(nvr,im,jm,km) :: w
!
!-------------------------------------------------------------------------------
!
! copy variables to a new array
!
    w(1:nqt,:,:,:) = u(1:nqt,:,:,:)
#ifdef FLUXCT
! interpolate cell centered components of the magnetic field
!
    call magtocen(w)
#endif /* FLUXCT */

! reset the increment array
!
    du(:,:,:,:) = 0.0

!  update along X-direction
!
    do k = 1, km
#if NDIMS == 3
#ifdef MHD
      km1 = max(k - 1, 1)
      kp1 = min(k + 1,km)
#endif /* MHD */
#endif /* NDIMS == 3 */
      do j = 1, jm
#ifdef MHD
        jm1 = max(j - 1, 1)
        jp1 = min(j + 1,jm)
#endif /* MHD */

! copy directional vectors of variables for the one dimensional solver
!
        do i = 1, im
          ux(idn,i) = w(idn,i,j,k)
          ux(imx,i) = w(imx,i,j,k)
          ux(imy,i) = w(imy,i,j,k)
          ux(imz,i) = w(imz,i,j,k)
#ifdef ADI
          ux(ien,i) = w(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          ux(ibx,i) = w(ibx,i,j,k)
          ux(iby,i) = w(iby,i,j,k)
          ux(ibz,i) = w(ibz,i,j,k)
#ifdef FLUXCT
          ux(icx,i) = w(icx,i,j,k)
          ux(icy,i) = w(icy,i,j,k)
          ux(icz,i) = w(icz,i,j,k)
#endif /* FLUXCT */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (im, ux, fx, .false.)
#endif /* HLL */
#ifdef HLLC
        call hllc(im, ux, fx, .false.)
#endif /* HLLC */

! update the arrays of increments
!
        do i = 1, im
          du(idn,i,j,k) = du(idn,i,j,k) + dxi * fx(idn,i)
          du(imx,i,j,k) = du(imx,i,j,k) + dxi * fx(imx,i)
          du(imy,i,j,k) = du(imy,i,j,k) + dxi * fx(imy,i)
          du(imz,i,j,k) = du(imz,i,j,k) + dxi * fx(imz,i)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dxi * fx(ien,i)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          im1 = max(i - 1, 1)
#ifdef FIELDCD
          ip1 = min(i + 1,im)

          du(iby,i,j,k) = du(iby,i,j,k) - 0.5 * dxi * (fx(iby,ip1) - fx(iby,im1))
#if NDIMS == 3
          du(ibz,i,j,k) = du(ibz,i,j,k) - 0.5 * dxi * (fx(ibz,ip1) - fx(ibz,im1))
#endif /* NDIMS == 3 */
#endif /* FIELDCD */
#ifdef FLUXCT
          du(ibx,i  ,jm1,k  ) = du(ibx,i  ,jm1,k  ) + dyi *  fx(iby,i)
          du(ibx,i  ,jp1,k  ) = du(ibx,i  ,jp1,k  ) - dyi *  fx(iby,i)
          du(iby,i  ,j  ,k  ) = du(iby,i  ,j  ,k  ) - dxi * (fx(iby,i) - fx(iby,im1))
          du(iby,i  ,jm1,k  ) = du(iby,i  ,jm1,k  ) - dxi * (fx(iby,i) - fx(iby,im1))
#if NDIMS == 3
          du(ibx,i  ,j  ,km1) = du(ibx,i  ,j  ,km1) + dzi *  fx(ibz,i)
          du(ibx,i  ,j  ,kp1) = du(ibx,i  ,j  ,kp1) - dzi *  fx(ibz,i)
          du(ibz,i  ,j  ,k  ) = du(ibz,i  ,j  ,k  ) - dxi * (fx(ibz,i) - fx(ibz,im1))
          du(ibz,i  ,j  ,km1) = du(ibz,i  ,j  ,km1) - dxi * (fx(ibz,i) - fx(ibz,im1))
#endif /* NDIMS == 3 */
#endif /* FLUXCT */
#endif /* MHD */
        end do
      end do
    end do

!  update along Y-direction
!
    do k = 1, km
#if NDIMS == 3
#ifdef MHD
      km1 = max(k - 1, 1)
      kp1 = min(k + 1,km)
#endif /* MHD */
#endif /* NDIMS == 3 */
      do i = 1, im
#ifdef MHD
        im1 = max(i - 1, 1)
        ip1 = min(i + 1,im)
#endif /* MHD */

! copy directional vectors of variables for the one dimensional solver
!
        do j = 1, jm
          uy(idn,j) = w(idn,i,j,k)
          uy(imx,j) = w(imy,i,j,k)
          uy(imy,j) = w(imz,i,j,k)
          uy(imz,j) = w(imx,i,j,k)
#ifdef ADI
          uy(ien,j) = w(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          uy(ibx,j) = w(iby,i,j,k)
          uy(iby,j) = w(ibz,i,j,k)
          uy(ibz,j) = w(ibx,i,j,k)
#ifdef FLUXCT
          uy(icx,j) = w(icy,i,j,k)
          uy(icy,j) = w(icz,i,j,k)
          uy(icz,j) = w(icx,i,j,k)
#endif /* FLUXCT */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (jm, uy, fy, .false.)
#endif /* HLL */
#ifdef HLLC
        call hllc(jm, uy, fy, .false.)
#endif /* HLLC */

! update the arrays of increments
!
        do j = 1, jm
          du(idn,i,j,k) = du(idn,i,j,k) + dyi * fy(idn,j)
          du(imx,i,j,k) = du(imx,i,j,k) + dyi * fy(imz,j)
          du(imy,i,j,k) = du(imy,i,j,k) + dyi * fy(imx,j)
          du(imz,i,j,k) = du(imz,i,j,k) + dyi * fy(imy,j)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dyi * fy(ien,j)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          jm1 = max(j - 1, 1)
#ifdef FIELDCD
          jp1 = min(j + 1,jm)

#if NDIMS == 3
          du(ibz,i,j,k) = du(ibz,i,j,k) - 0.5 * dyi * (fy(iby,jp1) - fy(iby,jm1))
#endif /* NDIMS == 3 */
          du(ibx,i,j,k) = du(ibx,i,j,k) - 0.5 * dyi * (fy(ibz,jp1) - fy(ibz,jm1))
#endif /* FIELDCD */
#ifdef FLUXCT
#if NDIMS == 3
          du(iby,i  ,j  ,km1) = du(iby,i  ,j  ,km1) + dzi *  fy(iby,j)
          du(iby,i  ,j  ,kp1) = du(iby,i  ,j  ,kp1) - dzi *  fy(iby,j)
          du(ibz,i  ,j  ,k  ) = du(ibz,i  ,j  ,k  ) - dyi * (fy(iby,j) - fy(iby,jm1))
          du(ibz,i  ,j  ,km1) = du(ibz,i  ,j  ,km1) - dyi * (fy(iby,j) - fy(iby,jm1))
#endif /* NDIMS == 3 */
          du(iby,im1,j  ,k  ) = du(iby,im1,j  ,k  ) + dxi *  fy(ibz,j)
          du(iby,ip1,j  ,k  ) = du(iby,ip1,j  ,k  ) - dxi *  fy(ibz,j)
          du(ibx,i  ,j  ,k  ) = du(ibx,i  ,j  ,k  ) - dyi * (fy(ibz,j) - fy(ibz,jm1))
          du(ibx,im1,j  ,k  ) = du(ibx,im1,j  ,k  ) - dyi * (fy(ibz,j) - fy(ibz,jm1))
#endif /* FLUXCT */
#endif /* MHD */
        end do
      end do
    end do

#if NDIMS == 3
!  update along Z-direction
!
    do j = 1, jm
#ifdef MHD
      jm1 = max(j - 1, 1)
      jp1 = min(j + 1,jm)
#endif /* MHD */
      do i = 1, im
#ifdef MHD
        im1 = max(i - 1, 1)
        ip1 = min(i + 1,im)
#endif /* MHD */

! copy directional vectors of variables for the one dimensional solver
!
        do k = 1, km
          uz(idn,k) = w(idn,i,j,k)
          uz(imx,k) = w(imz,i,j,k)
          uz(imy,k) = w(imx,i,j,k)
          uz(imz,k) = w(imy,i,j,k)
#ifdef ADI
          uz(ien,k) = w(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          uz(ibx,k) = w(ibz,i,j,k)
          uz(iby,k) = w(ibx,i,j,k)
          uz(ibz,k) = w(iby,i,j,k)
#ifdef FLUXCT
          uz(icx,k) = w(icz,i,j,k)
          uz(icy,k) = w(icx,i,j,k)
          uz(icz,k) = w(icy,i,j,k)
#endif /* FLUXCT */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (km, uz, fz, .false.)
#endif /* HLL */
#ifdef HLLC
        call hllc(km, uz, fz, .false.)
#endif /* HLLC */

! update the arrays of increments
!
        do k = 1, km
          du(idn,i,j,k) = du(idn,i,j,k) + dzi * fz(idn,k)
          du(imx,i,j,k) = du(imx,i,j,k) + dzi * fz(imy,k)
          du(imy,i,j,k) = du(imy,i,j,k) + dzi * fz(imz,k)
          du(imz,i,j,k) = du(imz,i,j,k) + dzi * fz(imx,k)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dzi * fz(ien,k)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          km1 = max(k - 1, 1)
#ifdef FIELDCD
          kp1 = min(k + 1,km)

          du(ibx,i,j,k) = du(ibx,i,j,k) - 0.5 * dzi * (fz(iby,kp1) - fz(iby,km1))
          du(iby,i,j,k) = du(iby,i,j,k) - 0.5 * dzi * (fz(ibz,kp1) - fz(ibz,km1))
#endif /* FIELDCD */
#ifdef FLUXCT
          du(ibz,im1,j  ,k  ) = du(ibz,im1,j  ,k  ) + dxi *  fz(iby,k)
          du(ibz,ip1,j  ,k  ) = du(ibz,ip1,j  ,k  ) - dxi *  fz(iby,k)
          du(ibx,i  ,j  ,k  ) = du(ibx,i  ,j  ,k  ) - dzi * (fz(iby,k) - fz(iby,km1))
          du(ibx,im1,j  ,k  ) = du(ibx,im1,j  ,k  ) - dzi * (fz(iby,k) - fz(iby,km1))

          du(ibz,i  ,jm1,k  ) = du(ibz,i  ,jm1,k  ) + dyi *  fz(ibz,k)
          du(ibz,i  ,jp1,k  ) = du(ibz,i  ,jp1,k  ) - dyi *  fz(ibz,k)
          du(iby,i  ,j  ,k  ) = du(iby,i  ,j  ,k  ) - dzi * (fz(ibz,k) - fz(ibz,km1))
          du(iby,i  ,jm1,k  ) = du(iby,i  ,jm1,k  ) - dzi * (fz(ibz,k) - fz(ibz,km1))
#endif /* FLUXCT */
#endif /* MHD */
        end do
      end do
    end do
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
  end subroutine update
#ifdef HLL
!
!===============================================================================
!
! hll: subroutine computes the approximated flux using HLL method
!
!===============================================================================
!
  subroutine hll(n, u, f, d)

    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : ivx, ivz
#ifdef MHD
    use variables    , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f
    logical               , intent(in)  :: d

! local variables
!
    integer                :: p, i
    real, dimension(nvr,n) :: q, ql, qr, ul, ur
    real, dimension(nqt,n) :: fl, fr, fn
    real, dimension(n)     :: cl, cr
    real                   :: al, ar, ap, div
!
!-------------------------------------------------------------------------------
!
! calculate the primitive variables
!
    call cons2prim(n, u, q)

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
      call reconstruct(n, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef MHD
#ifdef FLUXCT
! reconstruct the left and right states of the magnetic field components
!
    ql(ibx,:) = q(ibx,:)
    qr(ibx,:) = q(ibx,:)
    call reconstruct(n, q(icy,:), ql(iby,:), qr(iby,:))
    call reconstruct(n, q(icz,:), ql(ibz,:), qr(ibz,:))
    ql(icx:icz,:) = ql(ibx:ibz,:)
    qr(icx:icz,:) = qr(ibx:ibz,:)
#endif /* FLUXCT */
#endif /* MHD */

! calculate conservative variables at states
!
    call prim2cons(n, ql, ul)
    call prim2cons(n, qr, ur)

! calculate fluxes and speeds
!
    call fluxspeed(n, ql, ul, fl, cl)
    call fluxspeed(n, qr, ur, fr, cr)

! iterate over all points
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      al = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      ar = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate HLL flux
!
      if (al .ge. 0.0) then
        fn(:,i) = fl(:,i)
      else if (ar .le. 0.0) then
        fn(:,i) = fr(:,i)
      else
        ap  = ar * al
        div = 1.0 / (ar - al)

        fn(1:nqt,i) = div * (ar * fl(1:nqt,i) - al * fr(1:nqt,i)               &
                                             + ap * (ur(1:nqt,i) - ul(1:nqt,i)))
      end if
    end do

! calculate numerical flux
!
    if (d) then
      f(1:nqt,1:n) = fn(1:nqt,1:n)
    else
      f(  1:nfl,2:n) = - fn(  1:nfl,2:n) + fn(   1:nfl,1:n-1)
#ifdef MHD
#ifdef FIELDCD
      call emf(n, q(ivx:ivz,:), q(ibx:ibz,:), f(ibx:ibz,:))
#endif /* FIELDCD */
#ifdef FLUXCT
      f(ibx:ibz,1:n) = 0.25 * fn(ibx:ibz,1:n)
#endif /* FLUXCT */
#endif /* MHD */
    end if

!-------------------------------------------------------------------------------
!
  end subroutine hll
#endif /* HLL */
#ifdef HLLC
!===============================================================================
!
! hllc: subroutine to compute flux approximated by HLLC method (HYDRO only)
!         ([1] Batten et al., 1997, JSC, 18, 6, 1553)
!         ([2] Miyoshi & Kusano, 2005, JCP, 208, 315)
!
!===============================================================================
!
  subroutine hllc(n, u, f, d)

    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f
    logical               , intent(in)  :: d

! local variables
!
    integer                :: p, i
    real, dimension(nvr,n) :: ql, qr, qc, ul, ur
    real, dimension(nfl,n) :: fl, fr, fx
    real, dimension(n)     :: cl, cr, cm
    real                   :: sl, sr, sm, sml, smr, srmv, slmv, srmm, slmm     &
                            , smvl, smvr, div, pt
    real, dimension(nvr)   :: q1l, q1r, u1l, u1r
!
!-------------------------------------------------------------------------------
!
    f (:,:) = 0.0
    fx(:,:) = 0.0

! calculate primitive variables
!
    call cons2prim(n, uc, qc)

! reconstruct left and right states of primitive variables
!
    do p = 1, nfl
      call reconstruct(n, qc(p,:), ql(p,:), qr(p,:))
    enddo

! calculate conservative variables at states
!
    call prim2cons(n, ql, ul)
    call prim2cons(n, qr, ur)

! calculate fluxes and speeds
!
    call fluxspeed(n, ql, ul, fl, cl)
    call fluxspeed(n, qr, ur, fr, cr)

! calculate fluxes
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      sl = min(ql(2,i) - cl(i),qr(2,i) - cr(i))
      sr = max(ql(2,i) + cl(i),qr(2,i) + cr(i))

! all speeds >= 0, left side flux
!
      if (sl .ge. 0.0) then

        fx(:,i) = fl(:,i)

! all speeds <= 0, right side flux
!
      else if (sr .le. 0.0) then

        fx(:,i) = fr(:,i)

! intermediate states
!
      else  ! sl < 0 & sr > 0

! useful differences
!
        slmv = sl - ql(2,i)
        srmv = sr - qr(2,i)

! speed of contact discontinuity (eq. 34 [1], 14 [2])
!
        div =  srmv*qr(1,i) - slmv*ql(1,i)
        sml = (srmv*ur(2,i) - slmv*ul(2,i) - qr(5,i) + ql(5,i))/div
        div =  slmv*ql(1,i) - srmv*qr(1,i)
        smr = (slmv*ul(2,i) - srmv*ur(2,i) - ql(5,i) + qr(5,i))/div
        sm  = 0.5d0 * (sml + smr)

        if (sm .eq. 0.0d0) then

! calculate left intermediate state
!
          pt = ql(5,i) - ul(2,i)*slmv
          u1l(1) = ql(1,i)*slmv/sl
          u1l(2) = 0.0d0
          u1l(3) = u1l(1)*ql(3,i)
          u1l(4) = u1l(1)*ql(4,i)
          if (sl .eq. 0.0) then
            u1l(5) = ul(5,i)
          else
            u1l(5) = (slmv*ul(5,i) - ql(5,i)*ql(2,i))/sl
          endif

! calculate right intermediate state
!
            pt = qr(5,i) - ur(2,i)*srmv
            u1r(1) = qr(1,i)*srmv/sr
            u1r(2) = 0.0d0
            u1r(3) = u1r(1)*qr(3,i)
            u1r(4) = u1r(1)*qr(4,i)
            if (sr .eq. 0.0) then
              u1r(5) = ur(5,i)
            else
              u1r(5) = (srmv*ur(5,i) - qr(5,i)*qr(2,i))/sr
            endif

! calculate intermediate flux
!
            fx(:,i) = 0.5*(fl(:,i) + sl*(u1l(:) - ul(:,i)) &
                    +      fr(:,i) + sr*(u1r(:) - ur(:,i)))

        else

! useful differences
!
          slmm = sl - sm
          srmm = sr - sm
          smvl = sm - ql(2,i)
          smvr = sm - qr(2,i)

! intermediate discontinuities
!
          if (sm .gt. 0.0) then

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
            pt = ql(5,i) + ql(1,i)*slmv*smvl

! calculate left intermediate state
!
            u1l(1) = ql(1,i)*slmv/slmm                               ! eq. (35 [1], 17 [2])
            u1l(2) = u1l(1)*sm                                       ! eq. (37 [1])
            u1l(3) = u1l(1)*ql(3,i)                                  ! eq. (38 [1], 18 [2])
            u1l(4) = u1l(1)*ql(4,i)                                  ! eq. (39 [1], 19 [2])
            if (slmm .eq. 0.0) then
              u1l(5) = ul(5,i)
            else
              u1l(5) = (slmv*ul(5,i) - ql(5,i)*ql(2,i) + pt*sm)/slmm ! eq. (40 [1], 20 [2])
            endif

! calculate left intermediate flux
!
            fx(:,i) = fl(:,i) + sl*(u1l(:) - ul(:,i))    ! eq. (29 [1], 15 [2])

          else if (sm .lt. 0.0) then

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
            pt = qr(5,i) + qr(1,i)*srmv*smvr

! calculate right intermediate state
!
            u1r(1) = qr(1,i)*srmv/srmm                               ! eq. (35 [1], 17 [2])
            u1r(2) = u1r(1)*sm                                       ! eq. (37 [1])
            u1r(3) = u1r(1)*qr(3,i)                                  ! eq. (38 [1], 18 [2])
            u1r(4) = u1r(1)*qr(4,i)                                  ! eq. (39 [1], 19 [2])
            if (srmm .eq. 0.0) then
              u1r(5) = ur(5,i)
            else
              u1r(5) = (srmv*ur(5,i) - qr(5,i)*qr(2,i) + pt*sm)/srmm ! eq. (40 [1], 20 [2])
            endif

! calculate right intermediate flux
!
            fx(:,i) = fr(:,i) + sr*(u1r(:) - ur(:,i))    ! eq. (30 [1], 15 [2])

          endif
        endif
      endif

    enddo

! calculate numerical flux
!
    if (d) then
      f(:,:) = fx(:,:)
    else
      f(:,2:n) = - fx(:,2:n) + fx(:,1:n-1)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine hllc
#endif /* HLLC */
!
!===============================================================================
!
! fluxspeed: subroutine computes fluxes and speeds for a given set of equations
!
!===============================================================================
!
  subroutine fluxspeed(n, q, u, f, c)

    use config   , only : gamma, csnd, csnd2
    use variables, only : nvr, nqt
    use variables, only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr, ien
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real, dimension(nvr,n), intent(in)  :: q, u
    real, dimension(nqt,n), intent(out) :: f
    real, dimension(n)    , intent(out) :: c

! local variables
!
    integer :: i
    real    :: bb, pm, vb, cs, cb, ca
!
!-------------------------------------------------------------------------------
!
! sweep over all points
!
    do i = 1, n

! compute fluxes
!
      f(idn,i) = u(imx,i)
#ifdef ADI
      f(imx,i) = q(ivx,i) * u(imx,i) + q(ipr,i)
#endif /* ADI */
#ifdef ISO
      f(imx,i) = q(ivx,i) * u(imx,i) + q(idn,i) * csnd2
#endif /* ISO */
      f(imy,i) = q(ivx,i) * u(imy,i)
      f(imz,i) = q(ivx,i) * u(imz,i)
#ifdef ADI
      f(ien,i) = q(ivx,i) * (u(ien,i) + q(ipr,i))
#endif /* ADI */
#ifdef MHD
      bb       = sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      pm       = 0.5 * bb
      vb       = sum(q(ivx:ivz,i) * q(ibx:ibz,i))
      f(imx,i) = f(imx,i) - q(ibx,i) * q(ibx,i) + pm
      f(imy,i) = f(imy,i) - q(ibx,i) * q(iby,i)
      f(imz,i) = f(imz,i) - q(ibx,i) * q(ibz,i)
#ifdef ADI
      f(ien,i) = f(ien,i) + q(ivx,i) * pm - q(ibx,i) * vb
#endif /* ADI */
#ifdef FLUXCT
      f(ibx,i) = 0.0
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
#endif /* FLUXCT */
#endif /* MHD */

! compute speeds
!
#ifdef MHD
#ifdef ADI
      cs = gamma * q(ipr,i)
#endif /* ADI */
#ifdef ISO
      cs = csnd2 * q(idn,i)
#endif /* ISO */
      cb = cs + bb
      ca = q(ibx,i) * q(ibx,i)
      c(i) = sqrt(0.5 * (cb + sqrt(max(0.0, cb * cb - 4.0 * cs * ca))) / q(idn,i))
#else /* MHD */
#ifdef ADI
      c(i) = sqrt(gamma * q(ipr,i) / q(idn,i))
#endif /* ADI */
#ifdef ISO
      c(i) = csnd
#endif /* ISO */
#endif /* MHD */
    end do

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed
#ifdef MHD
!
!===============================================================================
!
! emf: subroutine computes magnetic fluxes (electromotive force)
!
!===============================================================================
!
  subroutine emf(n, v, b, f)

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: n
    real, dimension(3,n), intent(in)  :: v, b
    real, dimension(3,n), intent(out) :: f

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
! sweep over all points
!
    do i = 1, n
      f(1,i) = 0.0
      f(2,i) = v(1,i) * b(2,i) - b(1,i) * v(2,i)
      f(3,i) = v(1,i) * b(3,i) - b(1,i) * v(3,i)
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine emf
#endif /* MHD */
!
!===============================================================================
!
! cons2prim: subroutine converts primitive variables to conservative
!
!===============================================================================
!
  subroutine cons2prim(n, u, q)

    use config   , only : gammam1
    use variables, only : nvr
    use variables, only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr, ien
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nvr,n), intent(out) :: q

! local variables
!
    integer :: i
    real    :: dni, ei, ek, em
!
!-------------------------------------------------------------------------------
!
    do i = 1, n
      dni      = 1.0 / u(idn,i)

      q(idn,i) =       u(idn,i)
      q(ivx,i) = dni * u(imx,i)
      q(ivy,i) = dni * u(imy,i)
      q(ivz,i) = dni * u(imz,i)
#ifdef ADI
      ek       = 0.5 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      ei       = u(ien,i) - ek
#ifdef MHD
      em       = 0.5 * sum(u(icx:icz,i) * u(icx:icz,i))
      ei       = ei - em
#endif /* MHD */
      q(ipr,i) = gammam1 * ei
#endif /* ADI */
#ifdef MHD
      q(ibx,i) = u(ibx,i)
      q(iby,i) = u(iby,i)
      q(ibz,i) = u(ibz,i)
#ifdef FLUXCT
      q(icx,i) = u(icx,i)
      q(icy,i) = u(icy,i)
      q(icz,i) = u(icz,i)
#endif /* FLUXCT */
#endif /* MHD */
    end do

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim
!
!===============================================================================
!
! prim2cons: subroutine converts primitive variables to conservative
!
!===============================================================================
!
  subroutine prim2cons(n, q, u)

    use config   , only : gammam1i
    use variables, only : nvr
    use variables, only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr, ien
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real, dimension(nvr,n), intent(in)  :: q
    real, dimension(nvr,n), intent(out) :: u

! local variables
!
    integer :: i
    real    :: ei, ek, em
!
!-------------------------------------------------------------------------------
!
    do i = 1, n
      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
#ifdef ADI
      ei       = gammam1i * q(ipr,i)
      ek       = 0.5 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      u(ien,i) = ei + ek
#endif /* ADI */
#ifdef MHD
#ifdef ADI
      em       = 0.5 * sum(q(icx:icz,i) * q(icx:icz,i))
      u(ien,i) = u(ien,i) + em
#endif /* ADI */
      u(ibx,i) = q(ibx,i)
      u(iby,i) = q(iby,i)
      u(ibz,i) = q(ibz,i)
#ifdef FLUXCT
      u(icx,i) = q(icx,i)
      u(icy,i) = q(icy,i)
      u(icz,i) = q(icz,i)
#endif /* FLUXCT */
#endif /* MHD */
    end do

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons
!
!===============================================================================
!
! maxspeed: function to calculate maximum speed in the system
!
!===============================================================================
!
  function maxspeed(u)

    use config       , only : im, jm, km, ib, ie, jb, je, kb, ke
#ifdef ADI
    use config       , only : gamma
#endif /* ADI */
#ifdef ISO
    use config       , only : csnd, csnd2
#endif /* ISO */
#if defined MHD && defined FLUXCT
    use interpolation, only : magtocen
#endif /* MHD && FLUXCT */
    use variables    , only : nvr, nqt
    use variables    , only : idn, ivx, ivz
#ifdef ADI
    use variables    , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

    implicit none

! input arguments
!
    real, dimension(nqt,im,jm,km), intent(in)  :: u

! local variables
!
    integer :: i, j, k
    real    :: vv, v, c
#ifdef MHD
    real    :: bb
#endif /* MHD */
    real    :: maxspeed

! local arrays
!
    real, dimension(nvr,im,jm,km) :: w
    real, dimension(nvr,im)       :: q
!
!-------------------------------------------------------------------------------
!
    maxspeed = 0.0

! copy variables to a new array
!
    w(1:nqt,1:im,1:jm,1:km) = u(1:nqt,1:im,1:jm,1:km)
#ifdef FLUXCT
! interpolate cell centered components of the magnetic field
!
    call magtocen(w)
#endif /* FLUXCT */

! iterate over all points and calculate maximum speed
!
    do k = kb, ke
      do j = jb, je

        call cons2prim(im, w(1:nvr,1:im,j,k), q(1:nvr,1:im))

        do i = ib, ie

! calculate the velocity
!
          vv = sum(q(ivx:ivz,i)**2)
          v  = sqrt(vv)
#ifdef MHD
          bb = sum(q(icx:icz,i)**2)
#endif /* MHD */

! calculate the maximum characteristic speed
!
#ifdef MHD
#ifdef ADI
          c = sqrt((gamma * q(ipr,i) + bb) / q(idn,i))
#endif /* ADI */
#ifdef ISO
          c = sqrt(csnd2 + bb / q(idn,i))
#endif /* ISO */
#else /* MHD */
#ifdef ADI
          c = sqrt(gamma * q(ipr,i) / q(idn,i))
#endif /* ADI */
#ifdef ISO
          c = csnd
#endif /* ISO */
#endif /* MHD */

! calculate maximum of the speed
!
          maxspeed = max(maxspeed, v + c)
        end do
      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end function maxspeed

!===============================================================================
!
end module
