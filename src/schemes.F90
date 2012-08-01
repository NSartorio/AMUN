!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2012 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: SCHEMES
!!
!!  The module handles numerical schemes, subroutines to calculate variable
!!  increment and one dimensional Riemann solvers.
!!
!!******************************************************************************
!
module schemes

! module variables are not implicit by default
!
  implicit none

! by default everything is private
!
  private

! declare public subroutines
!
  public :: update_flux, update_increment

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
! subroutine UPDATE_FLUX:
! ----------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     idir - direction along which the flux is calculated;
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!
!===============================================================================
!
  subroutine update_flux(idir, dx, q, f)

! include external variables
!
    use coordinates   , only : im, jm, km
    use variables     , only : idn, ivx, ivy, ivz, imx, imy, imz
#ifdef ADI
    use variables     , only : ipr, ien
#endif /* ADI */
#ifdef MHD
    use variables     , only : ibx, iby, ibz
#ifdef GLM
    use variables     , only : iph
#endif /* GLM */
#endif /* MHD */
    use variables     , only : nt

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer                     , intent(in)    :: idir
    real                        , intent(in)    :: dx
    real, dimension(nt,im,jm,km), intent(in)    :: q
    real, dimension(nt,im,jm,km), intent(inout) :: f

! local variables
!
    integer                      :: i, j, k

! local temporary arrays
!
    real, dimension(nt,im)       :: qx, fx
    real, dimension(nt,jm)       :: qy, fy
#if NDIMS == 3
    real, dimension(nt,km)       :: qz, fz
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
! reset the flux array
!
    f(:,:,:,:) = 0.0d0

! select the directional flux to compute
!
    select case(idir)
    case(1)

!  calculate the flux along the X-direction
!
      do k = 1, km
        do j = 1, jm

! copy directional variable vectors to pass to the one dimensional solver
!
          qx(idn,1:im) = q(idn,1:im,j,k)
          qx(imx,1:im) = q(ivx,1:im,j,k)
          qx(imy,1:im) = q(ivy,1:im,j,k)
          qx(imz,1:im) = q(ivz,1:im,j,k)
#ifdef ADI
          qx(ien,1:im) = q(ien,1:im,j,k)
#endif /* ADI */
#ifdef MHD
          qx(ibx,1:im) = q(ibx,1:im,j,k)
          qx(iby,1:im) = q(iby,1:im,j,k)
          qx(ibz,1:im) = q(ibz,1:im,j,k)
#ifdef GLM
          qx(iph,1:im) = q(iph,1:im,j,k)
#endif /* GLM */
#endif /* MHD */

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(im, dx, qx(1:nt,1:im), fx(1:nt,1:im))

! update the array of fluxes
!
          f(idn,1:im,j,k) = fx(idn,1:im)
          f(imx,1:im,j,k) = fx(imx,1:im)
          f(imy,1:im,j,k) = fx(imy,1:im)
          f(imz,1:im,j,k) = fx(imz,1:im)
#ifdef ADI
          f(ien,1:im,j,k) = fx(ien,1:im)
#endif /* ADI */
#ifdef MHD
          f(ibx,1:im,j,k) = fx(ibx,1:im)
          f(iby,1:im,j,k) = fx(iby,1:im)
          f(ibz,1:im,j,k) = fx(ibz,1:im)
#ifdef GLM
          f(iph,1:im,j,k) = fx(iph,1:im)
#endif /* GLM */
#endif /* MHD */

        end do
      end do

    case(2)

!  calculate the flux along the Y direction
!
      do k = 1, km
        do i = 1, im

! copy directional variable vectors to pass to the one dimensional solver
!
          qy(idn,1:jm) = q(idn,i,1:jm,k)
          qy(ivx,1:jm) = q(ivy,i,1:jm,k)
          qy(ivy,1:jm) = q(ivz,i,1:jm,k)
          qy(ivz,1:jm) = q(ivx,i,1:jm,k)
#ifdef ADI
          qy(ien,1:jm) = q(ien,i,1:jm,k)
#endif /* ADI */
#ifdef MHD
          qy(ibx,1:jm) = q(iby,i,1:jm,k)
          qy(iby,1:jm) = q(ibz,i,1:jm,k)
          qy(ibz,1:jm) = q(ibx,i,1:jm,k)
#ifdef GLM
          qy(iph,1:jm) = q(iph,i,1:jm,k)
#endif /* GLM */
#endif /* MHD */

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(jm, dx, qy(1:nt,1:jm), fy(1:nt,1:jm))

! update the array of fluxes
!
          f(idn,i,1:jm,k) = fy(idn,1:jm)
          f(imx,i,1:jm,k) = fy(imz,1:jm)
          f(imy,i,1:jm,k) = fy(imx,1:jm)
          f(imz,i,1:jm,k) = fy(imy,1:jm)
#ifdef ADI
          f(ien,i,1:jm,k) = fy(ien,1:jm)
#endif /* ADI */
#ifdef MHD
          f(ibx,i,1:jm,k) = fy(ibz,1:jm)
          f(iby,i,1:jm,k) = fy(ibx,1:jm)
          f(ibz,i,1:jm,k) = fy(iby,1:jm)
#ifdef GLM
          f(iph,i,1:jm,k) = fy(iph,1:jm)
#endif /* GLM */
#endif /* MHD */

        end do
      end do

#if NDIMS == 3
    case(3)

!  calculate the flux along the Z direction
!
      do j = 1, jm
        do i = 1, im

! copy directional variable vectors to pass to the one dimensional solver
!
          qz(idn,1:km) = q(idn,i,j,1:km)
          qz(ivx,1:km) = q(ivz,i,j,1:km)
          qz(ivy,1:km) = q(ivx,i,j,1:km)
          qz(ivz,1:km) = q(ivy,i,j,1:km)
#ifdef ADI
          qz(ien,1:km) = q(ien,i,j,1:km)
#endif /* ADI */
#ifdef MHD
          qz(ibx,1:km) = q(ibz,i,j,1:km)
          qz(iby,1:km) = q(ibx,i,j,1:km)
          qz(ibz,1:km) = q(iby,i,j,1:km)
#ifdef GLM
          qz(iph,1:km) = q(iph,i,j,1:km)
#endif /* GLM */
#endif /* MHD */

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(km, dx, qz(1:nt,1:km), fz(1:nt,1:km))

! update the array of fluxes
!
          f(idn,i,j,1:km) = fz(idn,1:km)
          f(imx,i,j,1:km) = fz(imy,1:km)
          f(imy,i,j,1:km) = fz(imz,1:km)
          f(imz,i,j,1:km) = fz(imx,1:km)
#ifdef ADI
          f(ien,i,j,1:km) = fz(ien,1:km)
#endif /* ADI */
#ifdef MHD
          f(ibx,i,j,1:km) = fz(iby,1:km)
          f(iby,i,j,1:km) = fz(ibz,1:km)
          f(ibz,i,j,1:km) = fz(ibx,1:km)
#ifdef GLM
          f(iph,i,j,k) = fz(iph,k)
#endif /* GLM */
#endif /* MHD */

        end do
      end do
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine update_flux
!
!===============================================================================
!
! subroutine UPDATE_INCREMENT:
! ---------------------------
!
!   Subroutine calculate the conservative variable increment from the fluxes.
!
!   Arguments:
!
!     dh   - the ratio of the time step to the spatial step;
!     f    - the array of numerical fluxes;
!     du   - the array of variable increment;
!
!
!===============================================================================
!
  subroutine update_increment(dh, f, du)

! include external variables
!
    use coordinates   , only : im, jm, km
    use variables     , only : nfl

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real, dimension(3)                 , intent(in)    :: dh
    real, dimension(NDIMS,nfl,im,jm,km), intent(in)    :: f
    real, dimension(      nfl,im,jm,km), intent(inout) :: du

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! reset the increment array du
!
    du(:,:,:,:) = 0.0d0

! perform update along the X direction
!
    do i = 2, im
      du(:,i,:,:) = du(:,i,:,:) - dh(1) * (f(1,:,i,:,:) - f(1,:,i-1,:,:))
    end do
    du(:,1,:,:) = du(:,1,:,:) - dh(1) * f(1,:,1,:,:)

! perform update along the Y direction
!
    do j = 2, jm
      du(:,:,j,:) = du(:,:,j,:) - dh(2) * (f(2,:,:,j,:) - f(2,:,:,j-1,:))
    end do
    du(:,:,1,:) = du(:,:,1,:) - dh(2) * f(2,:,:,1,:)

#if NDIMS == 3
! perform update along the Z direction
!
    do k = 2, km
      du(:,:,:,k) = du(:,:,:,k) - dh(3) * (f(3,:,:,:,k) - f(3,:,:,:,k-1))
    end do
    du(:,:,:,1) = du(:,:,:,1) - dh(3) * f(3,:,:,:,1)
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine update_increment
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
#ifdef HYDRO
#ifdef HLL
!===============================================================================
!
! subroutine RIEMANN:
! ==================
!
!   Subroutine solves one dimensional Riemann problem using the HLL method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann(n, h, q, f)

! include external procedures
!
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct
#ifdef FIX_POSITIVITY
    use interpolations, only : fix_positivity
#endif /* FIX_POSITIVITY */

! include external variables
!
    use variables     , only : nt
    use variables     , only : ivx, ivy, ivz
#ifdef ADI
    use variables     , only : ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n
    real                 , intent(in)  :: h
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: f

! local variables
!
    integer               :: p, i
    real                  :: sl, sr, srl, srml

! local arrays to store the states
!
    real, dimension(nt,n) :: ql, qr, ul, ur, fl, fr
    real, dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states of primitive variables
!
    do p = 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
#ifdef ADI
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* ADI */
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else

! prepare coefficients for the intermediate state calculation
!
        srl  = sr * sl
        srml = sr - sl

! calculate the fluxes for the intermediate state
!
        f(:,i) = ((sr * fl(:,i) - sl * fr(:,i))                                &
                                           + srl * (ur(:,i) - ul(:,i))) / srml

      end if

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* HLL */
#ifdef HLLC
!
!===============================================================================
!
! subroutine RIEMANN:
! ------------------
!
!   Subroutine solves one dimensional Riemann problem using the HLLC method,
!   by Toro.  In the HLLC method the tangential components of the velocity are
!   discontinuous, which in the HLLCC method they are continuous and calculated
!   from the HLL average.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Toro, E. F., Spruce, M., & Speares, W.
!         "Restoration of the contact surface in the HLL-Riemann solver",
!         Shock Waves, 1994, Volume 4, Issue 1, pp. 25-34
!
!===============================================================================
!
  subroutine riemann(n, h, q, f)

! include external procedures
!
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct
#ifdef FIX_POSITIVITY
    use interpolations, only : fix_positivity
#endif /* FIX_POSITIVITY */

! include external variables
!
    use variables     , only : nt
    use variables     , only : idn, ivx, ivy, ivz, imx, imy, imz
#ifdef ADI
    use variables     , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n
    real                 , intent(in)  :: h
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: f

! local variables
!
    integer               :: p, i
    real                  :: dn, mx, pr
    real                  :: sl, sr, sm
    real                  :: slmv, srmv, slmm, srmm, smvl, smvr

! local arrays to store the states
!
    real, dimension(nt,n) :: ql, qr, ul, ur, fl, fr
    real, dimension(n)    :: cl, cr
    real, dimension(nt)   :: q1l, q1r, u1l, u1r
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states of primitive variables
!
    do p = 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLC flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! the speed of contact discontinuity (eq. 34 [1], 14 [2])
!
        dn = (sr * ur(idn,i) - sl * ul(idn,i)) - (fr(idn,i) - fl(idn,i))
        mx = (sr * ur(imx,i) - sl * ul(imx,i)) - (fr(imx,i) - fl(imx,i))
        sm = mx / dn

! compute speed differences
!
        slmv = sl - ql(ivx,i)
        srmv = sr - qr(ivx,i)
        slmm = sl - sm
        srmm = sr - sm
        smvl = sm - ql(ivx,i)
        smvr = sm - qr(ivx,i)

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
        pr = 0.5d0 * ((ql(ipr,i) + ql(idn,i) * slmv * smvl)                    &
                    + (qr(ipr,i) + qr(idn,i) * srmv * smvr))

! intermediate discontinuities
!
        if (sm > 0.0d0) then

! calculate the left intermediate state
!
          u1l(idn) = ql(idn,i) * slmv / slmm
          u1l(imx) = u1l(idn) * sm
          u1l(imy) = u1l(idn) * ql(ivy,i)
          u1l(imz) = u1l(idn) * ql(ivz,i)
          u1l(ien) = (slmv * ul(ien,i) - ql(ipr,i) * ql(ivx,i)                 &
                                                             + pr * sm) / slmm

! calculate the left intermediate flux
!
          f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

        else if (sm < 0.0d0) then

! calculate the right intermediate state
!
          u1r(idn) = qr(idn,i) * srmv / srmm
          u1r(imx) = u1r(idn) * sm
          u1r(imy) = u1r(idn) * qr(ivy,i)
          u1r(imz) = u1r(idn) * qr(ivz,i)
          u1r(ien) = (srmv * ur(ien,i) - qr(ipr,i) * qr(ivx,i)                 &
                                                             + pr * sm) / srmm

! calculate right intermediate flux
!
          f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

        else ! sm = 0

! calculate the left intermediate state
!
          u1l(idn) = ql(idn,i) * slmv / sl
          u1l(imx) = 0.0d0
          u1l(imy) = u1l(idn) * ql(ivy,i)
          u1l(imz) = u1l(idn) * ql(ivz,i)
          u1l(ien) = (slmv * ul(ien,i) - ql(ipr,i) * ql(ivx,i)) / sl

! calculate the left intermediate flux
!
          f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

! calculate the right intermediate state
!
          u1r(idn) = qr(idn,i) * srmv / sr
          u1r(imx) = 0.0d0
          u1r(imy) = u1r(idn) * qr(ivy,i)
          u1r(imz) = u1r(idn) * qr(ivz,i)
          u1r(ien) = (srmv * ur(ien,i) - qr(ipr,i) * qr(ivx,i)) / sr

! calculate right intermediate flux
!
          f(:,i) = 0.5d0 * (f(:,i) + (fr(:,i) + sr * (u1r(:) - ur(:,i))))

        end if ! sm = 0

      end if ! sl < 0.0 & sr > 0.0

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* HLLC */
#ifdef ROE
!
!===============================================================================
!
! subroutine RIEMANN:
! ------------------
!
!   Subroutine solves one dimensional Riemann problem using the HLLC method,
!   by Toro.  In the HLLC method the tangential components of the velocity are
!   discontinuous, which in the HLLCC method they are continuous and calculated
!   from the HLL average.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Roe, P. L., 1981, Journal of Computational Physics, 43, 357
!     [2] Toro, E. F., Spruce, M., & Speares, W.
!         "Restoration of the contact surface in the HLL-Riemann solver",
!         Shock Waves, 1994, Volume 4, Issue 1, pp. 25-34
!
!===============================================================================
!
  subroutine riemann(n, h, q, f)

! include external procedures
!
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct
#ifdef FIX_POSITIVITY
    use interpolations, only : fix_positivity
#endif /* FIX_POSITIVITY */

! include external variables
!
    use equations     , only : gamma
    use variables     , only : nt
    use variables     , only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables     , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n
    real                 , intent(in)  :: h
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: f

! local variables
!
    integer                :: p, i
    real                   :: al, ar, ap, div
    real                   :: sdl, sdr, sds, sfl, sfr

! local arrays to store the states
!
    real, dimension(nt,n)  :: ql, qr, ul, ur, fl, fr
    real, dimension(n)     :: cl, cr
    real, dimension(nt,nt) :: li, ri
    real, dimension(nt)    :: qi, ci, et, du
!
!-------------------------------------------------------------------------------
!
! reset the eigensystem values
!
    ci(:)   = 0.0d0
    li(:,:) = 0.0d0
    ri(:,:) = 0.0d0

! reconstruct the left and right states of primitive variables
!
    do p = 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! calculate conserved states difference
!
      du(:) = ur(:,i) - ul(:,i)

! calculate Roe variables for the eigenproblem solution
!
      sdl     = sqrt(ql(idn,i))
      sdr     = sqrt(qr(idn,i))
      sds     = sdl + sdr
      sfl     = sdl / sds
      sfr     = sdr / sds

! prepare the Roe intermediate state
!
      qi(idn) = sdl * sdr
      qi(ivx) = sfl * ql(ivx,i) + sfr * qr(ivx,i)
      qi(ivy) = sfl * ql(ivy,i) + sfr * qr(ivy,i)
      qi(ivz) = sfl * ql(ivz,i) + sfr * qr(ivz,i)
#ifdef ADI
      qi(ipr) = sfl * (ul(ien,i) + ql(ipr,i)) / ql(idn,i)                      &
              + sfr * (ur(ien,i) + qr(ipr,i)) / qr(idn,i)
#endif /* ADI */

! obtain eigenvalues and eigenvectors
!
      call eigensystem(qi(:), ci(:), ri(:,:), li(:,:))

! calculate vector (Ur - Ul).L
!
      et(:) = 0.0d0
      do p = 1, nt
        et(:) = et(:) + du(p) * li(p,:)
      end do

! calculate the flux average
!
      f(:,i) = 0.5d0 * (fl(:,i) + fr(:,i))

! correct the flux by adding the term abs(lambda) R.(Ur - Ul).L
!
      do p = 1, nt
        f(:,i) = f(:,i) - 0.5d0 * abs(ci(p)) * et(p) * ri(p,:)
      end do

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
!
!===============================================================================
!
! eigensystem: subroutine computes eigenvalues and eigenmatrices for a given
!              set of equations and input variables
!
!===============================================================================
!
#ifdef ADI
  subroutine eigensystem(q, c, r, l)

    use equations, only : gamma
    use variables, only : nqt
    use variables, only : idn, ivx, ivy, ivz
    use variables, only : ien

    implicit none

! input/output arguments
!
    real, dimension(nqt)    , intent(in)    :: q
    real, dimension(nqt)    , intent(inout) :: c
    real, dimension(nqt,nqt), intent(inout) :: l, r

! local variables
!
    real :: gm, vv, vh, c2, na, cc, vc, ng, nd, nv, nh, nc
!
!-------------------------------------------------------------------------------
!
! calculate characteristic speeds and useful variables
!
    gm = gamma - 1.0d0
    vv = sum(q(ivx:ivz)**2)
    vh = 0.5d0 * vv
    c2 = gm * (q(ien) - vh)
    na = 0.5d0 / c2
    cc = sqrt(c2)
    vc = q(ivx) * cc
    ng = na * gm
    nd = 2.0 * ng
    nv = na * vc
    nh = na * gm * vh
    nc = na * cc

! prepare eigenvalues
!
    c(1) = q(ivx) - cc
    c(2) = q(ivx)
    c(3) = q(ivx)
    c(4) = q(ivx)
    c(5) = q(ivx) + cc

! prepare the right eigenmatrix
!
    r(1,idn) = 1.0d0
    r(1,ivx) = q(ivx) - cc
    r(1,ivy) = q(ivy)
    r(1,ivz) = q(ivz)
    r(1,ien) = q(ien) - vc

    r(2,ivy) = 1.0d0
    r(2,ien) = q(ivy)

    r(3,ivz) = 1.0d0
    r(3,ien) = q(ivz)

    r(4,idn) = 1.0d0
    r(4,ivx) = q(ivx)
    r(4,ivy) = q(ivy)
    r(4,ivz) = q(ivz)
    r(4,ien) = vh

    r(5,idn) = 1.0d0
    r(5,ivx) = q(ivx) + cc
    r(5,ivy) = q(ivy)
    r(5,ivz) = q(ivz)
    r(5,ien) = q(ien) + vc

! prepare the left eigenmatrix
!
    l(idn,1) =   nh + nv
    l(ivx,1) = - ng * q(ivx) - nc
    l(ivy,1) = - ng * q(ivy)
    l(ivz,1) = - ng * q(ivz)
    l(ien,1) =   ng

    l(idn,2) = - q(ivy)
    l(ivy,2) = 1.0d0

    l(idn,3) = - q(ivz)
    l(ivz,3) = 1.0d0

    l(idn,4) = 1.0d0 - ng * vv
    l(ivx,4) =   nd * q(ivx)
    l(ivy,4) =   nd * q(ivy)
    l(ivz,4) =   nd * q(ivz)
    l(ien,4) = - nd

    l(idn,5) =   nh - nv
    l(ivx,5) = - ng * q(ivx) + nc
    l(ivy,5) = - ng * q(ivy)
    l(ivz,5) = - ng * q(ivz)
    l(ien,5) =   ng

!-------------------------------------------------------------------------------
!
  end subroutine eigensystem
#endif /* ADI */
#ifdef ISO
  subroutine eigensystem(q, c, r, l)

    use equations, only : csnd
    use variables, only : nqt
    use variables, only : idn, ivx, ivy, ivz

    implicit none

! input/output arguments
!
    real, dimension(nqt)    , intent(in)    :: q
    real, dimension(nqt)    , intent(inout) :: c
    real, dimension(nqt,nqt), intent(inout) :: l, r

! local variables
!
    real :: ch, vc
!
!-------------------------------------------------------------------------------
!
! calculate useful variables
!
    ch = 0.5d0 / csnd
    vc = ch * q(ivx)

! prepare eigenvalues
!
    c(1) = q(ivx) - csnd
    c(2) = q(ivx)
    c(3) = q(ivx)
    c(4) = q(ivx) + csnd

! prepare the right eigenmatrix
!
    r(1,idn) = 1.0d0
    r(1,ivx) = q(ivx) - csnd
    r(1,ivy) = q(ivy)
    r(1,ivz) = q(ivz)

    r(2,ivy) = 1.0d0

    r(3,ivz) = 1.0d0

    r(4,idn) = 1.0d0
    r(4,ivx) = q(ivx) + csnd
    r(4,ivy) = q(ivy)
    r(4,ivz) = q(ivz)

! prepare the left eigenmatrix
!
    l(idn,1) = 0.5d0 + vc
    l(ivx,1) = - ch

    l(idn,2) = - q(ivy)
    l(ivy,2) = 1.0d0

    l(idn,3) = - q(ivz)
    l(ivz,3) = 1.0d0

    l(idn,4) = 0.5d0 - vc
    l(ivx,4) =   ch

!-------------------------------------------------------------------------------
!
  end subroutine eigensystem
#endif /* ISO */
#endif /* ROE */
#endif /* HYDRO */
#ifdef MHD
#ifdef HLL
!===============================================================================
!
! subroutine RIEMANN:
! ==================
!
!   Subroutine solves one dimensional Riemann problem using the HLL method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann(n, h, q, f)

! include external procedures
!
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct
#ifdef FIX_POSITIVITY
    use interpolations, only : fix_positivity
#endif /* FIX_POSITIVITY */

! include external variables
!
    use variables     , only : nt
    use variables     , only : ivx, ivy, ivz
#ifdef ADI
    use variables     , only : ien
#endif /* ADI */
    use variables     , only : ibx, iby, ibz
#ifdef GLM
    use variables     , only : iph
    use variables     , only : cmax
#endif /* GLM */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer              , intent(in)  :: n
    real                 , intent(in)  :: h
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: f

! local variables
!
    integer               :: p, i
    real                  :: sl, sr, srl, srml

! local arrays to store the states
!
    real, dimension(nt,n) :: ql, qr, ul, ur, fl, fr
    real, dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states of primitive variables
!
    do p = 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! reconstruct the left and right states of the magnetic field components
!
    do p = ibx, ibz
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef GLM
! reconstruct the left and right states of the scalar potential
!
    call reconstruct(n, h, q(iph,:), ql(iph,:), qr(iph,:))

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d0 * ((qr(ibx,:) + ql(ibx,:)) - (qr(iph,:) - ql(iph,:)) / cmax)
    cr(:) = 0.5d0 * ((qr(iph,:) + ql(iph,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(iph,:) = cr(:)
    qr(iph,:) = cr(:)
#endif /* GLM */

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
#ifdef ADI
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* ADI */
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else

! prepare coefficients for the intermediate state calculation
!
        srl  = sr * sl
        srml = sr - sl

! calculate the fluxes for the intermediate state
!
        f(:,i) = ((sr * fl(:,i) - sl * fr(:,i))                                &
                                           + srl * (ur(:,i) - ul(:,i))) / srml

      end if

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* HLL */
#if defined HLLC || defined HLLCC || defined HLLCL
!
!===============================================================================
!
! subroutine RIEMANN:
! ==================
!
!   Subroutine solves one dimensional Riemann problem using the HLLC method
!   by Gurski or Li.  The HLLC, HLLCC, and HLLC-L differs by the definitions of
!   the tangential components of the velocity and magnetic field.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     b - the input vector of the normal magnetic field component;
!     f - the output array of fluxes;
!     s - the input array of shock indicators;
!
!   References:
!
!     [1] Gurski, K. F.,
!         "An HLLC-Type Approximate Riemann Solver for Ideal
!          Magnetohydrodynamics",
!         SIAM Journal on Scientific Computing, 2004, Volume 25, Issue 6,
!         pp. 2165–2187
!     [2] Li, S.,
!         "An HLLC Riemann solver for magneto-hydrodynamics",
!         Journal of Computational Physics, 2005, Volume 203, Issue 1,
!         pp. 344-357
!
!===============================================================================
!
  subroutine riemann(n, h, q, b, f, s)

! include external procedures and variables
!
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity
    use variables     , only : nt
    use variables     , only : idn, imx, imy, imz, ien
    use variables     , only : ivx, ivy, ivz, ipr
    use variables     , only : ibx, iby, ibz

! local variables are not implicit by default
!
    implicit none

! input / output arguments
!
    integer                 , intent(in)  :: n
    real                    , intent(in)  :: h
    real   , dimension(nt,n), intent(in)  :: q
    real   , dimension(   n), intent(in)  :: b
    real   , dimension(nt,n), intent(out) :: f
    logical, dimension(n)   , intent(in)  :: s

! local variables
!
    integer :: p, i
#ifdef SHOCK_DETECTION
    integer :: ip1
#endif /* SHOCK_DETECTION */
    real    :: dn, mx, vy, vz, by, bz, pt, pm
    real    :: sl, sr, srl, sm
    real    :: slmv, srmv, slmm, srmm, smvl, smvr, srml
    real    :: pml, pmr, ptl, ptr
    real    :: fc

! local arrays to store the states
!
    real, dimension(nt,n) :: ul, ur, ql, qr, fl, fr
    real, dimension(n)    :: cl, cr
    real, dimension(nt)   :: q1l, q1r, u1l, u1r
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states of primitive variables
!
    do p = 1, ibx - 1
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do
    do p = ibx + 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do

! copy normal component to the left and right states
!
    ql(ibx,:) = b(:)
    qr(ibx,:) = b(:)

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLC flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

#ifdef SHOCK_DETECTION
! calculate the next index
!
        ip1 = min(n, i + 1)

! check if we are in the vicinity of a shock
!
        if (s(i) .or. s(ip1)) then

! prepare coefficients for the intermediate state calculation
!
          srl  = sr * sl
          srml = sr - sl

! calculate the fluxes for the intermediate state
!
          f(:,i) = ((sr * fl(:,i) - sl * fr(:,i))                              &
                                           + srl * (ur(:,i) - ul(:,i))) / srml

        else
#endif /* SHOCK_DETECTION */

! the speed of contact discontinuity (eq. 34 [1], 14 [2])
!
          dn = (sr * ur(idn,i) - sl * ul(idn,i)) - (fr(idn,i) - fl(idn,i))
          mx = (sr * ur(imx,i) - sl * ul(imx,i)) - (fr(imx,i) - fl(imx,i))
          sm = mx / dn

! compute speed differences
!
          slmv = sl - ql(ivx,i)
          srmv = sr - qr(ivx,i)
          slmm = sl - sm
          srmm = sr - sm
          smvl = sm - ql(ivx,i)
          smvr = sm - qr(ivx,i)

! calculate magnetic pressure
!
          pml = 0.5d0 * sum(ql(ibx:ibz,i) * ql(ibx:ibz,i))
          pmr = 0.5d0 * sum(qr(ibx:ibz,i) * qr(ibx:ibz,i))

! calculate total pressures
!
          ptl = ql(ipr,i) + pml
          ptr = qr(ipr,i) + pmr

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
          pt = 0.5d0 * ((ptl + ql(idn,i) * slmv * smvl)                        &
                      + (ptr + qr(idn,i) * srmv * smvr))

#if defined HLLC || defined HLLCL
! separate cases when Bx = 0 and Bx ≠ 0
!
          if (b(i) == 0.0d0) then

! intermediate discontinuities
!
            if (sm > 0.0d0) then

! calculate the ratio of speed differences
!
              fc = slmv / slmm

! calculate density
!
              dn = ql(idn,i) * fc

! calculate the left intermediate state
!
              u1l(idn) = dn
              u1l(imx) = dn * sm
              u1l(imy) = dn * ql(ivy,i)
              u1l(imz) = dn * ql(ivz,i)
              u1l(ibx) = b(i)
              u1l(iby) = ql(iby,i) * fc
              u1l(ibz) = ql(ibz,i) * fc
              u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i) + pt * sm) / slmm

! calculate the left intermediate flux
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

            else if (sm < 0.0d0) then

! calculate the ratio of speed differences
!
              fc = srmv / srmm

! calculate density
!
              dn = qr(idn,i) * fc

! calculate the right intermediate state
!
              u1r(idn) = dn
              u1r(imx) = dn * sm
              u1r(imy) = dn * qr(ivy,i)
              u1r(imz) = dn * qr(ivz,i)
              u1r(ibx) = b(i)
              u1r(iby) = qr(iby,i) * fc
              u1r(ibz) = qr(ibz,i) * fc
              u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i) + pt * sm) / srmm

! calculate right intermediate flux
!
              f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

            else ! sm = 0

! calculate the ratio of speed differences
!
              fc = slmv / sl

! calculate density
!
              dn = ql(idn,i) * fc

! calculate the left intermediate state
!
              u1l(idn) = dn
              u1l(imx) = 0.0d0
              u1l(imy) = dn * ql(ivy,i)
              u1l(imz) = dn * ql(ivz,i)
              u1l(ibx) = b(i)
              u1l(iby) = ql(iby,i) * fc
              u1l(ibz) = ql(ibz,i) * fc
              u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i)) / sl

! calculate the left intermediate flux
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

! calculate the ratio of speed differences
!
              fc = srmv / sr

! calculate density
!
              dn = qr(idn,i) * fc

! calculate the right intermediate state
!
              u1r(idn) = dn
              u1r(imx) = 0.0d0
              u1r(imy) = dn * qr(ivy,i)
              u1r(imz) = dn * qr(ivz,i)
              u1r(ibx) = b(i)
              u1r(iby) = qr(iby,i) * fc
              u1r(ibz) = qr(ibz,i) * fc
              u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i)) / sr

! calculate right intermediate flux
!
              f(:,i) = 0.5d0 * (f(:,i) + (fr(:,i) + sr * (u1r(:) - ur(:,i))))

            end if ! sm = 0

          else ! Bx ≠ 0

! prepare the difference between the fastests speeds
!
            srml = sr - sl

! when Bx ≠ 0 the tangential components of velocity and magnetic field are
! calculated from the HLL state
!
#ifdef HLLC
            vy = ((sr * ur(imy,i) - sl * ul(imy,i))                            &
                                             - (fr(imy,i) - fl(imy,i))) / dn
            vz = ((sr * ur(imz,i) - sl * ul(imz,i))                            &
                                             - (fr(imz,i) - fl(imz,i))) / dn
#endif /* HLLC */
            by = ((sr * ur(iby,i) - sl * ul(iby,i))                            &
                                             - (fr(iby,i) - fl(iby,i))) / srml
            bz = ((sr * ur(ibz,i) - sl * ul(ibz,i))                            &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml

! intermediate discontinuities
!
            if (sm > 0.0d0) then

! calculate the ratio of speed differences
!
              fc = slmv / slmm

! calculate density
!
              dn = ql(idn,i) * fc

! calculate the left intermediate state
!
              u1l(idn) = dn
              u1l(imx) = dn * sm
#ifdef HLLC
              u1l(imy) = dn * vy
              u1l(imz) = dn * vz
#endif /* HLLC */
#ifdef HLLCL
              fc = b(i) / (ql(idn,i) * slmv)
              u1l(imy) = dn * (ql(ivy,i) + fc * (ql(iby,i) - by))
              u1l(imz) = dn * (ql(ivz,i) + fc * (ql(ibz,i) - bz))
#endif /* HLLCL */
              u1l(ibx) = b(i)
              u1l(iby) = by
              u1l(ibz) = bz
              u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i) + pt * sm) / slmm

! calculate the left intermediate flux
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

            else if (sm < 0.0d0) then

! calculate the ratio of speed differences
!
              fc = srmv / srmm

! calculate density
!
              dn = qr(idn,i) * fc

! calculate the right intermediate state
!
              u1r(idn) = dn
              u1r(imx) = dn * sm
#ifdef HLLC
              u1r(imy) = dn * vy
              u1r(imz) = dn * vz
#endif /* HLLC */
#ifdef HLLCL
              fc = b(i) / (qr(idn,i) * srmv)
              u1r(imy) = dn * (qr(ivy,i) + fc * (qr(iby,i) - by))
              u1r(imz) = dn * (qr(ivz,i) + fc * (qr(ibz,i) - bz))
#endif /* HLLCL */
              u1r(ibx) = b(i)
              u1r(iby) = by
              u1r(ibz) = bz
              u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i) + pt * sm) / srmm

! calculate right intermediate flux
!
              f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

            else ! sm = 0

! calculate the ratio of speed differences
!
              fc = slmv / sl

! calculate density
!
              dn = ql(idn,i) * fc

! calculate the left intermediate state
!
              u1l(idn) = dn
              u1l(imx) = 0.0d0
#ifdef HLLC
              u1l(imy) = dn * vy
              u1l(imz) = dn * vz
#endif /* HLLC */
#ifdef HLLCL
              fc = b(i) / (ql(idn,i) * slmv)
              u1l(imy) = dn * (ql(ivy,i) + fc * (ql(iby,i) - by))
              u1l(imz) = dn * (ql(ivz,i) + fc * (ql(ibz,i) - bz))
#endif /* HLLCL */
              u1l(ibx) = b(i)
              u1l(iby) = by
              u1l(ibz) = bz
              u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i)) / sl

! calculate the left intermediate flux
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

! calculate the ratio of speed differences
!
              fc = srmv / sr

! calculate density
!
              dn = qr(idn,i) * fc

! calculate the right intermediate state
!
              u1r(idn) = dn
              u1r(imx) = 0.0d0
#ifdef HLLC
              u1r(imy) = dn * vy
              u1r(imz) = dn * vz
#endif /* HLLC */
#ifdef HLLCL
              fc = b(i) / (qr(idn,i) * srmv)
              u1r(imy) = dn * (qr(ivy,i) + fc * (qr(iby,i) - by))
              u1r(imz) = dn * (qr(ivz,i) + fc * (qr(ibz,i) - bz))
#endif /* HLLCL */
              u1r(ibx) = b(i)
              u1r(iby) = by
              u1r(ibz) = bz
              u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i)) / sr

! calculate right intermediate flux
!
              f(:,i) = 0.5d0 * (f(:,i) + (fr(:,i) + sr * (u1r(:) - ur(:,i))))

            end if ! sm = 0

          end if ! Bx ≠ 0
#endif /* HLLC | HLLCL */
#ifdef HLLCC
! prepare the difference between the fastests speeds
!
          srml = sr - sl

! when Bx ≠ 0 the tangential components of velocity and magnetic field are
! calculated from the HLL state
!
          vy = ((sr * ur(imy,i) - sl * ul(imy,i))                              &
                                             - (fr(imy,i) - fl(imy,i))) / dn
          vz = ((sr * ur(imz,i) - sl * ul(imz,i))                              &
                                             - (fr(imz,i) - fl(imz,i))) / dn
          by = ((sr * ur(iby,i) - sl * ul(iby,i))                              &
                                             - (fr(iby,i) - fl(iby,i))) / srml
          bz = ((sr * ur(ibz,i) - sl * ul(ibz,i))                              &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml

! intermediate discontinuities
!
          if (sm > 0.0d0) then

! calculate the ratio of speed differences
!
            fc = slmv / slmm

! calculate density
!
            dn = ql(idn,i) * fc

! calculate the left intermediate state
!
            u1l(idn) = dn
            u1l(imx) = dn * sm
            u1l(imy) = dn * vy
            u1l(imz) = dn * vz
            u1l(ibx) = b(i)
            u1l(iby) = by
            u1l(ibz) = bz
            u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i) + pt * sm) / slmm

! calculate the left intermediate flux
!
            f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

          else if (sm < 0.0d0) then

! calculate the ratio of speed differences
!
            fc = srmv / srmm

! calculate density
!
            dn = qr(idn,i) * fc

! calculate the right intermediate state
!
            u1r(idn) = dn
            u1r(imx) = dn * sm
            u1r(imy) = dn * vy
            u1r(imz) = dn * vz
            u1r(ibx) = b(i)
            u1r(iby) = by
            u1r(ibz) = bz
            u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i) + pt * sm) / srmm

! calculate right intermediate flux
!
            f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

          else ! sm = 0

! calculate the ratio of speed differences
!
            fc = slmv / sl

! calculate density
!
            dn = ql(idn,i) * fc

! calculate the left intermediate state
!
            u1l(idn) = dn
            u1l(imx) = 0.0d0
            u1l(imy) = dn * vy
            u1l(imz) = dn * vz
            u1l(ibx) = b(i)
            u1l(iby) = by
            u1l(ibz) = bz
            u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i)) / sl

! calculate the left intermediate flux
!
            f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

! calculate the ratio of speed differences
!
            fc = srmv / sr

! calculate density
!
            dn = qr(idn,i) * fc

! calculate the right intermediate state
!
            u1r(idn) = dn
            u1r(imx) = 0.0d0
            u1r(imy) = dn * vy
            u1r(imz) = dn * vz
            u1r(ibx) = b(i)
            u1r(iby) = by
            u1r(ibz) = bz
            u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i)) / sr

! calculate right intermediate flux
!
            f(:,i) = 0.5d0 * (f(:,i) + (fr(:,i) + sr * (u1r(:) - ur(:,i))))

          end if ! sm = 0
#endif /* HLLCC */

#ifdef SHOCK_DETECTION
        end if ! shock vicinity
#endif /* SHOCK_DETECTION */

      end if ! sl < 0.0 & sr > 0.0

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* HLLC | HLLCC | HLLCL */
#ifdef HLLD
#ifdef ADI
!
!===============================================================================
!
! subroutine RIEMANN:
! ------------------
!
!   Subroutine solves one dimensional Riemann problem using the adiabatic HLLD
!   method described by Miyoshi & Kusano.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     b - the input vector of the normal magnetic field component;
!     f - the output array of fluxes;
!     s - the input array of shock indicators;
!
!   References:
!
!     [1] Miyoshi, T. & Kusano, K.,
!         "A multi-state HLL approximate Riemann solver for ideal
!          magnetohydrodynamics",
!         Journal of Computational Physics, 2005, 208, pp. 315-344
!
!===============================================================================
!
  subroutine riemann(n, h, q, b, f, s)

! include external procedures and variables
!
    use equations     , only : gamma
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity
    use variables     , only : nt
    use variables     , only : idn, imx, imy, imz, ien
    use variables     , only : ivx, ivy, ivz, ipr
    use variables     , only : ibx, iby, ibz

! local variables are not implicit by default
!
    implicit none

! input / output arguments
!
    integer                 , intent(in)  :: n
    real                    , intent(in)  :: h
    real   , dimension(nt,n), intent(in)  :: q
    real   , dimension(   n), intent(in)  :: b
    real   , dimension(nt,n), intent(out) :: f
    logical, dimension(n)   , intent(in)  :: s

! local variables
!
    integer :: p, i
#ifdef SHOCK_DETECTION
    integer :: ip1
#endif /* SHOCK_DETECTION */
    real    :: sl, sr, sm, sml, smr
    real    :: srl, srml, slmv, srmv, slmm, srmm, smvl, smvr
    real    :: dn, dnl, dnr, dlsq, drsq
    real    :: mx, dv, fc, b2, bs, ds, vbl, vbr, vb1l, vb1r, vb2
    real    :: pml, pmr, ptl, ptr, pt, pm

! local arrays to store the states
!
    real, dimension(nt,n) :: ul, ur, ql, qr, fl, fr
    real, dimension(n)    :: cl, cr
    real, dimension(nt)   :: q1l, q1r, u1l, u1r, q2, u2
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states of primitive variables
!
    do p = 1, ibx - 1
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do
    do p = ibx + 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do

! copy normal component to the left and right states
!
    ql(ibx,:) = b(:)
    qr(ibx,:) = b(:)

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))
#endif /* FIX_POSITIVITY */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLD flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

#ifdef SHOCK_DETECTION
! calculate the next index
!
        ip1 = min(n, i + 1)

! check if we are in the vicinity of a shock
!
        if (s(i) .or. s(ip1)) then

! prepare coefficients for the intermediate state calculation
!
          srl  = sr * sl
          srml = sr - sl

! calculate the fluxes for the intermediate state
!
          f(:,i) = ((sr * fl(:,i) - sl * fr(:,i))                              &
                                           + srl * (ur(:,i) - ul(:,i))) / srml

        else
#endif /* SHOCK_DETECTION */

! the speed of contact discontinuity (eq. 38)
!
          dn = (sr * ur(idn,i) - sl * ul(idn,i)) - (fr(idn,i) - fl(idn,i))
          mx = (sr * ur(imx,i) - sl * ul(imx,i)) - (fr(imx,i) - fl(imx,i))
          sm = mx / dn

! prepare speed differences
!
          slmv = sl - ql(ivx,i)
          srmv = sr - qr(ivx,i)
          slmm = sl - sm
          srmm = sr - sm
          smvl = sm - ql(ivx,i)
          smvr = sm - qr(ivx,i)
          srml = sr - sl

! prepare some coefficients
!
          dnl  = ql(idn,i) * slmv
          dnr  = qr(idn,i) * srmv
          b2   = b(i) * b(i)
          bs   = b2 / gamma

! calculate magnetic pressure
!
          pml = 0.5d0 * sum(ql(ibx:ibz,i) * ql(ibx:ibz,i))
          pmr = 0.5d0 * sum(qr(ibx:ibz,i) * qr(ibx:ibz,i))

! calculate total pressures
!
          ptl = ql(ipr,i) + pml
          ptr = qr(ipr,i) + pmr

! the pressure of the intermediate state (eq. 23)
!
          pt = 0.5d0 * ((ptl + dnl * smvl) + (ptr + dnr * smvr))

! calculate densities for the left and right intermediate states (eq. 43)
!
          q1l(idn) = dnl / slmm
          q1r(idn) = dnr / srmm

! substitute the normal components of velocity and magnetic field
!
          q1l(ivx) = sm
          q1r(ivx) = sm
          q1l(ibx) = b(i)
          q1r(ibx) = b(i)

! calculate transversal components of  magnetic field for the left state
! (eq. 45 & 47); take into account the B² > γp degeneracy
!
          dv = dnl * slmm - b2
          if (dv /= 0.0d0 .and. bs <= ql(ipr,i)) then

            fc       = (dnl * slmv - b2) / dv

            q1l(iby) = ql(iby,i) * fc
            q1l(ibz) = ql(ibz,i) * fc

          else

! in the case of degeneracy, calculate transversal components of magnetic field
! from the HLL states
!
            q1l(iby) = ((sr * ur(iby,i) - sl * ul(iby,i))                      &
                                             - (fr(iby,i) - fl(iby,i))) / srml
            q1l(ibz) = ((sr * ur(ibz,i) - sl * ul(ibz,i))                      &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml
          end if

! calculate transversal components of velocity (eq. 42)
!
          dv = b(i) / dnl
          q1l(ivy) = ql(ivy,i) + dv * (ql(iby,i) - q1l(iby))
          q1l(ivz) = ql(ivz,i) + dv * (ql(ibz,i) - q1l(ibz))

! calculate transversal components of  magnetic field for the left state
! (eq. 45 & 47); take into account the B² > γp degeneracy
!
          dv = dnr * srmm - b2
          if (dv /= 0.0d0 .and. bs <= qr(ipr,i)) then

            fc       = (dnr * srmv - b2) / dv

            q1r(iby) = qr(iby,i) * fc
            q1r(ibz) = qr(ibz,i) * fc

          else

! in the case of degeneracy, calculate transversal components of magnetic field
! from the HLL states
!
            q1r(iby) = ((sr * ur(iby,i) - sl * ul(iby,i))                      &
                                             - (fr(iby,i) - fl(iby,i))) / srml
            q1r(ibz) = ((sr * ur(ibz,i) - sl * ul(ibz,i))                      &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml

          end if

! calculate transversal components of velocity (eq. 42)
!
          dv = b(i) / dnr
          q1r(ivy) = qr(ivy,i) + dv * (qr(iby,i) - q1r(iby))
          q1r(ivz) = qr(ivz,i) + dv * (qr(ibz,i) - q1r(ibz))

! calculate scalar products of the velocity and magnetic field
!
          vbl  = sum(ql (ivx:ivz,i) * ql (ibx:ibz,i))
          vbr  = sum(qr (ivx:ivz,i) * qr (ibx:ibz,i))
          vb1l = sum(q1l(ivx:ivz)   * q1l(ibx:ibz))
          vb1r = sum(q1r(ivx:ivz)   * q1r(ibx:ibz))

! convert the left intermediate state to conservative form
!
          u1l(idn) = q1l(idn)
          u1l(imx) = q1l(idn) * q1l(ivx)
          u1l(imy) = q1l(idn) * q1l(ivy)
          u1l(imz) = q1l(idn) * q1l(ivz)
          u1l(ibx) = q1l(ibx)
          u1l(iby) = q1l(iby)
          u1l(ibz) = q1l(ibz)

! convert the right intermediate state to conservative form
!
          u1r(idn) = q1r(idn)
          u1r(imx) = q1r(idn) * q1r(ivx)
          u1r(imy) = q1r(idn) * q1r(ivy)
          u1r(imz) = q1r(idn) * q1r(ivz)
          u1r(ibx) = q1r(ibx)
          u1r(iby) = q1r(iby)
          u1r(ibz) = q1r(ibz)

! calculate the total energy of the left intermediate state (eq. 48)
!
          u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i) + pt * sm             &
                                                 + b(i) * (vbl - vb1l)) / slmm

! calculate the total energy of the right intermediate state (eq. 48)
!
          u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i) + pt * sm             &
                                                 + b(i) * (vbr - vb1r)) / srmm

! separate cases Bx ≠ 0 and Bx = 0
!
          if (abs(b(i)) > 0.0d0) then

! calculate the square root of the left and right intermediate densities
!
            dlsq = sqrt(q1l(idn))
            drsq = sqrt(q1r(idn))

! calculate the left and right going sound plus Alfvén speeds (eq. 51)
!
            sml = sm - abs(b(i)) / dlsq
            smr = sm + abs(b(i)) / drsq

! intermediate discontinuities
!
            if (sml > 0.0d0) then

! calculate the left intermediate flux, eq. (64)
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

            else if (smr < 0.0d0) then

! calculate the right intermediate flux, eq. (64)
!
              f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

            else ! sml ≤ 0 ≤ smr

! prepare the sign of the normal component of magnetic field
!
              if (b(i) >= 0.0d0) then
                bs =   1.0d0
              else
                bs = - 1.0d0
              end if

! calculate the sum and product of density root squares
!
              dv   = dlsq + drsq
              ds   = dlsq * drsq

! calculate velocity components, (eq. 39, and 59, 60)
!
              q2(ivx) = sm
              q2(ivy) = ((dlsq * q1l(ivy) + drsq * q1r(ivy))                   &
                                            + bs * (q1r(iby) - q1l(iby))) / dv
              q2(ivz) = ((dlsq * q1l(ivz) + drsq * q1r(ivz))                   &
                                            + bs * (q1r(ibz) - q1l(ibz))) / dv

! calculate magnetic field components, (eq. 61, 62)
!
              q2(ibx) = b(i)
              q2(iby) = ((dlsq * q1r(iby) + drsq * q1l(iby))                   &
                                       + bs * ds * (q1r(ivy) - q1l(ivy))) / dv
              q2(ibz) = ((dlsq * q1r(ibz) + drsq * q1l(ibz))                   &
                                       + bs * ds * (q1r(ivz) - q1l(ivz))) / dv

! calculate the scalar product of velocity and magnetic field
!
              vb2 = sum(q2(ivx:ivz) * q2(ibx:ibz))

! depending on the contact discontinuity speed chose the right Alfvén state
!
              if (sm > 0.0d0) then

! prepare the conservative variables for the left Alfvén intermediate state
!
                u2(idn) = u1l(idn)
                u2(imx) = u2(idn) * q2(ivx)
                u2(imy) = u2(idn) * q2(ivy)
                u2(imz) = u2(idn) * q2(ivz)
                u2(ibx) = q2(ibx)
                u2(iby) = q2(iby)
                u2(ibz) = q2(ibz)

! calculate the energy for the Alfvén intermediate state (eq. 63)
!
                u2(ien) = u1l(ien) - bs * dlsq * (vb1l - vb2)

! calculate the left Alfvén intermediate flux (eq. 65)
!
                f(:,i) = fl(:,i)                                               &
                            + sml * u2(:) - (sml - sl) * u1l(:) - sl * ul(:,i)

              else if (sm < 0.0d0) then

! prepare the conservative variables for the right Alfvén intermediate state
!
                u2(idn) = u1r(idn)
                u2(imx) = u2(idn) * q2(ivx)
                u2(imy) = u2(idn) * q2(ivy)
                u2(imz) = u2(idn) * q2(ivz)
                u2(ibx) = q2(ibx)
                u2(iby) = q2(iby)
                u2(ibz) = q2(ibz)

! calculate the energy for the Alfvén intermediate state (eq. 63)
!
                u2(ien) = u1r(ien) + bs * drsq * (vb1r - vb2)

! calculate the right Alfvén intermediate flux (eq. 65)
!
                f(:,i) = fr(:,i)                                               &
                            + smr * u2(:) - (smr - sr) * u1r(:) - sr * ur(:,i)

              else ! sm = 0

! prepare the conservative variables for the left Alfvén intermediate state
!
                u2(idn) = u1l(idn)
                u2(imx) = u2(idn) * q2(ivx)
                u2(imy) = u2(idn) * q2(ivy)
                u2(imz) = u2(idn) * q2(ivz)
                u2(ibx) = q2(ibx)
                u2(iby) = q2(iby)
                u2(ibz) = q2(ibz)

! calculate the energy for the Alfvén intermediate state (eq. 63)
!
                u2(ien) = u1l(ien) - bs * dlsq * (vb1l - vb2)

! calculate the left Alfvén intermediate flux (eq. 65)
!
                f(:,i) = fl(:,i)                                               &
                            + sml * u2(:) - (sml - sl) * u1l(:) - sl * ul(:,i)

! prepare the conservative variables for the right Alfvén intermediate state
!
                u2(idn) = u1r(idn)
                u2(imx) = u2(idn) * q2(ivx)
                u2(imy) = u2(idn) * q2(ivy)
                u2(imz) = u2(idn) * q2(ivz)
                u2(ibx) = q2(ibx)
                u2(iby) = q2(iby)
                u2(ibz) = q2(ibz)

! calculate the energy for the Alfvén intermediate state (eq. 63)
!
                u2(ien) = u1r(ien) + bs * drsq * (vb1r - vb2)

! calculate the right Alfvén intermediate flux (eq. 65)
!
                f(:,i) = 0.5d0 * (f(:,i) + (fr(:,i)                            &
                          + smr * u2(:) - (smr - sr) * u1r(:) - sr * ur(:,i)))

              end if ! sm = 0

            end if ! sml ≤ 0 and smr ≥ 0

          else ! Bx = 0

! intermediate state for the case of Bx = 0
!
            if (sm > 0.0d0) then

! calculate the left intermediate flux, eq. (64)
!
              f(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

            else if (sm < 0.0d0) then

! calculate the right intermediate flux, eq. (64)
!
              f(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

            else ! sm = 0

! average the left and right fluxes if both Sm = 0, and Bx = 0
!
              f(:,i) = 0.5d0 * ((fl(:,i) + sl * (u1l(:) - ul(:,i)))            &
                             +  (fr(:,i) + sr * (u1r(:) - ur(:,i))))

            end if ! sm = 0

          end if ! Bx = 0

#ifdef SHOCK_DETECTION
        end if ! shock vicinity
#endif /* SHOCK_DETECTION */

      end if ! sl < 0 < sr

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* ADI */
#ifdef ISO
!
!===============================================================================
!
! subroutine RIEMANN:
! ------------------
!
!   Subroutine solves one dimensional Riemann problem using the isothermal HLLD
!   method described by Mignone.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     b - the input vector of the normal magnetic field component;
!     f - the output array of fluxes;
!     s - the input array of shock indicators;
!
!   References:
!
!     [1] Mignone, A.,
!         "A simple and accurate Riemann solver for isothermal MHD",
!         Journal of Computational Physics, 2007, 225, pp. 1427-1441
!
!===============================================================================
!
  subroutine riemann(n, h, q, b, f, s)

! include external procedures and variables
!
    use equations     , only : csnd
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity
    use variables     , only : nt
    use variables     , only : idn, imx, imy, imz
    use variables     , only : ivx, ivy, ivz
    use variables     , only : ibx, iby, ibz

! local variables are not implicit by default
!
    implicit none

! input / output arguments
!
    integer                 , intent(in)  :: n
    real                    , intent(in)  :: h
    real   , dimension(nt,n), intent(in)  :: q
    real   , dimension(   n), intent(in)  :: b
    real   , dimension(nt,n), intent(out) :: f
    logical, dimension(n)   , intent(in)  :: s

! local variables
!
    integer :: p, i
#ifdef SHOCK_DETECTION
    integer :: ip1
#endif /* SHOCK_DETECTION */
    real    :: dv, dn, mx, sq, fc, bs, b2, ca, dnl, dnr
    real    :: sl, sr, srl, srml, sm, sml, smr, slmv, srmv, slmm, srmm

! local arrays to store the states
!
    real, dimension(nt,n) :: ul, ur, ql, qr, fl, fr
    real, dimension(n)    :: cl, cr
    real, dimension(nt)   :: u1l, u1r, u2
!
!-------------------------------------------------------------------------------
!
! reconstruct the left and right states from primitive variables
!
    do p = 1, ibx - 1
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do
    do p = ibx + 1, nt
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:), s(:))
    end do

! copy the normal component of magnetic field to the left and right states
!
    ql(ibx,:) = b(:)
    qr(ibx,:) = b(:)

#ifdef FIX_POSITIVITY
! check if the reconstruction doesn't give negative densities, if so, correct
! the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
#endif /* FIX_POSITIVITY */

! prepare conservative variables of both states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate physical fluxes and speeds for the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLD flux
!
      if (sl >= 0.0d0) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d0) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! prepare some coefficients
!
        srl  = sr * sl
        srml = sr - sl

#ifdef SHOCK_DETECTION
! calculate the next index
!
        ip1 = min(n, i + 1)

! check if we are in the vicinity of a shock
!
        if (s(i) .or. s(ip1)) then

! calculate the fluxes for the intermediate state
!
          f(:,i) = ((sr * fl(:,i) - sl * fr(:,i))                              &
                                           + srl * (ur(:,i) - ul(:,i))) / srml

        else
#endif /* SHOCK_DETECTION */

! density and normal momentum are constant acros the intermediate states
! (eq. 20 and 21)
!
          dn = ((sr * ur(idn,i) - sl * ul(idn,i))                              &
                                             - (fr(idn,i) - fl(idn,i))) / srml
          mx = ((sr * ur(imx,i) - sl * ul(imx,i))                              &
                                             - (fr(imx,i) - fl(imx,i))) / srml

! calculate corresponding fluxes for density and normal momentum for
! the intermediate states (eq. 22 and 23)
!
          f(idn,i) = ((sr * fl(idn,i) - sl * fr(idn,i))                        &
                                       + srl * (ur(idn,i) - ul(idn,i))) / srml
          f(imx,i) = ((sr * fl(imx,i) - sl * fr(imx,i))                        &
                                       + srl * (ur(imx,i) - ul(imx,i))) / srml
          f(ibx,i) = ((sr * fl(ibx,i) - sl * fr(ibx,i))                        &
                                       + srl * (ur(ibx,i) - ul(ibx,i))) / srml

! the speed of advection in the intermediate states (calculated from eq. 15
! and eq. 17)
!
          sm   = f(idn,i) / dn

! calculate speed differences
!
          slmv = sl - ql(ivx,i)
          srmv = sr - qr(ivx,i)
          slmm = sl - sm
          srmm = sr - sm

! prepare some coefficients
!
          dnl  = ql(idn,i) * slmv
          dnr  = qr(idn,i) * srmv
          sq   = sqrt(dn)
          b2   = b(i) * b(i)

! calculate Alfvén speed (eq. 29)
!
          ca   = abs(b(i)) / sq

! calculate speeds of left and right going Alfvén waves
!
          sml  = sm - ca
          smr  = sm + ca

! calculate the density and normal component of the momentum and magnetic
! field for the left and right intermediate states
!
          u1l(idn) = dn
          u1l(imx) = mx
          u1l(ibx) = b(i)
          u1r(idn) = dn
          u1r(imx) = mx
          u1r(ibx) = b(i)

! calculate tangential components of the magnetic field for the left and right
! intermediate states (like eq. 32-33, but use expressions as in the adabatic
! case and in the degeneracy takes place, revert tangential components to
! the HLL averages)
!
          dv = dnl * slmm - b2
          if (dv /= 0.0d0 .and. ca <= csnd) then
            fc       = (dnl * slmv - b2) / dv
            u1l(iby) = ul(iby,i) * fc
            u1l(ibz) = ul(ibz,i) * fc
          else
            u1l(iby) = ((sr * ur(iby,i) - sl * ul(iby,i))                      &
                                             - (fr(iby,i) - fl(iby,i))) / srml
            u1l(ibz) = ((sr * ur(ibz,i) - sl * ul(ibz,i))                      &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml
          end if

          dv = dnr * srmm - b2
          if (dv /= 0.0d0 .and. ca <= csnd) then
            fc       = (dnr * srmv - b2) / dv
            u1r(iby) = ur(iby,i) * fc
            u1r(ibz) = ur(ibz,i) * fc
          else
            u1r(iby) = ((sr * ur(iby,i) - sl * ul(iby,i))                      &
                                             - (fr(iby,i) - fl(iby,i))) / srml
            u1r(ibz) = ((sr * ur(ibz,i) - sl * ul(ibz,i))                      &
                                             - (fr(ibz,i) - fl(ibz,i))) / srml
          end if

! calculate tangential components of the momentum for the left and right
! intermediate states (like eqs. 30-31, but using already calculated magnetic
! field components for these states)
!
          u1l(imy) = (ul(imy,i) * slmv + b(i) * (ul(iby,i) - u1l(iby))) / slmm
          u1l(imz) = (ul(imz,i) * slmv + b(i) * (ul(ibz,i) - u1l(ibz))) / slmm
          u1r(imy) = (ur(imy,i) * srmv + b(i) * (ur(iby,i) - u1r(iby))) / srmm
          u1r(imz) = (ur(imz,i) * srmv + b(i) * (ur(ibz,i) - u1r(ibz))) / srmm

! calculate fluxes for the intermediate states
!
          if (sml > 0.0d0) then

! calculate the left intermediate flux, eq. (38)
!
            f(imy,i) = fl(imy,i) + sl * (u1l(imy) - ul(imy,i))
            f(imz,i) = fl(imz,i) + sl * (u1l(imz) - ul(imz,i))
            f(iby,i) = fl(iby,i) + sl * (u1l(iby) - ul(iby,i))
            f(ibz,i) = fl(ibz,i) + sl * (u1l(ibz) - ul(ibz,i))

          else if (smr < 0.0d0) then

! calculate the right intermediate flux, eq. (38)
!
            f(imy,i) = fr(imy,i) + sr * (u1r(imy) - ur(imy,i))
            f(imz,i) = fr(imz,i) + sr * (u1r(imz) - ur(imz,i))
            f(iby,i) = fr(iby,i) + sr * (u1r(iby) - ur(iby,i))
            f(ibz,i) = fr(ibz,i) + sr * (u1r(ibz) - ur(ibz,i))

          else ! sml ≤ 0 ≤ smr

! calculate the intermediate state (eq. 34-37)
!
            if (b(i) > 0.0d0) then

              u2(imy) = 0.5d0 * ((u1l(imy) + u1r(imy))                         &
                                                 + sq * (u1r(iby) - u1l(iby)))
              u2(imz) = 0.5d0 * ((u1l(imz) + u1r(imz))                         &
                                                 + sq * (u1r(ibz) - u1l(ibz)))

              u2(iby) = 0.5d0 * ((u1l(iby) + u1r(iby))                         &
                                                 + (u1r(imy) - u1l(imy)) / sq)
              u2(ibz) = 0.5d0 * ((u1l(ibz) + u1r(ibz))                         &
                                                 + (u1r(imz) - u1l(imz)) / sq)

            else if (b(i) < 0.0d0) then

              u2(imy) = 0.5d0 * ((u1l(imy) + u1r(imy))                         &
                                                 - sq * (u1r(iby) - u1l(iby)))
              u2(imz) = 0.5d0 * ((u1l(imz) + u1r(imz))                         &
                                                 - sq * (u1r(ibz) - u1l(ibz)))

              u2(iby) = 0.5d0 * ((u1l(iby) + u1r(iby))                         &
                                                 - (u1r(imy) - u1l(imy)) / sq)
              u2(ibz) = 0.5d0 * ((u1l(ibz) + u1r(ibz))                         &
                                                 - (u1r(imz) - u1l(imz)) / sq)
            else

              u2(imy) = 0.5d0 * (u1l(imy) + u1r(imy))
              u2(imz) = 0.5d0 * (u1l(imz) + u1r(imz))

              u2(iby) = 0.5d0 * (u1l(iby) + u1r(iby))
              u2(ibz) = 0.5d0 * (u1l(ibz) + u1r(ibz))

            end if

! calculate the intermediate flux (eq. 24)
!
            f(imy,i) = sm * u2(imy) - b(i) * u2(iby)
            f(imz,i) = sm * u2(imz) - b(i) * u2(ibz)
            f(iby,i) = sm * u2(iby) - b(i) * u2(imy) / dn
            f(ibz,i) = sm * u2(ibz) - b(i) * u2(imz) / dn

          end if ! sml ≤ 0 ≤ smr

#ifdef SHOCK_DETECTION
        end if ! shock vicinity
#endif /* SHOCK_DETECTION */

      end if ! sl < 0 < sr

    end do

!-------------------------------------------------------------------------------
!
  end subroutine riemann
#endif /* ISO */
#endif /* HLLD */
#endif /* MHD */

!===============================================================================
!
end module schemes
