!!******************************************************************************
!!
!! module: scheme - handling the actual solver of the set of equations
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
module scheme

  implicit none

! the maximal speed in the system
!
  real, save :: cmax

  contains
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
    use variables    , only : nvr, nqt, nfl
    use variables    , only : idn, imx, imy, imz
#ifdef ADI
    use variables    , only : ien
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    real, dimension(nqt,im,jm,km)      , intent(in)  :: u
    real, dimension(nqt,im,jm,km)      , intent(out) :: du
    real                               , intent(in)  :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k
    real    :: dx, dy, dz

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
!
!-------------------------------------------------------------------------------
!
! reset the increment array
!
    du(:,:,:,:) = 0.0

! prepare the spacial increment
!
    dx = 1.0d0 / dxi
    dy = 1.0d0 / dyi
#if NDIMS == 3
    dz = 1.0d0 / dzi
#endif /* NDIMS == 3 */

!  update along X-direction
!
    do k = 1, km
      do j = 1, jm

! copy directional vectors of variables for the one dimensional solver
!
        do i = 1, im
          ux(idn,i) = u(idn,i,j,k)
          ux(imx,i) = u(imx,i,j,k)
          ux(imy,i) = u(imy,i,j,k)
          ux(imz,i) = u(imz,i,j,k)
#ifdef ADI
          ux(ien,i) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          ux(ibx,i) = u(ibx,i,j,k)
          ux(iby,i) = u(iby,i,j,k)
          ux(ibz,i) = u(ibz,i,j,k)
#ifdef GLM
          ux(iph,i) = u(iph,i,j,k)
#endif /* GLM */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (im, dx, ux(:,:), fx(:,:))
#endif /* HLL */
#ifdef HLLC
        call hllc(im, dx, ux(:,:), fx(:,:))
#endif /* HLLC */
#ifdef HLLD
        call hlld(im, dx, ux(:,:), fx(:,:))
#endif /* HLLD */
#ifdef ROE
        call roe (im, dx, ux(:,:), fx(:,:))
#endif /* ROE */

! update the arrays of increments
!
        do i = 1, im

! update fluid variables
!
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
          du(ibx,i,j,k) = du(ibx,i,j,k) + dxi * fx(ibx,i)
          du(iby,i,j,k) = du(iby,i,j,k) + dxi * fx(iby,i)
          du(ibz,i,j,k) = du(ibz,i,j,k) + dxi * fx(ibz,i)
#ifdef GLM

! update scalar potential
!
          du(iph,i,j,k) = du(iph,i,j,k) + dxi * fx(iph,i)
#endif /* GLM */
#endif /* MHD */
        end do
      end do
    end do

!  update along Y-direction
!
    do k = 1, km
      do i = 1, im

! copy directional vectors of variables for the one dimensional solver
!
        do j = 1, jm
          uy(idn,j) = u(idn,i,j,k)
          uy(imx,j) = u(imy,i,j,k)
          uy(imy,j) = u(imz,i,j,k)
          uy(imz,j) = u(imx,i,j,k)
#ifdef ADI
          uy(ien,j) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          uy(ibx,j) = u(iby,i,j,k)
          uy(iby,j) = u(ibz,i,j,k)
          uy(ibz,j) = u(ibx,i,j,k)
#ifdef GLM
          uy(iph,j) = u(iph,i,j,k)
#endif /* GLM */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (jm, dy, uy(:,:), fy(:,:))
#endif /* HLL */
#ifdef HLLC
        call hllc(jm, dy, uy(:,:), fy(:,:))
#endif /* HLLC */
#ifdef HLLD
        call hlld(jm, dy, uy(:,:), fy(:,:))
#endif /* HLLD */
#ifdef ROE
        call roe (jm, dy, uy(:,:), fy(:,:))
#endif /* ROE */

! update the arrays of increments
!
        do j = 1, jm

! update fluid variables
!
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
          du(ibx,i,j,k) = du(ibx,i,j,k) + dyi * fy(ibz,j)
          du(iby,i,j,k) = du(iby,i,j,k) + dyi * fy(ibx,j)
          du(ibz,i,j,k) = du(ibz,i,j,k) + dyi * fy(iby,j)
#ifdef GLM

! update scalar potential
!
          du(iph,i,j,k) = du(iph,i,j,k) + dyi * fy(iph,j)
#endif /* GLM */
#endif /* MHD */
        end do
      end do
    end do

#if NDIMS == 3
!  update along Z-direction
!
    do j = 1, jm
      do i = 1, im

! copy directional vectors of variables for the one dimensional solver
!
        do k = 1, km
          uz(idn,k) = u(idn,i,j,k)
          uz(imx,k) = u(imz,i,j,k)
          uz(imy,k) = u(imx,i,j,k)
          uz(imz,k) = u(imy,i,j,k)
#ifdef ADI
          uz(ien,k) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          uz(ibx,k) = u(ibz,i,j,k)
          uz(iby,k) = u(ibx,i,j,k)
          uz(ibz,k) = u(iby,i,j,k)
#ifdef GLM
          uz(iph,k) = u(iph,i,j,k)
#endif /* GLM */
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll (km, dz, uz(:,:), fz(:,:))
#endif /* HLL */
#ifdef HLLC
        call hllc(km, dz, uz(:,:), fz(:,:))
#endif /* HLLC */
#ifdef HLLD
        call hlld(km, dz, uz(:,:), fz(:,:))
#endif /* HLLD */
#ifdef ROE
        call roe (km, dz, uz(:,:), fz(:,:))
#endif /* ROE */

! update the arrays of increments
!
        do k = 1, km

! update fluid variables
!
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
          du(ibx,i,j,k) = du(ibx,i,j,k) + dzi * fz(iby,k)
          du(iby,i,j,k) = du(iby,i,j,k) + dzi * fz(ibz,k)
          du(ibz,i,j,k) = du(ibz,i,j,k) + dzi * fz(ibx,k)
#ifdef GLM

! update scalar potential
!
          du(iph,i,j,k) = du(iph,i,j,k) + dzi * fz(iph,k)
#endif /* GLM */
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
! hll: subroutine computes the approximated flux using the HLL method
!
!===============================================================================
!
  subroutine hll(n, h, u, f)

#ifdef VISCOSITY
    use config       , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config       , only : ueta
#endif /* MHD & RESISTIVITY */
    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : ivx, ivy, ivz
#ifdef ADI
    use variables    , only : ien
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real                  , intent(in)  :: h
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f

! local variables
!
    integer                :: p, i, ip1
    real, dimension(nvr,n) :: q, ql, qr, ul, ur
    real, dimension(nqt,n) :: fl, fr, fn
    real, dimension(n)     :: cl, cr
    real                   :: al, ar, ap, div
#ifdef VISCOSITY
    real                   :: dvx, dvy, dvz, vih
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    real                   :: dbx, dby, dbz, ueh
#endif /* MHD & RESISTIVITY */
!
!-------------------------------------------------------------------------------
!
! usefull parameters
!
#ifdef VISCOSITY
    vih = visc / h
#endif /* VISCOSITY */
#ifdef RESISTIVITY
    ueh = ueta / h
#endif /* RESISTIVITY */

! calculate the primitive variables
!
    call cons2prim(n, u(:,:), q(:,:))

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef MHD
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
#endif /* MHD */

! calculate conservative variables at states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate fluxes and speeds
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      al = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      ar = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate HLL flux
!
      if (al .ge. 0.0d0) then
        fn(:,i) = fl(:,i)
      else if (ar .le. 0.0d0) then
        fn(:,i) = fr(:,i)
      else
        ap  = ar * al
        div = 1.0d0 / (ar - al)

        fn(:,i) = div * (ar * fl(:,i) - al * fr(:,i) + ap * (ur(:,i) - ul(:,i)))
      end if
    end do

#ifdef VISCOSITY
! add viscous term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dvx = vih * (q(ivx,ip1) - q(ivx,i))
      fn(ivx,i  ) = fn(ivx,i  ) - dvx

      dvy = vih * (q(ivy,ip1) - q(ivy,i))
      fn(ivy,i  ) = fn(ivy,i  ) - dvy

      dvz = vih * (q(ivz,ip1) - q(ivz,i))
      fn(ivz,i  ) = fn(ivz,i  ) - dvz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - 0.5d0 * (ql(ivx,i) + qr(ivx,i)) * dvx
#endif /* ADI */
    end do

#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
! add resistivity term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dbx = ueh * (q(ibx,ip1) - q(ibx,i))
      fn(ibx,i) = fn(ibx,i) - dbx

      dby = ueh * (q(iby,ip1) - q(iby,i))
      fn(iby,i) = fn(iby,i) - dby

      dbz = ueh * (q(ibz,ip1) - q(ibz,i))
      fn(ibz,i) = fn(ibz,i) - dbz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - ql(ibx,i) * dbx
#endif /* ADI */
    end do

#endif /* MHD & RESISTIVITY */
! calculate numerical flux
!
    f(  1:nfl,2:n) = - fn(  1:nfl,2:n) + fn(  1:nfl,1:n-1)
#ifdef MHD
#ifdef GLM
    f(ibx:ibz,2:n) = - fn(ibx:ibz,2:n) + fn(ibx:ibz,1:n-1)
    f(iph    ,2:n) = - fn(iph    ,2:n) + fn(iph    ,1:n-1)
#endif /* GLM */
#endif /* MHD */

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
  subroutine hllc(n, h, u, f)

#ifdef VISCOSITY
    use config       , only : visc
#endif /* VISCOSITY */
    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : idn, imx, imy, imz, ien, ivx, ivy, ivz, ipr

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real                  , intent(in)  :: h
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f

! local variables
!
    integer                :: p, i, ip1
    real, dimension(nvr,n) :: ql, qr, q, ul, ur
    real, dimension(nfl,n) :: fl, fr, fn
    real, dimension(n)     :: cl, cr, cm
    real                   :: sl, sr, sm, sml, smr, srmv, slmv, srmm, slmm     &
                            , smvl, smvr, div, pt
    real, dimension(nvr)   :: q1l, q1r, u1l, u1r
#ifdef VISCOSITY
    real                   :: dvx, dvy, dvz, vih
#endif /* VISCOSITY */
!
!-------------------------------------------------------------------------------
!
#ifdef VISCOSITY
! usefull parameters
!
    vih = visc / h

#endif /* VISCOSITY */
! obtain the primitive variables
!
    call cons2prim(n, u(:,:), q(:,:))

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the conservative variables at both states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate along the direction
!
    do i = 1, n

! calculate the minimum and maxximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! all speeds >= 0, left side flux
!
      if (sl .ge. 0.0) then

        fn(:,i) = fl(:,i)

! all speeds <= 0, right side flux
!
      else if (sr .le. 0.0) then

        fn(:,i) = fr(:,i)

! intermediate states
!
      else  ! sl < 0 & sr > 0

! useful differences
!
        slmv = sl - ql(ivx,i)
        srmv = sr - qr(ivx,i)

! speed of the contact discontinuity (eq. 34 [1], 14 [2])
!
        div =  srmv * qr(idn,i) - slmv * ql(idn,i)
        sml = (srmv * ur(imx,i) - slmv * ul(imx,i)                             &
                                - qr(ipr,i) + ql(ipr,i)) / div
        div =  slmv * ql(idn,i) - srmv * qr(idn,i)
        smr = (slmv * ul(imx,i) - srmv * ur(imx,i)                             &
                                - ql(ipr,i) + qr(ipr,i)) / div
        sm  = 0.5d0 * (sml + smr)

        if (sm .eq. 0.0d0) then

! calculate the left intermediate state
!
          pt = ql(ipr,i) - ul(imx,i) * slmv
          u1l(idn) = ql(idn,i) * slmv / sl
          u1l(imx) = 0.0d0
          u1l(imy) = u1l(idn) * ql(ivy,i)
          u1l(imz) = u1l(idn) * ql(ivz,i)
          if (sl .eq. 0.0d0) then
            u1l(ien) = ul(ien,i)
          else
            u1l(ien) = (slmv * ul(ien,i) - ql(ipr,i) * ql(ivx,i)) / sl
          end if

! calculate right intermediate state
!
          pt = qr(ipr,i) - ur(imx,i) * srmv
          u1r(idn) = qr(idn,i) * srmv / sr
          u1r(imx) = 0.0d0
          u1r(imy) = u1r(idn) * qr(ivy,i)
          u1r(imz) = u1r(idn) * qr(ivz,i)
          if (sr .eq. 0.0d0) then
            u1r(ien) = ur(ien,i)
          else
            u1r(ien) = (srmv * ur(ien,i) - qr(ipr,i) * qr(ivx,i)) / sr
          endif

! calculate intermediate flux
!
          fn(:,i) = 0.5d0 * (fl(:,i) + sl * (u1l(:) - ul(:,i))                 &
                          +  fr(:,i) + sr * (u1r(:) - ur(:,i)))

        else

! useful differences
!
          slmm = sl - sm
          srmm = sr - sm
          smvl = sm - ql(ivx,i)
          smvr = sm - qr(ivx,i)

! intermediate discontinuities
!
          if (sm .gt. 0.0d0) then

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
            pt = ql(ipr,i) + ql(idn,i) * slmv * smvl

! calculate the left intermediate state
!
            u1l(idn) = ql(idn,i) * slmv / slmm
            u1l(imx) = u1l(idn) * sm
            u1l(imy) = u1l(idn) * ql(ivy,i)
            u1l(imz) = u1l(idn) * ql(ivz,i)
            if (slmm .eq. 0.0d0) then
              u1l(ien) = ul(ien,i)
            else
              u1l(ien) = (slmv * ul(ien,i) - ql(ipr,i) * ql(ivx,i)             &
                                                            + pt * sm) / slmm
            end if

! calculate the left intermediate flux
!
            fn(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

          else if (sm .lt. 0.0) then

! pressure of intermediate states (eq. 36 [1], 16 [2])
!
            pt = qr(ipr,i) + qr(idn,i) * srmv * smvr

! calculate the right intermediate state
!
            u1r(idn) = qr(idn,i) * srmv / srmm
            u1r(imx) = u1r(idn) * sm
            u1r(imy) = u1r(idn) * qr(ivy,i)
            u1r(imz) = u1r(idn) * qr(ivz,i)
            if (srmm .eq. 0.0d0) then
              u1r(ien) = ur(ien,i)
            else
              u1r(ien) = (srmv * ur(ien,i) - qr(ipr,i) * qr(ivx,i)             &
                                                            + pt * sm) / srmm
            end if

! calculate the right intermediate flux
!
            fn(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

          end if
        end if
      end if

    end do

#ifdef VISCOSITY
! add viscous term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dvx = vih * (q(ivx,ip1) - q(ivx,i))
      fn(ivx,i  ) = fn(ivx,i  ) - dvx

      dvy = vih * (q(ivy,ip1) - q(ivy,i))
      fn(ivy,i  ) = fn(ivy,i  ) - dvy

      dvz = vih * (q(ivz,ip1) - q(ivz,i))
      fn(ivz,i  ) = fn(ivz,i  ) - dvz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - 0.5d0 * (ql(ivx,i) + qr(ivx,i)) * dvx
#endif /* ADI */
    end do

#endif /* VISCOSITY */
! calculate numerical flux
!
    f(:,2:n) = - fn(:,2:n) + fn(:,1:n-1)

!-------------------------------------------------------------------------------
!
  end subroutine hllc
#endif /* HLLC */
#ifdef MHD
#ifdef HLLD
#ifdef ISO
!
!===============================================================================
!
! hlld: subroutine computes the approximated flux using the HLLD method
!       for the isothermal equation of state
!
!===============================================================================
!
  subroutine hlld(n, h, u, f)

#ifdef VISCOSITY
    use config       , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config       , only : ueta
#endif /* MHD & RESISTIVITY */
    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : idn, imx, imy, imz, ivx, ivy, ivz
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real                  , intent(in)  :: h
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f

! local variables
!
    integer                :: p, i, ip1
    real, dimension(nvr,n) :: q, ql, qr, ul, ur
    real, dimension(nqt,n) :: fl, fr, fn
    real, dimension(n)     :: cl, cr
    real, dimension(nvr)   :: u1l, u1r, u2
    real                   :: sl, sr, srl, srml, sm, sml, smr
    real                   :: dnm, mxm, sqd, div, fac, bxs
#ifdef VISCOSITY
    real                   :: dvx, dvy, dvz, vih
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    real                   :: dbx, dby, dbz, ueh
#endif /* MHD & RESISTIVITY */
!
!-------------------------------------------------------------------------------
!
! usefull parameters
!
#ifdef VISCOSITY
    vih = visc / h
#endif /* VISCOSITY */
#ifdef RESISTIVITY
    ueh = ueta / h
#endif /* RESISTIVITY */

! calculate the primitive variables
!
    call cons2prim(n, u, q)

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
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

! calculate conservative variables at states
!
    call prim2cons(n, ql, ul)
    call prim2cons(n, qr, ur)

! calculate fluxes and speeds
!
    call fluxspeed(n, ql, ul, fl, cl)
    call fluxspeed(n, qr, ur, fr, cr)

! iterate over all points and calculate the HLLD flux
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! all speeds > 0, left side flux
!
      if (sl .ge. 0.0) then

        fn(:,i) = fl(:,i)

! all speeds < 0, right side flux
!
      else if (sr .le. 0.0) then

        fn(:,i) = fr(:,i)

! intermediate states
!
      else  ! sl < 0 & sr > 0

! product and difference of speeds
!
        srl  = sr * sl
        srml = sr - sl

! density of the intermediate state (eq. 20 and 21)
!
        dnm = (sr * ur(idn,i) - sl * ul(idn,i) - fr(idn,i) + fl(idn,i)) / srml
        mxm = (sr * ur(imx,i) - sl * ul(imx,i) - fr(imx,i) + fl(imx,i)) / srml
        sqd = sqrt(dnm)

! fluxes for density and x-momentum are the same for all intermediate states (eq. 22 and 23)
!
        fn(idn,i) = (sr * fl(idn,i) - sl * fr(idn,i)                           &
                                    + srl * (ur(idn,i) - ul(idn,i))) / srml
        fn(imx,i) = (sr * fl(imx,i) - sl * fr(imx,i)                           &
                                    + srl * (ur(imx,i) - ul(imx,i))) / srml

#ifdef GLM
! fluxes for parallel magnetic component and the scalar potential is the same
! as well
        fn(ibx,i) = (sr * fl(ibx,i) - sl * fr(ibx,i)                           &
                                    + srl * (ur(ibx,i) - ul(ibx,i))) / srml
        fn(iph,i) = (sr * fl(iph,i) - sl * fr(iph,i)                           &
                                    + srl * (ur(iph,i) - ul(iph,i))) / srml
#endif /* GLM */

! the speed of contact discontinuity (from eq. 15 and eq. 17)
!
        sm  = fn(idn,i) / dnm

! Alfven speeds (eq. 29)
!
        sml = sm - abs(ql(ibx,i)) / sqd
        smr = sm + abs(qr(ibx,i)) / sqd

! calculate the left intermediate state
!
        u1l(idn) = dnm
        u1l(imx) = mxm

        div    = (sl - sml) * (sl - smr)
        if (sm .eq. ql(ivx,i) .or. div .eq. 0.0 .or. ql(ibx,i) .eq. 0.0) then
          u1l(imy) = dnm * ql(ivy,i)
          u1l(imz) = dnm * ql(ivz,i)
          u1l(iby) = ql(iby,i)
          u1l(ibz) = ql(ibz,i)
        else
          fac      = ql(ibx,i) * (sm - ql(ivx,i)) / div
          u1l(imy) = dnm * ql(ivy,i) - ql(iby,i) * fac
          u1l(imz) = dnm * ql(ivz,i) - ql(ibz,i) * fac
          fac      = (ql(idn,i) * (sl - ql(ivx,i))**2 - ql(ibx,i)**2)          &
                                                                / (dnm * div)
          u1l(iby) = ql(iby,i) * fac
          u1l(ibz) = ql(ibz,i) * fac
        end if

! calculate the right intermediate state
!
        u1r(idn) = dnm
        u1r(imx) = mxm

        div    = (sr - sml) * (sr - smr)
        if (sm .eq. qr(ivx,i) .or. div .eq. 0.0 .or. qr(ibx,i) .eq. 0.0) then
          u1r(imy) = dnm * qr(ivy,i)
          u1r(imz) = dnm * qr(ivz,i)
          u1r(iby) = qr(iby,i)
          u1r(ibz) = qr(ibz,i)
        else
          fac      = qr(ibx,i) * (sm - qr(ivx,i)) / div
          u1r(imy) = dnm * qr(ivy,i) - qr(iby,i) * fac
          u1r(imz) = dnm * qr(ivz,i) - qr(ibz,i) * fac
          fac      = (qr(idn,i) * (sr - qr(ivx,i))**2 - qr(ibx,i)**2)          &
                                                                / (dnm * div)
          u1r(iby) = qr(iby,i) * fac
          u1r(ibz) = qr(ibz,i) * fac
        end if

! intermediate discontinuities
!
        if (sml .ge. 0.0) then

! calculate the left intermediate flux
!
          fn(imy,i) = fl(imy,i) + sl * (u1l(imy) - ul(imy,i))  ! eq. (38)
          fn(imz,i) = fl(imz,i) + sl * (u1l(imz) - ul(imz,i))
          fn(iby,i) = fl(iby,i) + sl * (u1l(iby) - ul(iby,i))
          fn(ibz,i) = fl(ibz,i) + sl * (u1l(ibz) - ul(ibz,i))

        else if (smr .le. 0.0) then

! calculate right intermediate flux
!
          fn(imy,i) = fr(imy,i) + sr * (u1r(imy) - ur(imy,i))  ! eq. (38)
          fn(imz,i) = fr(imz,i) + sr * (u1r(imz) - ur(imz,i))
          fn(iby,i) = fr(iby,i) + sr * (u1r(iby) - ur(iby,i))
          fn(ibz,i) = fr(ibz,i) + sr * (u1r(ibz) - ur(ibz,i))

        else ! sml < 0 & smr > 0

! normal component of magnetic field multiplied by sqrt(dnm)
!
          if (ql(ibx,i) .ge. 0.0) then
            bxs =   sqd
          else
            bxs = - sqd
          end if

! calculate the intermediate state (eq. 34-37)
!
          u2(imy) = 0.5d0 * (u1r(imy) + u1l(imy) + bxs * (u1r(iby) - u1l(iby)))
          u2(imz) = 0.5d0 * (u1r(imz) + u1l(imz) + bxs * (u1r(ibz) - u1l(ibz)))
          u2(iby) = 0.5d0 * (u1r(iby) + u1l(iby) + (u1r(imy) - u1l(imy)) / bxs)
          u2(ibz) = 0.5d0 * (u1r(ibz) + u1l(ibz) + (u1r(imz) - u1l(imz)) / bxs)

! calculate the intermediate flux (eq. 24)
!
          fn(imy,i) = sm * u2(imy) - ql(ibx,i) * u2(iby)
          fn(imz,i) = sm * u2(imz) - ql(ibx,i) * u2(ibz)
          fn(iby,i) = sm * u2(iby) - ql(ibx,i) * u2(imy) / dnm
          fn(ibz,i) = sm * u2(ibz) - ql(ibx,i) * u2(imz) / dnm

        end if
      end if

    end do

#ifdef VISCOSITY
! add viscous term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dvx = vih * (q(ivx,ip1) - q(ivx,i))
      fn(ivx,i  ) = fn(ivx,i  ) - dvx

      dvy = vih * (q(ivy,ip1) - q(ivy,i))
      fn(ivy,i  ) = fn(ivy,i  ) - dvy

      dvz = vih * (q(ivz,ip1) - q(ivz,i))
      fn(ivz,i  ) = fn(ivz,i  ) - dvz
    end do

#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
! add resistivity term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dbx = ueh * (q(ibx,ip1) - q(ibx,i))
      fn(ibx,i) = fn(ibx,i) - dbx

      dby = ueh * (q(iby,ip1) - q(iby,i))
      fn(iby,i) = fn(iby,i) - dby

      dbz = ueh * (q(ibz,ip1) - q(ibz,i))
      fn(ibz,i) = fn(ibz,i) - dbz
    end do

#endif /* MHD & RESISTIVITY */
! calculate numerical flux
!
    f(  1:nfl,2:n) = - fn(  1:nfl,2:n) + fn(   1:nfl,1:n-1)
    f(ibx:ibz,2:n) = - fn(ibx:ibz,2:n) + fn(ibx:ibz,1:n-1)
#ifdef GLM
    f(iph    ,2:n) = - fn(iph    ,2:n) + fn(iph    ,1:n-1)
#endif /* GLM */

!-------------------------------------------------------------------------------
!
  end subroutine hlld
#endif /* ISO */
#ifdef ADI
!
!===============================================================================
!
! hlld: subroutine computes the approximated flux using the HLLD method
!       for the adiabatic equation of state
!
!===============================================================================
!
  subroutine hlld(n, h, u, f)

    use config       , only : gamma
#ifdef VISCOSITY
    use config       , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config       , only : ueta
#endif /* MHD & RESISTIVITY */
    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : idn, imx, imy, imz, ien, ivx, ivy, ivz, ipr
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real                  , intent(in)  :: h
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f

! local variables
!
    integer                :: p, i, ip1
    real, dimension(nvr,n) :: q, ql, qr, ul, ur
    real, dimension(nqt,n) :: fl, fr, fn
    real, dimension(n)     :: cl, cr
    real, dimension(nvr)   :: u1l, u1r, u2, q1l, q1r, q2
    real                   :: sl, sr, slmv, srmv, slmm, srmm, sm, smvl, smvr   &
                            , sml, smr
    real                   :: ptl, ptr, pt, bx2, div, fac, bxs, dlsq, drsq
#ifdef VISCOSITY
    real                   :: dvx, dvy, dvz, vih
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    real                   :: dbx, dby, dbz, ueh
#endif /* MHD & RESISTIVITY */
!
!-------------------------------------------------------------------------------
!
! usefull parameters
!
#ifdef VISCOSITY
    vih = visc / h
#endif /* VISCOSITY */
#ifdef RESISTIVITY
    ueh = ueta / h
#endif /* RESISTIVITY */

! calculate the primitive variables
!
    call cons2prim(n, u, q)

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
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

! calculate conservative variables at states
!
    call prim2cons(n, ql, ul)
    call prim2cons(n, qr, ur)

! calculate fluxes and speeds
!
    call fluxspeed(n, ql, ul, fl, cl)
    call fluxspeed(n, qr, ur, fr, cr)

! iterate over all points and calculate the HLLD flux
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! all speeds > 0, left side flux
!
      if (sl .ge. 0.0) then

        fn(:,i) = fl(:,i)

! all speeds < 0, right side flux
!
      else if (sr .le. 0.0) then

        fn(:,i) = fr(:,i)

! intermediate states
!
      else  ! sl < 0 & sr > 0

! calculate the total left and right pressures
!
        ptl = ql(ipr,i) + 0.5d0 * sum(ql(ibx:ibz,i)**2)
        ptr = qr(ipr,i) + 0.5d0 * sum(qr(ibx:ibz,i)**2)

! useful speed differences
!
        slmv = sl - ql(ivx,i)
        srmv = sr - qr(ivx,i)

! the speed of contact discontinuity (eq. 38, average from the both states)
!
        div  =  slmv * ql(idn,i) - srmv * qr(idn,i)
        slmm = (slmv * ul(imx,i) - srmv * ur(imx,i) - ptl + ptr) / div
        div  =  srmv * qr(idn,i) - slmv * ql(idn,i)
        srmm = (srmv * ur(imx,i) - slmv * ul(imx,i) - ptr + ptl) / div
        sm   = 0.5d0 * (slmm + srmm)

! more useful speed differences
!
        slmm = sl - sm
        srmm = sr - sm
        smvl = sm - ql(ivx,i)
        smvr = sm - qr(ivx,i)
        bx2  = ql(ibx,i) * qr(ibx,i)

! pressure of the intermediate states (eq. 41)
!
        pt   = 0.5d0 * (ptl + ptr + ql(idn,i) * slmv * smvl                    &
                                  + qr(idn,i) * srmv * smvr)

! calculate the left intermediate state variables
!
        q1l(idn) = ql(idn,i) * slmv / slmm
        q1l(ivx) = sm
        q1l(ibx) = ql(ibx,i)
        div = ql(idn,i) * slmv * slmm - bx2
        if ((sm .eq. ql(ivx,i)) .or. (div .eq. 0.0)                            &
                                .or. (bx2 .ge. gamma * ql(ipr,i))              &
                                .or. (sl .eq. (ql(ivx,i) + cl(i)))             &
                                .or. (sl .eq. (ql(ivx,i) - cl(i)))) then
          q1l(ivy) = ql(ivy,i)
          q1l(ivz) = ql(ivz,i)
          q1l(iby) = ql(iby,i)
          q1l(ibz) = ql(ibz,i)
        else
          fac = ql(ibx,i) * smvl / div
          q1l(ivy) = ql(ivy,i) - ql(iby,i) * fac
          q1l(ivz) = ql(ivz,i) - ql(ibz,i) * fac
          fac = (ql(idn,i) * slmv**2 - bx2) / div
          q1l(iby) = ql(iby,i) * fac
          q1l(ibz) = ql(ibz,i) * fac
        end if

! convert the left intermediate state to the conservative form
!
        u1l(idn) = q1l(idn)
        u1l(imx) = q1l(idn) * q1l(ivx)
        u1l(imy) = q1l(idn) * q1l(ivy)
        u1l(imz) = q1l(idn) * q1l(ivz)

        if (slmm .ne. 0.0) then
          u1l(ien) = (slmv * ul(ien,i) - ptl * ql(ivx,i) + pt * sm             &
                   + ql(ibx,i) * (sum(ql(ivx:ivz,i) * ql(ibx:ibz,i))           &
                   - sum(q1l(ivx:ivz) * q1l(ibx:ibz)))) / slmm
        else
          u1l(ien) = ul(ien,i)
        end if
        u1l(ibx) = q1l(ibx)
        u1l(iby) = q1l(iby)
        u1l(ibz) = q1l(ibz)
#ifdef GLM
        u1l(iph) = ul(iph,i)
#endif /* GLM */

! calculate the right intermediate state variables
!
        q1r(idn) = qr(idn,i) * srmv / srmm
        q1r(ivx) = sm
        q1r(ibx) = qr(ibx,i)
        div = qr(idn,i) * srmv * srmm - bx2
        if ((sm .eq. qr(ivx,i)) .or. (div .eq. 0.0)                            &
                                .or. (bx2 .ge. gamma * qr(ipr,i))              &
                                .or. (sr .eq. (qr(ivx,i) + cr(i)))             &
                                .or. (sr .eq. (qr(ivx,i) - cr(i)))) then
          q1r(ivy) = qr(ivy,i)
          q1r(ivz) = qr(ivz,i)
          q1r(iby) = qr(iby,i)
          q1r(ibz) = qr(ibz,i)
        else
          fac = qr(ibx,i) * smvr / div
          q1r(ivy) = qr(ivy,i) - qr(iby,i) * fac
          q1r(ivz) = qr(ivz,i) - qr(ibz,i) * fac
          fac = (qr(idn,i) * srmv**2 - bx2) / div
          q1r(iby) = qr(iby,i) * fac
          q1r(ibz) = qr(ibz,i) * fac
        end if

! convert the right intermediate state to the conservative form
!
        u1r(idn) = q1r(idn)
        u1r(imx) = q1r(idn) * q1r(ivx)
        u1r(imy) = q1r(idn) * q1r(ivy)
        u1r(imz) = q1r(idn) * q1r(ivz)

        if (srmm .ne. 0.0) then
          u1r(ien) = (srmv * ur(ien,i) - ptr * qr(ivx,i) + pt * sm             &
                   + qr(ibx,i) * (sum(qr(ivx:ivz,i) * qr(ibx:ibz,i))           &
                   - sum(q1r(ivx:ivz) * q1r(ibx:ibz)))) / srmm
        else
          u1r(ien) = ur(ien,i)
        end if
        u1r(ibx) = q1r(ibx)
        u1r(iby) = q1r(iby)
        u1r(ibz) = q1r(ibz)
#ifdef GLM
        u1r(iph) = ur(iph,i)
#endif /* GLM */

! Alfven speeds (eq. 51)
!
        sml = sm - abs(ql(ibx,i)) / sqrt(q1l(idn))
        smr = sm + abs(qr(ibx,i)) / sqrt(q1r(idn))

! intermediate discontinuities
!
        if (sml .ge. 0.0d0) then

! calculate the left intermediate flux
!
          fn(:,i) = fl(:,i) + sl * (u1l(:) - ul(:,i))

        else if (smr .le. 0.0d0) then

! calculate the right intermediate flux
!
          fn(:,i) = fr(:,i) + sr * (u1r(:) - ur(:,i))

        else ! sml < 0 & smr > 0

! obtain the normal component of magnetic field
!
          if (ql(ibx,i) .gt. 0.0d0) then
            bxs =  1.0d0
          else if (ql(ibx,i) .lt. 0.0d0) then
            bxs = -1.0d0
          else
            bxs =  0.0d0
          end if

! compute the density root squares
!
          dlsq = sqrt(q1l(idn))
          drsq = sqrt(q1r(idn))
          div  = dlsq + drsq

! calculate the velocity components
!
          q2(ivx) = sm
          q2(ivy) = (dlsq * q1l(ivy) + drsq * q1r(ivy)                         &
                                     + (q1r(iby) - q1l(iby)) * bxs) / div
          q2(ivz) = (dlsq * q1l(ivz) + drsq * q1r(ivz)                         &
                                     + (q1r(ibz) - q1l(ibz)) * bxs) / div

! calculate the magnetic field components
!
          q2(ibx) = ql(ibx,i)
          q2(iby) = (dlsq * q1r(iby) + drsq * q1l(iby)                         &
                       + dlsq * drsq * (q1r(ivy) - q1l(ivy)) * bxs) / div
          q2(ibz) = (dlsq * q1r(ibz) + drsq * q1l(ibz)                         &
                       + dlsq * drsq * (q1r(ivz) - q1l(ivz)) * bxs) / div

          if (sm .ge. 0.0) then

! convert the left Alfven intermediate state to the conservative form
!
            u2(idn) = u1l(idn)
            u2(imx) = u1l(idn) * q2(ivx)
            u2(imy) = u1l(idn) * q2(ivy)
            u2(imz) = u1l(idn) * q2(ivz)
            u2(ien) = u1l(ien) - dlsq * (sum(q1l(ivx:ivz) * q1l(ibx:ibz))      &
                                       - sum(q2 (ivx:ivz) * q2 (ibx:ibz))) * bxs
            u2(ibx) = u1l(ibx)
            u2(iby) =  q2(iby)
            u2(ibz) =  q2(ibz)
#ifdef GLM
            u2(iph) = u1l(iph)
#endif /* GLM */

! calculate the numerical flux
!
            fn(:,i) = fl(:,i) + sml * u2(:) - (sml - sl) * u1l(:) - sl * ul(:,i)

          else ! sm < 0

! convert the right Alfven intermediate state to the conservative form
!
            u2(idn) = u1r(idn)
            u2(imx) = u1r(idn) * q2(ivx)
            u2(imy) = u1r(idn) * q2(ivy)
            u2(imz) = u1r(idn) * q2(ivz)
            u2(ien) = u1r(ien) + drsq * (sum(q1r(ivx:ivz) * q1r(ibx:ibz))      &
                                       - sum(q2 (ivx:ivz) * q2 (ibx:ibz))) * bxs
            u2(ibx) = u1r(ibx)
            u2(iby) =  q2(iby)
            u2(ibz) =  q2(ibz)
#ifdef GLM
            u2(iph) = u1r(iph)
#endif /* GLM */

! calculate the numerical flux
!
            fn(:,i) = fr(:,i) + smr * u2(:) - (smr - sr) * u1r(:) - sr * ur(:,i)

          end if
        end if

      end if

    end do

#ifdef VISCOSITY
! add viscous term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dvx = vih * (q(ivx,ip1) - q(ivx,i))
      fn(ivx,i  ) = fn(ivx,i  ) - dvx

      dvy = vih * (q(ivy,ip1) - q(ivy,i))
      fn(ivy,i  ) = fn(ivy,i  ) - dvy

      dvz = vih * (q(ivz,ip1) - q(ivz,i))
      fn(ivz,i  ) = fn(ivz,i  ) - dvz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - 0.5d0 * (ql(ivx,i) + qr(ivx,i)) * dvx
#endif /* ADI */
    end do

#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
! add resistivity term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dbx = ueh * (q(ibx,ip1) - q(ibx,i))
      fn(ibx,i) = fn(ibx,i) - dbx

      dby = ueh * (q(iby,ip1) - q(iby,i))
      fn(iby,i) = fn(iby,i) - dby

      dbz = ueh * (q(ibz,ip1) - q(ibz,i))
      fn(ibz,i) = fn(ibz,i) - dbz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - ql(ibx,i) * dbx
#endif /* ADI */
    end do

#endif /* MHD & RESISTIVITY */
! calculate the numerical flux derivative
!
    f(  1:nfl,2:n) = - fn(  1:nfl,2:n) + fn(  1:nfl,1:n-1)
    f(ibx:ibz,2:n) = - fn(ibx:ibz,2:n) + fn(ibx:ibz,1:n-1)
#ifdef GLM
    f(iph    ,2:n) = - fn(iph    ,2:n) + fn(iph    ,1:n-1)
#endif /* GLM */

!-------------------------------------------------------------------------------
!
  end subroutine hlld
#endif /* ADI */
#endif /* HLLD */
#endif /* MHD */
#ifdef ROE
!
!===============================================================================
!
! roe: subroutine computes the approximated flux using the ROE method
!
! references:
!
!   - Roe, P. L., 1981, Journal of Computational Physics, 43, 357
!
!===============================================================================
!
  subroutine roe(n, h, u, f)

    use config       , only : gamma
#ifdef VISCOSITY
    use config       , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config       , only : ueta
#endif /* MHD & RESISTIVITY */
    use interpolation, only : reconstruct
    use variables    , only : nvr, nfl, nqt
    use variables    , only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables    , only : ien, ipr
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input/output arguments
!
    integer               , intent(in)  :: n
    real                  , intent(in)  :: h
    real, dimension(nvr,n), intent(in)  :: u
    real, dimension(nqt,n), intent(out) :: f

! local variables
!
    integer                  :: p, i, ip1
    real, dimension(nvr,n)   :: q, ql, qr, ul, ur
    real, dimension(nqt,n)   :: fl, fr, fn
    real, dimension(n)       :: cl, cr
    real, dimension(nqt)     :: qi, ci, et, du
    real, dimension(nqt,nqt) :: li, ri
    real                     :: al, ar, ap, div
    real                     :: sdl, sdr, sds, sfl, sfr
#ifdef MHD
    real                     :: pbl, pbr, xfc, yfc
#endif /* MHD */
#ifdef VISCOSITY
    real                     :: dvx, dvy, dvz, vih
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    real                     :: dbx, dby, dbz, ueh
#endif /* MHD & RESISTIVITY */
!
!-------------------------------------------------------------------------------
!
! usefull parameters
!
#ifdef VISCOSITY
    vih = visc / h
#endif /* VISCOSITY */
#ifdef RESISTIVITY
    ueh = ueta / h
#endif /* RESISTIVITY */

! reset eigensystem values
!
    ci(:)   = 0.0d0
    li(:,:) = 0.0d0
    ri(:,:) = 0.0d0

! calculate the primitive variables
!
    call cons2prim(n, u(:,:), q(:,:))

! reconstruct the left and right states of the primitive variables
!
    do p = 1, nfl
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef MHD
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
#endif /* MHD */

! calculate conservative variables at states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate fluxes and speeds
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
      sdl = sqrt(ql(idn,i))
      sdr = sqrt(qr(idn,i))
      sds = sdl + sdr
      sfl     = sdl / sds
      sfr     = sdr / sds

      qi(idn) = sdl * sdr
      qi(ivx) = sfl * ql(ivx,i) + sfr * qr(ivx,i)
      qi(ivy) = sfl * ql(ivy,i) + sfr * qr(ivy,i)
      qi(ivz) = sfl * ql(ivz,i) + sfr * qr(ivz,i)
#ifdef HYDRO
#ifdef ADI
      qi(ipr) = sfl * (ul(ien,i) + ql(ipr,i)) / ql(idn,i)                      &
              + sfr * (ur(ien,i) + qr(ipr,i)) / qr(idn,i)
#endif /* ADI */
#endif /* HYDRO */
#ifdef MHD
      qi(ibx) = ql(ibx,i)
      qi(iby) = sfl * ql(iby,i) + sfr * qr(iby,i)
      qi(ibz) = sfl * ql(ibz,i) + sfr * qr(ibz,i)
#ifdef ADI
      pbl     = 0.5d0 * sum(ql(ibx:ibz,i)**2)
      pbr     = 0.5d0 * sum(qr(ibx:ibz,i)**2)
      qi(ipr) = sfl * (ul(ien,i) + ql(ipr,i) + pbl) / ql(idn,i)                &
              + sfr * (ur(ien,i) + qr(ipr,i) + pbr) / qr(idn,i)
#endif /* ADI */
      xfc     = 0.5d0 * sum(du(iby:ibz)**2) / sds**2
      yfc     = 0.5d0 * (ql(idn,i) + qr(idn,i)) / qi(idn)
#endif /* MHD */

! check if density and pressure are positive
!
#ifdef ADI
      if (qi(idn) .gt. 0.0d0 .and. qi(ipr) .gt. 0.0d0) then
#else /* ADI */
      if (qi(idn) .gt. 0.0d0) then
#endif /* ADI */

! obtain eigenvalues and eigenvectors
!
#ifdef HYDRO
        call eigensystem(qi(:), ci(:), ri(:,:), li(:,:))
#endif /* HYDRO */
#ifdef MHD
        call eigensystem(qi(:), ci(:), ri(:,:), li(:,:), xfc, yfc)
#endif /* MHD */

! calculate vector (Ur - Ul).L
!
        et(:) = 0.0d0
        do p = 1, nqt
          et(:) = et(:) + du(p) * li(p,:)
        end do

! calculate numerical flux
!
        fn(:,i) = 0.5d0 * (fl(:,i) + fr(:,i))

! add term abs(lambda) R.(Ur - Ul).L
!
        do p = 1, nqt
          fn(:,i) = fn(:,i) - 0.5d0 * abs(ci(p)) * et(p) * ri(p,:)
        end do

      else ! in the case when density or pressure of the intermediate state
           ! are negative, use the simplest and most robust HLL flux

! calculate min and max speeds
!
        al = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
        ar = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
        if (al .ge. 0.0d0) then
          fn(:,i) = fl(:,i)
        else if (ar .le. 0.0d0) then
          fn(:,i) = fr(:,i)
        else
          ap  = ar * al
          div = 1.0d0 / (ar - al)

          fn(:,i) = div * (ar * fl(:,i) - al * fr(:,i)                         &
                                                   + ap * (ur(:,i) - ul(:,i)))
        end if

      end if

    end do

#ifdef VISCOSITY
! add viscous term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dvx = vih * (q(ivx,ip1) - q(ivx,i))
      fn(ivx,i  ) = fn(ivx,i  ) - dvx

      dvy = vih * (q(ivy,ip1) - q(ivy,i))
      fn(ivy,i  ) = fn(ivy,i  ) - dvy

      dvz = vih * (q(ivz,ip1) - q(ivz,i))
      fn(ivz,i  ) = fn(ivz,i  ) - dvz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - 0.5d0 * (ql(ivx,i) + qr(ivx,i)) * dvx
#endif /* ADI */
    end do

#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
! add resistivity term to the left and right fluxes
!
    do i = 1, n - 1
      ip1 = i + 1

      dbx = ueh * (q(ibx,ip1) - q(ibx,i))
      fn(ibx,i) = fn(ibx,i) - dbx

      dby = ueh * (q(iby,ip1) - q(iby,i))
      fn(iby,i) = fn(iby,i) - dby

      dbz = ueh * (q(ibz,ip1) - q(ibz,i))
      fn(ibz,i) = fn(ibz,i) - dbz
#ifdef ADI
      fn(ien,i) = fn(ien,i) - ql(ibx,i) * dbx
#endif /* ADI */
    end do

#endif /* MHD & RESISTIVITY */
! calculate numerical flux
!
    f(  1:nfl,2:n) = - fn(  1:nfl,2:n) + fn(  1:nfl,1:n-1)
#ifdef MHD
#ifdef GLM
    f(ibx:ibz,2:n) = - fn(ibx:ibz,2:n) + fn(ibx:ibz,1:n-1)
    f(iph    ,2:n) = - fn(iph    ,2:n) + fn(iph    ,1:n-1)
#endif /* GLM */
#endif /* MHD */

!-------------------------------------------------------------------------------
!
  end subroutine roe
!
!===============================================================================
!
! eigensystem: subroutine computes eigenvalues and eigenmatrices for a given
!              set of equations and input variables
!
!===============================================================================
!
#ifdef HYDRO
#ifdef ADI
  subroutine eigensystem(q, c, r, l)

    use config   , only : gamma
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

    use config   , only : csnd
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
#endif /* HYDRO */
#ifdef MHD
#ifdef ADI
  subroutine eigensystem(q, c, r, l, x, y)

    use config   , only : gamma
    use variables, only : nqt
    use variables, only : idn, ivx, ivy, ivz, ibx, iby, ibz, ien

    implicit none

! input/output arguments
!
    real, dimension(nqt)    , intent(in)    :: q
    real, dimension(nqt)    , intent(inout) :: c
    real, dimension(nqt,nqt), intent(inout) :: l, r
    real                    , intent(in)    :: x, y

! local variables
!
    real :: gm1, gm2, v2, b2
    real :: ah2, bt2, aa2, cb2, cx2, ct2, ca2, cs2, cf2
    real :: cx , ct , cs, ca, cf
!
!-------------------------------------------------------------------------------
!
! calculate characteristic speeds
!
    gm1 = gamma - 1.0d0
    gm2 = gamma - 2.0d0
    v2  = sum(q(ivx:ivz)**2)
    b2  = sum(q(ibx:ibz)**2)
    ah2 = gm1 * (q(ien) - 0.5d0 * v2 - b2 / q(idn)) - gm2 * x
    bt2 = (gm1 - gm2 * y) * sum(q(iby:ibz)**2)
    cx2 = q(ibx) * q(ibx) / q(idn)
    ct2 = bt2 / q(idn)
    ca2 = cx2 + ct2
    aa2 = ah2 + ca2
    cb2 = sqrt(max(0.0d0, aa2 * aa2 - 4.0d0 * ah2 * cx2))
    cs2 = 0.5d0 * (aa2 - cb2)
    cf2 = 0.5d0 * (aa2 + cb2)
    af2 = (ah2  - cs2) / (cf2 - cs2)
    cx  = sqrt(cx2)
    ca  = sqrt(ca2)
    cs  = sqrt(cs2)
    cf  = sqrt(cf2)
    af  = sqrt(af2)

! prepare eigenvalues
!
    c(1) = q(ivx) - cf
    c(2) = q(ivx) - cx
    c(3) = q(ivx) - cs
    c(4) = q(ivx)
    c(5) = q(ivx) + cs
    c(6) = q(ivx) + cx
    c(7) = q(ivx) + cf

! prepare the right eigenmatrix
!
    r(1,idn) = af
    r(1,ivx) = af * (q(ivx) - cf)

! prepare the left eigenmatrix
!

!-------------------------------------------------------------------------------
!
  end subroutine eigensystem
#endif /* ADI */
#ifdef ISO
  subroutine eigensystem(q, c, r, l, x, y)

    use config   , only : csnd
    use variables, only : nqt
    use variables, only : idn, ivx, ivy, ivz, ibx, iby, ibz

    implicit none

! input/output arguments
!
    real, dimension(nqt)    , intent(in)    :: q
    real, dimension(nqt)    , intent(inout) :: c
    real, dimension(nqt,nqt), intent(inout) :: l, r
    real                    , intent(in)    :: x, y

! local variables
!
    real :: cs, ca, cx, cf
!
!-------------------------------------------------------------------------------
!
! calculate characteristic speeds
!
    ca = sqrt(sum(q(ibx:ibz)**2) / q(idn))
    cx = abs(q(ibx)) / sqrt(q(idn))

! prepare eigenvalues
!
    c(1) = q(ivx) - cf
    c(2) = q(ivx) - ca
    c(3) = q(ivx) - cs
    c(4) = q(ivx) + cs
    c(5) = q(ivx) + ca
    c(6) = q(ivx) + cf

! prepare the right eigenmatrix
!

! prepare the left eigenmatrix
!

!-------------------------------------------------------------------------------
!
  end subroutine eigensystem
#endif /* ISO */
#endif /* MHD */
#endif /* ROE */
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
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
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
      f(ibx,i) = 0.0d0
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
#ifdef GLM
      f(ibx,i) = q(iph,i)
      f(iph,i) = q(ibx,i)
#endif /* GLM */
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
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
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
      em       = 0.5 * sum(u(ibx:ibz,i) * u(ibx:ibz,i))
      ei       = ei - em
#endif /* MHD */
      q(ipr,i) = gammam1 * ei
#endif /* ADI */
#ifdef MHD
      q(ibx,i) = u(ibx,i)
      q(iby,i) = u(iby,i)
      q(ibz,i) = u(ibz,i)
#ifdef GLM
      q(iph,i) = u(iph,i)
#endif /* GLM */
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
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
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
      em       = 0.5 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      u(ien,i) = u(ien,i) + em
#endif /* ADI */
      u(ibx,i) = q(ibx,i)
      u(iby,i) = q(iby,i)
      u(ibz,i) = q(ibz,i)
#ifdef GLM
      u(iph,i) = q(iph,i)
#endif /* GLM */
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
    use variables    , only : nvr, nqt
    use variables    , only : idn, ivx, ivz
#ifdef ADI
    use variables    , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
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
    real, dimension(nvr,im)       :: q
!
!-------------------------------------------------------------------------------
!
    maxspeed = 0.0

! iterate over all points and calculate maximum speed
!
    do k = kb, ke
      do j = jb, je

        call cons2prim(im, u(1:nqt,1:im,j,k), q(1:nqt,1:im))

        do i = ib, ie

! calculate the velocity
!
          vv = sum(q(ivx:ivz,i)**2)
          v  = sqrt(vv)
#ifdef MHD
          bb = sum(q(ibx:ibz,i)**2)
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
