!!*****************************************************************************
!!
!! module: scheme - handling the actual solver of the set of equations
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
module scheme

  implicit none

  contains
!
!==============================================================================
!
! update: subroutine sweeps over all directions and integrates the directional
!         derivatives of the flux in order to get the increment of solution
!
!==============================================================================
!
  subroutine update(u, du, dxi, dyi, dzi)

    use blocks, only : nv => nvars
    use blocks, only : idn, imx, imy, imz, ien, ibx, iby, ibz, icx, icy, icz
    use config, only : im, jm, km, ngrids

    implicit none

! input arguments
!
    real, dimension(nv,im,jm,km), intent(in)  :: u
    real, dimension(nv,im,jm,km), intent(out) :: du
    real                        , intent(in)  :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k, im1, jm1, km1, ip1, jp1, kp1

! local arrays
!
    real, dimension(nv,ngrids) :: ul, fl
!
!-------------------------------------------------------------------------------
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
          ul(idn,i) = u(idn,i,j,k)
          ul(imx,i) = u(imx,i,j,k)
          ul(imy,i) = u(imy,i,j,k)
          ul(imz,i) = u(imz,i,j,k)
#ifdef ADI
          ul(ien,i) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          ul(ibx,i) = u(ibx,i,j,k)
          ul(iby,i) = u(iby,i,j,k)
          ul(ibz,i) = u(ibz,i,j,k)
          ul(icx,i) = u(icx,i,j,k)
          ul(icy,i) = u(icy,i,j,k)
          ul(icz,i) = u(icz,i,j,k)
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll(nv, im, ul, fl)
#endif /* HLL */
#ifdef HLLC
        call hllc(nv, im, ul, fl)
#endif /* HLLC */

! update the arrays of increments
!
        do i = 1, im
          du(idn,i,j,k) = du(idn,i,j,k) + dxi * fl(idn,i)
          du(imx,i,j,k) = du(imx,i,j,k) + dxi * fl(imx,i)
          du(imy,i,j,k) = du(imy,i,j,k) + dxi * fl(imy,i)
          du(imz,i,j,k) = du(imz,i,j,k) + dxi * fl(imz,i)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dxi * fl(ien,i)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          fl(ibx:ibz,i) = 0.25 * fl(ibx:ibz,i)

          im1 = max(i - 1, 1)

          du(ibx,i,jm1,k  ) = du(ibx,i,jm1,k  ) + dyi *  fl(iby,i)
          du(ibx,i,jp1,k  ) = du(ibx,i,jp1,k  ) - dyi *  fl(iby,i)
          du(iby,i,j  ,k  ) = du(iby,i,j  ,k  ) - dxi * (fl(iby,i) - fl(iby,im1))
          du(iby,i,jm1,k  ) = du(iby,i,jm1,k  ) - dxi * (fl(iby,i) - fl(iby,im1))
#if NDIMS == 3
          du(ibx,i,j  ,km1) = du(ibx,i,j  ,km1) + dzi *  fl(ibz,i)
          du(ibx,i,j  ,kp1) = du(ibx,i,j  ,kp1) - dzi *  fl(ibz,i)
          du(ibz,i,j  ,k  ) = du(ibz,i,j  ,k  ) - dxi * (fl(ibz,i) - fl(ibz,im1))
          du(ibz,i,j  ,km1) = du(ibz,i,j  ,km1) - dxi * (fl(ibz,i) - fl(ibz,im1))
#endif /* NDIMS == 3 */
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
          ul(idn,j) = u(idn,i,j,k)
          ul(imx,j) = u(imy,i,j,k)
          ul(imy,j) = u(imz,i,j,k)
          ul(imz,j) = u(imx,i,j,k)
#ifdef ADI
          ul(ien,j) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          ul(ibx,j) = u(iby,i,j,k)
          ul(iby,j) = u(ibz,i,j,k)
          ul(ibz,j) = u(ibx,i,j,k)
          ul(icx,j) = u(icy,i,j,k)
          ul(icy,j) = u(icz,i,j,k)
          ul(icz,j) = u(icx,i,j,k)
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll(nv, jm, ul, fl)
#endif /* HLL */
#ifdef HLLC
        call hllc(nv, jm, ul, fl)
#endif /* HLLC */

! update the arrays of increments
!
        do j = 1, jm
          du(idn,i,j,k) = du(idn,i,j,k) + dyi * fl(idn,j)
          du(imx,i,j,k) = du(imx,i,j,k) + dyi * fl(imz,j)
          du(imy,i,j,k) = du(imy,i,j,k) + dyi * fl(imx,j)
          du(imz,i,j,k) = du(imz,i,j,k) + dyi * fl(imy,j)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dyi * fl(ien,j)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          fl(ibx:ibz,j) = 0.25 * fl(ibx:ibz,j)

          jm1 = max(j - 1, 1)
#if NDIMS == 3
          du(iby,i  ,j,km1) = du(iby,i  ,j,km1) + dzi *  fl(iby,j)
          du(iby,i  ,j,kp1) = du(iby,i  ,j,kp1) - dzi *  fl(iby,j)
          du(ibz,i  ,j,k  ) = du(ibz,i  ,j,k  ) - dyi * (fl(iby,j) - fl(iby,jm1))
          du(ibz,i  ,j,km1) = du(ibz,i  ,j,km1) - dyi * (fl(iby,j) - fl(iby,jm1))
#endif /* NDIMS == 3 */
          du(iby,im1,j,k  ) = du(iby,im1,j,k  ) + dxi *  fl(ibz,j)
          du(iby,ip1,j,k  ) = du(iby,ip1,j,k  ) - dxi *  fl(ibz,j)
          du(ibx,i  ,j,k  ) = du(ibx,i  ,j,k  ) - dyi * (fl(ibz,j) - fl(ibz,jm1))
          du(ibx,im1,j,k  ) = du(ibx,im1,j,k  ) - dyi * (fl(ibz,j) - fl(ibz,jm1))
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
          ul(idn,k) = u(idn,i,j,k)
          ul(imx,k) = u(imz,i,j,k)
          ul(imy,k) = u(imx,i,j,k)
          ul(imz,k) = u(imy,i,j,k)
#ifdef ADI
          ul(ien,k) = u(ien,i,j,k)
#endif /* ADI */
#ifdef MHD
          ul(ibx,k) = u(ibz,i,j,k)
          ul(iby,k) = u(ibx,i,j,k)
          ul(ibz,k) = u(iby,i,j,k)
          ul(icx,k) = u(icz,i,j,k)
          ul(icy,k) = u(icx,i,j,k)
          ul(icz,k) = u(icy,i,j,k)
#endif /* MHD */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll(nv, km, ul, fl)
#endif /* HLL */
#ifdef HLLC
        call hllc(nv, km, ul, fl)
#endif /* HLLC */

! update the arrays of increments
!
        do k = 1, km
          du(idn,i,j,k) = du(idn,i,j,k) + dzi * fl(idn,k)
          du(imx,i,j,k) = du(imx,i,j,k) + dzi * fl(imy,k)
          du(imy,i,j,k) = du(imy,i,j,k) + dzi * fl(imz,k)
          du(imz,i,j,k) = du(imz,i,j,k) + dzi * fl(imx,k)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dzi * fl(ien,k)
#endif /* ADI */
#ifdef MHD
! update magnetic variables
!
          fl(ibx:ibz,k) = 0.25 * fl(ibx:ibz,k)

          km1 = max(k - 1, 1)

          du(ibz,im1,j  ,k) = du(ibz,im1,j  ,k) + dxi *  fl(iby,k)
          du(ibz,ip1,j  ,k) = du(ibz,ip1,j  ,k) - dxi *  fl(iby,k)
          du(ibx,i  ,j  ,k) = du(ibx,i  ,j  ,k) - dzi * (fl(iby,k) - fl(iby,km1))
          du(ibx,im1,j  ,k) = du(ibx,im1,j  ,k) - dzi * (fl(iby,k) - fl(iby,km1))

          du(ibz,i  ,jm1,k) = du(ibz,i  ,jm1,k) + dyi *  fl(ibz,k)
          du(ibz,i  ,jp1,k) = du(ibz,i  ,jp1,k) - dyi *  fl(ibz,k)
          du(iby,i  ,j  ,k) = du(iby,i  ,j  ,k) - dzi * (fl(ibz,k) - fl(ibz,km1))
          du(iby,i  ,jm1,k) = du(iby,i  ,jm1,k) - dzi * (fl(ibz,k) - fl(ibz,km1))
#endif /* MHD */
        end do
      end do
    end do
#endif /* NDIMS == 3 */

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
  subroutine hll(m, n, u, f)

    use blocks       , only : idn, imx, imy, imz, ivx, ivy, ivz, ipr, ien      &
                            , ibx, iby, ibz, icx, icy, icz, ifl, iqt
    use interpolation, only : reconstruct

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: u
    real, dimension(m,n), intent(out) :: f

! local variables
!
    integer              :: p, i
    real, dimension(m,n) :: ul, ur, ql, qr, q
    real, dimension(m,n) :: fl, fr, fx
    real, dimension(n)   :: cl, cr
    real                 :: al, ar, ap, div
!
!-------------------------------------------------------------------------------
!
#ifdef CONREC
! reconstruct left and right states of conserved variables
!
    do p = 1, ifl
      call reconstruct(n, u(p,:), ul(p,:), ur(p,:))
    end do

#ifdef MHD
! reconstruct left and right states of magnetic field components
!
    ul(ibx,:) = u(ibx,:)
    ur(ibx,:) = u(ibx,:)
    call reconstruct(n, u(icy,:), ul(iby,:), ur(iby,:))
    call reconstruct(n, u(icz,:), ul(ibz,:), ur(ibz,:))
    ul(icx:icz,:) = ul(ibx:ibz,:)
    ur(icx:icz,:) = ur(ibx:ibz,:)
#endif /* MHD */

! calculate primitive variables
!
    call cons2prim(m, n, ul, ql)
    call cons2prim(m, n, ur, qr)
#else /* CONREC */

! calculate primitive variables
!
    call cons2prim(m, n, u, q)

! reconstruct left and right states of primitive variables
!
    do p = 1, ifl
      call reconstruct(n, q(p,:), ql(p,:), qr(p,:))
    end do

#ifdef MHD
! reconstruct left and right states of magnetic field components
!
    ql(ibx,:) = q(ibx,:)
    qr(ibx,:) = q(ibx,:)
    call reconstruct(n, q(icy,:), ql(iby,:), qr(iby,:))
    call reconstruct(n, q(icz,:), ql(ibz,:), qr(ibz,:))
    ql(icx:icz,:) = ql(ibx:ibz,:)
    qr(icx:icz,:) = qr(ibx:ibz,:)
#endif /* MHD */

! calculate conservative variables at states
!
    call prim2cons(m, n, ql, ul)
    call prim2cons(m, n, qr, ur)
#endif /* CONREC */

! calculate fluxes and speeds
!
    call fluxspeed(m, n, ql, ul, fl, cl)
    call fluxspeed(m, n, qr, ur, fr, cr)

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
        fx(1:iqt,i) = fl(1:iqt,i)
      else if (ar .le. 0.0) then
        fx(1:iqt,i) = fr(1:iqt,i)
      else
        ap  = ar * al
        div = 1.0 / (ar - al)

        fx(1:iqt,i) = div * (ar * fl(1:iqt,i) - al * fr(1:iqt,i) + ap * (ur(1:iqt,i) - ul(1:iqt,i)))
      end if

    end do

! calculate numerical flux
!
    f(  1:ifl,2:n) = - fx(  1:ifl,2:n) + fx(   1:ifl,1:n-1)
#ifdef MHD
    f(ibx:ibz,1:n) =   fx(ibx:ibz,1:n)
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
  subroutine hllc(m, n, uc, f)

    use interpolation, only : reconstruct

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: uc
    real, dimension(m,n), intent(out) :: f

! local variables
!
    integer              :: p, i
    real, dimension(m,n) :: ql, qr, qc, ul, ur
    real, dimension(m,n) :: fl, fr, fx
    real, dimension(n)   :: cl, cr, cm
    real                 :: sl, sr, sm, sml, smr, srmv, slmv, srmm, slmm  &
                          , smvl, smvr, div, pt
    real, dimension(m)   :: q1l, q1r, u1l, u1r
!
!----------------------------------------------------------------------
!
    f (:,:) = 0.0
    fx(:,:) = 0.0

#ifdef CONREC
! reconstruct left and right states of conserved variables
!
    do p = 1, m
      call reconstruct(n,uc(p,:),ul(p,:),ur(p,:))
    enddo

! calculate primitive variables
!
    call cons2prim(m,n,ul,ql)
    call cons2prim(m,n,ur,qr)
#else /* CONREC */

! calculate primitive variables
!
    call cons2prim(m,n,uc,qc)

! reconstruct left and right states of primitive variables
!
    do p = 1, m
      call reconstruct(n,qc(p,:),ql(p,:),qr(p,:))
    enddo

! calculate conservative variables at states
!
    call prim2cons(m,n,ql,ul)
    call prim2cons(m,n,qr,ur)
#endif /* CONREC */

! calculate fluxes and speeds
!
    call fluxspeed(m,n,ql,ul,fl,cl)
    call fluxspeed(m,n,qr,ur,fr,cr)

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
    f(:,2:n) = - fx(:,2:n) + fx(:,1:n-1)

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
  subroutine fluxspeed(m, n, q, u, f, c)

    use blocks, only : idn, imx, imy, imz, ivx, ivy, ivz, ipr, ien, ibx, iby, ibz
    use config, only : gamma, csnd, csnd2

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: q, u
    real, dimension(m,n), intent(out) :: f
    real, dimension(n)  , intent(out) :: c

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
      f(ibx,i) = 0.0
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
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
  subroutine cons2prim(m, n, u, q)

    use blocks, only : idn, imx, imy, imz, ivx, ivy, ivz, ipr, ien, icx, icy, icz, ibx, iby, ibz
    use config, only : gammam1

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: u
    real, dimension(m,n), intent(out) :: q

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
      q(icx,i) = u(icx,i)
      q(icy,i) = u(icy,i)
      q(icz,i) = u(icz,i)
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
  subroutine prim2cons(m, n, q, u)

    use blocks, only : idn, imx, imy, imz, ivx, ivy, ivz, ipr, ien, icx, icy, icz, ibx, iby, ibz
    use config, only : gammam1i

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: q
    real, dimension(m,n), intent(out) :: u

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
      u(icx,i) = q(icx,i)
      u(icy,i) = q(icy,i)
      u(icz,i) = q(icz,i)
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

    use blocks, only : nv => nvars, idn, ivx, ivz, ipr, icx, icz
    use config, only : im, jm, km, ib, ie, jb, je, kb, ke, gamma

    implicit none

! input arguments
!
    real, dimension(nv,im,jm,km), intent(in)  :: u

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
    real, dimension(nv,im) :: q
!
!----------------------------------------------------------------------
!
    maxspeed = 0.0

! iterate over all points and calculate maximum speed
!
    do k = kb, ke
      do j = jb, je

        call cons2prim(nv, im, u(:,:,j,k), q(:,:))

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

!-------------------------------------------------------------------------------
!
  end function maxspeed

!===============================================================================
!
end module
