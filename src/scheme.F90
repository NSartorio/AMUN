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

    use blocks, only : nv => nvars, idn, imx, imy, imz, ien
    use config, only : im, jm, km, ngrids

    implicit none

! input arguments
!
    real, dimension(nv,im,jm,km), intent(in)  :: u
    real, dimension(nv,im,jm,km), intent(out) :: du
    real                        , intent(in)  :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k

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
      do j = 1, jm

! copy directional vectors of variables for the one dimensional solver
!
        do i = 1, im
          ul(1,i) = u(idn,i,j,k)
          ul(2,i) = u(imx,i,j,k)
          ul(3,i) = u(imy,i,j,k)
          ul(4,i) = u(imz,i,j,k)
#ifdef ADI
          ul(5,i) = u(ien,i,j,k)
#endif /* ADI */
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
          du(idn,i,j,k) = du(idn,i,j,k) + dxi * fl(1,i)
          du(imx,i,j,k) = du(imx,i,j,k) + dxi * fl(2,i)
          du(imy,i,j,k) = du(imy,i,j,k) + dxi * fl(3,i)
          du(imz,i,j,k) = du(imz,i,j,k) + dxi * fl(4,i)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dxi * fl(5,i)
#endif /* ADI */
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
          ul(1,j) = u(idn,i,j,k)
          ul(2,j) = u(imy,i,j,k)
          ul(3,j) = u(imz,i,j,k)
          ul(4,j) = u(imx,i,j,k)
#ifdef ADI
          ul(5,j) = u(ien,i,j,k)
#endif /* ADI */
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
          du(idn,i,j,k) = du(idn,i,j,k) + dyi * fl(1,j)
          du(imx,i,j,k) = du(imx,i,j,k) + dyi * fl(4,j)
          du(imy,i,j,k) = du(imy,i,j,k) + dyi * fl(2,j)
          du(imz,i,j,k) = du(imz,i,j,k) + dyi * fl(3,j)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dyi * fl(5,j)
#endif /* ADI */
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
          ul(1,k) = u(idn,i,j,k)
          ul(2,k) = u(imz,i,j,k)
          ul(3,k) = u(imx,i,j,k)
          ul(4,k) = u(imy,i,j,k)
#ifdef ADI
          ul(5,k) = u(ien,i,j,k)
#endif /* ADI */
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
          du(idn,i,j,k) = du(idn,i,j,k) + dzi * fl(1,k)
          du(imx,i,j,k) = du(imx,i,j,k) + dzi * fl(3,k)
          du(imy,i,j,k) = du(imy,i,j,k) + dzi * fl(4,k)
          du(imz,i,j,k) = du(imz,i,j,k) + dzi * fl(2,k)
#ifdef ADI
          du(ien,i,j,k) = du(ien,i,j,k) + dzi * fl(5,k)
#endif /* ADI */
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
  subroutine hll(m, n, uc, f)

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
    real, dimension(m,n) :: ul, ur, ql, qr, qc
    real, dimension(m,n) :: fl, fr, fx
    real, dimension(n)   :: cl, cr
    real                 :: al, ar, ap, div
!
!-------------------------------------------------------------------------------
!
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

! iterate over all points
!
    do i = 1, n

! calculate min and max and intermediate speeds: eq. (67)
!
      al = min(ql(2,i) - cl(i),qr(2,i) - cr(i))
      ar = max(ql(2,i) + cl(i),qr(2,i) + cr(i))

! calculate HLL flux
!
      if (al .ge. 0.0) then
        fx(:,i) = fl(:,i)
      else if (ar .le. 0.0) then
        fx(:,i) = fr(:,i)
      else
        ap  = ar * al
        div = 1.0/(ar - al)

        fx(:,i) = div*(ar*fl(:,i) - al*fr(:,i) + ap*(ur(:,i) - ul(:,i)))
      endif

    enddo

! calculate numerical flux
!
    f(:,2:n) = - fx(:,2:n) + fx(:,1:n-1)

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
    real    :: cs
!
!-------------------------------------------------------------------------------
!
! sweep over all points
!
    do i = 1, n

! compute fluxes
!
      f(1,i) = u(2,i)
#ifdef ADI
      f(2,i) = q(2,i)*u(2,i) + q(5,i)
#endif /* ADI */
#ifdef ISO
      f(2,i) = q(2,i)*u(2,i) + q(1,i)*csnd2
#endif /* ISO */
      f(3,i) = q(2,i)*u(3,i)
      f(4,i) = q(2,i)*u(4,i)
#ifdef ADI
      f(5,i) = q(2,i)*(u(5,i) + q(5,i))
#endif /* ADI */

! compute speeds
!
#ifdef ADI
      c(i) = sqrt(gamma*q(5,i)/q(1,i))
#endif /* ADI */
#ifdef ISO
      c(i) = csnd
#endif /* ISO */
    enddo

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

    use config   , only : gammam1

    implicit none

! input/output arguments
!
    integer             , intent(in)  :: m, n
    real, dimension(m,n), intent(in)  :: u
    real, dimension(m,n), intent(out) :: q

! local variables
!
    integer :: i
    real    :: dni
!
!-------------------------------------------------------------------------------
!
    do i = 1, n
      dni    = 1.0 / u(1,i)

      q(1,i) =     u(1,i)
      q(2,i) = dni*u(2,i)
      q(3,i) = dni*u(3,i)
      q(4,i) = dni*u(4,i)
#ifdef ADI
      q(5,i) = u(5,i) &
             - 0.5*(u(2,i)*q(2,i) + u(3,i)*q(3,i) + u(4,i)*q(4,i))
      q(5,i) = gammam1*q(5,i)
#endif /* ADI */
    enddo

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
!
!-------------------------------------------------------------------------------
!
    do i = 1, n
      u(1,i) = q(1,i)
      u(2,i) = q(1,i)*q(2,i)
      u(3,i) = q(1,i)*q(3,i)
      u(4,i) = q(1,i)*q(4,i)
#ifdef ADI
      u(5,i) = gammam1i*q(5,i) &
             + 0.5*(u(2,i)*q(2,i) + u(3,i)*q(3,i) + u(4,i)*q(4,i))
#endif /* ADI */
    enddo

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

    use blocks, only : nv => nvars, idn, ivx, ivy, ivz, ipr
    use config, only : im, jm, km, ib, ie, jb, je, kb, ke, gamma

    implicit none

! input arguments
!
    real, dimension(nv,im,jm,km), intent(in)  :: u

! local variables
!
    integer :: i, j, k
    real    :: vv, v, c
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

        call cons2prim(nv,im,u(:,:,j,k),q(:,:))

        do i = ib, ie

! calculate the velocity
!
          vv = sum(q(ivx:ivz,i)**2)
          v  = sqrt(vv)

! calculate the maximum characteristic speed
!
#ifdef ADI
          c = sqrt(gamma*q(ipr,i)/q(idn,i))
#endif /* ADI */
#ifdef ISO
          c = csnd
#endif /* ISO */

! calculate maximum of the speed
!
          maxspeed = max(maxspeed, v + c)
        enddo
      enddo
    enddo

!-------------------------------------------------------------------------------
!
  end function maxspeed

!===============================================================================
!
end module
