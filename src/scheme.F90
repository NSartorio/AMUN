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

    use blocks, only : nvars, idn, imx, imy, imz, ien
    use config, only : igrids, jgrids, kgrids, ngrids

    implicit none

! input arguments
!
    real, dimension(nvars,igrids,jgrids,kgrids), intent(in)  :: u
    real, dimension(nvars,igrids,jgrids,kgrids), intent(out) :: du
    real                                       , intent(in)  :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k

! local arrays
!
    real, dimension(nvars,ngrids) :: ul, fl
!
!-------------------------------------------------------------------------------
!
    du(:,:,:,:) = 0.0

!  update along X-direction
!
    do k = 1, kgrids
      do j = 1, jgrids

! copy directional vectors of variables for the one dimensional solver
!
        do i = 1, igrids
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
        call hll(nvars, igrids, ul, fl)
#endif /* HLL */

! update the arrays of increments
!
        do i = 1, igrids
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
    do k = 1, kgrids
      do i = 1, igrids

! copy directional vectors of variables for the one dimensional solver
!
        do j = 1, jgrids
          ul(1,j) = u(idn,i,j,k)
          ul(2,i) = u(imy,i,j,k)
          ul(3,j) = u(imz,i,j,k)
          ul(4,j) = u(imx,i,j,k)
#ifdef ADI
          ul(5,j) = u(ien,i,j,k)
#endif /* ADI */
        end do

! execute solver (returns fluxes for the update)
!
#ifdef HLL
        call hll(nvars, jgrids, ul, fl)
#endif /* HLL */

! update the arrays of increments
!
        do j = 1, jgrids
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
    do j = 1, jgrids
      do i = 1, igrids

! copy directional vectors of variables for the one dimensional solver
!
        do k = 1, kgrids
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
        call hll(nvars, kgrids, ul, fl)
#endif /* HLL */

! update the arrays of increments
!
        do k = 1, kgrids
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

    use blocks, only : nvars, idn, imx, imy, imz, ien
    use config, only : igrids, jgrids, kgrids, nghost, gamma

    implicit none

! input arguments
!
    real, dimension(nvars,igrids,jgrids,kgrids), intent(in)  :: u

! local variables
!
    integer :: i, j, k
    real    :: vv, v, c
    real    :: maxspeed

! local arrays
!
    real, dimension(nvars,igrids) :: q
!
!----------------------------------------------------------------------
!
    maxspeed = 0.0

! iterate over all points and calculate maximum speed
!
#if NDIMS == 3
    do k = nghost, kgrids - nghost
#else /* NDIMS == 3 */
    k = 1
#endif /* NDIMS == 3 */
      do j = nghost, jgrids - nghost

        call cons2prim(nvars,igrids,u(:,:,j,k),q(:,:))

        do i = nghost, igrids - nghost

! calculate the velocity
!
          vv = sum(q(imx:imz,i)**2)
          v  = sqrt(vv)

! calculate the maximum characteristic speed
!
#ifdef ADI
          c = sqrt(gamma*q(ien,i)/q(idn,i))
#endif /* ADI */
#ifdef ISO
          c = csnd
#endif /* ISO */

! calculate maximum of the speed
!
          maxspeed = max(maxspeed, v + c)
        enddo
      enddo
#if NDIMS == 3
    enddo
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end function maxspeed

!===============================================================================
!
end module
