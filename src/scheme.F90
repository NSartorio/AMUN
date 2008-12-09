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
!======================================================================
!
! update: subroutine sweeps over all directions and integrates
!         the increment of the solution
!
!======================================================================
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
!----------------------------------------------------------------------
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
        call hll(igrids, ul, fl)
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
        call hll(jgrids, ul, fl)
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
        call hll(kgrids, ul, fl)
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

!----------------------------------------------------------------------
!
  end subroutine update
#ifdef HLL
!
!======================================================================
!
! hll: subroutine to compute flux approximated by HLL method
!
!======================================================================
!
  subroutine hll(n, uc, f)

    use blocks, only : nvars
!     use interpolation, only : reconstruct, point2avg

    implicit none

! input/output arguments
!
    integer                 , intent(in)  :: n
    real, dimension(nvars,n), intent(in)  :: uc
    real, dimension(nvars,n), intent(out) :: f

! local variables
!
    integer                  :: p, i
    real, dimension(nvars,n) :: ul, ur, ql, qr, qc
    real, dimension(nvars,n) :: fl, fr, fx
    real, dimension(n)       :: cl, cr
    real                     :: al, ar, ap, div

    integer, parameter, dimension(6)    :: pos = (/ 1, 0, 0, 0, 1, 1 /)
!
!----------------------------------------------------------------------
!
#ifdef CONREC
! reconstruct left and right states of conserved variables
!
    do p = 1, nvars
!       call reconstruct(n,uc(p,:),ul(p,:),ur(p,:),pos(p))
    enddo

! calculate primitive variables
!
!     call cons2prim(n,ul,ql)
!     call cons2prim(n,ur,qr)
#else /* CONREC */

! calculate primitive variables
!
!     call cons2prim(n,uc,qc)

! reconstruct left and right states of primitive variables
!
    do p = 1, nvars
!       call reconstruct(n,qc(p,:),ql(p,:),qr(p,:),pos(p))
    enddo

! calculate conservative variables at states
!
!     call prim2cons(n,ql,ul)
!     call prim2cons(n,qr,ur)
#endif /* CONREC */

! calculate fluxes and speeds
!
!     call fluxspeed(n,ql,ul,fl,cl)
!     call fluxspeed(n,qr,ur,fr,cr)

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

!----------------------------------------------------------------------
!
  end subroutine hll
#endif /* HLL */

!======================================================================
!
end module
