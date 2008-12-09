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
!         call hll(igrids, ul, fl)
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
!         call hll(jgrids, ul, fl)
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
!         call hll(kgrids, ul, fl)
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

  end subroutine update

!======================================================================
!
end module
