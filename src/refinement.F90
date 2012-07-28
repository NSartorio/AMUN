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
!! module: REFINEMENT
!!
!!  This module handles the refinement criterion calculation.
!!
!!
!!******************************************************************************
!
module refinement

! module variables are not implicit by default
!
  implicit none

! refinement criterion parameters
!
  real             , save :: crefmin = 2.0d-1
  real             , save :: crefmax = 5.0d-1
  real             , save :: epsref  = 1.0d-2

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_refinement, check_refinement_criterion

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
! subroutine INITIALIZE_REFINEMENT:
! --------------------------------
!
!   Subroutine prepares module REFINEMENT.
!
!
!===============================================================================
!
  subroutine initialize_refinement()

! include external procedures and variables
!
    use parameters, only : get_parameter_real

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the refinement parameters
!
    call get_parameter_real  ("crefmin", crefmin)
    call get_parameter_real  ("crefmax", crefmax)
    call get_parameter_real  ("epsref" , epsref )

!-------------------------------------------------------------------------------
!
  end subroutine initialize_refinement
!
!===============================================================================
!
! function CHECK_REFINEMENT_CRITERION:
! -----------------------------------
!
!   Function scans the given data block and checks for the refinement
!   criterion.  It returns +1 if the criterion is met, which indicates that
!   the block needs to be refined, 0 if there is no need for the refinement, and
!   -1 if the block can be derefined.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  function check_refinement_criterion(pdata) result(criterion)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke
    use scheme     , only : cons2prim
    use variables  , only : nvr, nqt
    use variables  , only : idn
#ifdef ADI
    use variables  , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables  , only : ibx, iby, ibz
#endif /* MHD */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! return variable
!
    integer(kind=4) :: criterion

! local variables
!
    integer :: i, j, k, im2, jm2, km2, ip2, jp2, kp2
    real    :: cref, fl, fr, fc, fx, fy, fz

! local arrays
!
    real, dimension(nvr,im)   :: u, q
    real, dimension(im,jm,km) :: dn
#ifdef ADI
    real, dimension(im,jm,km) :: pr
#endif /* ADI */
#ifdef MHD
    real, dimension(im,jm,km) :: bx, by, bz
#endif /* MHD */
    real, parameter           :: cf  = 1.0d0 / NDIMS
!
!-------------------------------------------------------------------------------
!
! convert conserved variables to primitive ones
!
    do k = 1, km
      do j = 1, jm
        dn(1:im,j,k) = pdata%u(idn,1:im,j,k)
#ifdef ADI
        u(1:nqt,1:im) = pdata%u(1:nqt,1:im,j,k)

        call cons2prim(im, u(:,:), q(:,:))

        pr(1:im,j,k) = q(ipr,1:im)
#endif /* ADI */
#ifdef MHD
        bx(1:im,j,k) = pdata%u(ibx,1:im,j,k)
        by(1:im,j,k) = pdata%u(iby,1:im,j,k)
        bz(1:im,j,k) = pdata%u(ibz,1:im,j,k)
#endif /* MHD */
      end do
    end do

! reset indicators
!
    cref = 0.0d0

! find the maximum smoothness indicator over all cells
!
    do k = kb, ke
#if NDIMS == 3
      km2 = k - 2
      kp2 = k + 2
#endif /* NDIMS == 3 */
      do j = jb, je
        jm2 = j - 2
        jp2 = j + 2
        do i = ib, ie
          im2 = i - 2
          ip2 = i + 2

! density
!
          fr   = dn(ip2,j,k) - dn(i,j,k)
          fl   = dn(im2,j,k) - dn(i,j,k)
          fc   = dn(ip2,j,k) + dn(im2,j,k) + 2.0d0 * dn(i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = dn(i,jp2,k) - dn(i,j,k)
          fl   = dn(i,jm2,k) - dn(i,j,k)
          fc   = dn(i,jp2,k) + dn(i,jm2,k) + 2.0d0 * dn(i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = dn(i,j,kp2) - dn(i,j,k)
          fl   = dn(i,j,km2) - dn(i,j,k)
          fc   = dn(i,j,kp2) + dn(i,j,km2) + 2.0d0 * dn(i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#ifdef ADI
! pressure
!
          fr   = pr(ip2,j,k) - pr(i,j,k)
          fl   = pr(im2,j,k) - pr(i,j,k)
          fc   = pr(ip2,j,k) + pr(im2,j,k) + 2.0d0 * pr(i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = pr(i,jp2,k) - pr(i,j,k)
          fl   = pr(i,jm2,k) - pr(i,j,k)
          fc   = pr(i,jp2,k) + pr(i,jm2,k) + 2.0d0 * pr(i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pr(i,j,kp2) - pr(i,j,k)
          fl   = pr(i,j,km2) - pr(i,j,k)
          fc   = pr(i,j,kp2) + pr(i,j,km2) + 2.0d0 * pr(i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#endif /* ADI */
#ifdef MHD
! X magnetic component
!
          fr   = bx(ip2,j,k) - bx(i,j,k)
          fl   = bx(im2,j,k) - bx(i,j,k)
          fc   = abs(bx(ip2,j,k)) + abs(bx(im2,j,k)) + 2.0d0 * abs(bx(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = bx(i,jp2,k) - bx(i,j,k)
          fl   = bx(i,jm2,k) - bx(i,j,k)
          fc   = abs(bx(i,jp2,k)) + abs(bx(i,jm2,k)) + 2.0d0 * abs(bx(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = bx(i,j,kp2) - bx(i,j,k)
          fl   = bx(i,j,km2) - bx(i,j,k)
          fc   = abs(bx(i,j,kp2)) + abs(bx(i,j,km2)) + 2.0d0 * abs(bx(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Y magnetic component
!
          fr   = by(ip2,j,k) - by(i,j,k)
          fl   = by(im2,j,k) - by(i,j,k)
          fc   = abs(by(ip2,j,k)) + abs(by(im2,j,k)) + 2.0d0 * abs(by(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = by(i,jp2,k) - by(i,j,k)
          fl   = by(i,jm2,k) - by(i,j,k)
          fc   = abs(by(i,jp2,k)) + abs(by(i,jm2,k)) + 2.0d0 * abs(by(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = by(i,j,kp2) - by(i,j,k)
          fl   = by(i,j,km2) - by(i,j,k)
          fc   = abs(by(i,j,kp2)) + abs(by(i,j,km2)) + 2.0d0 * abs(by(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Z magnetic component
!
          fr   = bz(ip2,j,k) - bz(i,j,k)
          fl   = bz(im2,j,k) - bz(i,j,k)
          fc   = abs(bz(ip2,j,k)) + abs(bz(im2,j,k)) + 2.0d0 * abs(bz(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = bz(i,jp2,k) - bz(i,j,k)
          fl   = bz(i,jm2,k) - bz(i,j,k)
          fc   = abs(bz(i,jp2,k)) + abs(bz(i,jm2,k)) + 2.0d0 * abs(bz(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = bz(i,j,kp2) - bz(i,j,k)
          fl   = bz(i,j,km2) - bz(i,j,k)
          fc   = abs(bz(i,j,kp2)) + abs(bz(i,j,km2)) + 2.0d0 * abs(bz(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#endif /* MHD */
        end do
      end do
    end do

! return the refinement flag depending on the condition value
!
    criterion = 0

    if (cref .ge. crefmax) then
      criterion =  1
    end if
    if (cref .le. crefmin) then
      criterion = -1
    end if

    return

!-------------------------------------------------------------------------------
!
  end function check_refinement_criterion

!===============================================================================
!
end module refinement
