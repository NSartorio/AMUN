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
    real, parameter :: cf  = 1.0d0 / NDIMS
!
!-------------------------------------------------------------------------------
!
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
          fr   = pdata%q(idn,ip2,j,k) - pdata%q(idn,i,j,k)
          fl   = pdata%q(idn,im2,j,k) - pdata%q(idn,i,j,k)
          fc   = pdata%q(idn,ip2,j,k) + pdata%q(idn,im2,j,k)                   &
                                                  + 2.0d0 * pdata%q(idn,i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = pdata%q(idn,i,jp2,k) - pdata%q(idn,i,j,k)
          fl   = pdata%q(idn,i,jm2,k) - pdata%q(idn,i,j,k)
          fc   = pdata%q(idn,i,jp2,k) + pdata%q(idn,i,jm2,k)                   &
                                                  + 2.0d0 * pdata%q(idn,i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pdata%q(idn,i,j,kp2) - pdata%q(idn,i,j,k)
          fl   = pdata%q(idn,i,j,km2) - pdata%q(idn,i,j,k)
          fc   = pdata%q(idn,i,j,kp2) + pdata%q(idn,i,j,km2)                   &
                                                  + 2.0d0 * pdata%q(idn,i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#ifdef ADI
! pressure
!
          fr   = pdata%q(ipr,ip2,j,k) - pdata%q(ipr,i,j,k)
          fl   = pdata%q(ipr,im2,j,k) - pdata%q(ipr,i,j,k)
          fc   = pdata%q(ipr,ip2,j,k) + pdata%q(ipr,im2,j,k)                   &
                                                  + 2.0d0 * pdata%q(ipr,i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = pdata%q(ipr,i,jp2,k) - pdata%q(ipr,i,j,k)
          fl   = pdata%q(ipr,i,jm2,k) - pdata%q(ipr,i,j,k)
          fc   = pdata%q(ipr,i,jp2,k) + pdata%q(ipr,i,jm2,k)                   &
                                                  + 2.0d0 * pdata%q(ipr,i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pdata%q(ipr,i,j,kp2) - pdata%q(ipr,i,j,k)
          fl   = pdata%q(ipr,i,j,km2) - pdata%q(ipr,i,j,k)
          fc   = pdata%q(ipr,i,j,kp2) + pdata%q(ipr,i,j,km2)                   &
                                                  + 2.0d0 * pdata%q(ipr,i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#endif /* ADI */
#ifdef MHD
! X magnetic component
!
          fr   = pdata%q(ibx,ip2,j,k) - pdata%q(ibx,i,j,k)
          fl   = pdata%q(ibx,im2,j,k) - pdata%q(ibx,i,j,k)
          fc   = abs(pdata%q(ibx,ip2,j,k)) + abs(pdata%q(ibx,im2,j,k))         &
                                             + 2.0d0 * abs(pdata%q(ibx,i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = pdata%q(ibx,i,jp2,k) - pdata%q(ibx,i,j,k)
          fl   = pdata%q(ibx,i,jm2,k) - pdata%q(ibx,i,j,k)
          fc   = abs(pdata%q(ibx,i,jp2,k)) + abs(pdata%q(ibx,i,jm2,k))         &
                                             + 2.0d0 * abs(pdata%q(ibx,i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pdata%q(ibx,i,j,kp2) - pdata%q(ibx,i,j,k)
          fl   = pdata%q(ibx,i,j,km2) - pdata%q(ibx,i,j,k)
          fc   = abs(pdata%q(ibx,i,j,kp2)) + abs(pdata%q(ibx,i,j,km2))         &
                                             + 2.0d0 * abs(pdata%q(ibx,i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Y magnetic component
!
          fr   = pdata%q(iby,ip2,j,k) - pdata%q(iby,i,j,k)
          fl   = pdata%q(iby,im2,j,k) - pdata%q(iby,i,j,k)
          fc   = abs(pdata%q(iby,ip2,j,k)) + abs(pdata%q(iby,im2,j,k))         &
                                             + 2.0d0 * abs(pdata%q(iby,i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = pdata%q(iby,i,jp2,k) - pdata%q(iby,i,j,k)
          fl   = pdata%q(iby,i,jm2,k) - pdata%q(iby,i,j,k)
          fc   = abs(pdata%q(iby,i,jp2,k)) + abs(pdata%q(iby,i,jm2,k))         &
                                             + 2.0d0 * abs(pdata%q(iby,i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pdata%q(iby,i,j,kp2) - pdata%q(iby,i,j,k)
          fl   = pdata%q(iby,i,j,km2) - pdata%q(iby,i,j,k)
          fc   = abs(pdata%q(iby,i,j,kp2)) + abs(pdata%q(iby,i,j,km2))         &
                                             + 2.0d0 * abs(pdata%q(iby,i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Z magnetic component
!
          fr   = pdata%q(ibz,ip2,j,k) - pdata%q(ibz,i,j,k)
          fl   = pdata%q(ibz,im2,j,k) - pdata%q(ibz,i,j,k)
          fc   = abs(pdata%q(ibz,ip2,j,k)) + abs(pdata%q(ibz,im2,j,k))         &
                                             + 2.0d0 * abs(pdata%q(ibz,i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = pdata%q(ibz,i,jp2,k) - pdata%q(ibz,i,j,k)
          fl   = pdata%q(ibz,i,jm2,k) - pdata%q(ibz,i,j,k)
          fc   = abs(pdata%q(ibz,i,jp2,k)) + abs(pdata%q(ibz,i,jm2,k))         &
                                             + 2.0d0 * abs(pdata%q(ibz,i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pdata%q(ibz,i,j,kp2) - pdata%q(ibz,i,j,k)
          fl   = pdata%q(ibz,i,j,km2) - pdata%q(ibz,i,j,k)
          fc   = abs(pdata%q(ibz,i,j,kp2)) + abs(pdata%q(ibz,i,j,km2))         &
                                             + 2.0d0 * abs(pdata%q(ibz,i,j,k))
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

    if (cref >= crefmax) then
      criterion =  1
    end if
    if (cref <= crefmin) then
      criterion = -1
    end if

    return

!-------------------------------------------------------------------------------
!
  end function check_refinement_criterion

!===============================================================================
!
end module refinement