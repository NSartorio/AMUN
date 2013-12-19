!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!  This module handles the error estimation and refinement criterion
!!  determination.
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
  real(kind=8), save :: crefmin = 4.0d-02
  real(kind=8), save :: crefmax = 2.0d-01
  real(kind=8), save :: epsref  = 1.0d-02

! flags for variable included in the refinement criterion calculation
!
  logical     , save :: dens_ref  = .true.
  logical     , save :: pres_ref  = .true.
  logical     , save :: velo_ref  = .false.
  logical     , save :: magn_ref  = .false.

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
!   Subroutine initializes module REFINEMENT.
!
!
!===============================================================================
!
  subroutine initialize_refinement()

! include external procedures and variables
!
    use parameters, only : get_parameter_real, get_parameter_string

! local variables are not implicit by default
!
    implicit none

! local variables
!
    character(len=255) :: variables = "dens pres"
!
!-------------------------------------------------------------------------------
!
! get the refinement parameters
!
    call get_parameter_real("crefmin", crefmin)
    call get_parameter_real("crefmax", crefmax)
    call get_parameter_real("epsref" , epsref )

! get variables to include in the refinement criterion calculation
!
    call get_parameter_string("refinement_variables", variables)

! check if density should be take into account
!
    dens_ref = index(variables, 'dens') > 0
    pres_ref = index(variables, 'pres') > 0
    velo_ref = index(variables, 'velo') > 0
    magn_ref = index(variables, 'magn') > 0

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
!   the block needs to be refined, 0 if there is no need for the refinement,
!   and -1 if the block can be derefined.
!
!   Arguments:
!
!     pdata - pointer to the data block for which the refinement criterion
!             has to be determined;
!
!===============================================================================
!
  function check_refinement_criterion(pdata) result(criterion)

! variables and subroutines imported from other modules
!
    use blocks     , only : block_data
    use equations  , only : idn, ipr, ivx, ivy, ivz, ibx, iby, ibz

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(in) :: pdata

! return variable
!
    integer(kind=1)                       :: criterion

! local variables
!
    real(kind=4) :: cref
!
!-------------------------------------------------------------------------------
!
! reset indicators
!
    cref = 0.0e+00

! find the second derivative error for density
!
    if (dens_ref) cref = max(cref, second_derivative_error(idn, pdata))

! find the second derivative error for pressure
!
    if (pres_ref) cref = max(cref, second_derivative_error(ipr, pdata))

! find the second derivative error for velocity
!
    if (velo_ref) then

      cref = max(cref, second_derivative_error(ivx, pdata))
      cref = max(cref, second_derivative_error(ivy, pdata))
#if NDIMS == 3
      cref = max(cref, second_derivative_error(ivz, pdata))
#endif /* NDIMS == 3 */

    end if

! find the second derivative error for magnetic field
!
    if (magn_ref) then

      cref = max(cref, second_derivative_error(ibx, pdata))
      cref = max(cref, second_derivative_error(iby, pdata))
#if NDIMS == 3
      cref = max(cref, second_derivative_error(ibz, pdata))
#endif /* NDIMS == 3 */

    end if

! return the refinement flag depending on the condition value
!
    criterion = 0

    if (cref > crefmax) then
      criterion =  1
    end if
    if (cref < crefmin) then
      criterion = -1
    end if

! return the refinement flag
!
    return

!-------------------------------------------------------------------------------
!
  end function check_refinement_criterion
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
! function SECOND_DERIVATIVE_ERROR:
! --------------------------------
!
!   Function calculate the second derivative error for a given data block
!   and selected primitive variables.  The total error is returned then.
!
!   Arguments:
!
!     iqt   - the index of primitive variable;
!     pdata - pointer to the data block for which error is calculated;
!
!===============================================================================
!
  function second_derivative_error(iqt, pdata) result(error)

! variables and subroutines imported from other modules
!
    use blocks     , only : block_data
    use coordinates, only : ib, jb, kb, ie, je, ke

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                  , intent(in) :: iqt
    type(block_data), pointer, intent(in) :: pdata

! return variable
!
    real(kind=4)                          :: error

! local variables
!
    integer      :: i, im1, ip1
    integer      :: j, jm1, jp1
    integer      :: k, km1, kp1
    real(kind=4) :: fl, fr, fc, fx, fy, fz
    real(kind=4) :: eloc

! local parameters
!
    real(kind=8), parameter :: eps = epsilon(1.0d+00)
!
!-------------------------------------------------------------------------------
!
! reset indicators
!
    eloc = 0.0e+00

! calculate local refinement criterion for variable which exists
!
    if (iqt > 0) then

! find the maximum smoothness indicator over all cells
!
      do k = kb, ke
#if NDIMS == 3
        km1 = k - 1
        kp1 = k + 1
#endif /* NDIMS == 3 */
        do j = jb, je
          jm1 = j - 1
          jp1 = j + 1
          do i = ib, ie
            im1 = i - 1
            ip1 = i + 1

! calculate the second derivative error the X direction
!
            fr   = pdata%q(iqt,ip1,j,k) - pdata%q(iqt,i  ,j,k)
            fl   = pdata%q(iqt,im1,j,k) - pdata%q(iqt,i  ,j,k)
            fc   = abs(pdata%q(iqt,ip1,j,k)) + abs(pdata%q(iqt,im1,j,k))       &
                                           + 2.0d+00 * abs(pdata%q(iqt,i,j,k))
            fx   = (fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)

! calculate the second derivative error along the Y direction
!
            fr   = pdata%q(iqt,i,jp1,k) - pdata%q(iqt,i,j  ,k)
            fl   = pdata%q(iqt,i,jm1,k) - pdata%q(iqt,i,j  ,k)
            fc   = abs(pdata%q(iqt,i,jp1,k)) + abs(pdata%q(iqt,i,jm1,k))       &
                                           + 2.0d+00 * abs(pdata%q(iqt,i,j,k))
            fy   = (fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)

#if NDIMS == 3
! calculate the second derivative error along the Z direction
!
            fr   = pdata%q(iqt,i,j,kp1) - pdata%q(iqt,i,j,k  )
            fl   = pdata%q(iqt,i,j,km1) - pdata%q(iqt,i,j,k  )
            fc   = abs(pdata%q(iqt,i,j,kp1)) + abs(pdata%q(iqt,i,j,km1))       &
                                           + 2.0d+00 * abs(pdata%q(iqt,i,j,k))
            fz   = (fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)
#endif /* NDIMS == 3 */

! calculate the square of the second derivative error
!
#if NDIMS == 2
            eloc = max(eloc, fx * fx + fy * fy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            eloc = max(eloc, fx * fx + fy * fy + fz * fz)
#endif /* NDIMS == 3 */

          end do
        end do
      end do

    end if ! iqt > 0

! calculate the second derivative error
!
    error = sqrt(eloc)

! return the refinement flag
!
    return

!-------------------------------------------------------------------------------
!
  end function second_derivative_error

!===============================================================================
!
end module refinement
