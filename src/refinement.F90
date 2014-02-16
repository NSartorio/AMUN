!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2014 Grzegorz Kowal <grzegorz@amuncode.org>
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

#ifdef PROFILE
! import external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer            , save :: iri, irc
#endif /* PROFILE */

! refinement criterion parameters
!
  real(kind=8), save :: crefmin = 2.0d-01
  real(kind=8), save :: crefmax = 8.0d-01
  real(kind=8), save :: epsref  = 1.0d-02

! flags for variable included in the refinement criterion calculation
!
  logical, dimension(:), allocatable, save :: qvar_ref

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_refinement, finalize_refinement
  public :: check_refinement_criterion

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
!===============================================================================
!
! subroutine INITIALIZE_REFINEMENT:
! --------------------------------
!
!   Subroutine initializes module REFINEMENT.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_refinement(verbose, iret)

! import external procedures and variables
!
    use equations      , only : nv, pvars
    use parameters     , only : get_parameter_real, get_parameter_string

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    integer                :: p
    character(len=255)     :: variables = "dens pres"
    character(len=255)     :: rvars     = ""
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('refinement:: initialization', iri)
    call set_timer('refinement:: criterion'     , irc)

! start accounting time for module initialization/finalization
!
    call start_timer(iri)
#endif /* PROFILE */

! get the refinement parameters
!
    call get_parameter_real("crefmin", crefmin)
    call get_parameter_real("crefmax", crefmax)
    call get_parameter_real("epsref" , epsref )

! get variables to include in the refinement criterion calculation
!
    call get_parameter_string("refinement_variables", variables)

! allocate vector for indicators, which variables are taken into account in
! calculating the refinement criterion
!
    allocate(qvar_ref(nv))

! check which primitive variable is used to determine the refinement criterion
!
    do p = 1, nv
      qvar_ref(p) = index(variables, trim(pvars(p))) > 0
      if (qvar_ref(p)) rvars = adjustl(trim(rvars) // ' ' // trim(pvars(p)))
    end do ! p = 1, nv

! print information about the refinement criterion
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "refined variables      =", trim(rvars)
      write (*,"(4x,a,1x,1e9.2)") "derefinement threshold =", crefmin
      write (*,"(4x,a,1x,1e9.2)") "refinement threshold   =", crefmax

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_refinement
!
!===============================================================================
!
! subroutine FINALIZE_REFINEMENT:
! ------------------------------
!
!   Subroutine releases memory used by the module variables.
!
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_refinement(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(iri)
#endif /* PROFILE */

! deallocate refined variable indicators
!
    if (allocated(qvar_ref)) deallocate(qvar_ref)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(iri)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_refinement
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

! import external procedures and variables
!
    use blocks         , only : block_data
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(in) :: pdata

! return variable
!
    integer(kind=4)                       :: criterion

! local variables
!
    integer      :: p
    real(kind=4) :: cref
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the refinement criterion estimation
!
    call start_timer(irc)
#endif /* PROFILE */

! reset indicators
!
    cref = 0.0e+00

! calculate the second derivative error for selected primitive variables only
!
    do p = 1, nv
      if (qvar_ref(p)) cref = max(cref, second_derivative_error(p, pdata))
    end do ! p = 1, nv

! return the refinement flag depending on the condition value
!
    criterion = 0

    if (cref > crefmax) then
      criterion =  1
    end if
    if (cref < crefmin) then
      criterion = -1
    end if

#ifdef PROFILE
! stop accounting time for the refinement criterion estimation
!
    call stop_timer(irc)
#endif /* PROFILE */

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

! import external procedures and variables
!
    use blocks         , only : block_data
    use coordinates    , only : ib, jb, kb, ie, je, ke

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

! local parameters
!
    real(kind=8), parameter :: eps = epsilon(1.0d+00)
!
!-------------------------------------------------------------------------------
!
! reset indicators
!
    error = 0.0e+00

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
            fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)

! calculate the second derivative error along the Y direction
!
            fr   = pdata%q(iqt,i,jp1,k) - pdata%q(iqt,i,j  ,k)
            fl   = pdata%q(iqt,i,jm1,k) - pdata%q(iqt,i,j  ,k)
            fc   = abs(pdata%q(iqt,i,jp1,k)) + abs(pdata%q(iqt,i,jm1,k))       &
                                           + 2.0d+00 * abs(pdata%q(iqt,i,j,k))
            fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)

#if NDIMS == 3
! calculate the second derivative error along the Z direction
!
            fr   = pdata%q(iqt,i,j,kp1) - pdata%q(iqt,i,j,k  )
            fl   = pdata%q(iqt,i,j,km1) - pdata%q(iqt,i,j,k  )
            fc   = abs(pdata%q(iqt,i,j,kp1)) + abs(pdata%q(iqt,i,j,km1))       &
                                           + 2.0d+00 * abs(pdata%q(iqt,i,j,k))
            fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc + eps)
#endif /* NDIMS == 3 */

! take the maximum second derivative error
!
#if NDIMS == 2
            error = max(error, fx, fy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            error = max(error, fx, fy, fz)
#endif /* NDIMS == 3 */

          end do
        end do
      end do

    end if ! iqt > 0

! return the refinement flag
!
    return

!-------------------------------------------------------------------------------
!
  end function second_derivative_error

!===============================================================================
!
end module refinement
