!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2016 Grzegorz Kowal <grzegorz@amuncode.org>
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
  real(kind=8), save :: jabsmin = 1.0d-03
  real(kind=8), save :: jabsmax = 1.0d-01
  real(kind=8), save :: vortmin = 1.0d-03
  real(kind=8), save :: vortmax = 1.0d-01
  real(kind=8), save :: epsref  = 1.0d-02

! flags for variable included in the refinement criterion calculation
!
  logical, dimension(:), allocatable, save :: qvar_ref
  logical                           , save :: jabs_ref = .false.
  logical                           , save :: vort_ref = .false.

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
    use equations      , only : nv, pvars, ibx
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
    call get_parameter_real("jabsmin", jabsmin)
    call get_parameter_real("jabsmax", jabsmax)
    call get_parameter_real("vortmin", vortmin)
    call get_parameter_real("vortmax", vortmax)
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

! turn on refinement based on vorticity if specified
!
    vort_ref = index(variables, 'vort') > 0
    if (vort_ref) rvars = adjustl(trim(rvars) // ' vort')

! turn on refinement based on current density if specified
!
    if (ibx > 0) then
      jabs_ref = index(variables, 'jabs') > 0
      if (jabs_ref) rvars = adjustl(trim(rvars) // ' jabs')
    end if

! print information about the refinement criterion
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "refined variables      =", trim(rvars)
      write (*,"(4x,a,1x,2e9.2)") "2nd order error limits =", crefmin, crefmax
      if (vort_ref) &
        write (*,"(4x,a,1x,2e9.2)") "vorticity limits       =", vortmin, vortmax
      if (jabs_ref) &
        write (*,"(4x,a,1x,2e9.2)") "current density limits =", jabsmin, jabsmax

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

! reset the refinement criterion flag
!
    criterion = -1

! check the second derivative error from selected primitive variables
!
    do p = 1, nv
      if (qvar_ref(p)) then

        cref = second_derivative_error(p, pdata)

        if (cref > crefmin) criterion = max(criterion, 0)
        if (cref > crefmax) criterion = max(criterion, 1)

      end if
    end do ! p = 1, nv

! check vorticity criterion
!
    if (vort_ref) then
      cref = vorticity_magnitude(pdata)

      if (cref > vortmin) criterion = max(criterion, 0)
      if (cref > vortmax) criterion = max(criterion, 1)
    end if

! check current density criterion
!
    if (jabs_ref) then
      cref = current_density_magnitude(pdata)

      if (cref > jabsmin) criterion = max(criterion, 0)
      if (cref > jabsmax) criterion = max(criterion, 1)
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
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu

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
      do k = kbl, keu
#if NDIMS == 3
        km1 = k - 1
        kp1 = k + 1
#endif /* NDIMS == 3 */
        do j = jbl, jeu
          jm1 = j - 1
          jp1 = j + 1
          do i = ibl, ieu
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
!
!===============================================================================
!
! function VORTICITY_MAGNITUDE:
! ----------------------------
!
!   Function finds the maximum magnitude of vorticity in the block associated
!   with pdata.
!
!   Arguments:
!
!     pdata - pointer to the data block for which error is calculated;
!
!===============================================================================
!
  function vorticity_magnitude(pdata) result(wmax)

! import external procedures and variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use coordinates    , only : adx, ady, adz
    use equations      , only : inx, iny, inz
    use equations      , only : ivx, ivy, ivz
    use operators      , only : curl

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(in) :: pdata

! return variable
!
    real(kind=4)  :: wmax

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: vort

! local arrays
!
    real(kind=8), dimension(3)          :: dh = 1.0d+00
    real(kind=8), dimension(3,im,jm,km) :: wc
!
!-------------------------------------------------------------------------------
!
! reset indicators
!
    wmax = 0.0e+00

! calculate current density W = ∇xV
!
    call curl(dh(:), pdata%q(ivx:ivz,:,:,:), wc(inx:inz,:,:,:))

! find maximum current density
!
      do k = kbl, keu
        do j = jbl, jeu
          do i = ibl, ieu

! calculate the squared magnitude of vorticity
!
            vort = sum(wc(inx:inz,i,j,k)**2)

! find the maximum of squared vorticity
!
            wmax = max(wmax, real(vort, kind=4))

          end do ! i = ibl, ieu
        end do ! j = jbl, jeu
      end do ! kbl, keu

! return the maximum vorticity
!
      wmax = sqrt(wmax)

!-------------------------------------------------------------------------------
!
  end function vorticity_magnitude
!
!===============================================================================
!
! function CURRENT_DENSITY_MAGNITUDE:
! ----------------------------------
!
!   Function finds the maximum magnitude of current density from magnetic field
!   in the block associated with pdata.
!
!   Arguments:
!
!     pdata - pointer to the data block for which error is calculated;
!
!===============================================================================
!
  function current_density_magnitude(pdata) result(jmax)

! import external procedures and variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use coordinates    , only : adx, ady, adz
    use equations      , only : inx, iny, inz
    use equations      , only : ibx, iby, ibz
    use operators      , only : curl

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(in) :: pdata

! return variable
!
    real(kind=4)  :: jmax

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: jabs

! local arrays
!
    real(kind=8), dimension(3)          :: dh = 1.0d+00
    real(kind=8), dimension(3,im,jm,km) :: jc
!
!-------------------------------------------------------------------------------
!
! reset indicators
!
    jmax = 0.0e+00

! return if there is no magnetic field
!
    if (ibx <= 0) return

! calculate current density J = ∇xB
!
    call curl(dh(:), pdata%q(ibx:ibz,:,:,:), jc(inx:inz,:,:,:))

! find maximum current density
!
      do k = kbl, keu
        do j = jbl, jeu
          do i = ibl, ieu

! calculate the squared magnitude of current density
!
            jabs = sum(jc(inx:inz,i,j,k)**2)

! find the maximum of squared current density
!
            jmax = max(jmax, real(jabs, kind=4))

          end do ! i = ibl, ieu
        end do ! j = jbl, jeu
      end do ! kbl, keu

! return the maximum current density
!
      jmax = sqrt(jmax)

!-------------------------------------------------------------------------------
!
  end function current_density_magnitude

!===============================================================================
!
end module refinement
