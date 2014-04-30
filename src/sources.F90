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
!! module: SOURCES
!!
!!  This modules adds source terms.
!!
!!******************************************************************************
!
module sources

#ifdef PROFILE
! include external procedures
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer, save :: imi, imu
#endif /* PROFILE */

! gravitational acceleration coefficient
!
  real(kind=8), save :: gpoint = 0.0d+00

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_sources, finalize_sources
  public :: update_sources

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
! subroutine INITIALIZE_SOURCES:
! -----------------------------
!
!   Subroutine initializes module SOURCES.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_sources(verbose, iret)

! include external procedures and variables
!
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('sources:: initialize', imi)
    call set_timer('sources:: update'    , imu)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get acceleration coefficient
!
      call get_parameter_real("gpoint" , gpoint)

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,1e9.2)") "point mass constant    =", gpoint

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_sources
!
!===============================================================================
!
! subroutine FINALIZE_SOURCES:
! ---------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_sources(iret)

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
    call start_timer(imi)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_sources
!
!===============================================================================
!
! subroutine UPDATE_SOURCES:
! -------------------------
!
!   Subroutine add the source terms.
!
!   Arguments:
!
!     q    - the array of primitive variables;
!     du   - the array of variable increment;
!
!===============================================================================
!
  subroutine update_sources(pdata, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az
    use equations      , only : nv, idn, ivx, ivy, ivz, imx, imy, imz, ien

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: du

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: r2, r3, gc, gx, gy, gz

! local arrays
!
    real(kind=8), dimension(im) :: x
    real(kind=8), dimension(jm) :: y
    real(kind=8), dimension(km) :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for source terms
!
    call start_timer(imu)
#endif /* PROFILE */

! proceed only if the gravitational acceleration coefficient is not zero
!
    if (gpoint /= 0.0d+00) then

! prepare block coordinates
!
      x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
      y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
      z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#endif /* NDIMS == 3 */

! iterate over all positions in the YZ plane
!
      do k = 1, km
        do j = 1, jm
          do i = 1, jm

! calculate distance from the origin
!
#if NDIMS == 2
            r2 = x(i) * x(i) + y(j) * y(j)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            r2 = x(i) * x(i) + y(j) * y(j) + z(k) * z(k)
#endif /* NDIMS == 3 */
            r3 = r2 * sqrt(r2)

! calculate gravitational acceleration factors
!
            gc = gpoint * pdata%q(idn,i,j,k) / max(1.0d-16, r3)
            gx = gc * x(i)
            gy = gc * y(j)
#if NDIMS == 3
            gz = gc * z(k)
#endif /* NDIMS == 3 */

! add source terms to momentum equations
!
            du(imx,i,j,k) = du(imx,i,j,k) + gx
            du(imy,i,j,k) = du(imy,i,j,k) + gy
#if NDIMS == 3
            du(imz,i,j,k) = du(imz,i,j,k) + gz
#endif /* NDIMS == 3 */

! add source terms to total energy equation
!
          if (ien > 0) then

#if NDIMS == 2
              du(ien,i,j,k) = du(ien,i,j,k) + gx * pdata%q(ivx,i,j,k)          &
                                            + gy * pdata%q(ivy,i,j,k)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              du(ien,i,j,k) = du(ien,i,j,k) + gx * pdata%q(ivx,i,j,k)          &
                                            + gy * pdata%q(ivy,i,j,k)          &
                                            + gz * pdata%q(ivz,i,j,k)
#endif /* NDIMS == 3 */

            end if

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

    end if ! gpoint is not zero

#ifdef PROFILE
! stop accounting time for source terms
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_sources

!===============================================================================
!
end module sources
