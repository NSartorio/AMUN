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
!!*****************************************************************************
!!
!! module: OPERATORS
!!
!!  This module provides differential operators like gradient, divergence, or
!!  curl.
!!
!!*****************************************************************************
!
module operators

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
  integer     , save :: imi, imd, img, imc, iml, im1, im2
#endif /* PROFILE */

! by default everything is public
!
  private

! declare public subroutines
!
  public :: initialize_operators, finalize_operators
  public :: divergence, gradient, curl, laplace
  public :: derivative_1st, derivative_2nd

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
! subroutine INITIALIZE_OPERATORS:
! -------------------------------
!
!   Subroutine initializes the module structures, pointers and variables.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_operators(verbose, iret)

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
    call set_timer('operators:: initialization', imi)
    call set_timer('operators:: divergence'    , imd)
    call set_timer('operators:: gradient'      , img)
    call set_timer('operators:: curl'          , imc)
    call set_timer('operators:: laplace'       , iml)
    call set_timer('operators:: 1st derivative', im1)
    call set_timer('operators:: 2nd derivative', im2)

! start accounting time for the module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_operators
!
!===============================================================================
!
! subroutine FINALIZE_OPERATORS:
! -----------------------------
!
!   Subroutine releases the memory used by module variables and pointers.
!
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_operators(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout)    :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_operators
!
!===============================================================================
!
! subroutine DIVERGENCE:
! ---------------------
!
!   Subroutine calculates the cell centered divergence of the input vector
!   field.
!
!      div(U) = ∇.([Ux, Uy, Uz]) = ∂x Ux + ∂y Uy + ∂z Uz
!
!   Arguments:
!
!     dh  - the spacial intervals in all direction;
!     u   - the input vector field;
!     v   - the output divergence field;
!
!===============================================================================
!
  subroutine divergence(dh, u, v)

! local variables are not implicit by default
!
    implicit none

! input and output variables
!
    real(kind=8), dimension(3)      , intent(in)  :: dh
    real(kind=8), dimension(:,:,:,:), intent(in)  :: u
    real(kind=8), dimension(:,:,:)  , intent(out) :: v

! local arrays
!
    real(kind=8), dimension(:,:,:), allocatable :: w

! local variables
!
    integer :: dir
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the divergence operator calculation
!
    call start_timer(imd)
#endif /* PROFILE */

! allocate temporary array
!
    allocate(w(size(u,2), size(u,3), size(u,4)))

! reset the output array
!
    v(:,:,:) = 0.0d+00

! iterate over directions and update divergence with directional derivatives
!
    do dir = 1, NDIMS

! calculate derivative along the current direction
!
      call derivative_1st(dir, dh(dir), u(dir,:,:,:), w(:,:,:))

! update the divergence array
!
      v(:,:,:) = v(:,:,:) + w(:,:,:)

    end do ! directions

! deallocate temporary array
!
    deallocate(w)

#ifdef PROFILE
! stop accounting time for the divergence operator calculation
!
    call stop_timer(imd)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine divergence
!
!===============================================================================
!
! subroutine GRADIENT:
! -------------------
!
!   Subroutine calculates the cell centered gradient of the input scalar field.
!
!      grad(U) = ∇ U = [ ∂x U, ∂y U, ∂z U ]
!
!   Arguments:
!
!     dh  - the spacial intervals in all direction;
!     u   - the input scalar field;
!     v   - the output gradient vector field;
!
!===============================================================================
!
  subroutine gradient(dh, u, v)

! local variables are not implicit by default
!
    implicit none

! input and output variables
!
    real(kind=8), dimension(3)      , intent(in)  :: dh
    real(kind=8), dimension(:,:,:)  , intent(in)  :: u
    real(kind=8), dimension(:,:,:,:), intent(out) :: v

! local variables
!
    integer :: dir
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the gradient operator calculation
!
    call start_timer(img)
#endif /* PROFILE */

! reset the output array
!
    v(:,:,:,:) = 0.0d+00

! iterate over directions and calculate gradient components
!
    do dir = 1, NDIMS

! calculate derivative along the current direction
!
      call derivative_1st(dir, dh(dir), u(:,:,:), v(dir,:,:,:))

    end do ! directions

#ifdef PROFILE
! stop accounting time for the gradient operator calculation
!
    call stop_timer(img)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine gradient
!
!===============================================================================
!
! subroutine CURL:
! ---------------
!
!   Subroutine calculates the cell centered curl of the input vector field.
!
!      curl(U) = ∇x([Ux, Uy, Uz])
!                              = [∂y Uz - ∂z Uy, ∂z Ux - ∂x Uz, ∂x Uy - ∂y Ux]
!
!   Arguments:
!
!     dh  - the spacial intervals in all direction;
!     u   - the input vector field;
!     v   - the output divergence field;
!
!===============================================================================
!
  subroutine curl(dh, u, v)

! local variables are not implicit by default
!
    implicit none

! input and output variables
!
    real(kind=8), dimension(3)      , intent(in)  :: dh
    real(kind=8), dimension(:,:,:,:), intent(in)  :: u
    real(kind=8), dimension(:,:,:,:), intent(out) :: v

! local arrays
!
    real(kind=8), dimension(:,:,:), allocatable :: w
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the rotation operator calculation
!
    call start_timer(imc)
#endif /* PROFILE */

! allocate temporary array
!
    allocate(w(size(u,2), size(u,3), size(u,4)))

! === calculate Vx component ===
!
! contribution from the Y derivative of Uz
!
    call derivative_1st(2, dh(2), u(3,:,:,:), w(:,:,:))

! update Vx
!
    v(1,:,:,:) = w(:,:,:)

#if NDIMS == 3
! contribution from the Z derivative of Uy
!
    call derivative_1st(3, dh(3), u(2,:,:,:), w(:,:,:))

! update Vx
!
    v(1,:,:,:) = v(1,:,:,:) - w(:,:,:)
#endif /* NDIMS == 3 */

! === calculate Vy component ===
!
#if NDIMS == 3
! contribution from the Z derivative of Ux
!
    call derivative_1st(3, dh(3), u(1,:,:,:), w(:,:,:))

! update Vy
!
    v(2,:,:,:) = w(:,:,:)

! contribution from the X derivative of Uz
!
    call derivative_1st(1, dh(1), u(3,:,:,:), w(:,:,:))

! update Vy
!
    v(2,:,:,:) = v(2,:,:,:) - w(:,:,:)
#else /* NDIMS == 3 */
! contribution from the X derivative of Az
!
    call derivative_1st(1, dh(1), u(3,:,:,:), w(:,:,:))

! update Vy
!
    v(2,:,:,:) = - w(:,:,:)
#endif /* NDIMS == 3 */

! === calculate Vz component ===
!
! contribution from the X derivative of Uy
!
    call derivative_1st(1, dh(1), u(2,:,:,:), w(:,:,:))

! update Vz
!
    v(3,:,:,:) = w(:,:,:)

! contribution from the Y derivative of Ux
!
    call derivative_1st(2, dh(2), u(1,:,:,:), w(:,:,:))

! update Vz
!
    v(3,:,:,:) = v(3,:,:,:) - w(:,:,:)

! deallocate temporary array
!
    deallocate(w)

#ifdef PROFILE
! stop accounting time for the rotation operator calculation
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine curl
!
!===============================================================================
!
! subroutine LAPLACE:
! ------------------
!
!   Subroutine calculates the Laplace operator of the input scalar field.
!
!      laplace(U) = Δ(U) = ∂²x U + ∂²y U + ∂²z U
!
!   Arguments:
!
!     dh  - the spacial intervals in all direction;
!     u   - the input scalar field U;
!     v   - the output field representing the laplacian of u;
!
!===============================================================================
!
  subroutine laplace(dh, u, v)

! local variables are not implicit by default
!
    implicit none

! input and output variables
!
    real(kind=8), dimension(3)    , intent(in)  :: dh
    real(kind=8), dimension(:,:,:), intent(in)  :: u
    real(kind=8), dimension(:,:,:), intent(out) :: v

! local arrays
!
    real(kind=8), dimension(:,:,:), allocatable :: w
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the laplace operator calculation
!
    call start_timer(iml)
#endif /* PROFILE */

! allocate temporary array
!
    allocate(w(size(u,1), size(u,2), size(u,3)))

! calculate the second derivative of U along the X direction
!
    call derivative_2nd(1, dh(1), u(:,:,:), w(:,:,:))

! update V
!
    v(:,:,:) = w(:,:,:)

! calculate the second derivative of U along the X direction
!
    call derivative_2nd(2, dh(2), u(:,:,:), w(:,:,:))

! update V
!
    v(:,:,:) = v(:,:,:) + w(:,:,:)

#if NDIMS == 3
! calculate the second derivative of U along the Z direction
!
    call derivative_2nd(3, dh(3), u(:,:,:), w(:,:,:))

! update V
!
    v(:,:,:) = v(:,:,:) + w(:,:,:)
#endif /* NDIMS == 3 */

! deallocate temporary array
!
    deallocate(w)

#ifdef PROFILE
! stop accounting time for the laplace operator calculation
!
    call stop_timer(iml)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine laplace
!
!===============================================================================
!
! subroutine DERIVATIVE_1ST:
! -------------------------
!
!   Subroutine calculates the first order derivative of the input scalar field
!   along a given direction.
!
!   Arguments:
!
!     dir - the direction of derivative;
!     dh  - the spacial interval;
!     u   - the input scalar field;
!     v   - the first derivative of the input field;
!
!===============================================================================
!
  subroutine derivative_1st(dir, dh, u, v)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                       , intent(in)  :: dir
    real(kind=8)                  , intent(in)  :: dh
    real(kind=8), dimension(:,:,:), intent(in)  :: u
    real(kind=8), dimension(:,:,:), intent(out) :: v

! local variables
!
    integer :: m0, m1, m2
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! return if the direction is wrong
!
    if (dir < 1 .or. dir > NDIMS .or. dh == 0.0d+00) return
#endif /* DEBUG */

#ifdef PROFILE
! start accounting time for the 1st derivative calculation
!
    call start_timer(im1)
#endif /* PROFILE */

! prepare index limits
!
    m0 = size(u, dir)
    m1 = m0 - 1
    m2 = m0 - 2

! select the direction
!
    select case(dir)

! derivative along the X direction
!
    case(1)

      v(2:m1,:,:) = 0.5d+00 * (u(3:m0,:,:) - u(1:m2,:,:)) / dh
      v(1   ,:,:) =           (u(2   ,:,:) - u(1   ,:,:)) / dh
      v(  m0,:,:) =           (u(  m0,:,:) - u(  m1,:,:)) / dh

! derivative along the Y direction
!
    case(2)

      v(:,2:m1,:) = 0.5d+00 * (u(:,3:m0,:) - u(:,1:m2,:)) / dh
      v(:,1   ,:) =           (u(:,2   ,:) - u(:,1   ,:)) / dh
      v(:,  m0,:) =           (u(:,  m0,:) - u(:,  m1,:)) / dh

#if NDIMS == 3
! derivative along the Z direction
!
    case(3)

      v(:,:,2:m1) = 0.5d+00 * (u(:,:,3:m0) - u(:,:,1:m2)) / dh
      v(:,:,1   ) =           (u(:,:,2   ) - u(:,:,1   )) / dh
      v(:,:,  m0) =           (u(:,:,  m0) - u(:,:,  m1)) / dh
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for the 1st derivative calculation
!
    call stop_timer(im1)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine derivative_1st
!
!===============================================================================
!
! subroutine DERIVATIVE_2ND:
! -------------------------
!
!   Subroutine calculates the second order derivative of the input scalar field
!   along a given direction.
!
!   Arguments:
!
!     dir - the direction of derivative;
!     dh  - the spacial interval;
!     u   - the input scalar field;
!     v   - the output scalar field representing the second derivative of u;
!
!===============================================================================
!
  subroutine derivative_2nd(dir, dh, u, v)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                       , intent(in)  :: dir
    real(kind=8)                  , intent(in)  :: dh
    real(kind=8), dimension(:,:,:), intent(in)  :: u
    real(kind=8), dimension(:,:,:), intent(out) :: v

! local variables
!
    integer :: m0, m1, m2
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! return if the direction is wrong
!
    if (dir < 1 .or. dir > NDIMS .or. dh == 0.0d+00) return
#endif /* DEBUG */

#ifdef PROFILE
! start accounting time for the 2nd derivative calculation
!
    call start_timer(im2)
#endif /* PROFILE */

! prepare index limits
!
    m0 = size(u, dir)
    m1 = m0 - 1
    m2 = m0 - 2

! select the direction
!
    select case(dir)

! derivative along the X direction
!
    case(1)

      v(2:m1,:,:) = ((u(3:m0,:,:) + u(1:m2,:,:)) - 2.0d+00 * u(2:m1,:,:)) / dh**2
      v(1   ,:,:) = 0.0d+00
      v(  m0,:,:) = 0.0d+00

! derivative along the Y direction
!
    case(2)

      v(:,2:m1,:) = ((u(:,3:m0,:) + u(:,1:m2,:)) - 2.0d+00 * u(:,2:m1,:)) / dh**2
      v(:,1   ,:) = 0.0d+00
      v(:,  m0,:) = 0.0d+00

#if NDIMS == 3
! derivative along the Z direction
!
    case(3)

      v(:,:,2:m1) = ((u(:,:,3:m0) + u(:,:,1:m2)) - 2.0d+00 * u(:,:,2:m1)) / dh**2
      v(:,:,1   ) = 0.0d+00
      v(:,:,  m0) = 0.0d+00
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for the 2nd derivative calculation
!
    call stop_timer(im2)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine derivative_2nd

!===============================================================================
!
end module operators
