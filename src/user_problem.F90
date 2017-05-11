!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2017 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: USER_PROBLEM
!!
!!  This module provides subroutines to setup custom problem.
!!
!!*****************************************************************************
!
module user_problem

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
  integer, save :: imi, imp, ims, imu, img, imb
#endif /* PROFILE */

! default problem parameter values are defined here
!

! flag indicating if the gravitational source term is enabled
!
  logical     , save :: gravity_enabled_user = .false.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_USER_PROBLEM:
! ----------------------------------
!
!   Subroutine initializes user problem. It could read problem parameters which
!   are used in all subroutines defining this specific problem.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_user_problem(verbose, iret)

! include external procedures and variables
!
    use parameters, only : get_parameter_string

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64) :: problem_name   = "none"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('user_problem:: initialize'   , imi)
    call set_timer('user_problem:: problem setup', imp)
    call set_timer('user_problem:: shape'        , ims)
    call set_timer('user_problem:: sources'      , imu)
    call set_timer('user_problem:: gravity'      , img)
    call set_timer('user_problem:: boundaries'   , imb)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the problem name
!
    call get_parameter_string("problem", problem_name)

! print information about the user problem such as problem name, its
! parameters, etc.
!
    if (verbose) then

      write (*,"(4x,a14, 9x,'=',2x,a)") "problem name  ", trim(problem_name)

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_user_problem
!
!===============================================================================
!
! subroutine FINALIZE_USER_PROBLEM:
! --------------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_user_problem(iret)

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
  end subroutine finalize_user_problem
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_USER:
! -----------------------------
!
!   Subroutine sets the initial conditions for the user specific problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_user(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! local variables
!
    integer :: i, j, k

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imp)
#endif /* PROFILE */

! set the variables
!
    q(idn,:) = 1.0d+00
    if (ipr > 0) q(ipr,:) = 1.0d+00
    q(ivx,:) = 0.0d+00
    q(ivy,:) = 0.0d+00
    q(ivz,:) = 0.0d+00
    if (ibx > 0) then
      q(ibx,:) = 0.0d+00
      q(iby,:) = 0.0d+00
      q(ibz,:) = 0.0d+00
      q(ibp,:) = 0.0d+00
    end if

! convert the primitive variables to conservative ones
!
    call prim2cons(im, q(1:nv,1:im), u(1:nv,1:im))

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do j = 1, jm

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the problems setup
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_user
!
!===============================================================================
!
! subroutine UPDATE_SHAPES_USER:
! -----------------------------
!
!   Subroutine defines the regions updated by user.
!
!   Arguments:
!
!     pdata - pointer to the data block structure of the currently initialized
!             block;
!     time  - time at the moment of update;
!     dt    - time step since the last update;
!
!===============================================================================
!
  subroutine update_shapes_user(pdata, time, dt)

! include external procedures and variables
!
    use blocks         , only : block_data

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time, dt
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the shape update
!
    call start_timer(ims)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the shape update
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_shapes_user
!
!===============================================================================
!
! subroutine UPDATE_SOURCES_USER:
! ------------------------------
!
!   Subroutine adds the user defined source terms.
!
!   Arguments:
!
!     pdata - the pointer to a data block;
!     t, dt - the time and time increment;
!     du    - the array of variable increment;
!
!===============================================================================
!
  subroutine update_sources_user(pdata, t, dt, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8)                        , intent(in)    :: t, dt
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: du
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for source terms
!
    call start_timer(imu)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for source terms
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_sources_user
!
!===============================================================================
!
! subroutine GRAVITATIONAL_ACCELERATION_USER:
! ------------------------------------------
!
!   Subroutine returns the user defined gravitational acceleration.
!
!   Arguments:
!
!     t, dt   - time and the time increment;
!     x, y, z - rectangular coordinates;
!     acc     - vector of the gravitational acceleration;
!
!===============================================================================
!
  subroutine gravitational_acceleration_user(t, dt, x, y, z, acc)

! include external procedures and variables
!
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)              , intent(in)  :: t, dt
    real(kind=8)              , intent(in)  :: x, y, z
    real(kind=8), dimension(3), intent(out) :: acc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the gravitational acceleration calculation
!
    call start_timer(img)
#endif /* PROFILE */

! reset gravitational acceleration
!
    acc(:) = 0.0d+00

#ifdef PROFILE
! stop accounting time for the gravitational acceleration calculation
!
    call stop_timer(img)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine gravitational_acceleration_user
!
!===============================================================================
!
! subroutine BOUNDARY_USER_X:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the X direction.
!
!   Arguments:
!
!     ic      - the block side along the X direction for the ghost zone update;
!     jl, ju  - the cell index limits for the Y direction;
!     kl, ku  - the cell index limits for the Z direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_x(ic, jl, ju, kl, ku, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: ic
    integer                                     , intent(in)    :: jl, ju
    integer                                     , intent(in)    :: kl, ku
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_x
!
!===============================================================================
!
! subroutine BOUNDARY_USER_Y:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the Y direction.
!
!   Arguments:
!
!     jc      - the block side along the Y direction for the ghost zone update;
!     il, iu  - the cell index limits for the X direction;
!     kl, ku  - the cell index limits for the Z direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_y(jc, il, iu, kl, ku, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: jc
    integer                                     , intent(in)    :: il, iu
    integer                                     , intent(in)    :: kl, ku
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_y
!
!===============================================================================
!
! subroutine BOUNDARY_USER_Z:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the Z direction.
!
!   Arguments:
!
!     kc      - the block side along the Z direction for the ghost zone update;
!     il, iu  - the cell index limits for the X direction;
!     jl, ju  - the cell index limits for the Y direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_z(kc, il, iu, jl, ju, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: kc
    integer                                     , intent(in)    :: il, iu
    integer                                     , intent(in)    :: jl, ju
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_z

!===============================================================================
!
end module user_problem
