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

! default problem parameter values
!
  real(kind=8), save :: dens = 1.00d+00
  real(kind=8), save :: pres = 1.00d+00
  real(kind=8), save :: bamp = 1.00d+00
  real(kind=8), save :: bper = 0.00d+00
  real(kind=8), save :: bgui = 0.00d+00
  real(kind=8), save :: vper = 1.00d-02
  real(kind=8), save :: xcut = 4.00d-01
  real(kind=8), save :: ycut = 1.00d-02
  real(kind=8), save :: yth  = 1.00d-16
  real(kind=8), save :: pth  = 1.00d-02
  real(kind=8), save :: pmag = 5.00d-01

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
    use parameters, only : get_parameter_string, get_parameter_real

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

! get the reconnection problem parameters
!
    call get_parameter_real("dens", dens)
    call get_parameter_real("pres", pres)
    call get_parameter_real("bamp", bamp)
    call get_parameter_real("bper", bper)
    call get_parameter_real("bgui", bgui)
    call get_parameter_real("vper", vper)
    call get_parameter_real("xcut", xcut)
    call get_parameter_real("ycut", ycut)
    call get_parameter_real("yth" , yth )
    call get_parameter_real("pth" , pth )

! calculate the maximum magnetic pressure
!
    pmag = 0.5d+00 * (bamp**2 + bgui**2)

! print information about the user problem such as problem name, its
! parameters, etc.
!
    if (verbose) then

      write (*,*)
      write (*,"(1x,a)") "User problem:"
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
    use constants  , only : pi2
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, adx, ady, adz
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use equations  , only : csnd2
    use operators  , only : curl
    use random     , only : randomn

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: yt, yp

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
    real(kind=8), dimension(jm)    :: pm
    real(kind=8), dimension(3)     :: dh
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imp)
#endif /* PROFILE */

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)

! prepare cell sizes
!
    dh(1) = adx(pdata%meta%level)
    dh(2) = ady(pdata%meta%level)
    dh(3) = adz(pdata%meta%level)

! calculate the perturbation of magnetic field
!
    if (bper /= 0.0d+00) then

! initiate the vector potential (we use velocity components to store vector
! potential temporarily, and we store magnetic field perturbation in U)
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im
            yt = abs(y(j) / yth)
            yt = log(exp(yt) + exp(- yt))
            yp = y(j) / pth

            pdata%q(ivx,i,j,k) = 0.0d+00
            pdata%q(ivy,i,j,k) = 0.0d+00
            pdata%q(ivz,i,j,k) = bper * cos(pi2 * x(i)) * exp(- yp * yp) / pi2
          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

! calculate magnetic field components from vector potential
!
      call curl(dh(1:3), pdata%q(ivx:ivz,1:im,1:jm,1:km)                       &
                                            , pdata%q(ibx:ibz,1:im,1:jm,1:km))

    else

! reset magnetic field components
!
      pdata%q(ibx:ibz,:,:,:) = 0.0d+00

    end if ! bper /= 0.0

! iterate over all positions in the XZ plane
!
    do k = 1, km
      do i = 1, im

! if magnetic field is present, set its initial configuration
!
        if (ibx > 0) then

! set antiparallel magnetic field component
!
          do j = 1, jm
            q(ibx,j) = bamp * tanh(y(j) / yth)
          end do ! j = 1, jm

! set tangential magnetic field components
!
          q(iby,1:jm) = 0.0d+00
          q(ibz,1:jm) = bgui
          q(ibp,1:jm) = 0.0d+00

! calculate local magnetic pressure
!
          pm(1:jm)    = 0.5d+00 * sum(q(ibx:ibz,:) * q(ibx:ibz,:), 1)

! add magnetic field perturbation
!
          if (bper /= 0.0d+00) then
            q(ibx,1:jm) = q(ibx,1:jm) + pdata%q(ibx,i,1:jm,k)
            q(iby,1:jm) = q(iby,1:jm) + pdata%q(iby,i,1:jm,k)
            q(ibz,1:jm) = q(ibz,1:jm) + pdata%q(ibz,i,1:jm,k)
          end if ! bper /= 0.0

        end if ! ibx > 0

! set the uniform density and pressure
!
        if (ipr > 0) then
          q(idn,1:jm) = dens
          q(ipr,1:jm) = pres + (pmag - pm(1:jm))
        else
          q(idn,1:jm) = dens + (pmag - pm(1:jm)) / csnd2
        end if

! reset velocity components
!
        q(ivx,1:jm) = 0.0d+00
        q(ivy,1:jm) = 0.0d+00
        q(ivz,1:jm) = 0.0d+00

! set the random velocity field in a layer near current sheet
!
        if (abs(x(i)) <= xcut) then
          do j = 1, jm
            if (abs(y(j)) <= ycut) then

              q(ivx,j) = vper * randomn()
              q(ivy,j) = vper * randomn()
#if NDIMS == 3
              q(ivz,j) = vper * randomn()
#endif /* NDIMS == 3 */

            end if ! |y| < ycut
          end do ! j = 1, jm
        end if ! |x| < xcut

! convert the primitive variables to conservative ones
!
        call prim2cons(jm, q(1:nv,1:jm), u(1:nv,1:jm))

! copy the primitive variables to the current block
!
        pdata%q(1:nv,i,1:jm,k) = q(1:nv,1:jm)

! copy the conserved variables to the current block
!
        pdata%u(1:nv,i,1:jm,k) = u(1:nv,1:jm)

      end do ! i = 1, im
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
