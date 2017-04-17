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
  real(kind=8), save :: kper = 1.00d+00
  real(kind=8), save :: kvec = 1.00d+00
  real(kind=8), save :: xcut = 5.00d-01
  real(kind=8), save :: ycut = 1.00d-01
  real(kind=8), save :: xdec = 1.00d-01
  real(kind=8), save :: ydec = 1.00d-01
  real(kind=8), save :: yth  = 1.00d-16
  real(kind=8), save :: pth  = 1.00d-02
  real(kind=8), save :: pmag = 5.00d-01
  real(kind=8), save :: blim = 1.00d+00
  real(kind=8), save :: zeta = 0.00d+00
  integer     , save :: pert = 0
  integer     , save :: nper = 10

! flag indicating if the gravitational source term is enabled
!
  logical     , save :: gravity_enabled_user = .false.

! allocatable arrays for velocity perturbation
!
  real(kind=8), dimension(:), allocatable :: kx, ky, kz, ux, uy, uz, ph

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
    use constants      , only : pi, pi2
    use coordinates    , only : ng
    use coordinates    , only : ady
    use parameters     , only : get_parameter_string, get_parameter_real       &
                              , get_parameter_integer
    use random         , only : randomu, randomn

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
    character(len=64) :: perturbation   = "noise"
    integer           :: n
    real(kind=8)      :: thh, fc
#if NDIMS == 3
    real(kind=8)      :: thv, tx, ty, tz
#endif /* NDIMS == 3 */
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
    call get_parameter_string("problem"     , problem_name)
    call get_parameter_string("perturbation", perturbation)

! choose the perturbation type
!
    select case(perturbation)
    case('noise', 'random')
      pert = 0
    case('wave')
      pert = 1
    case('vel', 'velocity_mode')
      pert = 2
    case('mag', 'magnetic_mode')
      pert = 3
    case('both')
      pert = 4
    end select

! get the reconnection problem parameters
!
    call get_parameter_real("dens"  , dens)
    call get_parameter_real("pres"  , pres)
    call get_parameter_real("bamp"  , bamp)
    call get_parameter_real("bper"  , bper)
    call get_parameter_real("bgui"  , bgui)
    call get_parameter_real("vper"  , vper)
    call get_parameter_real("kper"  , kper)
    call get_parameter_real("xcut"  , xcut)
    call get_parameter_real("ycut"  , ycut)
    call get_parameter_real("xdec"  , xdec)
    call get_parameter_real("ydec"  , ydec)
    call get_parameter_real("yth"   , yth )
    call get_parameter_real("pth"   , pth )
    call get_parameter_real("blimit", blim)
    call get_parameter_real("zeta"  , zeta)
    call get_parameter_integer("nper", nper)

! calculate the maximum magnetic pressure
!
    pmag = 0.5d+00 * (bamp**2 + bgui**2)

! prepare the wave vector of the perturbation
!
    kvec = pi2 * kper

! upper limit for blim
!
    blim = max(blim, ng * ady(1))

! prepare the random perturbation of velocity
!
    if (pert == 1) then

! allocate phase and wave vector components
!
      allocate(kx(nper), ky(nper), kz(nper))
      allocate(ux(nper), uy(nper), uz(nper))
      allocate(ph(nper))

! choose random wave vector directions
!
      fc = 1.0d+00 / sqrt(1.0d+00 * nper)
      do n = 1, nper
        thh   = pi2 * randomu()
#if NDIMS == 3
        thv   = pi  * randomn()
        ux(n) = cos(thh) * cos(thv)
        uy(n) = sin(thh) * cos(thv)
        uz(n) =            sin(thv)
        kx(n) = pi2 * kper * ux(n)
        ky(n) = pi2 * kper * uy(n)
        kz(n) = pi2 * kper * uz(n)
        tx    = cos(thh) * cos(thv)
        ty    = sin(thh) * cos(thv)
        tz    = sin(thv)
        ux(n) = ty * kz(n) - tz * ky(n)
        uy(n) = tz * kx(n) - tx * kz(n)
        uz(n) = tx * ky(n) - ty * kx(n)
        tx    = sqrt(ux(n)**2 + uy(n)**2 + uz(n)**2)
        ux(n) = fc * ux(n) / tx
        uy(n) = fc * uy(n) / tx
        uz(n) = fc * uz(n) / tx
#else /* NDIMS == 3 */
        kx(n) = pi2 * kper * cos(thh)
        ky(n) = pi2 * kper * sin(thh)
        kz(n) = 0.0d+00
        ux(n) = fc * sin(thh)
        uy(n) = fc * cos(thh)
        uz(n) = 0.0d+00
#endif /* NDIMS == 3 */
        ph(n) = pi2 * randomu()
      end do

    end if

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

! deallocate wave vector components, random directions, and random phase
!
    if (allocated(kx)) deallocate(kx, ky, kz)
    if (allocated(ux)) deallocate(ux, uy, uz)
    if (allocated(ph)) deallocate(ph)

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
    use constants  , only : pi
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az, adx, ady, adz
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use equations  , only : csnd2
    use operators  , only : curl
    use random     , only : randomn, randomu

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! local variables
!
    integer      :: i, j, k, n
    real(kind=8) :: xp, yp, yt, yl, yu
    real(kind=8) :: yrat, itanh
    real(kind=8) :: vx, vy, vz, vv, kv, fv, va

! local arrays
!
    real(kind=8), dimension(3 ,im,jm,km) :: tmp
    real(kind=8), dimension(nv,im)       :: q, u
    real(kind=8), dimension(im)          :: x, fx
    real(kind=8), dimension(jm)          :: y, fy
    real(kind=8), dimension(km)          :: z
    real(kind=8), dimension(jm)          :: pm
    real(kind=8), dimension(3)           :: dh
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
#if NDIMS == 3
    z(1:km) = pdata%meta%ymin + az(pdata%meta%level,1:km)
#endif /* NDIMS == 3 */

! prepare cell sizes
!
    dh(1) = adx(pdata%meta%level)
    dh(2) = ady(pdata%meta%level)
    dh(3) = adz(pdata%meta%level)

! ratio of the current sheet thickness to the cell size
!
    yrat  = yth / dh(2)

! prepare decaying factors
!
    fv = 0.5d+00 * pi
    do i = 1, im
      xp    = fv * min(1.0d+00, max(0.0d+00, abs(x(i)) - xcut) / xdec)
      fx(i) = cos(xp)**2
    end do ! i = 1, im
    do j = 1, jm
      yp    = fv * min(1.0d+00, max(0.0d+00, abs(y(j)) - ycut) / ydec)
      fy(j) = cos(yp)**2
    end do ! i = 1, im

! reset the velocity components
!
    pdata%q(ivx:ivz,1:im,1:jm,1:km) = 0.0d+00

! prepare the random perturbation of velocity
!
    if (pert == 0) then

      if (vper /= 0.0d+00) then

! initiate the random velocity components
!
        do k = 1, km
          do j = 1, jm
            if (fy(j) > 0.0d+00) then
              do i = 1, im
                if (fx(i) > 0.0d+00) then

! calculate the velocity amplitude profile
!
                  va = vper * fx(i) * fy(j)

! get the random direction
!
                  vv = 0.0d+00
                  do while(vv == 0.0d+00)
                    vx = randomn()
                    vy = randomn()
#if NDIMS == 3
                    vz = randomn()
#endif /* NDIMS == 3 */

#if NDIMS == 3
                    vv = sqrt(vx * vx + vy * vy + vz * vz)
#else /* NDIMS == 3 */
                    vv = sqrt(vx * vx + vy * vy)
#endif /* NDIMS == 3 */
                  end do

                  pdata%q(ivx,i,j,k) = va * (vx / vv)
                  pdata%q(ivy,i,j,k) = va * (vy / vv)
#if NDIMS == 3
                  pdata%q(ivz,i,j,k) = va * (vz / vv)
#endif /* NDIMS == 3 */

                end if ! |x| < xcut
              end do ! i = 1, im
            end if ! |y| < ycut
          end do ! j = 1, jm
        end do ! k = 1, km

      end if ! vper /= 0.0

    end if ! pert == 0

! prepare the random perturbation of velocity
!
    if (pert == 1) then

      if (vper /= 0.0d+00) then

! iterate over the block position and initiate the velocity perturbation
!
        do k = 1, km
          do j = 1, jm
            if (fy(j) > 0.0d+00) then
              do i = 1, im
                if (fx(i) > 0.0d+00) then

! calculate the velocity amplitude profile
!
                  fv = vper * fx(i) * fy(j)

! add perturbation components
!
                  do n = 1, nper
                    kv = kx(n) * x(i) + ky(n) * y(j) + kz(n) * z(k) + ph(n)
                    va = fv * sin(kv)

                    pdata%q(ivx,i,j,k) = pdata%q(ivx,i,j,k) + va * ux(n)
                    pdata%q(ivy,i,j,k) = pdata%q(ivy,i,j,k) + va * uy(n)
#if NDIMS == 3
                    pdata%q(ivz,i,j,k) = pdata%q(ivz,i,j,k) + va * uz(n)
#endif /* NDIMS == 3 */
                  end do

                end if ! fx > 0.0
              end do ! i = 1, im
            end if ! fy > 0.0
          end do ! j = 1, jm
        end do ! k = 1, km

      end if ! vper /= 0.0

    end if ! pert == 1

! prepare the wave-like perturbation of velocity
!
    if (pert == 2 .or. pert == 4) then

      if (vper /= 0.0d+00) then

! initiate the potential to calculate the divergence-free velocity components
!
        do k = 1, km
          do j = 1, jm
            do i = 1, im
              xp = kvec * x(i)
              yp = y(j) / pth
              yt = - 0.5d+00 * yp**2

              tmp(1,i,j,k) = 0.0d+00
              tmp(2,i,j,k) = 0.0d+00
              tmp(3,i,j,k) = vper * sin(xp) * tanh(yp) * exp(yt) / kvec
            end do ! i = 1, im
          end do ! j = 1, jm
        end do ! k = 1, km

! calculate velocity components from the potential
!
        call curl(dh(1:3), tmp(1:3,1:im,1:jm,1:km)                             &
                                            , pdata%q(ivx:ivz,1:im,1:jm,1:km))

      end if ! vper /= 0.0

    end if ! pert == 2

! calculate the perturbation of magnetic field
!
    if (ibx > 0) then

! reset magnetic field components
!
      pdata%q(ibx:ibz,:,:,:) = 0.0d+00

! prepare the wave-like perturbation of magnetic field
!
      if (bper /= 0.0d+00 .and. pert >= 3) then

! initiate the vector potential to calculate the divergence-free magnetic field
! components
!
        do k = 1, km
          do j = 1, jm
            do i = 1, im
              xp = kvec * x(i)
              yp = - 0.5d+00 * (y(j) / pth)**2

              tmp(1,i,j,k) = 0.0d+00
              tmp(2,i,j,k) = 0.0d+00
              tmp(3,i,j,k) = vper * cos(xp) * exp(yp) / kvec
            end do ! i = 1, im
          end do ! j = 1, jm
        end do ! k = 1, km

! calculate magnetic field components from vector potential
!
        call curl(dh(1:3), tmp(1:3,1:im,1:jm,1:km)                             &
                                            , pdata%q(ibx:ibz,1:im,1:jm,1:km))

      end if ! bper /= 0.0

    end if ! ibx > 0

! iterate over all positions in the XZ plane
!
    do k = 1, km
      do i = 1, im

! if magnetic field is present, set its initial configuration
!
        if (ibx > 0) then

! set the magnetic field configuration
!
          do j = 1, jm
            yl    = (y(j) - 0.5d+00 * dh(2)) / yth
            yu    = (y(j) + 0.5d+00 * dh(2)) / yth

            itanh = (log_cosh(yu) - log_cosh(yl)) * yrat

            q(ibx,j) = bamp * sign(min(1.0d+00, abs(itanh)), itanh)
            q(iby,j) = 0.0d+00
            q(ibz,j) = (sqrt(2.0d+00 * pmag - q(ibx,j)**2) - bgui) * zeta + bgui
            q(ibp,j) = 0.0d+00
          end do ! j = 1, jm

! calculate the local magnetic pressure
!
          pm(1:jm)    = 0.5d+00 * sum(q(ibx:ibz,:) * q(ibx:ibz,:), 1)

! add the magnetic field perturbation
!
          q(ibx,1:jm) = q(ibx,1:jm) + pdata%q(ibx,i,1:jm,k)
          q(iby,1:jm) = q(iby,1:jm) + pdata%q(iby,i,1:jm,k)
          q(ibz,1:jm) = q(ibz,1:jm) + pdata%q(ibz,i,1:jm,k)

! set the density and pressure profiles
!
          if (ipr > 0) then
            q(idn,1:jm) = dens
            q(ipr,1:jm) = pres + (pmag - pm(1:jm))
          else
            q(idn,1:jm) = dens + (pmag - pm(1:jm)) / csnd2
          end if

        else ! ibx > 0

          if (ipr > 0) then
            q(idn,1:jm) = dens
            q(ipr,1:jm) = pres
          else
            q(idn,1:jm) = dens
          end if

        end if ! ibx > 0

! add the velocity field perturbation
!
        q(ivx,1:jm) = pdata%q(ivx,i,1:jm,k)
        q(ivy,1:jm) = pdata%q(ivy,i,1:jm,k)
        q(ivz,1:jm) = pdata%q(ivz,i,1:jm,k)

! convert the primitive variables to the conservative ones
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
! subroutine LOG_COSH:
! -------------------
!
!   Function calculates the logarithm of the hyperbolic cosine, which is
!   the result of the integration of tanh(x). Direct calculation using
!   Fortran intrinsic subroutines fails for large values of x, therefore
!   the logarithm of cosh is approximated as |x| + log(1/2) for
!   |x| > threshold.
!
!   Arguments:
!
!     x - function argument;
!
!===============================================================================
!
  function log_cosh(x) result(y)

! local variables are not implicit by default
!
    implicit none

! function arguments
!
    real(kind=8), intent(in) :: x
    real(kind=8)             :: y

! local parameters
!
    real(kind=8), parameter :: th = acosh(huge(x)), lh = log(0.5d+00)
!
!-------------------------------------------------------------------------------
!
    if (abs(x) < th) then
      y = log(cosh(x))
    else
      y = abs(x) + lh
    end if

!-------------------------------------------------------------------------------
!
  end function log_cosh
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
    use coordinates    , only : ib, ibl, ie, ieu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp

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

! local variables
!
    integer      :: im2, im1, i, ip1, ip2
    integer      :: jm2, jm1, j, jp1, jp2
    integer      :: km2, km1, k, kp1, kp2
    real(kind=8) :: dx, dy, dz, dxy, dxz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

! process case with magnetic field, otherwise revert to standard outflow
!
    if (ibx > 0) then

! get the cell sizes and their ratios
!
      dx  = x(2) - x(1)
      dy  = y(2) - y(1)
#if NDIMS == 3
      dz  = z(2) - z(1)
#endif /* NDIMS == 3 */
      dxy = dx / dy
#if NDIMS == 3
      dxz = dx / dz
#endif /* NDIMS == 3 */

! process left and right side boundary separatelly
!
      if (ic == 1) then

! iterate over left-side ghost layers
!
        do i = ibl, 1, -1

! calculate neighbor cell indices
!
          ip1 = min(im, i + 1)
          ip2 = min(im, i + 2)

! iterate over boundary layer
!
          do k = kl, ku
#if NDIMS == 3
            km2 = max( 1, k - 2)
            km1 = max( 1, k - 1)
            kp1 = min(km, k + 1)
            kp2 = min(km, k + 2)
#endif /* NDIMS == 3 */
            do j = jl, ju
              jm2 = max( 1, j - 2)
              jm1 = max( 1, j - 1)
              jp1 = min(jm, j + 1)
              jp2 = min(jm, j + 2)

! make the normal derivative zero
!
              qn(1:nv,i,j,k) = qn(1:nv,ib,j,k)

! prevent the inflow
!
              qn(ivx,i,j,k) = min(0.0d+00, qn(ivx,ib,j,k))

! update the normal component of magnetic field from divergence-free condition
!
              qn(ibx,i,j,k) = qn(ibx,ip2,j,k)                                  &
                            + (qn(iby,ip1,jp1,k) - qn(iby,ip1,jm1,k)) * dxy
#if NDIMS == 3
              qn(ibx,i,j,k) = qn(ibx,i  ,j,k)                                  &
                            + (qn(ibz,ip1,j,kp1) - qn(ibz,ip1,j,km1)) * dxz
#endif /* NDIMS == 3 */
              qn(ibp,i,j,k) = 0.0d+00
            end do ! j = jl, ju
          end do ! k = kl, ku
        end do ! i = ibl, 1, -1

      else ! ic == 1

! iterate over right-side ghost layers
!
        do i = ieu, im

! calculate neighbor cell indices
!
          im1 = max( 1, i - 1)
          im2 = max( 1, i - 2)

! iterate over boundary layer
!
          do k = kl, ku
#if NDIMS == 3
            km1 = max( 1, k - 1)
            kp1 = min(km, k + 1)
            km2 = max( 1, k - 2)
            kp2 = min(km, k + 2)
#endif /* NDIMS == 3 */
            do j = jl, ju
              jm1 = max( 1, j - 1)
              jp1 = min(jm, j + 1)
              jm2 = max( 1, j - 2)
              jp2 = min(jm, j + 2)

! make the normal derivative zero
!
              qn(1:nv,i,j,k) = qn(1:nv,ie,j,k)

! prevent the inflow
!
              qn(ivx,i,j,k) = max(0.0d+00, qn(ivx,ie,j,k))

! update the normal component of magnetic field from divergence-free condition
!
              qn(ibx,i,j,k) = qn(ibx,im2,j,k)                              &
                            + (qn(iby,im1,jm1,k) - qn(iby,im1,jp1,k)) * dxy
#if NDIMS == 3
              qn(ibx,i,j,k) = qn(ibx,i  ,j,k)                              &
                            + (qn(ibz,im1,j,km1) - qn(ibz,im1,j,kp1)) * dxz
#endif /* NDIMS == 3 */
              qn(ibp,i,j,k) = 0.0d+00
            end do ! j = jl, ju
          end do ! k = kl, ku
        end do ! i = ieu, im
      end if ! ic == 1
    else ! ibx > 0
      if (ic == 1) then
        do i = ibl, 1, -1
          qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ib,jl:ju,kl:ku)
          qn(ivx ,i,jl:ju,kl:ku) = min(0.0d+00, qn(ivx,ib,jl:ju,kl:ku))
        end do ! i = ibl, 1, -1
      else
        do i = ieu, im
          qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ie,jl:ju,kl:ku)
          qn(ivx ,i,jl:ju,kl:ku) = max(0.0d+00, qn(ivx,ie,jl:ju,kl:ku))
        end do ! i = ieu, im
      end if
    end if ! ibx > 0

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
    use coordinates    , only : jb, jbl, je, jeu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp

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

! local variables
!
    integer      :: im2, im1, i, ip1, ip2
    integer      :: jm2, jm1, j, jp1, jp2
    integer      :: km2, km1, k, kp1, kp2
    real(kind=8) :: dx, dy, dz, dyx, dyz
    real(kind=8) :: fl, fr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

! process case with magnetic field, otherwise revert to standard outflow
!
    if (ibx > 0) then

! get the cell sizes and their ratios
!
      dx  = x(2) - x(1)
      dy  = y(2) - y(1)
#if NDIMS == 3
      dz  = z(2) - z(1)
#endif /* NDIMS == 3 */
      dyx = dy / dx
#if NDIMS == 3
      dyz = dy / dz
#endif /* NDIMS == 3 */

! process left and right side boundary separatelly
!
      if (jc == 1) then

! iterate over left-side ghost layers
!
        do j = jbl, 1, -1

! calculate neighbor cell indices
!
          jp1 = min(jm, j + 1)
          jp2 = min(jm, j + 2)

! calculate variable decay coefficients
!
          fr  = (dy * (jb - j - 5.0d-01)) / blim
          fl  = 1.0d+00 - fr

! iterate over boundary layer
!
          do k = kl, ku
#if NDIMS == 3
            km1 = max( 1, k - 1)
            kp1 = min(km, k + 1)
#endif /* NDIMS == 3 */
            do i = il, iu
              im1 = max( 1, i - 1)
              ip1 = min(im, i + 1)

! make normal derivatives zero
!
              qn(1:nv,i,j,k) = qn(1:nv,i,jb,k)

! decay density and pressure to their limits
!
              qn(idn,i,j,k) = fl * qn(idn,i,jb,k) + fr * dens
              if (ipr > 0) qn(ipr,i,j,k) = fl * qn(ipr,i,jb,k) + fr * pres

! decay magnetic field to its limit
!
              qn(ibx,i,j,k) = fl * qn(ibx,i,jb,k) - fr * bamp
              qn(ibz,i,j,k) = fl * qn(ibz,i,jb,k) + fr * bgui

! update By from div(B)=0
!
              qn(iby,i,j,k) = qn(iby,i,jp2,k)                              &
                           + (qn(ibx,ip1,jp1,k) - qn(ibx,im1,jp1,k)) * dyx
#if NDIMS == 3
              qn(iby,i,j,k) = qn(iby,i,j  ,k)                              &
                           + (qn(ibz,i,jp1,kp1) - qn(ibz,i,jp1,km1)) * dyz
#endif /* NDIMS == 3 */
              qn(ibp,i,j,k) = 0.0d+00
            end do ! i = il, iu
          end do ! k = kl, ku
        end do ! j = jbl, 1, -1
      else ! jc = 1

! iterate over right-side ghost layers
!
        do j = jeu, jm

! calculate neighbor cell indices
!
          jm1 = max( 1, j - 1)
          jm2 = max( 1, j - 2)

! calculate variable decay coefficients
!
          fr  = (dy * (j - je - 5.0d-01)) / blim
          fl  = 1.0d+00 - fr

! iterate over boundary layer
!
          do k = kl, ku
#if NDIMS == 3
            km1 = max( 1, k - 1)
            kp1 = min(km, k + 1)
#endif /* NDIMS == 3 */
            do i = il, iu
              im1 = max( 1, i - 1)
              ip1 = min(im, i + 1)

! make normal derivatives zero
!
              qn(1:nv,i,j,k) = qn(1:nv,i,je,k)

! decay density and pressure to their limits
!
              qn(idn,i,j,k) = fl * qn(idn,i,je,k) + fr * dens
              if (ipr > 0) qn(ipr,i,j,k) = fl * qn(ipr,i,je,k) + fr * pres

! decay magnetic field to its limit
!
              qn(ibx,i,j,k) = fl * qn(ibx,i,je,k) + fr * bamp
              qn(ibz,i,j,k) = fl * qn(ibz,i,je,k) + fr * bgui

! update By from div(B)=0
!
              qn(iby,i,j,k) = qn(iby,i,jm2,k)                              &
                           + (qn(ibx,im1,jm1,k) - qn(ibx,ip1,jm1,k)) * dyx
#if NDIMS == 3
              qn(iby,i,j,k) = qn(iby,i,j  ,k)                              &
                           + (qn(ibz,i,jm1,km1) - qn(ibz,i,jm1,kp1)) * dyz
#endif /* NDIMS == 3 */
              qn(ibp,i,j,k) = 0.0d+00
            end do ! i = il, iu
          end do ! k = kl, ku
        end do ! j = jeu, jm
      end if ! jc = 1
    else ! ibx > 0
      if (jc == 1) then
        do j = jbl, 1, -1
          qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,jb,kl:ku)
          qn(ivy ,il:iu,j,kl:ku) = min(0.0d+00, qn(ivy,il:iu,jb,kl:ku))
        end do ! j = jbl, 1, -1
      else
        do j = jeu, jm
          qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,je,kl:ku)
          qn(ivy ,il:iu,j,kl:ku) = max(0.0d+00, qn(ivy,il:iu,je,kl:ku))
        end do ! j = jeu, jm
      end if
    end if ! ibx > 0

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
