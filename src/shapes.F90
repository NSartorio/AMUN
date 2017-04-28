!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2017 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: SHAPES
!!
!!  This modules handles the update of the shapes embedded in the domain. In
!!  such a region, the primitive and conservative variables have to be updated
!!  after each temporal integration substep.
!!
!!******************************************************************************
!
module shapes

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

! pointer to the shape update subroutine
!
  procedure(update_shapes_none), pointer, save :: update_shapes => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_shapes, finalize_shapes
  public :: update_shapes

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
! subroutine INITIALIZE_SHAPES:
! ----------------------------
!
!   Subroutine initializes module SHAPES.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_shapes(verbose, iret)

! include external procedures and variables
!
    use parameters     , only : get_parameter_string
    use user_problem   , only : update_shapes_user

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64) :: problem_name  = "blast"
    character(len=64) :: enable_shapes = "off"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('shapes:: initialize', imi)
    call set_timer('shapes:: update'    , imu)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the problem name
!
    call get_parameter_string("problem"      , problem_name )
    call get_parameter_string("enable_shapes", enable_shapes)

! set the correct procedure pointer if shapes are enabled
!
    select case(trim(enable_shapes))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")

! select the shape update subroutine depending on the problem
!
      select case(trim(problem_name))

! general test problems
!
      case("blast")
        update_shapes => update_shapes_blast

! no shape update
!
      case default
        update_shapes => update_shapes_user

      end select

    case default

! by default the shape update is turned off, so reset the procedure pointer
!
      update_shapes => update_shapes_none

    end select

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,a)") "embedded shapes        =", trim(enable_shapes)

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_shapes
!
!===============================================================================
!
! subroutine FINALIZE_SHAPES:
! --------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_shapes(iret)

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

! nullify procedure pointers
!
    nullify(update_shapes)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_shapes
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine UPDATE_SHAPES_NONE:
! -----------------------------
!
!   Subroutine does not do anything, but it is required to define the interface
!   for other specific shape update subroutines.
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
  subroutine update_shapes_none(pdata, time, dt)

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
    call start_timer(imu)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the shape update
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_shapes_none
!
!===============================================================================
!
! subroutine UPDATE_SHAPES_BLAST:
! ------------------------------
!
!   Subroutine resets the primitive and conserved variables within a defined
!   shape for the blast problem.
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
  subroutine update_shapes_blast(pdata, time, dt)

! include external procedures and variables
!
    use blocks         , only : block_data
    use constants      , only : d2r
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az, adx, ady, adz
    use equations      , only : prim2cons
    use equations      , only : gamma
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters     , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time, dt

! default parameter values
!
    real(kind=8), save :: dens   = 1.00d+00
    real(kind=8), save :: ratio  = 1.00d+02
    real(kind=8), save :: radius = 1.00d-01
    real(kind=8), save :: csnd   = 4.0824829046386301635d-01
    real(kind=8), save :: buni   = 1.00d+00
    real(kind=8), save :: angle  = 4.50d+01

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: r2
    real(kind=8), save :: dn_ovr, pr_ovr, bx_ovr, by_ovr

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: xl, yl, zl, xu, yu, zu, rl, ru
    real(kind=8)  :: xb, yb, xt, yt
    real(kind=8)  :: dx, dy, dz, dxh, dyh, dzh, daxy
    real(kind=8)  :: fc_amb, fc_ovr

! local arrays
!
    real(kind=8), dimension(nv)    :: qi
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the shape update
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("dens"  , dens  )
      call get_parameter_real("ratio" , ratio )
      call get_parameter_real("radius", radius)
      call get_parameter_real("csnd"  , csnd  )
      call get_parameter_real("buni"  , buni  )
      call get_parameter_real("angle" , angle )

! set the conditions inside the radius
!
      if (ipr > 0) then
        dn_ovr = dens
        pr_ovr = dens * ratio * csnd * csnd / gamma
      else
        dn_ovr = dens * ratio
      end if
      bx_ovr = buni * cos(d2r * angle)
      by_ovr = buni * sin(d2r * angle)

! calculate the square of radius
!
      r2    = radius * radius

! reset the first execution flag
!
      first = .false.

    end if ! first call

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
#endif /* NDIMS == 3 */

! calculate mesh intervals and areas
!
    dx   = adx(pdata%meta%level)
    dy   = ady(pdata%meta%level)
    dz   = adz(pdata%meta%level)
    dxh  = 0.5d+00 * dx
    dyh  = 0.5d+00 * dy
#if NDIMS == 3
    dzh  = 0.5d+00 * dz
#else /* NDIMS == 3 */
    dzh  = 1.0d+00
#endif /* NDIMS == 3 */
    daxy = dx * dy

! set the conditions inside the radius
!
    if (ipr > 0) then
      qi(idn) = dn_ovr
      qi(ipr) = pr_ovr
    else
      qi(idn) = dn_ovr
    end if
    qi(ivx) = 0.0d+00
    qi(ivy) = 0.0d+00
    qi(ivz) = 0.0d+00
    if (ibx > 0) then
      qi(ibx) = bx_ovr
      qi(iby) = by_ovr
      qi(ibz) = 0.0d+00
      qi(ibp) = 0.0d+00
    end if

! iterate over all positions in the YZ plane
!
    do k = 1, km

#if NDIMS == 3
! calculate the corner Z coordinates
!
      zl = abs(z(k)) - dzh
      zu = abs(z(k)) + dzh
#endif /* NDIMS == 3 */

      do j = 1, jm

! calculate the corner Y coordinates
!
        yl = abs(y(j)) - dyh
        yu = abs(y(j)) + dyh

! sweep along the X coordinate
!
        do i = 1, im

! calculate the corner X coordinates
!
          xl = abs(x(i)) - dxh
          xu = abs(x(i)) + dxh

! calculate the minimum and maximum corner distances from the origin
!
#if NDIMS == 3
          rl = xl * xl + yl * yl + zl * zl
          ru = xu * xu + yu * yu + zu * zu
#else /* NDIMS == 3 */
          rl = xl * xl + yl * yl
          ru = xu * xu + yu * yu
#endif /* NDIMS == 3 */

! set the initial density and pressure in cells laying completely within
! the blast radius
!
          if (ru <= r2) then

! set the overpressure region primitive variables
!
            q(1:nv,i) = qi(1:nv)

! variables in cells completely outside the given radius are not changed
!
          else if (rl >= r2) then

            q(1:nv,i) = pdata%q(1:nv,i,j,k)

! integrate density or pressure in cells which are crossed by the circle with
! the given radius
!
          else

#if NDIMS == 3
! in 3D simply set the ambient values since the integration is more complex
!
            q(1:nv,i) = pdata%q(1:nv,i,j,k)
#else /* NDIMS == 3 */
! calculate the bounds of area integration
!
            xb = max(xl, sqrt(max(0.0d+00, r2 - yu * yu)))
            xt = min(xu, sqrt(max(0.0d+00, r2 - yl * yl)))
            yb = max(yl, sqrt(max(0.0d+00, r2 - xu * xu)))
            yt = min(yu, sqrt(max(0.0d+00, r2 - xl * xl)))

! integrate the area below the circle within the current cell for both
! functions, y = f(x) and x = g(y), and then average them to be sure that we
! are getting the ideal symmetry
!
            fc_ovr = 0.5d+00 * (r2 * (asin(xt / radius) - asin(xb / radius))   &
                   + (xt * yb - xb * yt)) - yl * (xt - xb)
            fc_ovr = fc_ovr + (xb - xl) * dy

            fc_amb = 0.5d+00 * (r2 * (asin(yt / radius) - asin(yb / radius))   &
                   + (yt * xb - yb * xt)) - xl * (yt - yb)
            fc_amb = fc_amb + (yb - yl) * dx

            fc_ovr = 0.5d+00 * (fc_ovr + fc_amb)

! normalize coefficients
!
            fc_ovr = fc_ovr / daxy
            fc_amb = 1.0d+00 - fc_ovr

! integrate the primitive variables over the edge cells
!
            q(1:nv,i) = fc_ovr * qi(1:nv) + fc_amb * pdata%q(1:nv,i,j,k)
#endif /* NDIMS == 3 */

          end if

        end do ! i = 1, im

! convert the primitive variables to conservative ones
!
        call prim2cons(im, q(1:nv,1:im), u(1:nv,1:im))

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the shape update
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_shapes_blast

!===============================================================================
!
end module shapes
