!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2015 Grzegorz Kowal <grzegorz@amuncode.org>
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

! binaries problems
!
      case("binaries")
        update_shapes => update_shapes_binaries

! no shape update
!
      case default
        update_shapes => update_shapes_none

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
!
!===============================================================================
!
  subroutine update_shapes_none(pdata, time)

! include external procedures and variables
!
    use blocks         , only : block_data

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time
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
!
!===============================================================================
!
  subroutine update_shapes_blast(pdata, time)

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
    real(kind=8)             , intent(in)    :: time

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
!
!===============================================================================
!
! subroutine UPDATE_SHAPES_BINARIES:
! ---------------------------------
!
!   Subroutine resets the primitive and conserved variables within a defined
!   shapes for the binaries problem.
!
!   Arguments:
!
!     pdata - pointer to the data block structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine update_shapes_binaries(pdata, time)

! include external procedures and variables
!
    use blocks         , only : block_data
    use constants      , only : d2r, pi2
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az, adx, ady, adz, advol
    use equations      , only : prim2cons
    use equations      , only : gamma
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters     , only : get_parameter_real, get_parameter_integer

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time

! default parameter values
!
    real(kind=8), save :: mstar    = 2.00d+00
    real(kind=8), save :: mcomp    = 1.00d+00
    real(kind=8), save :: rstar    = 5.00d-02
    real(kind=8), save :: rcomp    = 5.00d-02
    real(kind=8), save :: dstar    = 1.00d+06
    real(kind=8), save :: dratio   = 1.00d+02
    real(kind=8), save :: msstar   = 3.00d+00
    real(kind=8), save :: mscomp   = 1.50d+01
    real(kind=8), save :: tstar    = 0.00d+00
    real(kind=8), save :: tcomp    = 0.00d+00
    real(kind=8), save :: omstar   = 0.00d+00
    real(kind=8), save :: omcomp   = 0.00d+00
    real(kind=8), save :: vstar    = 5.00d-01
    real(kind=8), save :: vcomp    = 2.50d+00
    real(kind=8), save :: dist     = 2.00d-01
    real(kind=8), save :: period   = 1.00d+02
    real(kind=8), save :: ecc      = 9.00d-01
    real(kind=8), save :: buni     = 1.00d-03
    real(kind=8), save :: tol      = 1.00d-14
    integer     , save :: maxit    = 20

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: ms, mc
    real(kind=8), save :: dcomp, pstar, pcomp
    real(kind=8), save :: r2star, r2comp
    real(kind=8), save :: acomp , bcomp
    real(kind=8), save :: omega, man, ean, vs, vc
    real(kind=8), save :: xps, yps, xpc, ypc, uxs, uys, uxc, uyc
    real(kind=8), save :: tprev

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: xs, ys, zs
    real(kind=8)  :: xc, yc, zc
    real(kind=8)  :: dv
    real(kind=8)  :: sn, cs, res, om
    real(kind=8)  :: rs2, rc2, rs, rc, rd
    real(kind=8)  :: dns, prs, vxs, vys, vzs
    real(kind=8)  :: dnc, prc, vxc, vyc, vzc

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
! start accounting time for the shape update
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! star masses, radia, densities and sonic Mach numbers
!
      call get_parameter_real("mstar"       , mstar )
      call get_parameter_real("mcomp"       , mcomp )
      call get_parameter_real("rstar"       , rstar )
      call get_parameter_real("rcomp"       , rcomp )
      call get_parameter_real("dstar"       , dstar )
      call get_parameter_real("dratio"      , dratio)
      call get_parameter_real("msstar"      , msstar)
      call get_parameter_real("mscomp"      , mscomp)

! wind speeds
!
      call get_parameter_real("vstar"       , vstar )
      call get_parameter_real("vcomp"       , vcomp )

! orbit parameters
!
      call get_parameter_real("distance"    , dist  )
      call get_parameter_real("period"      , period)
      call get_parameter_real("eccentricity", ecc   )

! star rotation periods
!
      call get_parameter_real("tstar"       , tstar )
      call get_parameter_real("tcomp"       , tcomp )

! magnetic field profile parameters
!
      call get_parameter_real("buni"        , buni  )

! Kepler's equation solver parameters
!
      call get_parameter_real("tolerance", tol   )
      call get_parameter_integer("maxit" , maxit )

! calculate mass fractions
!
      ms = mcomp / (mstar + mcomp)
      mc = mstar / (mstar + mcomp)

! calculate the square of radia
!
      r2star = rstar * rstar
      r2comp = rcomp * rcomp

! calculate densities and pressures
!
      dcomp = dstar / dratio
      if (ipr > 0) then
        pstar = (vstar / msstar)**2 * dstar / gamma
        pcomp = (vcomp / mscomp)**2 * dcomp / gamma
      end if

! calculate orbit parameters
!
      acomp  = dist / (1.0d+00 - ecc)
      bcomp  = acomp * sqrt(1.0d+00 - ecc * ecc)
      omega  = pi2 / period
      ean    = 0.0d+00

! calculate the rotation speed of the stars
!
      if (tstar > 0.0d+00) omstar = 1.0d+00 / tstar
      if (tcomp > 0.0d+00) omcomp = 1.0d+00 / tcomp

! set the previous time
!
      tprev  = -1.0d+00

! reset the first execution flag
!
      first = .false.

    end if ! first call

! if time changes, we have to recalculate the star positions
!
    if (time /= tprev) then

! solve the Kepler's equation to obtain the new value of eccentric anomaly
!
      man = omega * time
      res = 1.0d+00
      i   = 1
      do while(abs(res) >= tol .and. i <= maxit)
        res = (ean - ecc * sin(ean) - man) / (1.0d+00 - ecc * cos(ean))
        ean = ean - res
        i   = i + 1
      end do
      if (abs(res) >= tol) then
        print *, "Kepler's equations could not be solved!"
        print *, "it: ", i, ", E: ", ean, ", dE: ", res
        print *, ""
      end if

! calculate trigonometric functions of the true anomaly
!
      dv = 1.0d+00 - ecc * cos(ean)
      sn = sqrt(1.0d+00 - ecc**2) * sin(ean) / dv
      cs = (cos(ean) - ecc) / dv

! calculate the angular velocity
!
      om = omega / (1.0d+00 - ecc * cs)

! calculate the position and velocity of the companion star
!
      rd  =   acomp * dv
      rs  =   ms * rd
      rc  =   mc * rd
      vs  =   rs * om
      vc  =   rc * om
      xps = - rs * cs
      yps = - rs * sn
      xpc =   rc * cs
      ypc =   rc * sn
      uxs =   vs * sn
      uys = - vs * cs
      uxc = - vc * sn
      uyc =   vc * cs

! update tprev
!
      tprev = time

    end if ! time /= tprev

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
#endif /* NDIMS == 3 */

! iterate over all positions in the YZ plane
!
    do k = 1, km

! calculate the Z coordinates of the central and companion stars
!
      zs = z(k)
      zc = z(k)

      do j = 1, jm

! calculate the Y coordinates of the central and companion stars
!
        ys = y(j) - yps
        yc = y(j) - ypc

! copy the primitive variable vector
!
        q(1:nv,1:im) = pdata%q(1:nv,1:im,j,k)

! sweep along the X coordinate
!
        do i = 1, im

! calculate the X coordinates of the central and companion stars
!
          xs = x(i) - xps
          xc = x(i) - xpc

! calculate the distances from the centers of the central and companion stars
!
          rs2 = max(xs * xs + ys * ys + zs * zs, 1.0d-16)
          rc2 = max(xc * xc + yc * yc + zc * zc, 1.0d-16)
          rs  = sqrt(rs2)
          rc  = sqrt(rc2)

! calculate profiles from the central star
!
#if NDIMS == 3
          rd = max(r2star, rs2) / r2star
#else /* NDIMS == 3 */
          rd = max(rstar , rs ) / rstar
#endif /* NDIMS == 3 */

          dns = dstar / rd
          if (ipr > 0) prs = pstar / rd

! set the central star wind velocity
!
          vxs = vstar * xs / rs
          vys = vstar * ys / rs
          vzs = vstar * zs / rs

! calculate profiles from the companion star
!
#if NDIMS == 3
          rd = max(r2comp, rc2) / r2comp
#else /* NDIMS == 3 */
          rd = max(rcomp , rc ) / rcomp
#endif /* NDIMS == 3 */

          dnc = dcomp / rd
          if (ipr > 0) prc = pcomp / rd

! set the companion star wind velocity
!
          vxc = vcomp * xc / rc
          vyc = vcomp * yc / rc
          vzc = vcomp * zc / rc

! set the variables
!
          if (rs2 <= r2star) then
            q(idn,i) = dns
            q(ivx,i) = vxs - omstar * ys / rs + uxs
            q(ivy,i) = vys + omstar * xs / rs + uys
            q(ivz,i) = vzs
            if (ipr > 0) q(ipr,i) = prs
            if (ibx > 0) then
              q(ibx,i) = 0.0d+00
              q(iby,i) = 0.0d+00
              q(ibz,i) = buni
              q(ibp,i) = 0.0d+00
            end if
          else if (rc2 <= r2comp) then
            q(idn,i) = dnc
            q(ivx,i) = vxc - omcomp * yc + uxc
            q(ivy,i) = vyc + omcomp * xc + uyc
            q(ivz,i) = vzc
            if (ipr > 0) q(ipr,i) = prc
            if (ibx > 0) then
              q(ibx,i) = 0.0d+00
              q(iby,i) = 0.0d+00
              q(ibz,i) = buni
              q(ibp,i) = 0.0d+00
            end if
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
  end subroutine update_shapes_binaries

!===============================================================================
!
end module shapes
