!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2017-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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

  real(kind=8), save :: dcomp, pstar, pcomp
  real(kind=8), save :: r2star, r2comp
  real(kind=8), save :: acomp , bcomp
  real(kind=8), save :: xps, xpc, uys, uyc
  real(kind=8), save :: om

  real(kind=8), save :: ms, mc
  real(kind=8), save :: asemi, omega, ean
  real(kind=8), save :: yps, ypc, xrs, yrs, xrc, yrc
  real(kind=8), save :: uxs, uxc
  real(kind=8), save :: tprev

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
    use constants  , only : pi2
    use equations  , only : ipr
    use equations  , only : gamma
    use parameters , only : get_parameter_string, get_parameter_real           &
                          , get_parameter_integer

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
    om     = (pi2 / period) / (1.0d+00 - ecc)

! calculate initial positions and velocities
!
    xps    = - mcomp / (mstar + mcomp) * dist
    xpc    =   mstar / (mstar + mcomp) * dist
    uys    = - xps * om
    uyc    =   xpc * om

! calculate the rotation speeds
!
    if (tstar > 0.0d+00) omstar = 1.0d+00 / tstar
    if (tcomp > 0.0d+00) omcomp = 1.0d+00 / tcomp

! calculate mass fractions
!
    ms = mcomp / (mstar + mcomp)
    mc = mstar / (mstar + mcomp)

! calculate orbit parameters
!
    asemi  = dist / (1.0d+00 - ecc)
    omega  = pi2 / period
    ean    = 0.0d+00

! set the initial positions
!
    xrs = - ms * dist
    yrs =   0.0d+00
    xrc =   mc * dist
    yrc =   0.0d+00

! set the previous time
!
    tprev  = 0.0d+00

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
!   Subroutine sets the initial conditions for the wind interaction in
!   the binary star problem.
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
    use coordinates, only : ax, ay, az
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
    integer       :: i, j, k, ic, jc, kc
    real(kind=8)  :: xs, ys, zs
    real(kind=8)  :: xc, yc, zc
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
! start accounting time for the problem setup
!
    call start_timer(imp)
#endif /* PROFILE */

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
#endif /* NDIMS == 3 */

! set magnetic field components
!
    if (ibx > 0) then

      q(ibx,1:im) = 0.0d+00
      q(iby,1:im) = 0.0d+00
      q(ibz,1:im) = buni
      q(ibp,1:im) = 0.0d+00

    end if

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
        ys = y(j)
        yc = y(j)

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
            q(ivx,i) = vxs - omstar * ys
            q(ivy,i) = vys + omstar * xs + uys
            q(ivz,i) = vzs
            if (ipr > 0) q(ipr,i) = prs
          else if (rc2 <= r2comp) then
            q(idn,i) = dnc
            q(ivx,i) = vxc - omcomp * yc
            q(ivy,i) = vyc + omcomp * xc + uyc
            q(ivz,i) = vzc
            if (ipr > 0) q(ipr,i) = prc
          else
            q(idn,i) = min(dstar, dns + dnc)
            q(ivx,i) = vxs + vxc
            q(ivy,i) = vys + vyc
            q(ivz,i) = vzs + vzc
            if (ipr > 0) q(ipr,i) = min(pstar, prs + prc)
          end if

        end do ! i = 1, im

! convert the primitive variables to conservative ones
!
        call prim2cons(im, q(1:nv,1:im), u(1:nv,1:im))

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
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az!, adx, ady, adz, advol
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use equations      , only : prim2cons

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time, dt

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: xs, ys, zs
    real(kind=8)  :: xc, yc, zc
    real(kind=8)  :: dv, ds
    real(kind=8)  :: sn, cs, man, res
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
    call start_timer(ims)
#endif /* PROFILE */

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

! calculate the position and velocity of the companion star
!
      rd  =   asemi * dv
      rs  =   ms * rd
      rc  =   mc * rd
      xps = - rs * cs
      yps = - rs * sn
      xpc =   rc * cs
      ypc =   rc * sn
      ds  = time - tprev
      uxs = (xps - xrs) / ds
      uys = (yps - yrs) / ds
      uxc = (xpc - xrc) / ds
      uyc = (ypc - yrc) / ds

! update tprev, previous positions
!
      tprev = time
      xrs   = xps
      yrs   = yps
      xrc   = xpc
      yrc   = ypc

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

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

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
