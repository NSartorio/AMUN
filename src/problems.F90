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
!! module: PROBLEMS
!!
!!  This module handles the initialization of various test and research
!!  problems.
!!
!!
!!******************************************************************************
!
module problems

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

! pointer to the problem setup subroutine
!
  procedure(setup_problem_blast), pointer, save :: setup_problem => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_problems, finalize_problems
  public :: setup_problem

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
! subroutine INITIALIZE_PROBLEMS:
! ------------------------------
!
!   Subroutine prepares module PROBLEMS.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_problems(verbose, iret)

! include external procedures and variables
!
    use error          , only : print_error
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
    character(len=64)      :: problem_name  = "blast"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('problems:: initialize', imi)
    call set_timer('problems:: update'    , imu)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the problem name
!
    call get_parameter_string("problem", problem_name)

! associate the setup_problem pointer with the respective problem setup
! subroutine
!
    select case(trim(problem_name))

! general test problems
!
    case("blast")
      setup_problem => setup_problem_blast

    case("implosion")
      setup_problem => setup_problem_implosion

    case("kh", "kelvinhelmholtz", "kelvin-helmholtz")
      setup_problem => setup_problem_kelvin_helmholtz

    case("current_sheet")
      setup_problem => setup_problem_current_sheet

! research problems
!
    case("jet")
      setup_problem => setup_problem_jet

    case("binaries")
      setup_problem => setup_problem_binaries

    case default
      call print_error("problems::initialize_problems()"                       &
                     , "Setup subroutine is not implemented for this problem!")
      iret = 600
    end select

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_problems
!
!===============================================================================
!
! subroutine FINALIZE_PROBLEMS:
! ----------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_problems(iret)

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
    nullify(setup_problem)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_problems
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_BLAST:
! ------------------------------
!
!   Subroutine sets the initial conditions for the spherical blast problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!
!===============================================================================
!
  subroutine setup_problem_blast(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use constants  , only : d2r
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az, adx, ady, adz, advol
    use equations  , only : prim2cons
    use equations  , only : gamma
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real, get_parameter_integer

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real(kind=8), save :: dens     = 1.00d+00
    real(kind=8), save :: ratio    = 1.00d+02
    real(kind=8), save :: radius   = 1.00d-01
    real(kind=8), save :: csnd     = 4.0824829046386301635d-01
    real(kind=8), save :: buni     = 1.00d+00
    real(kind=8), save :: angle    = 4.50d+01
#if NDIMS == 3
    integer     , save :: nsubgrid = 10
#endif /* NDIMS == 3 */

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: dn_amb, dn_ovr
    real(kind=8), save :: pr_amb, pr_ovr
    real(kind=8), save :: bx, by
    real(kind=8), save :: r2

! local variables
!
    integer       :: i, j, k, ic, jc, kc
    real(kind=8)  :: xl, yl, zl, xu, yu, zu, rl, ru
    real(kind=8)  :: sn
#if NDIMS == 3
    real(kind=8)  :: xb, yb, zb
    real(kind=8)  :: xt, yt, zt
    real(kind=8)  :: fc_inc
#else /* NDIMS == 3 */
    real(kind=8)  :: rlu, rul
    real(kind=8)  :: xb, yb
    real(kind=8)  :: xt, yt
    real(kind=8)  :: ph
#endif /* NDIMS == 3 */
    real(kind=8)  :: dx, dy, dz, dxh, dyh, dzh, dvol
    real(kind=8)  :: fc_amb, fc_ovr

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z

#if NDIMS == 3
! allocatable arrays
!
    real(kind=8), dimension(:), allocatable :: xm, ym, zm
    real(kind=8), dimension(:), allocatable :: xp, yp, zp
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
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

#if NDIMS == 3
! get the fine grid resolution
!
      call get_parameter_integer("nsubgrid", nsubgrid)

! correct subgrid resolution if necessary
!
      nsubgrid = max(1, nsubgrid)
#endif /* NDIMS == 3 */

! calculate the overdense and ambient region densities
!
      dn_amb = dens
      if (ipr > 0) then
        dn_ovr = dn_amb

! calculate parallel and perpendicular pressures from sound speeds
!
        pr_amb = dens * csnd * csnd / gamma
        pr_ovr = pr_amb * ratio

      else
        dn_ovr = dn_amb * ratio
      end if

! calculate initial uniform field components
!
      if (ibx > 0) then
        sn = sin(d2r * angle)
        bx = buni * sqrt(1.0d+00 - sn * sn)
        by = buni * sn
      end if

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
    dvol = advol(pdata%meta%level)

#if NDIMS == 3
! allocate subgrid coordinates
!
    allocate(xm(nsubgrid), ym(nsubgrid), zm(nsubgrid))
    allocate(xp(nsubgrid), yp(nsubgrid), zp(nsubgrid))

! and generate them
!
    xm(:) = (1.0d+00 * (/(i, i = 0, nsubgrid - 1)/)) / nsubgrid
    ym(:) = xm(:)
    zm(:) = xm(:)
    xm(:) = xm(:) * dx
    ym(:) = ym(:) * dy
    zm(:) = zm(:) * dz
    xp(:) = (1.0d+00 * (/(i, i = 1, nsubgrid    )/)) / nsubgrid
    yp(:) = xp(:)
    zp(:) = xp(:)
    xp(:) = xp(:) * dx
    yp(:) = yp(:) * dy
    zp(:) = zp(:) * dz

! calculate the factor increment for the given subgrid
!
    fc_inc = dvol / nsubgrid**3
#endif /* NDIMS == 3 */

! set the ambient density and pressure
!
    q(idn,:) = dn_amb
    if (ipr > 0) q(ipr,:) = pr_amb

! reset velocity components
!
    q(ivx,:) = 0.0d+00
    q(ivy,:) = 0.0d+00
    q(ivz,:) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
    if (ibx > 0) then

      q(ibx,:) = bx
      q(iby,:) = by
      q(ibz,:) = 0.0d+00
      q(ibp,:) = 0.0d+00

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

! set the overpressure region density
!
            q(idn,i) = dn_ovr

! set the overpressure region pressure
!
            if (ipr > 0) q(ipr,i) = pr_ovr

! set the initial pressure in the cell completely outside the radius
!
          else if (rl >= r2) then

! set the ambient region density
!
            q(idn,i) = dn_amb

! set the ambient medium pressure
!
            if (ipr > 0) q(ipr,i) = pr_amb

! integrate density or pressure in cells which are crossed by the circule with
! the given radius
!
          else

#if NDIMS == 3
! interpolate the factor using subgrid
!
            fc_ovr = 0.0d+00
            do kc = 1, nsubgrid
              zb = (zl + zm(kc))**2
              zt = (zl + zp(kc))**2
              do jc = 1, nsubgrid
                yb = (yl + ym(jc))**2
                yt = (yl + yp(jc))**2
                do ic = 1, nsubgrid
                  xb = (xl + xm(ic))**2
                  xt = (xl + xp(ic))**2

! update the integration factor depending on the subcell position
!
                  if ((xt + yt + zt) <= r2) then
                    fc_ovr = fc_ovr + fc_inc
                  else if ((xb + yb + zb) < r2) then
                    fc_ovr = fc_ovr + 0.5d+00 * fc_inc
                  end if

                end do ! ic = 1, nsubgrid
              end do ! jc = 1, nsubgrid
            end do ! kc = 1, nsubgrid
#else /* NDIMS == 3 */

! calculate the distance of remaining corners
!
            rlu = xl * xl + yu * yu
            rul = xu * xu + yl * yl

! separate in the cases of which corners lay inside, and which outside
! the radius
!
            if (min(rlu, rul) >= r2) then

! only one cell corner inside the radius
!
! calculate middle coordinates of the radius-edge crossing point
!
              xb = sqrt(r2 - yl**2) - xl
              yb = sqrt(r2 - xl**2) - yl

! calculate the sin(½φ), φ, and sin(φ)
!
              sn = 0.5d+00 * sqrt(xb**2 + yb**2) / radius
              ph = 2.0d+00 * asin(sn)
              sn = sin(ph)

! calculate the area of cell intersection with the radius
!
              fc_ovr = 0.5d+00 * (xb * yb + (ph - sn) * r2)

            else if (rlu >= r2) then

! two lower corners inside the radius
!
! calculate middle coordinates of the radius-edge crossing point
!
              yb = sqrt(r2 - xl**2) - yl
              yt = sqrt(r2 - xu**2) - yl

! calculate the sin(½φ), φ, and sin(φ)
!
              sn = 0.5d+00 * sqrt(dx**2 + (yt - yb)**2) / radius
              ph = 2.0d+00 * asin(sn)
              sn = sin(ph)

! calculate the area of cell intersection with the radius
!
              fc_ovr = 0.5d+00 * ((yt + yb) * dx + (ph - sn) * r2)

            else if (rul >= r2) then

! two left corners inside the radius
!
! calculate middle coordinates of the radius-edge crossing point
!
              xb = sqrt(r2 - yl**2) - xl
              xt = sqrt(r2 - yu**2) - xl

! calculate the sin(½φ), φ, and sin(φ)
!
              sn = 0.5d+00 * sqrt((xt - xb)**2 + dy**2) / radius
              ph = 2.0d+00 * asin(sn)
              sn = sin(ph)

! calculate the area of cell intersection with the radius
!
              fc_ovr = 0.5d+00 * ((xt + xb) * dy + (ph - sn) * r2)

            else

! three corners inside the radius
!
! calculate middle coordinates of the radius-edge crossing point
!
              xt = xu - sqrt(r2 - yu**2)
              yt = yu - sqrt(r2 - xu**2)

! calculate the sin(½φ), φ, and sin(φ)
!
              sn = 0.5d+00 * sqrt(xt**2 + yt**2) / radius
              ph = 2.0d+00 * asin(sn)
              sn = sin(ph)

! calculate the area of cell intersection with the radius
!
              fc_ovr = dvol - 0.5d+00 * (xt * yt - (ph - sn) * r2)

            end if
#endif /* NDIMS == 3 */

! normalize coefficients
!
            fc_ovr = fc_ovr / dvol
            fc_amb = 1.0d+00 - fc_ovr

! integrate the density over the edge cells
!
            q(idn,i) = fc_ovr * dn_ovr + fc_amb * dn_amb

! integrate the pressure over the edge cells
!
            if (ipr > 0) q(ipr,i) = fc_ovr * pr_ovr + fc_amb * pr_amb

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

#if NDIMS == 3
! deallocate subgrid coordinates
!
    deallocate(xm, ym, zm)
    deallocate(xp, yp, zp)
#endif /* NDIMS == 3 */

#ifdef PROFILE
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_blast
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_IMPLOSION:
! ----------------------------------
!
!   Subroutine sets the initial conditions for the implosion problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_implosion(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data, ndims
    use constants  , only : d2r
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az, adx, ady, adz
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real(kind=8), save :: sline  = 1.50d-01
    real(kind=8), save :: adens  = 1.00d+00
    real(kind=8), save :: apres  = 1.00d+00
    real(kind=8), save :: drat   = 1.25d-01
    real(kind=8), save :: prat   = 1.40d-01
    real(kind=8), save :: buni   = 1.00d+00
    real(kind=8), save :: bgui   = 0.00d+00
    real(kind=8), save :: angle  = 0.00d+00

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: odens = 1.25d-01
    real(kind=8), save :: opres = 1.40d-01

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: rl, ru, dx, dy, dz, dxh, dyh, dzh, ds, dl, dr
    real(kind=8)  :: sn, cs

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x, xl, xu
    real(kind=8), dimension(jm)    :: y, yl, yu
    real(kind=8), dimension(km)    :: z, zl, zu
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("shock_line"      , sline )
      call get_parameter_real("ambient_density" , adens )
      call get_parameter_real("ambient_pressure", apres )
      call get_parameter_real("density_ratio"   , drat  )
      call get_parameter_real("pressure_ratio"  , prat  )
      call get_parameter_real("buni"            , buni  )
      call get_parameter_real("bgui"            , bgui  )
      call get_parameter_real("angle"           , angle )

! calculate parameters
!
      odens = drat * adens
      opres = prat * apres

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
#endif /* NDIMS == 3 */

! calculate mesh intervals and areas
!
    dx  = adx(pdata%meta%level)
    dy  = ady(pdata%meta%level)
#if NDIMS == 3
    dz  = adz(pdata%meta%level)
#endif /* NDIMS == 3 */
    dxh  = 0.5d+00 * dx
    dyh  = 0.5d+00 * dy
#if NDIMS == 3
    dzh  = 0.5d+00 * dz
#endif /* NDIMS == 3 */

! calculate edge coordinates
!
    xl(:) = abs(x(:)) - dxh
    xu(:) = abs(x(:)) + dxh
    yl(:) = abs(y(:)) - dyh
    yu(:) = abs(y(:)) + dyh
#if NDIMS == 3
    zl(:) = abs(z(:)) - dzh
    zu(:) = abs(z(:)) + dzh
#endif /* NDIMS == 3 */

! reset velocity components
!
    q(ivx,:) = 0.0d+00
    q(ivy,:) = 0.0d+00
    q(ivz,:) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
    if (ibx > 0) then

! calculate the orientation angles
!
      sn = sin(d2r * angle)
      cs = sqrt(1.0d+00 - sn * sn)

! set magnetic field components
!
      q(ibx,:) = buni * cs
      q(iby,:) = buni * sn
      q(ibz,:) = bgui
      q(ibp,:) = 0.0d+00

    end if

! iterate over all positions
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im

! calculate the distance from the origin
!
#if NDIMS == 3
          rl = xl(i) + yl(j) + zl(k)
          ru = xu(i) + yu(j) + zu(k)
#else /* NDIMS == 3 */
          rl = xl(i) + yl(j)
          ru = xu(i) + yu(j)
#endif /* NDIMS == 3 */

! initialize density and pressure
!
          if (ru <= sline) then
            q(idn,i) = odens
            if (ipr > 0) q(ipr,i) = opres
          else if (rl >= sline) then
            q(idn,i) = adens
            if (ipr > 0) q(ipr,i) = apres
          else
            ds = (sline - rl) / dx
            if (ds <= 1.0d+00) then
              dl = 5.0d-01 * ds**ndims
              dr = 1.0d+00 - dl
            else
              ds = (ru - sline) / dx
              dr = 5.0d-01 * ds**ndims
              dl = 1.0d+00 - dr
            end if

            q(idn,i) = adens * dl + odens * dr
            if (ipr > 0) q(ipr,i) = apres * dl + opres * dr
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
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_implosion
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_KELVIN_HELMHOLTZ:
! -----------------------------------------
!
!   Subroutine sets the initial conditions for the Kelvin-Helmholtz instability
!   problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_kelvin_helmholtz(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use constants  , only : d2r
    use coordinates, only : im, jm, km
    use coordinates, only : ay, ady
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real
    use random     , only : randomn

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real(kind=8), save :: ycut   = 2.50d-01
    real(kind=8), save :: dens   = 1.00d+00
    real(kind=8), save :: drat   = 2.00d+00
    real(kind=8), save :: pres   = 2.50d+00
    real(kind=8), save :: vamp   = 1.00d+00
    real(kind=8), save :: vper   = 1.00d-02
    real(kind=8), save :: buni   = 1.00d+00
    real(kind=8), save :: bgui   = 0.00d+00
    real(kind=8), save :: angle  = 0.00d+00

! local saved parameters
!
    logical     , save :: first = .true.

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: yl, yu, dy, dyh
    real(kind=8)  :: sn, cs

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
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("ycut"  , ycut  )
      call get_parameter_real("dens"  , dens  )
      call get_parameter_real("drat"  , drat  )
      call get_parameter_real("pres"  , pres  )
      call get_parameter_real("vamp"  , vamp  )
      call get_parameter_real("vper"  , vper  )
      call get_parameter_real("buni"  , buni  )
      call get_parameter_real("bgui"  , bgui  )
      call get_parameter_real("angle" , angle )

! reset the first execution flag
!
      first = .false.

    end if ! first call

! prepare block coordinates
!
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)

! calculate mesh intervals and areas
!
    dy   = ady(pdata%meta%level)
    dyh  = 0.5d+00 * dy

! set the ambient density and pressure
!
    q(idn,:) = dens
    if (ipr > 0) q(ipr,:) = pres

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
    if (ibx > 0) then

! calculate the orientation angles
!
      sn = sin(d2r * angle)
      cs = sqrt(1.0d+00 - sn * sn)

! set magnetic field components
!
      q(ibx,:) = buni * cs
      q(iby,:) = buni * sn
      q(ibz,:) = bgui
      q(ibp,:) = 0.0d+00

    end if

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do j = 1, jm

! calculate the corner Y coordinates
!
        yl = abs(y(j)) - dyh
        yu = abs(y(j)) + dyh

! set the primitive variables for two regions
!
        if (yu <= ycut) then
          q(idn,1:im) =   dens * drat
          q(ivx,1:im) =   vamp
        else if (yl >= ycut) then
          q(idn,1:im) =   dens
          q(ivx,1:im) = - vamp
        else
          q(idn,1:im) = dens * ((yu - ycut) * drat + (ycut - yl)) / dy
          q(ivx,1:im) = vamp * ((yu - ycut)        - (ycut - yl)) / dy
        end if

! reset remaining velocity components
!
        q(ivy,1:im) = 0.0d+00
        q(ivz,1:im) = 0.0d+00

! set the pressure
!
        if (ipr > 0) q(ipr,:) = pres

! add a random seed velocity component
!
        do i = 1, im
          q(ivx,i) = q(ivx,i) + vper * randomn()
          q(ivy,i) = q(ivy,i) + vper * randomn()
#if NDIMS == 3
          q(ivz,i) = q(ivz,i) + vper * randomn()
#endif /* NDIMS == 3 */
        end do

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
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_kelvin_helmholtz
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_CURRENT_SHEET:
! --------------------------------------
!
!   Subroutine sets the initial conditions for the current sheet test problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_current_sheet(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use constants  , only : pi2
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, adx
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real(kind=8), save :: xcut  = 2.50d-01
    real(kind=8), save :: dens  = 1.00d+00
    real(kind=8), save :: beta  = 1.00d-01
    real(kind=8), save :: vamp  = 1.00d-01
    real(kind=8), save :: buni  = 1.00d+00
    real(kind=8), save :: bgui  = 0.00d+00

! local saved parameters
!
    logical     , save :: first = .true.

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: xl, xu, dxh

! local arrays
!
    real(kind=8), dimension(nv,jm) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("xcut", xcut)
      call get_parameter_real("dens", dens)
      call get_parameter_real("beta", beta)
      call get_parameter_real("vamp", vamp)
      call get_parameter_real("buni", buni)
      call get_parameter_real("bgui", bgui)

! reset the first execution flag
!
      first = .false.

    end if ! first call

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)

! calculate mesh intervals and areas
!
    dxh  = 0.5d+00 * adx(pdata%meta%level)

! set the ambient density and pressure
!
    q(idn,1:jm) = dens
    if (ipr > 0) q(ipr,1:jm) = 0.5d+00 * beta

! set initial velocity
!
    q(ivx,1:jm) = vamp * sin(pi2 * y(1:jm))
    q(ivy,1:jm) = 0.0d+00
    q(ivz,1:jm) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
    if (ibx > 0) then

! set magnetic field components
!
      q(ibx,:) = 0.0d+00
      q(iby,:) = 0.0d+00
      q(ibz,:) = bgui
      q(ibp,:) = 0.0d+00

    end if

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do i = 1, im

! calculate the corner Y coordinates
!
        xl = abs(x(i)) - dxh
        xu = abs(x(i)) + dxh

! set two regions of magnetic field
!
        if (xu <= xcut) then
          if (iby > 0) q(iby,1:jm) =   buni
        else if (xl >= xcut) then
          if (iby > 0) q(iby,1:jm) = - buni
        else
          if (iby > 0) q(iby,1:jm) = 0.0d+00
        end if

! convert the primitive variables to conservative ones
!
        call prim2cons(jm, q(1:nv,1:jm), u(1:nv,1:jm))

! copy the conserved variables to the current block
!
        pdata%u(1:nv,i,1:jm,k) = u(1:nv,1:jm)

! copy the primitive variables to the current block
!
        pdata%q(1:nv,i,1:jm,k) = q(1:nv,1:jm)

      end do ! i = 1, im
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_current_sheet
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_JET:
! ----------------------------
!
!   Subroutine sets the initial conditions for the relativistic jet
!   problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_jet(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real(kind=8), save :: djet  = 1.00d-01
    real(kind=8), save :: damb  = 1.00d+01
    real(kind=8), save :: bamb  = 1.00d-08
    real(kind=8), save :: pres  = 1.00d-02
    real(kind=8), save :: vjet  = 0.99d+00
    real(kind=8), save :: bjet  = 1.00d-03
    real(kind=8), save :: ljet  = 1.00d-00
    real(kind=8), save :: rjet  = 1.00d+00
    real(kind=8), save :: rjet2 = 1.00d+00

! local saved parameters
!
    logical     , save :: first = .true.

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: dx, dy, dz, rm, rr

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(nv)    :: qj
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("djet"  , djet)
      call get_parameter_real("damb"  , damb)
      call get_parameter_real("pres"  , pres)
      call get_parameter_real("bamb"  , bamb)
      call get_parameter_real("bjet"  , bjet)
      call get_parameter_real("vjet"  , vjet)
      call get_parameter_real("ljet"  , ljet)
      call get_parameter_real("rjet"  , rjet)

! calculate Rjet²
!
      rjet2 = rjet * rjet

! reset the first execution flag
!
      first = .false.

    end if ! first call

! set the conditions inside the jet radius
!
    qj(idn) = djet
    if (ipr > 0) qj(ipr) = pres
    qj(ivx) = vjet
    qj(ivy) = 0.0d+00
    qj(ivz) = 0.0d+00
    if (ibx > 0) then
      qj(ibx) = 0.0d+00
      qj(iby) = 0.0d+00
      qj(ibz) = bjet
      qj(ibp) = 0.0d+00
    end if ! ibx > 0

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    dx = x(2) - x(1)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
    dy = y(2) - y(1)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
    dz = z(2) - z(1)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
    dz      = 0.0d+00
#endif /* NDIMS == 3 */
    rm = dy * dy + dz * dz

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do j = 1, jm

! calculate radius
!
        rr = y(j) * y(j) + z(k) * z(k)

! set the ambient density, pressure, and velocity
!
        q(idn,1:im) = damb
        if (ipr > 0) q(ipr,1:im) = pres
        q(ivx,1:im) = 0.0d+00
        q(ivy,1:im) = 0.0d+00
        q(ivz,1:im) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
        if (ibx > 0) then
          q(ibx,1:im) = 0.0d+00
          q(iby,1:im) = 0.0d+00
          q(ibz,1:im) = bamb
          q(ibp,1:im) = 0.0d+00
        end if ! ibx > 0

! set the jet injection
!
        if (rr <= max(rm, rjet2)) then
          do i = 1, im
            if (x(i) <= max(dx, ljet)) then
              q(1:nv,i) = qj(1:nv)
            end if
          end do ! i = 1, im
        end if ! R < Rjet

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
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_jet
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_BINARIES:
! ---------------------------------
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
  subroutine setup_problem_binaries(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use constants  , only : d2r, pi2
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az, adx, ady, adz, advol
    use equations  , only : prim2cons
    use equations  , only : gamma
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real, get_parameter_integer

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

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

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: dcomp, pstar, pcomp
    real(kind=8), save :: r2star, r2comp
    real(kind=8), save :: acomp , bcomp
    real(kind=8), save :: xps, xpc, uys, uyc
    real(kind=8), save :: om

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

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the problems setup
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_binaries

!===============================================================================
!
end module problems
