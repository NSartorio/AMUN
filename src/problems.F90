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
!! module: PROBLEMS
!!
!!  This module handles the initialization of various test and research
!!  problems.
!!
!!
!!******************************************************************************
!
module problems

! module variables are not implicit by default
!
  implicit none

! module variable to store the problem name
!
  character(len=32), save :: problem = "blast"

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_problems, setup_problem

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
! subroutine INITIALIZE_PROBLEMS:
! ------------------------------
!
!   Subroutine prepares module PROBLEMS.
!
!
!===============================================================================
!
  subroutine initialize_problems()

! include external procedures and variables
!
    use parameters , only : get_parameter_string

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the problem name
!
    call get_parameter_string("problem", problem)

!-------------------------------------------------------------------------------
!
  end subroutine initialize_problems
!
!===============================================================================
!
! subroutine SETUP_PROBLEM:
! ------------------------
!
!   Subroutine sets the initial conditions for selected problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!
!===============================================================================
!
  subroutine setup_problem(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use error      , only : print_error

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! select the setup subroutine depending on the problem name
!
    select case(problem)

! general test problems
!
    case("blast")
      call setup_problem_blast(pdata)

! specific problems
!
    case("reconnection")
      call setup_problem_reconnection(pdata)

    case default
      call print_error("problems::init_problem()"                              &
                     , "Setup subroutime is not implemented for this problem!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
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
    use coordinates, only : ax, ay, az, adx, ady, adz
    use equations  , only : prim2cons
    use equations  , only : gamma
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
    real(kind=8), save :: dens   = 1.00d+00
    real(kind=8), save :: ratio  = 1.00d+02
    real(kind=8), save :: radius = 1.00d-01
    real(kind=8), save :: csnd   = 4.0824829046386301635d-01
    real(kind=8), save :: buni   = 1.00d+00
    real(kind=8), save :: angle  = 4.50d+01

! local saved parameters
!
    logical     , save :: first = .true.
    real(kind=8), save :: dn_amb, dn_ovr
    real(kind=8), save :: pr_amb, pr_ovr
    real(kind=8), save :: r2

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: xl, yl, zl, xu, yu, zu, rl, ru
    real(kind=8)  :: xb, yb, xt, yt
    real(kind=8)  :: dx, dy, dz, dxh, dyh, dzh, daxy
    real(kind=8)  :: fc_amb, fc_ovr

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
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

      q(ibx,:) = buni * cos(d2r * angle)
      q(iby,:) = buni * sin(d2r * angle)
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
! in 3D simply set the ambient values since the integration is more complex
!

! set the ambient region density
!
            q(idn,i) = dn_amb

! set the ambient medium pressure
!
            if (ipr > 0) q(ipr,i) = pr_amb
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

! integrate the density over the edge cells
!
            q(idn,i) = fc_ovr * dn_ovr + fc_amb * dn_amb

! integrate the pressure over the edge cells
!
            if (ipr > 0) q(ipr,i) = fc_ovr * pr_ovr + fc_amb * pr_amb
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

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_blast
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_RECONNECTION:
! -------------------------------------
!
!   Subroutine sets the initial conditions for the reconnection problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_reconnection(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km
    use coordinates, only : ay
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
    real(kind=8), save :: dens = 1.00d+00
    real(kind=8), save :: pres = 1.00d+00
    real(kind=8), save :: bamp = 1.00d+00
    real(kind=8), save :: bgui = 0.00d+00
    real(kind=8), save :: vper = 1.00d-02
    real(kind=8), save :: ycut = 1.00d-01

! local saved parameters
!
    logical     , save :: first = .true.

! local variables
!
    integer            :: i, j, k

! local arrays
!
    real(kind=8), dimension(nv,jm) :: q, u
    real(kind=8), dimension(jm)    :: y
!
!-------------------------------------------------------------------------------
!
! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("dens"  , dens)
      call get_parameter_real("pres"  , pres)
      call get_parameter_real("brecon", bamp)
      call get_parameter_real("bguide", bgui)
      call get_parameter_real("vper"  , vper)
      call get_parameter_real("ycut"  , ycut)

! reset the first execution flag
!
      first = .false.

    end if ! first call

! prepare block coordinates
!
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)

! set the density and pressure
!
    q(idn,:) = dens
    if (ipr > 0) q(ipr,:) = pres

! reset velocity components
!
    q(ivx,:) = 0.0d+00
    q(ivy,:) = 0.0d+00
    q(ivz,:) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
    if (ibx > 0) then

      q(ibx,:) = 0.0d+00
      q(iby,:) = 0.0d+00
      q(ibz,:) = bgui
      q(ibp,:) = 0.0d+00

    end if

! iterate over all positions in the XZ plane
!
    do k = 1, km
      do i = 1, im

! sweep along the Y coordinate
!
        do j = 1, jm

! set the random velocity field in a layer near current sheet
!
          if (abs(y(j)) <= ycut) then

            q(ivx,j) = vper * randomn()
            q(ivy,j) = vper * randomn()
#if NDIMS == 3
            q(ivz,j) = vper * randomn()
#endif /* NDIMS == 3 */

          end if

! set antiparallel magnetic field component
!
          if (ibx > 0) q(ibx,j) = sign(bamp, y(j))

        end do ! j = 1, jm

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

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_reconnection

!===============================================================================
!
end module problems
