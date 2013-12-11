!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az, adx, ady, adz
    use equations  , only : prim2cons
    use equations  , only : gamma
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use parameters , only : get_parameter_real
    use equations  , only : nv

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real   , save :: dens = 1.0d0, ratio = 1.0e2, radius = 0.1d0
    real   , save :: csnd = 0.40824829046386301635d0
    logical, save :: first = .true.
    real   , save :: dn_amb, dn_ovr
    real   , save :: pr_amb, pr_ovr
    real   , save :: rad

! local variables
!
    integer       :: i, j, k
    real          :: xl, yl, zl, xu, yu, zu, rl, ru
    real          :: xb, yb, xt, yt
    real          :: dx, dy, dz, dxh, dyh, dzh, daxy
    real          :: fc_amb, fc_ovr

! local arrays
!
    real, dimension(nv,im) :: q, u
    real, dimension(im)    :: x
    real, dimension(jm)    :: y
    real, dimension(km)    :: z
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
      rad    = radius * radius

! reset the first execution flag
!
      first = .false.

    end if

! obtain block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d0
#endif /* NDIMS == 3 */

! calculate mesh intervals and areas
!
    dx   = adx(pdata%meta%level)
    dy   = ady(pdata%meta%level)
    dz   = adz(pdata%meta%level)
    dxh  = 0.5d0 * dx
    dyh  = 0.5d0 * dy
#if NDIMS == 3
    dzh  = 0.5d0 * dz
#else /* NDIMS == 3 */
    dzh  = 1.0d0
#endif /* NDIMS == 3 */
    daxy = dx * dy

! set the uniform primitive variables
!
    q(ivx,:) = 0.0d0
    q(ivy,:) = 0.0d0
    q(ivz,:) = 0.0d0

! iterate over all positions in the YZ plane
!
    do k = 1, km

#ifdef R3D
! calculate the corner Z coordinates
!
      zl = abs(z(k)) - dzh
      zu = abs(z(k)) + dzh
#endif /* R3D */

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
#ifdef R3D
          rl = xl * xl + yl * yl + zl * zl
          ru = xu * xu + yu * yu + zu * zu
#else /* R3D */
          rl = xl * xl + yl * yl
          ru = xu * xu + yu * yu
#endif /* R3D */

! set the initial density and pressure
!
          q(idn,i) = dn_amb
          if (ipr > 0) q(ipr,i) = pr_amb

! set the initial pressure in cells laying completely within the radius
!
          if (ru .le. rad) then

! set the overpressure region density
!
            q(idn,i) = dn_ovr

! set the overpressure region pressure
!
            if (ipr > 0) q(ipr,i) = pr_ovr

! set the initial pressure in the cell completely outside the radius
!
          else if (rl .ge. rad) then

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

#ifdef R3D
! in 3D simply set the ambient values since the integration is more complex
!

! set the ambient region density
!
            q(idn,i) = dn_amb

! set the ambient medium pressure
!
            if (ipr > 0) q(ipr,i) = pr_amb
#else /* R3D */
! calculate the bounds of area integration
!
            xb = max(xl, sqrt(max(0.0d0, rad - yu * yu)))
            xt = min(xu, sqrt(max(0.0d0, rad - yl * yl)))
            yb = max(yl, sqrt(max(0.0d0, rad - xu * xu)))
            yt = min(yu, sqrt(max(0.0d0, rad - xl * xl)))

! integrate the area below the circle within the current cell for both
! functions, y = f(x) and x = g(y), and then average them to be sure that we
! are getting the ideal symmetry
!
            fc_ovr = 0.5d0 * (rad * (asin(xt / radius) - asin(xb / radius))   &
                   + (xt * yb - xb * yt)) - yl * (xt - xb)
            fc_ovr = fc_ovr + (xb - xl) * dy

            fc_amb = 0.5d0 * (rad * (asin(yt / radius) - asin(yb / radius))   &
                   + (yt * xb - yb * xt)) - xl * (yt - yb)
            fc_amb = fc_amb + (yb - yl) * dx

            fc_ovr = 0.5d0 * (fc_ovr + fc_amb)

! normalize coefficients
!
            fc_ovr = fc_ovr / daxy
            fc_amb = 1.0d0 - fc_ovr

! integrate the density over the edge cells
!
            q(idn,i) = fc_ovr * dn_ovr + fc_amb * dn_amb

! integrate the pressure over the edge cells
!
            if (ipr > 0) q(ipr,i) = fc_ovr * pr_ovr + fc_amb * pr_amb
#endif /* R3D */

          end if

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

      end do
    end do

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_blast

!===============================================================================
!
end module problems
