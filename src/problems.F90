!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2012 Grzegorz Kowal <grzegorz@amuncode.org>
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
    use coordinates, only : ax, ay, az
    use equations  , only : prim2cons
    use equations  , only : gamma
    use parameters , only : get_parameter_real
    use variables  , only : nt
    use variables  , only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables  , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables  , only : ibx, iby, ibz
#ifdef GLM
    use variables  , only : iph
#endif /* GLM */
#endif /* MHD */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real   , save :: dens = 1.0d0, ratio = 1.0e2, radius = 0.1d0
#ifdef ADI
    real   , save :: csnd = 0.40824829046386301635d0
    real   , save :: pres = 1.0d0
#endif /* ADI */
    logical, save :: first = .true.

! local variables
!
    integer                :: i, j, k
    real                   :: r

! local arrays
!
    real, dimension(nt,im) :: q, u
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
#ifdef ADI
      call get_parameter_real("csnd"  , csnd  )

! calculate pressure
!
      pres = dens * csnd * csnd / gamma
#endif /* ADI */

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

! set the uniform variables
!
    q(idn,1:im) = dens
    q(ivx,1:im) = 0.0d0
    q(ivy,1:im) = 0.0d0
    q(ivz,1:im) = 0.0d0
#ifdef ADI
    q(ipr,1:im) = pres
#endif /* ADI */
#ifdef MHD
    q(ibx,1:im) = 1.0d0 / sqrt(2.0d0)
    q(iby,1:im) = 1.0d0 / sqrt(2.0d0)
    q(ibz,1:im) = 0.0d0
#ifdef GLM
    q(iph,1:im) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial star profile (density for ISO or pressure for ADI)
!
    do k = 1, km
      do j = 1, jm

        do i = 1, im

! calculate distance from the coordinate system origin
!
          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

! fill in the internal radius
!
          if (r < radius) then
#ifdef ADI
            q(ipr,i) = pres * ratio
#endif /* ADI */
#ifdef ISO
            q(idn,i) = dens * ratio
#endif /* ISO */
          else
#ifdef ADI
            q(ipr,i) = pres
#endif /* ADI */
#ifdef ISO
            q(idn,i) = dens
#endif /* ISO */
          end if

        end do

! convert the primitive variables to conserved ones
!
        call prim2cons(im, q(:,:), u(:,:))

! copy the conserved variables to the current block
!
        pdata%u(1:nt,1:im,j,k) = u(1:nt,1:im)

! copy the primitive variables to the current block
!
        pdata%q(1:nt,1:im,j,k) = q(1:nt,1:im)

      end do
    end do

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_blast

!===============================================================================
!
end module problems
