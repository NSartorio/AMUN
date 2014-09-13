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
!! module: EVOLUTION
!!
!!  This module provides an interface for temporal integration with
!!  the stability handling.  New integration methods can be added by
!!  implementing more evolve_* subroutines.
!!
!!******************************************************************************
!
module evolution

! module variables are not implicit by default
!
  implicit none

! pointer to the temporal integration subroutine
!
  procedure(evolve_euler), pointer, save :: evolve => null()

! evolution parameters
!
  real(kind=8), save :: cfl     = 5.0d-01
  integer     , save :: stages  = 2

! coefficient controlling the decay of scalar potential ѱ
!
  real(kind=8), save :: alpha   = 2.0d+00
  real(kind=8), save :: decay   = 1.0d+00

! time variables
!
  integer     , save :: step    = 0
  real(kind=8), save :: time    = 0.0d+00
  real(kind=8), save :: dt      = 1.0d+00
  real(kind=8), save :: dtn     = 1.0d+00

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_evolution, finalize_evolution
  public :: advance, new_time_step

! declare public variables
!
  public :: cfl, step, time, dt, dtn

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
! subroutine INITIALIZE_EVOLUTION:
! -------------------------------
!
!   Subroutine initializes module EVOLUTION by setting its parameters.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_evolution(verbose, iret)

! include external procedures
!
    use parameters, only : get_parameter_string, get_parameter_real            &
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
    character(len=255) :: integration = "rk2"
    character(len=255) :: name_int    = ""
!
!-------------------------------------------------------------------------------
!
! get the integration method and the value of the CFL coefficient
!
    call get_parameter_string ("time_advance", integration)
    call get_parameter_integer("stages"      , stages     )
    call get_parameter_real   ("cfl"         , cfl        )
    call get_parameter_real   ("alpha"       , alpha      )

! select the integration method, check the correctness of the integration
! parameters and adjust the CFL coefficient if necessary
!
    select case(trim(integration))

    case ("euler", "EULER")

      name_int =  "1st order Euler"
      evolve   => evolve_euler

    case ("rk2", "RK2")

      name_int =  "2nd order Runge-Kutta"
      evolve   => evolve_rk2

    case ("ssprk(m,2)", "SSPRK(m,2)")

      stages   = max(2, min(9, stages))
      write(name_int, "('2nd order SSPRK(',i1,',2)')") stages
      evolve   => evolve_ssprk2
      cfl      = (stages - 1) * cfl

    case ("rk3", "RK3")

      name_int =  "3rd order Runge-Kutta"
      evolve   => evolve_rk3

    case ("rk3.4", "RK3.4", "ssprk(4,3)", "SSPRK(4,3)")

      name_int =  "3rd order SSPRK(4,3)"
      evolve   => evolve_ssprk34
      cfl      = 2.0d+00 * cfl

    case ("rk3.5", "RK3.5", "ssprk(5,3)", "SSPRK(5,3)")

      name_int =  "3rd order SSPRK(5,3)"
      evolve   => evolve_ssprk35
      cfl      = 2.65062919143939d+00 * cfl

    case default

      if (verbose) then
        write (*,"(1x,a)") "The selected time advance method is not " //       &
                           "implemented: " // trim(integration)
        name_int =  "2nd order Runge-Kutta method"
        evolve   => evolve_rk2
      end if

    end select

! calculate the decay factor for magnetic field divergence scalar source term
!
    decay = exp(- alpha * cfl)

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "time advance           =", trim(name_int)

    end if

!-------------------------------------------------------------------------------
!
  end subroutine initialize_evolution
!
!===============================================================================
!
! subroutine FINALIZE_EVOLUTION:
! -----------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_evolution(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
    nullify(evolve)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_evolution
!
!===============================================================================
!
! subroutine ADVANCE:
! ------------------
!
!   Subroutine advances the solution by one time step using the selected time
!   integration method.
!
!   Arguments:
!
!     dtnext - next time step;
!
!===============================================================================
!
  subroutine advance(dtnext)

! include external procedures
!
    use blocks        , only : set_blocks_update
    use boundaries    , only : boundary_variables
    use mesh          , only : update_mesh

! include external variables
!
    use coordinates   , only : toplev

! local variables are not implicit by default
!
    implicit none

! input variables
!
    real(kind=8), intent(in) :: dtnext
!
!-------------------------------------------------------------------------------
!
! find new time step
!
    call new_time_step(dtnext)

! advance the solution using the selected method
!
    call evolve()

! chec if we need to perform the refinement step
!
    if (toplev > 1) then

! set all meta blocks to not be updated
!
      call set_blocks_update(.false.)

! check refinement and refine
!
      call update_mesh()

! update primitive variables
!
      call update_variables()

! update boundaries
!
      call boundary_variables()

#ifdef DEBUG
! check variables for NaNs
!
      call check_variables()
#endif /* DEBUG */

! set all meta blocks to be updated
!
      call set_blocks_update(.true.)

    end if ! toplev > 1

!-------------------------------------------------------------------------------
!
  end subroutine advance
!
!===============================================================================
!
! subroutine NEW_TIME_STEP:
! ------------------------
!
!   Subroutine estimates the new time step from the maximum speed in the system
!   and source term constraints.
!
!   Arguments:
!
!     dtnext - next time step;
!
!===============================================================================
!
  subroutine new_time_step(dtnext)

! include external procedures
!
    use equations     , only : maxspeed, cmax, cmax2
#ifdef MPI
    use mpitools      , only : reduce_maximum_real, reduce_maximum_integer
#endif /* MPI */

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : adx, ady, adz
    use coordinates   , only : toplev
    use sources       , only : viscosity, resistivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(in)  :: dtnext

! local pointers
!
    type(block_data), pointer :: pdata

! local variables
!
    integer                   :: iret, mlev
    real(kind=8)              :: cm, dx_min

! local parameters
!
    real(kind=8), parameter   :: eps = tiny(cmax)
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed, and the highest level
!
    cmax   = eps
    mlev   = 1

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks in order to find the maximum speed among them
! and the highest level which is required to obtain the minimum spacial step
!
    do while (associated(pdata))

! update the maximum level
!
      mlev = max(mlev, pdata%meta%level)

! get the maximum speed for the current block
!
      cm = maxspeed(pdata%q(:,:,:,:))

! update the maximum speed
!
      cmax = max(cmax, cm)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

#ifdef MPI
! reduce maximum speed and level over all processors
!
    call reduce_maximum_real   (cmax, iret)
    call reduce_maximum_integer(mlev, iret)
#endif /* MPI */

! calculate the squared cmax
!
    cmax2 = cmax * cmax

! find the smallest spacial step
!
#if NDIMS == 2
    dx_min = min(adx(mlev), ady(mlev))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    dx_min = min(adx(mlev), ady(mlev), adz(mlev))
#endif /* NDIMS == 3 */

! calculate the new time step
!
    dtn = cfl * dx_min / max(cmax                                              &
                        + 2.0d+00 * max(viscosity, resistivity) / dx_min, eps)

! round the time
!
    if (dtnext > 0.0d+00) dt = min(dtn, dtnext)

!-------------------------------------------------------------------------------
!
  end subroutine new_time_step
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine EVOLVE_EULER:
! -----------------------
!
!   Subroutine advances the solution by one time step using the 1st order
!   Euler integration method.
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!===============================================================================
!
  subroutine evolve_euler()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer            :: pdata

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)            &
                                     + dt * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,:,:,:) = decay * pdata%u(ibp,:,:,:)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_euler
!
!===============================================================================
!
! subroutine EVOLVE_RK2:
! ---------------------
!
!   Subroutine advances the solution by one time step using the 2nd order
!   Runge-Kutta time integration method.
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!===============================================================================
!
  subroutine evolve_rk2()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer    :: pblock

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! update fluxes for the first step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! calculate variable increment for the current block
!
      call update_increment(pblock, du(:,:,:,:))

! add source terms
!
      call update_sources(pblock, du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u1(1:nv,:,:,:) = pblock%u0(1:nv,:,:,:) + dt * du(1:nv,:,:,:)

! update the conservative variable pointer
!
      pblock%u => pblock%u1

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

! update fluxes for the second step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! calculate variable increment for the current block
!
      call update_increment(pblock, du(:,:,:,:))

! add source terms
!
      call update_sources(pblock, du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u0(1:nv,:,:,:) = 0.5d+00 * (pblock%u0(1:nv,:,:,:)                 &
                                + pblock%u1(1:nv,:,:,:) + dt * du(1:nv,:,:,:))

! update the conservative variable pointer
!
      pblock%u => pblock%u0

! update ψ by its source term
!
      if (ibp > 0) pblock%u(ibp,:,:,:) = decay * pblock%u(ibp,:,:,:)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk2
!
!===============================================================================
!
! subroutine EVOLVE_SSPRK2:
! ------------------------
!
!   Subroutine advances the solution by one time step using the 2nd order
!   m-stage Strong Stability Preserving Runge-Kutta time integration method.
!   Up to 9 stages are allowed, due to stability problems with more stages.
!
!   References:
!
!     [1] Gottlieb, S. and Gottlieb, L.-A., J.
!         "Strong Stability Preserving Properties of Runge-Kutta Time
!          Discretization Methods for Linear Constant Coefficient Operators",
!         Journal of Scientific Computing,
!         2003, vol. 18, no. 1, pp. 83-109
!
!===============================================================================
!
  subroutine evolve_ssprk2()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer :: pdata

! local variables
!
    integer      :: n
    real(kind=8) :: ds

! local saved variables
!
    logical     , save :: first = .true.
    real(kind=8), save :: ft, fl, fr

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! prepare things which don't change later
!
    if (first) then

! calculate integration coefficients
!
      ft = 1.0d+00 / (stages - 1)
      fl = 1.0d+00 / stages
      fr = 1.0d+00 - fl

! update first flag
!
      first = .false.

    end if

! calculate the fractional time step
!
    ds = ft * dt

!= 1st step: U(0) = U(n)
!
! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! copy conservative array u0 to u1
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u1

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

!= 2nd step: U(i) = [1 + dt/(m-1) L] U(i-1), for i = 1, ..., m-1
!
! integrate the intermediate steps
!
    do n = 1, stages - 1

! update fluxes
!
      call update_fluxes()

! assign pdata with the first block on the data block list
!
      pdata => list_data

! iterate over all data blocks
!
      do while (associated(pdata))

! calculate the variable increment
!
        call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
        call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
        pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)          &
                                       + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
        pdata => pdata%next

      end do ! over data blocks

! update primitive variables
!
      call update_variables()

! update boundaries
!
      call boundary_variables()

    end do ! n = 1, stages - 1

!= the final step: U(n+1) = 1/m U(0) + (m-1)/m [1 + dt/(m-1) L] U(m-1)
!
! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = fl *  pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + fr * (pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km))

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ by its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pointer to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_ssprk2
!
!===============================================================================
!
! subroutine EVOLVE_RK3:
! ---------------------
!
!   Subroutine advances the solution by one time step using the 3rd order
!   Runge-Kutta time integration method.
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!
!===============================================================================
!
  subroutine evolve_rk3()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer    :: pblock

! local variables
!
    real(kind=8)                         :: ds

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du

! local integration parameters
!
    real(kind=8), parameter :: f21 = 3.0d+00 / 4.0d+00, f22 = 1.0d+00 / 4.0d+00
    real(kind=8), parameter :: f31 = 1.0d+00 / 3.0d+00, f32 = 2.0d+00 / 3.0d+00
!
!-------------------------------------------------------------------------------
!
!! 1st substep of integration
!!
! prepare fractional time step
!
    ds = dt

! update fluxes for the first step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! calculate variable increment for the current block
!
      call update_increment(pblock, du(:,:,:,:))

! add source terms
!
      call update_sources(pblock, du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u1(1:nv,:,:,:) = pblock%u0(1:nv,:,:,:) + ds * du(1:nv,:,:,:)

! update the conservative variable pointer
!
      pblock%u => pblock%u1

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!! 2nd substep of integration
!!
! prepare fractional time step
!
    ds = f22 * dt

! update fluxes for the first step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! calculate variable increment for the current block
!
      call update_increment(pblock, du(:,:,:,:))

! add source terms
!
      call update_sources(pblock, du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u1(1:nv,:,:,:) = f21 * pblock%u0(1:nv,:,:,:)                      &
                            + f22 * pblock%u1(1:nv,:,:,:) + ds * du(1:nv,:,:,:)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!! 3rd substep of integration
!!
! prepare fractional time step
!
    ds = f32 * dt

! update fluxes for the second step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! calculate variable increment for the current block
!
      call update_increment(pblock, du(:,:,:,:))

! add source terms
!
      call update_sources(pblock, du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u0(1:nv,:,:,:) = f31 * pblock%u0(1:nv,:,:,:)                      &
                            + f32 * pblock%u1(1:nv,:,:,:) + ds * du(1:nv,:,:,:)

! update the conservative variable pointer
!
      pblock%u => pblock%u0

! update ψ by its source term
!
      if (ibp > 0) pblock%u(ibp,:,:,:) = decay * pblock%u(ibp,:,:,:)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk3
!
!===============================================================================
!
! subroutine EVOLVE_SSPRK34:
! -------------------------
!
!   Subroutine advances the solution by one time step using the 3rd order
!   4-stage Strong Stability Preserving Runge-Kutta time integration method.
!
!   References:
!
!     [1] Ruuth, S. J.,
!         "Global Optimization of Explicit Strong-Stability-Preserving
!          Runge-Kutta methods",
!         Mathematics of Computation,
!         2006, vol. 75, no. 253, pp. 183-207
!
!===============================================================================
!
  subroutine evolve_ssprk34()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer :: pdata

! local variables
!
    real(kind=8)                         :: ds

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du

! local integration parameters
!
    real(kind=8), parameter :: b1  = 1.0d+00 / 2.0d+00, b3  = 1.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 = 2.0d+00 / 3.0d+00, a33 = 1.0d+00 / 3.0d+00
!
!-------------------------------------------------------------------------------
!
!= 1st step: U(1) = U(n) + 1/2 dt F[U(n)]
!
! calculate the fractional time step
!
    ds = b1 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u1

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= 2nd step: U(2) = U(1) + 1/2 dt F[U(1)]
!
! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= 3rd step: U(3) = 2/3 U(n) + 1/3 U(2) + 1/6 dt F[U(2)]
!
! calculate the fractional time step
!
    ds = b3 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = a31 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + a33 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= the final step: U(n+1) = U(3) + 1/2 dt F[U(3)]
!
! calculate fractional time step
!
    ds = b1 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ by its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_ssprk34
!
!===============================================================================
!
! subroutine EVOLVE_SSPRK35:
! -------------------------
!
!   Subroutine advances the solution by one time step using the 3rd order
!   5-stage Strong Stability Preserving Runge-Kutta time integration method.
!
!   References:
!
!     [1] Ruuth, S. J.,
!         "Global Optimization of Explicit Strong-Stability-Preserving
!          Runge-Kutta methods",
!         Mathematics of Computation,
!         2006, vol. 75, no. 253, pp. 183-207
!
!===============================================================================
!
  subroutine evolve_ssprk35()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use sources       , only : update_sources

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im, jm, km
    use equations     , only : nv, ibp

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer :: pdata

! local variables
!
    real(kind=8)                         :: ds

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du

! local integration parameters
!
    real(kind=8), parameter :: b1  = 3.77268915331368d-01
    real(kind=8), parameter :: b3  = 2.42995220537396d-01
    real(kind=8), parameter :: b4  = 2.38458932846290d-01
    real(kind=8), parameter :: b5  = 2.87632146308408d-01
    real(kind=8), parameter :: a31 = 3.55909775063327d-01
    real(kind=8), parameter :: a33 = 6.44090224936674d-01
    real(kind=8), parameter :: a41 = 3.67933791638137d-01
    real(kind=8), parameter :: a44 = 6.32066208361863d-01
    real(kind=8), parameter :: a53 = 2.37593836598569d-01
    real(kind=8), parameter :: a55 = 7.62406163401431d-01
!
!-------------------------------------------------------------------------------
!
!= 1st step: U(1) = U(n) + b1 dt F[U(n)]
!
! calculate the fractional time step
!
    ds = b1 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u1

! assign pdata to the next block
!
      pdata => pdata%next

    end do

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= 2nd step: U(2) = U(1) + b1 dt F[U(1)]
!
! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= 3rd step: U(3) = a31 U(n) + a33 U(2) + b3 dt F[U(2)]
!
! calculate the fractional time step
!
    ds = b3 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = a31 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + a33 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= 4th step: U(4) = a41 U(n) + a44 U(3) + b4 dt F[U(3)]
!
! calculate the fractional time step
!
    ds = b4 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = a41 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + a44 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!= the final step: U(n+1) = a53 U(2) + a55 U(4) + b5 dt F[U(4)]
!
! calculate the fractional time step
!
    ds = b5 * dt

! update fluxes
!
    call update_fluxes()

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! calculate the variable increment
!
      call update_increment(pdata, du(1:nv,1:im,1:jm,1:km))

! add the source terms
!
      call update_sources(pdata, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = a53 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                    + a55 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! update ψ by its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables()

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine evolve_ssprk35
!
!===============================================================================
!
! subroutine UPDATE_FLUXES:
! ------------------------
!
!   Subroutine iterates over all data blocks and calculates the numerical
!   fluxes for each block.  After the fluxes are updated, they are corrected
!   for blocks which have neighbours at higher refinement level.
!
!
!===============================================================================
!
  subroutine update_fluxes()

! include external procedures
!
    use boundaries    , only : boundary_fluxes
    use schemes       , only : update_flux

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : adx, ady, adz

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer  :: pblock

! local vectors
!
    real(kind=8), dimension(3) :: dx

! local variables
!
    integer                    :: n
!
!-------------------------------------------------------------------------------
!
! iterate over all data blocks
!
    pblock => list_data
    do while (associated(pblock))

! obtain dx, dy, and dz for the current block
!
      dx(1) = adx(pblock%meta%level)
      dx(2) = ady(pblock%meta%level)
      dx(3) = adz(pblock%meta%level)

! update the flux for the current block
!
      do n = 1, NDIMS
        call update_flux(n, dx(n), pblock%q(:,:,:,:), pblock%f(n,:,:,:,:))
      end do

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! correct the numerical fluxes of the blocks which have neighbours at higher
! level
!
    call boundary_fluxes()

!-------------------------------------------------------------------------------
!
  end subroutine update_fluxes
!
!===============================================================================
!
! subroutine UPDATE_INCREMENT:
! ---------------------------
!
!   Subroutine calculate the conservative variable increment from the fluxes.
!
!   Arguments:
!
!     dh   - the ratio of the time step to the spatial step;
!     f    - the array of numerical fluxes;
!     du   - the array of variable increment;
!
!===============================================================================
!
  subroutine update_increment(pdata, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use coordinates    , only : adxi, adyi, adzi
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: du

! local variables
!
    integer      :: i  , j  , k
    real(kind=8) :: dxi, dyi, dzi
!
!-------------------------------------------------------------------------------
!
! reset the increment array du
!
    du(:,:,:,:) = 0.0d+00

! prepare coordinate intervals
!
    dxi = adxi(pdata%meta%level)
    dyi = adyi(pdata%meta%level)
#if NDIMS == 3
    dzi = adzi(pdata%meta%level)
#endif /* NDIMS == 3 */

! perform update along the X direction
!
    do i = ibl, ieu
      du(:,i,:,:) = du(:,i,:,:)                                                &
                           - dxi * (pdata%f(1,:,i,:,:) - pdata%f(1,:,i-1,:,:))
    end do

! perform update along the Y direction
!
    do j = jbl, jeu
      du(:,:,j,:) = du(:,:,j,:)                                                &
                           - dyi * (pdata%f(2,:,:,j,:) - pdata%f(2,:,:,j-1,:))
    end do

#if NDIMS == 3
! perform update along the Z direction
!
    do k = kbl, keu
      du(:,:,:,k) = du(:,:,:,k)                                                &
                           - dzi * (pdata%f(3,:,:,:,k) - pdata%f(3,:,:,:,k-1))
    end do
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine update_increment
!
!===============================================================================
!
! subroutine UPDATE_VARIABLES:
! ---------------------------
!
!   Subroutine iterates over all data blocks and converts the conservative
!   variables to their primitive representation.
!
!
!===============================================================================
!
  subroutine update_variables()

! include external procedures
!
    use equations     , only : update_primitive_variables
    use shapes        , only : update_shapes

! include external variables
!
    use blocks        , only : block_meta, list_meta
    use blocks        , only : block_data, list_data

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! associate the pointer with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! associate pmeta with the corresponding meta block
!
      pmeta => pdata%meta

! convert conserved variables to primitive ones for the current block and
! update shapes if necessary
!
      if (pmeta%update) then
        call update_primitive_variables(pdata%u, pdata%q)
        call update_shapes(pdata)
      end if

! assign pointer to the next block
!
      pdata => pdata%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine update_variables
#ifdef DEBUG
!
!===============================================================================
!
! subroutine CHECK_VARIABLES:
! --------------------------
!
!   Subroutine iterates over all data blocks and converts the conservative
!   variables to their primitive representation.
!
!
!===============================================================================
!
  subroutine check_variables()

! include external procedures
!
    use coordinates   , only : im, jm, km
    use equations     , only : nv, pvars, cvars
#ifdef IBM
    use, intrinsic :: ieee_arithmetic
#endif /* IBM */

! include external variables
!
    use blocks        , only : block_meta, list_meta
    use blocks        , only : block_data, list_data

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: i, j, k, p

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! associate the pointer with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! associate pmeta with the corresponding meta block
!
      pmeta => pdata%meta

! check if there are NaNs in primitive variables
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im
            do p = 1, nv
#ifdef IBM
              if (ieee_is_nan(pdata%u(p,i,j,k))) then
                print *, 'U NaN:', cvars(p), pdata%meta%id, i, j, k
              end if
              if (ieee_is_nan(pdata%q(p,i,j,k))) then
                print *, 'Q NaN:', pvars(p), pdata%meta%id, i, j, k
              end if
#else /* IBM */
              if (isnan(pdata%u(p,i,j,k))) then
                print *, 'U NaN:', cvars(p), pdata%meta%id, i, j, k
              end if
              if (isnan(pdata%q(p,i,j,k))) then
                print *, 'Q NaN:', pvars(p), pdata%meta%id, i, j, k
              end if
#endif /* IBM */
            end do ! p = 1, nv
          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

! assign pointer to the next block
!
      pdata => pdata%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine check_variables
#endif /* DEBUG */

!===============================================================================
!
end module evolution
