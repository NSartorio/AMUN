!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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

#ifdef PROFILE
! import external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer     , save :: imi, ima, imt, imu, imf, iui, imv
#endif /* PROFILE */

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
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('evolution:: initialization'  , imi)
    call set_timer('evolution:: solution advance', ima)
    call set_timer('evolution:: new time step'   , imt)
    call set_timer('evolution:: solution update' , imu)
    call set_timer('evolution:: flux update'     , imf)
    call set_timer('evolution:: increment update', iui)
    call set_timer('evolution:: variable update' , imv)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

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

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

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
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! nullify pointer to integration subroutine
!
    nullify(evolve)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

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
#ifdef PROFILE
! start accounting time for solution advance
!
    call start_timer(ima)
#endif /* PROFILE */

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
      call update_variables(time + dt, 0.0d+00)

! set all meta blocks to be updated
!
      call set_blocks_update(.true.)

    end if ! toplev > 1

#ifdef DEBUG
! check variables for NaNs
!
    call check_variables()
#endif /* DEBUG */

#ifdef PROFILE
! stop accounting time for solution advance
!
    call stop_timer(ima)
#endif /* PROFILE */

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
#ifdef PROFILE
! start accounting time for new time step estimation
!
    call start_timer(imt)
#endif /* PROFILE */

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

#ifdef PROFILE
! stop accounting time for new time step estimation
!
    call stop_timer(imt)
#endif /* PROFILE */

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
    real(kind=8) :: tm, dtm

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

! prepare times
!
    tm  = time + dt
    dtm = dt

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)            &
                                     + dt * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    real(kind=8) :: tm, dtm

! local arrays
!
    real(kind=8), dimension(nv,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

!= 1st step: U(1) = U(n) + dt * F[U(n)]
!
! prepare times
!
    tm  = time + dt
    dtm = dt

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u0(1:nv,1:im,1:jm,1:km)            &
                                     + dt * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u1

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

!= 2nd step: U(n+1) = 1/2 U(n) + 1/2 U(1) + 1/2 dt * F[U(1)]
!
! prepare times
!
    tm  = time + dt
    dtm = 0.5d+00 * dt

! update fluxes from the intermediate stage
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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = 0.5d+00 * (pdata%u0(1:nv,1:im,1:jm,1:km) &
                                               + pdata%u1(1:nv,1:im,1:jm,1:km) &
                                               +  dt * du(1:nv,1:im,1:jm,1:km))

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    real(kind=8) :: tm, dtm

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
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

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

! prepare times
!
      tm  = time + n * ds
      dtm = ds

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
        call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
      call update_variables(tm, dtm)

    end do ! n = 1, stages - 1

!= the final step: U(n+1) = 1/m U(0) + (m-1)/m [1 + dt/(m-1) L] U(m-1)
!
! prepare times
!
    tm  = time + dt
    dtm = fl * dt

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = fl *  pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + fr * (pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km))

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pointer to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    real(kind=8) :: ds
    real(kind=8) :: tm, dtm

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
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

!= 1st substep: U(1) = U(n) + dt F[U(n)]
!
! prepare the fractional time step
!
    ds = dt

! prepare times
!
    tm  = time + ds
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the first intermediate solution
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
    call update_variables(tm, dtm)

!= 2nd step: U(2) = 3/4 U(n) + 1/4 U(1) + 1/4 dt F[U(1)]
!
! prepare the fractional time step
!
    ds = f22 * dt

! prepare times
!
    tm  = time + 0.5d+00 * dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the second intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = f21 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + f22 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

!= 3rd step: U(n+1) = 1/3 U(n) + 2/3 U(2) + 2/3 dt F[U(2)]
!
! prepare the fractional time step
!
    ds = f32 * dt

! prepare times
!
    tm  = time + dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = f31 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                    + f32 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    real(kind=8) :: ds
    real(kind=8) :: tm, dtm

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
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

!= 1st step: U(1) = U(n) + 1/2 dt F[U(n)]
!
! calculate the fractional time step
!
    ds = b1 * dt

! prepare times
!
    tm  = time + ds
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the intermediate solution
!
      pdata%u1(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! prepare times
!
    tm  = time + dt
    dtm = ds

! update primitive variables
!
    call update_variables(tm, dtm)

!= 3rd step: U(3) = 2/3 U(n) + 1/3 U(2) + 1/6 dt F[U(2)]
!
! calculate the fractional time step
!
    ds = b3 * dt

! prepare times
!
    tm  = time + 0.5d+00 * dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

!= the final step: U(n+1) = U(3) + 1/2 dt F[U(3)]
!
! calculate fractional time step
!
    ds = b1 * dt

! prepare times
!
    tm  = time + dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = pdata%u1(1:nv,1:im,1:jm,1:km)            &
                                     + ds * du(1:nv,1:im,1:jm,1:km)

! update the conservative variable pointer
!
      pdata%u => pdata%u0

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    real(kind=8) :: ds
    real(kind=8) :: tm, dtm

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
    real(kind=8), parameter :: a33 = 6.44090224936673d-01
    real(kind=8), parameter :: a41 = 3.67933791638137d-01
    real(kind=8), parameter :: a44 = 6.32066208361863d-01
    real(kind=8), parameter :: a53 = 2.37593836598569d-01
    real(kind=8), parameter :: a55 = 7.62406163401431d-01
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for one step update
!
    call start_timer(imu)
#endif /* PROFILE */

!= 1st step: U(1) = U(n) + b1 dt F[U(n)]
!
! calculate the fractional time step
!
    ds = b1 * dt

! prepare times
!
    tm  = time + ds
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

!= 2nd step: U(2) = U(1) + b1 dt F[U(1)]
!
! prepare times
!
    tm  = time + 2.0d+00 * ds
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

!= 3rd step: U(3) = a31 U(n) + a33 U(2) + b3 dt F[U(2)]
!
! calculate the fractional time step
!
    ds = b3 * dt

! prepare times
!
    tm  = time + (2.0d+00 * a33 * b1 + b3) * dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

!= 4th step: U(4) = a41 U(n) + a44 U(3) + b4 dt F[U(3)]
!
! calculate the fractional time step
!
    ds = b4 * dt

! prepare times
!
    tm  = time + ((2.0d+00 * b1 * a33 + b3) * a44 + b4) * dt
    dtm = ds

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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

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
    call update_variables(tm, dtm)

!= the final step: U(n+1) = a53 U(2) + a55 U(4) + b5 dt F[U(4)]
!
! calculate the fractional time step
!
    ds = b5 * dt

! prepare times
!
    tm  = time + dt
    dtm = ds


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
      call update_sources(pdata, tm, dtm, du(1:nv,1:im,1:jm,1:km))

! update the final solution
!
      pdata%u0(1:nv,1:im,1:jm,1:km) = a53 * pdata%u1(1:nv,1:im,1:jm,1:km)      &
                                    + a55 * pdata%u0(1:nv,1:im,1:jm,1:km)      &
                                           + ds * du(1:nv,1:im,1:jm,1:km)

! update ψ with its source term
!
      if (ibp > 0) pdata%u(ibp,1:im,1:jm,1:km) =                               &
                                           decay * pdata%u(ibp,1:im,1:jm,1:km)

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! update primitive variables
!
    call update_variables(tm, dtm)

#ifdef PROFILE
! stop accounting time for one step update
!
    call stop_timer(imu)
#endif /* PROFILE */

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
    use coordinates   , only : im, jm, km
    use equations     , only : nv

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer  :: pdata

! local vectors
!
    real(kind=8), dimension(NDIMS) :: dx
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for fluxe update
!
    call start_timer(imf)
#endif /* PROFILE */

! assign pdata with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks
!
    do while (associated(pdata))

! obtain dx, dy, and dz for the current block
!
      dx(1) = adx(pdata%meta%level)
      dx(2) = ady(pdata%meta%level)
#if NDIMS == 3
      dx(3) = adz(pdata%meta%level)
#endif /* NDIMS == 3 */

! update fluxes for the current block
!
      call update_flux(dx(1:NDIMS), pdata%q(1:nv,1:im,1:jm,1:km)               &
                                       , pdata%f(1:NDIMS,1:nv,1:im,1:jm,1:km))

! assign pdata to the next block
!
      pdata => pdata%next

    end do ! over data blocks

! correct the numerical fluxes of the blocks which have neighbours at higher
! levels
!
    call boundary_fluxes()

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_fluxes
!
!===============================================================================
!
! subroutine UPDATE_INCREMENT:
! ---------------------------
!
!   Subroutine calculates the conservative variable increment from
!   directional fluxes.
!
!   Arguments:
!
!     pdata - the point to data block storing the directional fluxes;
!     du    - the array of the conservative variable increment;
!
!===============================================================================
!
  subroutine update_increment(pdata, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im , jm , km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use coordinates    , only : adxi, adyi, adzi
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(  out) :: du

! local variables
!
    integer      :: i  , j  , k
    integer      :: im1, jm1, km1
    real(kind=8) :: dxi, dyi, dzi
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the increment update
!
    call start_timer(iui)
#endif /* PROFILE */

! prepare the coordinate intervals
!
    dxi = adxi(pdata%meta%level)
    dyi = adyi(pdata%meta%level)
#if NDIMS == 3
    dzi = adzi(pdata%meta%level)
#endif /* NDIMS == 3 */

! initialize du
!
    du(1:nv,1:im,1:jm,1:km) = 0.0d+00

! calculate the variable update from the directional fluxes
!
    do k = kbl, keu
#if NDIMS == 3
      km1 = k - 1
#endif /* NDIMS == 3 */
      do j = jbl, jeu
        jm1 = j - 1
        do i = ibl, ieu
          im1 = i - 1

#if NDIMS == 3
          du(1:nv,i,j,k) =                                                     &
                     - dxi * (pdata%f(1,1:nv,i,j,k) - pdata%f(1,1:nv,im1,j,k)) &
                     - dyi * (pdata%f(2,1:nv,i,j,k) - pdata%f(2,1:nv,i,jm1,k)) &
                     - dzi * (pdata%f(3,1:nv,i,j,k) - pdata%f(3,1:nv,i,j,km1))
#else /* NDIMS == 3 */
          du(1:nv,i,j,k) = &
                     - dxi * (pdata%f(1,1:nv,i,j,k) - pdata%f(1,1:nv,im1,j,k)) &
                     - dyi * (pdata%f(2,1:nv,i,j,k) - pdata%f(2,1:nv,i,jm1,k))
#endif /* NDIMS == 3 */
        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
    end do ! k = kbl, keu

#ifdef PROFILE
! stop accounting time for the increment update
!
    call stop_timer(iui)
#endif /* PROFILE */

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
!   Arguments:
!
!     tm  - time at the moment of update;
!     dtm - time step since the last update;
!
!===============================================================================
!
  subroutine update_variables(tm, dtm)

! include external procedures
!
    use blocks        , only : set_neighbors_update
    use boundaries    , only : boundary_variables
    use equations     , only : update_primitive_variables
    use equations     , only : fix_unphysical_cells, correct_unphysical_states
    use shapes        , only : update_shapes

! include external variables
!
    use blocks        , only : block_meta, list_meta
    use blocks        , only : block_data, list_data

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(in) :: tm, dtm

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable update
!
    call start_timer(imv)
#endif /* PROFILE */

! update primitive variables in the changed blocks
!
    pdata => list_data
    do while (associated(pdata))
      pmeta => pdata%meta

      if (pmeta%update) call update_primitive_variables(pdata%u, pdata%q)

      pdata => pdata%next
    end do

! update boundaries
!
    call boundary_variables(tm, dtm)

! correct unphysical states if detected
!
    if (fix_unphysical_cells) then

! if an unphysical cell appeared in a block while updating its primitive
! variables it could be propagated to its neighbors through boundary update;
! mark all neighbors of such a block to be verified and corrected for
! unphysical cells too
!
      pmeta => list_meta
      do while (associated(pmeta))

        if (pmeta%update) call set_neighbors_update(pmeta)

        pmeta => pmeta%next
      end do

! verify and correct, if necessary, unphysical cells in recently updated blocks
!
      pdata => list_data
      do while (associated(pdata))
        pmeta => pdata%meta

        if (pmeta%update)                                                      &
              call correct_unphysical_states(step, pmeta%id, pdata%q, pdata%u)

        pdata => pdata%next
      end do
    end if

! apply shapes in blocks which need it
!
    pdata => list_data
    do while (associated(pdata))
      pmeta => pdata%meta

      if (pmeta%update) call update_shapes(pdata, tm, dtm)

      pdata => pdata%next
    end do

#ifdef PROFILE
! stop accounting time for variable update
!
    call stop_timer(imv)
#endif /* PROFILE */

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
    use coordinates    , only : im, jm, km
    use equations      , only : nv, pvars, cvars
    use ieee_arithmetic, only : ieee_is_nan

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
              if (ieee_is_nan(pdata%u(p,i,j,k))) then
                print *, 'U NaN:', cvars(p), pdata%meta%id, i, j, k
              end if
              if (ieee_is_nan(pdata%q(p,i,j,k))) then
                print *, 'Q NaN:', pvars(p), pdata%meta%id, i, j, k
              end if
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
