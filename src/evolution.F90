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
!! module: EVOLUTION
!!
!!  This module constains subroutines to perform the time integration using
!!  different methods preserving the conservation and estimates the new time
!!  step from the stability condition.
!!
!!
!!******************************************************************************
!
module evolution

! module variables are not implicit by default
!
  implicit none

! evolution parameters
!
  real   , save :: cfl     = 0.5d+0

! time variables
!
  integer, save :: n       = 0
  real   , save :: t       = 0.0d+0
  real   , save :: dt      = 1.0d+0
  real   , save :: dtn     = 1.0d+0

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_evolution, advance

! declare public variables
!
  public :: cfl, n, t, dt, dtn

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
! subroutine INITIALIZE_EVOLUTION:
! -------------------------------
!
!   Subroutine initializes module EVOLUTION by setting its parameters.
!
!
!===============================================================================
!
  subroutine initialize_evolution()

! include external procedures
!
    use parameters    , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the value of the stability coefficient
!
    call get_parameter_real("cfl", cfl)

! calculate the initial time step
!
    call new_time_step()

!-------------------------------------------------------------------------------
!
  end subroutine initialize_evolution
!
!===============================================================================
!
! subroutine ADVANCE:
! ------------------
!
!   Subroutine advances the solution by one time step using the selected time
!   integration method.
!
!
!===============================================================================
!
  subroutine advance()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use mesh          , only : update_mesh

! include external variables
!
    use coordinates   , only : toplev

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
#ifdef RK2
! advance the solution using the 2nd order Runge-Kutta method
    call advance_rk2()
#endif /* RK2 */

! chec if we need to perform the refinement step
!
    if (toplev > 1) then

! check refinement and refine
!
      call update_mesh()

! update boundaries
!
      call boundary_variables()

    end if ! toplev > 1

! update primitive variables
!
    call update_variables()

! find new time step
!
    call new_time_step()

!-------------------------------------------------------------------------------
!
  end subroutine advance
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
! subroutine ADVANCE_RK2:
! ----------------------
!
!   Subroutine advances the solution by one time step using the 2nd order
!   Runge-Kutta time integration method.
!
!
!===============================================================================
!
  subroutine advance_rk2()

! include external procedures
!
    use boundaries    , only : boundary_variables
    use schemes       , only : update_increment

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
    type(block_data), pointer    :: pblock

! local arrays
!
    real, dimension(3)           :: dh
    real, dimension(nv,im,jm,km) :: du
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

! obtain dx, dy, and dz for the current block
!
      dh(1) = dt / adx(pblock%meta%level)
      dh(2) = dt / ady(pblock%meta%level)
      dh(3) = dt / adz(pblock%meta%level)

! calculate variable increment for the current block
!
      call update_increment(dh(:), pblock%f(:,:,:,:,:), du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u1(1:nv,:,:,:) = pblock%u0(1:nv,:,:,:) + du(1:nv,:,:,:)

! update the conservative variable pointer
!
      pblock%u => pblock%u1

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update boundaries
!
    call boundary_variables()

! update primitive variables
!
    call update_variables()

! update fluxes for the second step of the RK2 integration
!
    call update_fluxes()

! update the solution using numerical fluxes stored in the data blocks
!
    pblock => list_data
    do while (associated(pblock))

! obtain dx, dy, and dz for the current block
!
      dh(1) = dt / adx(pblock%meta%level)
      dh(2) = dt / ady(pblock%meta%level)
      dh(3) = dt / adz(pblock%meta%level)

! calculate variable increment for the current block
!
      call update_increment(dh(:), pblock%f(:,:,:,:,:), du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u0(1:nv,:,:,:) = 0.5d0 * (pblock%u0(1:nv,:,:,:)                 &
                                   + pblock%u1(1:nv,:,:,:) + du(1:nv,:,:,:))

! update the conservative variable pointer
!
      pblock%u => pblock%u0

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update boundaries
!
    call boundary_variables()

!-------------------------------------------------------------------------------
!
  end subroutine advance_rk2
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
    type(block_data), pointer :: pblock

! local vectors
!
    real, dimension(3)        :: dx

! local variables
!
    integer                   :: n
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

! include external variables
!
    use blocks        , only : block_data, list_data

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer :: pblock
!
!-------------------------------------------------------------------------------
!
! iterate over all data blocks
!
    pblock => list_data
    do while (associated(pblock))

! convert conserved variables to primitive ones for the current block
!
      call update_primitive_variables(pblock%u, pblock%q)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine update_variables
!
!===============================================================================
!
! subroutine NEW_TIME_STEP:
! ------------------------
!
!   Subroutine estimates the new time step from the maximum speed in the system
!   and source term constraints.
!
!
!===============================================================================
!
  subroutine new_time_step()

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

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer :: pblock

! local variables
!
    integer                   :: iret
    integer(kind=4)           :: lev
    real                      :: cm, dx_min

! local parameters
!
    real, parameter           :: eps = tiny(cmax)
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed, and the highest level
!
    cmax   = eps
    lev    = 1

! iterate over all data blocks in order to find the maximum speed among them
! and the highest level which is required to obtain the spatial step
!
    pblock => list_data
    do while (associated(pblock))

! find the maximum level occupied by blocks (can be smaller than toplev)
!
      lev = max(lev, pblock%meta%level)

! obtain the maximum speed for the current block
!
      cm = maxspeed(pblock%q(:,:,:,:))

! compare global and local maximum speeds
!
      cmax = max(cmax, cm)

! assiociate the pointer with the next block
!
      pblock => pblock%next

    end do

#ifdef MPI
! find maximum speed in the system from all processors
!
    call reduce_maximum_real   (cmax, iret)
    call reduce_maximum_integer(lev , iret)
#endif /* MPI */

! calculate squared cmax
!
    cmax2 = cmax * cmax

! find the smallest spatial step
!
#if NDIMS == 2
    dx_min = min(adx(lev), ady(lev))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    dx_min = min(adx(lev), ady(lev), adz(lev))
#endif /* NDIMS == 3 */

! calcilate the new time step
!
    dtn = dx_min / max(cmax, eps)

!-------------------------------------------------------------------------------
!
  end subroutine new_time_step

!===============================================================================
!
end module evolution
