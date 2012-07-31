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
!! module: EVOLUTION
!!
!!  This module performs the time integration using different methods.
!!
!!******************************************************************************
!
module evolution

! module variables are not implicit by default
!
  implicit none

! evolution parameters
!
  real   , save :: cfl   = 0.5d0

#if defined MHD && defined GLM
! coefficient controlling the decay of scalar potential Psi
!
  real   , save :: alpha_p = 0.5d0
  real   , save :: decay   = 1.0d0
#endif /* MHD & GLM */

  integer, save :: n
  real   , save :: t, dt, dtn, dxmin

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_evolution, advance

! declare public variables
!
  public :: n, cfl, t, dt, dtn, dxmin

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
!   Subroutine initializes the EVOLUTION module by setting its parameters.
!
!
!===============================================================================
!
  subroutine initialize_evolution()

! include external procedures and variables
!
    use parameters, only : get_parameter_real

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the value of the stability coefficient
!
    call get_parameter_real("cfl", cfl)

#if defined MHD && defined GLM
! get the value of the stability coefficient
!
    call get_parameter_real("alpha_p", alpha_p)

! calculate decay factor for magnetic field divergence scalar
!
    decay = exp(- alpha_p * cfl)
#endif /* MHD & GLM */

! calculate the initial time step
!
    call find_new_timestep()

!-------------------------------------------------------------------------------
!
  end subroutine initialize_evolution
!
!===============================================================================
!
! subroutine ADVANCE:
! ------------------
!
!   Subroutine advances the solution by one time step using the selected
!   time integration method.
!
!
!===============================================================================
!
  subroutine advance()

! include external procedures and variables
!
    use blocks     , only : block_data, list_data
    use boundaries , only : boundary_variables
#ifdef REFINE
    use coordinates, only : toplev
#endif /* REFINE */
    use mesh       , only : update_mesh

! local variables are not implicit by default
!
    implicit none

! local variables
!
    type(block_data), pointer :: pblock
!
!-------------------------------------------------------------------------------
!
#ifdef RK2
! advance the solution using the 2nd order Runge-Kutta method
    call advance_rk2()
#endif /* RK2 */

#ifdef REFINE
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
#endif /* REFINE */

! update primitive variables
!
    call update_variables()

! find new time step
!
    call find_new_timestep()

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

! include external procedures and variables
!
    use blocks     , only : block_data, list_data
    use boundaries , only : boundary_variables
    use coordinates, only : im, jm, km
    use coordinates, only : adx, ady, adz
    use variables  , only : nfl

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_data), pointer    :: pblock

! local arrays
!
    real, dimension(3)            :: dh
    real, dimension(nfl,im,jm,km) :: du
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
      call update_interval(dh(:), pblock%f(:,:,:,:,:), du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u1(1:nfl,:,:,:) = pblock%u0(1:nfl,:,:,:) + du(1:nfl,:,:,:)

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
      call update_interval(dh(:), pblock%f(:,:,:,:,:), du(:,:,:,:))

! update the solution for the fluid variables
!
      pblock%u0(1:nfl,:,:,:) = 0.5d0 * (pblock%u0(1:nfl,:,:,:)                 &
                                   + pblock%u1(1:nfl,:,:,:) + du(1:nfl,:,:,:))

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

! include external procedures and variables
!
    use blocks     , only : block_data, list_data
    use boundaries , only : boundary_correct_fluxes
    use coordinates, only : adx, ady, adz
    use scheme     , only : update_flux

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
    call boundary_correct_fluxes()

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

! include external procedures and variables
!
    use blocks     , only : block_data, list_data
    use equations  , only : update_primitive_variables

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
! subroutine UPDATE_INTERVAL:
! --------------------------
!
!   Subroutine calculate the conservative variable interval from fluxes.
!
!
!===============================================================================
!
  subroutine update_interval(dh, f, du)

! include external procedures and variables
!
    use coordinates, only : im, jm, km
    use variables  , only : nfl

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real, dimension(3)                 , intent(in)    :: dh
    real, dimension(NDIMS,nfl,im,jm,km), intent(in)    :: f
    real, dimension(      nfl,im,jm,km), intent(inout) :: du

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! reset the increment array du
!
    du(:,:,:,:) = 0.0d0

! perform update along the X direction
!
    do i = 2, im
      du(:,i,:,:) = du(:,i,:,:) - dh(1) * (f(1,:,i,:,:) - f(1,:,i-1,:,:))
    end do
    du(:,1,:,:) = du(:,1,:,:) - dh(1) * f(1,:,1,:,:)

! perform update along the Y direction
!
    do j = 2, jm
      du(:,:,j,:) = du(:,:,j,:) - dh(2) * (f(2,:,:,j,:) - f(2,:,:,j-1,:))
    end do
    du(:,:,1,:) = du(:,:,1,:) - dh(2) * f(2,:,:,1,:)

#if NDIMS == 3
! perform update along the Z direction
!
    do k = 2, km
      du(:,:,:,k) = du(:,:,:,k) - dh(3) * (f(3,:,:,:,k) - f(3,:,:,:,k-1))
    end do
    du(:,:,:,1) = du(:,:,:,1) - dh(3) * f(3,:,:,:,1)
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine update_interval
!
!===============================================================================
!
! find_new_timestep: subroutine updates the maximum speed among the leafs and
!                    calculates new time step
!
!===============================================================================
!
  subroutine find_new_timestep()

    use blocks     , only : block_meta, block_data, list_meta, list_data
    use coordinates, only : toplev
    use coordinates, only : adx, ady, adz
#ifdef MPI
    use mpitools   , only : reduce_maximum_real
#endif /* MPI */
    use scheme  , only : maxspeed, cmax

    implicit none

! local variables
!
    integer                   :: iret
    real                      :: cm
    integer(kind=4)           :: lev

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed, and highest level
!
    cmax   = 1.0d-16
    lev    = 1

! if toplev > 1, find the highest level
!
    if (toplev .gt. 1) then

! iterate over all meta blocks and find the highest level with leafs
!
      pmeta => list_meta
      do while (associated(pmeta))

! check if the metablock is a leaf, if so obtaind the highest level
!
        if (pmeta%leaf) lev = max(lev, pmeta%level)

! associate the pointer with the next block
!
        pmeta => pmeta%next

      end do ! meta blocks

    end if ! toplev > 1

! find the smallest spacial step
!
#if NDIMS == 2
    dxmin = min(adx(lev), ady(lev))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    dxmin = min(adx(lev), ady(lev), adz(lev))
#endif /* NDIMS == 3 */

! iterate over all data blocks in order to find the maximum speed among them
!
    pdata => list_data
    do while (associated(pdata))

! check if this block is a leaf
!
      if (pdata%meta%leaf) then

! find the maximum level occupied by blocks (can be smaller than toplev)
!


! obtain the maximum speed for the current block
!
        cm = maxspeed(pdata%q(:,:,:,:))

! compare global and local maximum speeds
!
        cmax = max(cmax, cm)

      end if ! leaf

! assiociate the pointer with the next block
!
      pdata => pdata%next

    end do

#ifdef MPI
! find maximum speed in the system from all processors
!
    call reduce_maximum_real(cmax, iret)
#endif /* MPI */

! calcilate new time step
!
    dtn = dxmin / max(cmax, 1.0d-16)

!-------------------------------------------------------------------------------
!
  end subroutine find_new_timestep

!===============================================================================
!
end module
