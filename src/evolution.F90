!!******************************************************************************
!!
!! module: evolution - handling the time evolution of the block structure
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of Godunov-AMR.
!!
!!  Godunov-AMR is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Godunov-AMR is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!
module evolution

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! evolve: subroutine sweeps over all leaf blocks and performs one step time
!         evolution for each according to the selected integration scheme
!
!===============================================================================
!
  subroutine evolve()

    use blocks      , only : block_data, list_data
    use boundaries  , only : boundary_variables
    use mesh        , only : update_mesh
    use mesh        , only : dx_min
    use scheme      , only : cmax
    use timer       , only : start_timer, stop_timer

    implicit none

! local variables
!
    type(block_data), pointer :: pblock
    real                      :: cm
!
!-------------------------------------------------------------------------------
!
! iterate over all data blocks and perform one step of time evolution
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a leaf
!
      if (pblock%meta%leaf) &
#ifdef EULER
        call evolve_euler(pblock)
#endif /* EULER */
#ifdef RK2
        call evolve_rk2(pblock)
#endif /* RK2 */
#ifdef RK3
        call evolve_rk3(pblock)
#endif /* RK3 */

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update boundaries
!
    call start_timer(4)
    call boundary_variables()
    call stop_timer(4)

! check refinement and refine
!
    call start_timer(5)
    call update_mesh()
    call stop_timer(5)

! update boundaries
!
    call start_timer(4)
    call boundary_variables()
    call stop_timer(4)

! update the maximum speed
!
    call update_maximum_speed()

! get maximum time step
!
    dtn = dx_min / max(cmax, 1.0d-16)
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve
!
!===============================================================================
!
! update_maximum_speed: subroutine updates module variable cmax with the value
!                       corresponding to the maximum speed in the system
!
!===============================================================================
!
  subroutine update_maximum_speed()

    use blocks  , only : block_data, list_data
#ifdef MPI
    use mpitools, only : mallreducemaxr
#endif /* MPI */
    use scheme  , only : maxspeed, cmax

    implicit none

! local variables
!
    type(block_data), pointer :: pblock
    real                      :: cm
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed
!
    cmax = 1.0d-16

! iterate over all blocks in order to find the maximum speed
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a leaf
!
      if (pblock%meta%leaf) &
        cm = maxspeed(pblock%u)

! compare global and local maximum speeds
!
      cmax = max(cmax, cm)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

#ifdef MPI
! reduce the maximum speed over all processes
!
    call mallreducemaxr(cmax)
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
  end subroutine update_maximum_speed
#ifdef EULER
!
!===============================================================================
!
! evolve_euler: subroutine evolves the current block using Euler integration
!
!===============================================================================
!
  subroutine evolve_euler(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : alpha_p
    use mesh     , only : dx_min
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
#if defined MHD && defined GLM
    real    :: ch2, decay
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! 1st step of integration
!
    call update(pblock%meta%level, pblock%u(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    pblock%u(1:nfl,:,:,:) = pblock%u(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:)

#ifdef MHD
! update the solution for the magnetic variables
!
    pblock%u(ibx:ibz,:,:,:) = pblock%u(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:)

#ifdef GLM
! calculate c_h^2
!
    ch2 = cmax * cmax

! update the solution for the scalar potential Psi
!
    pblock%u(iph,:,:,:) = pblock%u(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)

! evolve Psi due to the source term
!
    decay = exp(- alpha_p * cmax * dt / dx_min)
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_euler
#endif /* EULER */
#ifdef RK2
!
!===============================================================================
!
! evolve_rk2: subroutine evolves the current block using the 2nd order
!             Runge-Kutta method
!
!===============================================================================
!
  subroutine evolve_rk2(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : alpha_p
    use mesh     , only : dx_min
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
#if defined MHD && defined GLM
    real    :: ch2, decay
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du
!
!-------------------------------------------------------------------------------
!
#if defined MHD && defined GLM
! calculate c_h^2
!
    ch2 = cmax * cmax

#endif /* MHD & GLM */
!! 1st step of integration
!!
    call update(pblock%meta%level, pblock%u(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    u1(1:nfl,:,:,:) = pblock%u(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:)

#ifdef MHD
! update the solution for the magnetic variables
!
    u1(ibx:ibz,:,:,:) = pblock%u(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:)

#ifdef GLM
! update the solution for the scalar potential Psi
!
    u1(iph,:,:,:) = pblock%u(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

! 2nd step of integration
!
    call update(pblock%meta%level, u1(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    pblock%u(1:nfl,:,:,:) = 0.5d0 * (pblock%u(1:nfl,:,:,:)                     &
                                    + u1(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:))

#ifdef MHD
! update the solution for the magnetic variables
!
    pblock%u(ibx:ibz,:,:,:) = 0.5d0 * (pblock%u(ibx:ibz,:,:,:)                 &
                                + u1(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:))

#ifdef GLM
! update the solution for the scalar potential Psi
!
    pblock%u(iph,:,:,:) = 0.5d0 * (pblock%u(iph,:,:,:)                         &
                                  + u1(iph,:,:,:) + ch2 * dt * du(iph,:,:,:))

! evolve Psi due to the source term
!
    decay = exp(- alpha_p * cmax * dt / dx_min)
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk2
#endif /* RK2 */
#ifdef RK3
!
!===============================================================================
!
! evolve_rk3: subroutine evolves the current block using the 3rd order
!             Runge-Kutta method
!
!===============================================================================
!
  subroutine evolve_rk3(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : alpha_p
    use mesh     , only : dx_min
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
#if defined MHD && defined GLM
    real    :: ch2, decay
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du

! parameters
!
    real, parameter :: f4 = 1.0d0 / 4.0d0, f3 = 1.0d0 / 3.0d0
!
!-------------------------------------------------------------------------------
!
#if defined MHD && defined GLM
! calculate c_h^2
!
    ch2 = cmax * cmax

#endif /* MHD & GLM */
!! 1st step of integration
!!
    call update(pblock%meta%level, pblock%u(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    u1(1:nfl,:,:,:) = pblock%u(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:)

#ifdef MHD
! update the solution for the magnetic variables
!
    u1(ibx:ibz,:,:,:) = pblock%u(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:)

#ifdef GLM
! update the solution for the scalar potential Psi
!
    u1(iph,:,:,:) = pblock%u(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

!! 2nd step of integration
!!
    call update(pblock%meta%level, u1(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    u1(1:nfl,:,:,:) = f4 * (3.0d0 * pblock%u(1:nfl,:,:,:)                      &
                                    + u1(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:))

#ifdef MHD
! update the solution for the magnetic variables
!
    u1(ibx:ibz,:,:,:) = f4 * (3.0d0 * pblock%u(ibx:ibz,:,:,:)                  &
                                + u1(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:))

#ifdef GLM
! update the solution for the scalar potential Psi
!
    u1(iph,:,:,:) = f4 * (3.0d0 * pblock%u(iph,:,:,:)                          &
                                  + u1(iph,:,:,:) + ch2 * dt * du(iph,:,:,:))
#endif /* GLM */
#endif /* MHD */

!! 3rd step of integration
!!
    call update(pblock%meta%level, u1(:,:,:,:), du(:,:,:,:))

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du(:,:,:,:))
#endif /* SHAPE */

! update the solution for the fluid variables
!
    pblock%u(1:nfl,:,:,:) = f3 * (pblock%u(1:nfl,:,:,:)                        &
                          + 2.0d0 * (u1(1:nfl,:,:,:) + dt * du(1:nfl,:,:,:)))

#ifdef MHD
! update the solution for the magnetic variables
!
    pblock%u(ibx:ibz,:,:,:) = f3 * (pblock%u(ibx:ibz,:,:,:)                    &
                      + 2.0d0 * (u1(ibx:ibz,:,:,:) + dt * du(ibx:ibz,:,:,:)))

#ifdef GLM
! update the solution for the scalar potential Psi
!
    pblock%u(iph,:,:,:) = f3 * (pblock%u(iph,:,:,:)                            &
                        + 2.0d0 * (u1(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)))

! evolve analytically Psi due to the source term
!
    decay = exp(- alpha_p * cmax * dt / dx_min)
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk3
#endif /* RK3 */
!
!===============================================================================
!
end module
