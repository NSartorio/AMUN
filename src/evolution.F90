!!******************************************************************************
!!
!! module: evolution - handling the time evolution of the block structure
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of AMUN.
!!
!!  This program is free software; you can redistribute it and/or
!!  modify it under the terms of the GNU General Public License
!!  as published by the Free Software Foundation; either version 2
!!  of the License, or (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program; if not, write to the Free Software Foundation,
!!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
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

    use blocks    , only : block_data, list_data
    use boundaries, only : boundary_variables
#ifdef CONSERVATIVE
    use boundaries, only : boundary_correct_fluxes
#endif /* CONSERVATIVE */
#ifdef REFINE
    use config    , only : maxlev
#endif /* REFINE */
#ifdef VISCOSITY
    use config    , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config    , only : ueta
#endif /* MHD & RESISTIVITY */
#ifdef FORCE
    use forcing   , only : fourier_transform, evolve_forcing
#endif /* FORCE */
    use mesh      , only : update_mesh
    use mesh      , only : dx_min
    use scheme    , only : cmax
    use timer     , only : start_timer, stop_timer
    use variables , only : idn, imz

    implicit none

! local variables
!
    type(block_data), pointer :: pblock
    real                      :: cm
!
!-------------------------------------------------------------------------------
!
#ifdef FORCE
! perform the Fourier transform of the velocity field
!
    pblock => list_data
    do while (associated(pblock))

      if (pblock%meta%leaf)                                                    &
        call fourier_transform(pblock%meta%level                               &
                       , pblock%meta%xmin, pblock%meta%ymin, pblock%meta%zmin  &
                       , pblock%u(idn:imz,:,:,:))

      pblock => pblock%next ! assign pointer to the next block

    end do

! evolve the forcing source terms by the time interval dt
!
    call evolve_forcing(dt)
#endif /* FORCE */

! iterate over all data blocks and perform one step of time evolution
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a leaf
!
#ifdef CONSERVATIVE
      if (pblock%meta%leaf) &
#ifdef EULER
        call flux_euler(pblock)
#endif /* EULER */
#ifdef RK2
        call flux_rk2(pblock)
#endif /* RK2 */
#ifdef RK3
        call flux_rk3(pblock)
#endif /* RK3 */
#else /* CONSERVATIVE */
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
#endif /* CONSERVATIVE */

! assign pointer to the next block
!
      pblock => pblock%next

    end do

#ifdef CONSERVATIVE
! correct the numerical fluxes between neighboring blocks which are at different
! levels
!
    call start_timer(4)
    if (maxlev .gt. 1) call boundary_correct_fluxes()
    call stop_timer(4)

! update solution using numerical fluxes stored in data blocks
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a  leaf and update its conserved variables using
! corrected numerical fluxes
!
      if (pblock%meta%leaf) call update_solution(pblock)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

#endif /* CONSERVATIVE */
! update boundaries
!
    call start_timer(4)
    call boundary_variables()
    call stop_timer(4)

#ifdef REFINE
! chec if we need to perform the refinement step
!
    if (maxlev .gt. 1) then

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

    end if ! maxlev > 1

#endif /* REFINE */
! update the maximum speed
!
    call update_maximum_speed()

! get maximum time step
!
    dtn = dx_min / max(cmax, 1.0d-16)
#ifdef VISCOSITY
    dtn = min(dtn, 0.5d0 * dx_min * dx_min / max(1.0d-16, visc))
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    dtn = min(dtn, 0.5d0 * dx_min * dx_min / max(1.0d-16, ueta))
#endif /* MHD & RESISTIVITY */
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
#ifdef CONSERVATIVE
!
!===============================================================================
!
! update_solution: subroutine performs an one step update of the conserved
!                  variables for the given data block using the integrated
!                  numerical fluxes stored in the same data block
!
!===============================================================================
!
  subroutine update_solution(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
    use mesh     , only : adxi, adyi, adzi
#if defined MHD && defined GLM
    use config   , only : alpha_p
    use mesh     , only : dx_min
    use scheme   , only : cmax
    use variables, only : iph
#endif /* MHD & GLM */
#ifdef FORCE
    use forcing  , only : real_forcing
    use variables, only : inx, iny, inz
    use variables, only : idn, imx, imy, imz
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi
#if defined MHD && defined GLM
    real    :: decay
#endif /* MHD & GLM */
#ifdef FORCE

! local arrays
!
    real, dimension(  3,im,jm,km) :: f
#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
#if NDIMS == 3
    dzi = adzi(pblock%meta%level)
#endif /* NDIMS == 3 */

! perform update of conserved variables of the given block using its fluxes
!
    call advance_solution(pblock%u(:,:,:,:), pblock%f(:,:,:,:,:), dxi, dyi, dzi)
#if defined MHD && defined GLM

! evolve Psi due to the source term
!
    decay = exp(- alpha_p * cmax * dt / dx_min)
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* MHD & GLM */
#ifdef FORCE

! obtain the forcing terms in real space
!
    call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin    &
                                       , pblock%meta%zmin, f(:,:,:,:))

! update momenta due to the forcing terms
!
    pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:)                                  &
                                          + pblock%u(idn,:,:,:) * f(inx,:,:,:)
    pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:)                                  &
                                          + pblock%u(idn,:,:,:) * f(iny,:,:,:)
    pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:)                                  &
                                          + pblock%u(idn,:,:,:) * f(inz,:,:,:)
#endif /* FORCE */

!-------------------------------------------------------------------------------
!
  end subroutine update_solution
!
!===============================================================================
!
! advance_solution: subroutine performs an one step update of the conserved
!                   variables array using the numerical fluxes passed as an
!                   argument
!
!===============================================================================
!
  subroutine advance_solution(u, f, dxi, dyi, dzi)

    use config   , only : im, jm, km
    use variables, only : nqt, nfl
    use variables, only : inx, iny, inz
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use scheme   , only : cmax
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

    implicit none

! input arguments
!
    real, dimension(      nqt,im,jm,km), intent(inout) :: u
    real, dimension(NDIMS,nqt,im,jm,km), intent(in)    :: f
    real                                               :: dxi, dyi, dzi

! local variables
!
    integer :: i, j, k, im1, jm1, km1
    real    :: dhx, dhy, dhz
#if defined MHD && defined GLM
    real    :: ch2
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dhx = dt * dxi
    dhy = dt * dyi
#if NDIMS == 3
    dhz = dt * dzi
#endif /* NDIMS == 3 */

! reset the increment array du
!
    du(:,:,:,:) = 0.0d0

! perform update along the X direction
!
    do i = 1, im
      im1 = max(1, i - 1)

      du(:,i,:,:) = du(:,i,:,:) + dhx * (f(inx,:,im1,:,:) - f(inx,:,i,:,:))
    end do

! perform update along the Y direction
!
    do j = 1, jm
      jm1 = max(1, j - 1)

      du(:,:,j,:) = du(:,:,j,:) + dhy * (f(iny,:,:,jm1,:) - f(iny,:,:,j,:))
    end do
#if NDIMS == 3

! perform update along the Z direction
!
    do k = 1, km
      km1 = max(1, k - 1)

      du(:,:,:,k) = du(:,:,:,k) + dhz * (f(inz,:,:,:,km1) - f(inz,:,:,:,k))
    end do
#endif /* NDIMS == 3 */

! update the solution for the fluid variables
!
    u(  1:nfl,:,:,:) = u(  1:nfl,:,:,:) + du(  1:nfl,:,:,:)

#ifdef MHD
! update the solution for the magnetic variables
!
    u(ibx:ibz,:,:,:) = u(ibx:ibz,:,:,:) + du(ibx:ibz,:,:,:)

#ifdef GLM
! calculate c_h^2
!
    ch2 = cmax * cmax

! update the solution for the scalar potential Psi
!
    u(iph,:,:,:) = u(iph,:,:,:) + ch2 * du(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

!-------------------------------------------------------------------------------
!
  end subroutine advance_solution
#ifdef EULER
!
!===============================================================================
!
! flux_euler: subroutine performs the first order integration of the numerical
!             flux
!
!===============================================================================
!
  subroutine flux_euler(pblock)

    use blocks   , only : block_data
    use mesh     , only : adxi, adyi, adzi
    use scheme   , only : update_flux

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! 1st step of integration
!
    call update_flux(pblock%u(:,:,:,:), pblock%f(:,:,:,:,:), dxi, dyi, dzi)

!-------------------------------------------------------------------------------
!
  end subroutine flux_euler
#endif /* EULER */
#ifdef RK2
!
!===============================================================================
!
! flux_rk2: subroutine performs integration of the numerical flux using
!           the second order Runge-Kutta method
!
!===============================================================================
!
  subroutine flux_rk2(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
    use mesh     , only : adxi, adyi, adzi
    use scheme   , only : update_flux
    use variables, only : nqt

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi

! local arrays
!
    real, dimension(      nqt,im,jm,km) :: u
    real, dimension(NDIMS,nqt,im,jm,km) :: f0, f1
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! calculate fluxes at the moment t
!
    call update_flux(u(:,:,:,:), f0(:,:,:,:,:), dxi, dyi, dzi)

! advance the solution to (t + dt) using computed fluxes in this substep
!
    call advance_solution(u(:,:,:,:), f0(:,:,:,:,:), dxi, dyi, dzi)

! calculate fluxes at the moment (t + dt)
!
    call update_flux(u(:,:,:,:), f1(:,:,:,:,:), dxi, dyi, dzi)

! calculate the time averaged flux
!
    pblock%f(:,:,:,:,:) = 0.5d0 * (f0(:,:,:,:,:) + f1(:,:,:,:,:))

!-------------------------------------------------------------------------------
!
  end subroutine flux_rk2
#endif /* RK2 */
#ifdef RK3
!
!===============================================================================
!
! flux_rk3: subroutine performs integration of the numerical flux using
!           the third order Runge-Kutta method
!
!===============================================================================
!
  subroutine flux_rk3(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km
    use mesh     , only : adxi, adyi, adzi
    use scheme   , only : update_flux
    use variables, only : nqt

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi

! local arrays
!
    real, dimension(      nqt,im,jm,km) :: u
    real, dimension(NDIMS,nqt,im,jm,km) :: f0, f1, f2
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! calculate fluxes at the moment t
!
    call update_flux(u(:,:,:,:), f0(:,:,:,:,:), dxi, dyi, dzi)

! advance the solution to (t + dt) using computed fluxes in this substep
!
    call advance_solution(u(:,:,:,:), f0(:,:,:,:,:), dxi, dyi, dzi)

! calculate fluxes at the moment (t + dt)
!
    call update_flux(u(:,:,:,:), f1(:,:,:,:,:), dxi, dyi, dzi)

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! average fluxes from t and t + dt and prepare for half step update
!
    f2(:,:,:,:,:) = 0.25d0 * (f0(:,:,:,:,:) + f1(:,:,:,:,:))

! advance the solution to (t + dt/2) using computed flux
!
    call advance_solution(u(:,:,:,:), f2(:,:,:,:,:), dxi, dyi, dzi)

! calculate fluxes at the moment (t + dt/2)
!
    call update_flux(u(:,:,:,:), f2(:,:,:,:,:), dxi, dyi, dzi)

! calculate the time averaged flux using Gauss formula
!
    pblock%f(:,:,:,:,:) = (f0(:,:,:,:,:) + f1(:,:,:,:,:)                       &
                                              + 4.0d0 * f2(:,:,:,:,:)) / 6.0d0

!-------------------------------------------------------------------------------
!
  end subroutine flux_rk3
#endif /* RK3 */
#else /* CONSERVATIVE */
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
#ifdef FORCE
    use forcing  , only : real_forcing
#endif /* FORCE */
    use mesh     , only : adxi, adyi, adzi
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
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi, ch2
#if defined MHD && defined GLM
    real    :: decay
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! 1st step of integration
!
    call update(pblock%u(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du)
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
#ifdef FORCE
! obtain the forcing terms in real space
!
    call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin    &
                                       , pblock%meta%zmin, f(:,:,:,:))

! update momenta due to the forcing terms
!
    pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:) + pblock%u(idn,:,:,:) * f(1,:,:,:)
    pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:) + pblock%u(idn,:,:,:) * f(2,:,:,:)
    pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:) + pblock%u(idn,:,:,:) * f(3,:,:,:)

#endif /* FORCE */
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
#ifdef FORCE
    use forcing  , only : real_forcing
#endif /* FORCE */
    use mesh     , only : adxi, adyi, adzi
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
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi, ch2
#if defined MHD && defined GLM
    real    :: decay
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

!! 1st step of integration
!!
    call update(pblock%u(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

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
! calculate c_h^2
!
    ch2 = cmax * cmax

! update the solution for the scalar potential Psi
!
    u1(iph,:,:,:) = pblock%u(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)

#endif /* GLM */
#endif /* MHD */
! 2nd step of integration
!
    call update(u1(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

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
#ifdef FORCE
! obtain the forcing terms in real space
!
    call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin    &
                                       , pblock%meta%zmin, f(:,:,:,:))

! update momenta due to the forcing terms
!
    pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:) + pblock%u(idn,:,:,:) * f(1,:,:,:)
    pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:) + pblock%u(idn,:,:,:) * f(2,:,:,:)
    pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:) + pblock%u(idn,:,:,:) * f(3,:,:,:)

#endif /* FORCE */
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
#ifdef FORCE
    use forcing  , only : real_forcing
#endif /* FORCE */
    use mesh     , only : adxi, adyi, adzi
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
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi
#if defined MHD && defined GLM
    real    :: decay, ch2
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#endif /* FORCE */

! parameters
!
    real, parameter :: f4 = 1.0d0 / 4.0d0, f3 = 1.0d0 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

!! 1st step of integration
!!
    call update(pblock%u(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

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
! calculate c_h^2
!
    ch2 = cmax * cmax

! update the solution for the scalar potential Psi
!
    u1(iph,:,:,:) = pblock%u(iph,:,:,:) + ch2 * dt * du(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

!! 2nd step of integration
!!
    call update(u1(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

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
    call update(u1(:,:,:,:), du(:,:,:,:), dxi, dyi, dzi)

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
#ifdef FORCE
! obtain the forcing terms in real space
!
    call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin    &
                                       , pblock%meta%zmin, f(:,:,:,:))

! update momenta due to the forcing terms
!
    pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:) + pblock%u(idn,:,:,:) * f(1,:,:,:)
    pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:) + pblock%u(idn,:,:,:) * f(2,:,:,:)
    pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:) + pblock%u(idn,:,:,:) * f(3,:,:,:)

#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk3
#endif /* RK3 */
#endif /* CONSERVATIVE */
!
!===============================================================================
!
end module
