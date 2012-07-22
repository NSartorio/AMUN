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
!! module: EVOLUTION - handling the time evolution of the block structure
!!
!!******************************************************************************
!
module evolution

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn, dxmin

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
    use config    , only : toplev
#endif /* REFINE */
    use mesh      , only : update_mesh
#ifdef FORCE
    use config    , only : tbfor
    use forcing   , only : fourier_transform, evolve_forcing
    use variables , only : idn, imz
#endif /* FORCE */

    implicit none

! local variables
!
    type(block_data), pointer :: pblock
    real                      :: cm
!
!-------------------------------------------------------------------------------
!
#ifdef FORCE
! perform forcing evolution only when t >= tbfor
!
    if (t .ge. tbfor) then

! perform the Fourier transform of the velocity field
!
      pblock => list_data
      do while (associated(pblock))

        if (pblock%meta%leaf)                                                  &
          call fourier_transform(pblock%meta%level                             &
                         , pblock%meta%xmin, pblock%meta%ymin                  &
                         , pblock%meta%zmin, pblock%u(idn:imz,:,:,:))

        pblock => pblock%next ! assign pointer to the next block

      end do

! evolve the forcing source terms by the time interval dt
!
      call evolve_forcing(t, dt)

    end if ! t >= tbfor
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
    call boundary_correct_fluxes()

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
    call boundary_variables()

#ifdef REFINE
! chec if we need to perform the refinement step
!
    if (toplev .gt. 1) then

! check refinement and refine
!
      call update_mesh()

! update boundaries
!
      call boundary_variables()

    end if ! toplev > 1
#endif /* REFINE */

! find new time step
!
    call find_new_timestep()

!-------------------------------------------------------------------------------
!
  end subroutine evolve
!
!===============================================================================
!
! find_new_timestep: subroutine updates the maximum speed among the leafs and
!                    calculates new time step
!
!===============================================================================
!
  subroutine find_new_timestep()

    use blocks  , only : block_meta, block_data, list_meta, list_data
    use config  , only : toplev
#ifdef MPI
    use mpitools, only : reduce_maximum_real
#endif /* MPI */
    use coords  , only : adx, ady, adz
    use scheme  , only : maxspeed, cmax
#ifdef VISCOSITY
    use config  , only : visc
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    use config  , only : ueta
#endif /* MHD & RESISTIVITY */

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
        cm = maxspeed(pdata%u(:,:,:,:))

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
#ifdef VISCOSITY
    dtn = min(dtn, 0.5d0 * dxmin * dxmin / max(1.0d-16, visc))
#endif /* VISCOSITY */
#if defined MHD && defined RESISTIVITY
    dtn = min(dtn, 0.5d0 * dxmin * dxmin / max(1.0d-16, ueta))
#endif /* MHD & RESISTIVITY */

!-------------------------------------------------------------------------------
!
  end subroutine find_new_timestep
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
    use coords   , only : adxi, adyi, adzi
#if defined MHD && defined GLM
    use config   , only : decay
    use scheme   , only : cmax
    use variables, only : iph
#endif /* MHD & GLM */
#ifdef FORCE
    use config   , only : tbfor
    use forcing  , only : real_forcing
    use variables, only : inx, iny, inz
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */
#endif /* FORCE */
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi
#ifdef FORCE

! local arrays
!
    real, dimension(  3,im,jm,km) :: f
#ifdef ADI
    real, dimension(    im,jm,km) :: ek, dek
#endif /* ADI */
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
    call advance_solution(pblock%u(:,:,:,:), pblock%f(:,:,:,:,:)               &
                                                          , dt, dxi, dyi, dzi)
#if defined MHD && defined GLM

! evolve Psi due to the source term
!
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* MHD & GLM */

#ifdef FORCE
! add forcing term only if t >= tbfor
!
    if (t .ge. tbfor) then

! obtain the forcing terms in real space
!
      call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin  &
                                         , pblock%meta%zmin, f(:,:,:,:))

#ifdef ADI
! calculate kinetic energy before adding the forcing term
!
      ek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2     &
                               + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:)
#endif /* ADI */

! update momenta due to the forcing terms
!
      pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:)                                &
                                          + pblock%u(idn,:,:,:) * f(inx,:,:,:)
      pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:)                                &
                                          + pblock%u(idn,:,:,:) * f(iny,:,:,:)
      pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:)                                &
                                          + pblock%u(idn,:,:,:) * f(inz,:,:,:)

#ifdef ADI
! calculate kinetic energy after adding the forcing term
!
      dek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2    &
                   + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:) - ek(:,:,:)

! update total energy with the injected one
!
      pblock%u(ien,:,:,:) = pblock%u(ien,:,:,:) + dek(:,:,:)
#endif /* ADI */

    end if ! t >= tbfor
#endif /* FORCE */

#ifdef SHAPE
! update solid shapes
!
    call update_shapes(pblock, t)
#endif /* SHAPE */

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
  subroutine advance_solution(u, f, dh, dxi, dyi, dzi)

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
    real                                               :: dh, dxi, dyi, dzi

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
    dhx = dh * dxi
    dhy = dh * dyi
#if NDIMS == 3
    dhz = dh * dzi
#endif /* NDIMS == 3 */

! reset the increment array du
!
    du(:,:,:,:) = 0.0d0

! perform update along the X direction
!
    do i = 2, im
      im1 = i - 1

      du(:,i,:,:) = du(:,i,:,:) - dhx * (f(inx,:,i,:,:) - f(inx,:,im1,:,:))
    end do
    du(:,1,:,:) = du(:,1,:,:) - dhx *  f(inx,:,1,:,:)

! perform update along the Y direction
!
    do j = 2, jm
      jm1 = j - 1

      du(:,:,j,:) = du(:,:,j,:) - dhy * (f(iny,:,:,j,:) - f(iny,:,:,jm1,:))
    end do
    du(:,:,1,:) = du(:,:,1,:) - dhy * f(iny,:,:,1,:)

#if NDIMS == 3
! perform update along the Z direction
!
    do k = 2, km
      km1 = k - 1

      du(:,:,:,k) = du(:,:,:,k) - dhz * (f(inz,:,:,:,k) - f(inz,:,:,:,km1))
    end do
    du(:,:,:,1) = du(:,:,:,1) - dhz * f(inz,:,:,:,1)
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
!
!===============================================================================
!
! advance_solution_1d: subroutine performs an one step update of the conserved
!                      variables array using the numerical fluxes passed as an
!                      argument along one selected direction only
!
!===============================================================================
!
  subroutine advance_solution_1d(idir, dh, dxi, u, f)

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
    integer                      , intent(in)    :: idir
    real                         , intent(in)    :: dh, dxi
    real, dimension(nqt,im,jm,km), intent(inout) :: u
    real, dimension(nqt,im,jm,km), intent(in)    :: f

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
! reset the increment array du
!
    du(:,:,:,:) = 0.0d0

! calculate the conserved variables increment
!
    select case(idir)
    case(1)

! prepare dxi
!
      dhx = dh * dxi

! perform update along the X direction
!
      do i = 2, im
        im1 = i - 1

        du(:,i,:,:) = du(:,i,:,:) - dhx * (f(:,i,:,:) - f(:,im1,:,:))
      end do
      du(:,1,:,:) = du(:,1,:,:) - dhx *  f(:,1,:,:)

    case(2)

! prepare dxi
!
      dhy = dh * dxi

! perform update along the Y direction
!
      do j = 2, jm
        jm1 = j - 1

        du(:,:,j,:) = du(:,:,j,:) - dhy * (f(:,:,j,:) - f(:,:,jm1,:))
      end do
      du(:,:,1,:) = du(:,:,1,:) - dhy * f(:,:,1,:)

#if NDIMS == 3
    case(3)

! prepare dxi
!
      dhz = dh * dxi

! perform update along the Z direction
!
      do k = 2, km
        km1 = k - 1

        du(:,:,:,k) = du(:,:,:,k) - dhz * (f(:,:,:,k) - f(:,:,:,km1))
      end do
      du(:,:,:,1) = du(:,:,:,1) - dhz * f(:,:,:,1)
#endif /* NDIMS == 3 */

    end select

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
  end subroutine advance_solution_1d
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
    use coords   , only : adx, ady, adz
    use scheme   , only : update_flux

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real :: dx, dy, dz
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dx = adx(pblock%meta%level)
    dy = ady(pblock%meta%level)
    dz = adz(pblock%meta%level)

! 1st step of integration
!
    call update_flux(1, dx, pblock%u(:,:,:,:), pblock%f(1,:,:,:,:))
    call update_flux(2, dy, pblock%u(:,:,:,:), pblock%f(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, pblock%u(:,:,:,:), pblock%f(3,:,:,:,:))
#endif /* NDIMS == 3 */

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
    use coords   , only : adx, ady, adz, adxi, adyi, adzi
    use scheme   , only : update_flux
    use variables, only : nqt

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real :: dx, dy, dz, dxi, dyi, dzi

! local arrays
!
    real, dimension(      nqt,im,jm,km) :: u
    real, dimension(NDIMS,nqt,im,jm,km) :: f0, f1
!
!-------------------------------------------------------------------------------
!
! prepare dx, dy, dz, dxi, dyi, and dzi
!
    dx  = adx (pblock%meta%level)
    dy  = ady (pblock%meta%level)
    dz  = adz (pblock%meta%level)
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! calculate fluxes at the moment t
!
    call update_flux(1, dx, u(:,:,:,:), f0(1,:,:,:,:))
    call update_flux(2, dy, u(:,:,:,:), f0(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, u(:,:,:,:), f0(3,:,:,:,:))
#endif /* NDIMS == 3 */

! advance the solution to (t + dt) using the computed fluxes
!
    call advance_solution(u(:,:,:,:), f0(:,:,:,:,:), dt, dxi, dyi, dzi)

! calculate fluxes at the moment (t + dt/2)
!
    call update_flux(1, dx, u(:,:,:,:), f1(1,:,:,:,:))
    call update_flux(2, dy, u(:,:,:,:), f1(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, u(:,:,:,:), f1(3,:,:,:,:))
#endif /* NDIMS == 3 */

! average the flux at the time )t + dt/2)
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
    use coords   , only : adx, ady, adz, adxi, adyi, adzi
    use scheme   , only : update_flux
    use variables, only : nqt

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer :: lev
    real    :: dth, dx, dy, dz, dxi, dyi, dzi

! local arrays
!
    real, dimension(      nqt,im,jm,km) :: u
    real, dimension(NDIMS,nqt,im,jm,km) :: f0, f1, f2
!
!-------------------------------------------------------------------------------
!
! obtain the block level
!
    lev    = pblock%meta%level

! calculate the half time step
!
    dth = 0.5d0 * dt

! prepare dx, dy, dz, dxi, dyi, and dzi
!
    dx     = adx (lev)
    dy     = ady (lev)
    dz     = adz (lev)
    dxi    = adxi(lev)
    dyi    = adyi(lev)
    dzi    = adzi(lev)

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! calculate fluxes at time t
!
    call update_flux(1, dx, u(:,:,:,:), f0(1,:,:,:,:))
    call update_flux(2, dy, u(:,:,:,:), f0(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, u(:,:,:,:), f0(3,:,:,:,:))
#endif /* NDIMS == 3 */

! advance the solution to (t + dt) using the computed fluxes
!
    call advance_solution(u(:,:,:,:), f0(:,:,:,:,:), dt, dxi, dyi, dzi)

! calculate fluxes at time (t + dt)
!
    call update_flux(1, dx, u(:,:,:,:), f1(1,:,:,:,:))
    call update_flux(2, dy, u(:,:,:,:), f1(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, u(:,:,:,:), f1(3,:,:,:,:))
#endif /* NDIMS == 3 */

! average fluxes at the time (t + dt / 2)
!
    f2(:,:,:,:,:) = 0.5d0 * (f0(:,:,:,:,:) + f1(:,:,:,:,:))

! copy the initial state to the local array u
!
    u(:,:,:,:) = pblock%u(:,:,:,:)

! advance the solution to (t + dt / 2) using the computed fluxes
!
    call advance_solution(u(:,:,:,:), f2(:,:,:,:,:), dth, dxi, dyi, dzi)

! calculate fluxes at time (t + dt / 2)
!
    call update_flux(1, dx, u(:,:,:,:), f2(1,:,:,:,:))
    call update_flux(2, dy, u(:,:,:,:), f2(2,:,:,:,:))
#if NDIMS == 3
    call update_flux(3, dz, u(:,:,:,:), f2(3,:,:,:,:))
#endif /* NDIMS == 3 */

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
    use config   , only : tbfor
    use forcing  , only : real_forcing
#endif /* FORCE */
    use coords   , only : adxi, adyi, adzi
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : decay
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi, ch2

! local arrays
!
    real, dimension(nqt,im,jm,km) :: du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#ifdef ADI
    real, dimension(    im,jm,km) :: ek, dek
#endif /* ADI */
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
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

#ifdef FORCE
! add the forcing term only if t >= tbfor
!
    if (t .ge. tbfor) then

! obtain the forcing terms in real space
!
      call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin  &
                                       , pblock%meta%zmin, f(:,:,:,:))

#ifdef ADI
! calculate kinetic energy before adding the forcing term
!
      ek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2     &
                               + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:)
#endif /* ADI */

! update momenta due to the forcing terms
!
      pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(1,:,:,:)
      pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(2,:,:,:)
      pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(3,:,:,:)

#ifdef ADI
! calculate kinetic energy after adding the forcing term
!
      dek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2    &
                   + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:) - ek(:,:,:)

! update total energy with the injected one
!
      pblock%u(ien,:,:,:) = pblock%u(ien,:,:,:) + dek(:,:,:)
#endif /* ADI */

    end if ! t >= tbfor
#endif /* FORCE */

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, t)
#endif /* SHAPE */

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
    use config   , only : tbfor
    use forcing  , only : real_forcing
#endif /* FORCE */
    use coords   , only : adxi, adyi, adzi
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : decay
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi, ch2

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#ifdef ADI
    real, dimension(    im,jm,km) :: ek, dek
#endif /* ADI */
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
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)

#endif /* GLM */
#endif /* MHD */

#ifdef FORCE
! add the forcing term only if t >= tbfor
!
    if (t .ge. tbfor) then

! obtain the forcing terms in real space
!
      call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin  &
                                       , pblock%meta%zmin, f(:,:,:,:))

#ifdef ADI
! calculate kinetic energy before adding the forcing term
!
      ek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2     &
                               + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:)
#endif /* ADI */

! update momenta due to the forcing terms
!
      pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(1,:,:,:)
      pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(2,:,:,:)
      pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(3,:,:,:)

#ifdef ADI
! calculate kinetic energy after adding the forcing term
!
      dek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2    &
                   + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:) - ek(:,:,:)

! update total energy with the injected one
!
      pblock%u(ien,:,:,:) = pblock%u(ien,:,:,:) + dek(:,:,:)

    end if ! t >= tbfor
#endif /* ADI */
#endif /* FORCE */

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, t)

#endif /* SHAPE */

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
    use config   , only : tbfor
    use forcing  , only : real_forcing
#endif /* FORCE */
    use coords   , only : adxi, adyi, adzi
#ifdef SHAPE
    use problem  , only : update_shapes
#endif /* SHAPE */
    use scheme   , only : update, cmax
    use variables, only : nqt, nfl
#ifdef MHD
    use variables, only : ibx, ibz
#ifdef GLM
    use config   , only : decay
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */
#ifdef FORCE
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */
#endif /* FORCE */

    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi
#if defined MHD && defined GLM
    real    :: ch2
#endif /* MHD & GLM */

! local arrays
!
    real, dimension(nqt,im,jm,km) :: u1, du
#ifdef FORCE
    real, dimension(  3,im,jm,km) :: f
#ifdef ADI
    real, dimension(    im,jm,km) :: ek, dek
#endif /* ADI */
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
    pblock%u(iph,:,:,:) = decay * pblock%u(iph,:,:,:)
#endif /* GLM */
#endif /* MHD */

#ifdef FORCE
! add the forcing term only if t >= tbfor
!
    if (t .ge. tbfor) then

! obtain the forcing terms in real space
!
      call real_forcing(pblock%meta%level, pblock%meta%xmin, pblock%meta%ymin  &
                                       , pblock%meta%zmin, f(:,:,:,:))

#ifdef ADI
! calculate kinetic energy before adding the forcing term
!
      ek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2     &
                               + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:)
#endif /* ADI */

! update momenta due to the forcing terms
!
      pblock%u(imx,:,:,:) = pblock%u(imx,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(1,:,:,:)
      pblock%u(imy,:,:,:) = pblock%u(imy,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(2,:,:,:)
      pblock%u(imz,:,:,:) = pblock%u(imz,:,:,:)                                &
                                            + pblock%u(idn,:,:,:) * f(3,:,:,:)

#ifdef ADI
! calculate kinetic energy after adding the forcing term
!
      dek(:,:,:) = 0.5d0 * (pblock%u(imx,:,:,:)**2 + pblock%u(imy,:,:,:)**2    &
                   + pblock%u(imz,:,:,:)**2) / pblock%u(idn,:,:,:) - ek(:,:,:)

! update total energy with the injected one
!
      pblock%u(ien,:,:,:) = pblock%u(ien,:,:,:) + dek(:,:,:)

    end if ! t >= tbfor
#endif /* ADI */
#endif /* FORCE */

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, t)

#endif /* SHAPE */

!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk3
#endif /* RK3 */
#endif /* CONSERVATIVE */
!
!===============================================================================
!
end module
