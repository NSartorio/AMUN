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
!! module: SOURCES
!!
!!  This modules adds source terms.
!!
!!******************************************************************************
!
module sources

#ifdef PROFILE
! include external procedures
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer, save :: imi, imu
#endif /* PROFILE */

! gravitational acceleration coefficient
!
  real(kind=8), save :: gpoint      = 0.0d+00

! viscosity coefficient
!
  real(kind=8), save :: viscosity   = 0.0d+00

! resistivity coefficient
!
  real(kind=8), save :: resistivity = 0.0d+00

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_sources, finalize_sources
  public :: update_sources
  public :: viscosity, resistivity

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
! subroutine INITIALIZE_SOURCES:
! -----------------------------
!
!   Subroutine initializes module SOURCES.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_sources(verbose, iret)

! include external procedures and variables
!
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('sources:: initialize', imi)
    call set_timer('sources:: update'    , imu)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get acceleration coefficient
!
      call get_parameter_real("gpoint"     , gpoint   )

! get viscosity coefficient
!
      call get_parameter_real("viscosity"  , viscosity)

! get resistivity coefficient
!
      call get_parameter_real("resistivity", resistivity)

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,1e9.2)") "point mass constant    =", gpoint
      write (*,"(4x,a,1x,1e9.2)") "viscosity              =", viscosity
      write (*,"(4x,a,1x,1e9.2)") "resistivity            =", resistivity

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_sources
!
!===============================================================================
!
! subroutine FINALIZE_SOURCES:
! ---------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_sources(iret)

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

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_sources
!
!===============================================================================
!
! subroutine UPDATE_SOURCES:
! -------------------------
!
!   Subroutine add the source terms.
!
!   Arguments:
!
!     q    - the array of primitive variables;
!     du   - the array of variable increment;
!
!===============================================================================
!
  subroutine update_sources(pdata, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az, adx, ady, adz, adxi, adyi, adzi
    use equations      , only : nv, inx, iny, inz
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ien
    use equations      , only : ibx, iby, ibz
    use operators      , only : divergence, curl

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: du

! local variables
!
    integer       :: i  , j  , k
    integer       :: im1, jm1, km1
    integer       :: ip1, jp1, kp1
    real(kind=8)  :: r2, r3, gc, gx, gy, gz
    real(kind=8)  :: dxi2, dyi2, dzi2, dvx, dvy, dvz, dbx, dby, dbz, dvv

! local arrays
!
    real(kind=8), dimension(3)  :: dh
    real(kind=8), dimension(im) :: x
    real(kind=8), dimension(jm) :: y
    real(kind=8), dimension(km) :: z
    real(kind=8), dimension(3,im,jm,km) :: jc, ejc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for source terms
!
    call start_timer(imu)
#endif /* PROFILE */

! proceed only if the gravitational acceleration coefficient is not zero
!
    if (gpoint /= 0.0d+00) then

! prepare block coordinates
!
      x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
      y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
#if NDIMS == 3
      z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
#endif /* NDIMS == 3 */

! iterate over all positions in the YZ plane
!
      do k = 1, km
        do j = 1, jm
          do i = 1, jm

! calculate distance from the origin
!
#if NDIMS == 2
            r2 = x(i) * x(i) + y(j) * y(j)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            r2 = x(i) * x(i) + y(j) * y(j) + z(k) * z(k)
#endif /* NDIMS == 3 */
            r3 = r2 * sqrt(r2)

! calculate gravitational acceleration factors
!
            gc = gpoint * pdata%q(idn,i,j,k) / max(1.0d-16, r3)
            gx = gc * x(i)
            gy = gc * y(j)
#if NDIMS == 3
            gz = gc * z(k)
#endif /* NDIMS == 3 */

! add source terms to momentum equations
!
            du(imx,i,j,k) = du(imx,i,j,k) + gx
            du(imy,i,j,k) = du(imy,i,j,k) + gy
#if NDIMS == 3
            du(imz,i,j,k) = du(imz,i,j,k) + gz
#endif /* NDIMS == 3 */

! add source terms to total energy equation
!
            if (ien > 0) then

#if NDIMS == 2
              du(ien,i,j,k) = du(ien,i,j,k) + gx * pdata%q(ivx,i,j,k)          &
                                            + gy * pdata%q(ivy,i,j,k)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              du(ien,i,j,k) = du(ien,i,j,k) + gx * pdata%q(ivx,i,j,k)          &
                                            + gy * pdata%q(ivy,i,j,k)          &
                                            + gz * pdata%q(ivz,i,j,k)
#endif /* NDIMS == 3 */

            end if

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

    end if ! gpoint is not zero

! proceed only if the viscosity coefficient is not zero
!
    if (viscosity > 0.0d+00) then

! prepare coordinate increments
!
      dxi2 = viscosity * adxi(pdata%meta%level)**2
      dyi2 = viscosity * adyi(pdata%meta%level)**2
#if NDIMS == 3
      dzi2 = viscosity * adzi(pdata%meta%level)**2
#endif /* NDIMS == 3 */

! iterate over all positions in the YZ plane
!
      do k = 1, km
#if NDIMS == 3
        km1 = max( 1, k - 1)
        kp1 = min(km, k + 1)
#endif /* NDIMS == 3 */
        do j = 1, jm
          jm1 = max( 1, j - 1)
          jp1 = min(jm, j + 1)
          do i = 1, jm
            im1 = max( 1, i - 1)
            ip1 = min(im, i + 1)

! calculate second order derivatives of Vx
!
            dvx = (pdata%q(ivx,ip1,j,k) + pdata%q(ivx,im1,j,k))                &
                                                - 2.0d+00 * pdata%q(ivx,i,j,k)
            dvy = (pdata%q(ivx,i,jp1,k) + pdata%q(ivx,i,jm1,k))                &
                                                - 2.0d+00 * pdata%q(ivx,i,j,k)
#if NDIMS == 3
            dvz = (pdata%q(ivx,i,j,kp1) + pdata%q(ivx,i,j,km1))                &
                                                - 2.0d+00 * pdata%q(ivx,i,j,k)
#endif /* NDIMS == 3 */

! calculate the source term for Vx
!
#if NDIMS == 2
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy + dzi2 * dvz)
#endif /* NDIMS == 3 */

! add viscous source terms to X-momentum equation
!
            du(imx,i,j,k) = du(imx,i,j,k) + dvv

! add viscous source term to total energy equation
!
            if (ien > 0) then
              du(ien,i,j,k) = du(ien,i,j,k) + pdata%q(ivx,i,j,k) * dvv
            end if

! calculate second order derivatives of Vy
!
            dvx = (pdata%q(ivy,ip1,j,k) + pdata%q(ivy,im1,j,k))                &
                                                - 2.0d+00 * pdata%q(ivy,i,j,k)
            dvy = (pdata%q(ivy,i,jp1,k) + pdata%q(ivy,i,jm1,k))                &
                                                - 2.0d+00 * pdata%q(ivy,i,j,k)
#if NDIMS == 3
            dvz = (pdata%q(ivy,i,j,kp1) + pdata%q(ivy,i,j,km1))                &
                                                - 2.0d+00 * pdata%q(ivy,i,j,k)
#endif /* NDIMS == 3 */

! calculate the source term for Vy
!
#if NDIMS == 2
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy + dzi2 * dvz)
#endif /* NDIMS == 3 */

! add viscous source terms to Y-momentum equation
!
            du(imy,i,j,k) = du(imy,i,j,k) + dvv

! add viscous source term to total energy equation
!
            if (ien > 0) then
              du(ien,i,j,k) = du(ien,i,j,k) + pdata%q(ivy,i,j,k) * dvv
            end if

! calculate second order derivatives of Vz
!
            dvx = (pdata%q(ivz,ip1,j,k) + pdata%q(ivz,im1,j,k))                &
                                                - 2.0d+00 * pdata%q(ivz,i,j,k)
            dvy = (pdata%q(ivz,i,jp1,k) + pdata%q(ivz,i,jm1,k))                &
                                                - 2.0d+00 * pdata%q(ivz,i,j,k)
#if NDIMS == 3
            dvz = (pdata%q(ivz,i,j,kp1) + pdata%q(ivz,i,j,km1))                &
                                                - 2.0d+00 * pdata%q(ivz,i,j,k)
#endif /* NDIMS == 3 */

! calculate the source term for Vz
!
#if NDIMS == 2
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            dvv = pdata%q(idn,i,j,k) * (dxi2 * dvx + dyi2 * dvy + dzi2 * dvz)
#endif /* NDIMS == 3 */

! add viscous source terms to Y-momentum equation
!
            du(imz,i,j,k) = du(imz,i,j,k) + dvv

! add viscous source term to total energy equation
!
            if (ien > 0) then
              du(ien,i,j,k) = du(ien,i,j,k) + pdata%q(ivz,i,j,k) * dvv
            end if

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

    end if ! viscosity is not zero

! proceed only if the resistivity coefficient is not zero
!
    if (resistivity > 0.0d+00 .and. ibx > 0) then

! prepare coordinate increments
!
      dh(1) = adx(pdata%meta%level)
      dh(2) = ady(pdata%meta%level)
      dh(3) = adz(pdata%meta%level)

! calculate current density J = ∇xB
!
      call curl(dh(:), pdata%q(ibx:ibz,:,:,:), jc(inx:inz,:,:,:))

! calculate (η J)
!
      ejc(inx:inz,:,:,:) = resistivity * jc(inx:inz,:,:,:)

! calculate curl of (η J)
!
      call curl(dh(:), ejc(inx:inz,:,:,:), jc(inx:inz,:,:,:))

! update magnetic field component increments
!
      du(ibx,:,:,:) = du(ibx,:,:,:) - jc(inx,:,:,:)
      du(iby,:,:,:) = du(iby,:,:,:) - jc(iny,:,:,:)
      du(ibz,:,:,:) = du(ibz,:,:,:) - jc(inz,:,:,:)

! update energy equation
!
      if (ien > 0) then

! calculate B x (η J)
!
        jc(inx,:,:,:) = pdata%q(iby,:,:,:) * ejc(inz,:,:,:)                    &
                                         - pdata%q(ibz,:,:,:) * ejc(iny,:,:,:)
        jc(iny,:,:,:) = pdata%q(ibz,:,:,:) * ejc(inx,:,:,:)                    &
                                         - pdata%q(ibx,:,:,:) * ejc(inz,:,:,:)
        jc(inz,:,:,:) = pdata%q(ibx,:,:,:) * ejc(iny,:,:,:)                    &
                                         - pdata%q(iby,:,:,:) * ejc(inx,:,:,:)

! calculate the divergence of B x (η J), i.e. ∇.[B x (η J)]
!
        call divergence(dh(:), jc(:,:,:,:), ejc(inx,:,:,:))

! add the resistive source term to the energy equation, i.e.
!     d/dt E + ∇.F = ∇.(B x (η J))
!
        du(ien,:,:,:) = du(ien,:,:,:) + ejc(inx,:,:,:)

      end if ! energy equation present

    end if ! resistivity is not zero

#ifdef PROFILE
! stop accounting time for source terms
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_sources

!===============================================================================
!
end module sources
