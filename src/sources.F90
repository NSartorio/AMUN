!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2017 Grzegorz Kowal <grzegorz@amuncode.org>
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

! GLM-MHD source terms type (1 - EGLM, 2 - HEGLM)
!
  integer     , save :: glm_type    = 0

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
    use parameters , only : get_parameter_string, get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=8)       :: tglm = "none"
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

! get the type of the GLM source terms
!
    call get_parameter_string("glm_source_terms", tglm)

! set the glm_type variable to correct value
!
    select case(trim(tglm))
    case("eglm", "EGLM")
      glm_type = 1

    case("heglm", "HEGLM")
      glm_type = 2

    case default
      glm_type = 0
    end select

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

      write (*,"(4x,a,1x,a)    ") "glm source terms       =", trim(tglm)
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
    use coordinates    , only : ax, ay, az, adx, ady, adz
    use equations      , only : nv, inx, iny, inz
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ien
    use equations      , only : ibx, iby, ibz, ibp
    use operators      , only : divergence, gradient, laplace, curl

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
    real(kind=8)  :: fc, gc
    real(kind=8)  :: r2, r3, gx, gy, gz
    real(kind=8)  :: dbx, dby, dbz
    real(kind=8)  :: dvxdx, dvxdy, dvxdz, divv
    real(kind=8)  :: dvydx, dvydy, dvydz
    real(kind=8)  :: dvzdx, dvzdy, dvzdz

! local arrays
!
    real(kind=8), dimension(3)  :: dh
    real(kind=8), dimension(im) :: x
    real(kind=8), dimension(jm) :: y
    real(kind=8), dimension(km) :: z
    real(kind=8), dimension(im,jm,km)     :: db
    real(kind=8), dimension(3,3,im,jm,km) :: tmp
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
          do i = 1, im

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
      dh(1) = adx(pdata%meta%level)
      dh(2) = ady(pdata%meta%level)
      dh(3) = adz(pdata%meta%level)

! calculate the velocity Jacobian
!
      call gradient(dh(:), pdata%q(ivx,1:im,1:jm,1:km)                         &
                                            , tmp(inx,inx:inz,1:im,1:jm,1:km))
      call gradient(dh(:), pdata%q(ivy,1:im,1:jm,1:km)                         &
                                            , tmp(iny,inx:inz,1:im,1:jm,1:km))
      call gradient(dh(:), pdata%q(ivz,1:im,1:jm,1:km)                         &
                                            , tmp(inz,inx:inz,1:im,1:jm,1:km))

! iterate over all cells
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im

! prepare the νρ factor
!
            gc    = viscosity * pdata%q(idn,i,j,k)
            fc    = 2.0d+00 * gc

! get the velocity Jacobian elements
!
            dvxdx = tmp(inx,inx,i,j,k)
            dvxdy = tmp(inx,iny,i,j,k)
            dvxdz = tmp(inx,inz,i,j,k)
            dvydx = tmp(iny,inx,i,j,k)
            dvydy = tmp(iny,iny,i,j,k)
            dvydz = tmp(iny,inz,i,j,k)
            dvzdx = tmp(inz,inx,i,j,k)
            dvzdy = tmp(inz,iny,i,j,k)
            dvzdz = tmp(inz,inz,i,j,k)
            divv  = (dvxdx + dvydy + dvzdz) / 3.0d+00

! calculate elements of the viscous stress tensor
!
            tmp(inx,inx,i,j,k) = fc * (dvxdx - divv)
            tmp(iny,iny,i,j,k) = fc * (dvydy - divv)
            tmp(inz,inz,i,j,k) = fc * (dvzdz - divv)
            tmp(inx,iny,i,j,k) = gc * (dvxdy + dvydx)
            tmp(inx,inz,i,j,k) = gc * (dvxdz + dvzdx)
            tmp(iny,inz,i,j,k) = gc * (dvydz + dvzdy)
            tmp(iny,inx,i,j,k) = tmp(inx,iny,i,j,k)
            tmp(inz,inx,i,j,k) = tmp(inx,inz,i,j,k)
            tmp(inz,iny,i,j,k) = tmp(iny,inz,i,j,k)

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

! calculate the divergence of the first tensor row
!
      call divergence(dh(:), tmp(inx,inx:inz,1:im,1:jm,1:km)                   &
                                                         , db(1:im,1:jm,1:km))

! add viscous source terms to the X momentum equation
!
      du(imx,1:im,1:jm,1:km) = du(imx,1:im,1:jm,1:km) + db(1:im,1:jm,1:km)

! calculate the divergence of the second tensor row
!
      call divergence(dh(:), tmp(iny,inx:inz,1:im,1:jm,1:km)                   &
                                                         , db(1:im,1:jm,1:km))

! add viscous source terms to the Y momentum equation
!
      du(imy,1:im,1:jm,1:km) = du(imy,1:im,1:jm,1:km) + db(1:im,1:jm,1:km)

! calculate the divergence of the third tensor row
!
      call divergence(dh(:), tmp(inz,inx:inz,1:im,1:jm,1:km)                   &
                                                         , db(1:im,1:jm,1:km))

! add viscous source terms to the Z momentum equation
!
      du(imz,1:im,1:jm,1:km) = du(imz,1:im,1:jm,1:km) + db(1:im,1:jm,1:km)

! add viscous source term to total energy equation
!
      if (ien > 0) then

! iterate over all cells
!
        do k = 1, km
          do j = 1, jm
            do i = 1, im

! calculate scalar product of v and viscous stress tensor τ
!
              gx = pdata%q(ivx,i,j,k) * tmp(inx,inx,i,j,k)                     &
                 + pdata%q(ivy,i,j,k) * tmp(inx,iny,i,j,k)                     &
                 + pdata%q(ivz,i,j,k) * tmp(inx,inz,i,j,k)
              gy = pdata%q(ivx,i,j,k) * tmp(iny,inx,i,j,k)                     &
                 + pdata%q(ivy,i,j,k) * tmp(iny,iny,i,j,k)                     &
                 + pdata%q(ivz,i,j,k) * tmp(iny,inz,i,j,k)
              gz = pdata%q(ivx,i,j,k) * tmp(inz,inx,i,j,k)                     &
                 + pdata%q(ivy,i,j,k) * tmp(inz,iny,i,j,k)                     &
                 + pdata%q(ivz,i,j,k) * tmp(inz,inz,i,j,k)

! update (v.τ), use the first row of the tensor tmp
!
              tmp(inx,inx,i,j,k) = gx
              tmp(inx,iny,i,j,k) = gy
              tmp(inx,inz,i,j,k) = gz

            end do ! i = 1, im
          end do ! j = 1, jm
        end do ! k = 1, km

! calculate the divergence of (v.τ)
!
        call divergence(dh(:), tmp(inx,inx:inz,1:im,1:jm,1:km)                 &
                                                         , db(1:im,1:jm,1:km))

! update the energy increment
!
        du(ien,1:im,1:jm,1:km) = du(ien,1:im,1:jm,1:km) + db(1:im,1:jm,1:km)

      end if ! ien > 0

    end if ! viscosity is not zero

!=== add magnetic field related source terms ===
!
    if (ibx > 0) then

! prepare coordinate increments
!
      dh(1) = adx(pdata%meta%level)
      dh(2) = ady(pdata%meta%level)
      dh(3) = adz(pdata%meta%level)

! add the EGLM-MHD source terms
!
      if (glm_type > 0) then

! calculate the magnetic field divergence
!
        call divergence(dh(:), pdata%q(ibx:ibz,:,:,:), db(:,:,:))

! update the momentum component increments, i.e.
!     d/dt (ρv) + ∇.F = - (∇.B)B
!
        du(imx,:,:,:) = du(imx,:,:,:) - db(:,:,:) * pdata%q(ibx,:,:,:)
        du(imy,:,:,:) = du(imy,:,:,:) - db(:,:,:) * pdata%q(iby,:,:,:)
        du(imz,:,:,:) = du(imz,:,:,:) - db(:,:,:) * pdata%q(ibz,:,:,:)

! update the energy equation
!
        if (ien > 0 .and. ibp > 0) then

! calculate the gradient of divergence potential
!
          call gradient(dh(:), pdata%q(ibp,:,:,:), tmp(inx:inz,inx,:,:,:))

! add the divergence potential source term to the energy equation, i.e.
!     d/dt E + ∇.F = - B.(∇ψ)
!
          du(ien,:,:,:) = du(ien,:,:,:)                                        &
                     - sum(pdata%q(ibx:ibz,:,:,:) * tmp(inx:inz,inx,:,:,:), 1)

        end if ! ien > 0

! add the HEGLM-MHD source terms
!
        if (glm_type > 1) then

! update magnetic field component increments, i.e.
!     d/dt B + ∇.F = - (∇.B)v
!
          du(ibx,:,:,:) = du(ibx,:,:,:) - db(:,:,:) * pdata%q(ivx,:,:,:)
          du(iby,:,:,:) = du(iby,:,:,:) - db(:,:,:) * pdata%q(ivy,:,:,:)
          du(ibz,:,:,:) = du(ibz,:,:,:) - db(:,:,:) * pdata%q(ivz,:,:,:)

! update the energy equation
!
          if (ien > 0) then

! calculate scalar product of velocity and magnetic field
!
            tmp(inx,inx,:,:,:) = sum(pdata%q(ivx:ivz,:,:,:)                    &
                                                  * pdata%q(ibx:ibz,:,:,:), 1)

! add the divergence potential source term to the energy equation, i.e.
!     d/dt E + ∇.F = - (∇.B) (v.B)
!
            du(ien,:,:,:) = du(ien,:,:,:) - db(:,:,:) * tmp(inx,inx,:,:,:)

          end if ! ien > 0

        end if ! glm_type > 1

      end if ! glmtype > 0

! proceed only if the resistivity coefficient is not zero
!
      if (resistivity > 0.0d+00) then

! calculate the Laplace operator of B, i.e. Δ(B)
!
        call laplace(dh(:), pdata%q(ibx,1:im,1:jm,1:km)                        &
                                                , tmp(inx,inx,1:im,1:jm,1:km))
        call laplace(dh(:), pdata%q(iby,1:im,1:jm,1:km)                        &
                                                , tmp(inx,iny,1:im,1:jm,1:km))
        call laplace(dh(:), pdata%q(ibz,1:im,1:jm,1:km)                        &
                                                , tmp(inx,inz,1:im,1:jm,1:km))

! multiply by the resistivity coefficient
!
        tmp(iny,inx:inz,1:im,1:jm,1:km) =                                      &
                                 resistivity * tmp(inx,inx:inz,1:im,1:jm,1:km)

! update magnetic field component increments
!
        du(ibx,1:im,1:jm,1:km) = du(ibx,1:im,1:jm,1:km)                        &
                                                 + tmp(iny,inx,1:im,1:jm,1:km)
        du(iby,1:im,1:jm,1:km) = du(iby,1:im,1:jm,1:km)                        &
                                                 + tmp(iny,iny,1:im,1:jm,1:km)
        du(ibz,1:im,1:jm,1:km) = du(ibz,1:im,1:jm,1:km)                        &
                                                 + tmp(iny,inz,1:im,1:jm,1:km)

! update energy equation
!
        if (ien > 0) then

! add the first resistive source term to the energy equation, i.e.
!     d/dt E + ∇.F = η B.[Δ(B)]
!
          du(ien,1:im,1:jm,1:km) = du(ien,1:im,1:jm,1:km)                      &
                 + (pdata%q(ibx,1:im,1:jm,1:km) * tmp(iny,inx,1:im,1:jm,1:km)  &
                  + pdata%q(iby,1:im,1:jm,1:km) * tmp(iny,iny,1:im,1:jm,1:km)  &
                  + pdata%q(ibz,1:im,1:jm,1:km) * tmp(iny,inz,1:im,1:jm,1:km))

! calculate current density J = ∇xB
!
          call curl(dh(:), pdata%q(ibx:ibz,1:im,1:jm,1:km)                     &
                                            , tmp(inz,inx:inz,1:im,1:jm,1:km))

! add the second resistive source term to the energy equation, i.e.
!     d/dt E + ∇.F = η J²
!
          du(ien,1:im,1:jm,1:km) = du(ien,1:im,1:jm,1:km)                      &
                    + resistivity * sum(tmp(inz,inx:inz,1:im,1:jm,1:km)**2, 2)

        end if ! energy equation present

      end if ! resistivity is not zero

    end if ! ibx > 0

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
