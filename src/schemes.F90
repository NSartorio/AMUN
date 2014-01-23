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
!! module: SCHEMES
!!
!!  The module provides and interface to numerical schemes, subroutines to
!!  calculate variable increment and one dimensional Riemann solvers.
!!
!!  If you implement a new set of equations, you have to add at least one
!!  corresponding update_flux subroutine, and one Riemann solver.
!!
!!******************************************************************************
!
module schemes

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
  integer            , save :: imi, imu, imf, imr
#endif /* PROFILE */

! pointer to the flux update procedure
!
  procedure(update_flux_hd_iso), pointer, save :: update_flux => null()

! pointer to the Riemann solver
!
  procedure(riemann_hd_iso_hll), pointer, save :: riemann     => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_schemes, finalize_schemes
  public :: update_flux, update_increment

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
! subroutine INITIALIZE_SCHEMES:
! -----------------------------
!
!   Subroutine initiate the module by setting module parameters and subroutine
!   pointers.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_schemes(verbose, iret)

! include external procedures and variables
!
    use equations , only : eqsys, eos
    use parameters, only : get_parameter_string

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64)      :: solver   = "HLL"
    character(len=255)     :: name_sol = ""
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('schemes:: initialization'  , imi)
    call set_timer('schemes:: increment update', imu)
    call set_timer('schemes:: flux update'     , imf)
    call set_timer('schemes:: Riemann solver'  , imr)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the Riemann solver
!
    call get_parameter_string("riemann_solver", solver)

! depending on the system of equations initialize the module variables
!
    select case(trim(eqsys))

!--- HYDRODYNAMICS ---
!
    case("hd", "HD", "hydro", "HYDRO", "hydrodynamic", "HYDRODYNAMIC")

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("iso", "ISO", "isothermal", "ISOTHERMAL")

! set pointers to subroutines
!
        update_flux => update_flux_hd_iso

! select the Riemann solver
!
        select case(trim(solver))


! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_hd_iso_hll

        end select

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! set pointers to subroutines
!
        update_flux => update_flux_hd_adi

! select the Riemann solver
!
        select case(trim(solver))

        case("hllc", "HLLC")

! set the solver name
!
          name_sol =  "HLLC"

! set pointers to subroutines
!
          riemann => riemann_hd_adi_hllc

! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_hd_adi_hll

        end select

      end select

!--- MAGNETOHYDRODYNAMICS ---
!
    case("mhd", "MHD", "magnetohydrodynamic", "MAGNETOHYDRODYNAMIC")

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("iso", "ISO", "isothermal", "ISOTHERMAL")

! set pointers to the subroutines
!
        update_flux => update_flux_mhd_iso

! select the Riemann solver
!
        select case(trim(solver))

        case ("hlld", "HLLD")

! set the solver name
!
          name_sol =  "HLLD"

! set the solver pointer
!
          riemann  => riemann_mhd_iso_hlld

        case ("hlldm", "hlld-m", "HLLDM", "HLLD-M")

! set the solver name
!
          name_sol =  "HLLD-M"

! set the solver pointer
!
          riemann  => riemann_mhd_iso_hlldm

! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_mhd_iso_hll

        end select

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! set pointers to subroutines
!
        update_flux => update_flux_mhd_adi

! select the Riemann solver
!
        select case(trim(solver))

        case ("hllc", "HLLC")

! set the solver name
!
          name_sol =  "HLLC"

! set the solver pointer
!
          riemann  => riemann_mhd_adi_hllc

        case ("hlld", "HLLD")

! set the solver name
!
          name_sol =  "HLLD"

! set the solver pointer
!
          riemann  => riemann_mhd_adi_hlld

! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_mhd_adi_hll

        end select

      end select

    end select

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "Riemann solver         =", trim(name_sol)

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_schemes
!
!===============================================================================
!
! subroutine FINALIZE_SCHEMES:
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
  subroutine finalize_schemes(iret)

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

! nullify procedure pointers
!
    nullify(update_flux)
    nullify(riemann)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_schemes
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
!
!===============================================================================
!
  subroutine update_increment(dh, f, du)

! include external variables
!
    use coordinates, only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations  , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), dimension(3)                , intent(in)    :: dh
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(in)    :: f
    real(kind=8), dimension(      nv,im,jm,km), intent(inout) :: du

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for increment update
!
    call start_timer(imu)
#endif /* PROFILE */

! reset the increment array du
!
    du(:,:,:,:) = 0.0d+00

! perform update along the X direction
!
    do i = ibl, ieu
      du(:,i,:,:) = du(:,i,:,:) - dh(1) * (f(1,:,i,:,:) - f(1,:,i-1,:,:))
    end do

! perform update along the Y direction
!
    do j = jbl, jeu
      du(:,:,j,:) = du(:,:,j,:) - dh(2) * (f(2,:,:,j,:) - f(2,:,:,j-1,:))
    end do

#if NDIMS == 3
! perform update along the Z direction
!
    do k = kbl, keu
      du(:,:,:,k) = du(:,:,:,k) - dh(3) * (f(3,:,:,:,k) - f(3,:,:,:,k-1))
    end do
#endif /* NDIMS == 3 */

#ifdef PROFILE
! stop accounting time for increment update
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_increment
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine UPDATE_FLUX_HD_ISO:
! -----------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     idir - direction along which the flux is calculated;
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!
!===============================================================================
!
  subroutine update_flux_hd_iso(idir, dx, q, f)

! include external variables
!
    use coordinates, only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, imx, imy, imz

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer                             , intent(in)    :: idir
    real(kind=8)                        , intent(in)    :: dx
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: q
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: f

! local variables
!
    integer                              :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im)       :: qx, fx
    real(kind=8), dimension(nv,jm)       :: qy, fy
#if NDIMS == 3
    real(kind=8), dimension(nv,km)       :: qz, fz
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! reset the flux array
!
    f(:,:,:,:) = 0.0d+00

! select the directional flux to compute
!
    select case(idir)
    case(1)

!  calculate the flux along the X-direction
!
      do k = kbl, keu
        do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
          qx(idn,1:im) = q(idn,1:im,j,k)
          qx(ivx,1:im) = q(ivx,1:im,j,k)
          qx(ivy,1:im) = q(ivy,1:im,j,k)
          qx(ivz,1:im) = q(ivz,1:im,j,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(im, dx, qx(1:nv,1:im), fx(1:nv,1:im))

! update the array of fluxes
!
          f(idn,1:im,j,k) = fx(idn,1:im)
          f(imx,1:im,j,k) = fx(imx,1:im)
          f(imy,1:im,j,k) = fx(imy,1:im)
          f(imz,1:im,j,k) = fx(imz,1:im)

        end do ! j = jbl, jeu
      end do ! k = kbl, keu

    case(2)

!  calculate the flux along the Y direction
!
      do k = kbl, keu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qy(idn,1:jm) = q(idn,i,1:jm,k)
          qy(ivx,1:jm) = q(ivy,i,1:jm,k)
          qy(ivy,1:jm) = q(ivz,i,1:jm,k)
          qy(ivz,1:jm) = q(ivx,i,1:jm,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(jm, dx, qy(1:nv,1:jm), fy(1:nv,1:jm))

! update the array of fluxes
!
          f(idn,i,1:jm,k) = fy(idn,1:jm)
          f(imx,i,1:jm,k) = fy(imz,1:jm)
          f(imy,i,1:jm,k) = fy(imx,1:jm)
          f(imz,i,1:jm,k) = fy(imy,1:jm)

        end do ! i = ibl, ieu
      end do ! k = kbl, keu

#if NDIMS == 3
    case(3)

!  calculate the flux along the Z direction
!
      do j = jbl, jeu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qz(idn,1:km) = q(idn,i,j,1:km)
          qz(ivx,1:km) = q(ivz,i,j,1:km)
          qz(ivy,1:km) = q(ivx,i,j,1:km)
          qz(ivz,1:km) = q(ivy,i,j,1:km)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(km, dx, qz(1:nv,1:km), fz(1:nv,1:km))

! update the array of fluxes
!
          f(idn,i,j,1:km) = fz(idn,1:km)
          f(imx,i,j,1:km) = fz(imy,1:km)
          f(imy,i,j,1:km) = fz(imz,1:km)
          f(imz,i,j,1:km) = fz(imx,1:km)

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_hd_iso
!
!===============================================================================
!
! subroutine UPDATE_FLUX_HD_ADI:
! -----------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     idir - direction along which the flux is calculated;
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!
!===============================================================================
!
  subroutine update_flux_hd_adi(idir, dx, q, f)

! include external variables
!
    use coordinates, only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer                             , intent(in)    :: idir
    real(kind=8)                        , intent(in)    :: dx
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: q
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: f

! local variables
!
    integer                              :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im)       :: qx, fx
    real(kind=8), dimension(nv,jm)       :: qy, fy
#if NDIMS == 3
    real(kind=8), dimension(nv,km)       :: qz, fz
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! reset the flux array
!
    f(:,:,:,:) = 0.0d+00

! select the directional flux to compute
!
    select case(idir)
    case(1)

!  calculate the flux along the X-direction
!
      do k = kbl, keu
        do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
          qx(idn,1:im) = q(idn,1:im,j,k)
          qx(ivx,1:im) = q(ivx,1:im,j,k)
          qx(ivy,1:im) = q(ivy,1:im,j,k)
          qx(ivz,1:im) = q(ivz,1:im,j,k)
          qx(ipr,1:im) = q(ipr,1:im,j,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(im, dx, qx(1:nv,1:im), fx(1:nv,1:im))

! update the array of fluxes
!
          f(idn,1:im,j,k) = fx(idn,1:im)
          f(imx,1:im,j,k) = fx(imx,1:im)
          f(imy,1:im,j,k) = fx(imy,1:im)
          f(imz,1:im,j,k) = fx(imz,1:im)
          f(ien,1:im,j,k) = fx(ien,1:im)

        end do ! j = jbl, jeu
      end do ! k = kbl, keu

    case(2)

!  calculate the flux along the Y direction
!
      do k = kbl, keu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qy(idn,1:jm) = q(idn,i,1:jm,k)
          qy(ivx,1:jm) = q(ivy,i,1:jm,k)
          qy(ivy,1:jm) = q(ivz,i,1:jm,k)
          qy(ivz,1:jm) = q(ivx,i,1:jm,k)
          qy(ipr,1:jm) = q(ipr,i,1:jm,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(jm, dx, qy(1:nv,1:jm), fy(1:nv,1:jm))

! update the array of fluxes
!
          f(idn,i,1:jm,k) = fy(idn,1:jm)
          f(imx,i,1:jm,k) = fy(imz,1:jm)
          f(imy,i,1:jm,k) = fy(imx,1:jm)
          f(imz,i,1:jm,k) = fy(imy,1:jm)
          f(ien,i,1:jm,k) = fy(ien,1:jm)

        end do ! i = ibl, ieu
      end do ! k = kbl, keu

#if NDIMS == 3
    case(3)

!  calculate the flux along the Z direction
!
      do j = jbl, jeu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qz(idn,1:km) = q(idn,i,j,1:km)
          qz(ivx,1:km) = q(ivz,i,j,1:km)
          qz(ivy,1:km) = q(ivx,i,j,1:km)
          qz(ivz,1:km) = q(ivy,i,j,1:km)
          qz(ipr,1:km) = q(ipr,i,j,1:km)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(km, dx, qz(1:nv,1:km), fz(1:nv,1:km))

! update the array of fluxes
!
          f(idn,i,j,1:km) = fz(idn,1:km)
          f(imx,i,j,1:km) = fz(imy,1:km)
          f(imy,i,j,1:km) = fz(imz,1:km)
          f(imz,i,j,1:km) = fz(imx,1:km)
          f(ien,i,j,1:km) = fz(ien,1:km)

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_hd_adi
!
!===============================================================================
!
! subroutine UPDATE_FLUX_MHD_ADI:
! ------------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     idir - direction along which the flux is calculated;
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!
!===============================================================================
!
  subroutine update_flux_mhd_iso(idir, dx, q, f)

! include external variables
!
    use coordinates, only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, imx, imy, imz
    use equations  , only : ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer                             , intent(in)    :: idir
    real(kind=8)                        , intent(in)    :: dx
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: q
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: f

! local variables
!
    integer                              :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im)       :: qx, fx
    real(kind=8), dimension(nv,jm)       :: qy, fy
#if NDIMS == 3
    real(kind=8), dimension(nv,km)       :: qz, fz
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! reset the flux array
!
    f(:,:,:,:) = 0.0d+00

! select the directional flux to compute
!
    select case(idir)
    case(1)

!  calculate the flux along the X-direction
!
      do k = kbl, keu
        do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
          qx(idn,1:im) = q(idn,1:im,j,k)
          qx(ivx,1:im) = q(ivx,1:im,j,k)
          qx(ivy,1:im) = q(ivy,1:im,j,k)
          qx(ivz,1:im) = q(ivz,1:im,j,k)
          qx(ibx,1:im) = q(ibx,1:im,j,k)
          qx(iby,1:im) = q(iby,1:im,j,k)
          qx(ibz,1:im) = q(ibz,1:im,j,k)
          qx(ibp,1:im) = q(ibp,1:im,j,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(im, dx, qx(1:nv,1:im), fx(1:nv,1:im))

! update the array of fluxes
!
          f(idn,1:im,j,k) = fx(idn,1:im)
          f(imx,1:im,j,k) = fx(imx,1:im)
          f(imy,1:im,j,k) = fx(imy,1:im)
          f(imz,1:im,j,k) = fx(imz,1:im)
          f(ibx,1:im,j,k) = fx(ibx,1:im)
          f(iby,1:im,j,k) = fx(iby,1:im)
          f(ibz,1:im,j,k) = fx(ibz,1:im)
          f(ibp,1:im,j,k) = fx(ibp,1:im)

        end do ! j = jbl, jeu
      end do ! k = kbl, keu

    case(2)

!  calculate the flux along the Y direction
!
      do k = kbl, keu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qy(idn,1:jm) = q(idn,i,1:jm,k)
          qy(ivx,1:jm) = q(ivy,i,1:jm,k)
          qy(ivy,1:jm) = q(ivz,i,1:jm,k)
          qy(ivz,1:jm) = q(ivx,i,1:jm,k)
          qy(ibx,1:jm) = q(iby,i,1:jm,k)
          qy(iby,1:jm) = q(ibz,i,1:jm,k)
          qy(ibz,1:jm) = q(ibx,i,1:jm,k)
          qy(ibp,1:jm) = q(ibp,i,1:jm,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(jm, dx, qy(1:nv,1:jm), fy(1:nv,1:jm))

! update the array of fluxes
!
          f(idn,i,1:jm,k) = fy(idn,1:jm)
          f(imx,i,1:jm,k) = fy(imz,1:jm)
          f(imy,i,1:jm,k) = fy(imx,1:jm)
          f(imz,i,1:jm,k) = fy(imy,1:jm)
          f(ibx,i,1:jm,k) = fy(ibz,1:jm)
          f(iby,i,1:jm,k) = fy(ibx,1:jm)
          f(ibz,i,1:jm,k) = fy(iby,1:jm)
          f(ibp,i,1:jm,k) = fy(ibp,1:jm)

        end do ! i = ibl, ieu
      end do ! k = kbl, keu

#if NDIMS == 3
    case(3)

!  calculate the flux along the Z direction
!
      do j = jbl, jeu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qz(idn,1:km) = q(idn,i,j,1:km)
          qz(ivx,1:km) = q(ivz,i,j,1:km)
          qz(ivy,1:km) = q(ivx,i,j,1:km)
          qz(ivz,1:km) = q(ivy,i,j,1:km)
          qz(ibx,1:km) = q(ibz,i,j,1:km)
          qz(iby,1:km) = q(ibx,i,j,1:km)
          qz(ibz,1:km) = q(iby,i,j,1:km)
          qz(ibp,1:km) = q(ibp,i,j,1:km)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(km, dx, qz(1:nv,1:km), fz(1:nv,1:km))

! update the array of fluxes
!
          f(idn,i,j,1:km) = fz(idn,1:km)
          f(imx,i,j,1:km) = fz(imy,1:km)
          f(imy,i,j,1:km) = fz(imz,1:km)
          f(imz,i,j,1:km) = fz(imx,1:km)
          f(ibx,i,j,1:km) = fz(iby,1:km)
          f(iby,i,j,1:km) = fz(ibz,1:km)
          f(ibz,i,j,1:km) = fz(ibx,1:km)
          f(ibp,i,j,1:km) = fz(ibp,1:km)

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_mhd_iso
!
!===============================================================================
!
! subroutine UPDATE_FLUX_MHD_ADI:
! ------------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     idir - direction along which the flux is calculated;
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!
!===============================================================================
!
  subroutine update_flux_mhd_adi(idir, dx, q, f)

! include external variables
!
    use coordinates, only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien
    use equations  , only : ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    integer                             , intent(in)    :: idir
    real(kind=8)                        , intent(in)    :: dx
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: q
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: f

! local variables
!
    integer                              :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im)       :: qx, fx
    real(kind=8), dimension(nv,jm)       :: qy, fy
#if NDIMS == 3
    real(kind=8), dimension(nv,km)       :: qz, fz
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! reset the flux array
!
    f(:,:,:,:) = 0.0d+00

! select the directional flux to compute
!
    select case(idir)
    case(1)

!  calculate the flux along the X-direction
!
      do k = kbl, keu
        do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
          qx(idn,1:im) = q(idn,1:im,j,k)
          qx(ivx,1:im) = q(ivx,1:im,j,k)
          qx(ivy,1:im) = q(ivy,1:im,j,k)
          qx(ivz,1:im) = q(ivz,1:im,j,k)
          qx(ibx,1:im) = q(ibx,1:im,j,k)
          qx(iby,1:im) = q(iby,1:im,j,k)
          qx(ibz,1:im) = q(ibz,1:im,j,k)
          qx(ibp,1:im) = q(ibp,1:im,j,k)
          qx(ipr,1:im) = q(ipr,1:im,j,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(im, dx, qx(1:nv,1:im), fx(1:nv,1:im))

! update the array of fluxes
!
          f(idn,1:im,j,k) = fx(idn,1:im)
          f(imx,1:im,j,k) = fx(imx,1:im)
          f(imy,1:im,j,k) = fx(imy,1:im)
          f(imz,1:im,j,k) = fx(imz,1:im)
          f(ibx,1:im,j,k) = fx(ibx,1:im)
          f(iby,1:im,j,k) = fx(iby,1:im)
          f(ibz,1:im,j,k) = fx(ibz,1:im)
          f(ibp,1:im,j,k) = fx(ibp,1:im)
          f(ien,1:im,j,k) = fx(ien,1:im)

        end do ! j = jbl, jeu
      end do ! k = kbl, keu

    case(2)

!  calculate the flux along the Y direction
!
      do k = kbl, keu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qy(idn,1:jm) = q(idn,i,1:jm,k)
          qy(ivx,1:jm) = q(ivy,i,1:jm,k)
          qy(ivy,1:jm) = q(ivz,i,1:jm,k)
          qy(ivz,1:jm) = q(ivx,i,1:jm,k)
          qy(ibx,1:jm) = q(iby,i,1:jm,k)
          qy(iby,1:jm) = q(ibz,i,1:jm,k)
          qy(ibz,1:jm) = q(ibx,i,1:jm,k)
          qy(ibp,1:jm) = q(ibp,i,1:jm,k)
          qy(ipr,1:jm) = q(ipr,i,1:jm,k)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(jm, dx, qy(1:nv,1:jm), fy(1:nv,1:jm))

! update the array of fluxes
!
          f(idn,i,1:jm,k) = fy(idn,1:jm)
          f(imx,i,1:jm,k) = fy(imz,1:jm)
          f(imy,i,1:jm,k) = fy(imx,1:jm)
          f(imz,i,1:jm,k) = fy(imy,1:jm)
          f(ibx,i,1:jm,k) = fy(ibz,1:jm)
          f(iby,i,1:jm,k) = fy(ibx,1:jm)
          f(ibz,i,1:jm,k) = fy(iby,1:jm)
          f(ibp,i,1:jm,k) = fy(ibp,1:jm)
          f(ien,i,1:jm,k) = fy(ien,1:jm)

        end do ! i = ibl, ieu
      end do ! k = kbl, keu

#if NDIMS == 3
    case(3)

!  calculate the flux along the Z direction
!
      do j = jbl, jeu
        do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
          qz(idn,1:km) = q(idn,i,j,1:km)
          qz(ivx,1:km) = q(ivz,i,j,1:km)
          qz(ivy,1:km) = q(ivx,i,j,1:km)
          qz(ivz,1:km) = q(ivy,i,j,1:km)
          qz(ibx,1:km) = q(ibz,i,j,1:km)
          qz(iby,1:km) = q(ibx,i,j,1:km)
          qz(ibz,1:km) = q(iby,i,j,1:km)
          qz(ibp,1:km) = q(ibp,i,j,1:km)
          qz(ipr,1:km) = q(ipr,i,j,1:km)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
          call riemann(km, dx, qz(1:nv,1:km), fz(1:nv,1:km))

! update the array of fluxes
!
          f(idn,i,j,1:km) = fz(idn,1:km)
          f(imx,i,j,1:km) = fz(imy,1:km)
          f(imy,i,j,1:km) = fz(imz,1:km)
          f(imz,i,j,1:km) = fz(imx,1:km)
          f(ibx,i,j,1:km) = fz(iby,1:km)
          f(iby,i,j,1:km) = fz(ibz,1:km)
          f(ibz,i,j,1:km) = fz(ibx,1:km)
          f(ibp,i,j,1:km) = fz(ibp,1:km)
          f(ien,i,j,1:km) = fz(ien,1:km)

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

    end select

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_mhd_adi
!
!===============================================================================
!
! subroutine RIEMANN_HD_ISO_HLL:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann_hd_iso_hll(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! calculate the fluxes for the intermediate state
!
        f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_hd_iso_hll
!
!===============================================================================
!
! subroutine RIEMANN_HD_ADI_HLL:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann_hd_adi_hll(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ipr, ivx
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! calculate the fluxes for the intermediate state
!
        f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_hd_adi_hll
!
!===============================================================================
!
! subroutine RIEMANN_HD_ADI_HLLC:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using the HLLC method,
!   by Toro.  In the HLLC method the tangential components of the velocity are
!   discontinuous actoss the contact dicontinuity.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Toro, E. F., Spruce, M., & Speares, W.
!         "Restoration of the contact surface in the HLL-Riemann solver",
!         Shock Waves, 1994, Volume 4, Issue 1, pp. 25-34
!
!===============================================================================
!
  subroutine riemann_hd_adi_hllc(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, ivy, ivz, ipr, imx, imy, imz, ien
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, sm
    real(kind=8)                  :: srml, slmm, srmm
    real(kind=8)                  :: dn, pr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, ui
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml  = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:) = sl * ul(:,i) - fl(:,i)
        wr(:) = sr * ur(:,i) - fr(:,i)

! the speed of contact discontinuity
!
        dn    =  wr(idn) - wl(idn)
        sm    = (wr(imx) - wl(imx)) / dn

! calculate the pressure if the intermediate state
!
        pr    = 0.5d+00 * ((wr(idn) * sm - wr(imx)) + (wl(idn) * sm - wl(imx)))

! separate intermediate states depending on the sign of the advection speed
!
        if (sm > 0.0d+00) then ! sm > 0

! the left speed difference
!
          slmm    =  sl - sm

! conservative variables for the left intermediate state
!
          ui(idn) =  wl(idn) / slmm
          ui(imx) =  ui(idn) * sm
          ui(imy) =  ui(idn) * ql(ivy,i)
          ui(imz) =  ui(idn) * ql(ivz,i)
          ui(ien) = (wl(ien) + sm * pr) / slmm

! the left intermediate flux
!
          f(:,i)  = sl * ui(:) - wl(:)

        else if (sm < 0.0d+00) then ! sm < 0

! the right speed difference
!
          srmm    =  sr - sm

! conservative variables for the right intermediate state
!
          ui(idn) =  wr(idn) / srmm
          ui(imx) =  ui(idn) * sm
          ui(imy) =  ui(idn) * qr(ivy,i)
          ui(imz) =  ui(idn) * qr(ivz,i)
          ui(ien) = (wr(ien) + sm * pr) / srmm

! the right intermediate flux
!
          f(:,i)  = sr * ui(:) - wr(:)

        else ! sm = 0

! intermediate flux is constant across the contact discontinuity and all except
! the parallel momentum flux are zero
!
          f(idn,i) =   0.0d+00
          f(imx,i) = - 0.5d+00 * (wl(imx) + wr(imx))
          f(imy,i) =   0.0d+00
          f(imz,i) =   0.0d+00
          f(ien,i) =   0.0d+00

        end if ! sm = 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_hd_adi_hllc
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ISO_HLL:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann_mhd_iso_hll(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, ibx, ibp
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! calculate the fluxes for the intermediate state
!
        f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_iso_hll
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ISO_HLLD:
! -------------------------------
!
!   Subroutine solves one dimensional Riemann problem using the isothermal HLLD
!   method with correct state averaging and degeneracies treatement.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Kowal, G.,
!         "An isothermal MHD Riemann solver with correct state averaging and
!          degeneracies treatement",
!         Journal of Computational Physics, in preparation
!
!===============================================================================
!
  subroutine riemann_mhd_iso_hlld(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, imx, imy, imz, ibx, iby, ibz, ibp
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, sm, sml, smr, srml, slmm, srmm
    real(kind=8)                  :: bx, b2, dn, dnl, dnr, dvl, dvr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! the advection speed in the intermediate states
!
        dn =  wr(idn) - wl(idn)
        sm = (wr(imx) - wl(imx)) / dn

! square of Bₓ, i.e. Bₓ²
!
        bx = ql(ibx,i)
        b2 = ql(ibx,i) * qr(ibx,i)

! speed differences
!
        slmm = sl - sm
        srmm = sr - sm

! left and right state densities
!
        dnl = wl(idn) / slmm
        dnr = wr(idn) / srmm

! if there is an Alvén wave, apply the full HLLD solver, otherwise revert to
! the HLL one
!
        if (b2 > 0.0d+00) then ! Bₓ² > 0

! left and right Alfvén speeds
!
          sml = sm - sqrt(b2 / dnl)
          smr = sm + sqrt(b2 / dnr)

! calculate denominators in order to check degeneracy
!
          dvl = slmm * wl(idn) - b2
          dvr = srmm * wr(idn) - b2

! check degeneracy Sl* -> Sl or Sr* -> Sr
!
          if (min(dvl, dvr) > 0.0d+00) then ! no degeneracy

! choose the correct state depending on the speed signs
!
            if (sml >= 0.0d+00) then ! sl* ≥ 0

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl *  sm
              ui(imy) = dnl * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dnl * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! the outer left intermediate flux
!
              f(:,i)  = sl * ui(:) - wl(:)

            else if (smr <= 0.0d+00) then ! sr* ≤ 0

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr *  sm
              ui(imy) = dnr * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dnr * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! the outer right intermediate flux
!
              f(:,i)  = sr * ui(:) - wr(:)

            else ! sl* < 0 < sr*

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl *  sm
              ui(imy) = dnl * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dnl * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! vector of the left-going Alfvén wave
!
              wcl(:)  = (sml - sl) * ui(:) + wl(:)

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr *  sm
              ui(imy) = dnr * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dnr * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! vector of the right-going Alfvén wave
!
              wcr(:)  = (smr - sr) * ui(:) + wr(:)

! the flux corresponding to the middle intermediate state
!
              f(:,i)  = (sml * wcr(:) - smr * wcl(:)) / (smr - sml)

            end if ! sl* < 0 < sr*

          else ! one degeneracy

! separate degeneracies
!
            if (dvl > 0.0d+00) then ! sr* > sr

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl *  sm
              ui(imy) = dnl * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dnl * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! choose the correct state depending on the speed signs
!
              if (sml >= 0.0d+00) then ! sl* ≥ 0

! the outer left intermediate flux
!
                f(:,i)  = sl * ui(:) - wl(:)

              else ! sl* < 0

! vector of the left-going Alfvén wave
!
                wcl(:)  = (sml - sl) * ui(:) + wl(:)

! the flux corresponding to the middle intermediate state
!
                f(:,i)  = (sml * wr(:) - sr * wcl(:)) / (sr - sml)

              end if ! sl* < 0

            else if (dvr > 0.0d+00) then ! sl* < sl

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr *  sm
              ui(imy) = dnr * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dnr * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! choose the correct state depending on the speed signs
!
              if (smr <= 0.0d+00) then ! sr* ≤ 0

! the outer right intermediate flux
!
                f(:,i)  = sr * ui(:) - wr(:)

              else ! sr* > 0

! vector of the right-going Alfvén wave
!
                wcr(:)  = (smr - sr) * ui(:) + wr(:)

! the flux corresponding to the middle intermediate state
!
                f(:,i)  = (sl * wcr(:) - smr * wl(:)) / (smr - sl)

              end if ! sr* > 0

            else ! sl* < sl & sr* > sr

! both outer states are degenerate so revert to the HLL flux
!
              f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

            end if ! sl* < sl & sr* > sr

          end if ! one degeneracy

        else ! Bₓ² = 0

! in the vase of vanishing Bₓ there is no Alfvén wave, density is constant, and
! still we can have a discontinuity in the perpendicular components
!
          dn   = dn / srml

! take the right flux depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! conservative variables for the left intermediate state
!
            ui(idn) = dn
            ui(imx) = dn *  sm
            ui(imy) = wl(imy) / slmm
            ui(imz) = wl(imz) / slmm
            ui(ibx) = 0.0d+00
            ui(iby) = wl(iby) / slmm
            ui(ibz) = wl(ibz) / slmm
            ui(ibp) = ql(ibp,i)

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the right intermediate state
!
            ui(idn) = dn
            ui(imx) = dn *  sm
            ui(imy) = wr(imy) / srmm
            ui(imz) = wr(imz) / srmm
            ui(ibx) = 0.0d+00
            ui(iby) = wr(iby) / srmm
            ui(ibz) = wr(ibz) / srmm
            ui(ibp) = qr(ibp,i)

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! the intermediate flux; since the advection speed is zero, perpendicular
! components do not change
!
            f(idn,i) = (sl * wr(idn) - sr * wl(idn)) / srml
            f(imx,i) = (sl * wr(imx) - sr * wl(imx)) / srml
            f(imy,i) = 0.0d+00
            f(imz,i) = 0.0d+00
            f(ibx,i) = (sl * wr(ibx) - sr * wl(ibx)) / srml
            f(iby,i) = 0.0d+00
            f(ibz,i) = 0.0d+00
            f(ibp,i) = (sl * wr(ibp) - sr * wl(ibp)) / srml

          end if ! sm = 0

        end if ! Bₓ² = 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_iso_hlld
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ISO_HLLDM:
! --------------------------------
!
!   Subroutine solves one dimensional Riemann problem using the isothermal HLLD
!   method described by Mignone.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Mignone, A.,
!         "A simple and accurate Riemann solver for isothermal MHD",
!         Journal of Computational Physics, 2007, 225, pp. 1427-1441
!
!===============================================================================
!
  subroutine riemann_mhd_iso_hlldm(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, imx, imy, imz, ibx, iby, ibz, ibp
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, sm, sml, smr, srml, slmm, srmm
    real(kind=8)                  :: bx, b2, dn, dnl, dnr, dvl, dvr, ca

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! the advection speed in the intermediate states
!
        dn =  wr(idn) - wl(idn)
        sm = (wr(imx) - wl(imx)) / dn

! square of Bₓ, i.e. Bₓ²
!
        bx = ql(ibx,i)
        b2 = ql(ibx,i) * qr(ibx,i)

! speed differences
!
        slmm = sl - sm
        srmm = sr - sm

! calculate density of the intermediate states
!
        dn = dn / srml

! if there is an Alvén wave, apply the full HLLD solver, otherwise revert to
! the HLL one
!
        if (b2 > 0.0d+00) then ! Bₓ² > 0

! left and right Alfvén speeds
!
          ca   = sqrt(b2 / dn)
          sml  = sm - ca
          smr  = sm + ca

! calculate denominators in order to check degeneracy
!
          dvl = slmm * slmm * dn - b2
          dvr = srmm * srmm * dn - b2

! check degeneracy Sl* -> Sl or Sr* -> Sr
!
          if (min(dvl, dvr) > 0.0d+00) then ! no degeneracy

! choose the correct state depending on the speed signs
!
            if (sml >= 0.0d+00) then ! sl* ≥ 0

! conservative variables for the outer left intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dn * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (slmm * dn * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (slmm * dn * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! the outer left intermediate flux
!
              f(:,i)  = sl * ui(:) - wl(:)

            else if (smr <= 0.0d+00) then ! sr* ≤ 0

! conservative variables for the outer right intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dn * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (srmm * dn * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (srmm * dn * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! the outer right intermediate flux
!
              f(:,i)  = sr * ui(:) - wr(:)

            else ! sl* < 0 < sr*

! conservative variables for the outer left intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dn * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (slmm * dn * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (slmm * dn * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! vector of the left-going Alfvén wave
!
              wcl(:)  = (sml - sl) * ui(:) + wl(:)

! conservative variables for the outer right intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dn * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (srmm * dn * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (srmm * dn * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! vector of the right-going Alfvén wave
!
              wcr(:)  = (smr - sr) * ui(:) + wr(:)

! the flux corresponding to the middle intermediate state
!
              f(:,i)  = (sml * wcr(:) - smr * wcl(:)) / (smr - sml)

            end if ! sl* < 0 < sr*

          else ! one degeneracy

! separate degeneracies
!
            if (dvl > 0.0d+00) then ! sr* > sr

! conservative variables for the outer left intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (slmm * wl(imy) - bx * wl(iby)) / dvl
              ui(imz) = dn * (slmm * wl(imz) - bx * wl(ibz)) / dvl
              ui(ibx) = bx
              ui(iby) = (slmm * dn * wl(iby) - bx * wl(imy)) / dvl
              ui(ibz) = (slmm * dn * wl(ibz) - bx * wl(imz)) / dvl
              ui(ibp) = ql(ibp,i)

! choose the correct state depending on the speed signs
!
              if (sml >= 0.0d+00) then ! sl* ≥ 0

! the outer left intermediate flux
!
                f(:,i)  = sl * ui(:) - wl(:)

              else ! sl* < 0

! vector of the left-going Alfvén wave
!
                wcl(:)  = (sml - sl) * ui(:) + wl(:)

! the flux corresponding to the middle intermediate state
!
                f(:,i)  = (sml * wr(:) - sr * wcl(:)) / (sr - sml)

              end if ! sl* < 0

            else if (dvr > 0.0d+00) then ! sl* < sl

! conservative variables for the outer right intermediate state
!
              ui(idn) = dn
              ui(imx) = dn *  sm
              ui(imy) = dn * (srmm * wr(imy) - bx * wr(iby)) / dvr
              ui(imz) = dn * (srmm * wr(imz) - bx * wr(ibz)) / dvr
              ui(ibx) = bx
              ui(iby) = (srmm * dn * wr(iby) - bx * wr(imy)) / dvr
              ui(ibz) = (srmm * dn * wr(ibz) - bx * wr(imz)) / dvr
              ui(ibp) = qr(ibp,i)

! choose the correct state depending on the speed signs
!
              if (smr <= 0.0d+00) then ! sr* ≤ 0

! the outer right intermediate flux
!
                f(:,i)  = sr * ui(:) - wr(:)

              else ! sr* > 0

! vector of the right-going Alfvén wave
!
                wcr(:)  = (smr - sr) * ui(:) + wr(:)

! the flux corresponding to the middle intermediate state
!
                f(:,i)  = (sl * wcr(:) - smr * wl(:)) / (smr - sl)

              end if ! sr* > 0

            else ! sl* < sl & sr* > sr

! both outer states are degenerate so revert to the HLL flux
!
              f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

            end if ! sl* < sl & sr* > sr

          end if ! one degeneracy

        else ! Bₓ² = 0

! in the vase of vanishing Bₓ there is no Alfvén wave, density is constant, and
! still we can have a discontinuity in the perpendicular components
!
! take the right flux depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! conservative variables for the left intermediate state
!
            ui(idn) = dn
            ui(imx) = dn *  sm
            ui(imy) = wl(imy) / slmm
            ui(imz) = wl(imz) / slmm
            ui(ibx) = 0.0d+00
            ui(iby) = wl(iby) / slmm
            ui(ibz) = wl(ibz) / slmm
            ui(ibp) = ql(ibp,i)

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the right intermediate state
!
            ui(idn) = dn
            ui(imx) = dn *  sm
            ui(imy) = wr(imy) / srmm
            ui(imz) = wr(imz) / srmm
            ui(ibx) = 0.0d+00
            ui(iby) = wr(iby) / srmm
            ui(ibz) = wr(ibz) / srmm
            ui(ibp) = qr(ibp,i)

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! the intermediate flux; since the advection speed is zero, perpendicular
! components do not change
!
            f(idn,i) = (sl * wr(idn) - sr * wl(idn)) / srml
            f(imx,i) = (sl * wr(imx) - sr * wl(imx)) / srml
            f(imy,i) = 0.0d+00
            f(imz,i) = 0.0d+00
            f(ibx,i) = (sl * wr(ibx) - sr * wl(ibx)) / srml
            f(iby,i) = 0.0d+00
            f(ibz,i) = 0.0d+00
            f(ibp,i) = (sl * wr(ibp) - sr * wl(ibp)) / srml

          end if ! sm = 0

        end if ! Bₓ² = 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_iso_hlldm
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ADI_HLL:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Harten, A., Lax, P. D. & Van Leer, B.
!         "On Upstream Differencing and Godunov-Type Schemes for Hyperbolic
!          Conservation Laws",
!         SIAM Review, 1983, Volume 25, Number 1, pp. 35-61
!
!===============================================================================
!
  subroutine riemann_mhd_adi_hll(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ipr, ivx, ibx, ibp
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! calculate speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! calculate the fluxes for the intermediate state
!
        f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_adi_hll
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ADI_HLLC:
! -------------------------------
!
!   Subroutine solves one dimensional Riemann problem using the HLLC
!   method by Gurski or Li.  The HLLC and HLLC-C differ by definitions of
!   the tangential components of the velocity and magnetic field.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Gurski, K. F.,
!         "An HLLC-Type Approximate Riemann Solver for Ideal
!          Magnetohydrodynamics",
!         SIAM Journal on Scientific Computing, 2004, Volume 25, Issue 6,
!         pp. 2165–2187
!     [2] Li, S.,
!         "An HLLC Riemann solver for magneto-hydrodynamics",
!         Journal of Computational Physics, 2005, Volume 203, Issue 1,
!         pp. 344-357
!
!===============================================================================
!
  subroutine riemann_mhd_adi_hllc(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp, ipr
    use equations     , only : imx, imy, imz, ien
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, sm, srml, slmm, srmm
    real(kind=8)                  :: dn, bx, b2, pt, vy, vz, by, bz, vb

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, ui
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLC flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! the speed of contact discontinuity
!
        dn =  wr(idn) - wl(idn)
        sm = (wr(imx) - wl(imx)) / dn

! square of Bx, i.e. Bₓ²
!
        bx = ql(ibx,i)
        b2 = ql(ibx,i) * qr(ibx,i)

! separate the cases when Bₓ = 0 and Bₓ ≠ 0
!
        if (b2 > 0.0d+00) then

! the total pressure, constant across the contact discontinuity
!
          pt = 0.5d+00 * ((wl(idn) * sm - wl(imx))                             &
                        + (wr(idn) * sm - wr(imx))) + b2

! constant intermediate state tangential components of velocity and
! magnetic field
!
          vy = (wr(imy) - wl(imy)) / dn
          vz = (wr(imz) - wl(imz)) / dn
          by = (wr(iby) - wl(iby)) / srml
          bz = (wr(ibz) - wl(ibz)) / srml

! scalar product of velocity and magnetic field for the intermediate states
!
          vb = sm * bx + vy * by + vz * bz

! separate intermediate states depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! the left speed difference
!
            slmm    =  sl - sm

! conservative variables for the left intermediate state
!
            ui(idn) =  wl(idn) / slmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * vy
            ui(imz) =  ui(idn) * vz
            ui(ibx) =  bx
            ui(iby) =  by
            ui(ibz) =  bz
            ui(ibp) =  ql(ibp,i)
            ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! the right speed difference
!
            srmm    =  sr - sm

! conservative variables for the right intermediate state
!
            ui(idn) =  wr(idn) / srmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * vy
            ui(imz) =  ui(idn) * vz
            ui(ibx) =  bx
            ui(iby) =  by
            ui(ibz) =  bz
            ui(ibp) =  qr(ibp,i)
            ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! when Sₘ = 0 all variables are continuous, therefore the flux reduces
! to the HLL one
!
            f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

          end if ! sm = 0

        else ! Bₓ = 0

! the total pressure, constant across the contact discontinuity
!
          pt = 0.5d+00 * ((wl(idn) * sm - wl(imx))                             &
                        + (wr(idn) * sm - wr(imx)))

! separate intermediate states depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! the left speed difference
!
            slmm    =  sl - sm

! conservative variables for the left intermediate state
!
            ui(idn) =  wl(idn) / slmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * ql(ivy,i)
            ui(imz) =  ui(idn) * ql(ivz,i)
            ui(ibx) =  0.0d+00
            ui(iby) =  wl(iby) / slmm
            ui(ibz) =  wl(ibz) / slmm
            ui(ibp) =  ql(ibp,i)
            ui(ien) = (wl(ien) + sm * pt) / slmm

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! the right speed difference
!
            srmm    =  sr - sm

! conservative variables for the right intermediate state
!
            ui(idn) =  wr(idn) / srmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * qr(ivy,i)
            ui(imz) =  ui(idn) * qr(ivz,i)
            ui(ibx) =  0.0d+00
            ui(iby) =  wr(iby) / srmm
            ui(ibz) =  wr(ibz) / srmm
            ui(ibp) =  qr(ibp,i)
            ui(ien) = (wr(ien) + sm * pt) / srmm

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! intermediate flux is constant across the contact discontinuity and all except
! the parallel momentum flux are zero
!
            f(idn,i) =   0.0d+00
            f(imx,i) = - 0.5d+00 * (wl(imx) + wr(imx))
            f(imy,i) =   0.0d+00
            f(imz,i) =   0.0d+00
            f(ibx,i) = fl(ibx,i)
            f(iby,i) =   0.0d+00
            f(ibz,i) =   0.0d+00
            f(ibp,i) = fl(ibp,i)
            f(ien,i) =   0.0d+00

          end if ! sm = 0

        end if ! Bₓ = 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_adi_hllc
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ADI_HLLD:
! -------------------------------
!
!   Subroutine solves one dimensional Riemann problem using the adiabatic HLLD
!   method described by Miyoshi & Kusano.
!
!   Arguments:
!
!     n - the length of input vectors;
!     h - the spatial step;
!     q - the input array of primitive variables;
!     f - the output array of fluxes;
!
!   References:
!
!     [1] Miyoshi, T. & Kusano, K.,
!         "A multi-state HLL approximate Riemann solver for ideal
!          magnetohydrodynamics",
!         Journal of Computational Physics, 2005, 208, pp. 315-344
!
!===============================================================================
!
  subroutine riemann_mhd_adi_hlld(n, h, q, f)

! include external procedures
!
    use equations     , only : nv
    use equations     , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp, ipr
    use equations     , only : imx, imy, imz, ien
    use equations     , only : cmax
    use equations     , only : prim2cons, fluxspeed
    use interpolations, only : reconstruct, fix_positivity

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)  :: n
    real(kind=8)                 , intent(in)  :: h
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: f

! local variables
!
    integer                       :: p, i
    real(kind=8)                  :: sl, sr, sm, srml, slmm, srmm
    real(kind=8)                  :: dn, bx, b2, pt, vy, vz, by, bz, vb, dv
    real(kind=8)                  :: dnl, dnr, cal, car, sml, smr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ql, qr, ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: cl, cr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the left and right states of primitive variables
!
    do p = 1, nv
      call reconstruct(n, h, q(p,:), ql(p,:), qr(p,:))
    end do

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    cl(:) = 0.5d+00 * ((qr(ibx,:) + ql(ibx,:)) - (qr(ibp,:) - ql(ibp,:)) / cmax)
    cr(:) = 0.5d+00 * ((qr(ibp,:) + ql(ibp,:)) - (qr(ibx,:) - ql(ibx,:)) * cmax)
    ql(ibx,:) = cl(:)
    qr(ibx,:) = cl(:)
    ql(ibp,:) = cr(:)
    qr(ibp,:) = cr(:)

! check if the reconstruction doesn't give the negative density or pressure,
! if so, correct the states
!
    call fix_positivity(n, q(idn,:), ql(idn,:), qr(idn,:))
    call fix_positivity(n, q(ipr,:), ql(ipr,:), qr(ipr,:))

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), cl(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), cr(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(ql(ivx,i) - cl(i), qr(ivx,i) - cr(i))
      sr = max(ql(ivx,i) + cl(i), qr(ivx,i) + cr(i))

! calculate the HLLC flux
!
      if (sl >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else ! sl < 0 < sr

! speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(:)  = sl * ul(:,i) - fl(:,i)
        wr(:)  = sr * ur(:,i) - fr(:,i)

! the speed of contact discontinuity
!
        dn =  wr(idn) - wl(idn)
        sm = (wr(imx) - wl(imx)) / dn

! square of Bₓ, i.e. Bₓ²
!
        bx = ql(ibx,i)
        b2 = ql(ibx,i) * qr(ibx,i)

! if there is an Alvén wave, apply the full HLLD solver, otherwise revert to
! the HLLC solver
!
        if (b2 > 0.0d+00) then ! Bₓ² > 0

! the total pressure, constant across the contact discontinuity and Alfvén waves
!
          pt = (wl(idn) * wr(imx) - wr(idn) * wl(imx)) / dn + b2

! speed differences
!
          slmm = sl - sm
          srmm = sr - sm

! left and right state densities
!
          dnl = wl(idn) / slmm
          dnr = wr(idn) / srmm

! left and right Alfvén speeds
!
          cal = - sqrt(b2 / dnl)
          car =   sqrt(b2 / dnr)
          sml = sm + cal
          smr = sm + car

! check degeneracy Sl* -> Sl or Sr* -> Sr
!
          if (min(sml - sl, sr - smr) > 0.0d+00) then ! no degeneracy

! choose the correct state depending on the speed signs
!
            if (sml >= 0.0d+00) then ! sl* ≥ 0

! primitive variables for the outer left intermediate state
!
              dv      =     slmm * wl(idn) - b2
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dv
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dv
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dv
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl * sm
              ui(imy) = dnl * vy
              ui(imz) = dnl * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = ql(ibp,i)
              ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! the outer left intermediate flux
!
              f(:,i)  = sl * ui(:) - wl(:)

          else if (smr <= 0.0d+00) then ! sr* ≤ 0

! primitive variables for the outer right intermediate state
!
              dv      =     srmm * wr(idn) - b2
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dv
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dv
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dv
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr * sm
              ui(imy) = dnr * vy
              ui(imz) = dnr * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = qr(ibp,i)
              ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! the outer right intermediate flux
!
              f(:,i)  = sr * ui(:) - wr(:)

            else ! sl* < 0 < sr*

! primitive variables for the outer left intermediate state
!
              dv      =     slmm * wl(idn) - b2
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dv
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dv
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dv
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl * sm
              ui(imy) = dnl * vy
              ui(imz) = dnl * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = ql(ibp,i)
              ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! vector of the left-going Alfvén wave
!
              wcl(:)  = (sml - sl) * ui(:) + wl(:)

! primitive variables for the outer right intermediate state
!
              dv      =     srmm * wr(idn) - b2
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dv
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dv
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dv
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr * sm
              ui(imy) = dnr * vy
              ui(imz) = dnr * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = qr(ibp,i)
              ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! vector of the right-going Alfvén wave
!
              wcr(:)  = (smr - sr) * ui(:) + wr(:)

! prepare constant primitive variables of the intermediate states
!
              dv      = car * dnr - cal * dnl
              vy      = (wcr(imy) - wcl(imy)) / dv
              vz      = (wcr(imz) - wcl(imz)) / dv
              dv      = car       - cal
              by      = (wcr(iby) - wcl(iby)) / dv
              bz      = (wcr(ibz) - wcl(ibz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! choose the correct state depending on the sign of contact discontinuity
! advection speed
!
              if (sm > 0.0d+00) then ! sm > 0

! conservative variables for the inmost left intermediate state
!
                ui(idn) = dnl
                ui(imx) = dnl * sm
                ui(imy) = dnl * vy
                ui(imz) = dnl * vz
                ui(ibx) = bx
                ui(iby) = by
                ui(ibz) = bz
                ui(ibp) = ql(ibp,i)
                ui(ien) = (wcl(ien) + sm * pt - bx * vb) / cal

! the inmost left intermediate flux
!
                f(:,i)  = sml * ui(:) - wcl(:)

              else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the inmost right intermediate state
!
                ui(idn) = dnr
                ui(imx) = dnr * sm
                ui(imy) = dnr * vy
                ui(imz) = dnr * vz
                ui(ibx) = bx
                ui(iby) = by
                ui(ibz) = bz
                ui(ibp) = qr(ibp,i)
                ui(ien) = (wcr(ien) + sm * pt - bx * vb) / car

! the inmost right intermediate flux
!
                f(:,i)  = smr * ui(:) - wcr(:)

              else ! sm = 0

! in the case when Sₘ = 0 and Bₓ² > 0, all variables are continuous, therefore
! the flux can be averaged from the Alfvén waves using the simple HLL formula
!
                f(:,i) = (sml * wcr(:) - smr * wcl(:)) / (smr - sml)

              end if ! sm = 0

            end if ! sl* < 0 < sr*

          else ! one degeneracy

! separate degeneracies
!
            if ((sml - sl) > 0.0d+00) then ! sr* > sr

! primitive variables for the outer left intermediate state
!
              dv      =     slmm * wl(idn) - b2
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dv
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dv
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dv
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer left intermediate state
!
              ui(idn) = dnl
              ui(imx) = dnl * sm
              ui(imy) = dnl * vy
              ui(imz) = dnl * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = ql(ibp,i)
              ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! choose the correct state depending on the speed signs
!
              if (sml >= 0.0d+00) then ! sl* ≥ 0

! the outer left intermediate flux
!
                f(:,i)  = sl * ui(:) - wl(:)

              else ! sl* < 0

! vector of the left-going Alfvén wave
!
                wcl(:)  = (sml - sl) * ui(:) + wl(:)

! primitive variables for the inner left intermediate state
!
                dv      = srmm * dnr - cal * dnl
                vy      = (wr(imy) - wcl(imy)) / dv
                vz      = (wr(imz) - wcl(imz)) / dv
                dv      = sr - sml
                by      = (wr(iby) - wcl(iby)) / dv
                bz      = (wr(ibz) - wcl(ibz)) / dv
                vb      = sm * bx + vy * by + vz * bz

! choose the correct state depending on the sign of contact discontinuity
! advection speed
!
                if (sm >= 0.0d+00) then ! sm ≥ 0

! conservative variables for the inner left intermediate state
!
                  ui(idn) = dnl
                  ui(imx) = dnl * sm
                  ui(imy) = dnl * vy
                  ui(imz) = dnl * vz
                  ui(ibx) = bx
                  ui(iby) = by
                  ui(ibz) = bz
                  ui(ibp) = ql(ibp,i)
                  ui(ien) = (wcl(ien) + sm * pt - bx * vb) / cal

! the inner left intermediate flux
!
                  f(:,i)  = sml * ui(:) - wcl(:)

                else ! sm < 0

! vector of the left-going Alfvén wave
!
                  wcl(:)  = (sm - sml) * ui(:) + wcl(:)

! calculate the average flux over the right inner intermediate state
!
                  f(:,i)  = (sm * wr(:) - sr * wcl(:)) / (sr - sm)

                end if ! sm < 0

              end if ! sl* < 0

            else if ((sr - smr) > 0.0d+00) then ! sl* < sl

! primitive variables for the outer right intermediate state
!
              dv      =     srmm * wr(idn) - b2
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dv
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dv
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dv
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dv
              vb      = sm * bx + vy * by + vz * bz

! conservative variables for the outer right intermediate state
!
              ui(idn) = dnr
              ui(imx) = dnr * sm
              ui(imy) = dnr * vy
              ui(imz) = dnr * vz
              ui(ibx) = bx
              ui(iby) = by
              ui(ibz) = bz
              ui(ibp) = qr(ibp,i)
              ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! choose the correct state depending on the speed signs
!
              if (smr <= 0.0d+00) then ! sr* ≤ 0

! the outer right intermediate flux
!
                f(:,i)  = sr * ui(:) - wr(:)

              else ! sr* > 0

! vector of the right-going Alfvén wave
!
                wcr(:)  = (smr - sr) * ui(:) + wr(:)

! primitive variables for the inner right intermediate state
!
                dv      = slmm * dnl - car * dnr
                vy      = (wl(imy) - wcr(imy)) / dv
                vz      = (wl(imz) - wcr(imz)) / dv
                dv      = sl - smr
                by      = (wl(iby) - wcr(iby)) / dv
                bz      = (wl(ibz) - wcr(ibz)) / dv
                vb      = sm * bx + vy * by + vz * bz

! choose the correct state depending on the sign of contact discontinuity
! advection speed
!
                if (sm <= 0.0d+00) then ! sm ≤ 0

! conservative variables for the inner left intermediate state
!
                  ui(idn) = dnr
                  ui(imx) = dnr * sm
                  ui(imy) = dnr * vy
                  ui(imz) = dnr * vz
                  ui(ibx) = bx
                  ui(iby) = by
                  ui(ibz) = bz
                  ui(ibp) = qr(ibp,i)
                  ui(ien) = (wcr(ien) + sm * pt - bx * vb) / car

! the inner right intermediate flux
!
                  f(:,i)  = smr * ui(:) - wcr(:)

                else ! sm > 0

! vector of the right-going Alfvén wave
!
                  wcr(:)  = (sm - smr) * ui(:) + wcr(:)

! calculate the average flux over the left inner intermediate state
!
                  f(:,i)  = (sm * wl(:) - sl * wcr(:)) / (sl - sm)

                end if ! sm > 0

              end if ! sr* > 0

            else ! sl* < sl & sr* > sr

! so far we revert to HLL flux in the case of degeneracies
!
              f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

            end if ! sl* < sl & sr* > sr

          end if ! one degeneracy

        else ! Bₓ² = 0

! the total pressure, constant across the contact discontinuity
!
          pt = (wl(idn) * wr(imx) - wr(idn) * wl(imx)) / dn

! separate intermediate states depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! the left speed difference
!
            slmm    =  sl - sm

! conservative variables for the left intermediate state
!
            ui(idn) =  wl(idn) / slmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * ql(ivy,i)
            ui(imz) =  ui(idn) * ql(ivz,i)
            ui(ibx) =  0.0d+00
            ui(iby) =  wl(iby) / slmm
            ui(ibz) =  wl(ibz) / slmm
            ui(ibp) =  ql(ibp,i)
            ui(ien) = (wl(ien) + sm * pt) / slmm

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! the right speed difference
!
            srmm    =  sr - sm

! conservative variables for the right intermediate state
!
            ui(idn) =  wr(idn) / srmm
            ui(imx) =  ui(idn) * sm
            ui(imy) =  ui(idn) * qr(ivy,i)
            ui(imz) =  ui(idn) * qr(ivz,i)
            ui(ibx) =  0.0d+00
            ui(iby) =  wr(iby) / srmm
            ui(ibz) =  wr(ibz) / srmm
            ui(ibp) =  qr(ibp,i)
            ui(ien) = (wr(ien) + sm * pt) / srmm

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! intermediate flux is constant across the contact discontinuity and all except
! the parallel momentum flux are zero
!
            f(idn,i) =   0.0d+00
            f(imx,i) = - 0.5d+00 * (wl(imx) + wr(imx))
            f(imy,i) =   0.0d+00
            f(imz,i) =   0.0d+00
            f(ibx,i) = fl(ibx,i)
            f(iby,i) =   0.0d+00
            f(ibz,i) =   0.0d+00
            f(ibp,i) = fl(ibp,i)
            f(ien,i) =   0.0d+00

          end if ! sm = 0

        end if ! Bₓ² = 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_adi_hlld

!===============================================================================
!
end module schemes
