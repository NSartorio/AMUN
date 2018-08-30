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
  integer            , save :: imi, imf, ims, imr
#endif /* PROFILE */

! 4-vector reconstruction flag
!
  logical            , save :: states_4vec = .false.

! interfaces for procedure pointers
!
  abstract interface
    subroutine update_flux_iface(dx, q, f)
      use coordinates    , only : im, jm, km
      use equations      , only : nv
      real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
      real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
      real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f
    end subroutine
    subroutine riemann_iface(n, ql, qr, f)
      use equations, only : nv
      integer                      , intent(in)    :: n
      real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
      real(kind=8), dimension(nv,n), intent(out)   :: f
    end subroutine
  end interface

! pointer to the flux update procedure
!
  procedure(update_flux_iface), pointer, save :: update_flux => null()

! pointer to the Riemann solver
!
  procedure(riemann_iface)    , pointer, save :: riemann     => null()

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_schemes, finalize_schemes
  public :: update_flux

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
    use equations      , only : eqsys, eos
    use parameters     , only : get_parameter_string

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
    character(len=64)      :: statev   = "primitive"
    character(len=255)     :: name_sol = ""
    character(len=255)     :: name_sts = "primitive"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('schemes:: initialization'  , imi)
    call set_timer('schemes:: flux update'     , imf)
    call set_timer('schemes:: Riemann states'  , ims)
    call set_timer('schemes:: Riemann solver'  , imr)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the Riemann solver
!
    call get_parameter_string("riemann_solver" , solver)
    call get_parameter_string("state_variables", statev)

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

        case("roe", "ROE")

! set the solver name
!
          name_sol =  "ROE"

! set pointers to subroutines
!
          riemann => riemann_hd_iso_roe

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

        case("roe", "ROE")

! set the solver name
!
          name_sol =  "ROE"

! set pointers to subroutines
!
          riemann => riemann_hd_adi_roe

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

        case("roe", "ROE")

! set the solver name
!
          name_sol =  "ROE"

! set pointers to subroutines
!
          riemann => riemann_mhd_iso_roe

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

        case("roe", "ROE")

! set the solver name
!
          name_sol =  "ROE"

! set pointers to subroutines
!
          riemann => riemann_mhd_adi_roe

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

!--- SPECIAL RELATIVITY HYDRODYNAMICS ---
!
    case("srhd", "SRHD")

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! set pointers to subroutines
!
        update_flux => update_flux_srhd_adi

! select the state reconstruction method
!
        select case(trim(statev))

        case("4vec", "4-vector", "4VEC", "4-VECTOR")

! set the state reconstruction name
!
          name_sts =  "4-vector"

! set 4-vector reconstruction flag
!
          states_4vec = .true.

! in the case of state variables, revert to primitive
!
        case default

! set the state reconstruction name
!
          name_sts =  "primitive"

        end select

! select the Riemann solver
!
        select case(trim(solver))

        case("hllc", "HLLC", "hllcm", "HLLCM", "hllc-m", "HLLC-M")

! set the solver name
!
          name_sol =  "HLLC (Mignone & Bodo 2005)"

! set pointers to subroutines
!
          riemann => riemann_srhd_adi_hllc

! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_srhd_adi_hll

        end select

      end select

!--- SPECIAL RELATIVITY MAGNETOHYDRODYNAMICS ---
!
    case("srmhd", "SRMHD")

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! set pointers to subroutines
!
        update_flux => update_flux_srmhd_adi

! select the state reconstruction method
!
        select case(trim(statev))

        case("4vec", "4-vector", "4VEC", "4-VECTOR")

! set the state reconstruction name
!
          name_sts =  "4-vector"

! set 4-vector reconstruction flag
!
          states_4vec = .true.

! in the case of state variables, revert to primitive
!
        case default

! set the state reconstruction name
!
          name_sts =  "primitive"

        end select

! select the Riemann solver
!
        select case(trim(solver))

        case("hllc", "HLLC", "hllcm", "HLLCM", "hllc-m", "HLLC-M")

! set the solver name
!
          name_sol =  "HLLC (Mignone & Bodo 2006)"

! set pointers to subroutines
!
          riemann => riemann_srmhd_adi_hllc

! in the case of unknown Riemann solver, revert to HLL
!
        case default

! set the solver name
!
          name_sol =  "HLL"

! set pointers to subroutines
!
          riemann => riemann_srmhd_adi_hll

        end select

      end select

    end select

! print information about the Riemann solver
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "Riemann solver         =", trim(name_sol)
      write (*,"(4x,a,1x,a)"    ) "state variables        =", trim(name_sts)

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
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
! subroutine RECONSTRUCT_INTERFACES:
! ---------------------------------
!
!   Subroutine reconstructs the Riemann states using 3D reconstruction methods.
!
!   Arguments:
!
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     qi   - the array of reconstructed states (2 in each direction);
!
!===============================================================================
!
  subroutine reconstruct_interfaces(dx, q, qi)

! include external procedures
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv, idn, ipr
    use interpolations , only : interfaces

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), dimension(NDIMS)              , intent(in)  :: dx
    real(kind=8), dimension(nv,im,jm,km)        , intent(in)  :: q
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS), intent(out) :: qi

! local variables
!
    logical                        :: positive
    integer                        :: p
!
!-------------------------------------------------------------------------------
!
! iterate over all variables
!
    do p = 1, nv

! determine if the variable is positive
!
      positive = (p == idn .or. p == ipr)

! interpolate interfaces
!
      call interfaces(positive, dx(1:NDIMS), q (p,1:im,1:jm,1:km),             &
                                             qi(p,1:im,1:jm,1:km,1:2,1:NDIMS))

    end do ! p = 1, nv

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_interfaces
!
!===============================================================================
!
!***** ISOTHERMAL HYDRODYNAMICS *****
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
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_hd_iso(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! reconstruct interfaces
!
    call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)                 &
                                     , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy states to directional lines with proper vector component ordering
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

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
! subroutine RIEMANN_HD_ISO_HLL:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_hd_iso_hll(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : ivx
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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
! subroutine RIEMANN_HD_ISO_ROE:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Roe's method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Roe, P. L.
!         "Approximate Riemann Solvers, Parameter Vectors, and Difference
!          Schemes",
!         Journal of Computational Physics, 1981, 43, pp. 357-372
!     [2] Toro, E. F.,
!         "Riemann Solvers and Numerical Methods for Fluid dynamics",
!         Springer-Verlag Berlin Heidelberg, 2009
!
!===============================================================================
!
  subroutine riemann_hd_iso_roe(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz
    use equations      , only : prim2cons, fluxspeed, eigensystem_roe

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                        :: p, i
    real(kind=8)                   :: sdl, sdr, sds
    real(kind=8)                   :: xx

! local arrays to store the states
!
    real(kind=8), dimension(nv,n)  :: ul, ur, fl, fr
    real(kind=8), dimension(nv)    :: qi, ci, al
    real(kind=8), dimension(nv,nv) :: li, ri
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:))

! iterate over all points
!
    do i = 1, n

! calculate variables for the Roe intermediate state
!
      sdl     = sqrt(ql(idn,i))
      sdr     = sqrt(qr(idn,i))
      sds     = sdl + sdr

! prepare the Roe intermediate state vector (eq. 11.60 in [2])
!
      qi(idn) =  sdl * sdr
      qi(ivx) = (sdl * ql(ivx,i) + sdr * qr(ivx,i)) / sds
      qi(ivy) = (sdl * ql(ivy,i) + sdr * qr(ivy,i)) / sds
      qi(ivz) = (sdl * ql(ivz,i) + sdr * qr(ivz,i)) / sds

! obtain eigenvalues and eigenvectors
!
      call eigensystem_roe(0.0d+00, 0.0d+00, qi(:), ci(:), ri(:,:), li(:,:))

! return upwind fluxes
!
      if (ci(1) >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (ci(nv) <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else

! calculate wave amplitudes α = L.ΔU
!
        al(1:nv) = 0.0d+00
        do p = 1, nv
          al(1:nv) = al(1:nv) + li(p,1:nv) * (ur(p,i) - ul(p,i))
        end do

! calculate the flux average
!
        f(1:nv,i) = 0.5d+00 * (fl(1:nv,i) + fr(1:nv,i))

! correct the flux by adding the characteristic wave contribution: ∑(α|λ|K)
!
        if (qi(ivx) >= 0.0d+00) then
          do p = 1, nv
            xx        = - 0.5d+00 * al(p) * abs(ci(p))
            f(1:nv,i) = f(1:nv,i) + xx * ri(p,1:nv)
          end do
        else
          do p = nv, 1, - 1
            xx        = - 0.5d+00 * al(p) * abs(ci(p))
            f(1:nv,i) = f(1:nv,i) + xx * ri(p,1:nv)
          end do
        end if

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_hd_iso_roe
!
!===============================================================================
!
!***** ADIABATIC HYDRODYNAMICS *****
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
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_hd_adi(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! reconstruct interfaces
!
    call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)                 &
                                     , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy states to directional lines with proper vector component ordering
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)
        qx(ipr,1:im,1:2) = qs(ipr,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
        f(1,ien,1:im,j,k) = fx(ien,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)
        qy(ipr,1:jm,1:2) = qs(ipr,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
        f(2,ien,i,1:jm,k) = fy(ien,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)
        qz(ipr,1:km,1:2) = qs(ipr,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
        f(3,ien,i,j,1:km) = fz(ien,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

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
! subroutine RIEMANN_HD_ADI_HLL:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_hd_adi_hll(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : ivx
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all position
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate the inverse of speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(1:nv)  = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv)  = sr * ur(1:nv,i) - fr(1:nv,i)

! calculate fluxes for the intermediate state
!
        f(1:nv,i) = (sl * wr(1:nv) - sr * wl(1:nv)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for the Riemann solver
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
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Toro, E. F., Spruce, M., & Speares, W.
!         "Restoration of the contact surface in the HLL-Riemann solver",
!         Shock Waves, 1994, Volume 4, Issue 1, pp. 25-34
!
!===============================================================================
!
  subroutine riemann_hd_adi_hllc(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, imx, imy, imz, ien
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, sm
    real(kind=8)                  :: slmm, srmm
    real(kind=8)                  :: dn, pr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, ui
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate vectors of the left and right-going waves
!
        wl(1:nv) = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv) = sr * ur(1:nv,i) - fr(1:nv,i)

! the speed of contact discontinuity
!
        dn =  wr(idn) - wl(idn)
        sm = (wr(imx) - wl(imx)) / dn

! calculate the pressure of the intermediate state
!
        pr = (wl(idn) * wr(imx) - wr(idn) * wl(imx)) / dn

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
          f(1:nv,i)  = sl * ui(1:nv) - wl(1:nv)

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
          f(1:nv,i)  = sr * ui(1:nv) - wr(1:nv)

        else ! sm = 0

! intermediate flux is constant across the contact discontinuity and all except
! the parallel momentum flux are zero
!
          f(idn,i) =   0.0d+00
          f(imx,i) = - wl(imx)
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
! subroutine RIEMANN_HD_ADI_ROE:
! -----------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Roe's method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Roe, P. L.
!         "Approximate Riemann Solvers, Parameter Vectors, and Difference
!          Schemes",
!         Journal of Computational Physics, 1981, 43, pp. 357-372
!     [2] Toro, E. F.,
!         "Riemann Solvers and Numerical Methods for Fluid dynamics",
!         Springer-Verlag Berlin Heidelberg, 2009
!
!===============================================================================
!
  subroutine riemann_hd_adi_roe(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ien
    use equations      , only : prim2cons, fluxspeed, eigensystem_roe

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                        :: p, i
    real(kind=8)                   :: sdl, sdr, sds
    real(kind=8)                   :: xx

! local arrays to store the states
!
    real(kind=8), dimension(nv,n)  :: ul, ur, fl, fr
    real(kind=8), dimension(nv)    :: qi, ci, al
    real(kind=8), dimension(nv,nv) :: li, ri
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:))

! iterate over all points
!
    do i = 1, n

! calculate variables for the Roe intermediate state
!
      sdl     = sqrt(ql(idn,i))
      sdr     = sqrt(qr(idn,i))
      sds     = sdl + sdr

! prepare the Roe intermediate state vector (eq. 11.60 in [2])
!
      qi(idn) =  sdl * sdr
      qi(ivx) = (sdl * ql(ivx,i) + sdr * qr(ivx,i)) / sds
      qi(ivy) = (sdl * ql(ivy,i) + sdr * qr(ivy,i)) / sds
      qi(ivz) = (sdl * ql(ivz,i) + sdr * qr(ivz,i)) / sds
      qi(ipr) = ((ul(ien,i) + ql(ipr,i)) / sdl                                 &
              +  (ur(ien,i) + qr(ipr,i)) / sdr) / sds

! obtain eigenvalues and eigenvectors
!
      call eigensystem_roe(0.0d+00, 0.0d+00, qi(:), ci(:), ri(:,:), li(:,:))

! return upwind fluxes
!
      if (ci(1) >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (ci(nv) <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else

! calculate wave amplitudes α = L.ΔU
!
        al(1:nv) = 0.0d+00
        do p = 1, nv
          al(1:nv) = al(1:nv) + li(p,1:nv) * (ur(p,i) - ul(p,i))
        end do

! calculate the flux average
!
        f(1:nv,i) = 0.5d+00 * (fl(1:nv,i) + fr(1:nv,i))

! correct the flux by adding the characteristic wave contribution: ∑(α|λ|K)
!
        if (qi(ivx) >= 0.0d+00) then
          do p = 1, nv
            xx        = - 0.5d+00 * al(p) * abs(ci(p))
            f(1:nv,i) = f(1:nv,i) + xx * ri(p,1:nv)
          end do
        else
          do p = nv, 1, - 1
            xx        = - 0.5d+00 * al(p) * abs(ci(p))
            f(1:nv,i) = f(1:nv,i) + xx * ri(p,1:nv)
          end do
        end if

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_hd_adi_roe
!
!===============================================================================
!
!***** ISOTHERMAL MAGNETOHYDRODYNAMICS *****
!
!===============================================================================
!
! subroutine UPDATE_FLUX_MHD_ISO:
! ------------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_mhd_iso(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz
    use equations      , only : ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! reconstruct interfaces
!
    call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)                 &
                                     , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)
        qx(ibx,1:im,1:2) = qs(ibx,1:im,j,k,1:2,1)
        qx(iby,1:im,1:2) = qs(iby,1:im,j,k,1:2,1)
        qx(ibz,1:im,1:2) = qs(ibz,1:im,j,k,1:2,1)
        qx(ibp,1:im,1:2) = qs(ibp,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
        f(1,ibx,1:im,j,k) = fx(ibx,1:im)
        f(1,iby,1:im,j,k) = fx(iby,1:im)
        f(1,ibz,1:im,j,k) = fx(ibz,1:im)
        f(1,ibp,1:im,j,k) = fx(ibp,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)
        qy(ibx,1:jm,1:2) = qs(iby,i,1:jm,k,1:2,2)
        qy(iby,1:jm,1:2) = qs(ibz,i,1:jm,k,1:2,2)
        qy(ibz,1:jm,1:2) = qs(ibx,i,1:jm,k,1:2,2)
        qy(ibp,1:jm,1:2) = qs(ibp,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
        f(2,ibx,i,1:jm,k) = fy(ibz,1:jm)
        f(2,iby,i,1:jm,k) = fy(ibx,1:jm)
        f(2,ibz,i,1:jm,k) = fy(iby,1:jm)
        f(2,ibp,i,1:jm,k) = fy(ibp,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)
        qz(ibx,1:km,1:2) = qs(ibz,i,j,1:km,1:2,3)
        qz(iby,1:km,1:2) = qs(ibx,i,j,1:km,1:2,3)
        qz(ibz,1:km,1:2) = qs(iby,i,j,1:km,1:2,3)
        qz(ibp,1:km,1:2) = qs(ibp,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
        f(3,ibx,i,j,1:km) = fz(iby,1:km)
        f(3,iby,i,j,1:km) = fz(ibz,1:km)
        f(3,ibz,i,j,1:km) = fz(ibx,1:km)
        f(3,ibp,i,j,1:km) = fz(ibp,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

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
! subroutine RIEMANN_MHD_ISO_HLL:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_mhd_iso_hll(n, ql, qr, f)

! include external variables and procedures
!
    use equations      , only : nv
    use equations      , only : ivx, ibx, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml
    real(kind=8)                  :: bx, bp

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_mhd_iso_hlld(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, imx, imy, imz, ibx, iby, ibz, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, sm, sml, smr, srml, slmm, srmm
    real(kind=8)                  :: bx, bp, b2, dn, dnl, dnr, dvl, dvr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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

! apply full HLLD solver only if the Alfvén wave is strong enough
!
        if (b2 > 1.0d-08 * max(dnl * sl**2, dnr * sr**2)) then ! Bₓ² > ε

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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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

        else if (b2 > 0.0d+00) then ! Bₓ² > 0

! the Alfvén wave is too weak, so ignore it in order to not introduce any
! numerical instabilities; since Bₓ is not zero, the perpendicular components
! of velocity and magnetic field are continuous; we have only one intermediate
! state, therefore simply apply the HLL solver
!
          f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

        else ! Bₓ² = 0

! take the right flux depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! conservative variables for the left intermediate state
!
            ui(idn) = dnl
            ui(imx) = dnl * sm
            ui(imy) = wl(imy) / slmm
            ui(imz) = wl(imz) / slmm
            ui(ibx) = 0.0d+00
            ui(iby) = wl(iby) / slmm
            ui(ibz) = wl(ibz) / slmm
            ui(ibp) = ul(ibp,i)

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the right intermediate state
!
            ui(idn) = dnr
            ui(imx) = dnr * sm
            ui(imy) = wr(imy) / srmm
            ui(imz) = wr(imz) / srmm
            ui(ibx) = 0.0d+00
            ui(iby) = wr(iby) / srmm
            ui(ibz) = wr(ibz) / srmm
            ui(ibp) = ur(ibp,i)

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! the intermediate flux; since the advection speed is zero, perpendicular
! components do not change, so revert it to HLL
!
            f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

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
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Mignone, A.,
!         "A simple and accurate Riemann solver for isothermal MHD",
!         Journal of Computational Physics, 2007, 225, pp. 1427-1441
!
!===============================================================================
!
  subroutine riemann_mhd_iso_hlldm(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, imx, imy, imz, ibx, iby, ibz, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, sm, sml, smr, srml, slmm, srmm
    real(kind=8)                  :: bx, bp, b2, dn, dnl, dnr, dvl, dvr, ca, ca2

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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

! get the square of the Alfvén speed
!
        ca2 = b2 / dn

! apply full HLLD solver only if the Alfvén wave is strong enough
!
        if (ca2 > 1.0d-08 * max(sl**2, sr**2)) then ! Bₓ² > ε

! left and right Alfvén speeds
!
          ca   = sqrt(ca2)
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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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
              ui(ibp) = ul(ibp,i)

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
              ui(ibp) = ur(ibp,i)

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

        else if (b2 > 0.0d+00) then ! Bₓ² > 0

! the Alfvén wave is very weak, so ignore it in order to not introduce any
! numerical instabilities; since Bₓ is not zero, the perpendicular components
! of velocity and magnetic field are continuous; we have only one intermediate
! state, therefore simply apply the HLL solver
!
          f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

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
            ui(ibp) = ul(ibp,i)

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
            ui(ibp) = ur(ibp,i)

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! both states are equal so revert to the HLL flux
!
            f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

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
! subroutine RIEMANN_MHD_ISO_ROE:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Roe's method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!     [2] Toro, E. F.,
!         "Riemann Solvers and Numerical Methods for Fluid dynamics",
!         Springer-Verlag Berlin Heidelberg, 2009
!
!===============================================================================
!
  subroutine riemann_mhd_iso_roe(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp
    use equations      , only : imx, imy, imz
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed, eigensystem_roe

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                        :: p, i
    real(kind=8)                   :: sdl, sdr, sds
    real(kind=8)                   :: bx, bp
    real(kind=8)                   :: xx, yy

! local arrays to store the states
!
    real(kind=8), dimension(nv,n)  :: ul, ur, fl, fr
    real(kind=8), dimension(nv)    :: qi, ci, al
    real(kind=8), dimension(nv,nv) :: li, ri
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:))

! iterate over all points
!
    do i = 1, n

! calculate variables for the Roe intermediate state
!
      sdl     = sqrt(ql(idn,i))
      sdr     = sqrt(qr(idn,i))
      sds     = sdl + sdr

! prepare the Roe intermediate state vector (eq. 11.60 in [2])
!
      qi(idn) =  sdl * sdr
      qi(ivx) = (sdl * ql(ivx,i) + sdr * qr(ivx,i)) / sds
      qi(ivy) = (sdl * ql(ivy,i) + sdr * qr(ivy,i)) / sds
      qi(ivz) = (sdl * ql(ivz,i) + sdr * qr(ivz,i)) / sds
      qi(ibx) = ql(ibx,i)
      qi(iby) = (sdr * ql(iby,i) + sdl * qr(iby,i)) / sds
      qi(ibz) = (sdr * ql(ibz,i) + sdl * qr(ibz,i)) / sds
      qi(ibp) = ql(ibp,i)

! prepare coefficients
!
      xx = 0.5d+00 * ((ql(iby,i) - qr(iby,i))**2                               &
                                        + (ql(ibz,i) - qr(ibz,i))**2) / sds**2
      yy = 0.5d+00 * (ql(idn,i) + qr(idn,i)) / qi(idn)

! obtain eigenvalues and eigenvectors
!
      call eigensystem_roe(xx, yy, qi(:), ci(:), ri(:,:), li(:,:))

! return upwind fluxes
!
      if (ci(1) >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (ci(nv) <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else

! prepare fluxes which do not change across the states
!
        f(ibx,i) = fl(ibx,i)
        f(ibp,i) = fl(ibp,i)

! calculate wave amplitudes α = L.ΔU
!
        al(1:nv) = 0.0d+00
        do p = 1, nv
          al(1:nv) = al(1:nv) + li(p,1:nv) * (ur(p,i) - ul(p,i))
        end do

! calculate the flux average
!
        f(idn,i) = 0.5d+00 * (fl(idn,i) + fr(idn,i))
        f(imx,i) = 0.5d+00 * (fl(imx,i) + fr(imx,i))
        f(imy,i) = 0.5d+00 * (fl(imy,i) + fr(imy,i))
        f(imz,i) = 0.5d+00 * (fl(imz,i) + fr(imz,i))
        f(iby,i) = 0.5d+00 * (fl(iby,i) + fr(iby,i))
        f(ibz,i) = 0.5d+00 * (fl(ibz,i) + fr(ibz,i))

! correct the flux by adding the characteristic wave contribution: ∑(α|λ|K)
!
        if (qi(ivx) >= 0.0d+00) then
          do p = 1, nv
            xx       = - 0.5d+00 * al(p) * abs(ci(p))

            f(idn,i) = f(idn,i) + xx * ri(p,idn)
            f(imx,i) = f(imx,i) + xx * ri(p,imx)
            f(imy,i) = f(imy,i) + xx * ri(p,imy)
            f(imz,i) = f(imz,i) + xx * ri(p,imz)
            f(iby,i) = f(iby,i) + xx * ri(p,iby)
            f(ibz,i) = f(ibz,i) + xx * ri(p,ibz)
          end do
        else
          do p = nv, 1, -1
            xx       = - 0.5d+00 * al(p) * abs(ci(p))

            f(idn,i) = f(idn,i) + xx * ri(p,idn)
            f(imx,i) = f(imx,i) + xx * ri(p,imx)
            f(imy,i) = f(imy,i) + xx * ri(p,imy)
            f(imz,i) = f(imz,i) + xx * ri(p,imz)
            f(iby,i) = f(iby,i) + xx * ri(p,iby)
            f(ibz,i) = f(ibz,i) + xx * ri(p,ibz)
          end do
        end if

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_iso_roe
!
!===============================================================================
!
!***** ADIABATIC MAGNETOHYDRODYNAMICS *****
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
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_mhd_adi(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien
    use equations      , only : ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! reconstruct interfaces
!
    call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)                 &
                                     , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)
        qx(ibx,1:im,1:2) = qs(ibx,1:im,j,k,1:2,1)
        qx(iby,1:im,1:2) = qs(iby,1:im,j,k,1:2,1)
        qx(ibz,1:im,1:2) = qs(ibz,1:im,j,k,1:2,1)
        qx(ibp,1:im,1:2) = qs(ibp,1:im,j,k,1:2,1)
        qx(ipr,1:im,1:2) = qs(ipr,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
        f(1,ibx,1:im,j,k) = fx(ibx,1:im)
        f(1,iby,1:im,j,k) = fx(iby,1:im)
        f(1,ibz,1:im,j,k) = fx(ibz,1:im)
        f(1,ibp,1:im,j,k) = fx(ibp,1:im)
        f(1,ien,1:im,j,k) = fx(ien,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)
        qy(ibx,1:jm,1:2) = qs(iby,i,1:jm,k,1:2,2)
        qy(iby,1:jm,1:2) = qs(ibz,i,1:jm,k,1:2,2)
        qy(ibz,1:jm,1:2) = qs(ibx,i,1:jm,k,1:2,2)
        qy(ibp,1:jm,1:2) = qs(ibp,i,1:jm,k,1:2,2)
        qy(ipr,1:jm,1:2) = qs(ipr,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
        f(2,ibx,i,1:jm,k) = fy(ibz,1:jm)
        f(2,iby,i,1:jm,k) = fy(ibx,1:jm)
        f(2,ibz,i,1:jm,k) = fy(iby,1:jm)
        f(2,ibp,i,1:jm,k) = fy(ibp,1:jm)
        f(2,ien,i,1:jm,k) = fy(ien,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)
        qz(ibx,1:km,1:2) = qs(ibz,i,j,1:km,1:2,3)
        qz(iby,1:km,1:2) = qs(ibx,i,j,1:km,1:2,3)
        qz(ibz,1:km,1:2) = qs(iby,i,j,1:km,1:2,3)
        qz(ibp,1:km,1:2) = qs(ibp,i,j,1:km,1:2,3)
        qz(ipr,1:km,1:2) = qs(ipr,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
        f(3,ibx,i,j,1:km) = fz(iby,1:km)
        f(3,iby,i,j,1:km) = fz(ibz,1:km)
        f(3,ibz,i,j,1:km) = fz(ibx,1:km)
        f(3,ibp,i,j,1:km) = fz(ibp,1:km)
        f(3,ien,i,j,1:km) = fz(ien,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

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
! subroutine RIEMANN_MHD_ADI_HLL:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_mhd_adi_hll(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : ivx, ibx, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml
    real(kind=8)                  :: bx, bp

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_mhd_adi_hllc(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp, ipr
    use equations      , only : imx, imy, imz, ien
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, sm, srml, slmm, srmm
    real(kind=8)                  :: dn, bx, bp, b2, pt, vy, vz, by, bz, vb

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, ui
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

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
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_mhd_adi_hlld(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp, ipr
    use equations      , only : imx, imy, imz, ien
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, sm, srml, slmm, srmm
    real(kind=8)                  :: dn, bx, bp, b2, pt, vy, vz, by, bz, vb
    real(kind=8)                  :: dnl, dnr, cal, car, sml, smr
    real(kind=8)                  :: dv, dvl, dvr

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr, wcl, wcr, ui
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all points
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLLD flux
!
      if (sl >= 0.0d+00) then ! sl ≥ 0

        f(:,i) = fl(:,i)

      else if (sr <= 0.0d+00) then ! sr ≤ 0

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

! speed differences
!
        slmm = sl - sm
        srmm = sr - sm

! left and right state densities
!
        dnl = wl(idn) / slmm
        dnr = wr(idn) / srmm

! square of Bₓ, i.e. Bₓ²
!
        bx = ql(ibx,i)
        b2 = ql(ibx,i) * qr(ibx,i)

! the total pressure, constant across the contact discontinuity and Alfvén waves
!
        pt = (wl(idn) * wr(imx) - wr(idn) * wl(imx)) / dn + b2

! if Alfvén wave is strong enough, apply the full HLLD solver, otherwise
! revert to the HLLC one
!
        if (b2 > 1.0d-08 * max(dnl * sl**2, dnr * sr**2)) then ! Bₓ² > ε

! left and right Alfvén speeds
!
          cal = sqrt(b2 / dnl)
          car = sqrt(b2 / dnr)
          sml = sm - cal
          smr = sm + car

! calculate division factors (also used to determine degeneracies)
!
          dvl = slmm * wl(idn) - b2
          dvr = srmm * wr(idn) - b2

! check degeneracy Sl* -> Sl or Sr* -> Sr
!
          if (dvl > 0.0d+00 .and. dvr > 0.0d+00) then ! no degeneracy

! choose the correct state depending on the speed signs
!
            if (sml >= 0.0d+00) then ! sl* ≥ 0

! primitive variables for the outer left intermediate state
!
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dvl
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dvl
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
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
              ui(ibp) = ul(ibp,i)
              ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! the outer left intermediate flux
!
              f(:,i)  = sl * ui(:) - wl(:)

            else if (smr <= 0.0d+00) then ! sr* ≤ 0

! primitive variables for the outer right intermediate state
!
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dvr
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dvr
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
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
              ui(ibp) = ur(ibp,i)
              ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! the outer right intermediate flux
!
              f(:,i)  = sr * ui(:) - wr(:)

            else ! sl* < 0 < sr*

! primitive variables for the outer left intermediate state
!
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dvl
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dvl
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
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
              ui(ibp) = ul(ibp,i)
              ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! vector of the left-going Alfvén wave
!
              wcl(:)  = (sml - sl) * ui(:) + wl(:)

! primitive variables for the outer right intermediate state
!
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dvr
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dvr
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
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
              ui(ibp) = ur(ibp,i)
              ui(ien) = (wr(ien) + sm * pt - bx * vb) / srmm

! vector of the right-going Alfvén wave
!
              wcr(:)  = (smr - sr) * ui(:) + wr(:)

! prepare constant primitive variables of the intermediate states
!
              dv      = car * dnr + cal * dnl
              vy      = (wcr(imy) - wcl(imy)) / dv
              vz      = (wcr(imz) - wcl(imz)) / dv
              dv      = car       + cal
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
                ui(ibp) = ul(ibp,i)
                ui(ien) = - (wcl(ien) + sm * pt - bx * vb) / cal

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
                ui(ibp) = ur(ibp,i)
                ui(ien) =   (wcr(ien) + sm * pt - bx * vb) / car

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
            if (dvl > 0.0d+00) then ! sr* > sr

! primitive variables for the outer left intermediate state
!
              vy      = (   slmm * wl(imy) - bx * wl(iby)) / dvl
              vz      = (   slmm * wl(imz) - bx * wl(ibz)) / dvl
              by      = (wl(idn) * wl(iby) - bx * wl(imy)) / dvl
              bz      = (wl(idn) * wl(ibz) - bx * wl(imz)) / dvl
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
              ui(ibp) = ul(ibp,i)
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
                dv      = srmm * dnr + cal * dnl
                vy      = (wr(imy) - wcl(imy)) / dv
                vz      = (wr(imz) - wcl(imz)) / dv
                dv      = sr - sml
                by      = (wr(iby) - wcl(iby)) / dv
                bz      = (wr(ibz) - wcl(ibz)) / dv
                vb      = sm * bx + vy * by + vz * bz

! conservative variables for the inner left intermediate state
!
                ui(idn) = dnl
                ui(imx) = dnl * sm
                ui(imy) = dnl * vy
                ui(imz) = dnl * vz
                ui(ibx) = bx
                ui(iby) = by
                ui(ibz) = bz
                ui(ibp) = ul(ibp,i)
                ui(ien) = - (wcl(ien) + sm * pt - bx * vb) / cal

! choose the correct state depending on the sign of contact discontinuity
! advection speed
!
                if (sm >= 0.0d+00) then ! sm ≥ 0

! the inner left intermediate flux
!
                  f(:,i)  = sml * ui(:) - wcl(:)

                else ! sm < 0

! vector of the left-going Alfvén wave
!
                  wcr(:)  = (sm - sml) * ui(:) + wcl(:)

! calculate the average flux over the right inner intermediate state
!
                  f(:,i)  = (sm * wr(:) - sr * wcr(:)) / srmm

                end if ! sm < 0

              end if ! sl* < 0

            else if (dvr > 0.0d+00) then ! sl* < sl

! primitive variables for the outer right intermediate state
!
              vy      = (   srmm * wr(imy) - bx * wr(iby)) / dvr
              vz      = (   srmm * wr(imz) - bx * wr(ibz)) / dvr
              by      = (wr(idn) * wr(iby) - bx * wr(imy)) / dvr
              bz      = (wr(idn) * wr(ibz) - bx * wr(imz)) / dvr
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
              ui(ibp) = ur(ibp,i)
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

! conservative variables for the inner left intermediate state
!
                ui(idn) = dnr
                ui(imx) = dnr * sm
                ui(imy) = dnr * vy
                ui(imz) = dnr * vz
                ui(ibx) = bx
                ui(iby) = by
                ui(ibz) = bz
                ui(ibp) = ur(ibp,i)
                ui(ien) =   (wcr(ien) + sm * pt - bx * vb) / car

! choose the correct state depending on the sign of contact discontinuity
! advection speed
!
                if (sm <= 0.0d+00) then ! sm ≤ 0

! the inner right intermediate flux
!
                  f(:,i)  = smr * ui(:) - wcr(:)

                else ! sm > 0

! vector of the right-going Alfvén wave
!
                  wcl(:)  = (sm - smr) * ui(:) + wcr(:)

! calculate the average flux over the left inner intermediate state
!
                  f(:,i)  = (sm * wl(:) - sl * wcl(:)) / slmm

                end if ! sm > 0

              end if ! sr* > 0

            else ! sl* < sl & sr* > sr

! so far we revert to HLL flux in the case of degeneracies
!
              f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

            end if ! sl* < sl & sr* > sr

          end if ! one degeneracy

        else if (b2 > 0.0d+00) then ! Bₓ² > 0

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

! conservative variables for the left intermediate state
!
            ui(idn) =  dnl
            ui(imx) =  dnl * sm
            ui(imy) =  dnl * vy
            ui(imz) =  dnl * vz
            ui(ibx) =  bx
            ui(iby) =  by
            ui(ibz) =  bz
            ui(ibp) =  ul(ibp,i)
            ui(ien) = (wl(ien) + sm * pt - bx * vb) / slmm

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the right intermediate state
!
            ui(idn) =  dnr
            ui(imx) =  dnr * sm
            ui(imy) =  dnr * vy
            ui(imz) =  dnr * vz
            ui(ibx) =  bx
            ui(iby) =  by
            ui(ibz) =  bz
            ui(ibp) =  ur(ibp,i)
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

        else ! Bₓ² = 0

! separate intermediate states depending on the sign of the advection speed
!
          if (sm > 0.0d+00) then ! sm > 0

! conservative variables for the left intermediate state
!
            ui(idn) =  dnl
            ui(imx) =  dnl * sm
            ui(imy) =  wl(imy) / slmm
            ui(imz) =  wl(imz) / slmm
            ui(ibx) =  0.0d+00
            ui(iby) =  wl(iby) / slmm
            ui(ibz) =  wl(ibz) / slmm
            ui(ibp) =  ul(ibp,i)
            ui(ien) = (wl(ien) + sm * pt) / slmm

! the left intermediate flux
!
            f(:,i)  = sl * ui(:) - wl(:)

          else if (sm < 0.0d+00) then ! sm < 0

! conservative variables for the right intermediate state
!
            ui(idn) =  dnr
            ui(imx) =  dnr * sm
            ui(imy) =  wr(imy) / srmm
            ui(imz) =  wr(imz) / srmm
            ui(ibx) =  0.0d+00
            ui(iby) =  wr(iby) / srmm
            ui(ibz) =  wr(ibz) / srmm
            ui(ibp) =  ur(ibp,i)
            ui(ien) = (wr(ien) + sm * pt) / srmm

! the right intermediate flux
!
            f(:,i)  = sr * ui(:) - wr(:)

          else ! sm = 0

! when Sₘ = 0 all variables are continuous, therefore the flux reduces
! to the HLL one
!
            f(:,i) = (sl * wr(:) - sr * wl(:)) / srml

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
!
!===============================================================================
!
! subroutine RIEMANN_MHD_ADI_ROE:
! ------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Roe's method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!     [2] Toro, E. F.,
!         "Riemann Solvers and Numerical Methods for Fluid dynamics",
!         Springer-Verlag Berlin Heidelberg, 2009
!
!===============================================================================
!
  subroutine riemann_mhd_adi_roe(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use equations      , only : imx, imy, imz, ien
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed, eigensystem_roe

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                        :: p, i
    real(kind=8)                   :: sdl, sdr, sds
    real(kind=8)                   :: bx, bp
    real(kind=8)                   :: pml, pmr
    real(kind=8)                   :: xx, yy

! local arrays to store the states
!
    real(kind=8), dimension(nv,n)  :: ul, ur, fl, fr
    real(kind=8), dimension(nv)    :: qi, ci, al
    real(kind=8), dimension(nv,nv) :: li, ri
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate corresponding conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at the states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:))

! iterate over all points
!
    do i = 1, n

! calculate variables for the Roe intermediate state
!
      sdl     = sqrt(ql(idn,i))
      sdr     = sqrt(qr(idn,i))
      sds     = sdl + sdr

! prepare magnetic pressures
!
      pml     = 0.5d+00 * sum(ql(ibx:ibz,i)**2)
      pmr     = 0.5d+00 * sum(qr(ibx:ibz,i)**2)

! prepare the Roe intermediate state vector (eq. 11.60 in [2])
!
      qi(idn) =  sdl * sdr
      qi(ivx) = (sdl * ql(ivx,i) + sdr * qr(ivx,i)) / sds
      qi(ivy) = (sdl * ql(ivy,i) + sdr * qr(ivy,i)) / sds
      qi(ivz) = (sdl * ql(ivz,i) + sdr * qr(ivz,i)) / sds
      qi(ipr) = ((ul(ien,i) + ql(ipr,i) + pml) / sdl                           &
              +  (ur(ien,i) + qr(ipr,i) + pmr) / sdr) / sds
      qi(ibx) = ql(ibx,i)
      qi(iby) = (sdr * ql(iby,i) + sdl * qr(iby,i)) / sds
      qi(ibz) = (sdr * ql(ibz,i) + sdl * qr(ibz,i)) / sds
      qi(ibp) = ql(ibp,i)

! prepare coefficients
!
      xx = 0.5d+00 * ((ql(iby,i) - qr(iby,i))**2                               &
                                        + (ql(ibz,i) - qr(ibz,i))**2) / sds**2
      yy = 0.5d+00 * (ql(idn,i) + qr(idn,i)) / qi(idn)

! obtain eigenvalues and eigenvectors
!
      call eigensystem_roe(xx, yy, qi(:), ci(:), ri(:,:), li(:,:))

! return upwind fluxes
!
      if (ci(1) >= 0.0d+00) then

        f(:,i) = fl(:,i)

      else if (ci(nv) <= 0.0d+00) then

        f(:,i) = fr(:,i)

      else

! prepare fluxes which do not change across the states
!
        f(ibx,i) = fl(ibx,i)
        f(ibp,i) = fl(ibp,i)

! calculate wave amplitudes α = L.ΔU
!
        al(1:nv) = 0.0d+00
        do p = 1, nv
          al(1:nv) = al(1:nv) + li(p,1:nv) * (ur(p,i) - ul(p,i))
        end do

! calculate the flux average
!
        f(idn,i) = 0.5d+00 * (fl(idn,i) + fr(idn,i))
        f(imx,i) = 0.5d+00 * (fl(imx,i) + fr(imx,i))
        f(imy,i) = 0.5d+00 * (fl(imy,i) + fr(imy,i))
        f(imz,i) = 0.5d+00 * (fl(imz,i) + fr(imz,i))
        f(ien,i) = 0.5d+00 * (fl(ien,i) + fr(ien,i))
        f(iby,i) = 0.5d+00 * (fl(iby,i) + fr(iby,i))
        f(ibz,i) = 0.5d+00 * (fl(ibz,i) + fr(ibz,i))

! correct the flux by adding the characteristic wave contribution: ∑(α|λ|K)
!
        if (qi(ivx) >= 0.0d+00) then
          do p = 1, nv
            xx       = - 0.5d+00 * al(p) * abs(ci(p))

            f(idn,i) = f(idn,i) + xx * ri(p,idn)
            f(imx,i) = f(imx,i) + xx * ri(p,imx)
            f(imy,i) = f(imy,i) + xx * ri(p,imy)
            f(imz,i) = f(imz,i) + xx * ri(p,imz)
            f(ien,i) = f(ien,i) + xx * ri(p,ien)
            f(iby,i) = f(iby,i) + xx * ri(p,iby)
            f(ibz,i) = f(ibz,i) + xx * ri(p,ibz)
          end do
        else
          do p = nv, 1, -1
            xx       = - 0.5d+00 * al(p) * abs(ci(p))

            f(idn,i) = f(idn,i) + xx * ri(p,idn)
            f(imx,i) = f(imx,i) + xx * ri(p,imx)
            f(imy,i) = f(imy,i) + xx * ri(p,imy)
            f(imz,i) = f(imz,i) + xx * ri(p,imz)
            f(ien,i) = f(ien,i) + xx * ri(p,ien)
            f(iby,i) = f(iby,i) + xx * ri(p,iby)
            f(ibz,i) = f(ibz,i) + xx * ri(p,ibz)
          end do
        end if

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_mhd_adi_roe
!
!===============================================================================
!
!***** ADIABATIC SPECIAL RELATIVITY HYDRODYNAMICS *****
!
!===============================================================================
!
! subroutine UPDATE_FLUX_SRHD_ADI:
! -------------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_srhd_adi(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k, l, p
    real(kind=8)                   :: vm

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km)         :: qq
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! apply reconstruction on variables using 4-vector or 3-vector
!
    if (states_4vec) then

! convert velocities to four-velocities for physical reconstruction
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im

            vm = sqrt(1.0d+00 - sum(q(ivx:ivz,i,j,k)**2))

            qq(idn,i,j,k) = q(idn,i,j,k) / vm
            qq(ivx,i,j,k) = q(ivx,i,j,k) / vm
            qq(ivy,i,j,k) = q(ivy,i,j,k) / vm
            qq(ivz,i,j,k) = q(ivz,i,j,k) / vm
            qq(ipr,i,j,k) = q(ipr,i,j,k)

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

! reconstruct interfaces
!
      call reconstruct_interfaces(dx(:), qq(1:nv,1:im,1:jm,1:km)               &
                                       , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

! convert state four-velocities back to velocities
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im

            do l = 1, 2
              do p = 1, NDIMS

                vm = sqrt(1.0d+00 + sum(qs(ivx:ivz,i,j,k,l,p)**2))

                qs(idn,i,j,k,l,p) = qs(idn,i,j,k,l,p) / vm
                qs(ivx,i,j,k,l,p) = qs(ivx,i,j,k,l,p) / vm
                qs(ivy,i,j,k,l,p) = qs(ivy,i,j,k,l,p) / vm
                qs(ivz,i,j,k,l,p) = qs(ivz,i,j,k,l,p) / vm

              end do ! p = 1, ndims
            end do ! l = 1, 2
          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

    else

! reconstruct interfaces
!
      call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)               &
                                       , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

    end if

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)
        qx(ipr,1:im,1:2) = qs(ipr,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
        f(1,ien,1:im,j,k) = fx(ien,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)
        qy(ipr,1:jm,1:2) = qs(ipr,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
        f(2,ien,i,1:jm,k) = fy(ien,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)
        qz(ipr,1:km,1:2) = qs(ipr,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
        f(3,ien,i,j,1:km) = fz(ien,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_srhd_adi
!
!===============================================================================
!
! subroutine RIEMANN_SRHD_ADI_HLL:
! -------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_srhd_adi_hll(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : ivx
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all position
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate the inverse of speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(1:nv)  = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv)  = sr * ur(1:nv,i) - fr(1:nv,i)

! calculate fluxes for the intermediate state
!
        f(1:nv,i) = (sl * wr(1:nv) - sr * wl(1:nv)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for the Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_srhd_adi_hll
!
!===============================================================================
!
! subroutine RIEMANN_SRHD_ADI_HLLC:
! --------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer method with contact discontinuity resolution (HLLC)
!   by Mignone & Bodo.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Mignone, A. & Bodo, G.
!         "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics",
!         Monthly Notices of the Royal Astronomical Society,
!         2005, Volume 364, Pages 126-136
!
!===============================================================================
!
  subroutine riemann_srhd_adi_hllc(n, ql, qr, f)

! include external procedures
!
    use algebra        , only : quadratic
    use equations      , only : nv
    use equations      , only : ivx, idn, imx, imy, imz, ien
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i, nr
    real(kind=8)                  :: sl, sr, srml, sm
    real(kind=8)                  :: pr, dv, fc

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: uh, us, fh, wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
    real(kind=8), dimension(3)    :: a
    real(kind=8), dimension(2)    :: x
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all position
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate the inverse of speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(1:nv)  = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv)  = sr * ur(1:nv,i) - fr(1:nv,i)

! calculate fluxes for the intermediate state
!
        uh(1:nv)  = (     wr(1:nv) -      wl(1:nv)) / srml
        fh(1:nv)  = (sl * wr(1:nv) - sr * wl(1:nv)) / srml

! correct the energy waves
!
        wl(ien)   = wl(ien) + wl(idn)
        wr(ien)   = wr(ien) + wr(idn)

! prepare the quadratic coefficients (eq. 18 in [1])
!
        a(1) = uh(imx)
        a(2) = - (fh(imx) + uh(ien) + uh(idn))
        a(3) = fh(ien) + fh(idn)

! solve the quadratic equation
!
        nr   = quadratic(a(1:3), x(1:2))

! if Δ < 0, just use the HLL flux
!
        if (nr < 1) then
          f(1:nv,i) = fh(1:nv)
        else

! get the contact dicontinuity speed
!
          sm = x(1)

! if the contact discontinuity speed exceeds the sonic speeds, use the HLL flux
!
          if ((sm <= sl) .or. (sm >= sr)) then
            f(1:nv,i) = fh(1:nv)
          else

! calculate total pressure (eq. 17 in [1])
!
            pr = fh(imx) - (fh(ien) + fh(idn)) * sm

! if the pressure is negative, use the HLL flux
!
            if (pr <= 0.0d+00) then
              f(1:nv,i) = fh(1:nv)
            else

! depending in the sign of the contact dicontinuity speed, calculate the proper
! state and corresponding flux
!
              if (sm > 0.0d+00) then

! calculate the conserved variable vector (eqs. 16 in [1])
!
                dv      = sl - sm
                us(idn) = wl(idn) / dv
                us(imy) = wl(imy) / dv
                us(imz) = wl(imz) / dv
                us(ien) = (wl(ien) + pr * sm) / dv
                us(imx) = (us(ien) + pr) * sm
                us(ien) = us(ien) - us(idn)

! calculate the flux (eq. 14 in [1])
!
                f(1:nv,i) = fl(1:nv,i) + sl * (us(1:nv) - ul(1:nv,i))

              else if (sm < 0.0d+00) then

! calculate the conserved variable vector (eqs. 16 in [1])
!
                dv      = sr - sm
                us(idn) = wr(idn) / dv
                us(imy) = wr(imy) / dv
                us(imz) = wr(imz) / dv
                us(ien) = (wr(ien) + pr * sm) / dv
                us(imx) = (us(ien) + pr) * sm
                us(ien) = us(ien) - us(idn)

! calculate the flux (eq. 14 in [1])
!
                f(1:nv,i) = fr(1:nv,i) + sr * (us(1:nv) - ur(1:nv,i))

              else

! intermediate flux is constant across the contact discontinuity and all fluxes
! except the parallel momentum one are zero
!
                f(idn,i) = 0.0d+00
                f(imx,i) = pr
                f(imy,i) = 0.0d+00
                f(imz,i) = 0.0d+00
                f(ien,i) = 0.0d+00

              end if ! sm == 0

            end if ! p* < 0

          end if ! sl < sm < sr

        end if ! nr < 1

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for the Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_srhd_adi_hllc
!
!===============================================================================
!
!***** ADIABATIC SPECIAL RELATIVITY MAGNETOHYDRODYNAMICS *****
!
!===============================================================================
!
! subroutine UPDATE_FLUX_SRMHD_ADI:
! --------------------------------
!
!   Subroutine solves the Riemann problem along each direction and calculates
!   the numerical fluxes, which are used later to calculate the conserved
!   variable increment.
!
!   Arguments:
!
!     dx   - the spatial step;
!     q    - the array of primitive variables;
!     f    - the array of numerical fluxes;
!
!===============================================================================
!
  subroutine update_flux_srmhd_adi(dx, q, f)

! include external variables
!
    use coordinates    , only : im, jm, km, ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, imx, imy, imz, ipr, ien
    use equations      , only : ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(NDIMS)            , intent(in)  :: dx
    real(kind=8), dimension(      nv,im,jm,km), intent(in)  :: q
    real(kind=8), dimension(NDIMS,nv,im,jm,km), intent(out) :: f

! local variables
!
    integer                        :: i, j, k, l, p
    real(kind=8)                   :: vm

! local temporary arrays
!
    real(kind=8), dimension(nv,im,jm,km)         :: qq
    real(kind=8), dimension(nv,im,jm,km,2,NDIMS) :: qs
    real(kind=8), dimension(nv,im,2)             :: qx
    real(kind=8), dimension(nv,jm,2)             :: qy
    real(kind=8), dimension(nv,km,2)             :: qz
    real(kind=8), dimension(nv,im)               :: fx
    real(kind=8), dimension(nv,jm)               :: fy
    real(kind=8), dimension(nv,km)               :: fz
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux update
!
    call start_timer(imf)
#endif /* PROFILE */

! initialize fluxes
!
    f(1:NDIMS,1:nv,1:im,1:jm,1:km) = 0.0d+00

! apply reconstruction on variables using 4-vector or 3-vector
!
    if (states_4vec) then

! convert velocities to four-velocities for physical reconstruction
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im

            vm = sqrt(1.0d+00 - sum(q(ivx:ivz,i,j,k)**2))

            qq(idn,i,j,k) = q(idn,i,j,k) / vm
            qq(ivx,i,j,k) = q(ivx,i,j,k) / vm
            qq(ivy,i,j,k) = q(ivy,i,j,k) / vm
            qq(ivz,i,j,k) = q(ivz,i,j,k) / vm
            qq(ipr,i,j,k) = q(ipr,i,j,k)
            qq(ibx,i,j,k) = q(ibx,i,j,k)
            qq(iby,i,j,k) = q(iby,i,j,k)
            qq(ibz,i,j,k) = q(ibz,i,j,k)
            qq(ibp,i,j,k) = q(ibp,i,j,k)
            qq(ipr,i,j,k) = q(ipr,i,j,k)

          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

! reconstruct interfaces
!
      call reconstruct_interfaces(dx(:), qq(1:nv,1:im,1:jm,1:km)               &
                                       , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

! convert state four-velocities back to velocities
!
      do k = 1, km
        do j = 1, jm
          do i = 1, im

            do l = 1, 2
              do p = 1, NDIMS

                vm = sqrt(1.0d+00 + sum(qs(ivx:ivz,i,j,k,l,p)**2))

                qs(idn,i,j,k,l,p) = qs(idn,i,j,k,l,p) / vm
                qs(ivx,i,j,k,l,p) = qs(ivx,i,j,k,l,p) / vm
                qs(ivy,i,j,k,l,p) = qs(ivy,i,j,k,l,p) / vm
                qs(ivz,i,j,k,l,p) = qs(ivz,i,j,k,l,p) / vm

              end do ! p = 1, ndims
            end do ! l = 1, 2
          end do ! i = 1, im
        end do ! j = 1, jm
      end do ! k = 1, km

    else

! reconstruct interfaces
!
      call reconstruct_interfaces(dx(:), q (1:nv,1:im,1:jm,1:km)               &
                                       , qs(1:nv,1:im,1:jm,1:km,1:2,1:NDIMS))

    end if

!  calculate the flux along the X-direction
!
    do k = kbl, keu
      do j = jbl, jeu

! copy directional variable vectors to pass to the one dimensional solver
!
        qx(idn,1:im,1:2) = qs(idn,1:im,j,k,1:2,1)
        qx(ivx,1:im,1:2) = qs(ivx,1:im,j,k,1:2,1)
        qx(ivy,1:im,1:2) = qs(ivy,1:im,j,k,1:2,1)
        qx(ivz,1:im,1:2) = qs(ivz,1:im,j,k,1:2,1)
        qx(ibx,1:im,1:2) = qs(ibx,1:im,j,k,1:2,1)
        qx(iby,1:im,1:2) = qs(iby,1:im,j,k,1:2,1)
        qx(ibz,1:im,1:2) = qs(ibz,1:im,j,k,1:2,1)
        qx(ibp,1:im,1:2) = qs(ibp,1:im,j,k,1:2,1)
        qx(ipr,1:im,1:2) = qs(ipr,1:im,j,k,1:2,1)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(im, qx(1:nv,1:im,1), qx(1:nv,1:im,2), fx(1:nv,1:im))

! update the array of fluxes
!
        f(1,idn,1:im,j,k) = fx(idn,1:im)
        f(1,imx,1:im,j,k) = fx(imx,1:im)
        f(1,imy,1:im,j,k) = fx(imy,1:im)
        f(1,imz,1:im,j,k) = fx(imz,1:im)
        f(1,ibx,1:im,j,k) = fx(ibx,1:im)
        f(1,iby,1:im,j,k) = fx(iby,1:im)
        f(1,ibz,1:im,j,k) = fx(ibz,1:im)
        f(1,ibp,1:im,j,k) = fx(ibp,1:im)
        f(1,ien,1:im,j,k) = fx(ien,1:im)

      end do ! j = jbl, jeu
    end do ! k = kbl, keu

!  calculate the flux along the Y direction
!
    do k = kbl, keu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qy(idn,1:jm,1:2) = qs(idn,i,1:jm,k,1:2,2)
        qy(ivx,1:jm,1:2) = qs(ivy,i,1:jm,k,1:2,2)
        qy(ivy,1:jm,1:2) = qs(ivz,i,1:jm,k,1:2,2)
        qy(ivz,1:jm,1:2) = qs(ivx,i,1:jm,k,1:2,2)
        qy(ibx,1:jm,1:2) = qs(iby,i,1:jm,k,1:2,2)
        qy(iby,1:jm,1:2) = qs(ibz,i,1:jm,k,1:2,2)
        qy(ibz,1:jm,1:2) = qs(ibx,i,1:jm,k,1:2,2)
        qy(ibp,1:jm,1:2) = qs(ibp,i,1:jm,k,1:2,2)
        qy(ipr,1:jm,1:2) = qs(ipr,i,1:jm,k,1:2,2)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(jm, qy(1:nv,1:jm,1), qy(1:nv,1:jm,2), fy(1:nv,1:jm))

! update the array of fluxes
!
        f(2,idn,i,1:jm,k) = fy(idn,1:jm)
        f(2,imx,i,1:jm,k) = fy(imz,1:jm)
        f(2,imy,i,1:jm,k) = fy(imx,1:jm)
        f(2,imz,i,1:jm,k) = fy(imy,1:jm)
        f(2,ibx,i,1:jm,k) = fy(ibz,1:jm)
        f(2,iby,i,1:jm,k) = fy(ibx,1:jm)
        f(2,ibz,i,1:jm,k) = fy(iby,1:jm)
        f(2,ibp,i,1:jm,k) = fy(ibp,1:jm)
        f(2,ien,i,1:jm,k) = fy(ien,1:jm)

      end do ! i = ibl, ieu
    end do ! k = kbl, keu

#if NDIMS == 3
!  calculate the flux along the Z direction
!
    do j = jbl, jeu
      do i = ibl, ieu

! copy directional variable vectors to pass to the one dimensional solver
!
        qz(idn,1:km,1:2) = qs(idn,i,j,1:km,1:2,3)
        qz(ivx,1:km,1:2) = qs(ivz,i,j,1:km,1:2,3)
        qz(ivy,1:km,1:2) = qs(ivx,i,j,1:km,1:2,3)
        qz(ivz,1:km,1:2) = qs(ivy,i,j,1:km,1:2,3)
        qz(ibx,1:km,1:2) = qs(ibz,i,j,1:km,1:2,3)
        qz(iby,1:km,1:2) = qs(ibx,i,j,1:km,1:2,3)
        qz(ibz,1:km,1:2) = qs(iby,i,j,1:km,1:2,3)
        qz(ibp,1:km,1:2) = qs(ibp,i,j,1:km,1:2,3)
        qz(ipr,1:km,1:2) = qs(ipr,i,j,1:km,1:2,3)

! call one dimensional Riemann solver in order to obtain numerical fluxes
!
        call riemann(km, qz(1:nv,1:km,1), qz(1:nv,1:km,2), fz(1:nv,1:km))

! update the array of fluxes
!
        f(3,idn,i,j,1:km) = fz(idn,1:km)
        f(3,imx,i,j,1:km) = fz(imy,1:km)
        f(3,imy,i,j,1:km) = fz(imz,1:km)
        f(3,imz,i,j,1:km) = fz(imx,1:km)
        f(3,ibx,i,j,1:km) = fz(iby,1:km)
        f(3,iby,i,j,1:km) = fz(ibz,1:km)
        f(3,ibz,i,j,1:km) = fz(ibx,1:km)
        f(3,ibp,i,j,1:km) = fz(ibp,1:km)
        f(3,ien,i,j,1:km) = fz(ien,1:km)

      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

#ifdef PROFILE
! stop accounting time for flux update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_flux_srmhd_adi
!
!===============================================================================
!
! subroutine RIEMANN_SRMHD_ADI_HLL:
! --------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer (HLL) method.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
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
  subroutine riemann_srmhd_adi_hll(n, ql, qr, f)

! include external procedures
!
    use equations      , only : nv
    use equations      , only : ivx, ibx, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i
    real(kind=8)                  :: sl, sr, srml
    real(kind=8)                  :: bx, bp

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all position
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate the inverse of speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(1:nv)  = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv)  = sr * ur(1:nv,i) - fr(1:nv,i)

! calculate fluxes for the intermediate state
!
        f(1:nv,i) = (sl * wr(1:nv) - sr * wl(1:nv)) / srml

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for the Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_srmhd_adi_hll
!
!===============================================================================
!
! subroutine RIEMANN_SRMHD_ADI_HLLC:
! ---------------------------------
!
!   Subroutine solves one dimensional Riemann problem using
!   the Harten-Lax-van Leer method with contact discontinuity resolution (HLLC)
!   by Mignone & Bodo.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     ql, qr - the array of primitive variables at the Riemann states;
!     f      - the output array of fluxes;
!
!   References:
!
!     [1] Mignone, A. & Bodo, G.
!         "An HLLC Riemann solver for relativistic flows - II.
!          Magnetohydrodynamics",
!         Monthly Notices of the Royal Astronomical Society,
!         2006, Volume 368, Pages 1040-1054
!
!===============================================================================
!
  subroutine riemann_srmhd_adi_hllc(n, ql, qr, f)

! include external procedures
!
    use algebra        , only : quadratic
    use equations      , only : nv
    use equations      , only : ivx, idn, imx, imy, imz, ien, ibx, iby, ibz, ibp
    use equations      , only : cmax
    use equations      , only : prim2cons, fluxspeed

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                      , intent(in)    :: n
    real(kind=8), dimension(nv,n), intent(inout) :: ql, qr
    real(kind=8), dimension(nv,n), intent(out)   :: f

! local variables
!
    integer                       :: i, nr
    real(kind=8)                  :: sl, sr, srml, sm
    real(kind=8)                  :: bx, bp, bx2
    real(kind=8)                  :: pt, dv, fc, uu, ff, uf
    real(kind=8)                  :: vv, vb, gi

! local arrays to store the states
!
    real(kind=8), dimension(nv,n) :: ul, ur, fl, fr
    real(kind=8), dimension(nv)   :: uh, us, fh, wl, wr
    real(kind=8), dimension(n)    :: clm, clp, crm, crp
    real(kind=8), dimension(3)    :: a, vs
    real(kind=8), dimension(2)    :: x
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the Riemann solver
!
    call start_timer(imr)
#endif /* PROFILE */

! obtain the state values for Bx and Psi for the GLM-MHD equations
!
    do i = 1, n

      bx        = 0.5d+00 * ((qr(ibx,i) + ql(ibx,i))                           &
                                             - (qr(ibp,i) - ql(ibp,i)) / cmax)
      bp        = 0.5d+00 * ((qr(ibp,i) + ql(ibp,i))                           &
                                             - (qr(ibx,i) - ql(ibx,i)) * cmax)

      ql(ibx,i) = bx
      qr(ibx,i) = bx
      ql(ibp,i) = bp
      qr(ibp,i) = bp

    end do ! i = 1, n

! calculate the conserved variables of the left and right states
!
    call prim2cons(n, ql(:,:), ul(:,:))
    call prim2cons(n, qr(:,:), ur(:,:))

! calculate the physical fluxes and speeds at both states
!
    call fluxspeed(n, ql(:,:), ul(:,:), fl(:,:), clm(:), clp(:))
    call fluxspeed(n, qr(:,:), ur(:,:), fr(:,:), crm(:), crp(:))

! iterate over all position
!
    do i = 1, n

! estimate the minimum and maximum speeds
!
      sl = min(clm(i), crm(i))
      sr = max(clp(i), crp(i))

! calculate the HLL flux
!
      if (sl >= 0.0d+00) then

        f(1:nv,i) = fl(1:nv,i)

      else if (sr <= 0.0d+00) then

        f(1:nv,i) = fr(1:nv,i)

      else ! sl < 0 < sr

! calculate the inverse of speed difference
!
        srml = sr - sl

! calculate vectors of the left and right-going waves
!
        wl(1:nv)  = sl * ul(1:nv,i) - fl(1:nv,i)
        wr(1:nv)  = sr * ur(1:nv,i) - fr(1:nv,i)

! calculate fluxes for the intermediate state
!
        uh(1:nv)  = (     wr(1:nv) -      wl(1:nv)) / srml
        fh(1:nv)  = (sl * wr(1:nv) - sr * wl(1:nv)) / srml

! correct the energy waves
!
        wl(ien)   = wl(ien) + wl(idn)
        wr(ien)   = wr(ien) + wr(idn)

! calculate Bₓ²
!
        bx2       = ql(ibx,i) * qr(ibx,i)

! calculate the contact dicontinuity speed solving eq. (41)
!
        if (bx2 >= 1.0d-16) then

! prepare the quadratic coefficients, (eq. 42 in [1])
!
          uu = sum(uh(iby:ibz) * uh(iby:ibz))
          ff = sum(fh(iby:ibz) * fh(iby:ibz))
          uf = sum(uh(iby:ibz) * fh(iby:ibz))

          a(1) = uh(imx) - uf
          a(2) = uu + ff - (fh(imx) + uh(ien) + uh(idn))
          a(3) = fh(ien) + fh(idn) - uf

! solve the quadratic equation
!
          nr   = quadratic(a(1:3), x(1:2))

! if Δ < 0, just use the HLL flux
!
          if (nr < 1) then
            f(1:nv,i) = fh(1:nv)
          else

! get the contact dicontinuity speed
!
            sm = x(1)

! if the contact discontinuity speed exceeds the sonic speeds, use the HLL flux
!
            if ((sm <= sl) .or. (sm >= sr)) then
              f(1:nv,i) = fh(1:nv)
            else

! substitute magnetic field components from the HLL state (eq. 37 in [1])
!
              us(ibx) = ql(ibx,i)
              us(iby) = uh(iby)
              us(ibz) = uh(ibz)
              us(ibp) = ql(ibp,i)

! calculate velocity components without Bₓ normalization (eq. 38 in [1])
!
              vs(1)   = sm
              vs(2)   = (us(iby) * sm - fh(iby)) / us(ibx)
              vs(3)   = (us(ibz) * sm - fh(ibz)) / us(ibx)

! calculate v² and v.B
!
              vv      = sum(vs(1:3) * vs(  1:  3))
              vb      = sum(vs(1:3) * us(ibx:ibz))

! calculate inverse gamma
!
              gi      = 1.0d+00 - vv

! calculate total pressure (eq. 40 in [1])
!
              pt = fh(imx) + gi * bx2 - (fh(ien) + fh(idn) - us(ibx) * vb) * sm

! if the pressure is negative, use the HLL flux
!
              if (pt <= 0.0d+00) then
                f(1:nv,i) = fh(1:nv)
              else

! depending in the sign of the contact dicontinuity speed, calculate the proper
! state and corresponding flux
!
                if (sm > 0.0d+00) then

! calculate the conserved variable vector (eqs. 43-46 in [1])
!
                  dv      = sl - sm
                  fc      = (sl - ql(ivx,i)) / dv
                  us(idn) = fc * ul(idn,i)
                  us(imy) = (sl * ul(imy,i) - fl(imy,i)                        &
                          - us(ibx) * (gi * us(iby) + vs(2) * vb)) / dv
                  us(imz) = (sl * ul(imz,i) - fl(imz,i)                        &
                          - us(ibx) * (gi * us(ibz) + vs(3) * vb)) / dv
                  us(ien) = (sl * (ul(ien,i) + ul(idn,i))                      &
                          - (fl(ien,i) + fl(idn,i))                            &
                          + pt * sm - us(ibx) * vb) / dv - us(idn)

! calculate normal component of momentum
!
                  us(imx) = (us(ien) + us(idn) + pt) * sm - us(ibx) * vb

! calculate the flux (eq. 14 in [1])
!
                  f(1:nv,i) = fl(1:nv,i) + sl * (us(1:nv) - ul(1:nv,i))

                else if (sm < 0.0d+00) then

! calculate the conserved variable vector (eqs. 43-46 in [1])
!
                  dv      = sr - sm
                  fc      = (sr - qr(ivx,i)) / dv
                  us(idn) = fc * ur(idn,i)
                  us(imy) = (sr * ur(imy,i) - fr(imy,i)                        &
                          - us(ibx) * (gi * us(iby) + vs(2) * vb)) / dv
                  us(imz) = (sr * ur(imz,i) - fr(imz,i)                        &
                          - us(ibx) * (gi * us(ibz) + vs(3) * vb)) / dv
                  us(ien) = (sr * (ur(ien,i) + ur(idn,i))                      &
                          - (fr(ien,i) + fr(idn,i))                            &
                          + pt * sm - us(ibx) * vb) / dv - us(idn)

! calculate normal component of momentum
!
                  us(imx) = (us(ien) + us(idn) + pt) * sm - us(ibx) * vb

! calculate the flux (eq. 14 in [1])
!
                  f(1:nv,i) = fr(1:nv,i) + sr * (us(1:nv) - ur(1:nv,i))

                else

! intermediate flux is constant across the contact discontinuity so use HLL flux
!
                  f(1:nv,i) = fh(1:nv)

                end if ! sm == 0

              end if ! p* < 0

            end if ! sl < sm < sr

          end if ! nr < 1

        else ! Bx² > 0

! prepare the quadratic coefficients (eq. 47 in [1])
!
          a(1) = uh(imx)
          a(2) = - (fh(imx) + uh(ien) + uh(idn))
          a(3) = fh(ien) + fh(idn)

! solve the quadratic equation
!
          nr   = quadratic(a(1:3), x(1:2))

! if Δ < 0, just use the HLL flux
!
          if (nr < 1) then
            f(1:nv,i) = fh(1:nv)
          else

! get the contact dicontinuity speed
!
            sm = x(1)

! if the contact discontinuity speed exceeds the sonic speeds, use the HLL flux
!
            if ((sm <= sl) .or. (sm >= sr)) then
              f(1:nv,i) = fh(1:nv)
            else

! calculate total pressure (eq. 48 in [1])
!
              pt = fh(imx) - (fh(ien) + fh(idn)) * sm

! if the pressure is negative, use the HLL flux
!
              if (pt <= 0.0d+00) then
                f(1:nv,i) = fh(1:nv)
              else

! depending in the sign of the contact dicontinuity speed, calculate the proper
! state and corresponding flux
!
                if (sm > 0.0d+00) then

! calculate the conserved variable vector (eqs. 43-46 in [1])
!
                  dv      = sl - sm
                  us(idn) = wl(idn) / dv
                  us(imy) = wl(imy) / dv
                  us(imz) = wl(imz) / dv
                  us(ien) = (wl(ien) + pt * sm) / dv
                  us(imx) = (us(ien) + pt) * sm
                  us(ien) = us(ien) - us(idn)
                  us(ibx) = 0.0d+00
                  us(iby) = wl(iby) / dv
                  us(ibz) = wl(ibz) / dv
                  us(ibp) = ql(ibp,i)

! calculate the flux (eq. 27 in [1])
!
                  f(1:nv,i) = fl(1:nv,i) + sl * (us(1:nv) - ul(1:nv,i))

                else if (sm < 0.0d+00) then

! calculate the conserved variable vector (eqs. 43-46 in [1])
!
                  dv      = sr - sm
                  us(idn) = wr(idn) / dv
                  us(imy) = wr(imy) / dv
                  us(imz) = wr(imz) / dv
                  us(ien) = (wr(ien) + pt * sm) / dv
                  us(imx) = (us(ien) + pt) * sm
                  us(ien) = us(ien) - us(idn)
                  us(ibx) = 0.0d+00
                  us(iby) = wr(iby) / dv
                  us(ibz) = wr(ibz) / dv
                  us(ibp) = qr(ibp,i)

! calculate the flux (eq. 27 in [1])
!
                  f(1:nv,i) = fr(1:nv,i) + sr * (us(1:nv) - ur(1:nv,i))

                else

! intermediate flux is constant across the contact discontinuity so use HLL flux
!
                  f(1:nv,i) = fh(1:nv)

                end if ! sm == 0

              end if ! p* < 0

            end if ! sl < sm < sr

          end if ! nr < 1

        end if ! Bx² == 0 or > 0

      end if ! sl < 0 < sr

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for the Riemann solver
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine riemann_srmhd_adi_hllc

!===============================================================================
!
end module schemes
