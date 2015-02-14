!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2015 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: EQUATIONS
!!
!!  This module provides interface for the systems of equations.  Any set of
!!  equations gives us some basic informations, such as the number of variables,
!!  the primitive and conservative variable definitions, the conversion between
!!  those variables, the flux and characteristic speeds defined in terms of
!!  primitive variables.  All this information is provided by this module.
!!
!!  In order to implement a new set of equations, we need to:
!!
!!  1) define the number of independent variables (or equations) nv;
!!  2) define the variable indices and names (both primitive and conservative);
!!  3) provide subroutines for primitive-conservative variable conversion and
!!     point them to corresponding pointers;
!!  4) provide a subroutine to calculate physical fluxes and characteristic
!!     speeds;
!!  5) provide a subroutine to calculate the maximum speed;
!!  6) optionally, define and read all physical constants related to a given
!!     system;
!!
!!******************************************************************************
!
module equations

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
  integer            , save :: imi, imc, imf, imm
#endif /* PROFILE */

! pointers to the conversion procedures
!
  procedure(prim2cons_hd_iso)  , pointer, save :: prim2cons       => null()
  procedure(cons2prim_hd_iso)  , pointer, save :: cons2prim       => null()

! pointer to the flux procedure
!
  procedure(fluxspeed_hd_iso)  , pointer, save :: fluxspeed       => null()

! pointer to the maxspeed procedure
!
  procedure(maxspeed_hd_iso)   , pointer, save :: maxspeed        => null()

! pointer to the Roe eigensystem procedure
!
  procedure(esystem_roe_hd_iso), pointer, save :: eigensystem_roe => null()

! pointer to the variable conversion method
!
  procedure(nr_iterate_srhd_adi_1dw), pointer, save :: nr_iterate => null()

! the system of equations and the equation of state
!
  character(len=32), save :: eqsys = "hydrodynamic"
  character(len=32), save :: eos   = "adiabatic"

! the variable conversion method
!
  character(len=32), save :: c2p   = "1Dw"

! the number of independent variables
!
  integer(kind=4)  , save :: nv  = 0

! direction indices
!
  integer(kind=4)  , save :: inx =  1, iny =  2, inz =  3

! variable indices
!
  integer(kind=4)  , save :: idn = -1
  integer(kind=4)  , save :: ivx = -1, ivy = -1, ivz = -1
  integer(kind=4)  , save :: imx = -1, imy = -1, imz = -1
  integer(kind=4)  , save :: ibx = -1, iby = -1, ibz = -1
  integer(kind=4)  , save :: ibp = -1
  integer(kind=4)  , save :: ipr = -1, ien = -1

! variable names
!
  character(len=4), dimension(:), allocatable, save :: pvars, cvars

! eigenvectors
!
  real(kind=8), dimension(:,:,:), allocatable, save :: evroe

! adiabatic heat ratio
!
  real(kind=8)     , save :: gamma   = 5.0d+00 / 3.0d+00

! additional adiabatic parameters
!
  real(kind=8)     , save :: gammam1 = 2.0d+00 / 3.0d+00, gammam1i = 1.5d+00
  real(kind=8)     , save :: gammaxi = 2.0d+00 / 5.0d+00

! isothermal speed of sound and its second power
!
  real(kind=8)     , save :: csnd    = 1.0d+00, csnd2   = 1.0d+00

! maximum speed in the system
!
  real(kind=8)     , save :: cmax    = 0.0d+00, cmax2   = 0.0d+00

! the lower limits for density and pressure to be treated as physical
!
  real(kind=8)     , save :: dmin    = 1.0d-08, pmin    = 1.0d-08

! the upper bound for the sonic Mach number
!
  real(kind=8)     , save :: msmax   = 1.0d+03
  real(kind=8)     , save :: msfac   = 3.0d-06 / 5.0d+00

! the tolerance for Newton-Raphson interative method, the maximum number of
! iterations and the number of extra iterations for polishing
!
  real(kind=8)     , save :: tol     = 1.0d-10
  integer          , save :: nmax    = 100
  integer          , save :: next    = 2

! flags for reconstruction corrections
!
  logical          , save :: positivity = .false.

! by default everything is private
!
  private

! declare public variables and subroutines
!
  public :: initialize_equations, finalize_equations
  public :: prim2cons, cons2prim
  public :: fluxspeed
  public :: maxspeed, reset_maxspeed, get_maxspeed
  public :: eigensystem_roe
  public :: update_primitive_variables
  public :: gamma
  public :: csnd, csnd2
  public :: cmax, cmax2
  public :: nv
  public :: inx, iny, inz
  public :: idn, ivx, ivy, ivz, imx, imy, imz
  public :: ibx, iby, ibz, ibp, ipr, ien
  public :: eqsys, eos
  public :: pvars, cvars

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
! subroutine INITIALIZE_EQUATIONS:
! -------------------------------
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
  subroutine initialize_equations(verbose, iret)

! include external procedures and variables
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
    logical                :: relativistic   = .false.
    character(len=255)     :: name_eqsys     = ""
    character(len=255)     :: name_eos       = ""
    character(len=255)     :: name_c2p       = ""
    character(len=255)     :: positivity_fix = "off"
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('equations:: initialization'     , imi)
    call set_timer('equations:: variable conversion', imc)
    call set_timer('equations:: flux calculation'   , imf)
    call set_timer('equations:: speed estimation'   , imm)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the system of equations
!
    call get_parameter_string("equation_system"   , eqsys)

! get the equation of state
!
    call get_parameter_string("equation_of_state" , eos  )

! get the primitive variable solver
!
    call get_parameter_string("primitive_solver"  , c2p  )

! get the tolerance
!
    call get_parameter_real   ("tolerance"         , tol )

! get the maximum number of Newton-Raphson method iterations
!
    call get_parameter_integer("maximum_iterations", nmax)
    call get_parameter_integer("extra_iterations"  , next)

! depending on the system of equations initialize the module variables
!
    select case(trim(eqsys))

!--- HYDRODYNAMICS ---
!
    case("hd", "HD", "hydro", "HYDRO", "hydrodynamic", "HYDRODYNAMIC")

! the name of equation system
!
      name_eqsys = "HD"

! initialize the number of variables (density + 3 components of velocity)
!
      nv  = 4

! initialize the variable indices
!
      idn = 1
      ivx = 2
      ivy = 3
      ivz = 4
      imx = 2
      imy = 3
      imz = 4

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("iso", "ISO", "isothermal", "ISOTHERMAL")

! the type of equation of state
!
        name_eos  = "isothermal"

! set pointers to subroutines
!
        prim2cons       => prim2cons_hd_iso
        cons2prim       => cons2prim_hd_iso
        fluxspeed       => fluxspeed_hd_iso
        maxspeed        => maxspeed_hd_iso
        eigensystem_roe => esystem_roe_hd_iso

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! the type of equation of state
!
        name_eos  = "adiabatic"

! include the pressure/energy in the number of variables
!
        nv  = nv + 1

! initialize the pressure and energy indices
!
        ipr = nv
        ien = nv

! set pointers to subroutines
!
        prim2cons       => prim2cons_hd_adi
        cons2prim       => cons2prim_hd_adi
        fluxspeed       => fluxspeed_hd_adi
        maxspeed        => maxspeed_hd_adi
        eigensystem_roe => esystem_roe_hd_adi

! warn about the unimplemented equation of state
!
      case default

        if (verbose) then
          write (*,"(1x,a)") "The selected equation of state is not " //       &
                             "implemented: " // trim(eos)
          write (*,*)
        end if
        iret = 110
        return

      end select

! allocate arrays for variable names
!
      allocate(pvars(nv), cvars(nv))

! fill in the primitive variable names
!
      pvars(idn) = 'dens'
      pvars(ivx) = 'velx'
      pvars(ivy) = 'vely'
      pvars(ivz) = 'velz'
      if (ipr > 0) pvars(ipr) = 'pres'

! fill in the conservative variable names
!
      cvars(idn) = 'dens'
      cvars(imx) = 'momx'
      cvars(imy) = 'momy'
      cvars(imz) = 'momz'
      if (ien > 0) cvars(ien) = 'ener'

!--- MAGNETOHYDRODYNAMICS ---
!
    case("mhd", "MHD", "magnetohydrodynamic", "MAGNETOHYDRODYNAMIC")

! the name of equation system
!
      name_eqsys = "MHD"

! initialize the number of variables (density + 3 components of velocity
!                                             + 3 components of magnetic field)
!
      nv  = 8

! initialize the variable indices
!
      idn = 1
      ivx = 2
      ivy = 3
      ivz = 4
      imx = 2
      imy = 3
      imz = 4
      ibx = 5
      iby = 6
      ibz = 7
      ibp = 8

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("iso", "ISO", "isothermal", "ISOTHERMAL")

! the type of equation of state
!
        name_eos  = "isothermal"

! set pointers to the subroutines
!
        prim2cons       => prim2cons_mhd_iso
        cons2prim       => cons2prim_mhd_iso
        fluxspeed       => fluxspeed_mhd_iso
        maxspeed        => maxspeed_mhd_iso
        eigensystem_roe => esystem_roe_mhd_iso

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! the type of equation of state
!
        name_eos  = "adiabatic"

! increase the number of variables by the pressure/energy
!
        nv  = nv + 1

! initialize the pressure and energy indices
!
        ipr = nv
        ien = nv

! set pointers to subroutines
!
        prim2cons       => prim2cons_mhd_adi
        cons2prim       => cons2prim_mhd_adi
        fluxspeed       => fluxspeed_mhd_adi
        maxspeed        => maxspeed_mhd_adi
        eigensystem_roe => esystem_roe_mhd_adi

      case default

        if (verbose) then
          write (*,"(1x,a)") "The selected equation of state is not " //       &
                             "implemented: " // trim(eos)
          write (*,*)
        end if
        iret = 110
        return

      end select

! allocate arrays for variable names
!
      allocate(pvars(nv), cvars(nv))

! fill in the primitive variable names
!
      pvars(idn) = 'dens'
      pvars(ivx) = 'velx'
      pvars(ivy) = 'vely'
      pvars(ivz) = 'velz'
      pvars(ibx) = 'magx'
      pvars(iby) = 'magy'
      pvars(ibz) = 'magz'
      pvars(ibp) = 'bpsi'
      if (ipr > 0) pvars(ipr) = 'pres'

! fill in the conservative variable names
!
      cvars(idn) = 'dens'
      cvars(imx) = 'momx'
      cvars(imy) = 'momy'
      cvars(imz) = 'momz'
      cvars(ibx) = 'magx'
      cvars(iby) = 'magy'
      cvars(ibz) = 'magz'
      cvars(ibp) = 'bpsi'
      if (ien > 0) cvars(ien) = 'ener'

!--- SPECIAL RELATIVITY HYDRODYNAMICS ---
!
    case("srhd", "SRHD")

! the name of equation system
!
      name_eqsys = "Special Relativity HD"

! set relativistic flag
!
      relativistic = .true.

! initialize the number of variables (density + 3 components of velocity)
!
      nv  = 4

! initialize the variable indices
!
      idn = 1
      ivx = 2
      ivy = 3
      ivz = 4
      imx = 2
      imy = 3
      imz = 4

! depending on the equation of state complete the initialization
!
      select case(trim(eos))

      case("adi", "ADI", "adiabatic", "ADIABATIC")

! the type of equation of state
!
        name_eos  = "adiabatic"

! include the pressure/energy in the number of variables
!
        nv  = nv + 1

! initialize the pressure and energy indices
!
        ipr = nv
        ien = nv

! set pointers to subroutines
!
        prim2cons       => prim2cons_srhd_adi
        cons2prim       => cons2prim_srhd_adi
        fluxspeed       => fluxspeed_srhd_adi
        maxspeed        => maxspeed_srhd_adi

! warn about the unimplemented equation of state
!
      case default

        if (verbose) then
          write (*,"(1x,a)") "The selected equation of state is not " //       &
                             "implemented: " // trim(eos)
          write (*,*)
        end if
        iret = 110
        return

      end select

! choose the conserved to primitive variable conversion method
!
      select case(trim(c2p))

      case("1Dw", "1dw", "1DW")

! the type of equation of state
!
        name_c2p  = "1Dw"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srhd_adi_1dw

      case("2D", "2d")

! the type of equation of state
!
        name_c2p  = "2D"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srhd_adi_2d

! warn about the unimplemented method
!
      case default

        if (verbose) then
          write (*,"(1x,a)") "The selected conversion method is not " //       &
                             "implemented: " // trim(c2p)
          write (*,*)
        end if
        iret = 110
        return

      end select

! allocate arrays for variable names
!
      allocate(pvars(nv), cvars(nv))

! fill in the primitive variable names
!
      pvars(idn) = 'dens'
      pvars(ivx) = 'velx'
      pvars(ivy) = 'vely'
      pvars(ivz) = 'velz'
      if (ipr > 0) pvars(ipr) = 'pres'

! fill in the conservative variable names
!
      cvars(idn) = 'dens'
      cvars(imx) = 'momx'
      cvars(imy) = 'momy'
      cvars(imz) = 'momz'
      if (ien > 0) cvars(ien) = 'ener'

!--- EQUATION SYSTEM NOT IMPLEMENTED ---
!
    case default

      if (verbose) then
        write (*,"(1x,a)") "The selected equation system is not " //           &
                           "implemented: " // trim(eqsys)
        write (*,*)
      end if
      iret = 100
      return

    end select

! obtain the adiabatic specific heat ratio
!
    call get_parameter_real("gamma" , gamma )

! calculate additional parameters
!
    gammam1  = gamma - 1.0d+00
    gammam1i = 1.0d+00 / gammam1
    gammaxi  = gammam1 / gamma

! obtain the isothermal sound speed
!
    call get_parameter_real("csnd"  , csnd  )

! calculate additional parameters
!
    csnd2 = csnd * csnd

! allocate space for Roe eigenvectors
!
    allocate(evroe(2,nv,nv))

! get the minimum allowed density and pressure in the system
!
    call get_parameter_real("dmin"  , dmin  )
    call get_parameter_real("pmin"  , pmin  )

! get the upper bound for the sonic Mach number
!
    call get_parameter_real("msmax" , msmax )

! calculate the sonic Mach number factor
!
    msfac = 1.0d+00 / (gamma * msmax**2)

! get the positivity fix flag
!
    call get_parameter_string("fix_positivity", positivity_fix )

! check additional reconstruction limiting
!
    select case(trim(positivity_fix))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      positivity = .true.
    case default
      positivity = .false.
    end select

! print information about the equation module
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "equation system        =", trim(name_eqsys)
      write (*,"(4x,a,1x,a)"    ) "equation of state      =", trim(name_eos)
      if (relativistic) then
        write (*,"(4x,a,1x,a)"    ) "variable conversion    =", trim(name_c2p)
      end if

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_equations
!
!===============================================================================
!
! subroutine FINALIZE_EQUATIONS:
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
  subroutine finalize_equations(iret)

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

! deallocate variable name arrays
!
    if (allocated(pvars)) deallocate(pvars)
    if (allocated(cvars)) deallocate(cvars)

! deallocate Roe eigenvectors
!
    if (allocated(evroe)) deallocate(evroe)

! release the procedure pointers
!
    nullify(prim2cons)
    nullify(cons2prim)
    nullify(fluxspeed)
    nullify(maxspeed)
    nullify(eigensystem_roe)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_equations
!
!===============================================================================
!
! subroutine RESET_MAXSPEED:
! -------------------------
!
!   Subroutine resets the maximum speed in the domain to zero.
!
!
!===============================================================================
!
  subroutine reset_maxspeed()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed
!
    cmax = 0.0d+00

!-------------------------------------------------------------------------------
!
  end subroutine reset_maxspeed
!
!===============================================================================
!
! function GET_MAXSPEED:
! -----------------
!
!   Function returns the maximum speed in the domain.
!
!
!===============================================================================
!
  real(kind=8) function get_maxspeed()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! return the maximum speed
!
    get_maxspeed = cmax

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function get_maxspeed
!
!===============================================================================
!
! subroutine UPDATE_PRIMITIVE_VARIABLES:
! -------------------------------------
!
!   Subroutine updates primitive variables from their conservative
!   representation.  This process is done once after advance of the conserved
!   variables due to their evolution in time.
!
!   Arguments:
!
!     uu - the input array of conservative variables;
!     qq - the output array of primitive variables;
!
!===============================================================================
!
  subroutine update_primitive_variables(uu, qq)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, in, jn, kn
    use coordinates, only : ib, jb, kb, ie, je, ke

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: uu
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: qq

! temporary variables
!
    integer      :: i, j, k
    real(kind=8) :: pmin
!
!-------------------------------------------------------------------------------
!
! update primitive variables
!
    do k = kb, ke
      do j = jb, je

! convert conserved variables to primitive ones
!
        call cons2prim(in, uu(1:nv,ib:ie,j,k), qq(1:nv,ib:ie,j,k))

      end do ! j = jb, je
    end do ! k = kb, ke

! fix negative pressure is desired
!
    if (positivity .and. ipr > 0) then

! iterate over block interior
!
      do k = kb, ke
        do j = jb, je
          do i = ib, ie

! fix the cells where pressure is negative
!
            if (qq(ipr,i,j,k) <= 0.0d+00) then

! calculate pressure from the sonic Mach number limit and local velocity
!
              pmin          = msfac * qq(idn,i,j,k) * sum(qq(ivx:ivz,i,j,k)**2)

! update total energy and pressure
!
              uu(ien,i,j,k) = uu(ien,i,j,k) + gammam1i * (pmin - qq(ipr,i,j,k))
              qq(ipr,i,j,k) = pmin

            end if

          end do ! i = ib, ie
        end do ! j = jb, je
      end do ! k = kb, ke

    end if

!-------------------------------------------------------------------------------
!
  end subroutine update_primitive_variables
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!*******************************************************************************
!
! ISOTHERMAL HYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_HD_ISO:
! ---------------------------
!
!   Subroutine converts primitive variables to their corresponding
!   conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons_hd_iso(n, q, u)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: u

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_hd_iso
!
!===============================================================================
!
! subroutine CONS2PRIM_HD_ISO:
! ---------------------------
!
!   Subroutine converts conservative variables to their corresponding
!   primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim_hd_iso(n, u, q)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: u
    real(kind=8), dimension(nv,n), intent(out) :: q

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_hd_iso
!
!===============================================================================
!
! subroutine FLUXSPEED_HD_ISO:
! ---------------------------
!
!   Subroutine calculates physical fluxes and characteristic speeds from a
!   given equation system.
!
!   Arguments:
!
!     n      - the length of input and output vectors;
!     q      - the input array of primitive variables;
!     u      - the input array of conservative variables;
!     f      - the output vector of fluxes;
!     cm, cp - the output vector of left- and right-going characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_hd_iso(n, q, u, f, cm, cp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                , intent(in)  :: n
    real(kind=8), dimension(nv,n)          , intent(in)  :: q, u
    real(kind=8), dimension(nv,n)          , intent(out) :: f
    real(kind=8), dimension(n)   , optional, intent(out) :: cm, cp

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux calculation
!
    call start_timer(imf)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! calculate the hydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i)
      f(imy,i) = q(ivx,i) * u(imy,i)
      f(imz,i) = q(ivx,i) * u(imz,i)
      f(imx,i) = f(imx,i) + csnd2 * q(idn,i)

! calculate characteristic speeds
!
      cm(i)    = q(ivx,i) - csnd
      cp(i)    = q(ivx,i) + csnd

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_hd_iso
!
!===============================================================================
!
! function MAXSPEED_HD_ISO:
! ------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!
!===============================================================================
!
  function maxspeed_hd_iso(qq) result(maxspeed)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in) :: qq

! return value
!
    real(kind=8) :: maxspeed

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: vv, v, c
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the maximum speed estimation
!
    call start_timer(imm)
#endif /* PROFILE */

! reset the maximum speed
!
    maxspeed = 0.0d+00

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(qq(ivx:ivz,i,j,k) * qq(ivx:ivz,i,j,k))
          v  = sqrt(vv)

! calculate the maximum speed
!
          maxspeed = max(maxspeed, v + csnd)

        end do ! i = ib, ie
      end do ! j = jb, je
    end do ! k = kb, ke

#ifdef PROFILE
! stop accounting time for the maximum speed estimation
!
    call stop_timer(imm)
#endif /* PROFILE */

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_hd_iso
!
!===============================================================================
!
! subroutine ESYSTEM_ROE_HD_ISO:
! -----------------------------
!
!   Subroutine computes eigenvalues and eigenvectors for a given set of
!   equations and input variables.
!
!   Arguments:
!
!     x - ratio of the perpendicular magnetic field component difference
!     y - ratio of the density
!     q - the intermediate Roe state vector;
!     c - the vector of eigenvalues;
!     r - the matrix of right eigenvectors;
!     l - the matrix of left eigenvectors;
!
!   References:
!
!     [1] Roe, P. L.
!         "Approximate Riemann Solvers, Parameter Vectors, and Difference
!          Schemes",
!         Journal of Computational Physics, 1981, 43, pp. 357-372
!     [2] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!
!===============================================================================
!
  subroutine esystem_roe_hd_iso(x, y, q, c, r, l)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)                  , intent(in)    :: x, y
    real(kind=8), dimension(nv)   , intent(in)    :: q
    real(kind=8), dimension(nv)   , intent(inout) :: c
    real(kind=8), dimension(nv,nv), intent(inout) :: l, r

! local variables
!
    logical     , save :: first = .true.
    real(kind=8), save :: ch
!
!-------------------------------------------------------------------------------
!
! prepare the internal arrays at the first run
!
    if (first) then

! prepare constants
!
      ch      = 0.5d+00 / csnd

! reset all elements
!
      evroe(:, : ,:) = 0.0d+00

! initiate the matrix of left eigenvectors
!
      evroe(1,ivx,1) = - ch
      evroe(1,ivy,2) = 1.0d+00
      evroe(1,ivz,3) = 1.0d+00
      evroe(1,ivx,4) =   ch

! initiate the matrix of right eigenvectors
!
      evroe(2,1,idn) = 1.0d+00
      evroe(2,2,ivy) = 1.0d+00
      evroe(2,3,ivz) = 1.0d+00
      evroe(2,4,idn) = 1.0d+00

! unset the first execution flag
!
      first = .false.

    end if ! first execution

! prepare eigenvalues
!
    c(1)           = q(ivx) - csnd
    c(2)           = q(ivx)
    c(3)           = q(ivx)
    c(4)           = q(ivx) + csnd

! update the varying elements of the matrix of left eigenvectors
!
    evroe(1,idn,1) =   ch * c(4)
    evroe(1,idn,2) = - q(ivy)
    evroe(1,idn,3) = - q(ivz)
    evroe(1,idn,4) = - ch * c(1)

! update the varying elements of the matrix of right eigenvectors
!
    evroe(2,1,ivx) = c(1)
    evroe(2,1,ivy) = q(ivy)
    evroe(2,1,ivz) = q(ivz)

    evroe(2,4,ivx) = c(4)
    evroe(2,4,ivy) = q(ivy)
    evroe(2,4,ivz) = q(ivz)

! copy matrices of eigenvectors
!
    l(1:nv,1:nv) = evroe(1,1:nv,1:nv)
    r(1:nv,1:nv) = evroe(2,1:nv,1:nv)

!-------------------------------------------------------------------------------
!
  end subroutine esystem_roe_hd_iso
!
!*******************************************************************************
!
! ADIABATIC HYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_HD_ADI:
! ---------------------------
!
!   Subroutine converts primitive variables to their corresponding
!   conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons_hd_adi(n, q, u)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: u

! local variables
!
    integer      :: i
    real(kind=8) :: ek, ei
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
      ek       = 0.5d+00 * (u(imx,i) * q(ivx,i) + u(imy,i) * q(ivy,i)          &
                                                + u(imz,i) * q(ivz,i))
      ei       = gammam1i * q(ipr,i)
      u(ien,i) = ei + ek

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_hd_adi
!
!===============================================================================
!
! subroutine CONS2PRIM_HD_ADI:
! ---------------------------
!
!   Subroutine converts conservative variables to their corresponding
!   primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim_hd_adi(n, u, q)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: u
    real(kind=8), dimension(nv,n), intent(out) :: q

! local variables
!
    integer      :: i
    real(kind=8) :: ek, ei
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
      ek       = 0.5d+00 * (u(imx,i) * q(ivx,i) + u(imy,i) * q(ivy,i)          &
                                                + u(imz,i) * q(ivz,i))
      ei       = u(ien,i) - ek
      q(ipr,i) = gammam1 * ei

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_hd_adi
!
!===============================================================================
!
! subroutine FLUXSPEED_HD_ADI:
! ---------------------------
!
!   Subroutine calculates physical fluxes and characteristic speeds from a
!   given equation system.
!
!   Arguments:
!
!     n      - the length of input and output vectors;
!     q      - the input array of primitive variables;
!     u      - the input array of conservative variables;
!     f      - the output vector of fluxes;
!     cm, cp - the output vector of left- and right-going characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_hd_adi(n, q, u, f, cm, cp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                , intent(in)  :: n
    real(kind=8), dimension(nv,n)          , intent(in)  :: q, u
    real(kind=8), dimension(nv,n)          , intent(out) :: f
    real(kind=8), dimension(n)   , optional, intent(out) :: cm, cp

! local variables
!
    integer      :: i
    real(kind=8) :: cs
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux calculation
!
    call start_timer(imf)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! calculate the hydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i)
      f(imy,i) = q(ivx,i) * u(imy,i)
      f(imz,i) = q(ivx,i) * u(imz,i)
      f(imx,i) = f(imx,i) + q(ipr,i)
      f(ien,i) = q(ivx,i) * (u(ien,i) + q(ipr,i))

! calculate characteristic speeds
!
      cs       = sqrt(gamma * q(ipr,i) / q(idn,i))
      cm(i)    = q(ivx,i) - cs
      cp(i)    = q(ivx,i) + cs

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_hd_adi
!
!===============================================================================
!
! function MAXSPEED_HD_ADI:
! ------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!
!===============================================================================
!
  function maxspeed_hd_adi(qq) result(maxspeed)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in) :: qq

! return value
!
    real(kind=8) :: maxspeed

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: vv, v, c
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the maximum speed estimation
!
    call start_timer(imm)
#endif /* PROFILE */

! reset the maximum speed
!
    maxspeed = 0.0d+00

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(qq(ivx:ivz,i,j,k) * qq(ivx:ivz,i,j,k))
          v  = sqrt(vv)

! calculate the adiabatic speed of sound
!
          c = sqrt(gamma * qq(ipr,i,j,k) / qq(idn,i,j,k))

! calculate the maximum speed
!
          maxspeed = max(maxspeed, v + c)

        end do ! i = ib, ie
      end do ! j = jb, je
    end do ! k = kb, ke

#ifdef PROFILE
! stop accounting time for the maximum speed estimation
!
    call stop_timer(imm)
#endif /* PROFILE */

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_hd_adi
!
!===============================================================================
!
! subroutine ESYSTEM_ROE_HD_ADI:
! -----------------------------
!
!   Subroutine computes eigenvalues and eigenvectors for a given set of
!   equations and input variables.
!
!   Arguments:
!
!     x - ratio of the perpendicular magnetic field component difference
!     y - ratio of the density
!     q - the intermediate Roe state vector;
!     c - the vector of eigenvalues;
!     r - the matrix of right eigenvectors;
!     l - the matrix of left eigenvectors;
!
!   References:
!
!     [1] Roe, P. L.
!         "Approximate Riemann Solvers, Parameter Vectors, and Difference
!          Schemes",
!         Journal of Computational Physics, 1981, 43, pp. 357-372
!     [2] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!
!===============================================================================
!
  subroutine esystem_roe_hd_adi(x, y, q, c, r, l)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)                  , intent(in)    :: x, y
    real(kind=8), dimension(nv)   , intent(in)    :: q
    real(kind=8), dimension(nv)   , intent(inout) :: c
    real(kind=8), dimension(nv,nv), intent(inout) :: l, r

! local variables
!
    logical, save :: first = .true.
    real(kind=8)  :: vv, vh, c2, na, cc, vc, ng, nd, nw, nh, nc
!
!-------------------------------------------------------------------------------
!
! prepare the internal arrays at the first run
!
    if (first) then

! reset all elements
!
      evroe(:, : ,:) = 0.0d+00

! initiate the matrix of left eigenvectors
!
      evroe(1,ivy,2) = 1.0d+00
      evroe(1,ivz,3) = 1.0d+00

! initiate the matrix of right eigenvectors
!
      evroe(2,1,idn) = 1.0d+00
      evroe(2,2,ivy) = 1.0d+00
      evroe(2,3,ivz) = 1.0d+00
      evroe(2,4,idn) = 1.0d+00
      evroe(2,5,idn) = 1.0d+00

! unset the first execution flag
!
      first = .false.

    end if ! first execution

! calculate characteristic speeds and useful variables
!
    vv = sum(q(ivx:ivz)**2)
    vh = 0.5d+00 * vv
    c2 = gammam1 * (q(ien) - vh)
    na = 0.5d+00 / c2
    cc = sqrt(c2)
    vc = q(ivx) * cc
    ng = na * gammam1
    nd = 2.0d+00 * ng
    nw = na * vc
    nh = na * gammam1 * vh
    nc = na * cc

! prepare eigenvalues
!
    c(1)           = q(ivx) - cc
    c(2)           = q(ivx)
    c(3)           = q(ivx)
    c(4)           = q(ivx)
    c(5)           = q(ivx) + cc

! update the varying elements of the matrix of left eigenvectors
!
    evroe(1,idn,1) =   nh + nw
    evroe(1,ivx,1) = - ng * q(ivx) - nc
    evroe(1,ivy,1) = - ng * q(ivy)
    evroe(1,ivz,1) = - ng * q(ivz)
    evroe(1,ien,1) =   ng

    evroe(1,idn,2) = - q(ivy)

    evroe(1,idn,3) = - q(ivz)

    evroe(1,idn,4) = 1.0d+00 - ng * vv
    evroe(1,ivx,4) =   nd * q(ivx)
    evroe(1,ivy,4) =   nd * q(ivy)
    evroe(1,ivz,4) =   nd * q(ivz)
    evroe(1,ien,4) = - nd

    evroe(1,idn,5) =   nh - nw
    evroe(1,ivx,5) = - ng * q(ivx) + nc
    evroe(1,ivy,5) = - ng * q(ivy)
    evroe(1,ivz,5) = - ng * q(ivz)
    evroe(1,ien,5) =   ng

! update the varying elements of the matrix of right eigenvectors
!
    evroe(2,1,ivx) = q(ivx) - cc
    evroe(2,1,ivy) = q(ivy)
    evroe(2,1,ivz) = q(ivz)
    evroe(2,1,ien) = q(ien) - vc

    evroe(2,2,ien) = q(ivy)

    evroe(2,3,ien) = q(ivz)

    evroe(2,4,ivx) = q(ivx)
    evroe(2,4,ivy) = q(ivy)
    evroe(2,4,ivz) = q(ivz)
    evroe(2,4,ien) = vh

    evroe(2,5,ivx) = q(ivx) + cc
    evroe(2,5,ivy) = q(ivy)
    evroe(2,5,ivz) = q(ivz)
    evroe(2,5,ien) = q(ien) + vc

! copy matrices of eigenvectors
!
    l(1:nv,1:nv) = evroe(1,1:nv,1:nv)
    r(1:nv,1:nv) = evroe(2,1:nv,1:nv)

!-------------------------------------------------------------------------------
!
  end subroutine esystem_roe_hd_adi
!
!*******************************************************************************
!
! ISOTHERMAL MAGNETOHYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_MHD_ISO:
! ----------------------------
!
!   Subroutine converts primitive variables to their corresponding
!   conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons_mhd_iso(n, q, u)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: u

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
      u(ibx,i) = q(ibx,i)
      u(iby,i) = q(iby,i)
      u(ibz,i) = q(ibz,i)
      u(ibp,i) = q(ibp,i)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_mhd_iso
!
!===============================================================================
!
! subroutine CONS2PRIM_MHD_ISO:
! ----------------------------
!
!   Subroutine converts conservative variables to their corresponding
!   primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim_mhd_iso(n, u, q)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: u
    real(kind=8), dimension(nv,n), intent(out) :: q

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
      q(ibx,i) = u(ibx,i)
      q(iby,i) = u(iby,i)
      q(ibz,i) = u(ibz,i)
      q(ibp,i) = u(ibp,i)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_mhd_iso
!
!===============================================================================
!
! subroutine FLUXSPEED_MHD_ISO:
! ----------------------------
!
!   Subroutine calculates physical fluxes and characteristic speeds from a
!   given equation system.
!
!   Arguments:
!
!     n      - the length of input and output vectors;
!     q      - the input array of primitive variables;
!     u      - the input array of conservative variables;
!     f      - the output vector of fluxes;
!     cm, cp - the output vector of left- and right-going characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_mhd_iso(n, q, u, f, cm, cp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                , intent(in)  :: n
    real(kind=8), dimension(nv,n)          , intent(in)  :: q, u
    real(kind=8), dimension(nv,n)          , intent(out) :: f
    real(kind=8), dimension(n)   , optional, intent(out) :: cm, cp

! local variables
!
    integer      :: i
    real(kind=8) :: bx2, by2, bz2, bb, pr, pt
    real(kind=8) :: fa, fb, fc, cf
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux calculation
!
    call start_timer(imf)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! prepare pressures and scalar product
!
      bx2 = q(ibx,i) * q(ibx,i)
      by2 = q(iby,i) * q(iby,i)
      bz2 = q(ibz,i) * q(ibz,i)
      bb  = bx2 + by2 + bz2
      pr  = csnd2 * q(idn,i)
      pt  = pr + 0.5d+00 * bb

! calculate the magnetohydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i) - bx2
      f(imy,i) = q(ivx,i) * u(imy,i) - q(ibx,i) * q(iby,i)
      f(imz,i) = q(ivx,i) * u(imz,i) - q(ibx,i) * q(ibz,i)
      f(imx,i) = f(imx,i) + pt
      f(ibx,i) = q(ibp,i)
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
      f(ibp,i) = cmax2 * q(ibx,i)

! calculate the fast magnetosonic speed
!
      fa       = csnd2 * q(idn,i)
      fb       = fa + bb
      fc       = fb * fb - 4.0d+00 * fa * bx2
      if (fc > 0.0d+00) then
        cf     = sqrt(0.5d+00 * (fb + sqrt(fc)) / q(idn,i))
      else
        cf     = sqrt(0.5d+00 *  fb             / q(idn,i))
      end if

! calculate characteristic speeds
!
      cm(i)    = q(ivx,i) - cf
      cp(i)    = q(ivx,i) + cf

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_mhd_iso
!
!===============================================================================
!
! function MAXSPEED_MHD_ISO:
! -------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!===============================================================================
!
  function maxspeed_mhd_iso(qq) result(maxspeed)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in) :: qq

! return value
!
    real(kind=8) :: maxspeed

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: vv, bb, v, c
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the maximum speed estimation
!
    call start_timer(imm)
#endif /* PROFILE */

! reset the maximum speed
!
    maxspeed = 0.0d+00

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(qq(ivx:ivz,i,j,k) * qq(ivx:ivz,i,j,k))
          v  = sqrt(vv)
          bb = sum(qq(ibx:ibz,i,j,k) * qq(ibx:ibz,i,j,k))

! calculate the fast magnetosonic speed
!
          c = sqrt(csnd2 + bb / qq(idn,i,j,k))

! calculate the maximum of speed
!
          maxspeed = max(maxspeed, v + c)

        end do ! i = ib, ie
      end do ! j = jb, je
    end do ! k = kb, ke

#ifdef PROFILE
! stop accounting time for the maximum speed estimation
!
    call stop_timer(imm)
#endif /* PROFILE */

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_mhd_iso
!
!===============================================================================
!
! subroutine ESYSTEM_ROE_MHD_ISO:
! ------------------------------
!
!   Subroutine computes eigenvalues and eigenvectors for a given set of
!   equations and input variables.
!
!   Arguments:
!
!     x - ratio of the perpendicular magnetic field component difference
!     y - ratio of the density
!     q - the intermediate Roe state vector;
!     c - the vector of eigenvalues;
!     r - the matrix of right eigenvectors;
!     l - the matrix of left eigenvectors;
!
!   References:
!
!     [1] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!     [2] Balsara, D. S.
!         "Linearized Formulation of the Riemann Problem for Adiabatic and
!          Isothermal Magnetohydrodynamics",
!         The Astrophysical Journal Suplement Series, 1998, 116, pp. 119-131
!
!===============================================================================
!
  subroutine esystem_roe_mhd_iso(x, y, q, c, r, l)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)                  , intent(in)    :: x, y
    real(kind=8), dimension(nv)   , intent(in)    :: q
    real(kind=8), dimension(nv)   , intent(inout) :: c
    real(kind=8), dimension(nv,nv), intent(inout) :: l, r

! saved variables
!
    logical     , save :: first = .true.

! local variables
!
    real(kind=8) :: di, btsq, bt_starsq, casq, twid_csq
    real(kind=8) :: ct2, tsum, tdif, cf2_cs2, cfsq, cf, cssq, cs, ca
    real(kind=8) :: bt, bt_star, bet2, bet3, bet2_star, bet3_star, bet_starsq
    real(kind=8) :: alpha_f, alpha_s
    real(kind=8) :: sqrtd, s, twid_c, qf, qs, af_prime, as_prime
    real(kind=8) :: norm, cff, css, af, as, afpb, aspb, q2_star, q3_star, vqstr
!
!-------------------------------------------------------------------------------
!
! prepare the internal arrays at the first run
!
    if (first) then

! reset all elements
!
      evroe(:, : ,:) = 0.0d+00

! unset the first execution flag
!
      first = .false.

    end if ! first execution

! prepare coefficients for eigenvalues
!
    di        = 1.0d+00 / q(idn)
    casq      = q(ibx) * q(ibx) * di
    ca        = sqrt(casq)
    btsq      = q(iby) * q(iby) + q(ibz) * q(ibz)
    bt_starsq = btsq * y
    twid_csq  = csnd2 + x
    ct2       = bt_starsq * di
    tsum      = casq + ct2 + twid_csq
    tdif      = casq + ct2 - twid_csq
    cf2_cs2   = sqrt(tdif * tdif + 4.0d+00 * twid_csq * ct2)
    cfsq      = 0.5d+00 * (tsum + cf2_cs2)
    cf        = sqrt(cfsq)
    cssq      = twid_csq * casq / cfsq
    cs        = sqrt(cssq)

! prepare eigenvalues
!
    c(1)      = q(ivx) - cf
    c(2)      = q(ivx) - ca
    c(3)      = q(ivx) - cs
    c(4)      = q(ivx)
    c(5)      = q(ivx) + cs
    c(6)      = q(ivx) + ca
    c(7)      = q(ivx) + cf
    c(8)      = c(7)

! calculate the eigenvectors only if the waves propagate in both direction
!
    if (c(1) >= 0.0d+00) return
    if (c(7) <= 0.0d+00) return

! prepare remaining coefficients for eigenvectors
!
    bt        = sqrt(btsq)
    bt_star   = sqrt(bt_starsq)
    if (bt == 0.0d+00) then
      bet2 = 1.0d+00
      bet3 = 0.0d+00
    else
      bet2 = q(iby) / bt
      bet3 = q(ibz) / bt
    end if
    bet2_star  = bet2 / sqrt(y)
    bet3_star  = bet3 / sqrt(y)
    bet_starsq = bet2_star * bet2_star + bet3_star * bet3_star

    if ((cfsq - cssq) == 0.0d+00) then
      alpha_f = 1.0d+00
      alpha_s = 0.0d+00
    else if ((twid_csq - cssq) <= 0.0d+00) then
      alpha_f = 0.0d+00
      alpha_s = 1.0d+00
    else if ((cfsq - twid_csq) <= 0.0d+00) then
      alpha_f = 1.0d+00
      alpha_s = 0.0d+00
    else
      alpha_f = sqrt((twid_csq - cssq) / (cfsq - cssq))
      alpha_s = sqrt((cfsq - twid_csq) / (cfsq - cssq))
    end if

    sqrtd    = sqrt(q(idn))
    s        = sign(1.0d+00, q(ibx))
    twid_c   = sqrt(twid_csq)
    qf       = cf * alpha_f * s
    qs       = cs * alpha_s * s
    af_prime = twid_c * alpha_f / sqrtd
    as_prime = twid_c * alpha_s / sqrtd

! update the varying elements of the matrix of right eigenvectors
!
! left-going fast wave
!
    evroe(2,1,idn) = alpha_f
    evroe(2,1,ivx) = alpha_f * c(1)
    evroe(2,1,ivy) = alpha_f * q(ivy) + qs * bet2_star
    evroe(2,1,ivz) = alpha_f * q(ivz) + qs * bet3_star
    evroe(2,1,iby) = as_prime * bet2_star
    evroe(2,1,ibz) = as_prime * bet3_star

! left-going Alfvèn wave
!
    evroe(2,2,ivy) = - bet3
    evroe(2,2,ivz) =   bet2
    evroe(2,2,iby) = - bet3 * s / sqrtd
    evroe(2,2,ibz) =   bet2 * s / sqrtd

! left-going slow wave
!
    evroe(2,3,idn) = alpha_s
    evroe(2,3,ivx) = alpha_s * c(3)
    evroe(2,3,ivy) = alpha_s * q(ivy) - qf * bet2_star
    evroe(2,3,ivz) = alpha_s * q(ivz) - qf * bet3_star
    evroe(2,3,iby) = - af_prime * bet2_star
    evroe(2,3,ibz) = - af_prime * bet3_star

! right-going slow wave
!
    evroe(2,5,idn) = alpha_s
    evroe(2,5,ivx) = alpha_s * c(5)
    evroe(2,5,ivy) = alpha_s * q(ivy) + qf * bet2_star
    evroe(2,5,ivz) = alpha_s * q(ivz) + qf * bet3_star
    evroe(2,5,iby) = evroe(2,3,iby)
    evroe(2,5,ibz) = evroe(2,3,ibz)

! right-going Alfvèn wave
!
    evroe(2,6,ivy) =   bet3
    evroe(2,6,ivz) = - bet2
    evroe(2,6,iby) = evroe(2,2,iby)
    evroe(2,6,ibz) = evroe(2,2,ibz)

! right-going fast wave
!
    evroe(2,7,idn) = alpha_f
    evroe(2,7,ivx) = alpha_f * c(7)
    evroe(2,7,ivy) = alpha_f * q(ivy) - qs * bet2_star
    evroe(2,7,ivz) = alpha_f * q(ivz) - qs * bet3_star
    evroe(2,7,iby) = evroe(2,1,iby)
    evroe(2,7,ibz) = evroe(2,1,ibz)

! update the varying elements of the matrix of left eigenvectors
!
    norm = 0.5d+00 / twid_csq
    cff  = norm * alpha_f * cf
    css  = norm * alpha_s * cs
    qf   = qf * norm
    qs   = qs * norm
    af   = norm * af_prime * q(idn)
    as   = norm * as_prime * q(idn)
    afpb = norm * af_prime * bt_star
    aspb = norm * as_prime * bt_star

    q2_star = bet2_star / bet_starsq
    q3_star = bet3_star / bet_starsq
    vqstr   = q(ivy) * q2_star + q(ivz) * q3_star

! left-going fast wave
!
    evroe(1,idn,1) = cff * c(7) - qs * vqstr - aspb
    evroe(1,ivx,1) = - cff
    evroe(1,ivy,1) = qs * q2_star
    evroe(1,ivz,1) = qs * q3_star
    evroe(1,iby,1) = as * q2_star
    evroe(1,ibz,1) = as * q3_star

! left-going Alfvèn wave
!
    evroe(1,idn,2) =   0.5d+00 * (q(ivy) * bet3 - q(ivz) * bet2)
    evroe(1,ivy,2) = - 0.5d+00 * bet3
    evroe(1,ivz,2) =   0.5d+00 * bet2
    evroe(1,iby,2) = - 0.5d+00 * sqrtd * bet3 * s
    evroe(1,ibz,2) =   0.5d+00 * sqrtd * bet2 * s

! left-going slow wave
!
    evroe(1,idn,3) =   css * c(5) + qf * vqstr + afpb
    evroe(1,ivx,3) = - css
    evroe(1,ivy,3) = - qf * q2_star
    evroe(1,ivz,3) = - qf * q3_star
    evroe(1,iby,3) = - af * q2_star
    evroe(1,ibz,3) = - af * q3_star

! right-going slow wave
!
    evroe(1,idn,5) = - css * c(3) - qf * vqstr + afpb
    evroe(1,ivx,5) =   css
    evroe(1,ivy,5) = - evroe(1,ivy,3)
    evroe(1,ivz,5) = - evroe(1,ivz,3)
    evroe(1,iby,5) =   evroe(1,iby,3)
    evroe(1,ibz,5) =   evroe(1,ibz,3)

! right-going Alfvèn wave
!
    evroe(1,idn,6) = - evroe(1,idn,2)
    evroe(1,ivy,6) = - evroe(1,ivy,2)
    evroe(1,ivz,6) = - evroe(1,ivz,2)
    evroe(1,iby,6) =   evroe(1,iby,2)
    evroe(1,ibz,6) =   evroe(1,ibz,2)

! right-going fast wave
!
    evroe(1,idn,7) = - cff * c(1) + qs * vqstr - aspb
    evroe(1,ivx,7) =   cff
    evroe(1,ivy,7) = - evroe(1,ivy,1)
    evroe(1,ivz,7) = - evroe(1,ivz,1)
    evroe(1,iby,7) =   evroe(1,iby,1)
    evroe(1,ibz,7) =   evroe(1,ibz,1)

! copy matrices of eigenvectors
!
    l(1:nv,1:nv) = evroe(1,1:nv,1:nv)
    r(1:nv,1:nv) = evroe(2,1:nv,1:nv)

!-------------------------------------------------------------------------------
!
  end subroutine esystem_roe_mhd_iso
!
!*******************************************************************************
!
! ADIABATIC MAGNETOHYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_MHD_ADI:
! ----------------------------
!
!   Subroutine converts primitive variables to their corresponding
!   conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons_mhd_adi(n, q, u)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: u

! local variables
!
    integer      :: i
    real(kind=8) :: ei, ek, em
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
      u(ibx,i) = q(ibx,i)
      u(iby,i) = q(iby,i)
      u(ibz,i) = q(ibz,i)
      u(ibp,i) = q(ibp,i)
      ei       = gammam1i * q(ipr,i)
      ek       = 0.5d+00 * (u(imx,i) * q(ivx,i) + u(imy,i) * q(ivy,i)          &
                                                + u(imz,i) * q(ivz,i))
      em       = 0.5d+00 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      u(ien,i) = ei + ek + em

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_mhd_adi
!
!===============================================================================
!
! subroutine CONS2PRIM_MHD_ADI:
! ----------------------------
!
!   Subroutine converts conservative variables to their corresponding
!   primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim_mhd_adi(n, u, q)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: u
    real(kind=8), dimension(nv,n), intent(out) :: q

! local variables
!
    integer      :: i
    real(kind=8) :: ei, ek, em
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
      q(ibx,i) = u(ibx,i)
      q(iby,i) = u(iby,i)
      q(ibz,i) = u(ibz,i)
      q(ibp,i) = u(ibp,i)
      ek       = 0.5d+00 * (u(imx,i) * q(ivx,i) + u(imy,i) * q(ivy,i)          &
                                                + u(imz,i) * q(ivz,i))
      em       = 0.5d+00 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      ei       = u(ien,i) - (ek + em)
      q(ipr,i) = gammam1 * ei

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_mhd_adi
!
!===============================================================================
!
! subroutine FLUXSPEED_MHD_ADI:
! ----------------------------
!
!   Subroutine calculates physical fluxes and characteristic speeds from a
!   given equation system.
!
!   Arguments:
!
!     n      - the length of input and output vectors;
!     q      - the input array of primitive variables;
!     u      - the input array of conservative variables;
!     f      - the output vector of fluxes;
!     cm, cp - the output vector of left- and right-going characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_mhd_adi(n, q, u, f, cm, cp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                , intent(in)  :: n
    real(kind=8), dimension(nv,n)          , intent(in)  :: q, u
    real(kind=8), dimension(nv,n)          , intent(out) :: f
    real(kind=8), dimension(n)   , optional, intent(out) :: cm, cp

! local variables
!
    integer      :: i
    real(kind=8) :: bx2, by2, bz2, bb, pr, pt
    real(kind=8) :: vb
    real(kind=8) :: fa, fb, fc, cf
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux calculation
!
    call start_timer(imf)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! prepare pressures and scalar product
!
      bx2 = q(ibx,i) * q(ibx,i)
      by2 = q(iby,i) * q(iby,i)
      bz2 = q(ibz,i) * q(ibz,i)
      bb  = bx2 + by2 + bz2
      vb  = sum(q(ivx:ivz,i) * q(ibx:ibz,i))
      pr  = q(ipr,i)
      pt  = pr + 0.5d+00 * bb

! calculate the magnetohydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i) - bx2
      f(imy,i) = q(ivx,i) * u(imy,i) - q(ibx,i) * q(iby,i)
      f(imz,i) = q(ivx,i) * u(imz,i) - q(ibx,i) * q(ibz,i)
      f(imx,i) = f(imx,i) + pt
      f(ibx,i) = q(ibp,i)
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
      f(ibp,i) = cmax2 * q(ibx,i)
      f(ien,i) = q(ivx,i) * (u(ien,i) + pt) - q(ibx,i) * vb

! calculate the fast magnetosonic speed
!
      fa       = gamma * q(ipr,i)
      fb       = fa + bb
      fc       = fb * fb - 4.0d+00 * fa * bx2
      if (fc > 0.0d+00) then
        cf     = sqrt(0.5d+00 * (fb + sqrt(fc)) / q(idn,i))
      else
        cf     = sqrt(0.5d+00 *  fb             / q(idn,i))
      end if

! calculate characteristic speeds
!
      cm(i)    = q(ivx,i) - cf
      cp(i)    = q(ivx,i) + cf

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_mhd_adi
!
!===============================================================================
!
! function MAXSPEED_MHD_ADI:
! -------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!===============================================================================
!
  function maxspeed_mhd_adi(qq) result(maxspeed)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in) :: qq

! return value
!
    real(kind=8) :: maxspeed

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: vv, bb, v, c
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the maximum speed estimation
!
    call start_timer(imm)
#endif /* PROFILE */

! reset the maximum speed
!
    maxspeed = 0.0d+00

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(qq(ivx:ivz,i,j,k) * qq(ivx:ivz,i,j,k))
          v  = sqrt(vv)
          bb = sum(qq(ibx:ibz,i,j,k) * qq(ibx:ibz,i,j,k))

! calculate the fast magnetosonic speed
!
          c = sqrt((gamma * qq(ipr,i,j,k) + bb) / qq(idn,i,j,k))

! calculate the maximum of speed
!
          maxspeed = max(maxspeed, v + c)

        end do ! i = ib, ie
      end do ! j = jb, je
    end do ! k = kb, ke

#ifdef PROFILE
! stop accounting time for the maximum speed estimation
!
    call stop_timer(imm)
#endif /* PROFILE */

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_mhd_adi
!
!===============================================================================
!
! subroutine ESYSTEM_ROE_MHD_ADI:
! ------------------------------
!
!   Subroutine computes eigenvalues and eigenvectors for a given set of
!   equations and input variables.
!
!   Arguments:
!
!     x - ratio of the perpendicular magnetic field component difference
!     y - ratio of the density
!     q - the intermediate Roe state vector;
!     c - the vector of eigenvalues;
!     r - the matrix of right eigenvectors;
!     l - the matrix of left eigenvectors;
!
!   References:
!
!     [1] Stone, J. M. & Gardiner, T. A.,
!         "ATHENA: A New Code for Astrophysical MHD",
!         The Astrophysical Journal Suplement Series, 2008, 178, pp. 137-177
!     [2] Balsara, D. S.
!         "Linearized Formulation of the Riemann Problem for Adiabatic and
!          Isothermal Magnetohydrodynamics",
!         The Astrophysical Journal Suplement Series, 1998, 116, pp. 119-131
!
!===============================================================================
!
  subroutine esystem_roe_mhd_adi(x, y, q, c, r, l)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)                  , intent(in)    :: x, y
    real(kind=8), dimension(nv)   , intent(in)    :: q
    real(kind=8), dimension(nv)   , intent(inout) :: c
    real(kind=8), dimension(nv,nv), intent(inout) :: l, r

! saved variables
!
    logical     , save :: first = .true.
    real(kind=8), save :: gammam2

! local variables
!
    real(kind=8) :: di, vsq, btsq, bt_starsq, casq, hp, twid_asq
    real(kind=8) :: ct2, tsum, tdif, cf2_cs2, cfsq, cf, cssq, cs, ca, bt
    real(kind=8) :: bt_star, bet2, bet3, bet2_star, bet3_star, bet_starsq, vbet
    real(kind=8) :: alpha_f, alpha_s, af_prime, as_prime
    real(kind=8) :: sqrtd, isqrtd, s, twid_a, qf, qs, afpbb, aspbb
    real(kind=8) :: qa, qb, qc, qd
    real(kind=8) :: norm, cff, css, af, as, afpb, aspb
    real(kind=8) :: q2_star, q3_star, vqstr

! local parameters
!
    real(kind=8), parameter :: eps = epsilon(di)
!
!-------------------------------------------------------------------------------
!
! prepare the internal arrays at the first run
!
    if (first) then

! prepare coefficients
!
      gammam2  = gamma - 2.0d+00

! reset all elements
!
      evroe(:, : ,:) = 0.0d+00

! initiate the matrix of right eigenvectors
!
      evroe(2,4,idn) = 1.0d+00

! unset the first execution flag
!
      first = .false.

    end if ! first execution

! prepare coefficients for eigenvalues
!
    di        = 1.0d+00 / q(idn)
    casq      = q(ibx) * q(ibx) * di
    ca        = sqrt(casq)
    vsq       = q(ivx) * q(ivx) + q(ivy) * q(ivy) + q(ivz) * q(ivz)
    btsq      = q(iby) * q(iby) + q(ibz) * q(ibz)
    bt_starsq = (gammam1 - gammam2 * y) * btsq
    hp        = q(ien) - (casq + btsq * di)
    twid_asq  = max(eps, (gammam1 * (hp - 0.5d+00 * vsq) - gammam2 * x))
    ct2       = bt_starsq * di
    tsum      = casq + ct2 + twid_asq
    tdif      = casq + ct2 - twid_asq
    cf2_cs2   = sqrt((tdif * tdif + 4.0d+00 * twid_asq * ct2))
    cfsq      = 0.5d+00 * (tsum + cf2_cs2)
    cf        = sqrt(cfsq)
    cssq      = twid_asq * casq / cfsq
    cs        = sqrt(cssq)

! prepare eigenvalues
!
    c(1)      = q(ivx) - cf
    c(2)      = q(ivx) - ca
    c(3)      = q(ivx) - cs
    c(4)      = q(ivx)
    c(5)      = q(ivx)
    c(6)      = q(ivx) + cs
    c(7)      = q(ivx) + ca
    c(8)      = q(ivx) + cf
    c(9)      = c(8)

! calculate the eigenvectors only if the waves propagate in both direction
!
    if (c(1) >= 0.0d+00) return
    if (c(8) <= 0.0d+00) return

! prepare remaining coefficients for eigenvectors
!
    bt        = sqrt(btsq)
    bt_star   = sqrt(bt_starsq)
    if (bt == 0.0d+00) then
      bet2 = 1.0d+00
      bet3 = 0.0d+00
    else
      bet2 = q(iby) / bt
      bet3 = q(ibz) / bt
    end if
    bet2_star  = bet2 / sqrt(gammam1 - gammam2 * y)
    bet3_star  = bet3 / sqrt(gammam1 - gammam2 * y)
    bet_starsq = bet2_star * bet2_star + bet3_star * bet3_star
    vbet       = q(ivy) * bet2_star + q(ivz) * bet3_star

    if ( (cfsq - cssq) == 0.0d+00 ) then
      alpha_f = 1.0d+00
      alpha_s = 0.0d+00
    else if ( (twid_asq - cssq) <= 0.0d+00 ) then
      alpha_f = 0.0d+00
      alpha_s = 1.0d+00
    else if ( (cfsq - twid_asq) <= 0.0d+00 ) then
      alpha_f = 1.0d+00
      alpha_s = 0.0d+00
    else
      alpha_f = sqrt((twid_asq - cssq) / (cfsq - cssq))
      alpha_s = sqrt((cfsq - twid_asq) / (cfsq - cssq))
    end if

    sqrtd    = sqrt(q(idn))
    isqrtd   = 1.0d+00 / sqrtd
    s        = sign(1.0d+00, q(ibx))
    twid_a   = sqrt(twid_asq)
    qf       = cf * alpha_f * s
    qs       = cs * alpha_s * s
    af_prime = twid_a * alpha_f * isqrtd
    as_prime = twid_a * alpha_s * isqrtd
    afpbb    = af_prime * bt_star * bet_starsq
    aspbb    = as_prime * bt_star * bet_starsq

! update the varying elements of the matrix of right eigenvectors
!
    evroe(2,1,idn) = alpha_f
    evroe(2,3,idn) = alpha_s
    evroe(2,6,idn) = alpha_s
    evroe(2,8,idn) = alpha_f

    evroe(2,1,ivx) = alpha_f * c(1)
    evroe(2,3,ivx) = alpha_s * c(3)
    evroe(2,4,ivx) = q(ivx)
    evroe(2,6,ivx) = alpha_s * c(6)
    evroe(2,8,ivx) = alpha_f * c(8)

    qa = alpha_f * q(ivy)
    qb = alpha_s * q(ivy)
    qc = qs * bet2_star
    qd = qf * bet2_star

    evroe(2,1,ivy) = qa + qc
    evroe(2,2,ivy) = - bet3
    evroe(2,3,ivy) = qb - qd
    evroe(2,4,ivy) = q(ivy)
    evroe(2,6,ivy) = qb + qd
    evroe(2,7,ivy) = bet3
    evroe(2,8,ivy) = qa - qc

    qa = alpha_f * q(ivz)
    qb = alpha_s * q(ivz)
    qc = qs * bet3_star
    qd = qf * bet3_star

    evroe(2,1,ivz) = qa + qc
    evroe(2,2,ivz) = bet2
    evroe(2,3,ivz) = qb - qd
    evroe(2,4,ivz) = q(ivz)
    evroe(2,6,ivz) = qb + qd
    evroe(2,7,ivz) = - bet2
    evroe(2,8,ivz) = qa - qc

    evroe(2,1,ipr) = alpha_f * (hp - q(ivx) * cf) + qs * vbet + aspbb
    evroe(2,2,ipr) = -(q(ivy) * bet3 - q(ivz) * bet2)
    evroe(2,3,ipr) = alpha_s * (hp - q(ivx) * cs) - qf * vbet - afpbb
    evroe(2,4,ipr) = 0.5d+00 * vsq + gammam2 * x / gammam1
    evroe(2,6,ipr) = alpha_s * (hp + q(ivx) * cs) + qf * vbet - afpbb
    evroe(2,7,ipr) = - evroe(2,2,ipr)
    evroe(2,8,ipr) = alpha_f * (hp + q(ivx) * cf) - qs * vbet + aspbb

    evroe(2,1,iby) = as_prime * bet2_star
    evroe(2,2,iby) = - bet3 * s * isqrtd
    evroe(2,3,iby) = - af_prime * bet2_star
    evroe(2,6,iby) = evroe(2,3,iby)
    evroe(2,7,iby) = evroe(2,2,iby)
    evroe(2,8,iby) = evroe(2,1,iby)

    evroe(2,1,ibz) = as_prime * bet3_star
    evroe(2,2,ibz) = bet2 * s * isqrtd
    evroe(2,3,ibz) = - af_prime * bet3_star
    evroe(2,6,ibz) = evroe(2,3,ibz)
    evroe(2,7,ibz) = evroe(2,2,ibz)
    evroe(2,8,ibz) = evroe(2,1,ibz)

! update the varying elements of the matrix of left eigenvectors
!
    norm = 0.5d+00 / twid_asq
    cff  = norm * alpha_f * cf
    css  = norm * alpha_s * cs
    qf   = qf * norm
    qs   = qs * norm
    af   = norm * af_prime * q(idn)
    as   = norm * as_prime * q(idn)
    afpb = norm * af_prime * bt_star
    aspb = norm * as_prime * bt_star

    norm    = norm * gammam1
    alpha_f = alpha_f * norm
    alpha_s = alpha_s * norm
    q2_star = bet2_star / bet_starsq
    q3_star = bet3_star / bet_starsq
    vqstr   = (q(ivy) * q2_star + q(ivz) * q3_star)
    norm    = 2.0d+00 * norm

! left-going fast wave
!
    evroe(1,idn,1) = alpha_f * (vsq - hp) + cff * (cf + q(ivx))                &
                                                           - qs * vqstr - aspb
    evroe(1,ivx,1) = - alpha_f * q(ivx) - cff
    evroe(1,ivy,1) = - alpha_f * q(ivy) + qs * q2_star
    evroe(1,ivz,1) = - alpha_f * q(ivz) + qs * q3_star
    evroe(1,ipr,1) = alpha_f
    evroe(1,iby,1) = as * q2_star - alpha_f * q(iby)
    evroe(1,ibz,1) = as * q3_star - alpha_f * q(ibz)

! left-going Alfvèn wave
!
    evroe(1,idn,2) = 0.5d+00 * (q(ivy) * bet3 - q(ivz) * bet2)
    evroe(1,ivy,2) = - 0.5d+00 * bet3
    evroe(1,ivz,2) = 0.5d+00 * bet2
    evroe(1,iby,2) = - 0.5d+00 * sqrtd * bet3 * s
    evroe(1,ibz,2) = 0.5d+00 * sqrtd * bet2 * s

! left-going slow wave
!
    evroe(1,idn,3) = alpha_s * (vsq - hp) + css * (cs + q(ivx))                &
                                                           + qf * vqstr + afpb
    evroe(1,ivx,3) = - alpha_s * q(ivx) - css
    evroe(1,ivy,3) = - alpha_s * q(ivy) - qf * q2_star
    evroe(1,ivz,3) = - alpha_s * q(ivz) - qf * q3_star
    evroe(1,ipr,3) = alpha_s
    evroe(1,iby,3) = - af * q2_star - alpha_s * q(iby)
    evroe(1,ibz,3) = - af * q3_star - alpha_s * q(ibz)

! entropy wave
!
    evroe(1,idn,4) = 1.0d+00 - norm * (0.5d+00 * vsq - gammam2 * x / gammam1)
    evroe(1,ivx,4) = norm * q(ivx)
    evroe(1,ivy,4) = norm * q(ivy)
    evroe(1,ivz,4) = norm * q(ivz)
    evroe(1,ipr,4) = - norm
    evroe(1,iby,4) = norm * q(iby)
    evroe(1,ibz,4) = norm * q(ibz)

! right-going slow wave
!
    evroe(1,idn,6) = alpha_s * (vsq - hp) + css * (cs - q(ivx))                &
                                                           - qf * vqstr + afpb
    evroe(1,ivx,6) = - alpha_s * q(ivx) + css
    evroe(1,ivy,6) = - alpha_s * q(ivy) + qf * q2_star
    evroe(1,ivz,6) = - alpha_s * q(ivz) + qf * q3_star
    evroe(1,ipr,6) = alpha_s
    evroe(1,iby,6) = evroe(1,iby,3)
    evroe(1,ibz,6) = evroe(1,ibz,3)

! right-going Alfvèn wave
!
    evroe(1,idn,7) = - evroe(1,idn,2)
    evroe(1,ivy,7) = - evroe(1,ivy,2)
    evroe(1,ivz,7) = - evroe(1,ivz,2)
    evroe(1,iby,7) =   evroe(1,iby,2)
    evroe(1,ibz,7) =   evroe(1,ibz,2)

! right-going fast wave
!
    evroe(1,idn,8) = alpha_f * (vsq - hp) + cff * (cf - q(ivx))                &
                                                           + qs * vqstr - aspb
    evroe(1,ivx,8) = - alpha_f * q(ivx) + cff
    evroe(1,ivy,8) = - alpha_f * q(ivy) - qs * q2_star
    evroe(1,ivz,8) = - alpha_f * q(ivz) - qs * q3_star
    evroe(1,ipr,8) = alpha_f
    evroe(1,iby,8) = evroe(1,iby,1)
    evroe(1,ibz,8) = evroe(1,ibz,1)

! copy matrices of eigenvectors
!
    l(1:nv,1:nv) = evroe(1,1:nv,1:nv)
    r(1:nv,1:nv) = evroe(2,1:nv,1:nv)

!-------------------------------------------------------------------------------
!
  end subroutine esystem_roe_mhd_adi
!
!*******************************************************************************
!
! ADIABATIC SPECIAL RELATIVITY HYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_SRHD_ADI:
! -----------------------------
!
!   Subroutine converts primitive variables to their corresponding
!   conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons_srhd_adi(n, q, u)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q
    real(kind=8), dimension(nv,n), intent(out) :: u

! local variables
!
    integer      :: i
    real(kind=8) :: vv, vm, vs, ww, ek, ei
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! calculate the square of velocity, the Lorentz factor and specific enthalpy
!
      vv  = sum(q(ivx:ivz,i) * q(ivx:ivz,i))
      vm  = 1.0d+00 - vv
      vs  = sqrt(vm)
      ww  = (q(idn,i) + q(ipr,i) / gammaxi) / vm

! calculate conservative variables
!
      u(idn,i) =      q(idn,i) / vs
      u(imx,i) = ww * q(ivx,i)
      u(imy,i) = ww * q(ivy,i)
      u(imz,i) = ww * q(ivz,i)
      u(ien,i) = ww - q(ipr,i) - u(idn,i)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_srhd_adi
!
!===============================================================================
!
! subroutine CONS2PRIM_SRHD_ADI:
! -----------------------------
!
!   Subroutine converts conservative variables to their corresponding
!   primitive representation using an interative method.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim_srhd_adi(n, u, q)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: u
    real(kind=8), dimension(nv,n), intent(out) :: q

! local variables
!
    integer      :: i
    real(kind=8) :: mm, en, dn
    real(kind=8) :: em, w, wm
    real(kind=8) :: vv, vm, vs
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable conversion
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! prepare variables which do not change during the Newton-Ralphson iterations
!
      mm  = sum(u(imx:imz,i) * u(imx:imz,i))
      en  = u(ien,i) + u(idn,i)
      dn  = u(idn,i)

! get the lower bound for W updated by the minimum pressure with assumed Γ = 1
!
      wm = sqrt(dn * dn + mm) + pmin / gammaxi

! calculate V² corresponding to the minimum W
!
      call w_to_vv_srhd_adi(mm, wm, vv)

! estimate the total energy corresponding to the minimum of W
!
      em = wm - dn - pmin

! set the initial W to the minimum value
!
      w  = wm

! find the exact W using an Newton-Ralphson interative method
!
      call nr_iterate(mm, en, dn, wm, w, vv)

! prepare coefficients
!
      vm = 1.0d+00 - vv
      vs = sqrt(vm)

! calculate the primitive variables
!
      q(idn,i) = dn * vs
      q(ivx,i) = u(imx,i) / w
      q(ivy,i) = u(imy,i) / w
      q(ivz,i) = u(imz,i) / w
      q(ipr,i) = gammaxi * (w - dn / vs) * vm

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_srhd_adi
!
!===============================================================================
!
! subroutine FLUXSPEED_SRHD_ADI:
! -----------------------------
!
!   Subroutine calculates physical fluxes and characteristic speeds from a
!   given equation system.
!
!   Arguments:
!
!     n      - the length of input and output vectors;
!     q      - the input array of primitive variables;
!     u      - the input array of conservative variables;
!     f      - the output vector of fluxes;
!     cm, cp - the output vector of left- and right-going characteristic speeds;
!
!   References:
!
!     [1] Mignone, A., Bodo, G.,
!         "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics",
!         Monthly Notices of the Royal Astronomical Society, 2005, 364, 126-136
!
!===============================================================================
!
  subroutine fluxspeed_srhd_adi(n, q, u, f, cm, cp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                , intent(in)  :: n
    real(kind=8), dimension(nv,n)          , intent(in)  :: q, u
    real(kind=8), dimension(nv,n)          , intent(out) :: f
    real(kind=8), dimension(n)   , optional, intent(out) :: cm, cp

! local variables
!
    integer      :: i
    real(kind=8) :: vv, ww, c2, ss, cc, fc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for flux calculation
!
    call start_timer(imf)
#endif /* PROFILE */

! iterate over all positions
!
    do i = 1, n

! calculate the relativistic hydrodynamic fluxes (eq. 2 in [1])
!
      f(idn,i) = u(idn,i) * q(ivx,i)
      f(imx,i) = u(imx,i) * q(ivx,i) + q(ipr,i)
      f(imy,i) = u(imy,i) * q(ivx,i)
      f(imz,i) = u(imz,i) * q(ivx,i)
      f(ien,i) = u(imx,i) - f(idn,i)

! calculate the relativistic speed of sound (eqs. 4, 22 and 23 in [1])
!
      ww       = q(idn,i) + q(ipr,i) / gammaxi
      c2       = gamma * q(ipr,i) / ww
      vv       = sum(q(ivx:ivz,i) * q(ivx:ivz,i))
      ss       = c2 * (1.0d+00 - vv) / (1.0d+00 - c2)
      fc       = 1.0d+00 + ss
      cc       = sqrt(ss * (fc - q(ivx,i)**2))

! calculate characteristic speeds (eq. 23 in [1])
!
      cm(i)    = (q(ivx,i) - cc) / fc
      cp(i)    = (q(ivx,i) + cc) / fc

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_srhd_adi
!
!===============================================================================
!
! function MAXSPEED_SRHD_ADI:
! --------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!===============================================================================
!
  function maxspeed_srhd_adi(qq) result(maxspeed)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in) :: qq

! return value
!
    real(kind=8) :: maxspeed

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: vv, v, cc, ww, c2, ss, fc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the maximum speed estimation
!
    call start_timer(imm)
#endif /* PROFILE */

! reset the maximum speed
!
    maxspeed = 0.0d+00

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(qq(ivx:ivz,i,j,k) * qq(ivx:ivz,i,j,k))
          v  = sqrt(vv)

! calculate the square of the sound speed
!
          ww = qq(idn,i,j,k) + qq(ipr,i,j,k) / gammaxi
          c2 = gamma * qq(ipr,i,j,k) / ww
          ss = c2 * (1.0d+00 - vv) / (1.0d+00 - c2)
          fc = 1.0d+00 + ss
          cc = sqrt(ss * (fc - vv))

! calculate the maximum of speed
!
          maxspeed = max(maxspeed, (v + cc) / fc)

        end do ! i = ib, ie
      end do ! j = jb, je
    end do ! k = kb, ke

#ifdef PROFILE
! stop accounting time for the maximum speed estimation
!
    call stop_timer(imm)
#endif /* PROFILE */

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_srhd_adi
!
!===============================================================================
!
! subroutine NR_FUNCTION_SRHD_ADI:
! -------------------------------
!
!   Subroutine calculates the value of function
!
!     F(W)     = W - P - E
!
!   and its derivative
!
!     dF(W)/dW = 1 - dP/dW
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E;
!     w      - input coefficients W;
!     f, df  - output values of the function F(W) and its derivative,
!              respectively;
!
!   References:
!
!     Noble, S. C., Gammie, C. F., McKinney, J. C, Del Zanna, L.,
!     "Primitive Variable Solvers for Conservative General Relativistic
!      Magnetohydrodynamics",
!     The Astrophysical Journal, 2006, vol. 641, pp. 626-637
!
!===============================================================================
!
  subroutine nr_function_srhd_adi(mm, en, dn, w, f, df)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, en, dn, w
    real(kind=8), intent(out)   :: f, df

! local variables
!
    real(kind=8) :: w2, tm, vv, vm, gm, pr, dv, dg, dp
!
!-------------------------------------------------------------------------------
!
! prepare W multiplications
!
    w2 = w * w

! calculate the velocity and its derivative
!
    call w_to_vv_srhd_adi(mm, w, vv, dv)

! calculate 1 - |V|²
!
    vm = 1.0d+00 - vv

! calculate Lorentz factor and its derivative
!
!  Γ(|V|²) = 1 / sqrt(1 - |V|²)
! dΓ/dW    = ½ Γ³ d|V|²/dW
!
    gm = 1.0d+00 / sqrt(vm)
    dg = 0.5d+00 * gm**3 * dv

! calculate the thermal pressure and its derivative
!
!  P(W) = (γ - 1)/γ (W - DΓ) (1 - |V|²)
! dP/dW = (γ - 1)/γ [(1 - D dΓ/dW) (1 - |V|²) - (W - DΓ) d|V|²/dW]
!
    tm = w - dn * gm
    pr = gammaxi * tm * vm
    dp = gammaxi * ((1.0d+00 - dn * dg) * vm - tm * dv)

! calculate F(W)
!
    f  = w - pr - en

! calculate dF(W)/dW
!
! dF(W)/dW = 1 - dP/dW
!
    df = 1.0d+00 - dp

!-------------------------------------------------------------------------------
!
  end subroutine nr_function_srhd_adi
!
!===============================================================================
!
! subroutine NR_ITERATE_SRHD_ADI_1DW:
! ----------------------------------
!
!   Subroutine finds a root W of equation
!
!     F(W) = W - P - E = 0
!
!   using the Newton-Raphson 1Dw iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     w, vv  - input/output coefficients W and |V|²;
!
!   References:
!
!     Noble, S. C., Gammie, C. F., McKinney, J. C, Del Zanna, L.,
!     "Primitive Variable Solvers for Conservative General Relativistic
!      Magnetohydrodynamics",
!     The Astrophysical Journal, 2006, vol. 641, pp. 626-637
!
!===============================================================================
!
  subroutine nr_iterate_srhd_adi_1dw(mm, en, dn, wm, w, vv)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, en, dn, wm
    real(kind=8), intent(inout) :: w, vv

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: f, df, dw
    real(kind=8) :: er
!
!-------------------------------------------------------------------------------
!
! initialize iteration parameters
!
    keep = .true.
    it   = 100
    cn   = 2

! iterate using the Newton-Raphson method in order to find a root w of the
! function
!
! F(W) = W - P - E = 0
!
    do while(keep)

! calculate F(W) and dF(W)/dW
!
      call nr_function_srhd_adi(mm, en, dn, w, f, df)

! calculate the increment dW
!
      dw = - f / df

! correct W
!
      w  = max(wm, w + dw)

! calculate the normalized error
!
      er = abs(dw / w)

! check the convergence
!
      if (er < tol) then
        if (cn <= 0) keep = .false.
        cn  = cn - 1
      end if

! break if the number of iterations exceeded the maximum value
!
      if (it <= 0) keep = .false.

! decrease the number of remaining iterations
!
      it = it - 1

    end do

! calculate |V|² from W
!
    call w_to_vv_srhd_adi(mm, w, vv)

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srhd_adi_1dw
!
!===============================================================================
!
! subroutine NR_ITERATE_SRHD_ADI_2D:
! ---------------------------------
!
!   Subroutine finds a root W of equation
!
!     F(W) = W - P - E = 0
!
!   using the Newton-Raphson 2D iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     w, vv  - input/output coefficients W and |V|²;
!
!   References:
!
!     Noble, S. C., Gammie, C. F., McKinney, J. C, Del Zanna, L.,
!     "Primitive Variable Solvers for Conservative General Relativistic
!      Magnetohydrodynamics",
!     The Astrophysical Journal, 2006, vol. 641, pp. 626-637
!
!===============================================================================
!
  subroutine nr_iterate_srhd_adi_2d(mm, en, dn, wm, w, vv)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, en, dn, wm
    real(kind=8), intent(inout) :: w, vv

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: ww, vm, gm
    real(kind=8) :: pr, dpw, dpv
    real(kind=8) :: f, dfw, dfv, df
    real(kind=8) :: g, dgw, dgv, dg
    real(kind=8) :: det, jfw, jfv, jgw, jgv
    real(kind=8) :: dw, dv
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
! initialize iteration parameters
!
    keep = .true.
    it   = nmax
    cn   = next

! calculate the initial |V|² from the guess of W
!
!  |V|²(W) = |M|² / W²
!
    vv = mm / (w * w)

! iterate using the Newton-Raphson method in order to find the roots W and |V|²
! of functions
!
! F(W,|V|²) = W - P - E = 0
! G(W,|V|²) = |V|² W² - |M|² = 0
!
    do while(keep)

! calculate W², (1 - |V|²), and the Lorentz factor
!
      ww  = w * w
      vm  = 1.0d+00 - vv
      gm  = 1.0d+00 / sqrt(vm)

! calculate the thermal pressure and its derivatives
!
!  P(W,|V|²) = (γ - 1)/γ (W - D Γ) (1 - |V|²)
! dP/dW      = (γ - 1)/γ (1 - |V|²)
! dP/d|V|²   = (γ - 1)/γ (- W - 1/2 D Γ)
!
      pr  = gammaxi * (w - dn * gm) * vm
      dpw = gammaxi * vm
      dpv = gammaxi * (- w + 0.5d+00 * dn * gm)

! calculate F(W,|V|²) and G(W,|V|²)
!
      f   = w - pr - en
      g   = vv * ww  - mm

! calculate dF(W,|V|²)/dW and dF(W,|V|²)/d|V|²
!
      dfw = 1.0d+00 - dpw
      dfv =         - dpv

! calculate dG(W,|V|²)/dW and dG(W,|V|²)/d|V|²
!
      dgw = 2.0d+00 * vv * w
      dgv = ww

! invert the Jacobian J = | dF/dW, dF/d|V|² |
!                         | dG/dW, dG/d|V|² |
!
      det = dfw * dgv - dfv * dgw

      jfw =   dgv / det
      jgw = - dfv / det
      jfv = - dgw / det
      jgv =   dfw / det

! calculate increments dW and d|V|²
!
      dw  = f * jfw + g * jgw
      dv  = f * jfv + g * jgv

! correct W and |V|²
!
      w   = w  - dw
      vv  = vv - dv

! calculate the normalized error
!
      err = max(abs(dw / w), abs(dv))

! check the convergence
!
      if (err < tol) then
        if (cn <= 0) keep = .false.
        cn  = cn - 1
      end if

! break if the number of iterations exceeded the maximum value
!
      if (it <= 0) keep = .false.

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! continue interations

! print information about failed convergence
!
    if (err >= tol) then
      print *, '[SRHD, 2D ] Convergence not reached: ', err
    end if
    if (vv  >= 1.0d+00) then
      print *, '[SRHD, 2D ] Unphysical speed: ', vv
    end if
    if (w   <= 0.0d+00) then
      print *, '[SRHD, 2D ] Unphysical enthalpy: ', w
    end if

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srhd_adi_2d
!
!===============================================================================
!
! subroutine W_TO_VV_SRHD_ADI:
! ---------------------------
!
!   Subroutine calculates the squared velocity and its W derivative from W
!   and other parameters.
!
!     |V|²(W) =     [|M|² W² + S² (2 W + |B|²)]            / [W (W + |B|²)]²
!    d|V|²/dW = - 2 {|M|² W³ + S² [3 W (W + |B|²) + |B|⁴]} / [W (W + |B|²)]³
!
!   Arguments:
!
!     mm - input coefficients for |M|²;
!     w  - input coefficient W;
!     vv - output value of |V|²;
!
!   References:
!
!     Noble, S. C., Gammie, C. F., McKinney, J. C, Del Zanna, L.,
!     "Primitive Variable Solvers for Conservative General Relativistic
!      Magnetohydrodynamics",
!     The Astrophysical Journal, 2006, vol. 641, pp. 626-637
!
!===============================================================================
!
  subroutine w_to_vv_srhd_adi(mm, w, vv, dv)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8)          , intent(in)  :: mm, w
    real(kind=8)          , intent(out) :: vv
    real(kind=8), optional, intent(out) :: dv
!
!-------------------------------------------------------------------------------
!
! calculate the squared velocity and its derivative
!
    vv = mm / w**2
    if (present(dv)) then
      dv = - 2.0d+00 * vv / w
    end if

!-------------------------------------------------------------------------------
!
  end subroutine w_to_vv_srhd_adi

!===============================================================================
!
end module equations
