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
  integer            , save :: imi, imc, imf, imm, imp, imb
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

! variable boundary values
!
  real(kind=8), dimension(:,:,:), allocatable, save :: qpbnd

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

! the upper limits for the Lorentz factor and corresponding |v|²
!
  real(kind=8)     , save :: lmax    = 1.0d+06
  real(kind=8)     , save :: vmax    = 0.999999999999d+00

! the upper bound for the sonic Mach number
!
  real(kind=8)     , save :: msmax   = 1.0d+03
  real(kind=8)     , save :: msfac   = 3.0d-06 / 5.0d+00

! the tolerance for Newton-Raphson interative method, the maximum number of
! iterations and the number of extra iterations for polishing
!
  real(kind=8)     , save :: tol     = 1.0d-10
  integer          , save :: nrmax   = 100
  integer          , save :: nrext   = 2

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
  public :: qpbnd

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
    integer                :: p
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
    call set_timer('equations:: variable solver'    , imp)
    call set_timer('equations:: flux calculation'   , imf)
    call set_timer('equations:: speed estimation'   , imm)
    call set_timer('equations:: initial brackets'   , imb)

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
    call get_parameter_real   ("tolerance"        , tol  )

! get the maximum number of Newton-Raphson method iterations
!
    call get_parameter_integer("nr_maxit"         , nrmax)
    call get_parameter_integer("nr_extra"         , nrext)

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

      case("1Dw", "1dw", "1DW", "1D(w)", "1D(W)")

! the type of equation of state
!
        name_c2p  = "1D(W)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srhd_adi_1dw

      case("2dwv", "2Dwv", "2D(w,v)", "2D(W,v)")

! the type of equation of state
!
        name_c2p  = "2D(W,v²)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srhd_adi_2dwv

      case("2dwu", "2Dwu", "2D(w,u)", "2D(W,u)")

! the type of equation of state
!
        name_c2p  = "2D(W,u²)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srhd_adi_2dwu

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

!--- SPECIAL RELATIVITY MAGNETOHYDRODYNAMICS ---
!
    case("srmhd", "SRMHD")

! the name of equation system
!
      name_eqsys = "Special Relativity MHD"

! set relativistic flag
!
      relativistic = .true.

! initialize the number of variables (density + 3 components of velocity
!                                             + 3 components of magnetic field
!                                             + magnetic divergence potential)
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
        prim2cons       => prim2cons_srmhd_adi
        cons2prim       => cons2prim_srmhd_adi
        fluxspeed       => fluxspeed_srmhd_adi
        maxspeed        => maxspeed_srmhd_adi

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

      case("1Dw", "1dw", "1DW", "1D(w)", "1D(W)")

! the type of equation of state
!
        name_c2p  = "1D(W)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srmhd_adi_1dw

      case("2dwv", "2Dwv", "2D(w,v)", "2D(W,v)")

! the type of equation of state
!
        name_c2p  = "2D(W,v²)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srmhd_adi_2dwv

      case("2dwu", "2Dwu", "2D(w,u)", "2D(W,u)")

! the type of equation of state
!
        name_c2p  = "2D(W,u²)"

! set pointer to the conversion method
!
        nr_iterate => nr_iterate_srmhd_adi_2dwu

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

! allocate array for the boundary values
!
    allocate(qpbnd(nv,3,2))

! set the boundary values
!
    do p = 1, nv

! set the initial boundary values (1.0 for density and pressure, 0.0 otherwise)
!
      if (pvars(p) == "dens" .or. pvars(p) == "pres") then
        qpbnd(p,:,:) = 1.0d+00
      else
        qpbnd(p,:,:) = 0.0d+00
      end if

! read the boundary values from the parameter file
!
      call get_parameter_real(pvars(p) // "_bnd_xl", qpbnd(p,1,1))
      call get_parameter_real(pvars(p) // "_bnd_xr", qpbnd(p,1,2))
      call get_parameter_real(pvars(p) // "_bnd_yl", qpbnd(p,2,1))
      call get_parameter_real(pvars(p) // "_bnd_yr", qpbnd(p,2,2))
      call get_parameter_real(pvars(p) // "_bnd_zl", qpbnd(p,3,1))
      call get_parameter_real(pvars(p) // "_bnd_zr", qpbnd(p,3,2))

    end do ! over all variables

! allocate space for Roe eigenvectors
!
    allocate(evroe(2,nv,nv))

! get the minimum allowed density and pressure in the system, and the maximum
! Lorentz factor for special relativity
!
    call get_parameter_real("dmin"  , dmin  )
    call get_parameter_real("pmin"  , pmin  )
    call get_parameter_real("lmax"  , lmax  )

! calculate the maximum speed corresponding to the maximum Lorentz factor
!
    vmax = 1.0d+00 - 1.0d+00 / lmax**2

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

! deallocate boundary values array
!
    if (allocated(qpbnd)) deallocate(qpbnd)

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
    real(kind=8) :: vv, vm, vs, ww
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
    logical      :: info
    integer      :: i
    real(kind=8) :: mm, bb, mb, en, dn
    real(kind=8) :: w , vv, vm, vs
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

! find the exact W using an Newton-Ralphson interative method
!
      call nr_iterate(mm, bb, mb, en, dn, w, vv, info)

! if info is .true., the solution was found
!
      if (info) then

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
        q(ipr,i) = w - en

      else ! unphysical state

        write(*,"(a,1x,a)"           ) "ERROR in"                              &
                                     , "EQUATIONS::cons2prim_srhd_adi()"
        write(*,"(a,5(1x,1e24.16e3))") "Unphysical state for U = ", u(1:nv,i)
        write(*,"(a,3(1x,1e24.16e3))") "            D, |m|², E = ", dn, mm, en
        stop

      end if ! unphysical state

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
! subroutine NR_FUNCTION_SRHD_ADI_1D:
! ----------------------------------
!
!   Subroutine calculate the value of function
!
!     F(W) = W - P(W) - E
!
!   for a given enthalpy W. It is used to estimate the initial guess.
!
!   The pressure is
!
!     P(W) = (γ - 1)/γ (W - D / sqrt(1 - |v|²(W))) (1 - |v|²(W))
!
!   and the squared velocity is
!
!     |v|²(W) = |m|² / W²
!
!   Optional derivative is returned
!
!     dF(W)/dW = 1 - dP(W)/dW
!
!   Arguments:
!
!     mm, en, dn, w - input coefficients for |M|² E, D, and W, respectively;
!     f             - the value of function F(W);
!     df            - optional derivative F'(W);
!
!===============================================================================
!
  subroutine nr_function_srhd_adi_1d(mm, en, dn, w, f, df)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8)          , intent(in)    :: mm, en, dn, w
    real(kind=8)          , intent(out)   :: f
    real(kind=8), optional, intent(out)   :: df

! local variables
!
    real(kind=8) :: vv, vm, vs, pr
!
!-------------------------------------------------------------------------------
!
    vv = mm / (w * w)
    vm = 1.0d+00 - vv
    vs = sqrt(vm)
    f  = (1.0d+00 - gammaxi * vm) * w + gammaxi * dn * vs - en
    if (present(df)) then
      df  = 1.0d+00 - gammaxi * (1.0d+00 + (1.0d+00 - dn / (vs * w)) * vv)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine nr_function_srhd_adi_1d
!
!===============================================================================
!
! subroutine NR_ITERATE_SRHD_ADI_1DW:
! ----------------------------------
!
!   Subroutine finds a root W of equation
!
!     F(W) = W - P(W) - E = 0
!
!   using the Newton-Raphson 1Dw iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w, vv  - input/output coefficients W and |v|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
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
  subroutine nr_iterate_srhd_adi_1dw(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: f, df, dw
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! prepare the initial brackets
!
    wl   = sqrt(mm + dn * dn) + gammaxi * pmin
    wu   = en + pmin

! make sure that the upper bracket is larger than the lower one
!
    keep = wl >= wu
    it   = nrmax
    do while(keep)
      wu   = 2.0d+00 * wu
      it   = it - 1
      keep = (wl >= wu) .and. it > 0
    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_1dw()"
      write(*,"(a)"     ) "Could not find the upper limit for enthalpy!"
      info = .false.
      return
    end if

! check if the brackets bound the root region, if not proceed until
! opposite function signs are found for the brackets
!
    call nr_function_srhd_adi_1d(mm, en, dn, wl, fl)
    call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

    keep = (fl * fu > 0.0d+00)
    it   = nrmax

    do while (keep)

      wl = wu
      fl = fu
      wu = 2.0d+00 * wu

      call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

      it   = it - 1

      keep = (fl * fu > 0.0d+00) .and. it > 0

    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_1dw()"
      write(*,"(a)"     ) "No initial brackets found!"
      info = .false.
      return
    end if

! estimate the value of enthalpy close to the root and corresponding v²
!
    w    = wl - fl * (wu - wl) / (fu - fl)

! initialize iteration parameters
!
    info = .true.
    keep = .true.
    it   = nrmax
    cn   = nrext

! iterate using the Newton-Raphson method in order to find a root w of the
! function
!
    do while(keep)

! calculate F(W) and dF(W)/dW
!
      call nr_function_srhd_adi_1d(mm, en, dn, w, f, df)

! update brackets
!
      if (f > fl .and. f < 0.0d+00) then
        wl = w
        fl = f
      end if
      if (f < fu .and. f > 0.0d+00) then
        wu = w
        fu = f
      end if

! calculate the increment dW
!
      dw  = f / df

! update the solution
!
      w   = w - dw

! calculate the error
!
      err = abs(dw / w)

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! if new W lays out of the brackets, use the bisection method to estimate
! the new guess
!
      if (w < wl .or. w > wu) then
        w = 0.5d+00 * (wl + wu)
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! calculate |V|² from W
!
    vv  = mm / (w * w)

! print information about failed convergence or unphysical variables
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srhd_adi_1dw()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srhd_adi_1dw
!
!===============================================================================
!
! subroutine NR_ITERATE_SRHD_ADI_2DWV:
! -----------------------------------
!
!   Subroutine finds a root (W,v²) of 2D equations
!
!     F(W,v²) = W - E - P(W,v²) = 0
!     G(W,v²) = W² v² - m²      = 0
!
!   using the Newton-Raphson 2D iterative method.
!
!   All evaluated equations incorporate already the pressure of the form
!
!     P(W,|v|²) = (γ - 1)/γ (W - Γ D) (1 - |v|²)
!
!   in order to optimize calculations.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w, vv  - input/output coefficients W and |v|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
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
  subroutine nr_iterate_srhd_adi_2dwv(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: ww, vm, vs, gd, gv
    real(kind=8) :: f, dfw, dfv, df
    real(kind=8) :: g, dgw, dgv, dg
    real(kind=8) :: det, jfw, jfv, jgw, jgv
    real(kind=8) :: dw, dv
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! prepare the initial brackets
!
    wl   = sqrt(mm + dn * dn) + gammaxi * pmin
    wu   = en + pmin

! make sure that the upper bracket is larger than the lower one
!
    keep = wl >= wu
    it   = nrmax
    do while(keep)
      wu   = 2.0d+00 * wu
      it   = it - 1
      keep = (wl >= wu) .and. it > 0
    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_1dw()"
      write(*,"(a)"     ) "Could not find the upper limit for enthalpy!"
      info = .false.
      return
    end if

! check if the brackets bound the root region, if not proceed until
! opposite function signs are found for the brackets
!
    call nr_function_srhd_adi_1d(mm, en, dn, wl, fl)
    call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

    keep = (fl * fu > 0.0d+00)
    it   = nrmax

    do while (keep)

      wl = wu
      fl = fu
      wu = 2.0d+00 * wu

      call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

      it   = it - 1

      keep = (fl * fu > 0.0d+00) .and. it > 0

    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_2dwv()"
      write(*,"(a)"     ) "No initial brackets found!"
      info = .false.
      return
    end if

! estimate the value of enthalpy close to the root and corresponding v²
!
    w  = wl - fl * (wu - wl) / (fu - fl)
    vv = mm / (w * w)

! initialize iteration parameters
!
    info = .true.
    keep = .true.
    it   = nrmax
    cn   = nrext

! iterate using the Newton-Raphson method in order to find the roots W and |v|²
! of functions
!
    do while(keep)

! calculate W², (1 - |v|²), and the Lorentz factor
!
      ww  = w * w
      vm  = 1.0d+00 - vv
      vs  = sqrt(vm)
      gd  = gammaxi * dn
      gv  = 1.0d+00 - gammaxi * vm

! calculate F(W,|v|²) and G(W,|v|²)
!
      f   = gv * w - en + gd * vs
      g   = vv * ww - mm

! calculate dF(W,|v|²)/dW and dF(W,|v|²)/d|v|²
!
      dfw = gv
      dfv = gammaxi * w - 0.5d+00 * gd / vs

! calculate dG(W,|v|²)/dW and dG(W,|v|²)/d|v|²
!
      dgw = 2.0d+00 * vv * w
      dgv = ww

! invert the Jacobian J = | dF/dW, dF/d|v|² |
!                         | dG/dW, dG/d|v|² |
!
      det = dfw * dgv - dfv * dgw

      jfw =   dgv / det
      jgw = - dfv / det
      jfv = - dgw / det
      jgv =   dfw / det

! calculate increments dW and d|v|²
!
      dw  = f * jfw + g * jgw
      dv  = f * jfv + g * jgv

! correct W and |v|²
!
      w   = w  - dw
      vv  = vv - dv

! check if the new enthalpy and velocity are physical
!
      if (w < wl) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srhd_adi_2dwv()"
        write(*,"(a,1x,2e24.16e3)") "Enthalpy smaller than the limit: ", w, wl
        info = .false.
        return
      end if
      if (vv < 0.0d+00 .or. vv >= 1.0d+00) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srhd_adi_2dwv()"
        write(*,"(a,1x,1e24.16e3)") "Unphysical speed |v|²: ", vv
        info = .false.
        return
      end if

! calculate the error
!
      err = max(abs(dw / w), abs(dv))

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! print information about failed convergence or unphysical variables
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srhd_adi_2dwv()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srhd_adi_2dwv
!
!===============================================================================
!
! subroutine NR_ITERATE_SRHD_ADI_2DWU:
! -----------------------------------
!
!   Subroutine finds a root (W,u²) of 2D equations
!
!     F(W,u²) = (W - E - P(W,u²)) (u² + 1) = 0
!     G(W,u²) =  W² u² - (u² + 1) m²       = 0
!
!   using the Newton-Raphson 2D iterative method.
!
!   All evaluated equations incorporate already the pressure of the form
!
!     P(W,|u|²) = (γ - 1)/γ (W - Γ D) / (1 + |u|²)
!
!   in order to optimize calculations.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w, vv  - input/output coefficients W and |v|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
!
!===============================================================================
!
  subroutine nr_iterate_srhd_adi_2dwu(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: ww, uu, up, gm, gd
    real(kind=8) :: f, dfw, dfu, df
    real(kind=8) :: g, dgw, dgu, dg
    real(kind=8) :: det, jfw, jfu, jgw, jgu
    real(kind=8) :: dw, du
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! prepare the initial brackets
!
    wl   = sqrt(mm + dn * dn) + gammaxi * pmin
    wu   = en + pmin

! make sure that the upper bracket is larger than the lower one
!
    keep = wl >= wu
    it   = nrmax
    do while(keep)
      wu   = 2.0d+00 * wu
      it   = it - 1
      keep = (wl >= wu) .and. it > 0
    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_1dw()"
      write(*,"(a)"     ) "Could not find the upper limit for enthalpy!"
      info = .false.
      return
    end if

! check if the brackets bound the root region, if not proceed until
! opposite function signs are found for the brackets
!
    call nr_function_srhd_adi_1d(mm, en, dn, wl, fl)
    call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

    keep = (fl * fu > 0.0d+00)
    it   = nrmax

    do while (keep)

      wl = wu
      fl = fu
      wu = 2.0d+00 * wu

      call nr_function_srhd_adi_1d(mm, en, dn, wu, fu)

      it   = it - 1

      keep = (fl * fu > 0.0d+00) .and. it > 0

    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srhd_adi_2dwu()"
      write(*,"(a)"     ) "No initial brackets found!"
      info = .false.
      return
    end if

! estimate the value of enthalpy close to the root and corresponding u²
!
    w  = wl - fl * (wu - wl) / (fu - fl)
    uu = mm / (w * w - mm)

! initialize iteration parameters
!
    info = .true.
    keep = .true.
    it   = nrmax
    cn   = nrext

! iterate using the Newton-Raphson method in order to find the roots W and |u|²
! of functions
!
    do while(keep)

! calculate W², (1 + |u|²), and the Lorentz factor, and some repeated
! expressions
!
      ww  = w * w
      up  = 1.0d+00 + uu
      gm  = sqrt(up)
      gd  = gammaxi * dn

! calculate F(W,|u|²) and G(W,|u|²)
!
      f   = (up - gammaxi) * w - up * en + gm * gd
      g   = uu * ww - up * mm

! calculate dF(W,|u|²)/dW and dF(W,|u|²)/d|u|²
!
      dfw = up - gammaxi
      dfu = w - en + 0.5d+00 * gd / gm

! calculate dG(W,|u|²)/dW and dG(W,|u|²)/d|u|²
!
      dgw = 2.0d+00 * uu * w
      dgu = ww - mm

! invert the Jacobian J = | dF/dW, dF/d|u|² |
!                         | dG/dW, dG/d|u|² |
!
      det = dfw * dgu - dfu * dgw

      jfw =   dgu / det
      jgw = - dfu / det
      jfu = - dgw / det
      jgu =   dfw / det

! calculate increments dW and d|u|²
!
      dw  = f * jfw + g * jgw
      du  = f * jfu + g * jgu

! correct W and |u|²
!
      w   = w  - dw
      uu  = uu - du

! check if the new enthalpy gives physical pressure and velocity
!
      if (w < wl) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srhd_adi_2dwu()"
        write(*,"(a,1x,2e24.16e3)") "Enthalpy smaller than the limit: ", w, wl
        info = .false.
        return
      end if
      if (uu < 0.0d+00) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srhd_adi_2dwu()"
        write(*,"(a,1x,1e24.16e3)") "Unphysical speed |u|²: ", uu
        info = .false.
        return
      end if

! calculate the error
!
      err = max(abs(dw / w), abs(du))

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! calculate |v|² from |u|²
!
    vv = uu / (1.0d+00 + uu)

! print information about failed convergence or unphysical variables
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srhd_adi_2dwu()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srhd_adi_2dwu
!
!*******************************************************************************
!
! ADIABATIC SPECIAL RELATIVITY MAGNETOHYDRODYNAMIC EQUATIONS
!
!*******************************************************************************
!
!===============================================================================
!
! subroutine PRIM2CONS_SRMHD_ADI:
! ------------------------------
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
  subroutine prim2cons_srmhd_adi(n, q, u)

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
    real(kind=8) :: vv, bb, vb, vm, vs, ww, wt
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

! calculate the square of velocity, the quare of magnetic field, the scalar
! product of velocity and magnetic field, the Lorentz factor, specific and
! total enthalpies
!
      vv  = sum(q(ivx:ivz,i) * q(ivx:ivz,i))
      bb  = sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      vb  = sum(q(ivx:ivz,i) * q(ibx:ibz,i))
      vm  = 1.0d+00 - vv
      vs  = sqrt(vm)
      ww  = (q(idn,i) + q(ipr,i) / gammaxi) / vm
      wt  = ww + bb

! calculate conservative variables
!
      u(idn,i) =      q(idn,i) / vs
      u(imx,i) = wt * q(ivx,i) - vb * q(ibx,i)
      u(imy,i) = wt * q(ivy,i) - vb * q(iby,i)
      u(imz,i) = wt * q(ivz,i) - vb * q(ibz,i)
      u(ibx,i) = q(ibx,i)
      u(iby,i) = q(iby,i)
      u(ibz,i) = q(ibz,i)
      u(ibp,i) = q(ibp,i)
      u(ien,i) = wt - q(ipr,i) - u(idn,i) - 0.5d+00 * (vm * bb + vb * vb)

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons_srmhd_adi
!
!===============================================================================
!
! subroutine CONS2PRIM_SRMHD_ADI:
! ------------------------------
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
  subroutine cons2prim_srmhd_adi(n, u, q)

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
    logical      :: info
    integer      :: i
    real(kind=8) :: mm, mb, bb, en, dn
    real(kind=8) :: w, wt, vv, vm, vs, vb, fc
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
! (|B|², |M|² and B.M and their multiplications)
!
      mm  = sum(u(imx:imz,i) * u(imx:imz,i))
      mb  = sum(u(imx:imz,i) * u(ibx:ibz,i))
      bb  = sum(u(ibx:ibz,i) * u(ibx:ibz,i))
      en  = u(ien,i) + u(idn,i)
      dn  = u(idn,i)

! find the exact W using an Newton-Ralphson interative method
!
      call nr_iterate(mm, bb, mb, en, dn, w, vv, info)

! if info is .true., the solution was found
!
      if (info) then

! prepare coefficients
!
        vm = 1.0d+00 - vv
        vs = sqrt(vm)
        wt = w + bb
        fc = mb / w

! calculate the primitive variables
!
        q(idn,i) = dn * vs
        q(ivx,i) = (u(imx,i) + fc * u(ibx,i)) / wt
        q(ivy,i) = (u(imy,i) + fc * u(iby,i)) / wt
        q(ivz,i) = (u(imz,i) + fc * u(ibz,i)) / wt
        q(ibx,i) = u(ibx,i)
        q(iby,i) = u(iby,i)
        q(ibz,i) = u(ibz,i)
        q(ibp,i) = u(ibp,i)
        q(ipr,i) = w - en + 0.5d+00 * (bb + (bb * mm - mb * mb) / wt**2)

! check if the pressure is positive, if not, print a warning and replace it
! with the minimum allowed value pmin
!
        if (q(ipr,i) <= 0.0d+00) then

          write(*,*)
          write(*,"(a,1x,a)"           ) "WARNING in"                          &
                                       , "EQUATIONS::cons2prim_srmhd_adi()"
          write(*,"(a,9(1x,1e24.16e3))") "Negative pressure for U = ", u(1:nv,i)
          write(*,"(a,6(1x,1e24.16e3))") " D, |m|², m.B, |B|², E, W = "        &
                                                       , dn, mm, mb, bb, en, w
          write(*,"(a,1(1x,1e24.16e3))") "Pressure corrected to ", pmin
          q(ipr,i) = pmin

        end if ! p <= 0

      else ! unphysical state

        write(*,*)
        write(*,"(a,1x,a)"           ) "ERROR in"                              &
                                     , "EQUATIONS::cons2prim_srmhd_adi()"
        write(*,"(a,9(1x,1e24.16e3))") "Unphysical state for U = ", u(1:nv,i)
        write(*,"(a,5(1x,1e24.16e3))") " D, |m|², m.B, |B|², E = ", dn, mm, mb &
                                                                  , bb, en

      end if ! unphysical state

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for variable conversion
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim_srmhd_adi
!
!===============================================================================
!
! subroutine FLUXSPEED_SRMHD_ADI:
! ------------------------------
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
!         "An HLLC Riemann solver for relativistic flows -
!          II. Magnetohydrodynamics",
!         Monthly Notices of the Royal Astronomical Society,
!         2006, Volume 368, Pages 1040-1054
!     [2] van der Holst, B., Keppens, R., Meliani, Z.
!         "A multidimentional grid-adaptive relativistic magnetofluid code",
!         Computer Physics Communications, 2008, Volume 179, Pages 617-627
!
!===============================================================================
!
  subroutine fluxspeed_srmhd_adi(n, q, u, f, cm, cp)

! include external procedures
!
    use algebra        , only : quadratic, quartic

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
    integer      :: i, nr
    real(kind=8) :: vv, bb, vb, vm, vs
    real(kind=8) :: bx, by, bz, b2, pm, pt
    real(kind=8) :: rh, v1, v2
    real(kind=8) :: ca, cc, c2, gn, rt, zm, zp
    real(kind=8) :: fa, fb, fc, fd, fe, ff, fg

! local arrays for characteristic speeds
!
    real(kind=8), dimension(5) :: a
    real(kind=8), dimension(4) :: x
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

! calculate the square of velocity, magnetic field and their scalar product
!
      vv  = sum(q(ivx:ivz,i) * q(ivx:ivz,i))
      bb  = sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      vb  = sum(q(ivx:ivz,i) * q(ibx:ibz,i))

! calculate (1 - |V|²)
!
      vm  = 1.0d+00 - vv

! calculate magnetic field components of the magnetic four-vector divided by
! the Lorentz factor (eq. 3 in [1])
!
      bx  = q(ibx,i) * vm + vb * q(ivx,i)
      by  = q(iby,i) * vm + vb * q(ivy,i)
      bz  = q(ibz,i) * vm + vb * q(ivz,i)

! calculate magnetic and total pressures (eq. 6 in [1])
!
      b2  = bb * vm + vb * vb
      pm  = 0.5d+00 * b2
      pt  = q(ipr,i) + pm

! calculate the relativistic hydrodynamic fluxes (eq. 13 in [1])
!
      f(idn,i) = u(idn,i) * q(ivx,i)
      f(imx,i) = u(imx,i) * q(ivx,i) - q(ibx,i) * bx + pt
      f(imy,i) = u(imy,i) * q(ivx,i) - q(ibx,i) * by
      f(imz,i) = u(imz,i) * q(ivx,i) - q(ibx,i) * bz
      f(ibx,i) = q(ibp,i)
      f(ibp,i) = cmax2 * q(ibx,i)
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)
      f(ien,i) = u(imx,i) - f(idn,i)

! calculate the fast magnetosonic speed
!
! check if the total velocity |V|² is larger than zero
!
      if (vv > 0.0d+00) then

! calculate additional coefficients
!
        rh  = q(idn,i) + q(ipr,i) / gammaxi
        vs  = sqrt(vm)

! check if the normal component of magnetic field Bₓ is larger than zero
!
        if (q(ibx,i) /= 0.0d+00) then ! Bₓ ≠ 0

! prepare parameters for this case
!
          c2 = gamma * q(ipr,i) / rh
          v1 = abs(q(ivx,i))
          v2 = v1 * v1

          fa = rh * (1.0d+00 - c2)
          fb = c2 * (rh - vb * vb) + b2
          fc = sign(1.0d+00, q(ivx,i)) * q(ibx,i) * vs
          fd = c2 * fc * fc
          fe = 1.0d+00 - v2
          ff = c2 * vb * fc
          fg = v1 * vs

! prepare polynomial coefficients
!
          a(5) = fa + fb * vm
          a(4) = 2.0d+00 * (ff * vm + fb * fg)
          a(3) = - fd * vm + 4.0d+00 * ff * fg - fb * fe
          a(2) = - 2.0d+00 * (fd * fg + fe * ff)
          a(1) = fd * fe

! call the quartic solver
!
          nr = quartic(a(1:5), x(1:4))

! convert eigenvalues to charasteristic speeds
!
          x(1:nr) = sign(1.0d+00, q(ivx,i)) * (abs(v1) + x(1:nr) * vs)

        else ! Bₓ ≠ 0

! special case when Bₓ = 0, then the quartic equation reduces to quadratic one
!
! prepare parameters for this case
!
          c2 = gamma * q(ipr,i) / rh
          cc = (1.0d+00 - c2) / vm
          gn = b2 - c2 * vb * vb

! prepare polynomial coefficients
!
          a(3) = rh * (c2 + cc) + gn
          a(2) = - 2.0d+00 * rh * cc * q(ivx,i)
          a(1) = rh * (cc * q(ivx,i)**2 - c2) - gn

! solve the quadratic equation
!
          nr = quadratic(a(1:3), x(1:2))

        end if ! Bx ≠ 0

      else ! |V|² > 0

! special case when |V|² = 0 (Γ = 1), then the quartic equation reduces to
! bi-quartic one
!
! prepare parameters for this case
!
        rh  = q(idn,i) + q(ipr,i) / gammaxi
        vs  = sqrt(vm)
        bx  = q(ibx,i) * vs + vb * q(ivx,i) / vs
        rt  = rh + b2
        c2 = gamma * q(ipr,i) / rh
        ca = bx * bx

! prepare polynomial coefficients
!
        a(3) = 1.0d+00
        a(2) = - ((rh + ca) * c2 + b2) / rt
        a(1) = c2 * ca / rt

! solve the bi-quartic equation
!
        nr = quadratic(a(1:3), x(1:2))

! compute the roots
!
        if (nr > 0) then

          zm = min(x(1), x(2))
          zp = max(x(1), x(2))

          if (zm >= 0.0d+00) then

            zm   = sqrt(zm)
            zp   = sqrt(zp)

            x(1) =   zp
            x(2) =   zm
            x(3) = - zm
            x(4) = - zp

            nr   = 4

          else

            if (zp >= 0.0d+00) then

              zp   = sqrt(zp)

              x(1) =   zp
              x(2) = - zp

              nr = 2

            else

              x(:) = 0.0d+00

              nr   = 0

            end if
          end if
        end if

      end if ! |V|² > 0

! find the minimum and maximum characteristic speeds
!
      if (nr > 1) then

        cm(i) = minval(x(1:nr))
        cp(i) = maxval(x(1:nr))

#ifdef DEBUG
        if (max(abs(cm(i)), abs(cp(i))) >= 1.0d+00) then
          write(*,*)
          write(*,*) 'Estimation returned unphysical speeds!'
          write(*,"('A = ',5(1pe24.16))") a(1:5)
          write(*,"('N = ',1i2)"        ) nr
          write(*,"('X = ',4(1pe24.16))") x(1:4)
          stop
        end if
#endif /* DEBUG */

      else

! speed estimation failed, so we substitute the minimum and maximum physical
! speeds equal the speed of light
!
        cm(i) = - 1.0d+00
        cp(i) =   1.0d+00

#ifdef DEBUG
        write(*,*)
        write(*,*) 'Speed estimation failed!'
        write(*,"('A = ',5(1pe24.16))") a(1:5)
        write(*,"('N = ',1i2)"        ) nr
        write(*,"('X = ',4(1pe24.16))") x(1:4)
        stop
#endif /* DEBUG */

      end if

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for flux calculation
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed_srmhd_adi
!
!===============================================================================
!
! function MAXSPEED_SRMHD_ADI:
! ---------------------------
!
!   Function scans the variable array and returns the maximum speed in within.
!
!   Arguments:
!
!     qq - the array of primitive variables;
!
!===============================================================================
!
  function maxspeed_srmhd_adi(qq) result(maxspeed)

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
    maxspeed = 1.0d+00

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
  end function maxspeed_srmhd_adi
!
!===============================================================================
!
! subroutine NR_POSITIVITY_SRMHD_ADI_1D:
! -------------------------------------
!
!   Subroutine calculates the function
!
!     Q(W)  = W² - D² - [|m|² W² + (2 W + |B|²) S²] / (W + |B|²)²
!
!   and its derivative
!
!     Q'(W) = 2 W [ 1 - (|m|² |B|² - S²) / (W + |B|²)³]
!
!   for a given enthalpy W.
!
!   This subroutine is used to find the minimum enthalpy for which the velocity
!   is physical and the pressure is positive.
!
!   Arguments:
!
!     mm, bb, mb, dn, w - input coefficients for |M|², |B|², M.B, D, and W,
!                         respectively;
!     q, dq             - the values for the function Q(W) and its
!                         derivative Q'(W);
!
!===============================================================================
!
  subroutine nr_positivity_srmhd_adi_1d(mm, bb, mb, dn, w, q, dq)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)  :: mm, bb, mb, dn, w
    real(kind=8), intent(out) :: q, dq

! local variables
!
    real(kind=8) :: dd, ss, ww, wt
!
!-------------------------------------------------------------------------------
!
! temporary variables
!
    dd = dn * dn
    ss = mb * mb
    ww = w * w
    wt = w + bb

! the function and its derivative
!
    q  = ww - dd - (mm * ww + (w + wt) * ss) / wt**2
    dq = 2.0d+00 * w * (1.0d+00 - (mm * bb - ss) / wt**3)

!-------------------------------------------------------------------------------
!
  end subroutine nr_positivity_srmhd_adi_1d
!
!===============================================================================
!
! subroutine NR_VELOCITY_SRMHD_ADI_1D:
! -----------------------------------
!
!   Subroutine calculates the squared velocity
!
!     |V|²(W)  = (|m|² W² + S² (2 W + |B|²)) / (W² (W² + |B|²)²)
!
!   and its derivative
!
!     |V|²'(W) = - 2 (|m|² W³ + 3 S² W² + 3 S² |B|² W + S² |B|⁴)
!                                                           / (W³ (W² + |B|²)³)
!
!   for a given enthalpy W.
!
!   Arguments:
!
!     mm, bb, mb, w - input coefficients for |M|², |B|², M.B, and W,
!                     respectively;
!     vv, dv        - the values of squared velocity |V|² and its derivative;
!
!===============================================================================
!
  subroutine nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv, dv)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8),           intent(in)  :: mm, bb, mb, w
    real(kind=8),           intent(out) :: vv
    real(kind=8), optional, intent(out) :: dv

! local variables
!
    real(kind=8) :: ss, ww, www, wt, wt2
!
!-------------------------------------------------------------------------------
!
! temporary variables
!
    ss  = mb * mb
    ww  = w * w
    www = ww * w
    wt  = w + bb
    wt2 = wt * wt

! the function and its derivative
!
    vv = (mm * ww + (w + wt) * ss) / (ww * wt2)
    if (present(dv)) then
      dv = - 2.0d+00 * (mm * www + (3.0d+00 * ww                               &
                           + (2.0d+00 * w + wt) * bb) * ss) / (www * wt2 * wt)
    end if

!-------------------------------------------------------------------------------
!
  end subroutine nr_velocity_srmhd_adi_1d
!
!===============================================================================
!
! subroutine NR_PRESSURE_SRMHD_ADI_1D:
! -----------------------------------
!
!   Subroutine calculates the pressure function
!
!     P(W)  = W  - E + ½ |B|² + ½ (|m|² |B|² - S²) / (W + |B|²)²
!
!   and its derivative
!
!     P'(W) = 1 - (|m|² |B|² - S²) / (W + |B|²)³
!
!   for a given enthalpy W.
!
!   This subroutine is used to find the minimum enthalpy for which the velocity
!   is physical and the pressure is positive.
!
!   Arguments:
!
!     mm, bb, mb, dn, en, w - input coefficients for |M|², |B|², M.B, D, E,
!                             and W, respectively;
!     p, dp                 - the values for the function P(W) and its
!                             derivative P'(W);
!
!===============================================================================
!
  subroutine nr_pressure_srmhd_adi_1d(mm, bb, mb, en, dn, w, p, dp)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8)          , intent(in)  :: mm, bb, mb, en, dn, w
    real(kind=8)          , intent(out) :: p
    real(kind=8), optional, intent(out) :: dp

! local variables
!
    real(kind=8) :: wt, wd, ss, fn
!
!-------------------------------------------------------------------------------
!
! temporary variables
!
    wt = w  + bb
    wd = wt * wt
    ss = mb * mb
    fn = (mm * bb - ss) / wd

! the pressure function and its derivative
!
    p  = w - en + 0.5d+00 * (bb + fn)
    if (present(dp)) dp = 1.0d+00 - fn / wt

!-------------------------------------------------------------------------------
!
  end subroutine nr_pressure_srmhd_adi_1d
!
!===============================================================================
!
! subroutine NR_FUNCTION_SRMHD_ADI_1D:
! -----------------------------------
!
!   Subroutine calculates the energy function
!
!     F(W)  = W - P(W) + ½ |B|² + ½ (|m|² |B|² - S²) / (W + |B|²)² - E
!
!   and its derivative
!
!     F'(W) = 1 - dP(W)/dW - (|m|² |B|² - S²) / (W + |B|²)³
!
!   for a given enthalpy W. It is used to estimate the initial guess.
!
!   Arguments:
!
!     mm, bb, mb, en, dn, w - input coefficients for |M|², |B|², M.B, E, D,
!                             and W, respectively;
!     f, df                 - the values of F(W) and its derivative;
!
!===============================================================================
!
  subroutine nr_function_srmhd_adi_1d(mm, bb, mb, en, dn, w, f, df)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8),           intent(in)  :: mm, bb, mb, en, dn, w
    real(kind=8),           intent(out) :: f
    real(kind=8), optional, intent(out) :: df

! local variables
!
    real(kind=8) :: pr, dp, gm2, gm, dg
    real(kind=8) :: ww, wt, wd, ss, sw, ws, fn, dv, ds
!
!-------------------------------------------------------------------------------
!
! temporary variables
!
    ww  = w * w
    wt  = w + bb
    wd  = wt * wt
    ss  = mb * mb
    sw  = ss / ww
    ws  = (w + wt) * sw
    fn  = (mm * bb - ss) / wd
    dv  = wd - mm - ws
    ds  = sqrt(dv)

! calculate the Lorentz factor
!
    gm2 = wd / dv
    gm  = wt / ds

! calculate the pressure P(W) and energy function F(W)
!
    pr  = gammaxi * (w - gm * dn) / gm2
    f   = w - pr - en + 0.5d+00 * (bb + fn)

! if desired, calculate the derivatives dP(W)/dW and dF(W)/dW
!
    if (present(df)) then

      dg  = (1.0d+00 - wt * (wt - sw + ws / w) / dv) / ds
      dp  = gammaxi * (1.0d+00 - (2.0d+00 * w / gm - dn) * dg) / gm2
      df  = 1.0d+00 - dp - fn / wt

    end if ! df present

!-------------------------------------------------------------------------------
!
  end subroutine nr_function_srmhd_adi_1d
!
!===============================================================================
!
! subroutine NR_INITIAL_BRACKETS_SRMHD_ADI:
! ----------------------------------------
!
!   Subroutine finds the initial brackets and initial guess from
!   the positivity condition
!
!     W³ + (5/2 |B|² - E) W² + 2 (|B|² - E) |B|² W
!                          + 1/2 [(|B|⁴ - 2 |B|² E + |m|²) |B|² - S²] > 0
!
!   coming from the energy equation and
!
!     W⁴ + 2 |B|² W³ - (|m|² + D² - |B|⁴) W²
!                    - (2 S² + D² |B|²) W - (S² + D² |B|²) |B|² > 0
!
!   coming from the equation of state
!
!   using analytical. It takes the maximum estimated root as the lower bracket.
!   If the analytical estimation fails, the Newton-Raphson iterative method
!   is used.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, mb - input coefficients for |B|² and m.B, respectively;
!     wl, wu - the lower and upper limits for the enthalpy;
!     wc     - the initial root guess;
!     info   - the flag is .true. if the initial brackets and guess were found,
!              otherwise it is .false.;
!
!===============================================================================
!
  subroutine nr_initial_brackets_srmhd_adi(mm, bb, mb, en, dn                  &
                                                     , wl, wu, wc, fl, fu, info)

! include external procedures
!
    use algebra        , only : cubic_normalized, quartic

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)  :: mm, bb, mb, en, dn
    real(kind=8), intent(out) :: wl, wu, wc, fl, fu
    logical     , intent(out) :: info

! local variables
!
    logical      :: keep
    integer      :: it, nr
    real(kind=8) :: dd, ss, ec
    real(kind=8) :: f , df
    real(kind=8) :: dw, err

! local vectors
!
    real(kind=8), dimension(5) :: a
    real(kind=8), dimension(4) :: x
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for initial bracket solver
!
    call start_timer(imb)
#endif /* PROFILE */

! calculate temporary variables
!
    dd   = dn * dn
    ss   = mb * mb
    ec   = en + pmin

! set the initial upper bracket
!
    wu   = en + pmin

! calculate the cubic equation coefficients for the positivity condition
! coming from the energy equation; the condition, in fact, finds the minimum
! enthalphy for which the pressure is equal to pmin
!
    a(3) = 2.5d+00 * bb - ec
    a(2) = 2.0d+00 * (bb - ec) * bb
    a(1) = 0.5d+00 * (((bb - 2.0d+00 * ec) * bb + mm) * bb - ss)

! solve the cubic equation
!
    nr = cubic_normalized(a(1:3), x(1:3))

! if solution was found, use the maximum root as the lower bracket
!
    if (nr > 0) then

      wl = x(nr)

! calculate the quartic equation coefficients for the positivity condition
! coming from the pressure equation
!
      a(5) = 1.0d+00
      a(4) = 2.0d+00 * bb
      a(3) = bb * bb - dd - mm
      a(2) = - 2.0d+00 * (ss + dd * bb)
      a(1) = - (ss + dd * bb) * bb

! solve the quartic equation
!
      nr = quartic(a(1:5), x(1:4))

! take the maximum ethalpy from both conditions to guarantee that the pressure
! obtains from any of those equations is positive
!
      if (nr > 0) wl = max(wl, x(nr))

    else ! nr = 0

! the root could not be found analytically, so use the iterative solver
! to find the lower bracket; as the initial guess use the initial upper bracket
!
      keep = .true.
      it   = nrmax
      wl   = wu

      do while(keep)

        call nr_pressure_srmhd_adi_1d(mm, bb, mb, ec, dn, wl, f, df)

        dw   = f / df
        wl   = wl - dw

        err  = abs(dw / wl)
        it   = it - 1
        keep = (err > tol) .and. it > 0

      end do

      if (it <= 0) then
        write(*,*)
        write(*,"(a,1x,a)") "ERROR in"                                         &
                          , "EQUATIONS::nr_initial_brackets_srmhd_adi()"
        write(*,"(a)"     ) "Could not find the lower limit for the enthalpy!"
        write(*,"(a,5(1x,1e24.16e3))") " D, |m|², m.B, |B|², E = "             &
                                                          , dn, mm, mb, bb, en
        info = .false.
        return
      end if

    end if ! nr > 0

! check if the energy function is negative for the lower limit
!
    call nr_function_srmhd_adi_1d(mm, bb, mb, en, dn, wl, fl)

    if (fl > 0.0d+00) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_initial_brackets_srmhd_adi()"
      write(*,"(a)"     ) "Lower limit positive!"
      write(*,"(a,6(1x,1e24.16e3))") " D, |m|², m.B, |B|², E, W = "            &
                                                      , dn, mm, mb, bb, en, wl
      info = .false.
      return
    end if

! make sure that the upper limit is larger than the lower one
!
    keep = wl >= wu
    it   = nrmax
    do while(keep)
      wu   = 2.0d+00 * wu
      it   = it - 1
      keep = (wl >= wu) .and. it > 0
    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a)"     ) "Could not find the upper limit for enthalpy!"
      write(*,"(a,6(1x,1e24.16e3))") " D, |m|², m.B, |B|², E, W = "            &
                                                      , dn, mm, mb, bb, en, wl
      info = .false.
      return
    end if

! check if the brackets bound the root region, if not proceed until
! opposite function signs are found for the brackets
!
    call nr_function_srmhd_adi_1d(mm, bb, mb, en, dn, wu, fu)

    keep = (fl * fu > 0.0d+00)
    it   = nrmax

    do while (keep)

      wl = wu
      fl = fu
      wu = 2.0d+00 * wu

      call nr_function_srmhd_adi_1d(mm, bb, mb, en, dn, wu, fu)

      it   = it - 1
      keep = (fl * fu > 0.0d+00) .and. it > 0

    end do
    if (it <= 0) then
      write(*,*)
      write(*,"(a,1x,a)") "ERROR in"                                           &
                        , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a)"     ) "No initial brackets found!"
      write(*,"(a,5(1x,1e24.16e3))") " D, |m|², m.B, |B|², E = "               &
                                                          , dn, mm, mb, bb, en
      info = .false.
      return
    end if

! estimate the enthalpy value close to the root
!
    wc   = wl - fl * (wu - wl) / (fu - fl)

#ifdef PROFILE
! stop accounting time for the initial brackets
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_initial_brackets_srmhd_adi
!
!===============================================================================
!
! subroutine NR_ITERATE_SRMHD_ADI_1DW:
! -----------------------------------
!
!   Subroutine finds a root W of equation
!
!     F(W) = W - P(W) + ½ [(1 + |V|²) |B|² - S² / W²] - E = 0
!
!   using the Newton-Raphson 1Dw iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w , vv - input/output coefficients W and |V|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
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
  subroutine nr_iterate_srmhd_adi_1dw(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: f , df, dw
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! find the initial brackets and estimate the initial enthalpy
!
    call nr_initial_brackets_srmhd_adi(mm, bb, mb, en, dn                      &
                                                     , wl, wu, w, fl, fu, info)

! if the brackets could not be found, return the lower bracket as the solution
!
    if (.not. info) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a,1x)"          ) "The solution lays in unphysical regime."
      write(*,"(a,1x,1e24.16e3)") "Using the lower bracket as solution: ", wl

! use the lower bracket, since it guarantees the positive pressure
!
      w = wl

! calculate |V|² from W
!
      call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)

      info = .true.
      return

    end if

! initialize iteration parameters
!
    info = .true.
    keep = .true.
    it   = nrmax
    cn   = nrext

! iterate using the Newton-Raphson method in order to find a root w of the
! function
!
    do while(keep)

! calculate F(W) and dF(W)/dW
!
      call nr_function_srmhd_adi_1d(mm, bb, mb, en, dn, w, f, df)

! update brackets
!
      if (f > fl .and. f < 0.0d+00) then
        wl = w
        fl = f
      end if
      if (f < fu .and. f > 0.0d+00) then
        wu = w
        fu = f
      end if

! calculate the increment dW
!
      dw  = f / df

! update the solution
!
      w   = w - dw

! calculate the error
!
      err = abs(dw / w)

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! if new W lays out of the brackets, use the bisection method to estimate
! the new guess
!
      if (w < wl .or. w > wu) then
        w = 0.5d+00 * (wl + wu)
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! calculate |V|² from W
!
    call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)

! let know the user if the convergence failed
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srmhd_adi_1dw
!
!===============================================================================
!
! subroutine NR_ITERATE_SRMHD_ADI_2DWV:
! ------------------------------------
!
!   Subroutine finds a root (W, |V|²) of equations
!
!     F(W,|V|²) = W - P + ½ [(1 + |V|²) |B|² - S² / W²] - E     = 0
!     G(W,|V|²) = |V|² (|B|² + W)² - S² (|B|² + 2W) / W² - |M|² = 0
!
!   using the Newton-Raphson 2D iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w , vv - input/output coefficients W and |V|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
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
  subroutine nr_iterate_srmhd_adi_2dwv(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: vm, vs, wt, mw, wt2
    real(kind=8) :: f, df, dfw, dfv
    real(kind=8) :: g, dg, dgw, dgv
    real(kind=8) :: det, jfw, jfv, jgw, jgv
    real(kind=8) :: dv, dw
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! find the initial brackets and estimate the initial enthalpy
!
    call nr_initial_brackets_srmhd_adi(mm, bb, mb, en, dn                      &
                                                     , wl, wu, w, fl, fu, info)

! if the brackets could not be found, return the lower bracket as the solution
!
    if (.not. info) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a,1x)"          ) "The solution lays in unphysical regime."
      write(*,"(a,1x,1e24.16e3)") "Using the lower bracket as solution: ", wl

! use the lower bracket, since it guarantees the positive pressure
!
      w = wl

! calculate |V|² from W
!
      call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)

      info = .true.
      return

    end if

! and the corresponding |V|²
!
    call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)

! initialize iteration parameters
!
    info  = .true.
    keep  = .true.
    it    = nrmax
    cn    = nrext

! find root with the help of the Newton-Raphson method
!
    do while(keep)

! calculate (S/W)², Wt, Wt²
!
      mw  = (mb / w)**2
      wt  = w + bb
      wt2 = wt * wt

! prepare (1 - |V|²) and sqrt(1 - |V|²)
!
      vm  = 1.0d+00 - vv
      vs  = sqrt(vm)

! calculate functions F(W,|V|²) and G(W,|V|²)
!
      f   = w - en - gammaxi * (w * vm - dn * vs)                              &
                                        + 0.5d+00 * ((1.0d+00 + vv) * bb - mw)
      g   = vv * wt2 - (wt + w) * mw  - mm

! calculate derivatives dF(W,|V|²)/dW and dF(W,|V|²)/d|V|²
!
      dfw = 1.0d+00 - gammaxi * vm + mw / w
      dfv = - gammaxi * (0.5d+00 * dn / vs - w) + 0.5d+00 * bb

! calculate derivatives dG(W,|V|²)/dW and dG(W,|V|²)/d|V|²
!
      dgw = 2.0d+00 * wt * (vv + mw / w)
      dgv = wt2

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

! check if the new enthalpy and velocity are physical
!
      if (w < wl) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srmhd_adi_2dwv()"
        write(*,"(a,1x,2e24.16e3)") "Enthalpy smaller than the limit: ", w, wl
        info = .false.
        return
      end if
      if (vv < 0.0d+00 .or. vv >= 1.0d+00) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srmhd_adi_2dwv()"
        write(*,"(a,1x,1e24.16e3)") "Unphysical speed |v|²: ", vv
        info = .false.
        return
      end if

! calculate the error
!
      err = max(abs(dw / w), abs(dv))

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! let know the user if the convergence failed
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_2dwv()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srmhd_adi_2dwv
!
!===============================================================================
!
! subroutine NR_ITERATE_SRMHD_ADI_2DWU:
! ------------------------------------
!
!   Subroutine finds a root (W, |u|²) of equations
!
!     F(W,|u|²) = W - E - P + ½ [(1 + |u|² / (1 + |u|²)) |B|² - S² / W²]     = 0
!     G(W,|u|²) = (|B|² + W)² |u|² / (1 + |u|²) - (2W + |B|²) S² / W² - |M|² = 0
!
!   using the Newton-Raphson 2D iterative method.
!
!   Arguments:
!
!     mm, en - input coefficients for |M|² and E, respectively;
!     bb, bm - input coefficients for |B|² and B.M, respectively;
!     w , vv - input/output coefficients W and |v|²;
!     info   - the flag is .true. if the solution was found, otherwise
!              it is .false.;
!
!===============================================================================
!
  subroutine nr_iterate_srmhd_adi_2dwu(mm, bb, mb, en, dn, w, vv, info)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), intent(in)    :: mm, bb, mb, en, dn
    real(kind=8), intent(inout) :: w, vv
    logical     , intent(out)   :: info

! local variables
!
    logical      :: keep
    integer      :: it, cn
    real(kind=8) :: wl, wu, fl, fu
    real(kind=8) :: uu, up, gm
    real(kind=8) :: ss, ww, wt, wt2, wd, wp
    real(kind=8) :: f, df, dfw, dfu
    real(kind=8) :: g, dg, dgw, dgu
    real(kind=8) :: det, jfw, jfu, jgw, jgu
    real(kind=8) :: dw, du
    real(kind=8) :: err
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable solver
!
    call start_timer(imp)
#endif /* PROFILE */

! find the initial brackets and estimate the initial enthalpy
!
    call nr_initial_brackets_srmhd_adi(mm, bb, mb, en, dn                      &
                                                     , wl, wu, w, fl, fu, info)

! if the brackets could not be found, return the lower bracket as the solution
!
    if (.not. info) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_1dw()"
      write(*,"(a,1x)"          ) "The solution lays in unphysical regime."
      write(*,"(a,1x,1e24.16e3)") "Using the lower bracket as solution: ", wl

! use the lower bracket, since it guarantees the positive pressure
!
      w = wl

! calculate |V|² from W
!
      call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)

      info = .true.
      return

    end if

! and the corresponding |u|²
!
    call nr_velocity_srmhd_adi_1d(mm, bb, mb, w, vv)
    uu = vv / (1.0d+00 - vv)

! initialize iteration parameters
!
    info  = .true.
    keep  = .true.
    it    = nrmax
    cn    = nrext

! find root with the help of the Newton-Raphson method
!
    do while(keep)

! prepare (1 + |u|²) and the Lorentz factor
!
      up  = 1.0d+00 + uu
      gm  = sqrt(up)

! calculate temporary variables
!
      ss  = mb * mb
      ww  = w  * w
      wt  = w + bb
      wt2 = wt * wt
      wd  = w - gm * dn
      wp  = wt / up

! calculate functions F(W,|u|²) and G(W,|u|²)
!
      f   = w - en - gammaxi * wd / up                                         &
                              + 0.5d+00 * (bb * (1.0d+00 + uu / up) - ss / ww)
      g   = wp * wt * uu - (w + wt) * ss / ww - mm

! calculate derivatives dF(W,|u|²)/dW and dF(W,|u|²)/d|u|²
!
      dfw = 1.0d+00 - gammaxi / up + ss / ww / w
      dfu = 0.5d+00 * (gammaxi * (w + wd) + bb) / up**2

! calculate derivatives dG(W,|u|²)/dW and dG(W,|u|²)/d|u|²
!
      dgw = 2.0d+00 * wt * (uu / up + ss / ww / w)
      dgu = wp * wp

! invert the Jacobian J = | dF/dW, dF/d|u|² |
!                         | dG/dW, dG/d|u|² |
!
      det = dfw * dgu - dfu * dgw

      jfw =   dgu / det
      jgw = - dfu / det
      jfu = - dgw / det
      jgu =   dfw / det

! calculate increments dW and d|u|²
!
      dw  = f * jfw + g * jgw
      du  = f * jfu + g * jgu

! correct W and |u|²
!
      w   = w  - dw
      uu  = uu - du

! check if the new enthalpy and velocity are physical
!
      if (w < wl) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srmhd_adi_2dwu()"
        write(*,"(a,1x,2e24.16e3)") "Enthalpy smaller than the limit: ", w, wl
        info = .false.
        return
      end if
      if (uu < 0.0d+00) then
        write(*,*)
        write(*,"(a,1x,a)"        ) "ERROR in"                                 &
                                  , "EQUATIONS::nr_iterate_srmhd_adi_2dwu()"
        write(*,"(a,1x,1e24.16e3)") "Unphysical speed |u|²: ", uu
        info = .false.
        return
      end if

! calculate the error
!
      err = max(abs(dw / w), abs(du))

! check the convergence, if the convergence is not reached, iterate until
! the maximum number of iteration is reached
!
      if (err < tol) then
        keep = cn > 0
        cn   = cn - 1
      else
        keep = it > 0
      end if

! decrease the number of remaining iterations
!
      it = it - 1

    end do ! NR iterations

! calculate |v|² from |u|²
!
    vv = uu / (1.0d+00 + uu)

! let know the user if the convergence failed
!
    if (err >= tol) then
      write(*,*)
      write(*,"(a,1x,a)"        ) "WARNING in"                                 &
                                , "EQUATIONS::nr_iterate_srmhd_adi_2dwu()"
      write(*,"(a,1x,1e24.16e3)") "Convergence not reached: ", err
    end if

#ifdef PROFILE
! stop accounting time for variable solver
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine nr_iterate_srmhd_adi_2dwu

!===============================================================================
!
end module equations
