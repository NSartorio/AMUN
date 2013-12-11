!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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

! module variables are not implicit by default
!
  implicit none

! pointers to the conversion procedures
!
  procedure(prim2cons_hd_iso), pointer, save :: prim2cons => null()
  procedure(cons2prim_hd_iso), pointer, save :: cons2prim => null()

! pointer to the flux procedure
!
  procedure(fluxspeed_hd_iso), pointer, save :: fluxspeed => null()

! pointer to the maxspeed procedure
!
  procedure(maxspeed_hd_iso) , pointer, save :: maxspeed  => null()


! the system of equations and the equation of state
!
  character(len=32), save :: eqsys = "hydrodynamic"
  character(len=32), save :: eos   = "adiabatic"

! the number of independent variables
!
  integer(kind=4)  , save :: nv  = 0

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

! adiabatic heat ratio
!
  real(kind=8)     , save :: gamma   = 5.0d+00 / 3.0d+00

! additional adiabatic parameters
!
  real(kind=8)     , save :: gammam1 = 2.0d+00 / 3.0d+00, gammam1i = 1.5d+00

! isothermal speed of sound and its second power
!
  real(kind=8)     , save :: csnd    = 1.0d+00, csnd2   = 1.0d+00

! maximum speed in the system
!
  real(kind=8)     , save :: cmax    = 0.0d+00, cmax2   = 0.0d+00

! by default everything is private
!
  private

! declare public variables and subroutines
!
  public :: initialize_equations, finalize_equations
  public :: prim2cons, cons2prim
  public :: fluxspeed
  public :: maxspeed, reset_maxspeed, get_maxspeed
  public :: update_primitive_variables
  public :: gamma
  public :: csnd
  public :: cmax, cmax2
  public :: nv
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
    use parameters, only : get_parameter_string, get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=255)     :: name_eqsys = ""
    character(len=255)     :: name_eos   = ""
!
!-------------------------------------------------------------------------------
!
! get the system of equations
!
    call get_parameter_string("equation_system"   , eqsys)

! get the equation of state
!
    call get_parameter_string("equation_of_state" , eos  )

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
        prim2cons => prim2cons_hd_iso
        cons2prim => cons2prim_hd_iso
        fluxspeed => fluxspeed_hd_iso
        maxspeed  => maxspeed_hd_iso

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
        prim2cons => prim2cons_hd_adi
        cons2prim => cons2prim_hd_adi
        fluxspeed => fluxspeed_hd_adi
        maxspeed  => maxspeed_hd_adi

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
        prim2cons => prim2cons_mhd_iso
        cons2prim => cons2prim_mhd_iso
        fluxspeed => fluxspeed_mhd_iso
        maxspeed  => maxspeed_mhd_iso

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
        prim2cons => prim2cons_mhd_adi
        cons2prim => cons2prim_mhd_adi
        fluxspeed => fluxspeed_mhd_adi
        maxspeed  => maxspeed_mhd_adi

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

! obtain the isothermal sound speed
!
    call get_parameter_real("csnd"  , csnd  )

! calculate additional parameters
!
    csnd2 = csnd * csnd

! print information about the equation module
!
    if (verbose) then

      write (*,"(4x,a,1x,a)"    ) "equation system        =", trim(name_eqsys)
      write (*,"(4x,a,1x,a)"    ) "equation of state      =", trim(name_eos)

    end if

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
! deallocate variable name arrays
!
    if (allocated(pvars)) deallocate(pvars)
    if (allocated(cvars)) deallocate(cvars)

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
    use coordinates, only : im, jm, km

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(nv,im,jm,km), intent(in)    :: uu
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: qq

! temporary variables
!
    integer                        :: j, k

! temporary array to store conserved variable vector
!
    real(kind=8), dimension(nv,im) :: u
!
!-------------------------------------------------------------------------------
!
! update primitive variables
!
    do k = 1, km
      do j = 1, jm

! copy variables to temporary array of conserved variables
!
        u(1:nv,1:im) = uu(1:nv,1:im,j,k)

! convert conserved variables to primitive ones
!
        call cons2prim(im, u(1:nv,1:im), qq(1:nv,1:im,j,k))

      end do ! j = 1, jm
    end do ! k = 1, km

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
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)

    end do ! i = 1, n

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
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)

    end do ! i = 1, n

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
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the input array of conservative variables;
!     f - the output vector of fluxes;
!     c - the output vector of characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_hd_iso(n, q, u, f, c)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q, u
    real(kind=8), dimension(nv,n), intent(out) :: f
    real(kind=8), dimension(n)   , intent(out) :: c

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

! calculate the hydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i)
      f(imy,i) = q(ivx,i) * u(imy,i)
      f(imz,i) = q(ivx,i) * u(imz,i)
      f(imx,i) = f(imx,i) + csnd2 * q(idn,i)

! calculate the speed of sound
!
      c(i) = csnd

! calculate the maximum speed
!
      cmax = max(cmax, abs(q(ivx,i)) + c(i))

    end do ! i = 1, n

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

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_hd_iso
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
    integer :: i
    real    :: ek, ei
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
      ek       = 0.5d+00 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      ei       = gammam1i * q(ipr,i)
      u(ien,i) = ei + ek

    end do ! i = 1, n

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
    integer :: i
    real    :: ek, ei
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
      ek       = 0.5d+00 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      ei       = u(ien,i) - ek
      q(ipr,i) = gammam1 * ei

    end do ! i = 1, n

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
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the input array of conservative variables;
!     f - the output vector of fluxes;
!     c - the output vector of characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_hd_adi(n, q, u, f, c)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q, u
    real(kind=8), dimension(nv,n), intent(out) :: f
    real(kind=8), dimension(n)   , intent(out) :: c

! local variables
!
    integer :: i
!
!-------------------------------------------------------------------------------
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

! calculate the speed of sound
!
      c(i) = sqrt(gamma * q(ipr,i) / q(idn,i))

! calculate the maximum speed
!
      cmax = max(cmax, abs(q(ivx,i)) + c(i))

    end do ! i = 1, n

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

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_hd_adi
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
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the input array of conservative variables;
!     f - the output vector of fluxes;
!     c - the output vector of characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_mhd_iso(n, q, u, f, c)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q, u
    real(kind=8), dimension(nv,n), intent(out) :: f
    real(kind=8), dimension(n)   , intent(out) :: c

! local variables
!
    integer      :: i
    real(kind=8) :: bx2, by2, bz2, bb, pr, pt
    real(kind=8) :: fa, fb, fc
!
!-------------------------------------------------------------------------------
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
      fa   = csnd2 * q(idn,i)
      fb   = fa + bb
      fc   = fb * fb - 4.0d+00 * fa * bx2
      if (fc > 0.0d+00) then
        c(i) = sqrt(0.5d+00 * (fb + sqrt(fc)) / q(idn,i))
      else
        c(i) = sqrt(0.5d+00 *  fb             / q(idn,i))
      end if

! calculate the maximum speed
!
      cmax = max(cmax, abs(q(ivx,i)) + c(i))

    end do ! i = 1, n

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

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_mhd_iso
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
      ek       = 0.5d+00 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      em       = 0.5d+00 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      u(ien,i) = ei + ek + em

    end do ! i = 1, n

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
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
      q(ibx,i) = u(ibx,i)
      q(iby,i) = u(iby,i)
      q(ibz,i) = u(ibz,i)
      q(ibp,i) = u(ibp,i)
      ek       = 0.5d+00 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      em       = 0.5d+00 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      ei       = u(ien,i) - (ek + em)
      q(ipr,i) = gammam1 * ei

    end do ! i = 1, n

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
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the input array of conservative variables;
!     f - the output vector of fluxes;
!     c - the output vector of characteristic speeds;
!
!===============================================================================
!
  subroutine fluxspeed_mhd_adi(n, q, u, f, c)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                      , intent(in)  :: n
    real(kind=8), dimension(nv,n), intent(in)  :: q, u
    real(kind=8), dimension(nv,n), intent(out) :: f
    real(kind=8), dimension(n)   , intent(out) :: c

! local variables
!
    integer      :: i
    real(kind=8) :: bx2, by2, bz2, bb, pr, pt
    real(kind=8) :: vb
    real(kind=8) :: fa, fb, fc
!
!-------------------------------------------------------------------------------
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
      fa   = gamma * q(ipr,i)
      fb   = fa + bb
      fc   = fb * fb - 4.0d+00 * fa * bx2
      if (fc > 0.0d+00) then
        c(i) = sqrt(0.5d+00 * (fb + sqrt(fc)) / q(idn,i))
      else
        c(i) = sqrt(0.5d+00 *  fb             / q(idn,i))
      end if

! calculate the maximum speed
!
      cmax = max(cmax, abs(q(ivx,i)) + c(i))

    end do ! i = 1, n

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

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed_mhd_adi

!===============================================================================
!
end module equations
