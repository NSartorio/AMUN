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
!! module: EQUATIONS
!!
!!  This module handles supported sets of equations, provides subroutines to
!!  convert between conserved and primitive variables, calculate flux and
!!  characteristic speeds.
!!
!!******************************************************************************
!
module equations

! module variables are not implicit by default
!
  implicit none

#ifdef ADI
! adiabatic heat ratio
!
  real, save :: gamma = 5.0d0 / 3.0d0

! additional adiabatic parameters
!
  real, save :: gammam1 = 2.0d0 / 3.0d0, gammam1i = 1.5d0
#endif /* ADI */

#ifdef ISO
! isothermal speed of sound and its second power
!
  real, save :: csnd  = 1.0d0, csnd2 = 1.0d0
#endif /* ISO */

! by default everything is public
!
  public

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_EQUATIONS:
! -------------------------------
!
!   Subroutine sets the default values of module parameters and obtains their
!   values from the PARAMETERS module.
!
!===============================================================================
!
  subroutine initialize_equations()

! include external procedures and variables
!
    use parameters, only : get_parameter_real

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
#ifdef ADI
! obtain the adiabatic specific heat ratio
!
    call get_parameter_real("gamma" , gamma )

! calculate additional parameters
!
    gammam1  = gamma - 1.0d0
    gammam1i = 1.0d0 / gammam1
#endif /* ADI */

#ifdef ISO
! obtain the isothermal sound speed
!
    call get_parameter_real("csnd"  , csnd  )

! calculate additional parameters
!
    csnd2 = csnd * csnd
#endif /* ISO */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_equations
#ifdef HYDRO
!
!===============================================================================
!
! subroutine CONS2PRIM:
! --------------------
!
!   Subroutine converts the conservative representation of variables to
!   its corresponding primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim(n, u, q)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: u
    real, dimension(nt,n), intent(out) :: q

! local variables
!
    integer :: i
#ifdef ADI
    real    :: ek, ei
#endif /* ADI */
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

      q(idn,i) = u(idn,i)
      q(ivx,i) = u(imx,i) / u(idn,i)
      q(ivy,i) = u(imy,i) / u(idn,i)
      q(ivz,i) = u(imz,i) / u(idn,i)
#ifdef ADI
      ek       = 0.5d0 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      ei       = u(ien,i) - ek
      q(ipr,i) = gammam1 * ei
#endif /* ADI */

    end do

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim
!
!===============================================================================
!
! subroutine PRIM2CONS:
! --------------------
!
!   Subroutine converts the primitive variable representation to its
!   corresponding conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons(n, q, u)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: u

! local variables
!
    integer :: i
#ifdef ADI
    real    :: ek, ei
#endif /* ADI */
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

      u(idn,i) = q(idn,i)
      u(imx,i) = q(idn,i) * q(ivx,i)
      u(imy,i) = q(idn,i) * q(ivy,i)
      u(imz,i) = q(idn,i) * q(ivz,i)
#ifdef ADI
      ek       = 0.5d0 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      ei       = gammam1i * q(ipr,i)
      u(ien,i) = ei + ek
#endif /* ADI */

    end do

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons
!
!===============================================================================
!
! subroutine FLUXSPEED:
! --------------------
!
!   Subroutine calculates fluxes and characteristic speeds from primitive and
!   conservative variables.
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
  subroutine fluxspeed(n, q, u, f, c)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: q, u
    real, dimension(nt,n), intent(out) :: f
    real, dimension(n)   , intent(out) :: c

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
#ifdef ISO
      f(imx,i) = f(imx,i) + csnd2 * q(idn,i)
#endif /* ISO */
#ifdef ADI
      f(imx,i) = f(imx,i) + q(ipr,i)
      f(ien,i) = q(ivx,i) * (u(ien,i) + q(ipr,i))
#endif /* ADI */

! calculate the speed of sound
!
#ifdef ADI
      c(i) = sqrt(gamma * q(ipr,i) / q(idn,i))
#endif /* ADI */
#ifdef ISO
      c(i) = csnd
#endif /* ISO */

    end do

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed
!
!===============================================================================
!
! function MAXSPEED:
! -----------------
!
!   Function scans the variable array and returns the maximum speed in the
!   system.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!
!===============================================================================
!
  function maxspeed(q)

! include external procedures and variables
!
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke
    use variables  , only : nt
    use variables  , only : ivx, ivz
#ifdef ADI
    use variables  , only : idn, ipr
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real, dimension(nt,im,jm,km), intent(in) :: q

! return variable
!
    real     :: maxspeed

! local variables
!
    integer  :: i, j, k
    real     :: vv, v, c
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed
!
    maxspeed = 0.0d0

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(q(ivx:ivz,i,j,k) * q(ivx:ivz,i,j,k))
          v  = sqrt(vv)

#ifdef ADI
! calculate the adiabatic speed of sound
!
          c = sqrt(gamma * q(ipr,i,j,k) / q(idn,i,j,k))
#endif /* ADI */

! calculate the maximum speed
!
#ifdef ISO
          maxspeed = max(maxspeed, v + csnd)
#endif /* ISO */
#ifdef ADI
          maxspeed = max(maxspeed, v + c)
#endif /* ADI */

        end do
      end do
    end do

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed
!
!===============================================================================
!
! subroutine UPDATE_PRIMITIVE_VARIABLES:
! -------------------------------------
!
!   Subroutine updates primitive variables from their conservative
!   representation.  This process is done ones after advance of the conserved
!   variables due to their evolution in time.
!
!===============================================================================
!
  subroutine update_primitive_variables(uu, qq)

! include external procedures and variables
!
    use coordinates, only : im, jm, km
#ifdef DEBUG
    use variables  , only : idn
#ifdef ADI
    use variables  , only : ipr
#endif /* ADI */
#endif /* DEBUG */
    use variables  , only : nt

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real, dimension(nt,im,jm,km), intent(in)    :: uu
    real, dimension(nt,im,jm,km), intent(inout) :: qq

! temporary variables
!
#ifdef DEBUG
    integer                :: i
#endif /* DEBUG */
    integer                :: j, k

! temporary array to store conserved variable vector
!
    real, dimension(nt,im) :: u
!
!-------------------------------------------------------------------------------
!
! update primitive variables
!
    do k = 1, km
      do j = 1, jm

! copy variables to temporary array of conserved variables
!
        u(1:nt,1:im) = uu(1:nt,1:im,j,k)

! convert conserved variables to primitive ones
!
        call cons2prim(im, u(1:nt,1:im), qq(1:nt,1:im,j,k))

#ifdef DEBUG
! check the positivity of density and pressure
!
        do i = 1, im
          if (qq(idn,i,j,k) <= 0.0d0) then
            write(*,*)
            write(*,*) 'Unphysical state: ρ ≤ 0', i, j, k
            write(*,"('Q = ',4(1pe24.16))") qq(1:4 ,i,j,k)
            write(*,"('    ',4(1pe24.16))") qq(5:nt,i,j,k)
            stop
          end if
#ifdef ADI
          if (qq(ipr,i,j,k) <= 0.0d0) then
            write(*,*)
            write(*,*) 'Unphysical state: p ≤ 0', i, j, k
            write(*,"('Q = ',4(1pe24.16))") qq(1:4 ,i,j,k)
            write(*,"('    ',4(1pe24.16))") qq(5:nt,i,j,k)
            stop
          end if
#endif /* ADI */
        end do
#endif /* DEBUG */

      end do
    end do

!-------------------------------------------------------------------------------
!
  end subroutine update_primitive_variables
#endif /* HYDRO */
#ifdef MHD
!
!===============================================================================
!
! subroutine CONS2PRIM:
! --------------------
!
!   Subroutine converts the conservative representation of variables to
!   its corresponding primitive representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     u - the input array of conservative variables;
!     q - the output array of primitive variables;
!
!===============================================================================
!
  subroutine cons2prim(n, u, q)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz, ibx, iby, ibz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: u
    real, dimension(nt,n), intent(out) :: q

! local variables
!
    integer :: i
#ifdef ADI
    real    :: ei, ek, em
#endif /* ADI */
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
#ifdef ADI
      ek       = 0.5d0 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      em       = 0.5d0 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      ei       = u(ien,i) - (ek + em)
      q(ipr,i) = gammam1 * ei
#endif /* ADI */

    end do

!-------------------------------------------------------------------------------
!
  end subroutine cons2prim
!
!===============================================================================
!
! subroutine PRIM2CONS:
! --------------------
!
!   Subroutine converts the primitive variable representation to its
!   corresponding conservative representation.
!
!   Arguments:
!
!     n - the length of input and output vectors;
!     q - the input array of primitive variables;
!     u - the output array of conservative variables;
!
!===============================================================================
!
  subroutine prim2cons(n, q, u)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz, ibx, iby, ibz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: q
    real, dimension(nt,n), intent(out) :: u

! local variables
!
    integer :: i
#ifdef ADI
    real    :: ei, ek, em
#endif /* ADI */
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
#ifdef ADI
      ei       = gammam1i * q(ipr,i)
      ek       = 0.5d0 * sum(u(imx:imz,i) * q(ivx:ivz,i))
      em       = 0.5d0 * sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      u(ien,i) = ei + ek + em
#endif /* ADI */

    end do

!-------------------------------------------------------------------------------
!
  end subroutine prim2cons
!
!===============================================================================
!
! subroutine FLUXSPEED:
! --------------------
!
!   Subroutine calculates fluxes and characteristic speeds from primitive and
!   conservative variables.
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
  subroutine fluxspeed(n, q, u, f, c)

! include external procedures and variables
!
    use variables , only : nt
    use variables , only : idn, imx, imy, imz, ivx, ivy, ivz
    use variables , only : ibx, iby, ibz
#ifdef ADI
    use variables , only : ipr, ien
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer              , intent(in)  :: n
    real, dimension(nt,n), intent(in)  :: q, u
    real, dimension(nt,n), intent(out) :: f
    real, dimension(n)   , intent(out) :: c

! local variables
!
    integer :: i
    real    :: bb, vb, pm, pt
    real    :: cs2, ca2, cx2, cf2, cl2
!
!-------------------------------------------------------------------------------
!
    do i = 1, n

! prepare pressures and scalar product
!
      bb = sum(q(ibx:ibz,i) * q(ibx:ibz,i))
      pm = 0.5d0 * bb
#ifdef ADI
      vb = sum(q(ivx:ivz,i) * q(ibx:ibz,i))
      pt = q(ipr,i)
#endif /* ADI */
#ifdef ISO
      pt = csnd2 * q(idn,i)
#endif /* ISO */
      pt = pt + pm

! calculate the magnetohydrodynamic fluxes
!
      f(idn,i) = u(imx,i)
      f(imx,i) = q(ivx,i) * u(imx,i) - q(ibx,i) * q(ibx,i)
      f(imy,i) = q(ivx,i) * u(imy,i) - q(ibx,i) * q(iby,i)
      f(imz,i) = q(ivx,i) * u(imz,i) - q(ibx,i) * q(ibz,i)
      f(imx,i) = f(imx,i) + pt
#ifdef ADI
      f(ien,i) = q(ivx,i) * (u(ien,i) + pt) - q(ibx,i) * vb
#endif /* ADI */
      f(ibx,i) = 0.0d0
      f(iby,i) = q(ivx,i) * q(iby,i) - q(ibx,i) * q(ivy,i)
      f(ibz,i) = q(ivx,i) * q(ibz,i) - q(ibx,i) * q(ivz,i)

! calculate the fast magnetosonic speed
!
#ifdef ADI
      cs2 = gamma * q(ipr,i) / q(idn,i)
#endif /* ADI */
#ifdef ISO
      cs2 = csnd2
#endif /* ISO */
      ca2 = sum(q(ibx:ibz,i) * q(ibx:ibz,i)) / q(idn,i)
      cx2 = q(ibx,i) * q(ibx,i) / q(idn,i)
      cf2 = cs2 + ca2
      cl2 = max(0.0d0, cf2**2 - 4.0d0 * cs2 * cx2)

      c(i) = sqrt(0.5d0 * (cf2 + sqrt(cl2)))

    end do

!-------------------------------------------------------------------------------
!
  end subroutine fluxspeed
!
!===============================================================================
!
! function MAXSPEED:
! -----------------
!
!   Function scans the variable array and returns the maximum speed in the
!   system.
!
!   Arguments:
!
!     q - the array of primitive variables;
!
!
!===============================================================================
!
  function maxspeed()

! include external procedures and variables
!
    use mesh      , only : im, jm, km, ib, ie, jb, je, kb, ke
    use variables , only : nt
    use variables , only : idn, ivx, ivz, ibx, ibz
#ifdef ADI
    use variables , only : ipr
#endif /* ADI */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real, dimension(nt,im,jm,km), intent(in) :: q

! return variable
!
    real     :: maxspeed

! local variables
!
    integer  :: i, j, k
    real     :: vv, bb, v, c
!
!-------------------------------------------------------------------------------
!
! reset the maximum speed
!
    maxspeed = 0.0d0

! iterate over all positions
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! calculate the velocity amplitude
!
          vv = sum(q(ivx:ivz,i,j,k) * q(ivx:ivz,i,j,k))
          v  = sqrt(vv)
          bb = sum(q(ibx:ibz,i,j,k) * q(ibx:ibz,i,j,k))

! calculate the fast magnetosonic speed
!
#ifdef ISO
          c = sqrt(csnd2 + bb / q(idn,i,j,k))
#endif /* ISO */
#ifdef ADI
          c = sqrt((gamma * q(ipr,i,j,k) + bb) / q(idn,i,j,k))
#endif /* ADI */

! calculate the maximum of speed
!
          maxspeed = max(maxspeed, v + c)

        end do
      end do
    end do

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function maxspeed
#endif /* MHD */

!===============================================================================
!
end module equations
