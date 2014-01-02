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
!! module: INTERPOLATIONS
!!
!!  This module provides subroutines to interpolate variables and reconstruct
!!  the Riemann states.
!!
!!
!!******************************************************************************
!
module interpolations

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
  integer            , save :: imi, imr, imf
#endif /* PROFILE */

! pointers to the reconstruction and limiter procedures
!
  procedure(reconstruct) , pointer, save :: reconstruct_states => null()
  procedure(limiter_zero), pointer, save :: limiter => null()

! module parameters
!
  real(kind=8), save :: eps        = epsilon(1.0d+00)
  real(kind=8), save :: rad        = 0.5d+00

! flags for reconstruction corrections
!
  logical     , save :: positivity = .false.

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations, finalize_interpolations
  public :: reconstruct, limiter
  public :: fix_positivity

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_INTERPOLATIONS:
! ------------------------------------
!
!   Subroutine initializes the interpolation module by reading the module
!   parameters.
!
!
!===============================================================================
!
  subroutine initialize_interpolations(verbose, iret)

! include external procedures
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
    character(len=255) :: sreconstruction = "tvd"
    character(len=255) :: slimiter        = "mm"
    character(len=255) :: positivity_fix  = "off"
    character(len=255) :: name_rec        = ""
    character(len=255) :: name_lim        = ""
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('interpolation initialization', imi)
    call set_timer('reconstruction'              , imr)
    call set_timer('fix positivity'              , imf)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! obtain the user defined interpolation methods and coefficients
!
    call get_parameter_string("reconstruction", sreconstruction)
    call get_parameter_string("limiter"       , slimiter       )
    call get_parameter_string("fix_positivity", positivity_fix )
    call get_parameter_real  ("eps"           , eps            )
    call get_parameter_real  ("limo3_rad"     , rad            )

! select the reconstruction method
!
    select case(trim(sreconstruction))
    case ("tvd", "TVD")
      name_rec           =  "2nd order TVD"
      reconstruct_states => reconstruct_tvd
    case ("weno3", "WENO3")
      name_rec           =  "3rd order WENO"
      reconstruct_states => reconstruct_weno3
    case ("limo3", "LIMO3", "LimO3")
      name_rec           =  "3rd order logarithmic limited"
      reconstruct_states => reconstruct_limo3
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected reconstruction method is not " //     &
                           "implemented: " // trim(sreconstruction)
        stop
      end if
    end select

! select the limiter
!
    select case(trim(slimiter))
    case ("mm", "minmod")
      name_lim           =  "minmod"
      limiter            => limiter_minmod
    case ("mc", "monotonized_central")
      name_lim           =  "monotonized central"
      limiter            => limiter_monotonized_central
    case ("sb", "superbee")
      name_lim           =  "superbee"
      limiter            => limiter_superbee
    case ("vl", "vanleer")
      name_lim           =  "van Leer"
      limiter            => limiter_vanleer
    case ("va", "vanalbada")
      name_lim           =  "van Albada"
      limiter            => limiter_vanalbada
    case default
      name_lim           =  "zero derivative"
      limiter            => limiter_zero
    end select

! check additional reconstruction limiting
!
    select case(trim(positivity_fix))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      positivity = .true.
    case default
      positivity = .false.
    end select

! print informations about the reconstruction methods and parameters
!
    if (verbose) then

      write (*,"(4x,a14,9x,'=',1x,a)") "reconstruction", trim(name_rec)
      write (*,"(4x,a14,9x,'=',1x,a)") "limiter       ", trim(name_lim)
      write (*,"(4x,a14,9x,'=',1x,a)") "fix positivity", trim(positivity_fix)

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_interpolations
!
!===============================================================================
!
! subroutine FINALIZE_INTERPOLATIONS:
! ----------------------------------
!
!   Subroutine finalizes the interpolation module by releasing all memory used
!   by its module variables.
!
!
!===============================================================================
!
  subroutine finalize_interpolations(iret)

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

! release the procedure pointers
!
    nullify(reconstruct_states)
    nullify(limiter)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_interpolations
!
!===============================================================================
!
! subroutine RECONSTRUCT:
! ----------------------
!
!   Subroutine calls a reconstruction procedure, depending on the compilation
!   flag SPACE, in order to interpolate the left and right states from their
!   cell integrals.  These states are required by any approximate Riemann
!   solver.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     h  - the spatial step; this is required for some reconstruction methods;
!     f  - the input vector of cell averaged values;
!     fl - the left side state reconstructed for location (i+1/2);
!     fr - the right side state reconstructed for location (i+1/2);
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for reconstruction
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the states using the selected subroutine
!
    call reconstruct_states(n, h, f(:), fl(:), fr(:))

#ifdef PROFILE
! stop accounting time for reconstruction
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
!===============================================================================
!
! subroutine RECONSTRUCT_TVD:
! --------------------------
!
!   Subroutine reconstructs the interface states using the second order TVD
!   method with a selected limiter.
!
!   Arguments are described in subroutine reconstruct().
!
!
!===============================================================================
!
  subroutine reconstruct_tvd(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      ::  i, im1, ip1
    real(kind=8) :: df, dfl, dfr
!
!-------------------------------------------------------------------------------
!
! calculate the left- and right-side interface interpolations
!
    do i = 1, n

! calculate left and right indices
!
      im1     = max(1, i - 1)
      ip1     = min(n, i + 1)

! calculate left and right side derivatives
!
      dfl     = f(i  ) - f(im1)
      dfr     = f(ip1) - f(i  )

! obtain the TVD limited derivative
!
      df      = limiter(0.5d+00, dfl, dfr)

! update the left and right-side interpolation states
!
      fl(i  ) = f(i) + df
      fr(im1) = f(i) - df

    end do ! i = 1, n

! update the interpolation of the first and last points
!
    fl(1) = f(1)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_tvd
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO3:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the third order
!   Weighted Essentially Non-Oscillatory (WENO) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Yamaleev & Carpenter, 2009, J. Comput. Phys., 228, 3025
!
!===============================================================================
!
  subroutine reconstruct_weno3(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: bp, bm, ap, am, wp, wm, ww
    real(kind=8) :: dfl, dfr, df, fp, fm, fc, h2

! selection weights
!
    real, parameter :: dp = 2.0d+00 / 3.0d+00, dm = 1.0d+00 / 3.0d+00
!
!-------------------------------------------------------------------------------
!
! prepare common parameters
!
    h2 = h * h

! iterate along the vector
!
    do i = 1, n

! prepare neighbour indices
!
      im1     = max(1, i - 1)
      ip1     = min(n, i + 1)

! calculate the left and right derivatives
!
      dfl     = f(i  ) - f(im1)
      dfr     = f(ip1) - f(i  )

! calculate coefficient omega
!
      ww      = (dfr - dfl)**2

! calculate corresponding betas
!
      bp      = dfr * dfr
      bm      = dfl * dfl

! calculate improved alphas
!
      ap      = 1.0d+00 + ww / (bp + h2)
      am      = 1.0d+00 + ww / (bm + h2)

! calculate weights
!
      wp      = dp * ap
      wm      = dm * am
      ww      = 2.0d+00 * (wp + wm)

! calculate central interpolation
!
      fp      =   f(i  ) +         f(ip1)

! calculate left side interpolation
!
      fm      = - f(im1) + 3.0d+00 * f(i  )

! calculate the left state
!
      fl( i ) = (wp * fp + wm * fm) / ww

! calculate weights
!
      wp      = dp * am
      wm      = dm * ap
      ww      = 2.0d+00 * (wp + wm)

! calculate central interpolation
!
      fp      =   f(i  ) +         f(im1)

! calculate right side interpolation
!
      fm      = - f(ip1) + 3.0d+00 * f(i  )

! calculate the right state
!
      fr(im1) = (wp * fp + wm * fm) / ww

    end do ! i = 1, n

! update the interpolation of the first and last points
!
    fl(1) = f (1)
    fr(n) = fl(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno3
!
!===============================================================================
!
! subroutine RECONSTRUCT_LIMO3:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the third order method
!   with a limiter function LimO3.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Cada, M. & Torrilhon, M.,
!         "Compact third-order limiter functions for finite volume methods",
!         Journal of Computational Physics, 2009, 228, 4118-4145
!     [2] Mignone, A., Tzeferacos, P., & Bodo, G.,
!         "High-order conservative finite divergence GLM-MHD schemes for
!          cell-centered MHD",
!         Journal of Computational Physics, 2010, 229, 5896-5920
!
!===============================================================================
!
  subroutine reconstruct_limo3(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: dfl, dfr
    real(kind=8) :: th, et, f1, f2, xl, xi, rdx, rdx2
!
!-------------------------------------------------------------------------------
!
! prepare parameters
!
    rdx   = rad * h
    rdx2  = rdx * rdx

! iterate over positions and interpolate states
!
    do i = 1, n

! prepare neighbour indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! prepare left and right differences
!
      dfl = f(i  ) - f(im1)
      dfr = f(ip1) - f(i  )

! calculate the indicator function (eq. 3.17 in [1])
!
      et = (dfl * dfl + dfr * dfr) / rdx2

! the switching function (embedded in eq. 3.22 in [1], eq. 32 in [2])
!
      xi = max(0.0d+00, 0.5d+00 * min(2.0d+00, 1.0d+00 + (et - 1.0d+00) / eps))
      xl = 1.0d+00 - xi

! calculate values at i + ½
!
      if (dfr == 0.0d+00) then

        fl(i) = f(i)

      else

! calculate the slope ratio (eq. 2.8 in [1])
!
        th = dfl / dfr

! calculate the quadratic reconstruction (eq. 3.8 in [1], divided by 2)
!
        f1 = (2.0d+00 + th) / 6.0d+00

! calculate the third order limiter (eq. 3.13 in [1], cofficients divided by 2)
!
        if (th >= 0.0d+00) then
          f2 = max(0.0d+00, min(f1, th, 0.8d+00))
        else
          f2 = max(0.0d+00, min(f1, - 0.25d+00 * th))
        end if

! interpolate the left state (eq. 3.5 in [1], eq. 30 in [2])
!
        fl(i) = f(i) + dfr * (xl * f1 + xi * f2)

      end if

! calculate values at i - ½
!
      if (dfl == 0.0d+00) then

        fr(im1) = f(i)

      else

! calculate the slope ratio (eq. 2.8 in [1])
!
        th = dfr / dfl

! calculate the quadratic reconstruction (eq. 3.8 in [1], divided by 2)
!
        f1 = (2.0d+00 + th) / 6.0d+00

! calculate the third order limiter (eq. 3.13 in [1], cofficients divided by 2)
!
        if (th >= 0.0d+00) then
          f2 = max(0.0d+00, min(f1, th, 0.8d+00))
        else
          f2 = max(0.0d+00, min(f1, - 0.25d+00 * th))
        end if

! interpolate the right state (eq. 3.5 in [1], eq. 30 in [2])
!
        fr(im1) = f(i) - dfl * (xl * f1 + xi * f2)

      end if

    end do ! i = 1, n

! update the interpolation of the first and last points
!
    fl(1) = f (1)
    fr(n) = fl(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_limo3
!
!===============================================================================
!
! function LIMITER_ZERO:
! ---------------------
!
!   Function returns zero.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_zero(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.0d+00

!-------------------------------------------------------------------------------
!
  end function limiter_zero
!
!===============================================================================
!
! function LIMITER_MINMOD:
! -----------------------
!
!   Function returns the minimum module value among two arguments using
!   minmod limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_minmod(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.5d+00 * (sign(x, a) + sign(x, b)) * min(abs(a), abs(b))

!-------------------------------------------------------------------------------
!
  end function limiter_minmod
!
!===============================================================================
!
! function LIMITER_MONOTONIZED_CENTRAL:
! ------------------------------------
!
!   Function returns the minimum module value among two arguments using
!   the monotonized central TVD limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_monotonized_central(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = (sign(x, a) + sign(x, b)) * min(abs(a), abs(b), 2.5d-01 * abs(a + b))

!-------------------------------------------------------------------------------
!
  end function limiter_monotonized_central
!
!===============================================================================
!
! function LIMITER_SUPERBEE:
! -------------------------
!
!   Function returns the minimum module value among two arguments using
!   superbee limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_superbee(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.5d+00 * (sign(x, a) + sign(x, b))                                    &
           * max(min(2.0d+00 * abs(a), abs(b)), min(abs(a), 2.0d+00 * abs(b)))

!-------------------------------------------------------------------------------
!
  end function limiter_superbee
!
!===============================================================================
!
! function LIMITER_VANLEER:
! ------------------------
!
!   Function returns the minimum module value among two arguments using
!   van Leer's limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_vanleer(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = a * b
    if (c > 0.0d+00) then
      c = 2.0d+00 * x * c / (a + b)
    else
      c = 0.0d+00
    end if

!-------------------------------------------------------------------------------
!
  end function limiter_vanleer
!
!===============================================================================
!
! function LIMITER_VANALBADA:
! --------------------------
!
!   Function returns the minimum module value among two arguments using
!   van Albada's limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_vanalbada(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = x * a * b * (a + b) / max(eps, a * a + b * b)

!-------------------------------------------------------------------------------
!
  end function limiter_vanalbada
!
!===============================================================================
!
! subroutine FIX_POSITIVITY:
! -------------------------
!
!   Subroutine scans the input arrays of the left and right states fl(:) and
!   fr(:) for negative values.  If it finds a negative value, it repeates the
!   state reconstruction from f(:) using the zeroth order interpolation.
!
!===============================================================================
!
  subroutine fix_positivity(n, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(in)    :: f
    real(kind=8), dimension(n), intent(inout) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: fmn, fmx
!
!------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for positivity fix
!
    call start_timer(imf)
#endif /* PROFILE */

! check positivity only if desired
!
    if (positivity) then

! look for negative values in the states along the vector
!
      do i = 1, n

! check if the left state has a negative value
!
        if (fl(i) <= 0.0d+00) then

! calculate the left neighbour index
!
          im1 = max(1, i - 1)

! limit the states using the zeroth-order reconstruction
!
          fl(i  ) = f(i)
          fr(im1) = f(i)

        end if ! fl ≤ 0

! check if the right state has a negative value
!
        if (fr(i) <= 0.0d+00) then

! calculate the right neighbour index
!
          ip1 = min(n, i + 1)

! limit the states using the zeroth-order reconstruction
!
          fl(ip1) = f(ip1)
          fr(i  ) = f(ip1)

        end if ! fr ≤ 0

      end do ! i = 1, n

    end if ! positivity == .true.

#ifdef PROFILE
! stop accounting time for positivity fix
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fix_positivity

!===============================================================================
!
end module interpolations
