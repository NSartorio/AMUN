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
  integer            , save :: imi, imr, imf, imc
#endif /* PROFILE */

! pointers to the reconstruction and limiter procedures
!
  procedure(reconstruct)       , pointer, save :: reconstruct_states => null()
  procedure(limiter_zero)      , pointer, save :: limiter_tvd        => null()
  procedure(limiter_zero)      , pointer, save :: limiter_prol       => null()
  procedure(limiter_zero)      , pointer, save :: limiter_clip       => null()

! module parameters
!
  real(kind=8), save :: eps        = epsilon(1.0d+00)
  real(kind=8), save :: rad        = 0.5d+00

! monotonicity preserving reconstruction coefficients
!
  real(kind=8), save :: kappa      = 1.0d+00
  real(kind=8), save :: kbeta      = 1.0d+00

! number of ghost zones (required for compact schemes)
!
  integer     , save :: ng         = 2

! flags for reconstruction corrections
!
  logical     , save :: positivity = .false.
  logical     , save :: clip       = .false.

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations, finalize_interpolations
  public :: reconstruct, limiter_prol
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
    use error     , only : print_warning
    use parameters, only : get_parameter_string, get_parameter_integer         &
                         , get_parameter_real

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
    character(len=255) :: tlimiter        = "mm"
    character(len=255) :: plimiter        = "mm"
    character(len=255) :: climiter        = "mm"
    character(len=255) :: positivity_fix  = "off"
    character(len=255) :: clip_extrema    = "off"
    character(len=255) :: name_rec        = ""
    character(len=255) :: name_tlim       = ""
    character(len=255) :: name_plim       = ""
    character(len=255) :: name_clim       = ""
    real(kind=8)       :: cfl             = 0.5d+00
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('interpolations:: initialization', imi)
    call set_timer('interpolations:: reconstruction', imr)
    call set_timer('interpolations:: fix positivity', imf)
    call set_timer('interpolations:: clip extrema'  , imc)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! obtain the user defined interpolation methods and coefficients
!
    call get_parameter_string ("reconstruction"      , sreconstruction)
    call get_parameter_string ("limiter"             , tlimiter       )
    call get_parameter_string ("fix_positivity"      , positivity_fix )
    call get_parameter_string ("clip_extrema"        , clip_extrema   )
    call get_parameter_string ("extrema_limiter"     , climiter       )
    call get_parameter_string ("prolongation_limiter", plimiter       )
    call get_parameter_integer("nghosts"             , ng             )
    call get_parameter_real   ("eps"                 , eps            )
    call get_parameter_real   ("limo3_rad"           , rad            )
    call get_parameter_real   ("kappa"               , kappa          )
    call get_parameter_real   ("kbeta"               , kbeta          )
    call get_parameter_real   ("cfl"                 , cfl            )

! calculate κ = (1 - ν) / ν
!
    kappa = min(kappa, (1.0d+00 - cfl) / cfl)

! select the reconstruction method
!
    select case(trim(sreconstruction))
    case ("tvd", "TVD")
      name_rec           =  "2nd order TVD"
      reconstruct_states => reconstruct_tvd
      if (verbose .and. ng < 2)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 2).")
    case ("weno3", "WENO3")
      name_rec           =  "3rd order WENO"
      reconstruct_states => reconstruct_weno3
      if (verbose .and. ng < 2)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 2).")
    case ("limo3", "LIMO3", "LimO3")
      name_rec           =  "3rd order logarithmic limited"
      reconstruct_states => reconstruct_limo3
      if (verbose .and. ng < 2)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 2).")
      eps = max(1.0d-12, eps)
    case ("weno5z", "weno5-z", "WENO5Z", "WENO5-Z")
      name_rec           =  "5th order WENO-Z (Borges et al. 2008)"
      reconstruct_states => reconstruct_weno5z
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("weno5yc", "weno5-yc", "WENO5YC", "WENO5-YC")
      name_rec           =  "5th order WENO-YC (Yamaleev & Carpenter 2009)"
      reconstruct_states => reconstruct_weno5yc
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("weno5ns", "weno5-ns", "WENO5NS", "WENO5-NS")
      name_rec           =  "5th order WENO-NS (Ha et al. 2013)"
      reconstruct_states => reconstruct_weno5ns
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("crweno5z", "crweno5-z", "CRWENO5Z", "CRWENO5-Z")
      name_rec           =  "5th order Compact WENO-Z"
      reconstruct_states => reconstruct_crweno5z
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("crweno5yc", "crweno5-yc", "CRWENO5YC", "CRWENO5-YC")
      name_rec           =  "5th order Compact WENO-YC"
      reconstruct_states => reconstruct_crweno5yc
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("crweno5ns", "crweno5-ns", "CRWENO5NS", "CRWENO5-NS")
      name_rec           =  "5th order Compact WENO-NS"
      reconstruct_states => reconstruct_crweno5ns
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("mp5", "MP5")
      name_rec           =  "5th order Monotonicity Preserving"
      reconstruct_states => reconstruct_mp5
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case ("crmp5", "CRMP5")
      name_rec           =  "5th order Compact Monotonicity Preserving"
      reconstruct_states => reconstruct_crmp5
      if (verbose .and. ng < 4)                                                &
                  call print_warning("interpolations:initialize_interpolation" &
                         , "Increase the number of ghost cells (at least 4).")
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected reconstruction method is not " //     &
                           "implemented: " // trim(sreconstruction)
        stop
      end if
    end select

! select the TVD limiter
!
    select case(trim(tlimiter))
    case ("mm", "minmod")
      name_tlim          =  "minmod"
      limiter_tvd        => limiter_minmod
    case ("mc", "monotonized_central")
      name_tlim          =  "monotonized central"
      limiter_tvd        => limiter_monotonized_central
    case ("sb", "superbee")
      name_tlim          =  "superbee"
      limiter_tvd        => limiter_superbee
    case ("vl", "vanleer")
      name_tlim          =  "van Leer"
      limiter_tvd        => limiter_vanleer
    case ("va", "vanalbada")
      name_tlim          =  "van Albada"
      limiter_tvd        => limiter_vanalbada
    case default
      name_tlim          =  "zero derivative"
      limiter_tvd        => limiter_zero
    end select

! select the prolongation limiter
!
    select case(trim(plimiter))
    case ("mm", "minmod")
      name_plim          =  "minmod"
      limiter_prol       => limiter_minmod
    case ("mc", "monotonized_central")
      name_plim          =  "monotonized central"
      limiter_prol       => limiter_monotonized_central
    case ("sb", "superbee")
      name_plim          =  "superbee"
      limiter_prol       => limiter_superbee
    case ("vl", "vanleer")
      name_plim          =  "van Leer"
      limiter_prol       => limiter_vanleer
    case default
      name_plim          =  "zero derivative"
      limiter_prol       => limiter_zero
    end select

! select the clipping limiter
!
    select case(trim(climiter))
    case ("mm", "minmod")
      name_clim          =  "minmod"
      limiter_clip       => limiter_minmod
    case ("mc", "monotonized_central")
      name_clim          =  "monotonized central"
      limiter_clip       => limiter_monotonized_central
    case ("sb", "superbee")
      name_clim          =  "superbee"
      limiter_clip       => limiter_superbee
    case ("vl", "vanleer")
      name_clim          =  "van Leer"
      limiter_clip       => limiter_vanleer
    case default
      name_clim          =  "zero derivative"
      limiter_clip       => limiter_zero
    end select

! check additional reconstruction limiting
!
    select case(trim(positivity_fix))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      positivity = .true.
    case default
      positivity = .false.
    end select
    select case(trim(clip_extrema))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      clip = .true.
    case default
      clip = .false.
    end select

! print informations about the reconstruction methods and parameters
!
    if (verbose) then

      write (*,"(4x,a14, 9x,'=',1x,a)") "reconstruction"      , trim(name_rec)
      write (*,"(4x,a11,12x,'=',1x,a)") "TVD limiter"         , trim(name_tlim)
      write (*,"(4x,a20, 3x,'=',1x,a)") "prolongation limiter", trim(name_plim)
      write (*,"(4x,a14, 9x,'=',1x,a)") "fix positivity"      , trim(positivity_fix)
      write (*,"(4x,a12,11x,'=',1x,a)") "clip extrema"        , trim(clip_extrema)
      if (clip) then
        write (*,"(4x,a15,8x,'=',1x,a)") "extrema limiter", trim(name_clim)
      end if

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
    nullify(limiter_tvd)
    nullify(limiter_prol)
    nullify(limiter_clip)

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

! correct the reconstruction near extrema by clipping them in order to improve
! the stability of scheme
!
    if (clip) call clip_extrema(n, f(:), fl(:), fr(:))

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
    do i = 2, n - 1

! calculate left and right indices
!
      im1     = i - 1
      ip1     = i + 1

! calculate left and right side derivatives
!
      dfl     = f(i  ) - f(im1)
      dfr     = f(ip1) - f(i  )

! obtain the TVD limited derivative
!
      df      = limiter_tvd(0.5d+00, dfl, dfr)

! update the left and right-side interpolation states
!
      fl(i  ) = f(i) + df
      fr(im1) = f(i) - df

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
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
    real(kind=8), parameter :: dp = 2.0d+00 / 3.0d+00, dm = 1.0d+00 / 3.0d+00
!
!-------------------------------------------------------------------------------
!
! prepare common parameters
!
    h2 = h * h

! iterate along the vector
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1     = i - 1
      ip1     = i + 1

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

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

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
    do i = 2, n - 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

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
      if (abs(dfr) > eps) then

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

      else

        fl(i) = f(i)

      end if

! calculate values at i - ½
!
      if (abs(dfl) > eps) then

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

      else

        fr(im1) = f(i)

      end if

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_limo3
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5Z:
! -----------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with
!   stencil weights by Borges et al. (2008).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Borges, R., Carmona, M., Costa, B., & Don, W.-S.,
!         "An improved weighted essentially non-oscillatory scheme for
!          hyperbolic conservation laws"
!         Journal of Computational Physics,
!         2008, vol. 227, pp. 3191-3211,
!         http://dx.doi.org/10.1016/j.jcp.2007.11.038
!
!===============================================================================
!
  subroutine reconstruct_weno5z(n, h, f, fl, fr)

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: al, ac, ar
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! iterate along the vector
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate βₖ (eqs. 9-11 in [1])
!
      bl  = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc  = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br  = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 25 in [1])
!
      tt  = abs(br - bl)

! calculate αₖ (eq. 28 in [1])
!
      al  = 1.0d+00 + tt / (bl + eps)
      ac  = 1.0d+00 + tt / (bc + eps)
      ar  = 1.0d+00 + tt / (br + eps)

! calculate weights
!
      wl  = dl * al
      wc  = dc * ac
      wr  = dr * ar
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! normalize weights
!
      wl  = dl * ar
      wc  = dc * ac
      wr  = dr * al
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5z
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5YC:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with
!   stencil weights by Yamaleev & Carpenter (2009).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Yamaleev, N. K. & Carpenter, H. C.,
!         "A Systematic Methodology for Constructing High-Order Energy Stable
!          WENO Schemes"
!         Journal of Computational Physics,
!         2009, vol. 228, pp. 4248-4272,
!         http://dx.doi.org/10.1016/j.jcp.2009.03.002
!
!===============================================================================
!
  subroutine reconstruct_weno5yc(n, h, f, fl, fr)

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: al, ac, ar
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! iterate along the vector
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate βₖ (eq. 19 in [1])
!
      bl  = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc  = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br  = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 20 in [1])
!
      tt  = (6.0d+00 * f(i) - 4.0d+00 * (f(im1) + f(ip1))                      &
                                                       + (f(im2) + f(ip2)))**2

! calculate αₖ (eqs. 18 or 58 in [1])
!
      al  = 1.0d+00 + tt / (bl + eps)
      ac  = 1.0d+00 + tt / (bc + eps)
      ar  = 1.0d+00 + tt / (br + eps)

! calculate weights
!
      wl  = dl * al
      wc  = dc * ac
      wr  = dr * ar
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! normalize weights
!
      wl  = dl * ar
      wc  = dc * ac
      wr  = dr * al
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5yc
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5NS:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with new
!   smoothness indicators and stencil weights by Ha et al. (2013).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ha, Y., Kim, C. H., Lee, Y. J., & Yoon, J.,
!         "An improved weighted essentially non-oscillatory scheme with a new
!          smoothness indicator",
!         Journal of Computational Physics,
!         2013, vol. 232, pp. 68-86
!         http://dx.doi.org/10.1016/j.jcp.2012.06.016
!
!===============================================================================
!
  subroutine reconstruct_weno5ns(n, h, f, fl, fr)

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br
    real(kind=8) :: al, ac, ar, aa
    real(kind=8) :: wl, wc, wr
    real(kind=8) :: df, lq, l3, zt
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00

! the free parameter for smoothness indicators (see Eq. 3.6 in [1])
!
    real(kind=8), parameter :: xi  =   4.0d-01
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = 0.5d+00 * abs(dfp(:) - dfm(:))

! iterate along the vector
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate βₖ (eq. 3.6 in [1])
!
      df  = abs(dfp(i))
      lq  = xi * df
      bl  = df2(im1) + xi * abs(2.0d+00 * dfm(i) - dfm(im1))
      bc  = df2(i  ) + lq
      br  = df2(ip1) + lq

! calculate ζ (below eq. 3.6 in [1])
!
      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ (eq. 3.9 in [4])
!
      al  = dl * (1.0d+00 + zt / (bl + eps)**2)
      ac  = dc * (1.0d+00 + zt / (bc + eps)**2)
      ar  = dr * (1.0d+00 + zt / (br + eps)**2)

! calculate weights
!
      aa  = (al + ar) + ac
      wl  = al / aa
      wr  = ar / aa
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! calculate βₖ (eq. 3.6 in [1])
!
      df  = abs(dfm(i))
      lq  = xi * df
      bl  = df2(ip1) + xi * abs(2.0d+00 * dfp(i) - dfp(ip1))
      bc  = df2(i  ) + lq
      br  = df2(im1) + lq

! calculate ζ (below eq. 3.6 in [1])

      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ (eq. 3.9 in [4])
!
      al  = dl * (1.0d+00 + zt / (bl + eps)**2)
      ac  = dc * (1.0d+00 + zt / (bc + eps)**2)
      ar  = dr * (1.0d+00 + zt / (br + eps)**2)

! normalize weights
!
      aa  = (al + ar) + ac
      wl  = al / aa
      wr  = ar / aa
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5ns
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5Z:
! -------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method and smoothness indicators by Borges et al. (2008).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Borges, R., Carmona, M., Costa, B., & Don, W.-S.,
!         "An improved weighted essentially non-oscillatory scheme for
!          hyperbolic conservation laws"
!         Journal of Computational Physics,
!         2008, vol. 227, pp. 3191-3211,
!         http://dx.doi.org/10.1016/j.jcp.2007.11.038
!
!===============================================================================
!
  subroutine reconstruct_crweno5z(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n)   :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! prepare smoothness indicators
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1   = i - 1
      ip1   = i + 1

! calculate βₖ (eqs. 9-11 in [1])
!
      bl    = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc    = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br    = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 25 in [1])
!
      tt    = abs(br - bl)

! calculate αₖ (eq. 28 in [1])
!
      al(i) = 1.0d+00 + tt / (bl + eps)
      ac(i) = 1.0d+00 + tt / (bc + eps)
      ar(i) = 1.0d+00 + tt / (br + eps)

    end do ! i = 2, n - 1

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng + 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! calculate weights
!
      wl  = cl * al(i)
      wc  = cc * ac(i)
      wr  = cr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i)
      wc  = cc * ac(i)
      wr  = cr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 2, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 2, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n - 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n))

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n))

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5z
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5YC:
! --------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method and smoothness indicators by Yamaleev & Carpenter (2009).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Yamaleev, N. K. & Carpenter, H. C.,
!         "A Systematic Methodology for Constructing High-Order Energy Stable
!          WENO Schemes"
!         Journal of Computational Physics,
!         2009, vol. 228, pp. 4248-4272,
!         http://dx.doi.org/10.1016/j.jcp.2009.03.002
!
!===============================================================================
!
  subroutine reconstruct_crweno5yc(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n)   :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! prepare smoothness indicators
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im2   = max(1, i - 2)
      im1   = i - 1
      ip1   = i + 1
      ip2   = min(n, i + 2)

! calculate βₖ (eqs. 9-11 in [1])
!
      bl    = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc    = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br    = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 64 in [3])
!
      tt  = (6.0d+00 * f(i) + (f(im2) + f(ip2))                                &
                                             - 4.0d+00 * (f(im1) + f(ip1)))**2

! calculate αₖ (eq. 28 in [1])
!
      al(i) = 1.0d+00 + tt / (bl + eps)
      ac(i) = 1.0d+00 + tt / (bc + eps)
      ar(i) = 1.0d+00 + tt / (br + eps)

    end do ! i = 2, n - 1

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng + 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! calculate weights
!
      wl  = cl * al(i)
      wc  = cc * ac(i)
      wr  = cr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i)
      wc  = cc * ac(i)
      wr  = cr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 2, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 2, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n))

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n))

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5yc
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5NS:
! --------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method combined with the smoothness indicators by Ha et al. (2013).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Ha, Y., Kim, C. H., Lee, Y. J., & Yoon, J.,
!         "An improved weighted essentially non-oscillatory scheme with a new
!          smoothness indicator",
!         Journal of Computational Physics,
!         2013, vol. 232, pp. 68-86
!         http://dx.doi.org/10.1016/j.jcp.2012.06.016
!
!===============================================================================
!
  subroutine reconstruct_crweno5ns(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: df, lq, l3, zt
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n,2) :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! the free parameter for smoothness indicators (see eq. 3.6 in [3])
!
    real(kind=8), parameter :: xi  =   4.0d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! 3rd order interpolation coefficients for three stencils
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = 0.5d+00 * abs(dfp(:) - dfm(:))

! prepare smoothness indicators
!
    do i = 1, n

! prepare neighbour indices
!
      im1  = max(1, i - 1)
      ip1  = min(n, i + 1)

! calculate βₖ
!
      df  = abs(dfp(i))
      lq  = xi * df
      bl  = df2(im1) + xi * abs(2.0d+00 * dfm(i) - dfm(im1))
      bc  = df2(i  ) + lq
      br  = df2(ip1) + lq

! calculate ζ
!
      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ
!
      al(i,1) = 1.0d+00 + zt / (bl + eps)**2
      ac(i,1) = 1.0d+00 + zt / (bc + eps)**2
      ar(i,1) = 1.0d+00 + zt / (br + eps)**2

! calculate βₖ
!
      df  = abs(dfm(i))
      lq  = xi * df
      bl  = df2(im1) + lq
      bc  = df2(i  ) + lq
      br  = df2(ip1) + xi * abs(2.0d+00 * dfp(i) - dfp(ip1))

! calculate ζ

      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ
!
      al(i,2) = 1.0d+00 + zt / (bl + eps)**2
      ac(i,2) = 1.0d+00 + zt / (bc + eps)**2
      ar(i,2) = 1.0d+00 + zt / (br + eps)**2

    end do ! i = 1, n

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng

! prepare neighbour indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! calculate weights
!
      wl  = cl * al(i,1)
      wc  = cc * ac(i,1)
      wr  = cr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i,2)
      wc  = cc * ac(i,2)
      wr  = cr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = 1, n

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 1, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i,1)
      wc  = dc * ac(i,1)
      wr  = dr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 1, ng

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i,1)
      wc  = dc * ac(i,1)
      wr  = dr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 1, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i,2)
      wc  = dc * ac(i,2)
      wr  = dr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 1, ng + 1

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i,2)
      wc  = dc * ac(i,2)
      wr  = dr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 1, ng + 1

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n))

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n))

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    fl(1) = fr(1)
    fr(n) = fl(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5ns
!
!===============================================================================
!
! subroutine RECONSTRUCT_MP5:
! --------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Monotonicity Preserving (MP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_mp5(n, h, f, fl, fr)

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: df, ds, dc0, dc4, dm1, dp1, dml, dmr
    real(kind=8) :: flc, fmd, fmp, fmn, fmx, ful
    real(kind=8) :: sigma

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! obtain the face values using high order interpolation
!
    do i = 2, n - 1

      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

      fl(i) = (4.7d+01 * f(i  ) + (2.7d+01 * f(ip1) - 1.3d+01 * f(im1))        &
                                - (3.0d+00 * f(ip2) - 2.0d+00 * f(im2)))       &
                                                                     / 6.0d+01
      fr(i) = (4.7d+01 * f(i  ) + (2.7d+01 * f(im1) - 1.3d+01 * f(ip1))        &
                                - (3.0d+00 * f(im2) - 2.0d+00 * f(ip2)))       &
                                                                     / 6.0d+01

    end do ! i = 2, n - 1

! apply monotonicity preserving limiting
!
    do i = 2, n - 1

      im1 = i - 1
      ip1 = i + 1

      if (dfm(i) * dfp(i) >= 0.0d+00) then
        sigma = kappa
      else
        sigma = kbeta
      end if

! get the limiting condition for the left state
!
      df    = sigma * dfm(i)
      fmp   = f(i) + minmod(dfp(i), df)
      ds    = (fl(i) - f(i)) * (fl(i) - fmp)

! limit the left state
!
      if (ds > eps) then

        dm1   = dfp(im1) - dfm(im1)
        dc0   = dfp(i  ) - dfm(i  )
        dp1   = dfp(ip1) - dfm(ip1)
        dc4   = 4.0d+00 * dc0

        dml   = 0.5d+00 * minmod4(dc4 - dm1, 4.0d+00 * dm1 - dc0, dc0, dm1)
        dmr   = 0.5d+00 * minmod4(dc4 - dp1, 4.0d+00 * dp1 - dc0, dc0, dp1)

        fmd   = f(i) + 0.5d+00 * dfp(i) - dmr
        ful   = f(i) +           df
        flc   = f(i) + 0.5d+00 * df     + dml

        fmx   = max(min(f(i), f(ip1), fmd), min(f(i), ful, flc))
        fmn   = min(max(f(i), f(ip1), fmd), max(f(i), ful, flc))

        fl(i) = median(fl(i), fmn, fmx)

      end if

! get the limiting condition for the right state
!
      df    = sigma * dfp(i)
      fmp   = f(i) - minmod(dfm(i), df)
      ds    = (fr(i) - f(i)) * (fr(i) - fmp)

! limit the right state
!
      if (ds > eps) then

        dm1 = dfp(im1) - dfm(im1)
        dc0 = dfp(i  ) - dfm(i  )
        dp1 = dfp(ip1) - dfm(ip1)
        dc4 = 4.0d+00 * dc0

        dml = 0.5d+00 * minmod4(dc4 - dm1, 4.0d+00 * dm1 - dc0, dc0, dm1)
        dmr = 0.5d+00 * minmod4(dc4 - dp1, 4.0d+00 * dp1 - dc0, dc0, dp1)

        fmd   = f(i) - 0.5d+00 * dfm(i) - dml
        ful   = f(i) -           df
        flc   = f(i) - 0.5d+00 * df     + dmr

        fmx   = max(min(f(i), f(im1), fmd), min(f(i), ful, flc))
        fmn   = min(max(f(i), f(im1), fmd), max(f(i), ful, flc))

        fr(i) = median(fr(i), fmn, fmx)

      end if

! shift the right state
!
      fr(im1) = fr(i)

    end do ! n = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_mp5
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRMP5:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact Reconstruction Monotonicity Preserving (CRMP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_crmp5(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

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
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: df, ds, dc0, dc4, dm1, dp1, dml, dmr
    real(kind=8) :: flc, fmd, fmp, fmn, fmx, ful
    real(kind=8) :: sigma

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! prepare the tridiagonal system coefficients for the interior
!
    do i = ng, n - ng + 1

      im1 = i - 1
      ip1 = i + 1

      a(i,1) = 3.0d-01
      b(i,1) = 6.0d-01
      c(i,1) = 1.0d-01

      a(i,2) = 1.0d-01
      b(i,2) = 6.0d-01
      c(i,2) = 3.0d-01

      r(i,1) = (f(im1) + 1.9d+01 * f(i  ) + 1.0d+01 * f(ip1)) / 3.0d+01
      r(i,2) = (f(ip1) + 1.9d+01 * f(i  ) + 1.0d+01 * f(im1)) / 3.0d+01

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit method (left-side reconstruction)
!
    do i = 2, ng

      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00

      r(i,1) = (4.7d+01 * f(i  ) + (2.7d+01 * f(ip1) - 1.3d+01 * f(im1))       &
                                 - (3.0d+00 * f(ip2) - 2.0d+00 * f(im2)))      &
                                                                     / 6.0d+01

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

    do i = n - ng, n - 1

      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00

      r(i,1) = (4.7d+01 * f(i  ) + (2.7d+01 * f(ip1) - 1.3d+01 * f(im1))       &
                                 - (3.0d+00 * f(ip2) - 2.0d+00 * f(im2)))      &
                                                                     / 6.0d+01

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit method (right-side reconstruction)
!
    do i = 2, ng + 1

      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00

      r(i,2) = (4.7d+01 * f(i  ) + (2.7d+01 * f(im1) - 1.3d+01 * f(ip1))       &
                                 - (3.0d+00 * f(im2) - 2.0d+00 * f(ip2)))      &
                                                                     / 6.0d+01

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

    do i = n - ng + 1, n - 1

      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00

      r(i,2) = (4.7d+01 * f(i  ) + (2.7d+01 * f(im1) - 1.3d+01 * f(ip1))       &
                                 - (3.0d+00 * f(im2) - 2.0d+00 * f(ip2)))      &
                                                                     / 6.0d+01

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n))

! apply the monotonicity preserving limiting
!
    do i = 2, n - 1

      im1 = i - 1
      ip1 = i + 1

      if (dfm(i) * dfp(i) >= 0.0d+00) then
        sigma = kappa
      else
        sigma = kbeta
      end if

      df    = sigma * dfm(i)
      fmp   = f(i) + minmod(dfp(i), df)
      ds    = (u(i) - f(i)) * (u(i) - fmp)

      if (ds <= eps) then

        fl(i) = u(i)

      else

        dm1   = dfp(im1) - dfm(im1)
        dc0   = dfp(i  ) - dfm(i  )
        dp1   = dfp(ip1) - dfm(ip1)
        dc4   = 4.0d+00 * dc0

        dml   = 0.5d+00 * minmod4(dc4 - dm1, 4.0d+00 * dm1 - dc0, dc0, dm1)
        dmr   = 0.5d+00 * minmod4(dc4 - dp1, 4.0d+00 * dp1 - dc0, dc0, dp1)

        fmd   = f(i) + 0.5d+00 * dfp(i) - dmr
        ful   = f(i) +           df
        flc   = f(i) + 0.5d+00 * df     + dml

        fmx   = max(min(f(i), f(ip1), fmd), min(f(i), ful, flc))
        fmn   = min(max(f(i), f(ip1), fmd), max(f(i), ful, flc))

        fl(i) = median(u(i), fmn, fmx)

      end if

    end do ! i = 2, n - 1

! solve the tridiagonal system of equations for the right-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n))

! apply the monotonicity preserving limiting
!
    do i = 2, n - 1

      im1 = i - 1
      ip1 = i + 1

      if (dfm(i) * dfp(i) >= 0.0d+00) then
        sigma = kappa
      else
        sigma = kbeta
      end if

      df    = sigma * dfp(i)
      fmp   = f(i) - minmod(dfm(i), df)

      ds    = (u(i) - f(i)) * (u(i) - fmp)

      if (ds <= eps) then

        fr(i) = u(i)

      else

        dm1 = dfp(im1) - dfm(im1)
        dc0 = dfp(i  ) - dfm(i  )
        dp1 = dfp(ip1) - dfm(ip1)
        dc4 = 4.0d+00 * dc0

        dml = 0.5d+00 * minmod4(dc4 - dm1, 4.0d+00 * dm1 - dc0, dc0, dm1)
        dmr = 0.5d+00 * minmod4(dc4 - dp1, 4.0d+00 * dp1 - dc0, dc0, dp1)

        fmd   = f(i) - 0.5d+00 * dfm(i) - dml
        ful   = f(i) -           df
        flc   = f(i) - 0.5d+00 * df     + dmr

        fmx   = max(min(f(i), f(im1), fmd), min(f(i), ful, flc))
        fmn   = min(max(f(i), f(im1), fmd), max(f(i), ful, flc))

        fr(i) = median(u(i), fmn, fmx)

      end if

! shift the right state
!
      fr(im1) = fr(i)

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crmp5
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
! function MINMOD:
! ===============
!
!   Function returns the minimum module value among two arguments.
!
!   Arguments:
!
!     a, b - the input values;
!
!===============================================================================
!
  real(kind=8) function minmod(a, b)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b
!
!-------------------------------------------------------------------------------
!
    minmod = (sign(0.5d+00, a) + sign(0.5d+00, b)) * min(abs(a), abs(b))

    return

!-------------------------------------------------------------------------------
!
  end function minmod
!
!===============================================================================
!
! function MINMOD4:
! ================
!
!   Function returns the minimum module value among four arguments.
!
!   Arguments:
!
!     a, b, c, d - the input values;
!
!===============================================================================
!
  real(kind=8) function minmod4(a, b, c, d)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b, c, d
!
!-------------------------------------------------------------------------------
!
    minmod4 = minmod(minmod(a, b), minmod(c, d))

    return

!-------------------------------------------------------------------------------
!
  end function minmod4
!
!===============================================================================
!
! function MEDIAN:
! ===============
!
!   Function returns the median of three argument values.
!
!   Arguments:
!
!     a, b, c - the input values;
!
!===============================================================================
!
  real(kind=8) function median(a, b, c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b, c
!
!-------------------------------------------------------------------------------
!
    median = a + minmod(b - a, c - a)

    return

  end function median
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
!
!===============================================================================
!
! subroutine CLIP_EXTREMA:
! -----------------------
!
!   Subroutine scans the reconstructed states and check if they didn't leave
!   the allowed limits.  In the case where the limits where exceeded,
!   the states are limited using constant reconstruction.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     f      - the cell centered integrals of variable;
!     fl, fr - the left and right states of variable;
!
!===============================================================================
!
  subroutine clip_extrema(n, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(in)    :: f
    real(kind=8), dimension(n), intent(inout) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, ip2
    real(kind=8) :: fmn, fmx
    real(kind=8) :: dfl, dfr, df
!
!------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for extrema clipping
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all points
!
    do i = 1, n

! calculate indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! estimate the bounds of the allowed interval for reconstructed states
!
      fmn = min(f(i), f(ip1))
      fmx = max(f(i), f(ip1))

! check if the left state lays in the allowed range
!
      if (fl(i) < fmn .or. fl(i) > fmx) then

! calculate the left and right derivatives
!
        dfl = f(i  ) - f(im1)
        dfr = f(ip1) - f(i  )

! get the limited slope
!
        df  = limiter_clip(0.5d+00, dfl, dfr)

! calculate new states
!
        fl(i  ) = f(i  ) + df
        fr(im1) = f(i  ) - df

      end if

! check if the right state lays in the allowed range
!
      if (fr(i) < fmn .or. fr(i) > fmx) then

! calculate the missing index
!
        ip2 = min(n, i + 2)

! calculate the left and right derivatives
!
        dfl = f(ip1) - f(i  )
        dfr = f(ip2) - f(ip1)

! get the limited slope
!
        df  = limiter_clip(0.5d+00, dfl, dfr)

! calculate new states
!
        fl(ip1) = f(ip1) + df
        fr(i  ) = f(ip1) - df

      end if

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for extrema clipping
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine clip_extrema

!===============================================================================
!
end module interpolations
