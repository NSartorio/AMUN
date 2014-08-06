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
  procedure(reconstruct)       , pointer, save :: reconstruct_states => null()
  procedure(stencil_weights_js), pointer, save :: stencil_weights    => null()
  procedure(limiter_zero)      , pointer, save :: limiter            => null()

! module parameters
!
  real(kind=8), save :: eps        = epsilon(1.0d+00)
  real(kind=8), save :: rad        = 0.5d+00

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
    character(len=255) :: sweights        = "yc"
    character(len=255) :: slimiter        = "mm"
    character(len=255) :: positivity_fix  = "off"
    character(len=255) :: clip_extrema    = "off"
    character(len=255) :: name_rec        = ""
    character(len=255) :: name_wei        = ""
    character(len=255) :: name_lim        = ""
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('interpolations:: initialization', imi)
    call set_timer('interpolations:: reconstruction', imr)
    call set_timer('interpolations:: fix positivity', imf)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! obtain the user defined interpolation methods and coefficients
!
    call get_parameter_string("reconstruction" , sreconstruction)
    call get_parameter_string("stencil_weights", sweights       )
    call get_parameter_string("limiter"        , slimiter       )
    call get_parameter_string("fix_positivity" , positivity_fix )
    call get_parameter_string("clip_extrema"   , clip_extrema   )
    call get_parameter_real  ("eps"            , eps            )
    call get_parameter_real  ("limo3_rad"      , rad            )

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
    case ("weno5", "WENO5")
      name_rec           =  "5th order WENO"
      reconstruct_states => reconstruct_weno5
    case ("weno5yc", "weno5-yc", "WENO5YC", "WENO5-YC")
      name_rec           =  "5th order WENO-YC (Yamaleev & Carpenter 2013)"
      reconstruct_states => reconstruct_weno5yc
    case ("weno5ns", "weno5-ns", "WENO5NS", "WENO5-NS")
      name_rec           =  "5th order WENO-NS (Ha et al. 2013)"
      reconstruct_states => reconstruct_weno5ns
    case ("crweno5", "CRWENO5")
      name_rec           =  "5th order Compact WENO"
      reconstruct_states => reconstruct_crweno5
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected reconstruction method is not " //     &
                           "implemented: " // trim(sreconstruction)
        stop
      end if
    end select

! select the stencil weights
!
    select case(trim(sweights))
    case ("js", "JS")
      name_wei        =  "Jiang-Shu"
      stencil_weights => stencil_weights_js
    case ("z", "Z")
      name_wei        =  "Borges et al."
      stencil_weights => stencil_weights_z
    case ("yc", "YC")
      name_wei        =  "Yamaleev-Carpenter"
      stencil_weights => stencil_weights_yc
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected stencil weight method is not " //     &
                           "implemented: " // trim(sweights)
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
    select case(trim(clip_extrema))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      clip = .true.
    case default
      clip = .false.
    end select

! print informations about the reconstruction methods and parameters
!
    if (verbose) then

      write (*,"(4x,a15,8x,'=',1x,a)") "reconstruction ", trim(name_rec)
      select case(trim(sreconstruction))
      case ("weno5", "WENO5", "crweno5", "CRWENO5")
        write (*,"(4x,a15,8x,'=',1x,a)") "stencil weights", trim(name_wei)
      case default
      end select
      write (*,"(4x,a15,8x,'=',1x,a)") "limiter        ", trim(name_lim)
      write (*,"(4x,a15,8x,'=',1x,a)") "fix positivity ", trim(positivity_fix)
      write (*,"(4x,a15,8x,'=',1x,a)") "clip extrema   ", trim(clip_extrema)

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
    real(kind=8), parameter :: dp = 2.0d+00 / 3.0d+00, dm = 1.0d+00 / 3.0d+00
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
! subroutine RECONSTRUCT_WENO5:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Reconstruction Weighted Essentially Non-Oscillatory (WENO5)
!   method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Jiang, G.-S., Shu, C.-W.,
!         "Efficient Implementation of Weighted ENO Schemes"
!         Journal of Computational Physics,
!         1996, vol. 126, pp. 202-228,
!         http://dx.doi.org/10.1006/jcph.1996.0130
!
!===============================================================================
!
  subroutine reconstruct_weno5(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local arrays
!
    real(kind=8), dimension(1:n)   :: wl, wc, wr
    real(kind=8), dimension(1:n)   :: u
!
!-------------------------------------------------------------------------------
!
! calculate stencil weights
!
    call stencil_weights(n, f(1:n), wl(1:n), wc(1:n), wr(1:n))

! find the left state interpolation implicitelly
!
    call weno5_explicit(n, f(1:n), wl(1:n), wc(1:n), wr(1:n), u(1:n))

! substitute the left state
!
    fl(1:n) = u(1:n)

! find the right state interpolation implicitelly
!
    call weno5_explicit(n, f(n:1:-1), wr(n:1:-1), wc(n:1:-1), wl(n:1:-1)       &
                                                                  , u(n:1:-1))

! substitute the right state
!
    fr(1:n-1) = u (2:n)
    fr(  n  ) = fl(  n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5YC:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with new
!   smoothness indicators and stencil weights by Yamaleev & Carpenter (2009).
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
!     [2] Arshed, G. M. & Hoffmann, K. A.,
!         "Minimizing errors from linear and nonlinear weights of WENO scheme
!          for broadband applications with shock waves",
!         Journal of Computational Physics,
!         2013, vol. 246, pp. 58-77
!         http://dx.doi.org/10.1016/j.jcp.2013.03.037
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

! improved weight coefficients (Table 1 in [2])
!
    real(kind=8), parameter :: dl = 1.235341937d-01, dr = 3.699651429d-01      &
                             , dc = 5.065006634d-01

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
    do i = 1, n

! prepare neighbour indices
!
      im1 = max(1, i - 1)
      im2 = max(1, i - 2)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! calculate βₖ (eq. 19 in [1])
!
      bl  = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc  = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br  = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 64 in [1])
!
      tt  = (6.0d+00 * f(i) + (f(im2) + f(ip2))                                &
                                             - 4.0d+00 * (f(im1) + f(ip1)))**2

! calculate αₖ (eq. 58 in [1])
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

! calculate the interpolations of the left state (eq. 15 in [1])
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

! calculate the interpolations of the right state (eq. 15 in [1])
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 1, n

! update the interpolation of the first and last points
!
    fl(1) = fr(1)
    fr(n) = fl(n)

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
!   smoothness indicators and stencil weights by He et al. (2013).
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
!     [2] Arshed, G. M. & Hoffmann, K. A.,
!         "Minimizing errors from linear and nonlinear weights of WENO scheme
!          for broadband applications with shock waves",
!         Journal of Computational Physics,
!         2013, vol. 246, pp. 58-77
!         http://dx.doi.org/10.1016/j.jcp.2013.03.037
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

! improved weight coefficients (Table 1 in [2])
!
    real(kind=8), parameter :: dl = 1.235341937d-01, dr = 3.699651429d-01      &
                             , dc = 5.065006634d-01

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
    do i = 1, n

! prepare neighbour indices
!
      im1 = max(1, i - 1)
      im2 = max(1, i - 2)
      ip1 = min(n, i + 1)
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

! calculate the interpolations of the left state (eq. 15 in [1])
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

! calculate the interpolations of the right state (eq. 15 in [1])
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 1, n

! update the interpolation of the first and last points
!
    fl(1) = fr(1)
    fr(n) = fl(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5ns
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method.
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
!
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!
!===============================================================================
!
  subroutine reconstruct_crweno5(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local arrays
!
    real(kind=8), dimension(1:n)   :: wl, wc, wr
    real(kind=8), dimension(1:n)   :: u
!
!-------------------------------------------------------------------------------
!
! calculate stencil weights
!
    call stencil_weights(n, f(1:n), wl(1:n), wc(1:n), wr(1:n))

! find the left state interpolation implicitelly
!
    call crweno5_implicit(n, f(1:n), wl(1:n), wc(1:n), wr(1:n), u(1:n))

! substitute the left state
!
    fl(1:n) = u(1:n)

! find the right state interpolation implicitelly
!
    call crweno5_implicit(n, f(n:1:-1), wr(n:1:-1), wc(n:1:-1), wl(n:1:-1)     &
                                                                  , u(n:1:-1))

! substitute the right state
!
    fr(1:n-1) = u (2:n)
    fr(  n  ) = fl(  n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5
!
!===============================================================================
!
! subroutine WENO5_EXPLICIT:
! -------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Reconstruction Weighted Essentially Non-Oscillatory (WENO5)
!   method (see [1]). This subroutine uses improved weights from [2].
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Jiang, G.-S., Shu, C.-W.,
!         "Efficient Implementation of Weighted ENO Schemes"
!         Journal of Computational Physics,
!         1996, vol. 126, pp. 202-228,
!         http://dx.doi.org/10.1006/jcph.1996.0130
!     [2] Arshed, G. M. & Hoffmann, K. A.,
!         "Minimizing errors from linear and nonlinear weights of WENO scheme
!          for broadband applications with shock waves",
!         Journal of Computational Physics,
!         2013, 246, 58-77
!         http://dx.doi.org/10.1016/j.jcp.2013.03.037
!
!===============================================================================
!
  subroutine weno5_explicit(n, f, wl, wc, wr, u)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(in)  :: wl, wc, wr
    real(kind=8), dimension(n), intent(out) :: u

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: xl, xc, xr, xx
    real(kind=8) :: ql, qc, qr

! improved weight coefficients (Table 1 in [2])
!
    real(kind=8), parameter :: dl = 1.235341937d-01, dr = 3.699651429d-01      &
                             , dc = 5.065006634d-01

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
! prepare coefficients
!
    do i = 1, n

! prepare indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! normalize weigths
!
      xc = dc * wc(i)
      xl = dl * wl(i)
      xr = dr * wr(i)
      xx = (xl + xr) + xc
      xl = xl / xx
      xr = xr / xx
      xc = 1.0d+00 - (xl + xr)

! calculate the interpolations of the left state (eq. 15 in [1])
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the interpolation of the left state
!
      u(i) = (xl * ql + xr * qr) + xc * qc

    end do ! i = 1, n

!-------------------------------------------------------------------------------
!
  end subroutine weno5_explicit
!
!===============================================================================
!
! subroutine CRWENO5_IMPLICIT:
! ---------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method.
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
!
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!
!===============================================================================
!
  subroutine crweno5_implicit(n, f, wl, wc, wr, u)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(in)  :: wl, wc, wr
    real(kind=8), dimension(n), intent(out) :: u

! local variables
!
    integer      :: im2, im1, i  , ip1, ip2
    integer      :: nm1, nm2, nm3, nm4
    real(kind=8) :: xl, xc, xr, xx
    real(kind=8) :: bl, bc, br, bt

! local arrays
!
    real(kind=8), dimension(1:n)   :: fm, fp
    real(kind=8), dimension(1:n)   :: a, b, c
    real(kind=8), dimension(1:n)   :: r, g

! local constants
!
    real(kind=8), parameter :: al = 1.0d-01, ac = 6.0d-01, ar = 3.0d-01
    real(kind=8), parameter :: cl = 2.0d-01, cc = 5.0d-01, cr = 3.0d-01
    real(kind=8), parameter :: dh = 5.0d-01, ds = 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate indices
!
    nm1 = n - 1
    nm2 = n - 2
    nm3 = n - 3
    nm4 = n - 4

! prepare coefficients
!
    do i = 4, nm3

! prepare indices
!
      im1 = i - 1
      ip1 = i + 1

! normalize weigths
!
      xc = cc * wc(i)
      xl = cl * wl(i)
      xr = cr * wr(i)
      xx = xc + (xl + xr)
      bl = xl / xx
      br = xr / xx
      bc = 1.0d+00 - (bl + br)

! calculate tridiagonal matrix coefficients
!
      a(i) = 2.0d+00 * bl +            bc
      b(i) =           bl + 2.0d+00 * (bc + br)
      c(i) =                                br

! prepare right hand side of tridiagonal equation
!
      r(i) = ( bl                        * f(im1)                              &
           +  (5.0d+00 * (bl + bc) + br) * f(i)                                &
           +  (bc + 5.0d+00 * br)        * f(ip1)) * dh

    end do

! at the left boundaries, we apply 5th order explicit WENO interpolation with
! different weights
!
! normalize weigths
!
    xc = ac * wc(3)
    xl = al * wl(3)
    xr = ar * wr(3)
    xx = xc + (xl + xr)
    bl = xl / xx
    br = xr / xx
    bc = 1.0d+00 - (bl + br)

! prepare right hand side of tridiagonal equation
!
    r(1) = dh * (f(1) + f(2  ))
    r(2) = f(2) + limiter(dh, f(2) - f(1), f(3) - f(2))
    r(3) = (bl * (2.0d+00 * f(1) - 7.0d+00 * f(2) + 1.1d+01 * f(3))            &
         +  bc * (        - f(2) + 5.0d+00 * f(3) + 2.0d+00 * f(4))            &
         +  br * (2.0d+00 * f(3) + 5.0d+00 * f(4) -           f(5))) * ds

! at the right boundaries, we do the similar thing
!
! normalize weigths
!
    xc = ac * wc(nm2)
    xl = al * wl(nm2)
    xr = ar * wr(nm2)
    xx = xc + (xl + xr)
    bl = xl / xx
    br = xr / xx
    bc = 1.0d+00 - (bl + br)

! prepare right hand side of tridiagonal equation
!
    r(nm2) = (bl * (2.0d+00 * f(nm4) - 7.0d+00 * f(nm3) + 1.1d+01 * f(nm2))    &
           +  bc * (        - f(nm3) + 5.0d+00 * f(nm2) + 2.0d+00 * f(nm1))    &
           +  br * (2.0d+00 * f(nm2) + 5.0d+00 * f(nm1) -           f(n  )))   &
            * ds

    r(nm1) = f(nm1) + limiter(dh, f(nm1) - f(nm2), f(n) - f(nm1))
    r(n  ) = f(n  ) + dh * (f(n) - f(nm1))

! correct matrix coefficients for boundaries with explicit interpolation
!
    do i = 1, 3
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do

    do i = nm2, n
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do

! solve the tridiagonal system of equations
!
    bt   = b(1)
    u(1) = r(1) / bt
    do i = 2, n
      im1  = i - 1
      g(i) = c(im1) / bt
      bt   = b(i) - a(i) * g(i)
      u(i) = (r(i) - a(i) * u(im1)) / bt
    end do
    do i = nm1, 1, -1
      ip1  = i + 1
      u(i) = u(i) - g(ip1) * u(ip1)
    end do


!-------------------------------------------------------------------------------
!
  end subroutine crweno5_implicit
!
!===============================================================================
!
! subroutine SMOOTHNESS_INDICATORS_JS:
! -----------------------------------
!
!   Subroutine calculates Jiang-Shu smoothness indicators for a given vector
!   of variable values.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     f  - the input vector of cell averaged values;
!     bl - the smoothness indicators for the left stencil;
!     bc - the smoothness indicators for the central stencil;
!     br - the smoothness indicators for the right stencil;
!
!   References:
!
!     [1] Jiang, G.-S., Shu, C.-W.,
!         "Efficient Implementation of Weighted ENO Schemes"
!         Journal of Computational Physics,
!         1996, vol. 126, pp. 202-228,
!         http://dx.doi.org/10.1006/jcph.1996.0130
!
!===============================================================================
!
  subroutine smoothness_indicators_js(n, f, bl, bc, br)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: bl, bc, br

! local variables
!
    integer :: nm1, np1
    integer :: i  , ip1

! local arrays
!
    real(kind=8), dimension(0:n+1) :: df2
    real(kind=8), dimension(0:n  ) :: dfm
    real(kind=8), dimension(1:n+1) :: dfp

! local constants
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 1.0d+00 / 4.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate indices
!
    np1 = n + 1
    nm1 = n - 1

! calculate the left and right first order derivatives
!
    do i = 1, nm1
      ip1 = i + 1

      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1  ) = dfp(1)
    dfm(0  ) = dfm(1)
    dfp(n  ) = dfm(n)
    dfp(np1) = dfp(n)

! the second order derivative
!
    df2(1:n  ) = c1 * (dfp(1:n) - dfm(1:n))**2
    df2(0    ) = df2(1)
    df2(  np1) = df2(n)

! calculate the left, central and right smoothness indicators
!
    bl(1:n) = df2(0:nm1) + c2 * (3.0d+00 * dfm(1:n) - dfm(0:nm1))**2
    bc(1:n) = df2(1:n  ) + c2 * (          dfp(1:n) + dfm(1:n  ))**2
    br(1:n) = df2(2:np1) + c2 * (3.0d+00 * dfp(1:n) - dfp(2:np1))**2

!-------------------------------------------------------------------------------
!
  end subroutine smoothness_indicators_js
!
!===============================================================================
!
! subroutine STENCIL_WEIGHTS_JS:
! -----------------------------
!
!   Subroutine calculate the stencil weights using Jiang-Shu method.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     f  - the input vector of cell averaged values;
!     wl - the weights the left stencil;
!     wc - the weights for the central stencil;
!     wr - the weights for the right stencil;
!
!   References:
!
!     [1] Jiang, G.-S., Shu, C.-W.,
!         "Efficient Implementation of Weighted ENO Schemes"
!         Journal of Computational Physics,
!         1996, vol. 126, pp. 202-228,
!         http://dx.doi.org/10.1006/jcph.1996.0130
!
!===============================================================================
!
  subroutine stencil_weights_js(n, f, wl, wc, wr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: wl, wc, wr

! local arrays
!
    real(kind=8), dimension(n) :: bl, bc, br
!
!-------------------------------------------------------------------------------
!
! calculate smoothness indicators according to Jiang-Sho
!
    call smoothness_indicators_js(n, f(1:n), bl(1:n), bc(1:n), br(1:n))

! calculate the weights
!
    wl(1:n) = 1.0d+00 / (bl(1:n) + eps)**2
    wc(1:n) = 1.0d+00 / (bc(1:n) + eps)**2
    wr(1:n) = 1.0d+00 / (br(1:n) + eps)**2

!-------------------------------------------------------------------------------
!
  end subroutine stencil_weights_js
!
!===============================================================================
!
! subroutine STENCIL_WEIGHTS_Z:
! ----------------------------
!
!   Subroutine calculate the stencil weights using Borges et al. method.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     f  - the input vector of cell averaged values;
!     wl - the weights the left stencil;
!     wc - the weights for the central stencil;
!     wr - the weights for the right stencil;
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
  subroutine stencil_weights_z(n, f, wl, wc, wr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: wl, wc, wr

! local arrays
!
    real(kind=8), dimension(n) :: bl, bc, br
    real(kind=8), dimension(n) :: tt
!
!-------------------------------------------------------------------------------
!
! calculate smoothness indicators according to Jiang-Sho
!
    call smoothness_indicators_js(n, f(1:n), bl(1:n), bc(1:n), br(1:n))

! calculate the factor τ
!
    tt(1:n) = abs(bl(1:n) - br(1:n))

! calculate the weights
!
    wl(1:n) = 1.0d+00 + (tt(1:n) / (bl(1:n) + eps))**2
    wc(1:n) = 1.0d+00 + (tt(1:n) / (bc(1:n) + eps))**2
    wr(1:n) = 1.0d+00 + (tt(1:n) / (br(1:n) + eps))**2

!-------------------------------------------------------------------------------
!
  end subroutine stencil_weights_z
!
!===============================================================================
!
! subroutine STENCIL_WEIGHTS_YC:
! -----------------------------
!
!   Subroutine calculate the stencil weights using Yalamleev-Carpenter method.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     f  - the input vector of cell averaged values;
!     wl - the weights the left stencil;
!     wc - the weights for the central stencil;
!     wr - the weights for the right stencil;
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
  subroutine stencil_weights_yc(n, f, wl, wc, wr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: wl, wc, wr

! local variables
!
    integer      :: im2, im1, i  , ip1, ip2

! local arrays
!
    real(kind=8), dimension(n) :: bl, bc, br
    real(kind=8), dimension(n) :: tt
!
!-------------------------------------------------------------------------------
!
! calculate smoothness indicators according to Jiang-Sho
!
    call smoothness_indicators_js(n, f(1:n), bl(1:n), bc(1:n), br(1:n))

! calculate the factor τ
!
    do i = 1, n
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

      tt(i) = (6.0d+00 * f(i) + (f(im2) + f(ip2))                              &
                                             - 4.0d+00 * (f(im1) + f(ip1)))**2

    end do

! calculate the weights
!
    wl(1:n) = 1.0d+00 + (tt(1:n) / (bl(1:n) + eps))**2
    wc(1:n) = 1.0d+00 + (tt(1:n) / (bc(1:n) + eps))**2
    wr(1:n) = 1.0d+00 + (tt(1:n) / (br(1:n) + eps))**2

!-------------------------------------------------------------------------------
!
  end subroutine stencil_weights_yc
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
    integer      :: i, im1, ip1
    real(kind=8) :: fmn, fmx
!
!------------------------------------------------------------------------------
!
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

! calculate new states
!
        fl(i  ) = f(i  )
        fr(im1) = f(i  )

      end if

! check if the right state lays in the allowed range
!
      if (fr(i) < fmn .or. fr(i) > fmx) then

! calculate new states
!
        fl(ip1) = f(ip1)
        fr(i  ) = f(ip1)

      end if

    end do ! i = 1, n

!-------------------------------------------------------------------------------
!
  end subroutine clip_extrema

!===============================================================================
!
end module interpolations
