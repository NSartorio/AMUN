!!******************************************************************************
!!
!! module: INTERPOLATION - subroutines for different kinds of interpolation
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!******************************************************************************
!!
!!  This file is part of the AMUN code.
!!
!!  This program is free software; you can redistribute it and/or
!!  modify it under the terms of the GNU General Public License
!!  as published by the Free Software Foundation; either version 2
!!  of the License, or (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program; if not, write to the Free Software Foundation,
!!  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!!
!!******************************************************************************
!!
!
module interpolation

  implicit none

  contains
!
#ifdef TVD
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i
    real               :: df, ds
    real, dimension(n) :: dfl, dfr
!
!-------------------------------------------------------------------------------
!
!! second order TVD interpolation
!!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      dfr(i  ) = f(i+1) - f(i)
      dfl(i+1) = dfr(i)
    end do
    dfl(1) = dfr(1)
    dfr(n) = dfl(n)

! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n
#ifdef MINMOD
      df = 0.5d0 * minmod(dfr(i), dfl(i))

      fl(i) = f(i) + df
      fr(i) = f(i) - df
#endif /* MINMOD */
#ifdef LF
      ds = dfr(i) * dfl(i)

      if (ds .gt. eps) then
        df  = ds / (dfr(i) + dfl(i))

        fl(i) = f(i) + df
        fr(i) = f(i) - df
      else
        fl(i) = f(i)
        fr(i) = f(i)
      end if
#endif /* LF */
    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* TVD */
#ifdef CENO3
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer               :: i, im1, ip1, im2, ip2
    integer, dimension(1) :: loc
    real                  :: dfl, dfr, df
    real, dimension(3)    :: ql, qr, dl

! interpolation coefficients
!
    real, parameter    :: a11 =   0.375d0, a12 = - 1.250d0, a13 =   1.875d0    &
                        , a21 = - 0.125d0, a22 =   0.750d0, a23 =   0.375d0    &
                        , a31 =   0.375d0, a32 =   0.750d0, a33 = - 0.125d0
!
!-------------------------------------------------------------------------------
!
!! third order CENO interpolation
!!
! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n

! prepare indices
!
      im2   = max(1, i - 2)
      im1   = max(1, i - 1)
      ip1   = min(n, i + 1)
      ip2   = min(n, i + 2)

! calculate left and right derivatives
!
      dfr   = f(ip1) - f(i  )
      dfl   = f(i  ) - f(im1)

! calculate limited derivative
!
      df    = 0.5d0 * minmod(dfr, dfl)

! obtain the states with the linear interpolation
!
      fl(i) = f(i) + df
      fr(i) = f(i) - df

! calculate the second order interpolations of the left and right states
!
      ql(1) = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      ql(2) = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      ql(3) = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

      qr(1) = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qr(2) = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr(3) = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! select the second order interpolation the closest to the linear one
!
      dl(:) = ql(:) - fl(i)
      if(minval(dl) * maxval(dl) .gt. eps) then
        dl(2) = 0.7d0 * dl(2)
        loc   = minloc(abs(dl(:)))
        fl(i) = ql(loc(1))
      end if

      dl(:) = qr(:) - fr(i)
      if(minval(dl) * maxval(dl) .gt. eps) then
        dl(2) = 0.7d0 * dl(2)
        loc   = minloc(abs(dl(:)))
        fr(i) = qr(loc(1))
      end if

    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* CENO3 */
#ifdef WENO3
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation;
!              the subroutine implements the improved 3rd order WENO scheme as
!              described in Yamaleev & Carpenter, 2000, Journal of Computational
!              Physics, 228, 3025
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer         :: i, im1, ip1
    real            :: h2, dfp, dfm, fp, fm, bp, bm, ap, am, wp, wm, ww

! weight coefficients
!
    real, parameter :: dp = 2.0d0 / 3.0d0, dm = 1.0d0 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! prepare fixed parameters
!
    h2 = h * h

!! third order WENO interpolation
!!
! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n

! prepare indices
!
      im1   = max(1, i - 1)
      ip1   = min(n, i + 1)

! calculate left and right derivatives
!
      dfp    = f(ip1) - f(i  )
      dfm    = f(i  ) - f(im1)
      ww    = abs(dfp - dfm)**2

! calculate corresponding betas
!
      bp    = dfp * dfp
      bm    = dfm * dfm

! calculate improved alphas
!
      ap    = (1.0d0 + ww / (bp + h2))
      am    = (1.0d0 + ww / (bm + h2))

! calculate weights
!
      wp    = dp * ap
      wm    = dm * am
      ww    = wp + wm
      wp    = wp / ww
      wm    = wm / ww

! calculate right and left sides interpolations
!
      fp    = 0.5d0 * (  f(i  ) +         f(ip1))
      fm    = 0.5d0 * (- f(im1) + 3.0d0 * f(i  ))

! calculate the left state
!
      fl(i) = wp * fp + wm * fm

! calculate weights
!
      wp    = dp * am
      wm    = dm * ap
      ww    = wp + wm
      wp    = wp / ww
      wm    = wm / ww

! calculate right and left sides interpolations
!
      fp    = 0.5d0 * (  f(i  ) +         f(im1))
      fm    = 0.5d0 * (- f(ip1) + 3.0d0 * f(i  ))

! calculate the left state
!
      fr(i) = wp * fp + wm * fm

    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* WENO3 */
#ifdef WENOZ
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation;
!              the subroutine implements the improved 5rd order WENO scheme
!              called WENO-Z as described by Borges et al. 2008, Journal of
!              Computational Physics, 227, 3191
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i, im1, ip1, im2, ip2
    real               :: bl, bc, br, al, ac, ar, ww, wl, wc, wr
    real, dimension(3) :: ql, qr
    real, dimension(n) :: dfl, dfr, df2

! weight coefficients
!
    real, parameter :: c1 = 13.0d0 / 12.0d0, c2 = 0.25d0
    real, parameter :: dl = 0.1d0, dc = 0.6d0, dr = 0.3d0

! interpolation coefficients
!
    real, parameter :: a11 =  2.0d0/6.0d0, a12 = -7.0d0/6.0d0, a13 = 11.0d0/6.0d0 &
                     , a21 = -1.0d0/6.0d0, a22 =  5.0d0/6.0d0, a23 =  2.0d0/6.0d0 &
                     , a31 =  2.0d0/6.0d0, a32 =  5.0d0/6.0d0, a33 = -1.0d0/6.0d0
!
!-------------------------------------------------------------------------------
!
!! fifth order WENO-Z interpolation
!!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      dfr(i  ) = f(i+1) - f(i)
      dfl(i+1) = dfr(i)
    end do
    dfl(1) = dfr(1)
    dfr(n) = dfl(n)
    df2(:) = dfr(:) - dfl(:)

! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n

! prepare indices
!
      im2   = max(1, i - 2)
      im1   = max(1, i - 1)
      ip1   = min(n, i + 1)
      ip2   = min(n, i + 2)

! calculate betas
!
      bl    = c1 * df2(im1)**2 + c2 * (3.0d0 * dfl(i  ) - dfl(im1))**2
      bc    = c1 * df2(i  )**2 + c2 * (        dfr(i  ) + dfl(i  ))**2
      br    = c1 * df2(ip1)**2 + c2 * (3.0d0 * dfr(i  ) - dfr(ip1))**2
      ww    = abs(bl - br)

! calculate alphas
!
      al    = 1.0d0 + ww / (bl + eps)
      ac    = 1.0d0 + ww / (bc + eps)
      ar    = 1.0d0 + ww / (br + eps)

! calculate weights
!
      wl    = dl * al
      wc    = dc * ac
      wr    = dr * ar
      ww    = wl + wc + wr
      wl    = wl / ww
      wc    = wc / ww
      wr    = wr / ww

! calculate the second order interpolations of the left state
!
      ql(1) = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      ql(2) = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      ql(3) = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = wl * ql(1) + wc * ql(2) + wr * ql(3)

! calculate weights
!
      wl    = dl * ar
      wc    = dc * ac
      wr    = dr * al
      ww    = wl + wc + wr
      wl    = wl / ww
      wc    = wc / ww
      wr    = wr / ww

! calculate the second order interpolations of the right state
!
      qr(1) = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qr(2) = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr(3) = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = wl * qr(1) + wc * qr(2) + wr * qr(3)

    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* WENOZ */
#ifdef LIMO3
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps
    use config, only : rad

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i
    real               :: df, ds
    integer            :: im1, ip1
    real               :: rdx, rdx2, dfr, dfl, th, et, f0, f1, f2, ft, xi
!
!-------------------------------------------------------------------------------
!
!! third order interpolation
!!
! prepare parameters
!
    rdx  = rad * h
    rdx2 = rdx * rdx

    do i = 1, n

! prepare indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! prepare differences
!
      dfr = f(ip1) - f(i  )
      dfl = f(i  ) - f(im1)

      et = (dfl * dfl + dfr * dfr) / rdx2
      xi = max(0.0d0, min(1.0d0, 0.5d0 + 0.5d0 * (et - 1.0d0) / eps))

! calculate values at i+1/2
!
      if (dfr .ne. 0.0d0) then
        th = dfl / dfr
      else
        th = dfl / eps
      end if

      f1 = (2.0d0 + th) / 3.0d0

      if (th .ge. 0.0d0) then
        f2 = max(0.0d0, min(f1, 2.0d0 * th, 1.6d0))
      else
        f2 = max(0.0d0, min(f1, - 0.5d0 * th))
      end if

      f0 = f1 + xi * (f2 - f1)

      fl(i) = f(i) + 0.5d0 * dfr * f0

! calculate values at i-1/2
!
      if (dfl .ne. 0.0d0) then
        th = dfr / dfl
      else
        th = dfr / eps
      end if

      f1 = (2.0d0 + th) / 3.0d0

      if (th .ge. 0.0d0) then
        f2 = max(0.0d0, min(f1, 2.0d0 * th, 1.6d0))
      else
        f2 = max(0.0d0, min(f1, - 0.5d0 * th))
      end if

      f0 = f1 + xi * (f2 - f1)

      fr(im1) = f(i) - 0.5d0 * dfl * f0

    end do

! update boundaries
!
    fl(1) = f(1)
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* LIMO3 */
#ifdef MP
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

    use config, only : eps
    use config, only : alpha

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i
    real               :: df, ds
    integer            :: im2, im1, ip1, ip2
#ifndef MP5
    integer            :: im3, ip3
#endif /* !MP5 */
#ifdef MP9
    integer            :: im4, ip4
#endif /* MP9 */
    real               :: fh, fmp, fav, fmd, ful, flc, fmn, fmx
    real               :: dl, dr, dm1, dc0, dp1, dml, dmr

! parameters
!
    real, parameter    :: ac =   4.0d0 / 3.0d0
#ifdef MP5
    real, parameter    :: a1 =    2.0d0 / 60.0d0, a2 = - 13.0d0 / 60.0d0       &
                        , a3 =   47.0d0 / 60.0d0, a4 =   27.0d0 / 60.0d0       &
                        , a5 = -  3.0d0 / 60.0d0
#endif /* MP5 */
#ifdef MP7
    real, parameter    :: a1 = -  3.0d0 / 42.0d1, a2 =   25.0d0 / 42.0d1       &
                        , a3 = - 10.1d1 / 42.0d1, a4 =   31.9d1 / 42.0d1       &
                        , a5 =   21.4d1 / 42.0d1, a6 = - 38.0d0 / 42.0d1       &
                        , a7 =    4.0d0 / 42.0d1
#endif /* MP7 */
#ifdef MP9
    real, parameter    :: a1 =    4.0d0 / 25.2d2, a2 = - 41.0d0 / 25.2d2       &
                        , a3 =   19.9d1 / 25.2d2, a4 = - 64.1d1 / 25.2d2       &
                        , a5 =  187.9d1 / 25.2d2, a6 =  137.5d1 / 25.2d2       &
                        , a7 = - 30.5d1 / 25.2d2, a8 =   55.0d0 / 25.2d2       &
                        , a9 = -  5.0d0 / 25.2d2
#endif /* MP9 */
!
!-------------------------------------------------------------------------------
!
!! fifth or higher order monotonicity preserving interpolation
!!
! iterate over all positions
!
    do i = 1, n
#ifdef MP9
      im4 = max(1, i-4)
#endif /* MP9 */
#ifndef MP5
      im3 = max(1, i-3)
#endif /* !MP5 */
      im2 = max(1, i-2)
      im1 = max(1, i-1)
      ip1 = min(n, i+1)
      ip2 = min(n, i+2)
#ifndef MP5
      ip3 = min(n, i+3)
#endif /* !MP5 */
#ifdef MP9
      ip4 = min(n, i+4)
#endif /* MP9 */

! calculate values at i+1/2
!
#ifdef MP5
      fh = a1 * f(im2) + a2 * f(im1) + a3 * f(i  ) + a4 * f(ip1) + a5 * f(ip2)
#endif /* MP5 */
#ifdef MP7
      fh = a1 * f(im3) + a2 * f(im2) + a3 * f(im1) + a4 * f(i  ) + a5 * f(ip1) &
         + a5 * f(ip2) + a5 * f(ip3)
#endif /* MP7 */
#ifdef MP9
      fh = a1 * f(im4) + a2 * f(im3) + a3 * f(im2) + a4 * f(im1) + a5 * f(i  ) &
         + a6 * f(ip1) + a7 * f(ip2) + a8 * f(ip3) + a9 * f(ip4)
#endif /* MP9 */

      dl  = f(i  ) - f(im1)
      dr  = f(ip1) - f(i  )
      fmp = f(i) + minmod(dr, alpha * dl)
      ds = (fh - f(i)) * (fh - fmp)
      if (ds .le. eps) then
        fl(i) = fh
      else
        dm1 = f(i) + f(im2) - 2.0d0 * f(im1)
        dc0 = dr - dl
        dp1 = f(i) + f(ip2) - 2.0d0 * f(ip1)

        dml = minmod4(4.0d0 * dm1 - dc0, 4.0d0 * dc0 - dm1, dm1, dc0)
        dmr = minmod4(4.0d0 * dc0 - dp1, 4.0d0 * dp1 - dc0, dc0, dp1)

        ful = f(i) + alpha * dl
        fav = 0.5d0 * (f(i) + f(ip1))
        fmd = fav - 0.5d0 * dmr
        flc = f(i) + 0.5d0 * dl + ac * dml
        fmn = max(min(f(i),f(ip1),fmd), min(f(i),ful,flc))
        fmx = min(max(f(i),f(ip1),fmd), max(f(i),ful,flc))
        fl(i) = median(fh, fmn, fmx)
      end if

! calculate values at i-1/2
!
#ifdef MP5
      fh = a1 * f(ip2) + a2 * f(ip1) + a3 * f(i  ) + a4 * f(im1) + a5 * f(im2)
#endif /* MP5 */
#ifdef MP7
      fh = a1 * f(ip3) + a2 * f(ip2) + a3 * f(ip1) + a4 * f(i  ) + a5 * f(im1) &
         + a5 * f(im2) + a5 * f(im3)
#endif /* MP7 */
#ifdef MP9
      fh = a1 * f(ip4) + a2 * f(ip3) + a3 * f(ip2) + a4 * f(ip1) + a5 * f(i  ) &
         + a6 * f(im1) + a7 * f(im2) + a8 * f(im3) + a9 * f(im4)
#endif /* MP9 */

      dl  = f(i  ) - f(ip1)
      dr  = f(im1) - f(i  )
      fmp = f(i) + minmod(dr, alpha * dl)
      ds = (fh - f(i)) * (fh - fmp)
      if (ds .le. eps) then
        fr(i) = fh
      else
        dm1 = f(i) + f(ip2) - 2.0d0 * f(ip1)
        dc0 = dr - dl
        dp1 = f(i) + f(im2) - 2.0d0 * f(im1)

        dml = minmod4(4.0d0 * dm1 - dc0, 4.0d0 * dc0 - dm1, dm1, dc0)
        dmr = minmod4(4.0d0 * dc0 - dp1, 4.0d0 * dp1 - dc0, dc0, dp1)

        ful = f(i) + alpha * dl
        fav = 0.5d0 * (f(i) + f(im1))
        fmd = fav - 0.5d0 * dmr
        flc = f(i) + 0.5d0 * dl + ac * dml
        fmn = max(min(f(i),f(im1),fmd), min(f(i),ful,flc))
        fmx = min(max(f(i),f(im1),fmd), max(f(i),ful,flc))
        fr(i) = median(fh, fmn, fmx)
      end if

! if the left state creates a local extremum, limit both states to the zeroth
! interpolation
!
      if ((f(ip1) - fl(i)) * (fl(i) - f(i)) .le. eps) then
        fl(i) = f(i)
        fr(i) = f(i)
      end if

! if the right state creates a local extremum, limit both states to the zeroth
! interpolation
!
      if ((f(i) - fr(i)) * (fr(i) - f(im1)) .le. eps) then
        fl(i) = f(i)
        fr(i) = f(i)
      end if
    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
#endif /* MP */
!
!===============================================================================
!
! minmod: function returns the minimum module value among two arguments
!
!===============================================================================
!
  real function minmod(a, b)

    implicit none

! input arguments
!
    real, intent(in) :: a, b
!
!-------------------------------------------------------------------------------
!
    minmod = (sign(0.5d0, a) + sign(0.5d0, b)) * min(abs(a), abs(b))
    return
  end function minmod
!
!===============================================================================
!
! minmod4: function returns the minimum module value among four arguments
!
!===============================================================================
!
  real function minmod4(a, b, c, d)

    implicit none

! input arguments
!
    real, intent(in) :: a, b, c, d

! local variables
!
    real :: sb, sc, sd, s
!
!-------------------------------------------------------------------------------
!
    sb      = sign(1.0d0, a) + sign(1.0d0, b)
    sc      = sign(1.0d0, a) + sign(1.0d0, c)
    sd      = sign(1.0d0, a) + sign(1.0d0, d)
    s       = 0.125d0 * sb * abs(sc * sd)
    minmod4 = s * min(abs(a), abs(b), abs(c), abs(d))
    return
  end function minmod4
!
!===============================================================================
!
! median: function returns the median of three numbers
!
!===============================================================================
!
  real function median(a, b, c)

    implicit none

! input arguments
!
    real, intent(in) :: a, b, c
!
!-------------------------------------------------------------------------------
!
    median = a + minmod(b - a, c - a)

    return
  end function median
!
!===============================================================================
!
! limiter: function returns the minimum module value among two arguments
!
!===============================================================================
!
  real function limiter(a, b)

    implicit none

! input arguments
!
    real, intent(in) :: a, b

! local variables
!
    real             :: ds
!
!-------------------------------------------------------------------------------
!
    limiter = 0.0d0

    ds = a * b

    if (ds .gt. 0.0d0) then
#ifdef MINMOD
      limiter = sign(1.0d0, a) * min(abs(a), abs(b))
#else /* MINMOD */
      limiter = 2.0d0 * ds / (a + b)
#endif /* MINMOD */
    end if

    return
  end function limiter
!
!===============================================================================
!
! expand: expands a multi-dimentional array similar using only TVD or high order
!         interpolation
!
!===============================================================================
!
  subroutine expand(cm, dm, ng, u, v)

    implicit none

! input parameters
!
    integer, dimension(3)                     , intent(in)  :: cm, dm
    integer                                   , intent(in)  :: ng
    real        , dimension(cm(1),cm(2),cm(3)), intent(in)  :: u
    real        , dimension(dm(1),dm(2),dm(3)), intent(out) :: v

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real, dimension(:)    , allocatable :: x, y
    real, dimension(:,:,:), allocatable :: w, z
!
!-------------------------------------------------------------------------------
!
! allocate temporary arrays
!
    allocate(x(maxval(cm)))
    allocate(y(maxval(dm)))
    allocate(w(dm(1),cm(2),cm(3)))
    allocate(z(dm(1),dm(2),cm(3)))

! expand in X direction
!
    do k = 1, cm(3)
      do j = 1, cm(2)
        x(1:cm(1)) = u(1:cm(1),j,k)

        call expand_1d_const(cm(1), dm(1), ng, x(1:cm(1)), y(1:dm(1)))

        w(1:dm(1),j,k) = y(1:dm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, cm(3)
      do i = 1, dm(1)
        x(1:cm(2)) = w(i,1:cm(2),k)

        call expand_1d_const(cm(2), dm(2), ng, x(1:cm(2)), y(1:dm(2)))

        z(i,1:dm(2),k) = y(1:dm(2))
      end do
    end do

! expand in Z-direction
!
    if (cm(3) .gt. 1) then
      do j = 1, dm(2)
        do i = 1, dm(1)
          x(1:cm(3)) = z(i,j,1:cm(3))

          call expand_1d_const(cm(3), dm(3), ng, x(1:cm(3)), y(1:dm(3)))

          v(i,j,1:dm(3)) = y(1:dm(3))
        end do
      end do
    else
      v(:,:,:) = z(:,:,:)
    endif

! deallocate temporary arrays
!
    deallocate(w)
    deallocate(z)
    deallocate(x)
    deallocate(y)
!
!-------------------------------------------------------------------------------
!
  end subroutine expand
!
!===============================================================================
!
! shrink: shrinks multi-dimentional array using different interpolations
!
!===============================================================================
!
  subroutine shrink(dm, fm, u, v, xflag, yflag, zflag)

    implicit none

! input parameters
!
    integer, dimension(3)         , intent(in)  :: dm, fm
    real        , dimension(:,:,:), intent(in)  :: u
    real        , dimension(:,:,:), intent(out) :: v
    character                     , intent(in)  :: xflag, yflag, zflag

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real, dimension(:)    , allocatable :: x, y
    real, dimension(:,:,:), allocatable :: w, z
!
!-------------------------------------------------------------------------------
!
    v(:,:,:) = 0.0

    allocate(x(maxval(dm)))
    allocate(y(maxval(fm)))
    allocate(w(fm(1),dm(2),dm(3)))
    allocate(z(fm(1),fm(2),dm(3)))

! expand in X direction
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        x(1:dm(1)) = u(1:dm(1),j,k)

        call shrink_1d(dm(1), fm(1), x(1:dm(1)), y(1:fm(1)), xflag)

        w(1:fm(1),j,k) = y(1:fm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        x(1:dm(2)) = w(i,1:dm(2),k)

        call shrink_1d(dm(2), fm(2), x(1:dm(2)), y(1:fm(2)), yflag)

        z(i,1:fm(2),k) = y(1:fm(2))
      end do
    end do

! expand in Z-direction
!
    if (dm(3) .gt. 1) then
      do j = 1, fm(2)
        do i = 1, fm(1)
          x(1:dm(3)) = z(i,j,1:dm(3))

          call shrink_1d(dm(3), fm(3), x(1:dm(3)), y(1:fm(3)), zflag)

          v(i,j,1:fm(3)) = y(1:fm(3))
        end do
      end do
    else
      do k = 1, fm(3)
        do j = 1, fm(2)
          do i = 1, fm(1)
            v(i,j,k) = z(i,j,k)
          end do
        end do
      end do
    endif

    deallocate(w)
    deallocate(z)
    deallocate(x)
    deallocate(y)

!-------------------------------------------------------------------------------
!
  end subroutine shrink
!
!===============================================================================
!
! expand_1d: one dimensional expansion using the TVD or high ordred
!            interpolation
!
!===============================================================================
!
  subroutine expand_1d(nu, nv, ng, u, v)

    implicit none

! input parameters
!
    integer                    , intent(in)    :: nu, nv, ng
    real        , dimension(nu), intent(in)    :: u
    real        , dimension(nv), intent(inout) :: v

! local variables
!
    integer  :: i, ib, ie, il, ir
    real     :: dum, dup, ds, du
!
!-------------------------------------------------------------------------------
!
    ib =  1 + ng
    ie = nu - ng

    do i = ib, ie

      ir = 2 * (i - ng)
      il = ir - 1

#ifdef TVD
      dup = u(i+1) - u(i)
      dum = u(i) - u(i-1)
      ds  = dup * dum

      if (ds .gt. 0.0) then
#ifdef MINMOD
        du  = sign(0.25, dup) * min(abs(dup), abs(dum))
#endif /* MINMOD */
#ifdef LF
        du    = 0.5 * ds / (dup + dum)
#endif /* LF */
        v(il) = u(i) - du
        v(ir) = u(i) + du
      else
        v(il) = u(i)
        v(ir) = u(i)
      end if
#endif /* TVD */
#if defined LIMO3 || defined MP
      dup = u(i+1) - u(i)
      dum = u(i) - u(i-1)
      ds  = dup * dum

      if (ds .gt. 0.0) then
        du  = sign(0.25, dup) * min(abs(dup), abs(dum))
        v(il) = u(i) - du
        v(ir) = u(i) + du
      else
        v(il) = u(i)
        v(ir) = u(i)
      end if
#endif /* LIMO3 | MP */

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_1d
!
!===============================================================================
!
! expand_1d_const: one dimensional expansion using the constant interpolation
!
!===============================================================================
!
  subroutine expand_1d_const(nu, nv, ng, u, v)

    implicit none

! input parameters
!
    integer                    , intent(in)    :: nu, nv, ng
    real        , dimension(nu), intent(in)    :: u
    real        , dimension(nv), intent(inout) :: v

! local variables
!
    integer  :: i, ib, ie, il, ir
!
!-------------------------------------------------------------------------------
!
    ib =  1 + ng
    ie = nu - ng

    do i = ib, ie

      ir = 2 * (i - ng)
      il = ir - 1

      v(il) = u(i)
      v(ir) = u(i)

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_1d_const
!
!===============================================================================
!
! shrink_1d: one dimensional shrinking using different interpolations for
!            different location of the new points
!
!===============================================================================
!
  subroutine shrink_1d(n, m, u, v, flag)

    implicit none

! input parameters
!
    integer                   , intent(in)  :: n, m
    real        , dimension(n), intent(in)  :: u
    real        , dimension(m), intent(out) :: v
    character                 , intent(in)  :: flag

! local variables
!
    integer      :: i, il, ir
!
!-------------------------------------------------------------------------------
!
    v(:) = 0.0

    select case(flag)
    case('c')
      do i = 1, m
        ir   = 2 * i
        v(i) = u(ir)
      end do
    case default
      do i = 1, m
        ir   = 2 * i
        il   = ir - 1

        v(i) = 0.5 * (u(il) + u(ir))
      end do
    end select
!
!-------------------------------------------------------------------------------
!
  end subroutine shrink_1d
#ifdef MHD
!
!===============================================================================
!
! divergence: function calculates divergence of the input staggered field
!             [Bx(i+1/2,j,k), By(i,j+1/2,k), Bz(i,j,k+1/2)] -> div B
!
!===============================================================================
!
  subroutine divergence(b, db, dx, dy, dz)

    use config   , only : im, jm, km
    use variables, only : nvr

    implicit none

! input and output variables
!
    real, dimension(3,im,jm,km), intent(in)  :: b
    real, dimension(  im,jm,km), intent(out) :: db
    real, optional             , intent(in)  :: dx, dy, dz

! local  variables
!
    integer :: i, j, k
    real    :: dxi, dyi, dzi
!
!------------------------------------------------------------------------------
!
! check optional arguments
!
    dxi = 1.0
    dyi = 1.0
    dzi = 1.0

    if (present(dx)) dxi = 1.0 / dx
    if (present(dy)) dyi = 1.0 / dy
    if (present(dz)) dzi = 1.0 / dz

! reset output array
!
    db(:,:,:) = 0.0

! iterate over all points
!
#if NDIMS == 3
    do k = 2, km
      do j = 2, jm
        do i = 2, im
          db(i,j,k) = dxi * (b(1,i,j,k) - b(1,i-1,j,k))                        &
                    + dyi * (b(2,i,j,k) - b(2,i,j-1,k))                        &
                    + dzi * (b(3,i,j,k) - b(3,i,j,k-1))
#else /* NDIMS == 3 */
    do k = 1, km
      do j = 2, jm
        do i = 2, im
          db(i,j,k) = dxi * (b(1,i,j,k) - b(1,i-1,j,k))                        &
                    + dyi * (b(2,i,j,k) - b(2,i,j-1,k))
#endif /* NDIMS == 3 */
        end do
      end do
    end do
!
  end subroutine divergence
#endif /* MHD */
!
!===============================================================================
!
end module
