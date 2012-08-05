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
!! module: INTERPOLATIONS
!!
!!  This module provides subroutines to interpolate variables and reconstruct
!!  the Riemann states.
!!
!!
!!******************************************************************************
!
module interpolations

! module variables are not implicit by default
!
  implicit none

! module parameters
!
  real, save :: kappa = 1.0d0
  real, save :: rad   = 1.0d0
  real, save :: eps   = epsilon(rad)

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations, reconstruct, minmod

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
  subroutine initialize_interpolations()

! include external procedures
!
    use parameters, only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! local variables
!
    real :: cfl = 0.5d0
!
!-------------------------------------------------------------------------------
!
! obtain the interpolation coefficients
!
    call get_parameter_real("rad", rad)
    call get_parameter_real("eps", eps)
    call get_parameter_real("cfl", cfl)

! calculate κ = (1 - ν) / ν
!
    kappa = (1.0d0 - cfl) / cfl

!-------------------------------------------------------------------------------
!
  end subroutine initialize_interpolations
#ifdef TVD
!
!===============================================================================
!
! subroutine RECONSTRUCT:
! ----------------------
!
!   Subroutine reconstructs the interface states using the second order TVD
!   method with a limiter selected through a compilation flag LIMITER.
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

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i, im1, ip1
    real               :: df, ds, dm, dp
!
!-------------------------------------------------------------------------------
!
! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n

! prepare indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! calculate derivatives
!
      dm = f(i  ) - f(im1)
      dp = f(ip1) - f(i  )

! obtain the limited derivative
!
#ifdef MINMOD
      df = 0.5d0 * minmod(dm, dp)

! interpolate the states
!
      fl(i  ) = f(i) + df
      fr(im1) = f(i) - df
#endif /* MINMOD */

#ifdef LF
! obtain the sign change detector
!
      ds = dm * dp

! depending on the sign change choose the proper slope
!
      if (ds > 0.0d0) then

! calculate derivative
!
        df  = ds / (dm + dp)

! interpolate the states
!
        fl(i  ) = f(i) + df
        fr(im1) = f(i) - df

      else

! copy the states
!
        fl(i  ) = f(i)
        fr(im1) = f(i)

      end if
#endif /* LF */

    end do

! prepare the last point
!
    fl(1) = f(1)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
#endif /* TVD */
#ifdef WENO3
!
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation;
!              the subroutine implements the improved 3rd order WENO scheme as
!              described in Yamaleev & Carpenter, 2009, Journal of Computational
!              Physics, 228, 3025
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

! local variables are not implicit by default
!
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
    h2    = h * h

! prepare initial left derivative
!
    dfp   = f(2) - f(1)
    fp    = f(1)
    fr(n) = f(n)

!! third order WENO interpolation
!!
! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n

! prepare indices
!
      im1     = max(1, i - 1)
      ip1     = min(n, i + 1)

! calculate left and right derivatives
!
      dfm     = dfp
      dfp     = f(ip1) - f(i  )
      ww      = (dfp - dfm)**2

! calculate corresponding betas
!
      bp      = dfp * dfp
      bm      = dfm * dfm

! calculate improved alphas
!
      ap      = 1.0d0 + ww / (bp + h2)
      am      = 1.0d0 + ww / (bm + h2)

! calculate weights
!
      wp      = dp * am
      wm      = dm * ap
      ww      = 2.0d0 * (wp + wm)

! calculate right and left sides interpolations
!
      fm      = - f(ip1) + 3.0d0 * f(i  )

! calculate the left state
!
      fr(im1) = (wp * fp + wm * fm) / ww

! calculate weights
!
      wp      = dp * ap
      wm      = dm * am
      ww      = 2.0d0 * (wp + wm)

! calculate right and left sides interpolations
!
      fp      =   f(i  ) +         f(ip1)
      fm      = - f(im1) + 3.0d0 * f(i  )

! calculate the left state
!
      fl(i  ) = (wp * fp + wm * fm) / ww

    end do

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* WENO3 */
#ifdef LIMO3
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

! local variables are not implicit by default
!
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

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
#endif /* LIMO3 */
!
!===============================================================================
!
! subroutine MINMOD:
! -----------------
!
!   Function returns the minimum module value among two arguments.
!
!   Arguments:
!
!     a, b - two real values;
!
!
!===============================================================================
!
  real function minmod(a, b)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real, intent(in) :: a, b
!
!-------------------------------------------------------------------------------
!
! calculate the minimal module value
!
    minmod = (sign(0.5d0, a) + sign(0.5d0, b)) * min(abs(a), abs(b))

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function minmod
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
!   Arguments:
!
!     n  - the length of vectors f, fl, and fr;
!     f  - the cell centered vector representation;
!     fl - the left state reconstruction;
!     fr - the right state reconstruction;
!
!
!===============================================================================
!
  subroutine fix_positivity(n, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer           , intent(in)    :: n
    real, dimension(n), intent(in)    :: f
    real, dimension(n), intent(inout) :: fl, fr

! local variables
!
    integer :: i, im1, ip1
    real    :: fmn, fmx
!
!------------------------------------------------------------------------------
!
! look for negative values in the states along the vector
!
    do i = 1, n

! check if there is a negative value in the left or right states
!
      if (min(fl(i), fr(i)) .le. 0.0d0) then

! calculate the left and right indices
!
        im1 = max(1, i - 1)
        ip1 = min(n, i + 1)

! prepare the allowed interval
!
        fmn = min(f(i), f(ip1))
        fmx = max(f(i), f(ip1))

! limit the states using the zeroth-order reconstruction
!
        if ((fl(i) .lt. fmn) .or. (fl(i) .gt. fmx)) then
          fl(i  ) = f(i)
          fr(im1) = f(i)
        end if

        if ((fr(i) .lt. fmn) .or. (fr(i) .gt. fmx)) then
          fl(ip1) = f(ip1)
          fr(i  ) = f(ip1)
        end if
      end if
    end do

!-------------------------------------------------------------------------------
!
  end subroutine fix_positivity
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

    use coordinates, only : im, jm, km
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

!===============================================================================
!
end module interpolations
