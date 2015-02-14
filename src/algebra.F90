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
!! module: ALGEBRA
!!
!!  This module provides all sort of algebra subroutines.
!!
!!******************************************************************************
!
module algebra

! module variables are not implicit by default
!
  implicit none

! by default everything is private
!
  private

! declare public subroutines
!
  public :: quadratic, quadratic_normalized

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! function QUADRATIC:
! ------------------
!
!   Function finds real roots of a quadratic equation:
!
!     f(x) = a₂ x² + a₁ x + a₀ = 0
!
!   Remark:
!
!     The coefficients are renormalized in order to avoid overflow and
!     unnecessary underflow.  The catastrophic cancellation is avoided as well.
!
!   Arguments:
!
!     a₂, a₁, a₀ - the quadratic equation coefficients;
!     x          - the root array;
!
!   Return value:
!
!     The number of roots found.
!
!===============================================================================
!
  function quadratic(a, x) result(nr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(3), intent(in)  :: a
    real(kind=8), dimension(2), intent(out) :: x
    integer                                 :: nr

! local variables
!
    real(kind=8), dimension(3) :: b
    real(kind=8)               :: bh, dl, dr, tm
!
!-------------------------------------------------------------------------------
!
! in order to avoid overflow or needless underflow, divide all coefficients by
! the maximum among them, but only if the maximum value is really large
!
    b(1:3) = a(1:3) / maxval(abs(a(1:3)))

! check if the coefficient a₂ is not zero
!
    if (b(3) /= 0.0d+00) then

! check if the coefficient a₀ is not zero
!
      if (b(1) /= 0.0d+00) then

! coefficients a₂ and a₀ are nonzero, so we solve the quadratic equation
!
!   a₂ x² + a₁ x + a₀ = 0
!
! calculate half of a₁
!
        bh = 0.5d+00 * b(2)

! calculate Δ = a₁² - 4 a₂ a₀
!
        dl = bh * bh - b(3) * b(1)

! check if Δ is larger then zero
!
        if (dl > 0.0d+00) then

! calculate quare root of Δ
!
          dr = sqrt(dl)

! Δ > 0, so the quadratic has two real roots
!
          if (b(2) /= 0.0d+00) then
            tm   = - (bh + sign(dr, bh))
            x(1) =   tm / b(3)
            x(2) = b(1) / tm
          else
            x(2) = dr / b(3)
            x(1) = - x(2)
          end if

! update the number of roots
!
          nr   = 2

! sort roots
!
          if (x(1) > x(2)) then
            tm   = x(1)
            x(1) = x(2)
            x(2) = tm
          end if

        else if (dl == 0.0d+00) then ! Δ = 0

! Δ = 0, so the quadratic has two identical real roots
!
          x(1) = - 0.5d+00 * b(2) / b(3)
          x(2) = x(1)

! update the number of roots
!
          nr   = 2

        else ! Δ < 0

! Δ < 0, so the quadratic does not have any real roots
!
          x(1) = 0.0d+00
          x(2) = 0.0d+00

! update the number of roots
!
          nr   = 0

        end if ! Δ < 0

      else ! a₀ = 0

! since a₀ is zero, we have one zero root and the second root is calculated from
! the linear formula
!
!   (a₂ x + a₁) x = 0
!
        tm   = - b(2) / b(3)
        if (tm < 0.0d+00) then
          x(1) = tm
          x(2) = 0.0d+00
        else
          x(1) = 0.0d+00
          x(2) = tm
        end if

! update the number of roots
!
        nr   = 2

      end if ! a₀ = 0

    else ! a₂ = 0

! since coefficient a₂ is zero, the quadratic equation is reduced to a linear
! one; one root is undetermined in this case
!
!   a₁ x + a₀ = 0
!
      x(2) = 0.0d+00

! find the unique root
!
      if (b(2) /= 0.0d+00) then

! the coefficient a₁ is not zero, so the quadratic has only one root
!
        x(1) = - b(1) / b(2)

! update the number of roots
!
        nr   = 1

      else ! a₁ = 0

! the coefficient a₁ is zero, so the quadratic has no roots
!
!   a₀ = 0
!
        x(1) = 0.0d+00

! update the number of roots
!
        nr   = 0

      end if ! a₁ = 0

    end if ! a₂ = 0

!-------------------------------------------------------------------------------
!
  end function quadratic
!
!===============================================================================
!
! function QUADRATIC_NORMALIZED:
! -----------------------------
!
!   Function finds real roots of the normalized quadratic equation:
!
!     f(x) = x² + a₁ x + a₀ = 0
!
!   Remark:
!
!     The coefficients are renormalized in order to avoid overflow and
!     unnecessary underflow.  The catastrophic cancellation is avoided as well.
!
!   Arguments:
!
!     a₁, a₀ - the quadratic equation coefficients;
!     x     - the root array;
!
!   Return value:
!
!     The number of roots found.
!
!===============================================================================
!
  function quadratic_normalized(a, x) result(nr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(2), intent(in)  :: a
    real(kind=8), dimension(2), intent(out) :: x
    integer                                 :: nr

! local variables
!
    real(kind=8) :: bh, dl, dr, tm
!
!-------------------------------------------------------------------------------
!
! check if the coefficient a₀ is not zero
!
    if (a(1) /= 0.0d+00) then

! coefficients a₀ is nonzero, so we solve the quadratic equation
!
!   x² + a₁ x + a₀ = 0
!
! calculate half of a₁
!
      bh = 0.5d+00 * a(2)

! calculate Δ = ¼ a₁² - a₀
!
      dl = bh * bh - a(1)

! check if Δ is larger then zero
!
      if (dl > 0.0d+00) then

! calculate quare root of Δ
!
        dr = sqrt(dl)

! Δ > 0, so the quadratic has two real roots, x₁ = - ½ a₁ ± √Δ
!
        if (a(2) > 0.0d+00) then
          tm   = bh + dr
          x(1) = - tm
          x(2) = - a(1) / tm
        else if (a(2) < 0.0d+00) then
          tm   = bh - dr
          x(1) = - a(1) / tm
          x(2) = - tm
        else
          x(1) = - dr
          x(2) =   dr
        end if

! update the number of roots
!
        nr   = 2

      else if (dl == 0.0d+00) then

! Δ = 0, so the quadratic has two identical real roots
!
        x(:) = - bh

! update the number of roots
!
        nr   = 2

      else

! Δ < 0, so the quadratic does not have any real root
!
        x(1) = 0.0d+00
        x(2) = 0.0d+00

! update the number of roots
!
        nr   = 0

      end if

    else

! since a₀ is zero, we have one zero root and the second root is calculated from
! the linear formula
!
!   (x + a₁) x = 0
!
        x(1) = 0.0d+00
        x(2) = - a(2)

! update the number of roots
!
        nr   = 2

    end if ! a₀ = 0

!-------------------------------------------------------------------------------
!
  end function quadratic_normalized

!===============================================================================
!
end module algebra
