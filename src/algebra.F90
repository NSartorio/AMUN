!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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

! include external procedures
!
  use iso_fortran_env, only : real32, real64, real128

! module variables are not implicit by default
!
  implicit none

! maximum real kind
!
  integer, parameter :: max_real_kind = max(real32, real64, real128)

! by default everything is private
!
  private

! declare public subroutines
!
  public :: max_real_kind
  public :: quadratic, quadratic_normalized
  public :: cubic, cubic_normalized
  public :: quartic
  public :: tridiag, invert

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
!     x          - the root array; if there are two roots, x(1) corresponds to
!                  the one with '-' and x(2) to the one with '+';
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
          if (b(2) > 0.0d+00) then
            tm   = - bh - dr
            x(1) =   tm / b(3)
            x(2) = b(1) / tm
          else if (b(2) < 0.0d+00) then
            tm   = - bh + dr
            x(1) = b(1) / tm
            x(2) =   tm / b(3)
          else
            x(2) = dr / b(3)
            x(1) = - x(2)
          end if

! update the number of roots
!
          nr   = 2

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
        if (tm >= 0.0d+00) then
          x(1) = 0.0d+00
          x(2) = tm
        else
          x(1) = tm
          x(2) = 0.0d+00
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
          x(1) = - a(1) / tm
          x(2) = - tm
        else if (a(2) < 0.0d+00) then
          tm   = bh - dr
          x(1) = - tm
          x(2) = - a(1) / tm
        else
          x(1) =   dr
          x(2) = - dr
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
!
!===============================================================================
!
! function CUBIC:
! --------------
!
!   Function finds real roots of a cubic equation:
!
!     f(x) = a₃ x³ + a₂ x² + a₁ x + a₀ = 0
!
!   Remark:
!
!     The coefficients are renormalized in order to avoid overflow and
!     unnecessary underflow.
!
!   Arguments:
!
!     a(:) - the cubic equation coefficients vector (in increasing moments);
!     x(:) - the root vector;
!
!   Return value:
!
!     The number of roots found.
!
!===============================================================================
!
  function cubic(a, x) result(nr)

! include external constants
!
    use constants, only : pi2

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(4), intent(in)  :: a
    real(kind=8), dimension(3), intent(out) :: x
    integer                                 :: nr

! local variables
!
    real(kind=8)            :: ac, bc, qc, cc, dc, ec, rc, th, cs, a2, b3

! local parameters
!
    real(kind=8), parameter :: c0 = 1.0d+00 / 3.0d+00, eps = epsilon(1.0d+00)
!
!-------------------------------------------------------------------------------
!
! check if the coefficient a₃ is not zero
!
    if (a(4) /= 0.0d+00) then

! check if the coefficient a₀ is not zero
!
      if (a(1) /= 0.0d+00) then

! both coefficients a₃ and a₀ are not zero, so we have to solve full cubic
! equation
!
!   a₃ x³ + a₂ x² + a₁ x + a₀ = 0
!
! compute coefficient A = 2 a₂³ - 9 a₃ a₂ a₁ + 27 a₃² a₀
!
        ac = 2.0d+00 * a(3)**3 - 9.0d+00 * a(4) * a(3) * a(2)                  &
                                                   + 27.0d+00 * a(4)**2 * a(1)

! compute coefficient B = a₂² - 3 a₃ a₁
!
        bc = a(3)**2 - 3.0d+00 * a(4) * a(2)

! calculate A² and B³
!
        a2 = ac**2
        b3 = bc**3

! compute coefficient Q² = A² - 4 B³
!
        qc = a2 - 4.0d+00 * b3

! check the condition for any root existance
!
        if (qc > 0.0d+00) then

! the case Q > 0, in which only one root is real
!
! calculate C³ = ½ (Q + A)
!
          cc = 0.5d+00 * (sqrt(qc) + ac)
          dc = sign(1.0d+00, cc) * abs(cc)**c0
          if (cc == 0.0d+00) then
            ec = b3 / a2
            ec = 2.0d+00 * ec * (1.0d+00 / sqrt(1.0d+00 - 4.0d+00 * ec)        &
                                                               + 1.0d+00 / ac)
          else
            ec = b3 / cc
          end if
          ec = sign(1.0d+00, ec) * abs(ec)**c0

! calculate the unique real root
!
          x(1) = - c0 * (a(3) + dc + ec) / a(4)
          x(2) = 0.0d+00
          x(3) = 0.0d+00

! update the number of roots
!
          nr = 1

        else if (qc == 0.0d+00) then

          if (bc /= 0.0d+00) then

! the case Q = 0 and B ≠ 0, in which all three roots are real and one root
! is double
!
            x(1) = 0.5d+00 * (9.0d+00 * a(4) * a(1) - a(3) * a(2)) / bc
            x(2) = x(1)
            x(3) = (- 9.0d+00 * a(4) * a(1) + 4.0d+00 * a(3) * a(2)            &
                                                        - a(3)**3 / a(4)) / bc

! update the number of roots
!
            nr = 3

          else

! the case Q = 0 and B = 0, in which all three roots are real and equal
!
            x(:) = - c0 * a(3) / a(4)

! update the number of roots
!
            nr = 3

          end if

        else

! Q < 0, so there are three distinct and real roots
!
! compute Q and R
!
          th = 3.0d+00 * a(4)
          qc = bc / th**2
          rc = 0.5d+00 * ac / th**3

! compute cosine
!
          cs = sign(0.5d+00, ac * bc) * sqrt(ac**2 / abs(bc)**3)

! check condition
!
          if (cs <= 1.0d+00) then

! there are three real roots; prepare coefficients to calculate them
!
            th = acos(cs)
            ac = - 2.0d+00 * sqrt(qc)
            bc = - c0 * (a(3) / a(4))

! calculate roots
!
            x(1) = ac * cos(c0 *  th       ) + bc
            x(2) = ac * cos(c0 * (th + pi2)) + bc
            x(3) = ac * cos(c0 * (th - pi2)) + bc

! update the number of roots
!
            nr = 3

          else

! there is only one real root; prepare coefficients to calculate it
!
            ac = - sign(1.0d+00, rc) * (abs(rc) + sqrt(rc**2 - qc**3))**c0
            if (ac /= 0.0d+00) then
              bc = qc / ac
            else
              bc = 0.0d+00
           end if

! calculate the root
!
            x(1) = (ac + bc) - c0 * (a(3) / a(4))

! reset remaining roots
!
            x(2) = 0.0d+00
            x(3) = 0.0d+00

! update the number of roots
!
            nr = 1

          end if

        end if

      else

! since the coefficient a₀ is zero, there is one zero root, and the remaining
! roots are found from a quadratic equation
!
!   (a₃ x² + a₂ x + a₁) x = 0
!
        x(3) = 0.0d+00

! call function quadratic() to find the remaining roots
!
        nr = quadratic(a(2:4), x(1:2))

! increase the number of roots
!
        nr = nr + 1

      end if

    else

! since the coefficient a₃ is zero, the equation reduces to a quadratic one,
! and one root is undetermined
!
!   a₂ x² + a₁ x + a₀ = 0
!
      x(3) = 0.0d+00

! call function quadratic() to find the remaining roots
!
      nr = quadratic(a(1:3), x(1:2))

    end if

! sort roots in the increasing order
!
    call sort(nr, x(:))

!-------------------------------------------------------------------------------
!
  end function cubic
!
!===============================================================================
!
! function CUBIC_NORMALIZED:
! -------------------------
!
!   Function finds real roots of the normalized cubic equation:
!
!     f(x) = x³ + a₂ x² + a₁ x + a₀ = 0
!
!   Arguments:
!
!     a(:) - the cubic equation coefficients vector (in increasing moments);
!     x(:) - the root vector;
!
!   Return value:
!
!     The number of roots found.
!
!===============================================================================
!
  function cubic_normalized(a, x) result(nr)

! include external constants
!
    use constants, only : pi2

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(3), intent(in)  :: a
    real(kind=8), dimension(3), intent(out) :: x
    integer                                 :: nr

! local variables
!
    real(kind=8)            :: ac, bc, qc, a3, b2, aa, bb, sq, cs, ph

! local parameters
!
    real(kind=8), parameter :: c0 = 1.0d+00 / 3.0d+00, c1 = 1.0d+00 / 2.7d+01
    real(kind=8), parameter :: c2 = 1.0d+00 / 4.0d+00
!
!-------------------------------------------------------------------------------
!
! check if the coefficient a₀ is not zero
!
    if (a(1) /= 0.0d+00) then

! both coefficients a₃ and a₀ are not zero, so we have to solve full cubic
! equation
!
!   x³ + a₂ x² + a₁ x + a₀ = 0
!
! compute coefficient A = a₂² - 3 a₁
!
      ac = c0 * (a(3)**2 - 3.0d+00 * a(2))

! compute coefficient B = (2 a₂³ - 9 a₂ a₁ + 27 a₀) / 27
!
      bc = c1 * (2.0d+00 * a(3)**3 - 9.0d+00 * a(3) * a(2) + 2.7d+01 * a(1))

! calculate powers of A and B
!
      a3 = c1 * ac**3
      b2 = c2 * bc**2

! compute coefficient Q² = B² / 4 + A³ / 27
!
      qc = b2 - a3

! check the condition for any root existance
!
      if (qc > 0.0d+00) then

! the case Q > 0, in which only one root is real
!
! calculate A and B
!
        sq = sqrt(qc)
        ph = - 0.5d+00 * bc + sq
        aa = sign(1.0d+00, ph) * abs(ph)**c0
        ph = - 0.5d+00 * bc - sq
        bb = sign(1.0d+00, ph) * abs(ph)**c0

! calculate the unique real root
!
        x(1) = aa + bb - (c0 * a(3))
        x(2) = 0.0d+00
        x(3) = 0.0d+00

! update the number of roots
!
        nr = 1

      else if (qc < 0.0d+00) then

! calculate angle φ
!
        if (bc > 0.0d+00) then
          cs = - sqrt(b2 / a3)
        else if (bc < 0.0d+00) then
          cs =   sqrt(b2 / a3)
        else ! B = 0
          cs = 0.0d+00
        end if

! there are three real roots; prepare coefficients to calculate them
!
        ph = acos(cs)
        aa = 2.0d+00 * sqrt(c0 * ac)

! calculate all three roots
!
        x(1) = aa * cos(c0 * (ph - pi2))
        x(2) = aa * cos(c0 *  ph       )
        x(3) = aa * cos(c0 * (ph + pi2))

! return to the original variables
!
        x(:) = x(:) - (c0 * a(3))

! update the number of roots
!
        nr = 3

      else ! Q = 0

! the case Q = 0, for which all three roots are real and two roots are equal
!
        if (bc > 0.0d+00) then
          aa   = sqrt(c0 * ac)
          x(1) = - 2.0d+00 * aa
          x(2) =             aa
          x(3) =             aa
        else if (bc < 0.0d+00) then
          aa   = sqrt(c0 * ac)
          x(1) =   2.0d+00 * aa
          x(2) = -           aa
          x(3) = -           aa
        else
          x(:) = 0.0d+00
        end if

! return to the original variables
!
        x(:) = x(:) - (c0 * a(3))

! update the number of roots
!
        nr = 3

      end if ! Q = 0

    else ! a₀ = 0

! since the coefficient a₀ is zero, there is one zero root, and the remaining
! roots are found by solving the quadratic equation
!
!   (x² + a₂ x + a₁) x = 0
!
      x(3) = 0.0d+00

! call function quadratic() to find the remaining roots
!
      nr = quadratic_normalized(a(2:3), x(1:2))

! increase the number of roots
!
      nr = nr + 1

    end if ! a₀ = 0

! sort roots in the increasing order
!
    call sort(nr, x(:))

!-------------------------------------------------------------------------------
!
  end function cubic_normalized
!
!===============================================================================
!
! subroutine QUARTIC:
! ------------------
!
!   Subroutine finds real roots of a quartic equation:
!
!     f(x) = a₄ * x⁴ + a₃ * x³ + a₂ * x² + a₁ * x + a₀ = 0
!
!   Arguments:
!
!     a(:) - the quartic equation coefficients vector (in increasing moments);
!     x(:) - the root vector;
!
!   References:
!
!     Hacke, Jr. J. E., "A Simple Solution of the General Quartic",
!     The American Mathematical Monthly, 1941, vol. 48, no. 5, pp. 327-328
!
!===============================================================================
!
  function quartic(a, x) result(nr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    real(kind=8), dimension(5), intent(in)  :: a
    real(kind=8), dimension(4), intent(out) :: x
    integer                                 :: nr

! local variables
!
    integer                    :: n, m
    real(kind=8)               :: pc, qc, rc
    real(kind=8)               :: zm, k2, l2, kc, lc

! local arrays
!
    real(kind=8), dimension(4) :: b
    real(kind=8), dimension(2) :: y
    real(kind=8), dimension(3) :: c, z
!
!-------------------------------------------------------------------------------
!
! check if the coefficient a₄ is not zero
!
    if (a(5) /= 0.0d+00) then

! check if the coefficient a₀ is not zero
!
      if (a(1) /= 0.0d+00) then

! both coefficients a₄ and a₀ are not zero, so we have to solve full quartic
! equation
!
!   a₄ * x⁴ + a₃ * x³ + a₂ * x² + a₁ * x + a₀ = 0
!
! calculate coefficients of P, Q, and R
!
        pc = (- 3.0d+00 * a(4)**2 + 8.0d+00 * a(5) * a(3)) / (8.0d+00 * a(5)**2)
        qc = (a(4)**3 - 4.0d+00 * a(5) * a(4) * a(3)                           &
                             + 8.0d+00 * a(5)**2 * a(2)) / (8.0d+00 * a(5)**3)
        rc = (- 3.0d+00 * a(4)**4 + 16.0d+00 * a(5) * a(4)**2 * a(3)           &
              - 64.0d+00 * a(5)**2 * a(4) * a(2) + 256.0d+00 * a(5)**3 * a(1)) &
                                                       / (256.0d+00 * a(5)**4)

! prepare coefficients for the cubic equation
!
        b(4) =   8.0d+00
        b(3) = - 4.0d+00 * pc
        b(2) = - 8.0d+00 * rc
        b(1) =   4.0d+00 * pc * rc - qc**2

! find roots of the cubic equation
!
        nr = cubic(b(1:4), z(1:3))

! check if the cubic solver returned any roots
!
        if (nr > 0) then

! take the maximum root
!
          zm = z(1)

! calculate coefficients to find quartic roots
!
          k2 = 2.0d+00 * zm - pc
          l2 =      zm * zm - rc

! check if both coefficients are larger then zero
!
          if (k2 > 0.0d+00 .or. l2 > 0.0d+00) then

! prepare coefficients K and L
!
            if (k2 > l2) then
              kc = sqrt(k2)
              lc = 0.5d+00 * qc / kc
            else if (l2 > k2) then
              lc = sqrt(l2)
              kc = 0.5d+00 * qc / lc
            else
              kc = sqrt(k2)
              lc = sqrt(l2)
            end if

! prepare coefficients for the first quadratic equation
!
            c(3) = 1.0d+00
            c(2) = - kc
            c(1) = zm + lc

! find the first pair of roots
!
            m = quadratic(c(1:3), y(1:2))

! update the roots
!
            do n = 1, m
              x(n) = y(n) - 0.25d+00 * a(4) / a(5)
            end do

! update the number of roots
!
            nr = m

! prepare coefficients for the second quadratic equation
!
            c(3) = 1.0d+00
            c(2) = kc
            c(1) = zm - lc

! find the second pair of roots
!
            m = quadratic(c(1:3), y(1:2))

! update the roots
!
            do n = 1, m
              x(nr + n) = y(n) - 0.25d+00 * a(4) / a(5)
            end do

! update the number of roots
!
            nr = nr + m

! reset the remainings roots
!
            do n = nr + 1, 4
              x(n) = 0.0d+00
            end do

          else

! both K² and L² are negative, which signifies that the quartic equation has
! no real roots; reset them
!
            x(:) = 0.0d+00

! update the number of roots
!
            nr   = 0

          end if

        else

! weird, no cubic roots... it might mean that this is not real quartic equation
! or something is wrong; reset the roots
!
          x(:) = 0.0d+00

        end if

      else

! since the coefficient a₀ is zero, there is one zero root, and the remaining
! roots are found from a cubic equation
!
!   (a₄ * x³ + a₃ * x² + a₂ * x + a₁) * x = 0
!
        x(4) = 0.0d+00

! call function cubic() to find the remaining roots
!
        nr = cubic(a(2:5), x(1:3))

! increase the number of roots
!
        nr = nr + 1

      end if

    else

! since the coefficient a₄ is zero, the equation reduces to a cubic one, and
! one root is undetermined
!
!   a₃ * x³ + a₂ * x² + a₁ * x + a₀ = 0
!
      x(4) = 0.0d+00

! call function cubic() to find the remaining roots
!
      nr = cubic(a(1:4), x(1:3))

    end if

! sort roots in the increasing order
!
    call sort(nr, x(:))

!-------------------------------------------------------------------------------
!
  end function quartic
!
!===============================================================================
!
! subroutine SORT:
! ---------------
!
!   Subroutine sorts the input vector in the ascending order by straight
!   insertion.
!
!   Arguments:
!
!     n    - the number of elements to sort;
!     x(:) - the vector which needs to be sorted;
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!===============================================================================
!
  subroutine sort(n, x)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(:), intent(inout) :: x

! local variables
!
    integer      :: k, l
    real(kind=8) :: y
!
!-------------------------------------------------------------------------------
!
! return if only one element
!
    if (n < 2) return

! iterate over all elements
!
    do k = 2, n

! pick the element
!
      y = x(k)

! go back to find the right place
!
      do l = k - 1, 1, -1
        if (x(l) <= y) goto 10
        x(l+1) = x(l)
      end do

! reset the index
!
      l = 0

! insert the element
!
10    x(l+1) = y

    end do ! k = 2, n

!-------------------------------------------------------------------------------
!
  end subroutine sort
!
!===============================================================================
!
! subroutine TRIDIAG:
! ------------------
!
!   Subroutine solves the system of tridiagonal linear equations.
!
!   Arguments:
!
!     n    - the size of the matrix;
!     l(:) - the vector of the lower off-diagonal values;
!     d(:) - the vector of the diagonal values;
!     u(:) - the vector of the upper off-diagonal values;
!     r(:) - the vector of the right side values;
!     x(:) - the solution vector;
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!===============================================================================
!
  subroutine tridiag(n, l, d, u, r, x, iret)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: l, d, u, r
    real(kind=8), dimension(n), intent(out) :: x
    integer                   , intent(out) :: iret

! local variables
!
    integer                    :: i, j
    real(kind=8)               :: t
    real(kind=8), dimension(n) :: g

! local parameters
!
    character(len=*), parameter :: loc = 'ALGEBRA::tridiag()'
!
!-------------------------------------------------------------------------------
!
! decomposition and forward substitution
!
    t    = d(1)
    x(1) = r(1) / t
    do i = 2, n
      j  = i - 1
      g(i) = u(j) / t
      t    = d(i) - l(i) * g(i)
      if (t == 0.0d+00) then
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Could not find solution!"
        iret = 1
        return
      end if
      x(i) = (r(i) - l(i) * x(j)) / t
    end do

! backsubstitution
!
    do i = n - 1, 1, -1
      j    = i + 1
      x(i) = x(i) - g(j) * x(j)
    end do

! set return value to success
!
    iret = 0

!-------------------------------------------------------------------------------
!
  end subroutine tridiag
!
!===============================================================================
!
! subroutine INVERT:
! -----------------
!
!   Subroutine inverts the squared matrix.
!
!   Arguments:
!
!     n      - the size of the matrix;
!     m(:,:) - the matrix elements;
!     r(:,:) - the inverted matrix elements;
!     f      - the solution flag;
!
!   References:
!
!     [1] Press, W. H, Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.,
!         "Numerical Recipes in Fortran",
!         Cambridge University Press, Cambridge, 1992
!
!===============================================================================
!
  subroutine invert(n, m, r, f)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                                 , intent(in)  :: n
    real(kind=max_real_kind), dimension(n,n), intent(in)  :: m
    real(kind=max_real_kind), dimension(n,n), intent(out) :: r
    logical                                 , intent(out) :: f

! local variables
!
    logical                                    :: flag = .true.
    integer                                    :: i, j, k, l
    real(kind=max_real_kind)                   :: t
    real(kind=max_real_kind), dimension(n,2*n) :: g
!
!-------------------------------------------------------------------------------
!
! augment input matrix with an identity matrix
!
    do i = 1, n
      do j = 1, 2 * n
        if (j <= n) then
          g(i,j) = m(i,j)
        else if ((i+n) == j) then
          g(i,j) = 1.0d+00
        else
          g(i,j) = 0.0d+00
        end if
      end do
    end do

! reduce augmented matrix to upper traingular form
!
    do k = 1, n - 1
      if (g(k,k) == 0.0d+00) then
        flag = .false.
        do i = k + 1, n
          if (g(i,k) /= 0.0d+00) then
            do j = 1, 2 * n
              g(k,j) = g(k,j) + g(i,j)
            end do
            flag = .true.
            exit
          end if

          if (flag .eqv. .false.) then
            r(:,:) = 0.0d+00
            f      = .false.
            return
          end if
        end do
      end if
      do j = k + 1, n
        t = g(j,k) / g(k,k)
        do i = k, 2 * n
          g(j,i) = g(j,i) - t * g(k,i)
        end do
      end do
    end do

! test for invertibility
!
    do i = 1, n
      if (g(i,i) == 0.0d+00) then
        r(:,:) = 0.0d+00
        f = .false.
        return
      end if
    end do

! make diagonal elements as 1
!
    do i = 1, n
      t = g(i,i)
      do j = i , 2 * n
        g(i,j) = g(i,j) / t
      end do
    end do

! reduced right side half of augmented matrix to identity matrix
!
    do k = n - 1, 1, -1
      do i = 1, k
        t = g(i,k+1)
        do j = k, 2 * n
          g(i,j) = g(i,j) - g(k+1,j) * t
        end do
      end do
    end do

! store answer
!
    do i = 1, n
      do j = 1, n
        r(i,j) = g(i,j+n)
      end do
    end do

    f = .true.

!-------------------------------------------------------------------------------
!
  end subroutine invert

!===============================================================================
!
end module algebra
