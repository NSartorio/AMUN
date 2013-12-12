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

! pointers to the reconstruction and limiter procedures
!
  procedure(reconstruct), pointer, save :: reconstruct_states => null()
  procedure(limit_mm)   , pointer, save :: limit_derivatives  => null()

! module parameters
!
  real(kind=8), save :: eps        = epsilon(1.0d+00)

! flags for reconstruction corrections
!
  logical     , save :: positivity = .false.

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations, finalize_interpolations
  public :: reconstruct
  public :: fix_positivity
  public :: minmod, minmod3

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
    character(len=255) :: reconstruction = "tvd"
    character(len=255) :: limiter        = "mm"
    character(len=255) :: positivity_fix = "off"
    character(len=255) :: name_rec       = ""
    character(len=255) :: name_lim       = ""
!
!-------------------------------------------------------------------------------
!
! obtain the user defined interpolation methods and coefficients
!
    call get_parameter_string("reconstruction", reconstruction)
    call get_parameter_string("limiter"       , limiter       )
    call get_parameter_string("fix_positivity", positivity_fix)
    call get_parameter_real  ("eps"           , eps           )

! select the reconstruction method
!
    select case(trim(reconstruction))
    case ("tvd", "TVD")
      name_rec           =  "2nd order TVD"
      reconstruct_states => reconstruct_tvd
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected reconstruction method is not " //     &
                           "implemented: " // trim(reconstruction)
        stop
      end if
    end select

! select the limiter
!
    select case(trim(limiter))
    case ("mm")
      name_lim           =  "MinMod"
      limit_derivatives  => limit_mm
    case ("mc")
      name_lim           =  "McCormic"
      limit_derivatives  => limit_mc
    case ("lf")
      name_lim           =  "Lax-Friedrichs"
      limit_derivatives  => limit_lf
    case default
      if (verbose) then
        write (*,"(1x,a)") "The selected limiter is not implemented: " //      &
                           trim(limiter)
        stop
      end if
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

      write (*,"(4x,a,1x,a)"    ) "reconstruction         =", trim(name_rec)
      write (*,"(4x,a,1x,a)"    ) "limiter                =", trim(name_lim)
      write (*,"(4x,a,1x,a)"    ) "fix positivity         =", trim(positivity_fix)

    end if

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
! release the procedure pointers
!
    nullify(reconstruct_states)
    nullify(limit_derivatives)

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
! reconstruct the states using the selected subroutine
!
    call reconstruct_states(n, h, f(:), fl(:), fr(:))

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
    integer                    :: i, im1
    real(kind=8), dimension(n) :: df
!
!-------------------------------------------------------------------------------
!
! calculate limited derivative
!
    call limit_derivatives(n, f(:), df(:))

! calculate the left- and right-side interface interpolations
!
    do i = 2, n

! update the left and right-side interpolation states
!
      fl(i  ) = f(i) + df(i)
      fr(i-1) = f(i) - df(i)

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
! subroutine LIMIT_MM:
! -------------------
!
!   Subroutine calculates the local derivative applying the minmod TVD limiter.
!
!===============================================================================
!
  subroutine limit_mm(n, f, df)

    implicit none

! input/output arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: df

! local variables
!
    integer                    :: i, ip1
    real(kind=8)               :: ds
    real(kind=8), dimension(n) :: dfl, dfr
!
!------------------------------------------------------------------------------
!
! calculate the left- and right-side derivatives
!
    do i = 1, n - 1
      ip1      = i + 1

      dfr(i  ) = f(ip1) - f(i)
      dfl(ip1) = dfr(i)
    end do
    dfl(1) = 0.0d+00
    dfr(n) = 0.0d+00

! second order interpolation
!
    do i = 1, n

      df(i) = (sign(2.5d-01, dfl(i)) + sign(2.5d-01, dfr(i)))                  &
                                               * min(abs(dfl(i)), abs(dfr(i)))

    end do ! i = 1, n

!-------------------------------------------------------------------------------
!
  end subroutine limit_mm
!
!===============================================================================
!
! subroutine LIMIT_MC:
! -------------------
!
!   Subroutine calculates the local derivative applying the McCormic TVD
! limiter.
!
!===============================================================================
!
  subroutine limit_mc(n, f, df)

    implicit none

! input/output arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: df

! local variables
!
    integer                    :: i, ip1
    real(kind=8)               :: ds
    real(kind=8), dimension(n) :: dfl, dfr
!
!------------------------------------------------------------------------------
!
! calculate the left- and right-side derivatives
!
    do i = 1, n - 1
      ip1      = i + 1

      dfr(i  ) = f(ip1) - f(i)
      dfl(ip1) = dfr(i)
    end do
    dfl(1) = 0.0d+00
    dfr(n) = 0.0d+00

! second order interpolation
!
    do i = 1, n

      df(i) = (sign(5.0d-01, dfr(i)) + sign(5.0d-01, dfl(i)))                  &
                                      * min(abs(dfr(i)), abs(dfl(i))           &
                                             , 2.5d-01 * abs(dfr(i) + dfl(i)))

    end do ! i = 1, n

!-------------------------------------------------------------------------------
!
  end subroutine limit_mc
!
!===============================================================================
!
! subroutine LIMIT_LF:
! -------------------
!
!   Subroutine calculates the local derivative applying the Lax-Friendrich TVD
!   limiter.
!
!===============================================================================
!
  subroutine limit_lf(n, f, df)

    implicit none

! input/output arguments
!
    integer                   , intent(in)  :: n
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: df

! local variables
!
    integer                    :: i, ip1
    real(kind=8)               :: ds
    real(kind=8), dimension(n) :: dfl, dfr
!
!------------------------------------------------------------------------------
!
! calculate the left- and right-side derivatives
!
    do i = 1, n - 1
      ip1      = i + 1

      dfr(i  ) = f(ip1) - f(i)
      dfl(ip1) = dfr(i)
    end do
    dfl(1) = 0.0d+00
    dfr(n) = 0.0d+00

! second order interpolation
!
    do i = 1, n

! calculate the monotonicity indicator
!
      ds = dfr(i) * dfl(i)

! if the region is monotonic calculate limited derivative, otherwise set zero
!
      if (ds > 0.0d+00) then

! use selected limiter
!
        df(i) = ds / (dfr(i) + dfl(i))
      else
        df(i) = 0.0d+00
      end if

    end do ! i = 1, n

!-------------------------------------------------------------------------------
!
  end subroutine limit_lf
!
!===============================================================================
!
! function MINMOD:
! ---------------
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
! subroutine MINMOD3:
! ------------------
!
!   Function returns the minimum module value among two arguments and their
!   average.
!
!   Arguments:
!
!     a, b - two real values;
!
!
!===============================================================================
!
  real(kind=8) function minmod3(a, b)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b
!
!-------------------------------------------------------------------------------
!
! calculate the minimal module value
!
    minmod3 = (sign(1.0d+00, a) + sign(1.0d+00, b))                            &
                                   * min(abs(a), abs(b), 2.5d-01 * abs(a + b))

! return the value
!
    return

!-------------------------------------------------------------------------------
!
  end function minmod3
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

!-------------------------------------------------------------------------------
!
  end subroutine fix_positivity

!===============================================================================
!
end module interpolations
