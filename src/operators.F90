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
!!*****************************************************************************
!!
!! module: OPERATORS
!!
!!  This module provides differential operators like gradient, divergence, or
!!  curl.
!!
!!*****************************************************************************
!
module operators

! module variables are not implicit by default
!
  implicit none

! by default everything is public
!
  private

! declare public subroutines
!
  public :: initialize_operators, finalize_operators

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
!===============================================================================
!
! subroutine INITIALIZE_OPERATORS:
! -------------------------------
!
!   Subroutine initializes the module structures, pointers and variables.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_operators(verbose, iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!

!-------------------------------------------------------------------------------
!
  end subroutine initialize_operators
!
!===============================================================================
!
! subroutine FINALIZE_OPERATORS:
! -----------------------------
!
!   Subroutine releases the memory used by module variables and pointers.
!
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_operators(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout)    :: iret
!
!-------------------------------------------------------------------------------
!

!-------------------------------------------------------------------------------
!
  end subroutine finalize_operators

!===============================================================================
!
end module operators
