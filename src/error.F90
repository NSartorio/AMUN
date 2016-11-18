!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2016 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: ERROR
!!
!!  This module provides subroutines to print errors, warnings and
!!  notifications in a unified format.
!!
!!******************************************************************************
!
module error

! module variables are not implicit by default
!
  implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine PRINT_ERROR:
! ----------------------
!
!   Subroutine prints an error message with the module/subroutine where
!   the error occurred.
!
!   Arguments:
!
!     loc - string informing about the place where the error occurred
!           (for example 'module name::subroutine name:[line]');
!     msg - string of the actual error message;
!
!===============================================================================
!
  subroutine print_error(loc, msg)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in) :: loc, msg
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[ERROR in ', a, ']: ', a)") trim(loc), trim(msg)

!-------------------------------------------------------------------------------
!
  end subroutine print_error
!
!===============================================================================
!
! subroutine PRINT_WARNING:
! ------------------------
!
!   Subroutine prints a warning message with the module/subroutine where
!   the warning occurred.
!
!   Arguments:
!
!     loc - string informing about the place where the warning occurred
!           (for example 'module name::subroutine name:[line]');
!     msg - string of the actual warning message;
!
!===============================================================================
!
  subroutine print_warning(loc, msg)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in) :: loc, msg
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[WARNING in ', a, ']: ', a)") trim(loc), trim(msg)

!-------------------------------------------------------------------------------
!
  end subroutine print_warning
!
!===============================================================================
!
! subroutine PRINT_NOTIFICATION:
! -----------------------------
!
!   Subroutine prints a notification message with the module/subroutine where
!   the notification occurred.
!
!   Arguments:
!
!     loc - string informing about the place where the notification occurred
!           (for example 'module name::subroutine name:[line]');
!     msg - string of the actual notification message;
!
!===============================================================================
!
  subroutine print_notification(loc, msg)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in) :: loc, msg
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[NOTIFICATION in ', a, ']: ', a)") trim(loc), trim(msg)

!-------------------------------------------------------------------------------
!
  end subroutine print_notification

!===============================================================================
!
end module
