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
!! module: ERROR - handling errors
!!
!!******************************************************************************
!
module error

  implicit none

  contains
!
!===============================================================================
!
! print_error: subroutine prints error
!
!===============================================================================
!
  subroutine print_error(position, text)

    implicit none

! input arguments
!
    character(len=*), intent(in) :: position, text
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[error in ', a, ']: ', a)") trim(position), trim(text)
    write(*,*)
    stop

!-------------------------------------------------------------------------------
!
  end subroutine print_error
!
!===============================================================================
!
! print_warning: subroutine prints warning
!
!===============================================================================
!
  subroutine print_warning(position, text)

    implicit none

! input arguments
!
    character(len=*), intent(in) :: position, text
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[warning in ', a, ']: ', a)") trim(position), trim(text)
    write(*,*)

!-------------------------------------------------------------------------------
!
  end subroutine print_warning
!
!===============================================================================
!
! print_info: subroutine prints information
!
!===============================================================================
!
  subroutine print_info(position, text)

    implicit none

! input arguments
!
    character(len=*), intent(in) :: position, text
!
!-------------------------------------------------------------------------------
!
    write(*,*)
    write(*,"('[info in ', a, ']: ', a)") trim(position), trim(text)
    write(*,*)

!-------------------------------------------------------------------------------
!
  end subroutine print_info

!===============================================================================
!
end module
