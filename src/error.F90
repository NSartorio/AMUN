!!******************************************************************************
!!
!! module: error - handling errors
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of AMUN.
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
    write(*,"('  [error in ', a, ']: ', a)") trim(position), trim(text)
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
    write(*,"('[warning in ', a, ']: ', a)") trim(position), trim(text)

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
    write(*,"('   [info in ', a, ']: ', a)") trim(position), trim(text)

!-------------------------------------------------------------------------------
!
  end subroutine print_info

!===============================================================================
!
end module
