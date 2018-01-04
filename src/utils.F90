!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2015-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: UTILS
!!
!!  This module provides utility subroutines like assertion.
!!
!!*****************************************************************************
!
module utils

! module variables are not implicit by default
!
  implicit none

! subroutine interfaces
!
  interface assert
    module procedure assert_scalar
    module procedure assert_array_1d
    module procedure assert_array_2d
    module procedure assert_array_3d
  end interface

! by default everything is public
!
  public :: assert

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! ASSERTION SUBROUTINES
!
!===============================================================================
!
! subroutine ASSERT_SCALAR:
! ------------------------
!
!   Subroutine check assertion of the scalar condition.
!
!   Arguments:
!
!     cond - the condition;
!     tag  - the condition's tag;
!
!===============================================================================
!
  subroutine assert_scalar(cond, tag)

! import procedures and variables from other modules
!
    use error          , only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical         , intent(in) :: cond
    character(len=*), intent(in) :: tag
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! check if the condition was asserted
!
    if (.not. cond) then
      call print_error(tag, 'failed scalar assertion!')
      stop 'program terminated by assert_scalar()'
    end if
#endif /* DEBUG */

!-------------------------------------------------------------------------------
!
  end subroutine assert_scalar
!
!===============================================================================
!
! subroutine ASSERT_ARRAY_1D:
! --------------------------
!
!   Subroutine check assertion of the vector condition.
!
!   Arguments:
!
!     cond - the vector of conditions;
!     tag  - the condition's tag;
!
!===============================================================================
!
  subroutine assert_array_1d(cond, tag)

! import procedures and variables from other modules
!
    use error, only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical         , dimension(:), intent(in) :: cond
    character(len=*)              , intent(in) :: tag

! local variables
!
    integer :: i, im
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! check if the condition was asserted
!
    im = size(cond)
    do i = 1, im
      if (.not. cond(i)) then
        call print_error(tag, 'failed assertion!')
        stop 'program terminated by assert_array_1d()'
      end if
    end do ! i = 1, im
#endif /* DEBUG */

!-------------------------------------------------------------------------------
!
  end subroutine assert_array_1d
!
!===============================================================================
!
! subroutine ASSERT_ARRAY_2D:
! --------------------------
!
!   Subroutine check assertion of the 2D array condition.
!
!   Arguments:
!
!     cond - the array of conditions;
!     tag  - the condition's tag;
!
!===============================================================================
!
  subroutine assert_array_2d(cond, tag)

! import procedures and variables from other modules
!
    use error, only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical         , dimension(:,:), intent(in) :: cond
    character(len=*)                , intent(in) :: tag

! local variables
!
    integer :: i, j, im, jm
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! check if the condition was asserted
!
    im = size(cond, 1)
    jm = size(cond, 2)
    do j = 1, jm
      do i = 1, im
        if (.not. cond(i,j)) then
          call print_error(tag, 'failed assertion!')
          stop 'program terminated by assert_array_2d()'
        end if
      end do ! i = 1, im
    end do ! j = 1, jm
#endif /* DEBUG */

!-------------------------------------------------------------------------------
!
  end subroutine assert_array_2d
!
!===============================================================================
!
! subroutine ASSERT_ARRAY_3D:
! --------------------------
!
!   Subroutine check assertion of the 3D array condition.
!
!   Arguments:
!
!     cond - the array of conditions;
!     tag  - the condition's tag;
!
!===============================================================================
!
  subroutine assert_array_3d(cond, tag)

! import procedures and variables from other modules
!
    use error, only : print_error

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical         , dimension(:,:,:), intent(in) :: cond
    character(len=*)                  , intent(in) :: tag

! local variables
!
    integer :: i, j, k, im, jm, km
!
!-------------------------------------------------------------------------------
!
#ifdef DEBUG
! check if the condition was asserted
!
    im = size(cond, 1)
    jm = size(cond, 2)
    km = size(cond, 3)
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          if (.not. cond(i,j,k)) then
            call print_error(tag, 'failed assertion!')
            stop 'program terminated by assert_array_3d()'
          end if
        end do ! i = 1, im
      end do ! j = 1, jm
    end do ! k = 1, km
#endif /* DEBUG */

!-------------------------------------------------------------------------------
!
  end subroutine assert_array_3d

!===============================================================================
!
end module utils
