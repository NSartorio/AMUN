!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2017 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: CONSTANTS
!!
!!  This module provides mathematical, physical, and unit conversion
!!  constants and factors.
!!
!!*****************************************************************************
!
module constants

! module variables are not implicit by default
!
  implicit none

! π and its multiplications
!
  real(kind=8), parameter :: pi   = 3.14159265358979323846264338327950d+00
  real(kind=8), parameter :: pi2  = 6.28318530717958647692528676655900d+00
  real(kind=8), parameter :: pi4  = 1.25663706143591729538505735331180d+01

! exp(-1/2)
!
  real(kind=8), parameter :: ehi  = 0.60653065971263342360379953499118d+00

! conversion between angular units (degree to radian and radian to degree)
!
  real(kind=8), parameter :: d2r  = 1.74532925199432954743716805978693d-02
  real(kind=8), parameter :: r2d  = 5.72957795130823228646477218717337d+01

! by default everything is public
!
  public

!===============================================================================
!
end module constants
