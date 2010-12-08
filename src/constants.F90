!!*****************************************************************************
!!
!! module: constants - mathematical and physical constants
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!*****************************************************************************
!!
!!  This file is part of Godunov-AMR.
!!
!!  Godunov-AMR is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  Godunov-AMR is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!*****************************************************************************
!!
!
module constants

  implicit none

! number Pi and its multiplications
!
  real, parameter ::  pi =  3.1415926535897931159979634685442d0
  real, parameter :: dpi =  6.2831853071795862319959269370884d0
  real, parameter :: qpi = 12.5663706143591724639918538741770d0

end module constants