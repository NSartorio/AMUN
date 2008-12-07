!!*****************************************************************************
!!
!! program: Godunov-AMR
!!
!! Copyright (C) 2008 Grzegorz Kowal <kowal@astro.wisc.edu>
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
program godunov

! modules
!
  use config, only : read_config
  use io    , only : write_data
  use mesh  , only : init_mesh, clear_mesh
!
!----------------------------------------------------------------------
!
! local variables
!
!
!----------------------------------------------------------------------
!
! print info message
!
  write (*,"(1x,78('-'))")
  write (*,"(1x,18('='),4x,a,4x,19('='))") '      Godunov-AMR algorithm      '
  write (*,"(1x,18('='),4x,a,4x,19('='))") 'Copyright (C) 2008 Grzegorz Kowal'
  write (*,"(1x,78('-'))")
  write (*,*)

! read configuration file
!
  call read_config

! initialize our adaptive mesh, refine that mesh to the desired level
! according to the initialized problem
!
  call init_mesh

! write down the initial state
!
  call write_data('r', 0, 0)

! TODO: main loop, perform one step evolution of the system, do refinement/derefinement
! TODO: get new time step, dump data, print info about the progress

! deallocate and reset mesh
!
  call clear_mesh

!----------------------------------------------------------------------
!
end program
