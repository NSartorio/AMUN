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
  use blocks, only : init_blocks, clear_blocks
  use config, only : read_config
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

! setup adaptive mesh structure, allocate first block, initialize its mesh and variables
!

! initialize block structure
!
  call init_blocks

! fill the first block with initial conditions, check refinement, refine if necessary by creating more blocks, this should create the initial structure of domain
!

! write down the initial state
!

! finalize blocks structure
!
  call clear_blocks

!----------------------------------------------------------------------
!
end program
