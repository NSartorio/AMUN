!!*****************************************************************************
!!
!! module: evolution - handling the time evolution of the block structure
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
module evolution

  implicit none

  integer, save :: n
  real   , save :: t, dt

  contains
!
!======================================================================
!
! evolve: subroutine sweeps over all leaf blocks and performes one step
!         time evolution according to the selected integration scheme
!
!======================================================================
!
  subroutine evolve

    use blocks, only : block, plist

    implicit none

! local variables
!
    type(block), pointer :: pblock
!
!----------------------------------------------------------------------
!
! iterate over all blocks and perform one step of time evolution
!
      pblock => plist
      do while (associated(pblock))

! check if this block is a leaf
!
        if (pblock%leaf .eq. 'T') &
          call evolve_rk2(pblock)

! assign pointer to the next block
!
        pblock => pblock%next

      end do

! TODO: boundary conditions
!

!----------------------------------------------------------------------
!
  end subroutine evolve
!
!======================================================================
!
! evolve_rk2: subroutine evolves the current block using RK2 method
!
!======================================================================
!
  subroutine evolve_rk2(pblock)

    use blocks, only : block

    implicit none

! input arguments
!
    type(block), pointer, intent(in) :: pblock

! local variables
!
!
!----------------------------------------------------------------------
!

!----------------------------------------------------------------------
!
  end subroutine evolve_rk2

!======================================================================
!
end module
