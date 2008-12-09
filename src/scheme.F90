!!*****************************************************************************
!!
!! module: scheme - handling the actual solver of the set of equations
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
module scheme

  implicit none

  contains
!
!======================================================================
!
! update: subroutine sweeps over all directions and integrated
!         the increment to the solution
!
!======================================================================
!
  subroutine update(u, du)

    use blocks, only : nvars
    use config, only : igrids, jgrids, kgrids

    implicit none

! input arguments
!
    real, dimension(nvars,igrids,jgrids,kgrids), intent(in)  :: u
    real, dimension(nvars,igrids,jgrids,kgrids), intent(out) :: du
!
!----------------------------------------------------------------------
!
    du(:,:,:,:) = 0.0

  end subroutine update

!======================================================================
!
end module
