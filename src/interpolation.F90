!!*****************************************************************************
!!
!! module: interpolation - subroutines for different kinds of interpolation
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
module interpolation

  implicit none

  contains
!
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, vx, vl, vr)

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real, dimension(n), intent(in)  :: vx
    real, dimension(n), intent(out) :: vl, vr

! local variables
!
    integer            :: i
    real               :: dv
    real, dimension(n) :: ds, dvl, dvr
!
!------------------------------------------------------------------------------
!
! second order interpolation
!
    dvl(1) = 0.0
    dvr(n) = 0.0

    do i = 1, n-1
      dvr(i  ) = vx(i+1) - vx(i)
      dvl(i+1) = dvr(i)
    enddo

    do i = 1, n
      ds (i) = dvr(i) * dvl(i)

      if (ds(i) .gt. 0.0) then
        dv  = ds(i) / (dvr(i) + dvl(i))

        vl(i) = vx(i) + dv
        vr(i) = vx(i) - dv
      else
        vl(i) = vx(i)
        vr(i) = vx(i)
      endif
    enddo

    do i = 1, n-1
      vr(i) = vr(i+1)
    enddo
    vr(n) = vx(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct

!===============================================================================
!
end module
