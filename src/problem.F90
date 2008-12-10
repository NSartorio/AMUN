!!*****************************************************************************
!!
!! module: problem - handling the initial problem definition
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
module problem

  implicit none

  contains
!
!======================================================================
!
! init_problem: subroutine initializes the variables according to
!               the studied problem
!
!======================================================================
!
  subroutine init_problem(pblock)

    use blocks, only : block, idn, imx, imy, imz, ien
    use config, only : ncells, ngrids, igrids, jgrids, kgrids, nghost &
                     , dens, pres, gammam1i

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: r, dx, dy, dz, en, enamb

! local arrays
!
    real, dimension(:), allocatable :: x, y, z
!
!----------------------------------------------------------------------
!
! calculate parameters
!
    enamb = gammam1i * pres
    en    = 100.0 * enamb

! get dimensions
!
    dm(1) = igrids
    dm(2) = jgrids
    dm(3) = kgrids

! allocate coordinates
!
    allocate(x(dm(1)))
    allocate(y(dm(2)))
    allocate(z(dm(3)))

! calculate cell sizes
!
    dx = (pblock%xmax - pblock%xmin) / ncells
    dy = (pblock%ymax - pblock%ymin) / ncells
#if NDIMS == 3
    dz = (pblock%zmax - pblock%zmin) / ncells
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, dm(1))/) - nghost - 0.5) * dx + pblock%xmin
    y(:) = ((/(j, j = 1, dm(2))/) - nghost - 0.5) * dy + pblock%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, dm(3))/) - nghost - 0.5) * dz + pblock%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0
#endif /* NDIMS == 3 */

! set variables
!
    pblock%u(idn,:,:,:) = dens
    pblock%u(imx,:,:,:) = 0.0d0
    pblock%u(imy,:,:,:) = 0.0d0
    pblock%u(imz,:,:,:) = 0.0d0

! set initial pressure
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        do i = 1, dm(1)

          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

          if (r .le. 0.1) then
            pblock%u(ien,i,j,k) = en
          else
            pblock%u(ien,i,j,k) = enamb
          endif

        end do
      end do
    end do

! deallocate coordinates
!
    deallocate(x)
    deallocate(y)
    deallocate(z)

!----------------------------------------------------------------------
!
  end subroutine init_problem
!
!======================================================================
!
! check_ref: function checks refinement criterium and returns +1 if
!            the criterium fullfiled and block is selected for
!            refinement, 0 there is no need for refinement, and -1 if
!            block is selected for refinement
!
!======================================================================
!
  function check_ref(pblock)

    use blocks, only : block, ien
    use config, only : ncells, ngrids, igrids, jgrids, kgrids, nghost

! input arguments
!
    type(block), pointer, intent(in) :: pblock

! return variable
!
    integer(kind=1) :: check_ref

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: dpmax, dpdx, dpdy, dpdz
!
!----------------------------------------------------------------------
!
! get dimensions
!
    dm(1) = igrids
    dm(2) = jgrids
    dm(3) = kgrids

! check gradient of pressure
!
    dpmax = 0.0d0

    do k = 1, dm(3)
      do j = 2, dm(2)-1
        do i = 2, dm(1)-1
          dpdx = abs(pblock%u(ien,i+1,j,k) - 2 * pblock%u(ien,i,j,k) + pblock%u(ien,i-1,j,k)) / &
                 max(1.0e-8, (abs(pblock%u(ien,i+1,j,k) - pblock%u(ien,i,j,k)) + abs(pblock%u(ien,i,j,k) - pblock%u(ien,i-1,j,k)) + 1.0e-2 * (abs(pblock%u(ien,i+1,j,k)) - 2 * abs(pblock%u(ien,i,j,k)) + abs(pblock%u(ien,i-1,j,k)))))
          dpdy = abs(pblock%u(ien,i,j+1,k) - 2 * pblock%u(ien,i,j,k) + pblock%u(ien,i,j-1,k)) / &
                 max(1.e-8, (abs(pblock%u(ien,i,j+1,k) - pblock%u(ien,i,j,k)) + abs(pblock%u(ien,i,j,k) - pblock%u(ien,i,j-1,k)) + 1.0e-2 * (abs(pblock%u(ien,i,j+1,k)) - 2 * abs(pblock%u(ien,i,j,k)) + abs(pblock%u(ien,i,j-1,k)))))

          dpmax = max(dpmax, dpdx, dpdy)
        end do
      end do
    end do

! check condition
!
    check_ref = 0

    if (dpmax .gt. 0.7) then
      check_ref =  1
    endif
    if (dpmax .lt. 0.3) then
      check_ref = -1
    endif

    return

!----------------------------------------------------------------------
!
  end function check_ref

!======================================================================
!
end module
