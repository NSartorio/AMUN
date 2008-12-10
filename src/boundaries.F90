!!*****************************************************************************
!!
!! module: boundaries - routines for handling the boundary conditions
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
module boundaries

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! boundary: subroutine sweeps over all leaf blocks and performs the boundary
!           update
!
!===============================================================================
!
  subroutine boundary

    use blocks, only : block, plist, ndims
    use error , only : print_error

    implicit none

! local variables
!
    integer :: i, j, k, dl

! local pointers
!
    type(block), pointer :: pblock, pneigh
!
!-------------------------------------------------------------------------------
!
! iterate over all blocks and perform boundary update
!
    pblock => plist
    do while (associated(pblock))

! if the current block is a leaf...
!
      if (pblock%leaf .eq. 'T') then

! iterate over all neighbor blocks
!
        do i = 1, ndims
          do j = 1, 2
            do k = 1, 2

              pneigh => pblock%pneigh(i,j,k)%p

! check if neighbor is associated
!
              if (associated(pneigh)) then

! neighbor is associated, which means the periodic boundary conditions
! or interior of the domain
!
                if (pneigh%leaf .eq. 'T') then

! calculate the difference of current and neighbor levels
!
                  dl = pblock%level - pneigh%level

! depending on the level difference
!
                  select case(dl)
                  case(-1)  ! restriction

                  case(0)   ! the same level, copying
                    call bnd_copy(pblock%u, pneigh%u, i, j, k)
                  case(1)   ! prolongation

                  case default
                    call print_error("boundaries::boundary", "Level difference unsupported!")
                  end select

! perform copying, prolongation or restriction
!

                else
                  print *, pneigh%id, 'is not a leaf'
                endif
              else

! neighbor is not associated, it means that we have non periodic boundary here
!

              endif

            end do
          end do
        end do

      endif

! assign pointer to the next block
!
      pblock => pblock%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine boundary
!
!===============================================================================
!
! bnd_copy: subroutine copies the interior of neighbor to update the boundaries
!           of current block
!
!===============================================================================
!
  subroutine bnd_copy(u, b, id, is, ip)

    use blocks, only : nvars
    use config, only : igrids, jgrids, kgrids, nghost, ncells
    use error , only : print_warning

    implicit none

! arguments
!
    real, dimension(nvars,igrids,jgrids,kgrids), intent(inout) :: u
    real, dimension(nvars,igrids,jgrids,kgrids), intent(in)    :: b
    integer                                    , intent(in)    :: id, is, ip

! local variables
!
    integer :: ii
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * is

! perform update according to the flag
!
    select case(ii)
    case(110)
      u(:,1:nghost,:,:) = b(:,ncells:ncells+nghost,:,:)
    case(120)
      u(:,igrids-nghost:igrids,:,:) = b(:,nghost+1:2*nghost,:,:)
    case(210)
      u(:,:,1:nghost,:) = b(:,:,ncells:ncells+nghost,:)
    case(220)
      u(:,:,jgrids-nghost:jgrids,:) = b(:,:,nghost+1:2*nghost,:)
    case(310)
      u(:,:,:,1:nghost) = b(:,:,:,ncells:ncells+nghost)
    case(320)
      u(:,:,:,kgrids-nghost:kgrids) = b(:,:,:,nghost+1:2*nghost)
    case default
      call print_warning("boundaries::bnd_copy", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_copy

!===============================================================================
!
end module
