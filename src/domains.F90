!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2012 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!******************************************************************************
!!
!! module: DOMAINS
!!
!!  This module handles the initialization of the problem domains.
!!
!!
!!******************************************************************************
!
module domains

! module variables are not implicit by default
!
  implicit none

! module variable to store the problem name
!
  character(len=32), save :: problem = "blast"

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_domains, setup_domain

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!!
!!***  PUBLIC SUBROUTINES  *****************************************************
!!
!===============================================================================
!
! subroutine INITIALIZE_DOMAINS:
! -----------------------------
!
!   Subroutine prepares module DOMAINS.
!
!
!===============================================================================
!
  subroutine initialize_domains()

! include external procedures and variables
!
    use parameters, only : get_parameter_string

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the problem name
!
    call get_parameter_string("problem", problem)

!-------------------------------------------------------------------------------
!
  end subroutine initialize_domains
!
!===============================================================================
!
! subroutine SETUP_DOMAIN:
! -----------------------
!
!   Subroutine sets up the domain for selected problem.  If there is no special
!   domain required, sets up the default domain.
!
!
!===============================================================================
!
  subroutine setup_domain()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! select the domain setup depending on the problem name
!
    select case(problem)
    case default
      call setup_domain_default()
    end select

!-------------------------------------------------------------------------------
!
  end subroutine setup_domain
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
! subroutine SETUP_DOMAIN_DEFAULT:
! -------------------------------
!
!   Subroutine sets the default domain of N₁xN₂xN₃ blocks in the right
!   configuration.
!
!
!===============================================================================
!
  subroutine setup_domain_default()

! include external procedures and variables
!
    use blocks     , only : pointer_meta, block_meta, block_data               &
                          , append_metablock, append_datablock                 &
                          , associate_blocks, metablock_set_leaf               &
                          , metablock_set_config, metablock_set_level          &
                          , metablock_set_coord, metablock_set_bounds
    use blocks     , only : nsides, nfaces
    use boundaries , only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry
    use coordinates, only : xmin, xmax, ymin, ymax, zmin, zmax
    use coordinates, only : ir, jr, kr, res

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: i, j, k, n, p, il, jl, kl
    real    :: xl, xmn, xmx, yl, ymn, ymx, zl, zmn, zmx

! local arrays
!
    integer, dimension(3) :: loc, del

! local pointers
!
    type(block_meta), pointer :: pmeta, pnext
    type(block_data), pointer :: pdata

! allocatable arrays
!
    integer, dimension(:,:,:), allocatable :: cfg

! local pointer array
!
    type(pointer_meta), dimension(:,:,:), allocatable :: block_array
!
!-------------------------------------------------------------------------------
!
! obtain the number of blocks
!
    n = ir * jr * kr

!! PREPARE BLOCK CONFIGURATION ARRAY
!!
! allocate the configuration array
!
    allocate(cfg(ir,jr,kr))

! set the block configurations
!
    cfg(1:ir,1:jr:2,1:kr:2) = 12

    if (jr .gt. 1) then
      cfg(1:ir,2:jr:2,1:kr:2) = 43
      cfg(  ir,1:jr  ,1:kr:2) = 13
    end if

    if (kr .gt. 1) then
      cfg(1:ir,1:jr:2,2:kr:2) = 65
      if (jr .gt. 1) then
        cfg(1:ir,2:jr:2,2:kr:2) = 78
        cfg(  ir,1:jr  ,2:kr:2) = 75
      end if
      if (ir .eq. 1 .or. mod(jr,2) .eq. 1) then
        cfg(  ir,  jr  ,1:kr  ) = 15
      else
        cfg(  1       ,  jr  ,1:kr  ) = 48
      end if
    end if

!! ALLOCATE AND GENERATE META BLOCK CHAIN AND SET BLOCK CONFIGURATIONS
!!
! allocate the block pointer array
!
    allocate(block_array(ir,jr,kr))

! generate the gray code for a given configuration and link the block in
! the proper order
!
    loc(:) = (/ 0, 0, 0 /)
    del(:) = (/ 1, 1, 1 /)

    p = 1
    do k = 1, kr
      if (del(3) .eq. 1) loc(3) = loc(3) + del(3)
      do j = 1, jr
        if (del(2) .eq. 1) loc(2) = loc(2) + del(2)
        do i = 1, ir
          if (del(1) .eq. 1) loc(1) = loc(1) + del(1)

! append a new metablock
!
          call append_metablock(block_array(loc(1),loc(2),loc(3))%ptr)

! set the configuration type
!
          call metablock_set_config(block_array(loc(1),loc(2),loc(3))%ptr      &
                                         , cfg(loc(1),loc(2),loc(3)))

! increase the block number
!
          p = p + 1

          if (del(1) .eq. -1) loc(1) = loc(1) + del(1)
        end do
        if (del(2) .eq. -1) loc(2) = loc(2) + del(2)
        del(1) = - del(1)
      end do
      if (del(3) .eq. -1) loc(3) = loc(3) + del(3)
      del(2) = - del(2)
    end do

! deallocate the configuration array
!
    deallocate(cfg)

!! FILL OUT THE REMAINING FIELDS AND ALLOCATE AND ASSOCIATE DATA BLOCKS
!!
! calculate block sizes
!
    xl = (xmax - xmin) / ir
    yl = (ymax - ymin) / jr
    zl = (zmax - zmin) / kr

! fill out block structure fields
!
    do k = 1, kr

! claculate the block position along Z
!
      kl  = (k - 1) * res(1,3)

! calculate the Z bounds
!
      zmn = zmin + (k - 1) * zl
      zmx = zmin +  k      * zl

      do j = 1, jr

! claculate the block position along Y
!
        jl  = (j - 1) * res(1,2)

! calculate the Y bounds
!
        ymn = ymin + (j - 1) * yl
        ymx = ymin +  j      * yl

        do i = 1, ir

! claculate the block position along Y
!
          il  = (i - 1) * res(1,1)

! calculate the Z bounds
!
          xmn = xmin + (i - 1) * xl
          xmx = xmin +  i      * xl

! assign a pointer
!
          pmeta => block_array(i,j,k)%ptr

! mark it as the leaf
!
          call metablock_set_leaf(pmeta)

! set the level
!
          call metablock_set_level(pmeta, 1)

! create a new data block
!
          call append_datablock(pdata)

! associate meta and data blocks
!
          call associate_blocks(pmeta, pdata)

! set block coordinates
!
          call metablock_set_coord(pmeta, il, jl, kl)

! set the bounds
!
          call metablock_set_bounds(pmeta, xmn, xmx, ymn, ymx, zmn, zmx)
        end do
      end do
    end do

!! ASSIGN THE BLOCK NEIGHBORS
!!
! assign boundaries along the X direction
!
    do k = 1, kr
      do j = 1, jr
        do i = 1, ir - 1

! assign a pointer
!
          pmeta => block_array(i  ,j,k)%ptr

! assign neighbor
!
          pnext => block_array(i+1,j,k)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(1,2,p)%ptr => pnext
            pnext%neigh(1,1,p)%ptr => pmeta
          end do

        end do
      end do
    end do

! if periodic boundary conditions set edge block neighbors
!
    if (xlbndry .eq. 'periodic' .and. xubndry .eq. 'periodic') then
      do k = 1, kr
        do j = 1, jr

! assign pointers
!
          pmeta => block_array(      1 ,j,k)%ptr
          pnext => block_array(ir,j,k)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(1,1,p)%ptr => pnext
            pnext%neigh(1,2,p)%ptr => pmeta
          end do
        end do
      end do
    end if

! assign boundaries along the Y direction
!
    do k = 1, kr
      do j = 1, jr - 1
        do i = 1, ir

! assign a pointer
!
          pmeta => block_array(i,j  ,k)%ptr

! assign neighbor
!
          pnext => block_array(i,j+1,k)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(2,2,p)%ptr => pnext
            pnext%neigh(2,1,p)%ptr => pmeta
          end do

        end do
      end do
    end do

! if periodic boundary conditions set edge block neighbors
!
    if (ylbndry .eq. 'periodic' .and. yubndry .eq. 'periodic') then
      do k = 1, kr
        do i = 1, ir

! assign pointers
!
          pmeta => block_array(i,      1 ,k)%ptr
          pnext => block_array(i,jr,k)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(2,1,p)%ptr => pnext
            pnext%neigh(2,2,p)%ptr => pmeta
          end do
        end do
      end do
    end if
#if NDIMS == 3

! assign boundaries along the Z direction
!
    do k = 1, kr - 1
      do j = 1, jr
        do i = 1, ir

! assign a pointer
!
          pmeta => block_array(i,j,k  )%ptr

! assign neighbor
!
          pnext => block_array(i,j,k+1)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(3,2,p)%ptr => pnext
            pnext%neigh(3,1,p)%ptr => pmeta
          end do

        end do
      end do
    end do

! if periodic boundary conditions set edge block neighbors
!
    if (zlbndry .eq. 'periodic' .and. zubndry .eq. 'periodic') then
      do j = 1, jr
        do i = 1, ir

! assign pointers
!
          pmeta => block_array(i,j,      1 )%ptr
          pnext => block_array(i,j,kr)%ptr

! assign their neighbor pointers
!
          do p = 1, nfaces
            pmeta%neigh(3,1,p)%ptr => pnext
            pnext%neigh(3,2,p)%ptr => pmeta
          end do
        end do
      end do
    end if
#endif /* NDIMS == 3 */

! deallocate the block pointer array
!
    deallocate(block_array)

!-------------------------------------------------------------------------------
!
  end subroutine setup_domain_default

!===============================================================================
!
end module domains
