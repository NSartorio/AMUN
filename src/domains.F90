!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2016 Grzegorz Kowal <grzegorz@amuncode.org>
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
    use blocks         , only : pointer_meta, block_meta, block_data
    use blocks         , only : append_metablock, append_datablock, link_blocks
    use blocks         , only : metablock_set_leaf, metablock_set_level
    use blocks         , only : metablock_set_configuration
    use blocks         , only : metablock_set_coordinates, metablock_set_bounds
    use blocks         , only : nsides
    use boundaries     , only : bnd_type, bnd_periodic
    use coordinates    , only : xmin, ymin, zmin, xlen, ylen, zlen
    use coordinates    , only : ir, jr, kr

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer      :: i, j, k, n, p, ic, jc, kc
    real(kind=8) :: xl, xmn, xmx, yl, ymn, ymx, zl, zmn, zmx

! local arrays
!
    integer, dimension(3) :: loc, del

! local pointers
!
    type(block_meta), pointer :: pmeta, pnext

! allocatable arrays
!
    integer, dimension(:,:,:), allocatable :: cfg
    integer, dimension(:)    , allocatable :: im, jm, km
    integer, dimension(:)    , allocatable :: ip, jp, kp

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

    if (jr > 1) then
      cfg(1:ir,2:jr:2,1:kr:2) = 43
      cfg(  ir,1:jr  ,1:kr:2) = 13
    end if

    if (kr > 1) then
      cfg(1:ir,1:jr:2,2:kr:2) = 65
      if (jr > 1) then
        cfg(1:ir,2:jr:2,2:kr:2) = 78
        cfg(  ir,1:jr  ,2:kr:2) = 75
      end if
      if (ir == 1 .or. mod(jr,2) == 1) then
        cfg(  ir,  jr  ,1:kr  ) = 15
      else
        cfg(  1 ,  jr  ,1:kr  ) = 48
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
      if (del(3) == 1) loc(3) = loc(3) + del(3)
      do j = 1, jr
        if (del(2) == 1) loc(2) = loc(2) + del(2)
        do i = 1, ir
          if (del(1) == 1) loc(1) = loc(1) + del(1)

! append a new metablock
!
          call append_metablock(block_array(loc(1),loc(2),loc(3))%ptr)

! set the configuration type
!
          call metablock_set_configuration(                                    &
                                        block_array(loc(1),loc(2),loc(3))%ptr  &
                                                  , cfg(loc(1),loc(2),loc(3)))

! increase the block number
!
          p = p + 1

          if (del(1) == -1) loc(1) = loc(1) + del(1)
        end do
        if (del(2) == -1) loc(2) = loc(2) + del(2)
        del(1) = - del(1)
      end do
      if (del(3) == -1) loc(3) = loc(3) + del(3)
      del(2) = - del(2)
    end do

! deallocate the configuration array
!
    deallocate(cfg)

!! FILL OUT THE REMAINING FIELDS AND ALLOCATE AND ASSOCIATE DATA BLOCKS
!!
! calculate block sizes
!
    xl = xlen / ir
    yl = ylen / jr
    zl = zlen / kr

! fill out block structure fields
!
    do k = 1, kr

! claculate the block position along Z
!
      kc  = k - 1

! calculate the Z bounds
!
      zmn = zmin + kc * zl
      zmx = zmin + k  * zl

      do j = 1, jr

! claculate the block position along Y
!
        jc  = j - 1

! calculate the Y bounds
!
        ymn = ymin + jc * yl
        ymx = ymin + j  * yl

        do i = 1, ir

! claculate the block position along Y
!
          ic  = i - 1

! calculate the Z bounds
!
          xmn = xmin + ic * xl
          xmx = xmin + i  * xl

! assign a pointer
!
          pmeta => block_array(i,j,k)%ptr

! mark it as the leaf
!
          call metablock_set_leaf(pmeta)

! set the level
!
          call metablock_set_level(pmeta, 1)

! set block coordinates
!
          call metablock_set_coordinates(pmeta, ic, jc, kc)

! set the bounds
!
          call metablock_set_bounds(pmeta, xmn, xmx, ymn, ymx, zmn, zmx)
        end do
      end do
    end do

!! ASSIGN THE BLOCK NEIGHBORS
!!
! allocate indices
!
   allocate(im(ir), ip(ir))
   allocate(jm(jr), jp(jr))
#if NDIMS == 3
   allocate(km(kr), kp(kr))
#endif /* NDIMS == 3 */

! generate indices
!
   im(:) = cshift((/(i, i = 1, ir)/),-1)
   ip(:) = cshift((/(i, i = 1, ir)/), 1)
   jm(:) = cshift((/(j, j = 1, jr)/),-1)
   jp(:) = cshift((/(j, j = 1, jr)/), 1)
#if NDIMS == 3
   km(:) = cshift((/(k, k = 1, kr)/),-1)
   kp(:) = cshift((/(k, k = 1, kr)/), 1)
#endif /* NDIMS == 3 */

! check periodicity and reset the edge indices if box is not periodic
!
   if (bnd_type(1,1) /= bnd_periodic .or. bnd_type(1,2) /= bnd_periodic) then
     im( 1) = 0
     ip(ir) = 0
   end if
   if (bnd_type(2,1) /= bnd_periodic .or. bnd_type(2,2) /= bnd_periodic) then
     jm( 1) = 0
     jp(jr) = 0
   end if
#if NDIMS == 3
   if (bnd_type(3,1) /= bnd_periodic .or. bnd_type(3,2) /= bnd_periodic) then
     km( 1) = 0
     kp(kr) = 0
   end if
#endif /* NDIMS == 3 */

! iterate over all initial blocks
!
    do k = 1, kr
      do j = 1, jr
        do i = 1, ir

! assign pmeta with the current block
!
          pmeta => block_array(i,j,k)%ptr

#if NDIMS == 3
! assign face neighbor pointers
!
          if (im(i) > 0) then
            pmeta%faces(1,1,1,1)%ptr => block_array(im(i),j,k)%ptr
            pmeta%faces(1,2,1,1)%ptr => block_array(im(i),j,k)%ptr
            pmeta%faces(1,1,2,1)%ptr => block_array(im(i),j,k)%ptr
            pmeta%faces(1,2,2,1)%ptr => block_array(im(i),j,k)%ptr
          end if
          if (ip(i) > 0) then
            pmeta%faces(2,1,1,1)%ptr => block_array(ip(i),j,k)%ptr
            pmeta%faces(2,2,1,1)%ptr => block_array(ip(i),j,k)%ptr
            pmeta%faces(2,1,2,1)%ptr => block_array(ip(i),j,k)%ptr
            pmeta%faces(2,2,2,1)%ptr => block_array(ip(i),j,k)%ptr
          end if

          if (jm(j) > 0) then
            pmeta%faces(1,1,1,2)%ptr => block_array(i,jm(j),k)%ptr
            pmeta%faces(2,1,1,2)%ptr => block_array(i,jm(j),k)%ptr
            pmeta%faces(1,1,2,2)%ptr => block_array(i,jm(j),k)%ptr
            pmeta%faces(2,1,2,2)%ptr => block_array(i,jm(j),k)%ptr
          end if
          if (jp(j) > 0) then
            pmeta%faces(1,2,1,2)%ptr => block_array(i,jp(j),k)%ptr
            pmeta%faces(2,2,1,2)%ptr => block_array(i,jp(j),k)%ptr
            pmeta%faces(1,2,2,2)%ptr => block_array(i,jp(j),k)%ptr
            pmeta%faces(2,2,2,2)%ptr => block_array(i,jp(j),k)%ptr
          end if

          if (km(k) > 0) then
            pmeta%faces(1,1,1,3)%ptr => block_array(i,j,km(k))%ptr
            pmeta%faces(2,1,1,3)%ptr => block_array(i,j,km(k))%ptr
            pmeta%faces(1,2,1,3)%ptr => block_array(i,j,km(k))%ptr
            pmeta%faces(2,2,1,3)%ptr => block_array(i,j,km(k))%ptr
          end if
          if (kp(k) > 0) then
            pmeta%faces(1,1,2,3)%ptr => block_array(i,j,kp(k))%ptr
            pmeta%faces(2,1,2,3)%ptr => block_array(i,j,kp(k))%ptr
            pmeta%faces(1,2,2,3)%ptr => block_array(i,j,kp(k))%ptr
            pmeta%faces(2,2,2,3)%ptr => block_array(i,j,kp(k))%ptr
          end if
#endif /* NDIMS == 3 */

! assign edge neighbor pointers
!
#if NDIMS == 2
          if (im(i) > 0) then
            pmeta%edges(1,1,2)%ptr => block_array(im(i),j,k)%ptr
            pmeta%edges(1,2,2)%ptr => block_array(im(i),j,k)%ptr
          end if
          if (ip(i) > 0) then
            pmeta%edges(2,1,2)%ptr => block_array(ip(i),j,k)%ptr
            pmeta%edges(2,2,2)%ptr => block_array(ip(i),j,k)%ptr
          end if
          if (jm(j) > 0) then
            pmeta%edges(1,1,1)%ptr => block_array(i,jm(j),k)%ptr
            pmeta%edges(2,1,1)%ptr => block_array(i,jm(j),k)%ptr
          end if
          if (jp(j) > 0) then
            pmeta%edges(1,2,1)%ptr => block_array(i,jp(j),k)%ptr
            pmeta%edges(2,2,1)%ptr => block_array(i,jp(j),k)%ptr
          end if
#endif /* NDIMS == 2 */
#if NDIMS == 3
          if (jm(j) > 0 .and. km(k) > 0) then
            pmeta%edges(1,1,1,1)%ptr => block_array(i,jm(j),km(k))%ptr
            pmeta%edges(2,1,1,1)%ptr => block_array(i,jm(j),km(k))%ptr
          end if
          if (jp(j) > 0 .and. km(k) > 0) then
            pmeta%edges(1,2,1,1)%ptr => block_array(i,jp(j),km(k))%ptr
            pmeta%edges(2,2,1,1)%ptr => block_array(i,jp(j),km(k))%ptr
          end if
          if (jm(j) > 0 .and. kp(k) > 0) then
            pmeta%edges(1,1,2,1)%ptr => block_array(i,jm(j),kp(k))%ptr
            pmeta%edges(2,1,2,1)%ptr => block_array(i,jm(j),kp(k))%ptr
          end if
          if (jp(j) > 0 .and. kp(k) > 0) then
            pmeta%edges(1,2,2,1)%ptr => block_array(i,jp(j),kp(k))%ptr
            pmeta%edges(2,2,2,1)%ptr => block_array(i,jp(j),kp(k))%ptr
          end if

          if (im(i) > 0 .and. km(k) > 0) then
            pmeta%edges(1,1,1,2)%ptr => block_array(im(i),j,km(k))%ptr
            pmeta%edges(1,2,1,2)%ptr => block_array(im(i),j,km(k))%ptr
          end if
          if (ip(i) > 0 .and. km(k) > 0) then
            pmeta%edges(2,1,1,2)%ptr => block_array(ip(i),j,km(k))%ptr
            pmeta%edges(2,2,1,2)%ptr => block_array(ip(i),j,km(k))%ptr
          end if
          if (im(i) > 0 .and. kp(k) > 0) then
            pmeta%edges(1,1,2,2)%ptr => block_array(im(i),j,kp(k))%ptr
            pmeta%edges(1,2,2,2)%ptr => block_array(im(i),j,kp(k))%ptr
          end if
          if (ip(i) > 0 .and. kp(k) > 0) then
            pmeta%edges(2,1,2,2)%ptr => block_array(ip(i),j,kp(k))%ptr
            pmeta%edges(2,2,2,2)%ptr => block_array(ip(i),j,kp(k))%ptr
          end if

          if (im(i) > 0 .and. jm(j) > 0) then
            pmeta%edges(1,1,1,3)%ptr => block_array(im(i),jm(j),k)%ptr
            pmeta%edges(1,1,2,3)%ptr => block_array(im(i),jm(j),k)%ptr
          end if
          if (ip(i) > 0 .and. jm(j) > 0) then
            pmeta%edges(2,1,1,3)%ptr => block_array(ip(i),jm(j),k)%ptr
            pmeta%edges(2,1,2,3)%ptr => block_array(ip(i),jm(j),k)%ptr
          end if
          if (im(i) > 0 .and. jp(j) > 0) then
            pmeta%edges(1,2,1,3)%ptr => block_array(im(i),jp(j),k)%ptr
            pmeta%edges(1,2,2,3)%ptr => block_array(im(i),jp(j),k)%ptr
          end if
          if (ip(i) > 0 .and. jp(j) > 0) then
            pmeta%edges(2,2,1,3)%ptr => block_array(ip(i),jp(j),k)%ptr
            pmeta%edges(2,2,2,3)%ptr => block_array(ip(i),jp(j),k)%ptr
          end if
#endif /* NDIMS == 3 */

! assign corner neighbor pointers
!
#if NDIMS == 2
          if (im(i) > 0 .and. jm(j) > 0)                                       &
                      pmeta%corners(1,1)%ptr => block_array(im(i),jm(j),k)%ptr
          if (ip(i) > 0 .and. jm(j) > 0)                                       &
                      pmeta%corners(2,1)%ptr => block_array(ip(i),jm(j),k)%ptr
          if (im(i) > 0 .and. jp(j) > 0)                                       &
                      pmeta%corners(1,2)%ptr => block_array(im(i),jp(j),k)%ptr
          if (ip(i) > 0 .and. jp(j) > 0)                                       &
                      pmeta%corners(2,2)%ptr => block_array(ip(i),jp(j),k)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
          if (im(i) > 0 .and. jm(j) > 0 .and. km(k) > 0)                       &
                pmeta%corners(1,1,1)%ptr => block_array(im(i),jm(j),km(k))%ptr
          if (ip(i) > 0 .and. jm(j) > 0 .and. km(k) > 0)                       &
                pmeta%corners(2,1,1)%ptr => block_array(ip(i),jm(j),km(k))%ptr
          if (im(i) > 0 .and. jp(j) > 0 .and. km(k) > 0)                       &
                pmeta%corners(1,2,1)%ptr => block_array(im(i),jp(j),km(k))%ptr
          if (ip(i) > 0 .and. jp(j) > 0 .and. km(k) > 0)                       &
                pmeta%corners(2,2,1)%ptr => block_array(ip(i),jp(j),km(k))%ptr
          if (im(i) > 0 .and. jm(j) > 0 .and. kp(k) > 0)                       &
                pmeta%corners(1,1,2)%ptr => block_array(im(i),jm(j),kp(k))%ptr
          if (ip(i) > 0 .and. jm(j) > 0 .and. kp(k) > 0)                       &
                pmeta%corners(2,1,2)%ptr => block_array(ip(i),jm(j),kp(k))%ptr
          if (im(i) > 0 .and. jp(j) > 0 .and. kp(k) > 0)                       &
                pmeta%corners(1,2,2)%ptr => block_array(im(i),jp(j),kp(k))%ptr
          if (ip(i) > 0 .and. jp(j) > 0 .and. kp(k) > 0)                       &
                pmeta%corners(2,2,2)%ptr => block_array(ip(i),jp(j),kp(k))%ptr
#endif /* NDIMS == 3 */
        end do ! over i
      end do ! over j
    end do ! over k

! deallocate indices
!
    deallocate(im, ip)
    deallocate(jm, jp)
#if NDIMS == 3
    deallocate(km, kp)
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
