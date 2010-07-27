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
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! evolve: subroutine sweeps over all leaf blocks and performs one step time
!         evolution for each according to the selected integration scheme
!
!===============================================================================
!
  subroutine evolve

    use blocks      , only : block_data, list_data
    use boundaries  , only : boundary
    use mesh        , only : update_mesh
    use mesh        , only : dx_min
#ifdef MPI
    use mpitools    , only : mallreduceminr
#endif /* MPI */
    use scheme      , only : maxspeed
    use timer       , only : start_timer, stop_timer

    implicit none

! local variables
!
    type(block_data), pointer :: pblock
    real                      :: cmax, cm
!
!-------------------------------------------------------------------------------
!
! iterate over all data blocks and perform one step of time evolution
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a leaf
!
      if (pblock%meta%leaf) &
#ifdef EULER
        call evolve_euler(pblock)
#endif /* EULER */
#ifdef RK2
        call evolve_rk2(pblock)
#endif /* RK2 */

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! update boundaries
!
    call start_timer(4)
    call boundary
    call stop_timer(4)

! check refinement and refine
!
    call start_timer(5)
    call update_mesh
    call stop_timer(5)

! update boundaries
!
    call start_timer(4)
    call boundary
    call stop_timer(4)

! reset maximum speed
!
    cmax = 1.0e-8

! iterate over all blocks in order to find the maximum speed
!
    pblock => list_data
    do while (associated(pblock))

! check if this block is a leaf
!
      if (pblock%meta%leaf) &
        cm = maxspeed(pblock%u)

! compare global and local maximum speeds
!
      cmax = max(cmax, cm)

! assign pointer to the next block
!
      pblock => pblock%next

    end do

! get maximum time step
!
    dtn = dx_min / max(cmax, 1.e-8)

#ifdef MPI
! reduce new time step over all processes
!
    call mallreduceminr(dtn)
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve
#ifdef EULER
!
!===============================================================================
!
! evolve_euler: subroutine evolves the current block using Euler integration
!
!===============================================================================
!
  subroutine evolve_euler(pblock)

    use blocks       , only : block_data, nqt, nfl
    use config       , only : im, jm, km
    use mesh         , only : adxi, adyi, adzi
#ifdef SHAPE
    use problem      , only : update_shapes
#endif /* SHAPE */
    use scheme       , only : update

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi

! local arrays
!
    real, dimension(NDIMS,nfl,im,jm,km) :: f
    real, dimension(      nqt,im,jm,km) :: du
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! 1st step of integration
!
    call update(pblock%u, f, du, dxi, dyi, dzi)

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du)
#endif /* SHAPE */

! update solution
!
    pblock%u(:,:,:,:) = pblock%u(:,:,:,:) + dt * du(:,:,:,:)
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_euler
#endif /* EULER */
#ifdef RK2
!
!===============================================================================
!
! evolve_rk2: subroutine evolves the current block using RK2 integration
!
!===============================================================================
!
  subroutine evolve_rk2(pblock)

    use blocks       , only : block_data, nqt, nfl
    use config       , only : im, jm, km
    use mesh         , only : adxi, adyi, adzi
#ifdef SHAPE
    use problem      , only : update_shapes
#endif /* SHAPE */
    use scheme       , only : update

    implicit none

! input arguments
!
    type(block_data), intent(inout) :: pblock

! local variables
!
    real    :: dxi, dyi, dzi

! local arrays
!
    real, dimension(      nqt,im,jm,km) :: u1, du
    real, dimension(NDIMS,nfl,im,jm,km) :: f
!
!-------------------------------------------------------------------------------
!
! prepare dxi, dyi, and dzi
!
    dxi = adxi(pblock%meta%level)
    dyi = adyi(pblock%meta%level)
    dzi = adzi(pblock%meta%level)

! 1st step of integration
!
    call update(pblock%u, f, du, dxi, dyi, dzi)

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du)
#endif /* SHAPE */

! update solution
!
    u1(:,:,:,:) = pblock%u(:,:,:,:) + dt * du(:,:,:,:)

! 2nd step of integration
!
    call update(u1, f, du, dxi, dyi, dzi)

#ifdef SHAPE
! restrict update in a defined shape
!
    call update_shapes(pblock, du)
#endif /* SHAPE */

! update solution
!
    pblock%u(:,:,:,:) = 0.5 * (pblock%u(:,:,:,:) + u1(:,:,:,:) + dt * du(:,:,:,:))
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_rk2
#endif /* RK2 */
!
!===============================================================================
!
end module
