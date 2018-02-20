!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: BOUNDARIES
!!
!!  This module handles the boundary synchronization.
!!
!!******************************************************************************
!
module boundaries

! import external subroutines
!
#ifdef MPI
  use blocks, only : pointer_info
#endif /* MPI */
#ifdef PROFILE
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer                , save :: imi, imv, imf, ims, imu
  integer                , save :: ifc, ifr, ifp, iec, ier, iep, icc, icr, icp
#endif /* PROFILE */

! parameters corresponding to the boundary type
!
  integer, parameter            :: bnd_periodic   = 0
  integer, parameter            :: bnd_open       = 1
  integer, parameter            :: bnd_outflow    = 2
  integer, parameter            :: bnd_reflective = 3
  integer, parameter            :: bnd_gravity    = 4
  integer, parameter            :: bnd_user       = 5

! variable to store boundary type flags
!
  integer, dimension(3,2), save :: bnd_type       = bnd_periodic

#ifdef MPI
! arrays to store information about blocks which need to be exchange between
! processes
!
    type(pointer_info), dimension(:,:), allocatable, save :: barray
    integer           , dimension(:,:), allocatable, save :: bcount
#endif /* MPI */

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_boundaries, finalize_boundaries
  public :: boundary_variables, boundary_fluxes
  public :: bnd_type, bnd_periodic

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
! subroutine INITIALIZE_BOUNDARIES:
! --------------------------------
!
!   Subroutine initializes the module BOUNDARIES by setting its parameters.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_boundaries(verbose, iret)

! import external procedures and variables
!
    use coordinates    , only : periodic
#ifdef MPI
    use mpitools       , only : npmax
#endif /* MPI */
    use parameters     , only : get_parameter_string

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! module parameters for the boundary update order and boundary type
!
    character(len = 32)    :: xlbndry = "periodic"
    character(len = 32)    :: xubndry = "periodic"
    character(len = 32)    :: ylbndry = "periodic"
    character(len = 32)    :: yubndry = "periodic"
    character(len = 32)    :: zlbndry = "periodic"
    character(len = 32)    :: zubndry = "periodic"

! local variables
!
    integer                :: n
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('boundaries:: initialization' , imi)
    call set_timer('boundaries:: variables'      , imv)
    call set_timer('boundaries:: fluxes'         , imf)
    call set_timer('boundaries:: specific'       , ims)
    call set_timer('boundaries:: face copy'      , ifc)
    call set_timer('boundaries:: face restrict'  , ifr)
    call set_timer('boundaries:: face prolong'   , ifp)
    call set_timer('boundaries:: edge copy'      , iec)
    call set_timer('boundaries:: edge restrict'  , ier)
    call set_timer('boundaries:: edge prolong'   , iep)
    call set_timer('boundaries:: corner copy'    , icc)
    call set_timer('boundaries:: corner restrict', icr)
    call set_timer('boundaries:: corner prolong' , icp)
    call set_timer('boundaries:: update ghosts'  , imu)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get runtime values for the boundary types
!
    call get_parameter_string ("xlbndry" , xlbndry)
    call get_parameter_string ("xubndry" , xubndry)
    call get_parameter_string ("ylbndry" , ylbndry)
    call get_parameter_string ("yubndry" , yubndry)
    call get_parameter_string ("zlbndry" , zlbndry)
    call get_parameter_string ("zubndry" , zubndry)

! fill the boundary type flags
!
    select case(xlbndry)
    case("open")
      bnd_type(1,1) = bnd_open
    case("outflow", "out")
      bnd_type(1,1) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(1,1) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(1,1) = bnd_gravity
    case("user", "custom")
      bnd_type(1,1) = bnd_user
    case default
      bnd_type(1,1) = bnd_periodic
    end select

    select case(xubndry)
    case("open")
      bnd_type(1,2) = bnd_open
    case("outflow", "out")
      bnd_type(1,2) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(1,2) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(1,2) = bnd_gravity
    case("user", "custom")
      bnd_type(1,2) = bnd_user
    case default
      bnd_type(1,2) = bnd_periodic
    end select

    select case(ylbndry)
    case("open")
      bnd_type(2,1) = bnd_open
    case("outflow", "out")
      bnd_type(2,1) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(2,1) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(2,1) = bnd_gravity
    case("user", "custom")
      bnd_type(2,1) = bnd_user
    case default
      bnd_type(2,1) = bnd_periodic
    end select

    select case(yubndry)
    case("open")
      bnd_type(2,2) = bnd_open
    case("outflow", "out")
      bnd_type(2,2) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(2,2) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(2,2) = bnd_gravity
    case("user", "custom")
      bnd_type(2,2) = bnd_user
    case default
      bnd_type(2,2) = bnd_periodic
    end select

    select case(zlbndry)
    case("open")
      bnd_type(3,1) = bnd_open
    case("outflow", "out")
      bnd_type(3,1) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(3,1) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(3,1) = bnd_gravity
    case("user", "custom")
      bnd_type(3,1) = bnd_user
    case default
      bnd_type(3,1) = bnd_periodic
    end select

    select case(zubndry)
    case("open")
      bnd_type(3,2) = bnd_open
    case("outflow", "out")
      bnd_type(3,2) = bnd_outflow
    case("reflective", "reflecting", "reflect")
      bnd_type(3,2) = bnd_reflective
    case("hydrostatic", "gravity")
      bnd_type(3,2) = bnd_gravity
    case("user", "custom")
      bnd_type(3,2) = bnd_user
    case default
      bnd_type(3,2) = bnd_periodic
    end select

! set domain periodicity
!
    do n = 1, NDIMS
      periodic(n) = (bnd_type(n,1) == bnd_periodic) .and.                      &
                    (bnd_type(n,2) == bnd_periodic)
    end do

#ifdef MPI
! allocate the exchange arrays
!
    allocate(barray(0:npmax,0:npmax), bcount(0:npmax,0:npmax))

! prepare the exchange arrays
!
    call prepare_exchange_array()
#endif /* MPI */

! print information about the boundary conditions
!
    if (verbose) then

      write (*,"(4x,a10,13x,'=',2(1x,a))") "x-boundary"                        &
                                                , trim(xlbndry), trim(xubndry)
      write (*,"(4x,a10,13x,'=',2(1x,a))") "y-boundary"                        &
                                                , trim(ylbndry), trim(yubndry)
#if NDIMS == 3
      write (*,"(4x,a10,13x,'=',2(1x,a))") "z-boundary"                        &
                                                , trim(zlbndry), trim(zubndry)
#endif /* NDIMS == 3 */

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_boundaries
!
!===============================================================================
!
! subroutine FINALIZE_BOUNDARIES:
! ------------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_boundaries(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

#ifdef MPI
! release the exchange arrays
!
    call release_exchange_array()

! deallocate the exchange arrays
!
    deallocate(barray, bcount)
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_boundaries
!
!===============================================================================
!
! subroutine BOUNDARY_VARIABLES:
! -----------------------------
!
!   Subroutine updates the ghost zones of the data blocks from their neighbors
!   or applies the specific boundary conditions.
!
!   Arguments:
!
!     t, dt   - time and time increment;
!
!===============================================================================
!
  subroutine boundary_variables(t, dt)

! import external procedures and variables
!
    use blocks         , only : ndims
    use coordinates    , only : minlev, maxlev

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(in) :: t, dt

! local variables
!
    integer :: idir
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for variable boundary update
!
    call start_timer(imv)
#endif /* PROFILE */

#if NDIMS == 3
! update face boundaries between blocks at the same levels
!
    do idir = 1, ndims
      call boundaries_face_copy(idir)
    end do ! idir
#endif /* NDIMS == 3 */

! update edge boundaries between blocks at the same levels
!
    do idir = 1, ndims
      call boundaries_edge_copy(idir)
    end do ! idir

! update corner boundaries between blocks at the same levels
!
    call boundaries_corner_copy()

! do prolongation and restriction only if blocks are at different levels
!
    if (minlev /= maxlev) then

#if NDIMS == 3
! restrict face boundaries from higher level blocks
!
      do idir = 1, ndims
        call boundaries_face_restrict(idir)
      end do ! idir
#endif /* NDIMS == 3 */

! restricts edge boundaries from block at higher level
!
      do idir = 1, ndims
        call boundaries_edge_restrict(idir)
      end do ! idir

! restricts corner boundaries from blocks at higher levels
!
      call boundaries_corner_restrict()

! update specific boundaries
!
      call boundaries_specific(t, dt)

#if NDIMS == 3
! prolong face boundaries from lower level blocks
!
       do idir = 1, ndims
        call boundaries_face_prolong(idir)
       end do ! idir
#endif /* NDIMS == 3 */

! prolongs edge boundaries from block at lower level
!
      do idir = 1, ndims
        call boundaries_edge_prolong(idir)
      end do ! idir

! prolong corner boundaries from blocks at lower levels
!
      call boundaries_corner_prolong()

    end if ! minlev /= maxlev

! update specific boundaries
!
    call boundaries_specific(t, dt)

! convert updated primitive variables to conservative ones in all ghost cells
!
    call update_ghost_cells()

#ifdef PROFILE
! stop accounting time for variable boundary update
!
    call stop_timer(imv)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_variables
!
!===============================================================================
!
! subroutine BOUNDARY_FLUXES:
! --------------------------
!
!   Subroutine updates the numerical fluxes from neighbors which lay on
!   higher level. At higher levels the numerical fluxes are calculated with
!   smaller error, since the resolution is higher, therefore we take those
!   fluxes and restrict them down to the level of the updated block.
!
!
!===============================================================================
!
  subroutine boundary_fluxes()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
#ifdef MPI
    use blocks         , only : block_info, pointer_info
#endif /* MPI */
    use blocks         , only : ndims, nsides
    use coordinates    , only : minlev, maxlev
    use coordinates    , only : in, jn, kn
    use coordinates    , only : ib, ie, ibl
    use coordinates    , only : jb, je, jbl
    use coordinates    , only : kb, ke, kbl
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: n, m
    integer :: i, is, it, il, iu, ih
    integer :: j, js, jt, jl, ju, jh
    integer :: k, ks, kt, kl, ku, kh
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable  :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! quit if all blocks are at the same level
!
    if (minlev == maxlev) return

#ifdef PROFILE
! start accounting time for the flux boundary update
!
    call start_timer(imf)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
#if NDIMS == 2
    kh = kn
    kl = 1
    ku = kn
#endif /* NDIMS == 2 */
#if NDIMS == 3
    kh = kn / 2
#endif /* NDIMS == 3 */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update the fluxes between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! iterate over all dimensions
!
      do n = 1, ndims
#if NDIMS == 2
        m = 3 - n
#endif /* NDIMS == 2 */

! iterate over all corners
!
#if NDIMS == 3
        do k = 1, nsides
#endif /* NDIMS == 3 */
          do j = 1, nsides
            do i = 1, nsides

! associate pneigh with the neighbor
!
#if NDIMS == 2
              pneigh => pmeta%edges(i,j,m)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pneigh => pmeta%faces(i,j,k,n)%ptr
#endif /* NDIMS == 3 */

! process only if the neighbor is associated
!
              if (associated(pneigh)) then

! check if the neighbor lays at higher level
!
                if (pneigh%level > pmeta%level) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! update the flux depending on the direction
!
                      select case(n)
                      case(1)

! prepare the boundary layer indices for X-direction flux
!
                        if (i == 1) then
                          is = ie
                          it = ibl
                        else
                          is = ibl
                          it = ie
                        end if
                        if (j == 1) then
                          jl = jb
                          ju = jb + jh - 1
                        else
                          jl = je - jh + 1
                          ju = je
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = kb
                          ku = kb + kh - 1
                        else
                          kl = ke - kh + 1
                          ku = ke
                        end if
#endif /* NDIMS == 3 */

! update the flux at the X-face of the block
!
                        call block_update_flux(i, j, k, n                      &
                                       , pneigh%data%f(n,1:nv,is,jb:je,kb:ke)  &
                                       ,  pmeta%data%f(n,1:nv,it,jl:ju,kl:ku))

                      case(2)

! prepare the boundary layer indices for Y-direction flux
!
                        if (i == 1) then
                          il = ib
                          iu = ib + ih - 1
                        else
                          il = ie - ih + 1
                          iu = ie
                        end if
                        if (j == 1) then
                          js = je
                          jt = jbl
                        else
                         js = jbl
                          jt = je
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = kb
                          ku = kb + kh - 1
                        else
                          kl = ke - kh + 1
                          ku = ke
                        end if
#endif /* NDIMS == 3 */

! update the flux at the Y-face of the block
!
                        call block_update_flux(i, j, k, n                      &
                                       , pneigh%data%f(n,1:nv,ib:ie,js,kb:ke)  &
                                       ,  pmeta%data%f(n,1:nv,il:iu,jt,kl:ku))

#if NDIMS == 3
                      case(3)

! prepare the boundary layer indices for Z-direction flux
!
                        if (i == 1) then
                          il = ib
                          iu = ib + ih - 1
                        else
                          il = ie - ih + 1
                          iu = ie
                        end if
                        if (j == 1) then
                          jl = jb
                          ju = jb + jh - 1
                        else
                          jl = je - jh + 1
                          ju = je
                        end if
                        if (k == 1) then
                          ks = ke
                          kt = kbl
                        else
                          ks = kbl
                          kt = ke
                        end if

! update the flux at the Z-face of the block
!
                        call block_update_flux(i, j, k, n                      &
                                       , pneigh%data%f(n,1:nv,ib:ie,jb:je,ks)  &
                                       ,  pmeta%data%f(n,1:nv,il:iu,jl:ju,kt))
#endif /* NDIMS == 3 */

                      end select

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! pneigh%proc /= pmeta%proc

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, n, i, j, k)

                  end if ! pneigh%proc /= pmeta%proc
#endif /* MPI */
                end if ! pmeta level < pneigh level

              end if ! pneigh associated

            end do ! i = 1, nsides
          end do ! j = 1, nsides
#if NDIMS == 3
        end do ! k = 1, nsides
#endif /* NDIMS == 3 */
      end do ! n = 1, ndims

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
! update flux boundaries between neighbors laying on different processes
!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 1
          rtag = 16 * (sproc * nprocs + rproc) + 1

! allocate buffers for variable exchange
!
          allocate(sbuf(scount,nv,ih,kh))
          allocate(rbuf(rcount,nv,ih,kh))

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan all blocks on the list
!
          do while(associated(pinfo))

! increase the block count
!
            l = l + 1

! associate pneigh pointer
!
            pneigh => pinfo%neigh

! get neighbor direction and corner coordinates
!
            n = pinfo%direction
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! update directional flux from the neighbor
!
            select case(n)
            case(1)

! prepare the boundary layer index depending on the side
!
              if (i == 1) then
                is = ie
              else
                is = ibl
              end if

! update the flux edge from the neighbor at higher level
!
              call block_update_flux(i, j, k, n                                &
                                 , pneigh%data%f(n,1:nv,is,jb:je,kb:ke)        &
                                 ,          sbuf(l,1:nv,1:jh,1:kh))

            case(2)

! prepare the boundary layer index depending on the side
!
              if (j == 1) then
                js = je
              else
                js = jbl
              end if

! update the flux edge from the neighbor at higher level
!
              call block_update_flux(i, j, k, n                                &
                                 , pneigh%data%f(n,1:nv,ib:ie,js,kb:ke)        &
                                 ,          sbuf(l,1:nv,1:ih,1:kh))

#if NDIMS == 3
            case(3)

! prepare the boundary layer index depending on the side
!
              if (k == 1) then
                ks = ke
              else
                ks = kbl
              end if

! update the flux edge from the neighbor at higher level
!
              call block_update_flux(i, j, k, n                                &
                                 , pneigh%data%f(n,1:nv,ib:ie,jb:je,ks)        &
                                 ,          sbuf(l,1:nv,1:ih,1:jh))
#endif /* NDIMS == 3 */

            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr blocks

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! scan all blocks on the list
!
          do while(associated(pinfo))

! increase the block count
!
            l = l + 1

! associate pmeta pointer
!
            pmeta  => pinfo%block

! get neighbor direction and corner indices
!
            n = pinfo%direction
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! update directional flux from the neighbor
!
            select case(n)
            case(1)

! prepare the boundary layer indices depending on the corner position
!
              if (i == 1) then
                it = ibl
              else
                it = ie
              end if
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
#if NDIMS == 3
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
#endif /* NDIMS == 3 */

! update the flux edge from the neighbor at higher level
!
              pmeta%data%f(n,1:nv,it,jl:ju,kl:ku) = rbuf(l,1:nv,1:jh,1:kh)

            case(2)

! prepare the boundary layer indices depending on the corner position
!
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
              if (j == 1) then
                jt = jbl
              else
                jt = je
              end if
#if NDIMS == 3
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
#endif /* NDIMS == 3 */

! update the flux edge from the neighbor at higher level
!
              pmeta%data%f(n,1:nv,il:iu,jt,kl:ku) = rbuf(l,1:nv,1:ih,1:kh)

#if NDIMS == 3
            case(3)

! prepare the boundary layer indices depending on the corner position
!
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
              if (k == 1) then
                kt = kbl
              else
                kt = ke
              end if

! update the flux edge from the neighbor at higher level
!
              pmeta%data%f(n,1:nv,il:iu,jl:ju,kt) = rbuf(l,1:nv,1:ih,1:jh)
#endif /* NDIMS == 3 */

            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr blocks

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the flux boundary update
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_fluxes
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
!===============================================================================
!
!  DOMAIN SPECIFIC BOUNDARY SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BOUNDARIES_SPECIFIC:
! ------------------------------
!
!   Subroutine scans over all leaf blocks in order to find blocks without
!   neighbors and update the corresponding boundaries for the selected
!   boundary type.
!
!   Arguments:
!
!     t, dt   - time and time increment;
!
!===============================================================================
!
  subroutine boundaries_specific(t, dt)

! import external procedures and variables
!
    use blocks         , only : block_meta, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : ndims, nsides
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az
    use coordinates    , only : periodic
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), intent(in) :: t, dt

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf

! local variables
!
    integer                   :: i, j, k, n, m

! local arrays
!
    real(kind=8), dimension(im) :: x
    real(kind=8), dimension(jm) :: y
    real(kind=8), dimension(km) :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the specific boundary update
!
    call start_timer(ims)
#endif /* PROFILE */

! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! process only if this block is marked for update
!
      if (pmeta%update) then

#ifdef MPI
! check if the block belongs to the local process
!
        if (pmeta%process == nproc) then
#endif /* MPI */

! prepare block coordinates
!
        x(1:im) = pmeta%xmin + ax(pmeta%level,1:im)
        y(1:jm) = pmeta%ymin + ay(pmeta%level,1:jm)
#if NDIMS == 3
        z(1:km) = pmeta%zmin + az(pmeta%level,1:km)
#else /* NDIMS == 3 */
        z(1:km) = 0.0d+00
#endif /* NDIMS == 3 */

#if NDIMS == 2
! iterate over all directions
!
          do n = 1, ndims

! process boundaries only if they are not periodic along the given direction
!
            if (.not. periodic(n)) then

! calculate the edge direction (in 2D we don't have face neighbors, so we have
! to use edge neighbors)
!
              m = 3 - n

! iterate over all corners
!
              do j = 1, nsides
                do i = 1, nsides

! if the face neighbor is not associated, apply specific boundaries
!
                  if (.not. associated(pmeta%edges(i,j,m)%ptr))                &
                            call block_boundary_specific(i, j, k, n            &
                                          , t, dt, x(:), y(:), z(:)            &
                                          , pmeta%data%q(1:nv,1:im,1:jm,1:km))

                end do ! i = 1, sides
              end do ! j = 1, sides

            end if ! not periodic

          end do ! n = 1, ndims
#endif /* NDIMS == 2 */
#if NDIMS == 3
! iterate over all directions
!
          do n = 1, ndims

! process boundaries only if they are not periodic along the given direction
!
            if (.not. periodic(n)) then

! iterate over all corners
!
              do k = 1, nsides
                do j = 1, nsides
                  do i = 1, nsides

! if the face neighbor is not associated, apply specific boundaries
!
                    if (.not. associated(pmeta%faces(i,j,k,n)%ptr))            &
                            call block_boundary_specific(i, j, k, n            &
                                          , t, dt, x(:), y(:), z(:)            &
                                          , pmeta%data%q(1:nv,1:im,1:jm,1:km))

                  end do ! i = 1, sides
                end do ! j = 1, sides
              end do ! k = 1, sides

            end if ! not periodic

          end do ! n = 1, ndims
#endif /* NDIMS == 3 */

#ifdef MPI
        end if ! block belongs to the local process
#endif /* MPI */

      end if ! if pmeta marked for update

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef PROFILE
! stop accounting time for the specific boundary update
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_specific
#if NDIMS == 3
!
!===============================================================================
!
!  DOMAIN FACE BOUNDARY UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BOUNDARIES_FACE_COPY:
! -------------------------------
!
!   Subroutine updates the face boundaries between blocks on the same level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_face_copy(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : faces_gc, faces_dc
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: is, js, ks
    integer :: it, jt, kt
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the face boundary update by copying
!
    call start_timer(ifc)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
    kh = kn / 2

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
      do k = 1, nsides
        do j = 1, nsides
          do i = 1, nsides

! associate pneigh with the current neighbor
!
            pneigh => pmeta%faces(i,j,k,idir)%ptr

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at the same level
!
              if (pneigh%level == pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare region indices of the block and its neighbor for the face boundary
! update
!
                      il = faces_gc(i,j,k,idir)%l(1)
                      jl = faces_gc(i,j,k,idir)%l(2)
                      kl = faces_gc(i,j,k,idir)%l(3)
                      iu = faces_gc(i,j,k,idir)%u(1)
                      ju = faces_gc(i,j,k,idir)%u(2)
                      ku = faces_gc(i,j,k,idir)%u(3)

                      is = faces_dc(i,j,k,idir)%l(1)
                      js = faces_dc(i,j,k,idir)%l(2)
                      ks = faces_dc(i,j,k,idir)%l(3)
                      it = faces_dc(i,j,k,idir)%u(1)
                      jt = faces_dc(i,j,k,idir)%u(2)
                      kt = faces_dc(i,j,k,idir)%u(3)

! copy the corresponding face region from the neighbor to the current data
! block
!
                      pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                   &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
      end do ! k = 1, nsides

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 2
          rtag = 16 * (sproc * nprocs + rproc) + 2

! allocate data buffer for variables to exchange
!
          select case(idir)
          case(1)
            allocate(sbuf(scount,nv,ng,jh,kh))
            allocate(rbuf(rcount,nv,ng,jh,kh))
          case(2)
            allocate(sbuf(scount,nv,ih,ng,kh))
            allocate(rbuf(rcount,nv,ih,ng,kh))
          case(3)
            allocate(sbuf(scount,nv,ih,jh,ng))
            allocate(rbuf(rcount,nv,ih,jh,ng))
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! prepare region indices for the face boundary update
!
            is = faces_dc(i,j,k,idir)%l(1)
            js = faces_dc(i,j,k,idir)%l(2)
            ks = faces_dc(i,j,k,idir)%l(3)
            it = faces_dc(i,j,k,idir)%u(1)
            jt = faces_dc(i,j,k,idir)%u(2)
            kt = faces_dc(i,j,k,idir)%u(3)

! copy the corresponding face region from the neighbor and insert it to
! the buffer
!
            select case(idir)
            case(1)
              sbuf(l,1:nv,1:ng,1:jh,1:kh) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
            case(2)
              sbuf(l,1:nv,1:ih,1:ng,1:kh) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
            case(3)
              sbuf(l,1:nv,1:ih,1:jh,1:ng) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! prepare region indices for the face boundary update
!
            il = faces_gc(i,j,k,idir)%l(1)
            jl = faces_gc(i,j,k,idir)%l(2)
            kl = faces_gc(i,j,k,idir)%l(3)
            iu = faces_gc(i,j,k,idir)%u(1)
            ju = faces_gc(i,j,k,idir)%u(2)
            ku = faces_gc(i,j,k,idir)%u(3)

! update the corresponding face region of the current block
!
            select case(idir)
            case(1)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:jh,1:kh)
            case(2)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ih,1:ng,1:kh)
            case(3)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ih,1:jh,1:ng)
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the face boundary update by copying
!
    call stop_timer(ifc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_face_copy
!
!===============================================================================
!
! subroutine BOUNDARIES_FACE_RESTRICT:
! -----------------------------------
!
!   Subroutine updates the face boundaries from blocks on higher level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_face_restrict(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : faces_gr
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the face boundary update by restriction
!
    call start_timer(ifr)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
    kh = kn / 2

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
      do k = 1, nsides
        do j = 1, nsides
          do i = 1, nsides

! associate pneigh with the current neighbor
!
            pneigh => pmeta%faces(i,j,k,idir)%ptr

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at higher level
!
              if (pneigh%level > pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare region indices of the block and its neighbor for the face boundary
! update
!
                      il = faces_gr(i,j,k,idir)%l(1)
                      jl = faces_gr(i,j,k,idir)%l(2)
                      kl = faces_gr(i,j,k,idir)%l(3)
                      iu = faces_gr(i,j,k,idir)%u(1)
                      ju = faces_gr(i,j,k,idir)%u(2)
                      ku = faces_gr(i,j,k,idir)%u(3)

! extract the corresponding face region from the neighbor and insert it in
! the current data block
!
                      call block_face_restrict(idir, i, j, k                   &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
      end do ! k = 1, nsides

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 3
          rtag = 16 * (sproc * nprocs + rproc) + 3

! allocate data buffer for variables to exchange
!
          select case(idir)
          case(1)
            allocate(sbuf(scount,nv,ng,jh,kh))
            allocate(rbuf(rcount,nv,ng,jh,kh))
          case(2)
            allocate(sbuf(scount,nv,ih,ng,kh))
            allocate(rbuf(rcount,nv,ih,ng,kh))
          case(3)
            allocate(sbuf(scount,nv,ih,jh,ng))
            allocate(rbuf(rcount,nv,ih,jh,ng))
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! restrict the corresponding face region from the neighbor and insert it
! to the buffer
!
            select case(idir)
            case(1)
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:kh))
            case(2)
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:kh))
            case(3)
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:jh,1:ng))
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! prepare region indices for the face boundary update
!
            il = faces_gr(i,j,k,idir)%l(1)
            jl = faces_gr(i,j,k,idir)%l(2)
            kl = faces_gr(i,j,k,idir)%l(3)
            iu = faces_gr(i,j,k,idir)%u(1)
            ju = faces_gr(i,j,k,idir)%u(2)
            ku = faces_gr(i,j,k,idir)%u(3)

! update the corresponding face region of the current block
!
            select case(idir)
            case(1)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:kh)
            case(2)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:kh)
            case(3)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:jh,1:ng)
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the face boundary update by restriction
!
    call stop_timer(ifr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_face_restrict
!
!===============================================================================
!
! subroutine BOUNDARIES_FACE_PROLONG:
! ----------------------------------
!
!   Subroutine updates the face boundaries from blocks on lower level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_face_prolong(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : faces_gp
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ic, jc, kc
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the face boundary update by prolongation
!
    call start_timer(ifp)
#endif /* PROFILE */

! calculate the sizes
!
    ih = in + ng
    jh = jn + ng
    kh = kn + ng

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
      do k = 1, nsides
        kc = k
        do j = 1, nsides
          jc = j
          do i = 1, nsides
            ic = i

! associate pneigh with the current neighbor
!
            pneigh => pmeta%faces(i,j,k,idir)%ptr

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor lays at lower level
!
              if (pneigh%level < pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare indices of the region in which the boundaries should be updated
!
                      select case(idir)
                      case(1)
                        jc = pmeta%pos(2) + 1
                        kc = pmeta%pos(3) + 1
                        il = faces_gp(i ,jc,kc,idir)%l(1)
                        jl = faces_gp(i ,jc,kc,idir)%l(2)
                        kl = faces_gp(i ,jc,kc,idir)%l(3)
                        iu = faces_gp(i ,jc,kc,idir)%u(1)
                        ju = faces_gp(i ,jc,kc,idir)%u(2)
                        ku = faces_gp(i ,jc,kc,idir)%u(3)

                      case(2)
                        ic = pmeta%pos(1) + 1
                        kc = pmeta%pos(3) + 1
                        il = faces_gp(ic,j ,kc,idir)%l(1)
                        jl = faces_gp(ic,j ,kc,idir)%l(2)
                        kl = faces_gp(ic,j ,kc,idir)%l(3)
                        iu = faces_gp(ic,j ,kc,idir)%u(1)
                        ju = faces_gp(ic,j ,kc,idir)%u(2)
                        ku = faces_gp(ic,j ,kc,idir)%u(3)

                      case(3)
                        ic = pmeta%pos(1) + 1
                        jc = pmeta%pos(2) + 1
                        il = faces_gp(ic,jc,k ,idir)%l(1)
                        jl = faces_gp(ic,jc,k ,idir)%l(2)
                        kl = faces_gp(ic,jc,k ,idir)%l(3)
                        iu = faces_gp(ic,jc,k ,idir)%u(1)
                        ju = faces_gp(ic,jc,k ,idir)%u(2)
                        ku = faces_gp(ic,jc,k ,idir)%u(3)
                      end select

! take the neighbor volume, extract the corresponding face region and insert
! it in the current data block
!
                      call block_face_prolong(idir, ic, jc, kc                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */
                end if ! pmeta and pneigh marked for update

              end if ! neighbor at lower level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
      end do ! k = 1, nsides

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 4
          rtag = 16 * (sproc * nprocs + rproc) + 4

! allocate data buffer for variables to exchange
!
          select case(idir)
          case(1)
            allocate(sbuf(scount,nv,ng,jh,kh))
            allocate(rbuf(rcount,nv,ng,jh,kh))
          case(2)
            allocate(sbuf(scount,nv,ih,ng,kh))
            allocate(rbuf(rcount,nv,ih,ng,kh))
          case(3)
            allocate(sbuf(scount,nv,ih,jh,ng))
            allocate(rbuf(rcount,nv,ih,jh,ng))
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pmeta and pneigh to the right blocks
!
            pmeta  => pinfo%block
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! extract the corresponding face region from the neighbor and insert it
! to the buffer
!
            select case(idir)
            case(1)
              j = pmeta%pos(2) + 1
              k = pmeta%pos(3) + 1
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:kh))
            case(2)
              i = pmeta%pos(1) + 1
              k = pmeta%pos(3) + 1
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:kh))
            case(3)
              i = pmeta%pos(1) + 1
              j = pmeta%pos(2) + 1
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:jh,1:ng))
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
            k = pinfo%corner(3)

! update the corresponding face region of the current block
!
            select case(idir)
            case(1)

              jc = pmeta%pos(2) + 1
              kc = pmeta%pos(3) + 1

              il = faces_gp(i ,jc,kc,idir)%l(1)
              jl = faces_gp(i ,jc,kc,idir)%l(2)
              kl = faces_gp(i ,jc,kc,idir)%l(3)
              iu = faces_gp(i ,jc,kc,idir)%u(1)
              ju = faces_gp(i ,jc,kc,idir)%u(2)
              ku = faces_gp(i ,jc,kc,idir)%u(3)

              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:jh,1:kh)

            case(2)

              ic = pmeta%pos(1) + 1
              kc = pmeta%pos(3) + 1

              il = faces_gp(ic,j ,kc,idir)%l(1)
              jl = faces_gp(ic,j ,kc,idir)%l(2)
              kl = faces_gp(ic,j ,kc,idir)%l(3)
              iu = faces_gp(ic,j ,kc,idir)%u(1)
              ju = faces_gp(ic,j ,kc,idir)%u(2)
              ku = faces_gp(ic,j ,kc,idir)%u(3)

              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ih,1:ng,1:kh)

            case(3)

              ic = pmeta%pos(1) + 1
              jc = pmeta%pos(2) + 1

              il = faces_gp(ic,jc,k ,idir)%l(1)
              jl = faces_gp(ic,jc,k ,idir)%l(2)
              kl = faces_gp(ic,jc,k ,idir)%l(3)
              iu = faces_gp(ic,jc,k ,idir)%u(1)
              ju = faces_gp(ic,jc,k ,idir)%u(2)
              ku = faces_gp(ic,jc,k ,idir)%u(3)

              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ih,1:jh,1:ng)

            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the face boundary update by prolongation
!
    call stop_timer(ifp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_face_prolong
#endif /* NDIMS == 3 */
!
!===============================================================================
!
!  DOMAIN EDGE BOUNDARY UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BOUNDARIES_EDGE_COPY:
! -------------------------------
!
!   Subroutine updates the edge boundaries from blocks on the same level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_edge_copy(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in, jn, kn
    use coordinates    , only : im, jm, km
    use coordinates    , only : edges_gc, edges_dc
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: is, js, ks
    integer :: it, jt, kt
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the edge boundary update by copying
!
    call start_timer(iec)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
#if NDIMS == 3
    kh = kn / 2
#endif /* NDIMS == 3 */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
#endif /* NDIMS == 3 */
        do j = 1, nsides
          do i = 1, nsides

! associate pneigh with the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%edges(i,j,idir)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%edges(i,j,k,idir)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at the same level
!
              if (pneigh%level == pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare region indices of the block and its neighbor for the edge boundary
! update
#if NDIMS == 2
                      il = edges_gc(i,j  ,idir)%l(1)
                      jl = edges_gc(i,j  ,idir)%l(2)
                      iu = edges_gc(i,j  ,idir)%u(1)
                      ju = edges_gc(i,j  ,idir)%u(2)

                      is = edges_dc(i,j  ,idir)%l(1)
                      js = edges_dc(i,j  ,idir)%l(2)
                      it = edges_dc(i,j  ,idir)%u(1)
                      jt = edges_dc(i,j  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      il = edges_gc(i,j,k,idir)%l(1)
                      jl = edges_gc(i,j,k,idir)%l(2)
                      kl = edges_gc(i,j,k,idir)%l(3)
                      iu = edges_gc(i,j,k,idir)%u(1)
                      ju = edges_gc(i,j,k,idir)%u(2)
                      ku = edges_gc(i,j,k,idir)%u(3)

                      is = edges_dc(i,j,k,idir)%l(1)
                      js = edges_dc(i,j,k,idir)%l(2)
                      ks = edges_dc(i,j,k,idir)%l(3)
                      it = edges_dc(i,j,k,idir)%u(1)
                      jt = edges_dc(i,j,k,idir)%u(2)
                      kt = edges_dc(i,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

! copy the corresponding edge region from the neighbor and insert it in
! the current data block
!
#if NDIMS == 2
                      pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                   &
                                         pneigh%data%q(1:nv,is:it,js:jt, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                   &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 5
          rtag = 16 * (sproc * nprocs + rproc) + 5

! allocate buffers for variable exchange
!
          select case(idir)
#if NDIMS == 2
          case(1)
            allocate(sbuf(scount,nv,ih,ng,km))
            allocate(rbuf(rcount,nv,ih,ng,km))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,km))
            allocate(rbuf(rcount,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          case(1)
            allocate(sbuf(scount,nv,ih,ng,ng))
            allocate(rbuf(rcount,nv,ih,ng,ng))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,ng))
            allocate(rbuf(rcount,nv,ng,jh,ng))
          case(3)
            allocate(sbuf(scount,nv,ng,ng,kh))
            allocate(rbuf(rcount,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare indices of the region for edge boundary update
!
#if NDIMS == 2
            is = edges_dc(i,j  ,idir)%l(1)
            js = edges_dc(i,j  ,idir)%l(2)
            it = edges_dc(i,j  ,idir)%u(1)
            jt = edges_dc(i,j  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            is = edges_dc(i,j,k,idir)%l(1)
            js = edges_dc(i,j,k,idir)%l(2)
            ks = edges_dc(i,j,k,idir)%l(3)
            it = edges_dc(i,j,k,idir)%u(1)
            jt = edges_dc(i,j,k,idir)%u(2)
            kt = edges_dc(i,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

! copy the corresponding edge region from the neighbor and insert it in
! the buffer
!
            select case(idir)
            case(1)
#if NDIMS == 2
              sbuf(l,1:nv,1:ih,1:ng,1:km) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              sbuf(l,1:nv,1:ih,1:ng,1:ng) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              sbuf(l,1:nv,1:ng,1:jh,1:km) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              sbuf(l,1:nv,1:ng,1:jh,1:ng) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              sbuf(l,1:nv,1:ng,1:ng,1:kh) =                                    &
                                         pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare indices of the region for the edge update
!
#if NDIMS == 2
            il = edges_gc(i,j  ,idir)%l(1)
            jl = edges_gc(i,j  ,idir)%l(2)
            iu = edges_gc(i,j  ,idir)%u(1)
            ju = edges_gc(i,j  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            il = edges_gc(i,j,k,idir)%l(1)
            jl = edges_gc(i,j,k,idir)%l(2)
            kl = edges_gc(i,j,k,idir)%l(3)
            iu = edges_gc(i,j,k,idir)%u(1)
            ju = edges_gc(i,j,k,idir)%u(2)
            ku = edges_gc(i,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

! update the corresponding edge region of the current block
!
            select case(idir)
            case(1)
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) = rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) = rbuf(l,1:nv,1:ng,1:jh,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:jh,1:ng)
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the edge boundary update by copying
!
    call stop_timer(iec)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_edge_copy
!
!===============================================================================
!
! subroutine BOUNDARIES_EDGE_RESTRICT:
! -----------------------------------
!
!   Subroutine updates the edge boundaries from blocks on higher level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_edge_restrict(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : edges_gr
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the edge boundary update by restriction
!
    call start_timer(ier)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
#if NDIMS == 3
    kh = kn / 2
#endif /* NDIMS == 3 */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
#endif /* NDIMS == 3 */
        do j = 1, nsides
          do i = 1, nsides

! assign pneigh to the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%edges(i,j,idir)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%edges(i,j,k,idir)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at higher level
!
              if (pneigh%level > pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare the region indices for edge boundary update
!
#if NDIMS == 2
                      il = edges_gr(i,j,  idir)%l(1)
                      jl = edges_gr(i,j,  idir)%l(2)
                      iu = edges_gr(i,j,  idir)%u(1)
                      ju = edges_gr(i,j,  idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      il = edges_gr(i,j,k,idir)%l(1)
                      jl = edges_gr(i,j,k,idir)%l(2)
                      kl = edges_gr(i,j,k,idir)%l(3)
                      iu = edges_gr(i,j,k,idir)%u(1)
                      ju = edges_gr(i,j,k,idir)%u(2)
                      ku = edges_gr(i,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor to the current
! data block
!
#if NDIMS == 2
                      call block_edge_restrict(idir, i, j, k                   &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      call block_edge_restrict(idir, i, j, k                   &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 6
          rtag = 16 * (sproc * nprocs + rproc) + 6

! allocate buffers for variable exchange
!
          select case(idir)
#if NDIMS == 2
          case(1)
            allocate(sbuf(scount,nv,ih,ng,km))
            allocate(rbuf(rcount,nv,ih,ng,km))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,km))
            allocate(rbuf(rcount,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          case(1)
            allocate(sbuf(scount,nv,ih,ng,ng))
            allocate(rbuf(rcount,nv,ih,ng,ng))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,ng))
            allocate(rbuf(rcount,nv,ng,jh,ng))
          case(3)
            allocate(sbuf(scount,nv,ng,ng,kh))
            allocate(rbuf(rcount,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor and insert it
! to the buffer
!
            select case(idir)
            case(1)
#if NDIMS == 2
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:ng))
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:ng))
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:kh))
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare the region indices for edge boundary update
!
#if NDIMS == 2
            il = edges_gr(i,j,  idir)%l(1)
            jl = edges_gr(i,j,  idir)%l(2)
            iu = edges_gr(i,j,  idir)%u(1)
            ju = edges_gr(i,j,  idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            il = edges_gr(i,j,k,idir)%l(1)
            jl = edges_gr(i,j,k,idir)%l(2)
            kl = edges_gr(i,j,k,idir)%l(3)
            iu = edges_gr(i,j,k,idir)%u(1)
            ju = edges_gr(i,j,k,idir)%u(2)
            ku = edges_gr(i,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
            select case(idir)
            case(1)
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:ng)
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the edge boundary update by restriction
!
    call stop_timer(ier)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_edge_restrict
!
!===============================================================================
!
! subroutine BOUNDARIES_EDGE_PROLONG:
! ----------------------------------
!
!   Subroutine updates the edge boundaries from blocks on lower level.
!
!   Arguments:
!
!     idir - the direction to be processed;
!
!===============================================================================
!
  subroutine boundaries_edge_prolong(idir)

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : edges_gp
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(in) :: idir

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ic, jc, kc
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the edge boundary update by prolongation
!
    call start_timer(iep)
#endif /* PROFILE */

! calculate the sizes
!
    ih = in + ng
    jh = jn + ng
#if NDIMS == 3
    kh = kn + ng
#endif /* NDIMS == 3 */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
        kc = k
#endif /* NDIMS == 3 */
        do j = 1, nsides
          jc = j
          do i = 1, nsides
            ic = i

! assign pneigh to the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%edges(i,j,idir)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%edges(i,j,k,idir)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor lays at lower level
!
              if (pneigh%level < pmeta%level) then

! process only blocks and neighbors which are marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare the region indices for edge boundary update
!
                      select case(idir)
                      case(1)

                        ic = pmeta%pos(1) + 1

#if NDIMS == 2
                        il = edges_gp(ic,j  ,idir)%l(1)
                        iu = edges_gp(ic,j  ,idir)%u(1)
                        jl = edges_gp(i ,j  ,idir)%l(2)
                        ju = edges_gp(i ,j  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        il = edges_gp(ic,j,k,idir)%l(1)
                        iu = edges_gp(ic,j,k,idir)%u(1)
                        jl = edges_gp(i ,j,k,idir)%l(2)
                        ju = edges_gp(i ,j,k,idir)%u(2)
                        kl = edges_gp(i ,j,k,idir)%l(3)
                        ku = edges_gp(i ,j,k,idir)%u(3)
#endif /* NDIMS == 3 */
                      case(2)

                        jc = pmeta%pos(2) + 1

#if NDIMS == 2
                        il = edges_gp(i,j   ,idir)%l(1)
                        iu = edges_gp(i,j   ,idir)%u(1)
                        jl = edges_gp(i,jc  ,idir)%l(2)
                        ju = edges_gp(i,jc  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        il = edges_gp(i,j ,k,idir)%l(1)
                        iu = edges_gp(i,j ,k,idir)%u(1)
                        jl = edges_gp(i,jc,k,idir)%l(2)
                        ju = edges_gp(i,jc,k,idir)%u(2)
                        kl = edges_gp(i,j ,k,idir)%l(3)
                        ku = edges_gp(i,j ,k,idir)%u(3)
                      case(3)

                        kc = pmeta%pos(3) + 1

                        il = edges_gp(i,j,k ,idir)%l(1)
                        iu = edges_gp(i,j,k ,idir)%u(1)
                        jl = edges_gp(i,j,k ,idir)%l(2)
                        ju = edges_gp(i,j,k ,idir)%u(2)
                        kl = edges_gp(i,j,kc,idir)%l(3)
                        ku = edges_gp(i,j,kc,idir)%u(3)
#endif /* NDIMS == 3 */
                      end select

! extract the corresponding edge region from the neighbor and insert it in
! the current data block
!
#if NDIMS == 2
                      call block_edge_prolong(idir, ic, jc, kc                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      call block_edge_prolong(idir, ic, jc, kc                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, idir, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at lower level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 7
          rtag = 16 * (sproc * nprocs + rproc) + 7

! allocate buffers for variable exchange
!
          select case(idir)
#if NDIMS == 2
          case(1)
            allocate(sbuf(scount,nv,ih,ng,km))
            allocate(rbuf(rcount,nv,ih,ng,km))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,km))
            allocate(rbuf(rcount,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          case(1)
            allocate(sbuf(scount,nv,ih,ng,ng))
            allocate(rbuf(rcount,nv,ih,ng,ng))
          case(2)
            allocate(sbuf(scount,nv,ng,jh,ng))
            allocate(rbuf(rcount,nv,ng,jh,ng))
          case(3)
            allocate(sbuf(scount,nv,ng,ng,kh))
            allocate(rbuf(rcount,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
          end select

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pmeta and pneigh to the associated blocks
!
            pmeta  => pinfo%block
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor and insert it
! to the buffer
!
            select case(idir)
            case(1)
              i = pmeta%pos(1) + 1
#if NDIMS == 2
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ih,1:ng,1:ng))
#endif /* NDIMS == 3 */
            case(2)
              j = pmeta%pos(2) + 1
#if NDIMS == 2
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:jh,1:ng))
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              k = pmeta%pos(3) + 1
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:kh))
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
            select case(idir)
            case(1)
              ic = pmeta%pos(1) + 1

#if NDIMS == 2
              il = edges_gp(ic,j  ,idir)%l(1)
              iu = edges_gp(ic,j  ,idir)%u(1)
              jl = edges_gp(i ,j  ,idir)%l(2)
              ju = edges_gp(i ,j  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              il = edges_gp(ic,j,k,idir)%l(1)
              iu = edges_gp(ic,j,k,idir)%u(1)
              jl = edges_gp(i ,j,k,idir)%l(2)
              ju = edges_gp(i ,j,k,idir)%u(2)
              kl = edges_gp(i ,j,k,idir)%l(3)
              ku = edges_gp(i ,j,k,idir)%u(3)
#endif /* NDIMS == 3 */

#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */

            case(2)

              jc = pmeta%pos(2) + 1

#if NDIMS == 2
              il = edges_gp(i,j   ,idir)%l(1)
              iu = edges_gp(i,j   ,idir)%u(1)
              jl = edges_gp(i,jc  ,idir)%l(2)
              ju = edges_gp(i,jc  ,idir)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              il = edges_gp(i,j ,k,idir)%l(1)
              iu = edges_gp(i,j ,k,idir)%u(1)
              jl = edges_gp(i,jc,k,idir)%l(2)
              ju = edges_gp(i,jc,k,idir)%u(2)
              kl = edges_gp(i,j ,k,idir)%l(3)
              ku = edges_gp(i,j ,k,idir)%u(3)
#endif /* NDIMS == 3 */

#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:ng)

            case(3)
              kc = pmeta%pos(3) + 1

              il = edges_gp(i,j,k ,idir)%l(1)
              iu = edges_gp(i,j,k ,idir)%u(1)
              jl = edges_gp(i,j,k ,idir)%l(2)
              ju = edges_gp(i,j,k ,idir)%u(2)
              kl = edges_gp(i,j,kc,idir)%l(3)
              ku = edges_gp(i,j,kc,idir)%u(3)

              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the edge boundary update by prolongation
!
    call stop_timer(iep)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_edge_prolong
!
!===============================================================================
!
!  DOMAIN CORNER BOUNDARY UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BOUNDARIES_CORNER_COPY:
! ---------------------------------
!
!   Subroutine updates the corner boundaries from blocks on the same level.
!
!
!===============================================================================
!
  subroutine boundaries_corner_copy()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im, jm, km
    use coordinates    , only : corners_gc, corners_dc
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: is, js, ks
    integer :: it, jt, kt
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the corner boundary update by copying
!
    call start_timer(icc)
#endif /* PROFILE */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
#endif /* NDIMS == 3 */
        do j = 1, nsides
          do i = 1, nsides

! assign pneigh to the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%corners(i,j)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%corners(i,j,k)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at the same level
!
              if (pneigh%level == pmeta%level) then

! skip if the block and its neighbor are not marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare region indices of the block and its neighbor for the corner boundary
! update
#if NDIMS == 2
                      il = corners_gc(i,j  )%l(1)
                      jl = corners_gc(i,j  )%l(2)
                      iu = corners_gc(i,j  )%u(1)
                      ju = corners_gc(i,j  )%u(2)

                      is = corners_dc(i,j  )%l(1)
                      js = corners_dc(i,j  )%l(2)
                      it = corners_dc(i,j  )%u(1)
                      jt = corners_dc(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      il = corners_gc(i,j,k)%l(1)
                      jl = corners_gc(i,j,k)%l(2)
                      kl = corners_gc(i,j,k)%l(3)
                      iu = corners_gc(i,j,k)%u(1)
                      ju = corners_gc(i,j,k)%u(2)
                      ku = corners_gc(i,j,k)%u(3)

                      is = corners_dc(i,j,k)%l(1)
                      js = corners_dc(i,j,k)%l(2)
                      ks = corners_dc(i,j,k)%l(3)
                      it = corners_dc(i,j,k)%u(1)
                      jt = corners_dc(i,j,k)%u(2)
                      kt = corners_dc(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! copy the corresponding corner region from the neighbor to the current
! data block
!
#if NDIMS == 2
                      pmeta%data%q(1:nv,il:iu,jl:ju, 1:km)                     &
                                       = pneigh%data%q(1:nv,is:it,js:jt, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku)                     &
                                       = pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! pneigh on the current process

                  else ! block and neighbor belong to different processes

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, -1, i, j, k)

                  end if ! block and neighbor belong to different processes
#endif /* MPI */

                end if ! pmeta and pneigh marked for update

              end if ! neighbor at the same level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 8
          rtag = 16 * (sproc * nprocs + rproc) + 8

! allocate buffers for variable exchange
!
#if NDIMS == 2
          allocate(sbuf(scount,nv,ng,ng,km))
          allocate(rbuf(rcount,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          allocate(sbuf(scount,nv,ng,ng,ng))
          allocate(rbuf(rcount,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare the corner region indices for the neighbor
!
#if NDIMS == 2
            is = corners_dc(i,j  )%l(1)
            js = corners_dc(i,j  )%l(2)
            it = corners_dc(i,j  )%u(1)
            jt = corners_dc(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            is = corners_dc(i,j,k)%l(1)
            js = corners_dc(i,j,k)%l(2)
            ks = corners_dc(i,j,k)%l(3)
            it = corners_dc(i,j,k)%u(1)
            jt = corners_dc(i,j,k)%u(2)
            kt = corners_dc(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! copy the corresponding corner region from the neighbor to the buffer
!
#if NDIMS == 2
            sbuf(l,1:nv,1:ng,1:ng,1:km) = pneigh%data%q(1:nv,is:it,js:jt, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            sbuf(l,1:nv,1:ng,1:ng,1:ng) = pneigh%data%q(1:nv,is:it,js:jt,ks:kt)
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare the corner region indices for the block
!
#if NDIMS == 2
            il = corners_gc(i,j  )%l(1)
            jl = corners_gc(i,j  )%l(2)
            iu = corners_gc(i,j  )%u(1)
            ju = corners_gc(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            il = corners_gc(i,j,k)%l(1)
            jl = corners_gc(i,j,k)%l(2)
            kl = corners_gc(i,j,k)%l(3)
            iu = corners_gc(i,j,k)%u(1)
            ju = corners_gc(i,j,k)%u(2)
            ku = corners_gc(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
#if NDIMS == 2
            pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) = rbuf(l,1:nv,1:ng,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:ng,1:ng)
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the corner boundary update by copying
!
    call stop_timer(icc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_corner_copy
!
!===============================================================================
!
! subroutine BOUNDARIES_CORNER_RESTRICT:
! -------------------------------------
!
!   Subroutine updates the corner boundaries from blocks on higher level.
!
!
!===============================================================================
!
  subroutine boundaries_corner_restrict()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : corners_gr
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the corner boundary update by restriction
!
    call start_timer(icr)
#endif /* PROFILE */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
#endif /* NDIMS == 3 */
        do j = 1, nsides
          do i = 1, nsides

! assign pneigh to the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%corners(i,j)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%corners(i,j,k)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor is at higher level
!
              if (pneigh%level > pmeta%level) then

! skip if the block and its neighbor are not marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare the region indices for corner boundary update
!
#if NDIMS == 2
                      il = corners_gr(i,j  )%l(1)
                      jl = corners_gr(i,j  )%l(2)
                      iu = corners_gr(i,j  )%u(1)
                      ju = corners_gr(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      il = corners_gr(i,j,k)%l(1)
                      jl = corners_gr(i,j,k)%l(2)
                      kl = corners_gr(i,j,k)%l(3)
                      iu = corners_gr(i,j,k)%u(1)
                      ju = corners_gr(i,j,k)%u(2)
                      ku = corners_gr(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! extract and restrict the corresponding corner region from the neighbor and
! insert it in the current data block
!
#if NDIMS == 2
                      call block_corner_restrict(i, j, k                       &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      call block_corner_restrict(i, j, k                       &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! block on the current processor

                  else ! block and neighbor on different processors

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, -1, i, j, k)

                  end if ! block and neighbor on different processors
#endif /* MPI */
                end if ! pmeta and pneigh marked for update

              end if ! neighbor at higher level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 9
          rtag = 16 * (sproc * nprocs + rproc) + 9

! allocate buffers for variable exchange
!
#if NDIMS == 2
          allocate(sbuf(scount,nv,ng,ng,km))
          allocate(rbuf(rcount,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          allocate(sbuf(scount,nv,ng,ng,ng))
          allocate(rbuf(rcount,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! restrict and extract the corresponding corner region from the neighbor and
! insert it to the buffer
!
#if NDIMS == 2
            call block_corner_restrict(i, j, k                                 &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
            call block_corner_restrict(i, j, k                                 &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:ng))
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare the region indices for corner boundary update
!
#if NDIMS == 2
            il = corners_gr(i,j  )%l(1)
            jl = corners_gr(i,j  )%l(2)
            iu = corners_gr(i,j  )%u(1)
            ju = corners_gr(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            il = corners_gr(i,j,k)%l(1)
            jl = corners_gr(i,j,k)%l(2)
            kl = corners_gr(i,j,k)%l(3)
            iu = corners_gr(i,j,k)%u(1)
            ju = corners_gr(i,j,k)%u(2)
            ku = corners_gr(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
#if NDIMS == 2
            pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) = rbuf(l,1:nv,1:ng,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:ng,1:ng)
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the corner boundary update by restriction
!
    call stop_timer(icr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_corner_restrict
!
!===============================================================================
!
! subroutine BOUNDARIES_CORNER_PROLONG:
! ------------------------------------
!
!   Subroutine updates the corner boundaries from blocks on lower level.
!
!
!===============================================================================
!
  subroutine boundaries_corner_prolong()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data, block_leaf
    use blocks         , only : list_meta, list_leaf
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : corners_gp
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc, nprocs, npairs, pairs
    use mpitools       , only : exchange_real_arrays
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
    type(block_leaf), pointer :: pleaf
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
#ifdef MPI
    integer :: sproc, scount, stag
    integer :: rproc, rcount, rtag
    integer :: l, p, iret

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: sbuf, rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the corner boundary update by prolongation
!
    call start_timer(icp)
#endif /* PROFILE */

#ifdef MPI
! prepare the block exchange structures
!
    call prepare_exchange_array()
#endif /* MPI */

! update boundaries between blocks on the same process
!
! associate pleaf with the first block on the leaf list
!
    pleaf => list_leaf

! scan all leaf meta blocks in the list
!
    do while(associated(pleaf))

! get the associated meta block
!
      pmeta => pleaf%meta

! scan over all block corners
!
#if NDIMS == 3
      do k = 1, nsides
#endif /* NDIMS == 3 */
        do j = 1, nsides
          do i = 1, nsides

! assign pneigh to the current neighbor
!
#if NDIMS == 2
            pneigh => pmeta%corners(i,j)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pneigh => pmeta%corners(i,j,k)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
            if (associated(pneigh)) then

! check if the neighbor lays at lower level
!
              if (pneigh%level < pmeta%level) then

! skip if the block and its neighbor are not marked for update
!
                if (pmeta%update .or. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                  if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                    if (pneigh%process == nproc) then
#endif /* MPI */

! prepare the region indices for corner boundary update
!
#if NDIMS == 2
                      il = corners_gp(i,j  )%l(1)
                      jl = corners_gp(i,j  )%l(2)
                      iu = corners_gp(i,j  )%u(1)
                      ju = corners_gp(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      il = corners_gp(i,j,k)%l(1)
                      jl = corners_gp(i,j,k)%l(2)
                      kl = corners_gp(i,j,k)%l(3)
                      iu = corners_gp(i,j,k)%u(1)
                      ju = corners_gp(i,j,k)%u(2)
                      ku = corners_gp(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! restrict and extract the corresponding corner region from the neighbor and
! insert it in the current data block
!
#if NDIMS == 2
                      call block_corner_prolong(i, j, k                        &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                      call block_corner_prolong(i, j, k                        &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                    end if ! block on the current processor

                  else ! block and neighbor on different processors

! append the block to the exchange list
!
                    call append_exchange_block(pmeta, pneigh, -1, i, j, k)

                  end if ! block and neighbor on different processors
#endif /* MPI */
                end if ! pmeta and pneigh marked for update

              end if ! neighbor at lower level

            end if ! neighbor associated

          end do ! i = 1, nsides
        end do ! j = 1, nsides
#if NDIMS == 3
      end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pleaf with the next leaf on the list
!
      pleaf => pleaf%next

    end do ! over leaf blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! process only pairs related to this process
!
      if (pairs(p,1) == nproc .or. pairs(p,2) == nproc) then

! get sending and receiving process identifiers (depending on pair member)
!
        if (pairs(p,1) == nproc) then
          sproc = pairs(p,1)
          rproc = pairs(p,2)
        end if
        if (pairs(p,2) == nproc) then
          sproc = pairs(p,2)
          rproc = pairs(p,1)
        end if

! get the number of blocks to exchange
!
        scount = bcount(sproc,rproc)
        rcount = bcount(rproc,sproc)

! process only pairs which have anything to exchange
!
        if ((scount + rcount) > 0) then

! prepare the tag for communication
!
          stag = 16 * (rproc * nprocs + sproc) + 10
          rtag = 16 * (sproc * nprocs + rproc) + 10

! allocate buffers for variable exchange
!
#if NDIMS == 2
          allocate(sbuf(scount,nv,ng,ng,km))
          allocate(rbuf(rcount,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          allocate(sbuf(scount,nv,ng,ng,ng))
          allocate(rbuf(rcount,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

!! PREPARE BLOCKS FOR SENDING
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(sproc,rproc)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign pneigh to the associated neighbor block
!
            pneigh => pinfo%neigh

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prolong the corresponding corner region from the neighbor and insert it in
! the buffer
!
#if NDIMS == 2
            call block_corner_prolong(i, j, k                                  &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
            call block_corner_prolong(i, j, k                                  &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        sbuf(l,1:nv,1:ng,1:ng,1:ng))
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

!! SEND PREPARED BLOCKS AND RECEIVE NEW ONES
!!
! exchange data
!
          call exchange_real_arrays(rproc, stag, size(sbuf), sbuf              &
                                  , rproc, rtag, size(rbuf), rbuf, iret)

!! PROCESS RECEIVED BLOCKS
!!
! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => barray(rproc,sproc)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! assign a pointer to the associated data block
!
            pmeta => pinfo%block

! get the corner coordinates
!
            i = pinfo%corner(1)
            j = pinfo%corner(2)
#if NDIMS == 3
            k = pinfo%corner(3)
#endif /* NDIMS == 3 */

! prepare the region indices for corner boundary update
!
#if NDIMS == 2
            il = corners_gp(i,j  )%l(1)
            jl = corners_gp(i,j  )%l(2)
            iu = corners_gp(i,j  )%u(1)
            ju = corners_gp(i,j  )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            il = corners_gp(i,j,k)%l(1)
            jl = corners_gp(i,j,k)%l(2)
            kl = corners_gp(i,j,k)%l(3)
            iu = corners_gp(i,j,k)%u(1)
            ju = corners_gp(i,j,k)%u(2)
            ku = corners_gp(i,j,k)%u(3)
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
#if NDIMS == 2
            pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) = rbuf(l,1:nv,1:ng,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) = rbuf(l,1:nv,1:ng,1:ng,1:ng)
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! deallocate data buffer
!
          deallocate(sbuf, rbuf)

        end if ! (scount + rcount) > 0

      end if ! pairs(p,1) == nproc || pairs(p,2) == nproc

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for the corner boundary update by prolongation
!
    call stop_timer(icp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundaries_corner_prolong
!
!===============================================================================
!
!  BLOCK SPECIFIC BOUNDARY SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BLOCK_BOUNDARY_SPECIFIC:
! ----------------------------------
!
!   Subroutine applies specific boundary conditions to the pointed data block.
!
!   Arguments:
!
!     nc         - the edge direction;
!     ic, jc, kc - the corner position;
!     t, dt      - time and time increment;
!     x, y, z    - the block coordinates;
!     qn         - the variable array;
!
!===============================================================================
!
  subroutine block_boundary_specific(ic, jc, kc, nc, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im , jm , km , ng
    use coordinates    , only : ib , jb , kb , ie , je , ke
    use coordinates    , only : ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ipr, ivx, ivy, ivz, ibx, iby, ibz, ibp
    use equations      , only : csnd2
    use error          , only : print_error, print_warning
    use gravity        , only : gravitational_acceleration
    use user_problem   , only : boundary_user_x, boundary_user_y               &
                              , boundary_user_z

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: ic, jc, kc
    integer                                     , intent(in)    :: nc
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(inout) :: x
    real(kind=8), dimension(1:jm)               , intent(inout) :: y
    real(kind=8), dimension(1:km)               , intent(inout) :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn

! local variables
!
    integer      :: i, il, iu, is, it, im1, ip1
    integer      :: j, jl, ju, js, jt, jm1, jp1
    integer      :: k, kl, ku, ks, kt, km1, kp1
    real(kind=8) :: dx, dy, dz, dxh, dyh, dzh, xi, yi, zi

! local vectors
!
    real(kind=8), dimension(3) :: ga
!
!-------------------------------------------------------------------------------
!
! apply specific boundaries depending on the direction
!
    select case(nc)
    case(1)

! prepare indices for the boundaries
!
      if (jc == 1) then
        jl = 1
        ju = jm / 2 - 1
      else
        jl = jm / 2
        ju = jm
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = 1
        ku = km / 2 - 1
      else
        kl = km / 2
        ku = km
      end if
#else /* NDIMS == 3 */
        kl = 1
        ku = km
#endif /* NDIMS == 3 */

! apply selected boundary condition
!
      select case(bnd_type(nc,ic))

! "open" boundary conditions
!
      case(bnd_open)

        if (ic == 1) then
          do i = ibl, 1, -1
            qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ib,jl:ju,kl:ku)
          end do
        else
          do i = ieu, im
            qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ie,jl:ju,kl:ku)
          end do
        end if

! "outflow" boundary conditions
!
      case(bnd_outflow)

        if (ic == 1) then
          do i = ibl, 1, -1
            qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ib,jl:ju,kl:ku)
            qn(ivx ,i,jl:ju,kl:ku) = min(0.0d+00, qn(ivx,ib,jl:ju,kl:ku))
          end do ! i = ibl, 1, -1
        else
          do i = ieu, im
            qn(1:nv,i,jl:ju,kl:ku) = qn(1:nv,ie,jl:ju,kl:ku)
            qn(ivx ,i,jl:ju,kl:ku) = max(0.0d+00, qn(ivx,ie,jl:ju,kl:ku))
          end do ! i = ieu, im
        end if

! "reflective" boundary conditions
!
      case(bnd_reflective)

        if (ic == 1) then
          do i = 1, ng
            it = ib  - i
            is = ibl + i

            qn(1:nv,it,jl:ju,kl:ku) =   qn(1:nv,is,jl:ju,kl:ku)
            qn(ivx ,it,jl:ju,kl:ku) = - qn(ivx ,is,jl:ju,kl:ku)
            if (ibx > 0) then
              qn(ibx ,it,jl:ju,kl:ku) = - qn(ibx ,is,jl:ju,kl:ku)
            end if
          end do
        else
          do i = 1, ng
            it = ie  + i
            is = ieu - i

            qn(1:nv,it,jl:ju,kl:ku) =   qn(1:nv,is,jl:ju,kl:ku)
            qn(ivx ,it,jl:ju,kl:ku) = - qn(ivx ,is,jl:ju,kl:ku)
            if (ibx > 0) then
              qn(ibx ,it,jl:ju,kl:ku) = - qn(ibx ,is,jl:ju,kl:ku)
            end if
          end do
        end if

! "gravity" or "hydrostatic" boundary conditions
!
      case(bnd_gravity)

        dx  = x(ib) - x(ibl)
        dxh = 0.5d+00 * dx

        if (ipr > 0) then
          if (ic == 1) then
            do i = ibl, 1, -1
              ip1 = i + 1
              xi  = x(i) + dxh
              do k = kl, ku
                do j = jl, ju
                  qn(1:nv,i,j,k) = qn(1:nv,ib,j,k)

                  call gravitational_acceleration(t, dt, xi, y(j), z(k), ga(:))

                  qn(ipr,i,j,k) = qn(ipr,ip1,j,k)                              &
                             - (qn(idn,ip1,j,k) + qn(idn,i,j,k)) * ga(1) * dxh
                end do
              end do
            end do
          else
            do i = ieu, im
              im1 = i - 1
              xi  = x(i) - dxh
              do k = kl, ku
                do j = jl, ju
                  qn(1:nv,i,j,k) = qn(1:nv,ie,j,k)

                  call gravitational_acceleration(t, dt, xi, y(j), z(k), ga(:))

                  qn(ipr,i,j,k) = qn(ipr,im1,j,k)                              &
                             + (qn(idn,im1,j,k) + qn(idn,i,j,k)) * ga(1) * dxh
                end do
              end do
            end do
          end if
        else
          if (ic == 1) then
            do i = ibl, 1, -1
              ip1 = i + 1
              xi  = x(i) + dxh
              do k = kl, ku
                do j = jl, ju
                  qn(1:nv,i,j,k) = qn(1:nv,ib,j,k)

                  call gravitational_acceleration(t, dt, xi, y(j), z(k), ga(:))

                  qn(idn,i,j,k) = qn(idn,ip1,j,k) * exp(- ga(1) * dx / csnd2)
                end do
              end do
            end do
          else
            do i = ieu, im
              im1 = i - 1
              xi  = x(i) - dxh
              do k = kl, ku
                do j = jl, ju
                  qn(1:nv,i,j,k) = qn(1:nv,ie,j,k)

                  call gravitational_acceleration(t, dt, xi, y(j), z(k), ga(:))

                  qn(idn,i,j,k) = qn(idn,im1,j,k) * exp(  ga(1) * dx / csnd2)
                end do
              end do
            end do
          end if
        end if

! user specific boundary conditions
!
      case(bnd_user)

        call boundary_user_x(ic, jl, ju, kl, ku                                &
                           , t, dt, x(1:im), y(1:jm), z(1:km)                  &
                           , qn(1:nv,1:im,1:jm,1:km))

! wrong boundary conditions
!
      case default

        if (ic == 1) then
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong left X boundary type!")
        else
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong right X boundary type!")
        end if

      end select

    case(2)

! prepare indices for the boundaries
!
      if (ic == 1) then
        il = 1
        iu = im / 2 - 1
      else
        il = im / 2
        iu = im
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = 1
        ku = km / 2 - 1
      else
        kl = km / 2
        ku = km
      end if
#else /* NDIMS == 3 */
        kl = 1
        ku = km
#endif /* NDIMS == 3 */

! apply selected boundary condition
!
      select case(bnd_type(nc,jc))

! "open" boundary conditions
!
      case(bnd_open)

        if (jc == 1) then
          do j = jbl, 1, -1
            qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,jb,kl:ku)
          end do
        else
          do j = jeu, jm
            qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,je,kl:ku)
          end do
        end if

! "outflow" boundary conditions
!
      case(bnd_outflow)

        if (jc == 1) then
          do j = jbl, 1, -1
            qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,jb,kl:ku)
            qn(ivy ,il:iu,j,kl:ku) = min(0.0d+00, qn(ivy,il:iu,jb,kl:ku))
          end do ! j = jbl, 1, -1
        else
          do j = jeu, jm
            qn(1:nv,il:iu,j,kl:ku) = qn(1:nv,il:iu,je,kl:ku)
            qn(ivy ,il:iu,j,kl:ku) = max(0.0d+00, qn(ivy,il:iu,je,kl:ku))
          end do ! j = jeu, jm
        end if

! "reflective" boundary conditions
!
      case(bnd_reflective)

        if (jc == 1) then
          do j = 1, ng
            jt = jb  - j
            js = jbl + j

            qn(1:nv,il:iu,jt,kl:ku) =   qn(1:nv,il:iu,js,kl:ku)
            qn(ivy ,il:iu,jt,kl:ku) = - qn(ivy ,il:iu,js,kl:ku)
            if (iby > 0) then
              qn(iby ,il:iu,jt,kl:ku) = - qn(iby ,il:iu,js,kl:ku)
            end if
          end do
        else
          do j = 1, ng
            jt = je  + j
            js = jeu - j

            qn(1:nv,il:iu,jt,kl:ku) =   qn(1:nv,il:iu,js,kl:ku)
            qn(ivy ,il:iu,jt,kl:ku) = - qn(ivy ,il:iu,js,kl:ku)
            if (iby > 0) then
              qn(iby ,il:iu,jt,kl:ku) = - qn(iby ,il:iu,js,kl:ku)
            end if
          end do
        end if

! "gravity" or "hydrostatic" boundary conditions
!
      case(bnd_gravity)

        dy  = y(jb) - y(jbl)
        dyh = 0.5d+00 * dy

        if (ipr > 0) then
          if (jc == 1) then
            do j = jbl, 1, -1
              jp1 = j + 1
              yi  = y(j) + dyh
              do k = kl, ku
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,jb,k)

                  call gravitational_acceleration(t, dt, x(i), yi, z(k), ga(:))

                  qn(ipr,i,j,k) = qn(ipr,i,jp1,k)                              &
                             - (qn(idn,i,jp1,k) + qn(idn,i,j,k)) * ga(2) * dyh
                end do
              end do
            end do
          else
            do j = jeu, jm
              jm1 = j - 1
              yi  = y(j) - dyh
              do k = kl, ku
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,je,k)

                  call gravitational_acceleration(t, dt, x(i), yi, z(k), ga(:))

                  qn(ipr,i,j,k) = qn(ipr,i,jm1,k)                              &
                             + (qn(idn,i,jm1,k) + qn(idn,i,j,k)) * ga(2) * dyh
                end do
              end do
            end do
          end if
        else
          if (jc == 1) then
            do j = jbl, 1, -1
              jp1 = j + 1
              yi  = y(j) + dyh
              do k = kl, ku
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,jb,k)

                  call gravitational_acceleration(t, dt, x(i), yi, z(k), ga(:))

                  qn(idn,i,j,k) = qn(idn,i,jp1,k) * exp(- ga(2) * dy / csnd2)
                end do
              end do
            end do
          else
            do j = jeu, jm
              jm1 = j - 1
              yi  = y(j) - dyh
              do k = kl, ku
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,je,k)

                  call gravitational_acceleration(t, dt, x(i), yi, z(k), ga(:))

                  qn(idn,i,j,k) = qn(idn,i,jm1,k) * exp(  ga(2) * dy / csnd2)
                end do
              end do
            end do
          end if
        end if

! user specific boundary conditions
!
      case(bnd_user)

        call boundary_user_y(jc, il, iu, kl, ku                                &
                           , t, dt, x(1:im), y(1:jm), z(1:km)                  &
                           , qn(1:nv,1:im,1:jm,1:km))

! wrong boundary conditions
!
      case default

        if (jc == 1) then
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong left Y boundary type!")
        else
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong right Y boundary type!")
        end if

      end select

#if NDIMS == 3
    case(3)

! prepare indices for the boundaries
!
      if (ic == 1) then
        il = 1
        iu = im / 2 - 1
      else
        il = im / 2
        iu = im
      end if
      if (jc == 1) then
        jl = 1
        ju = jm / 2 - 1
      else
        jl = jm / 2
        ju = jm
      end if

! apply selected boundary condition
!
      select case(bnd_type(nc,kc))

! "open" boundary conditions
!
      case(bnd_open)

        if (kc == 1) then
          do k = kbl, 1, -1
            qn(1:nv,il:iu,jl:ju,k) = qn(1:nv,il:iu,jl:ju,kb)
          end do
        else
          do k = keu, km
            qn(1:nv,il:iu,jl:ju,k) = qn(1:nv,il:iu,jl:ju,ke)
          end do
        end if

! "outflow" boundary conditions
!
      case(bnd_outflow)

        if (kc == 1) then
          do k = kbl, 1, -1
            qn(1:nv,il:iu,jl:ju,k) = qn(1:nv,il:iu,jl:ju,kb)
            qn(ivz ,il:iu,jl:ju,k) = min(0.0d+00, qn(ivz,il:iu,jl:ju,kb))
          end do ! k = kbl, 1, -1
        else
          do k = keu, km
            qn(1:nv,il:iu,jl:ju,k) = qn(1:nv,il:iu,jl:ju,ke)
            qn(ivz ,il:iu,jl:ju,k) = max(0.0d+00, qn(ivz,il:iu,jl:ju,ke))
          end do ! k = keu, km
        end if

! "reflective" boundary conditions
!
      case(bnd_reflective)

        if (kc == 1) then
          do k = 1, ng
            kt = kb  - k
            ks = kbl + k

            qn(1:nv,il:iu,jl:ju,kt) =   qn(1:nv,il:iu,jl:ju,ks)
            qn(ivz ,il:iu,jl:ju,kt) = - qn(ivz ,il:iu,jl:ju,ks)
            if (ibz > 0) then
              qn(ibz ,il:iu,jl:ju,kt) = - qn(ibz ,il:iu,jl:ju,ks)
            end if
          end do
        else
          do k = 1, ng
            kt = ke  + k
            ks = keu - k

            qn(1:nv,il:iu,jl:ju,kt) =   qn(1:nv,il:iu,jl:ju,ks)
            qn(ivz ,il:iu,jl:ju,kt) = - qn(ivz ,il:iu,jl:ju,ks)
            if (ibz > 0) then
              qn(ibz ,il:iu,jl:ju,kt) = - qn(ibz ,il:iu,jl:ju,ks)
            end if
          end do
        end if

! "gravity" or "hydrostatic" boundary conditions
!
      case(bnd_gravity)

        dz  = z(kb) - z(kbl)
        dzh = 0.5d+00 * dz

        if (ipr > 0) then
          if (kc == 1) then
            do k = kbl, 1, -1
              kp1 = k + 1
              zi  = z(k) + dzh
              do j = jl, ju
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,j,kb)

                  call gravitational_acceleration(t, dt, x(i), y(j), zi, ga(:))

                  qn(ipr,i,j,k) = qn(ipr,i,j,kp1)                              &
                             - (qn(idn,i,j,kp1) + qn(idn,i,j,k)) * ga(3) * dzh
                end do
              end do
            end do
          else
            do k = keu, km
              km1 = k - 1
              zi  = z(k) - dzh
              do j = jl, ju
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,j,ke)

                  call gravitational_acceleration(t, dt, x(i), y(j), zi, ga(:))

                  qn(ipr,i,j,k) = qn(ipr,i,j,km1)                              &
                             + (qn(idn,i,j,km1) + qn(idn,i,j,k)) * ga(3) * dzh
                end do
              end do
            end do
          end if
        else
          if (kc == 1) then
            do k = kbl, 1, -1
              kp1 = k + 1
              zi  = z(k) + dzh
              do j = jl, ju
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,j,kb)

                  call gravitational_acceleration(t, dt, x(i), y(j), zi, ga(:))

                  qn(idn,i,j,k) = qn(idn,i,j,kp1) * exp(- ga(3) * dz / csnd2)
                end do
              end do
            end do
          else
            do k = keu, km
              km1 = k - 1
              zi  = z(k) - dzh
              do j = jl, ju
                do i = il, iu
                  qn(1:nv,i,j,k) = qn(1:nv,i,j,ke)

                  call gravitational_acceleration(t, dt, x(i), y(j), zi, ga(:))

                  qn(idn,i,j,k) = qn(idn,i,j,km1) * exp(  ga(3) * dz / csnd2)
                end do
              end do
            end do
          end if
        end if

! user specific boundary conditions
!
      case(bnd_user)

        call boundary_user_z(kc, il, iu, jl, ju                                &
                           , t, dt, x(1:im), y(1:jm), z(1:km)                  &
                           , qn(1:nv,1:im,1:jm,1:km))

! wrong boundary conditions
!
      case default

        if (kc == 1) then
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong left Z boundary type!")
        else
          call print_error("boundaries:boundary_specific()"                    &
                                              , "Wrong right Z boundary type!")
        end if

      end select

#endif /* NDIMS == 3 */
    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_boundary_specific
#if NDIMS == 3
!
!===============================================================================
!
!  BLOCK FACE UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BLOCK_FACE_RESTRICT:
! ------------------------------
!
!   Subroutine returns the face boundary region restricted from the provided
!   input variable array.
!
!   Arguments:
!
!     nc         - the face direction;
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output face boundary array;
!
!===============================================================================
!
  subroutine block_face_restrict(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nd
    use coordinates    , only : in
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : faces_dr
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local variables
!
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: ip, jp, kp
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! prepare indices for the face region
!
    il = faces_dr(ic,jc,kc,nc)%l(1)
    jl = faces_dr(ic,jc,kc,nc)%l(2)
    kl = faces_dr(ic,jc,kc,nc)%l(3)
    ip = il + 1
    jp = jl + 1
    kp = kl + 1
    iu = faces_dr(ic,jc,kc,nc)%u(1)
    ju = faces_dr(ic,jc,kc,nc)%u(2)
    ku = faces_dr(ic,jc,kc,nc)%u(3)

! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half sizes
!
      jh = jn / 2
      kh = kn / 2

! restrict the face region to the output array
!
      qb(1:nv,1:ng,1:jh,1:kh) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))

    case(2)

! calculate half sizes
!
      ih = in / 2
      kh = kn / 2

! restrict the face region to the output array
!
      qb(1:nv,1:ih,1:ng,1:kh) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))

    case(3)

! calculate half sizes
!
      ih = in / 2
      jh = jn / 2

! restrict the face region to the output array
!
      qb(1:nv,1:ih,1:jh,1:ng) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))

    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_face_restrict
!
!===============================================================================
!
! subroutine BLOCK_FACE_PROLONG:
! -----------------------------
!
!   Subroutine returns the face boundary region prolongated from the provided
!   input variable array.
!
!   Arguments:
!
!     nc         - the face direction;
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output face boundary array;
!
!===============================================================================
!
  subroutine block_face_prolong(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nh
    use coordinates    , only : in
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : faces_dp
    use equations      , only : nv, idn, ipr
    use interpolations , only : limiter_prol

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local variables
!
    integer      :: i, j, k, p
    integer      :: ih, jh, kh
    integer      :: il, jl, kl
    integer      :: iu, ju, ku
    integer      :: is, js, ks
    integer      :: it, jt, kt
    integer      :: im1, jm1, km1
    integer      :: ip1, jp1, kp1
    real(kind=8) :: dql, dqr
    real(kind=8) :: dq1, dq2, dq3, dq4

! local arrays
!
    real(kind=8), dimension(3) :: dq
!
!-------------------------------------------------------------------------------
!
! prepare subregion indices
!
    il = faces_dp(ic,jc,kc,nc)%l(1)
    jl = faces_dp(ic,jc,kc,nc)%l(2)
    kl = faces_dp(ic,jc,kc,nc)%l(3)
    iu = faces_dp(ic,jc,kc,nc)%u(1)
    ju = faces_dp(ic,jc,kc,nc)%u(2)
    ku = faces_dp(ic,jc,kc,nc)%u(3)

! iterate over all face region cells
!
    do k = kl, ku
      km1 = k - 1
      kp1 = k + 1
      ks = 2 * (k - kl) + 1
      kt = ks + 1
      do j = jl, ju
        jm1 = j - 1
        jp1 = j + 1
        js = 2 * (j - jl) + 1
        jt = js + 1
        do i = il, iu
          im1 = i - 1
          ip1 = i + 1
          is = 2 * (i - il) + 1
          it = is + 1

! iterate over all variables
!
          do p = 1, nv

! calculate limited derivatives in all directions
!
            dql   = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr   = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dq(1) = limiter_prol(0.5d+00, dql, dqr)

            dql   = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr   = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dq(2) = limiter_prol(0.5d+00, dql, dqr)

            dql   = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr   = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dq(3) = limiter_prol(0.5d+00, dql, dqr)

            if (p == idn .or. p == ipr) then
              do while (qn(p,i,j,k) <= sum(abs(dq(1:NDIMS))))
                dq(:) = 0.5d+00 * dq(:)
              end do
            end if

            dq(:) = 0.5d+00 * dq(:)

! calculate the derivative terms
!
            dq1 = dq(1) + dq(2) + dq(3)
            dq2 = dq(1) - dq(2) - dq(3)
            dq3 = dq(1) - dq(2) + dq(3)
            dq4 = dq(1) + dq(2) - dq(3)

! prolong the face region to the output array
!
            qb(p,is,js,ks) = qn(p,i,j,k) - dq1
            qb(p,it,js,ks) = qn(p,i,j,k) + dq2
            qb(p,is,jt,ks) = qn(p,i,j,k) - dq3
            qb(p,it,jt,ks) = qn(p,i,j,k) + dq4
            qb(p,is,js,kt) = qn(p,i,j,k) - dq4
            qb(p,it,js,kt) = qn(p,i,j,k) + dq3
            qb(p,is,jt,kt) = qn(p,i,j,k) - dq2
            qb(p,it,jt,kt) = qn(p,i,j,k) + dq1

          end do ! q = 1, nv

        end do ! i = il, iu
      end do ! j = jl, ju
    end do ! k = kl, ku

!-------------------------------------------------------------------------------
!
  end subroutine block_face_prolong
#endif /* NDIMS == 3 */
!
!===============================================================================
!
!  BLOCK EDGE UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BLOCK_EDGE_RESTRICT:
! ------------------------------
!
!   Subroutine returns the edge boundary region by restricting the corresponding
!   region from the provided input variable array.
!
!   Arguments:
!
!     nc         - the edge direction;
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output edge boundary array;
!
!===============================================================================
!
  subroutine block_edge_restrict(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nd
    use coordinates    , only : in
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : edges_dr
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local variables
!
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: ip, jp, kp
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! prepare indices for the edge region
!
#if NDIMS == 2
    il = edges_dr(ic,jc,nc)%l(1)
    jl = edges_dr(ic,jc,nc)%l(2)
    ip = il + 1
    jp = jl + 1
    iu = edges_dr(ic,jc,nc)%u(1)
    ju = edges_dr(ic,jc,nc)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
    il = edges_dr(ic,jc,kc,nc)%l(1)
    jl = edges_dr(ic,jc,kc,nc)%l(2)
    kl = edges_dr(ic,jc,kc,nc)%l(3)
    ip = il + 1
    jp = jl + 1
    kp = kl + 1
    iu = edges_dr(ic,jc,kc,nc)%u(1)
    ju = edges_dr(ic,jc,kc,nc)%u(2)
    ku = edges_dr(ic,jc,kc,nc)%u(3)
#endif /* NDIMS == 3 */

! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half size
!
      ih = in / 2

! restrict the edge region to the output array
!
#if NDIMS == 2
      qb(1:nv,1:ih,1:ng,1:km) =                                                &
                           2.50d-01 *  ((qn(1:nv,il:iu:2,jl:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2, 1:km  ))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2, 1:km  )))
#endif /* NDIMS == 2 */
#if NDIMS == 3
      qb(1:nv,1:ih,1:ng,1:ng) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

    case(2)

! calculate half size
!
      jh = jn / 2

! restrict the edge region to the output array
!
#if NDIMS == 2
      qb(1:nv,1:ng,1:jh,1:km) =                                                &
                           2.50d-01 *  ((qn(1:nv,il:iu:2,jl:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2, 1:km  ))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2, 1:km  )))
#endif /* NDIMS == 2 */
#if NDIMS == 3
      qb(1:nv,1:ng,1:jh,1:ng) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

! calculate half size
!
      kh = kn / 2

! restrict the edge region to the output array
!
      qb(1:nv,1:ng,1:ng,1:kh) =                                                &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_edge_restrict
!
!===============================================================================
!
! subroutine BLOCK_EDGE_PROLONG:
! -----------------------------
!
!   Subroutine returns the edge boundary region by prolongating
!   the corresponding region from the provided input variable array.
!
!   Arguments:
!
!     nc         - the edge direction;
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output edge boundary array;
!
!===============================================================================
!
  subroutine block_edge_prolong(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nh
    use coordinates    , only : in
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : edges_dp
    use equations      , only : nv, idn, ipr
    use interpolations , only : limiter_prol

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local variables
!
    integer      :: i, j, k, p
    integer      :: il, jl, kl
    integer      :: iu, ju, ku
    integer      :: is, js, ks
    integer      :: it, jt, kt
    integer      :: im1, jm1, km1
    integer      :: ip1, jp1, kp1
    real(kind=8) :: dql, dqr
    real(kind=8) :: dq1, dq2, dq3, dq4

! local arrays
!
    real(kind=8), dimension(NDIMS) :: dq
!
!-------------------------------------------------------------------------------
!
! prepare indices for the edge region
!
#if NDIMS == 2
    il = edges_dp(ic,jc   ,nc)%l(1)
    jl = edges_dp(ic,jc   ,nc)%l(2)
    iu = edges_dp(ic,jc   ,nc)%u(1)
    ju = edges_dp(ic,jc   ,nc)%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
    il = edges_dp(ic,jc,kc,nc)%l(1)
    jl = edges_dp(ic,jc,kc,nc)%l(2)
    kl = edges_dp(ic,jc,kc,nc)%l(3)
    iu = edges_dp(ic,jc,kc,nc)%u(1)
    ju = edges_dp(ic,jc,kc,nc)%u(2)
    ku = edges_dp(ic,jc,kc,nc)%u(3)
#endif /* NDIMS == 3 */

! iterate over all edge region cells
!
#if NDIMS == 2
    do k = 1, km
      kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = kl, ku
      km1 = k - 1
      kp1 = k + 1
      ks = 2 * (k - kl) + 1
      kt = ks + 1
#endif /* NDIMS == 3 */
      do j = jl, ju
        jm1 = j - 1
        jp1 = j + 1
        js = 2 * (j - jl) + 1
        jt = js + 1
        do i = il, iu
          im1 = i - 1
          ip1 = i + 1
          is = 2 * (i - il) + 1
          it = is + 1

! iterate over all variables
!
          do p = 1, nv

! calculate limited derivatives in all directions
!
            dql   = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr   = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dq(1) = limiter_prol(0.5d+00, dql, dqr)

            dql   = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr   = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dq(2) = limiter_prol(0.5d+00, dql, dqr)

#if NDIMS == 3
            dql   = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr   = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dq(3) = limiter_prol(0.5d+00, dql, dqr)
#endif /* NDIMS == 3 */

            if (p == idn .or. p == ipr) then
              do while (qn(p,i,j,k) <= sum(abs(dq(1:NDIMS))))
                dq(:) = 0.5d+00 * dq(:)
              end do
            end if

            dq(:) = 0.5d+00 * dq(:)

#if NDIMS == 2
! calculate the derivative terms
!
            dq1 = dq(1) + dq(2)
            dq2 = dq(1) - dq(2)

! prolong the edge region to the output array
!
            qb(p,is,js,k ) = qn(p,i,j,k) - dq1
            qb(p,it,js,k ) = qn(p,i,j,k) + dq2
            qb(p,is,jt,k ) = qn(p,i,j,k) - dq2
            qb(p,it,jt,k ) = qn(p,i,j,k) + dq1
#endif /* NDIMS == 2 */
#if NDIMS == 3
! calculate the derivative terms
!
            dq1 = dq(1) + dq(2) + dq(3)
            dq2 = dq(1) - dq(2) - dq(3)
            dq3 = dq(1) - dq(2) + dq(3)
            dq4 = dq(1) + dq(2) - dq(3)

! prolong the edge region to the output array
!
            qb(p,is,js,ks) = qn(p,i,j,k) - dq1
            qb(p,it,js,ks) = qn(p,i,j,k) + dq2
            qb(p,is,jt,ks) = qn(p,i,j,k) - dq3
            qb(p,it,jt,ks) = qn(p,i,j,k) + dq4
            qb(p,is,js,kt) = qn(p,i,j,k) - dq4
            qb(p,it,js,kt) = qn(p,i,j,k) + dq3
            qb(p,is,jt,kt) = qn(p,i,j,k) - dq2
            qb(p,it,jt,kt) = qn(p,i,j,k) + dq1
#endif /* NDIMS == 3 */

          end do ! q = 1, nv

        end do ! i = il, iu
      end do ! j = jl, ju
    end do ! k = kl, ku

!-------------------------------------------------------------------------------
!
  end subroutine block_edge_prolong
!
!===============================================================================
!
!  BLOCK CORNER UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BLOCK_CORNER_RESTRICT:
! --------------------------------
!
!   Subroutine returns the corner boundary region by restricting
!   the corresponding region from the provided input variable array.
!
!   Arguments:
!
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output corner boundary array;
!
!===============================================================================
!
  subroutine block_corner_restrict(ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nd
    use coordinates    , only : im , jm , km
    use coordinates    , only : corners_dr
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
#if NDIMS == 2
    real(kind=8), dimension(1:nv,1:ng,1:ng,1:km), intent(out) :: qb
#endif /* NDIMS == 2 */
#if NDIMS == 3
    real(kind=8), dimension(1:nv,1:ng,1:ng,1:ng), intent(out) :: qb
#endif /* NDIMS == 3 */

! local variables
!
    integer :: il, jl, kl
    integer :: ip, jp, kp
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! prepare indices for the corner region
!
#if NDIMS == 2
    il = corners_dr(ic,jc   )%l(1)
    jl = corners_dr(ic,jc   )%l(2)
    ip = il + 1
    jp = jl + 1
    iu = corners_dr(ic,jc   )%u(1)
    ju = corners_dr(ic,jc   )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
    il = corners_dr(ic,jc,kc)%l(1)
    jl = corners_dr(ic,jc,kc)%l(2)
    kl = corners_dr(ic,jc,kc)%l(3)
    ip = il + 1
    jp = jl + 1
    kp = kl + 1
    iu = corners_dr(ic,jc,kc)%u(1)
    ju = corners_dr(ic,jc,kc)%u(2)
    ku = corners_dr(ic,jc,kc)%u(3)
#endif /* NDIMS == 3 */

! restrict the corner region to the output array
!
#if NDIMS == 2
    qb(1:nv,1:ng,1:ng,1:km) =                                                  &
                           2.50d-01 *  ((qn(1:nv,il:iu:2,jl:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2, 1:km  ))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2, 1:km  )      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2, 1:km  )))
#endif /* NDIMS == 2 */
#if NDIMS == 3
    qb(1:nv,1:ng,1:ng,1:ng) =                                                  &
                           1.25d-01 * (((qn(1:nv,il:iu:2,jl:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kp:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jl:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jp:ju:2,kl:ku:2)))    &
                                    +  ((qn(1:nv,il:iu:2,jp:ju:2,kp:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kl:ku:2))     &
                                    +   (qn(1:nv,il:iu:2,jp:ju:2,kl:ku:2)      &
                                    +    qn(1:nv,ip:iu:2,jl:ju:2,kp:ku:2))))
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine block_corner_restrict
!
!===============================================================================
!
! subroutine BLOCK_CORNER_PROLONG:
! -------------------------------
!
!   Subroutine returns the corner boundary region by prolongating
!   the corresponding region from the provided input variable array.
!
!   Arguments:
!
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output corner boundary array;
!
!===============================================================================
!
  subroutine block_corner_prolong(ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng, nh
    use coordinates    , only : im , jm , km
    use coordinates    , only : corners_dp
    use equations      , only : nv, idn, ipr
    use interpolations , only : limiter_prol

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
#if NDIMS == 2
    real(kind=8), dimension(1:nv,1:ng,1:ng,1:km), intent(out) :: qb
#endif /* NDIMS == 2 */
#if NDIMS == 3
    real(kind=8), dimension(1:nv,1:ng,1:ng,1:ng), intent(out) :: qb
#endif /* NDIMS == 3 */

! local variables
!
    integer      :: i, j, k, p
    integer      :: il, jl, kl
    integer      :: iu, ju, ku
    integer      :: is, js, ks
    integer      :: it, jt, kt
    integer      :: im1, jm1, km1
    integer      :: ip1, jp1, kp1
    real(kind=8) :: dql, dqr
    real(kind=8) :: dq1, dq2, dq3, dq4

! local arrays
!
    real(kind=8), dimension(NDIMS) :: dq
!
!-------------------------------------------------------------------------------
!
! prepare indices for the corner region
!
#if NDIMS == 2
    il = corners_dp(ic,jc   )%l(1)
    jl = corners_dp(ic,jc   )%l(2)
    iu = corners_dp(ic,jc   )%u(1)
    ju = corners_dp(ic,jc   )%u(2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
    il = corners_dp(ic,jc,kc)%l(1)
    jl = corners_dp(ic,jc,kc)%l(2)
    kl = corners_dp(ic,jc,kc)%l(3)
    iu = corners_dp(ic,jc,kc)%u(1)
    ju = corners_dp(ic,jc,kc)%u(2)
    ku = corners_dp(ic,jc,kc)%u(3)
#endif /* NDIMS == 3 */

! iterate over all corner region cells
!
#if NDIMS == 2
    do k = 1, km
      kt = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = kl, ku
      km1 = k - 1
      kp1 = k + 1
      ks = 2 * (k - kl) + 1
      kt = ks + 1
#endif /* NDIMS == 3 */
      do j = jl, ju
        jm1 = j - 1
        jp1 = j + 1
        js = 2 * (j - jl) + 1
        jt = js + 1
        do i = il, iu
          im1 = i - 1
          ip1 = i + 1
          is = 2 * (i - il) + 1
          it = is + 1

! iterate over all variables
!
          do p = 1, nv

! calculate limited derivatives in all directions
!
            dql   = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr   = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dq(1) = limiter_prol(0.5d+00, dql, dqr)

            dql   = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr   = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dq(2) = limiter_prol(0.5d+00, dql, dqr)

#if NDIMS == 3
            dql   = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr   = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dq(3) = limiter_prol(0.5d+00, dql, dqr)
#endif /* NDIMS == 3 */

            if (p == idn .or. p == ipr) then
              do while (qn(p,i,j,k) <= sum(abs(dq(1:NDIMS))))
                dq(:) = 0.5d+00 * dq(:)
              end do
            end if

            dq(:) = 0.5d+00 * dq(:)

#if NDIMS == 2
! calculate the derivative terms
!
            dq1 = dq(1) + dq(2)
            dq2 = dq(1) - dq(2)

! prolong the corner region to the output array
!
            qb(p,is,js,k ) = qn(p,i,j,k) - dq1
            qb(p,it,js,k ) = qn(p,i,j,k) + dq2
            qb(p,is,jt,k ) = qn(p,i,j,k) - dq2
            qb(p,it,jt,k ) = qn(p,i,j,k) + dq1
#endif /* NDIMS == 2 */
#if NDIMS == 3
! calculate the derivative terms
!
            dq1 = dq(1) + dq(2) + dq(3)
            dq2 = dq(1) - dq(2) - dq(3)
            dq3 = dq(1) - dq(2) + dq(3)
            dq4 = dq(1) + dq(2) - dq(3)

! prolong the corner region to the output array
!
            qb(p,is,js,ks) = qn(p,i,j,k) - dq1
            qb(p,it,js,ks) = qn(p,i,j,k) + dq2
            qb(p,is,jt,ks) = qn(p,i,j,k) - dq3
            qb(p,it,jt,ks) = qn(p,i,j,k) + dq4
            qb(p,is,js,kt) = qn(p,i,j,k) - dq4
            qb(p,it,js,kt) = qn(p,i,j,k) + dq3
            qb(p,is,jt,kt) = qn(p,i,j,k) - dq2
            qb(p,it,jt,kt) = qn(p,i,j,k) + dq1
#endif /* NDIMS == 3 */

          end do ! q = 1, nv

        end do ! i = il, iu
      end do ! j = jl, ju
    end do ! k = kl, ku

!-------------------------------------------------------------------------------
!
  end subroutine block_corner_prolong
!
!===============================================================================
!
!  BLOCK FLUX UPDATE SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine BLOCK_UPDATE_FLUX:
! ----------------------------
!
!   Subroutine updates the boundary flux from the provided flux array.
!
!   Arguments:
!
!     nc         - the edge direction;
!     ic, jc, kc - the corner position;
!     fn         - the correcting flux array;
!     fb         - the corrected flux array;
!
!===============================================================================
!
  subroutine block_update_flux(nc, ic, jc, kc, fn, fb)

! import external procedures and variables
!
    use blocks         , only : block_data
    use coordinates    , only : in, jn, kn
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                       , intent(in)    :: ic, jc, kc
    integer                       , intent(in)    :: nc
    real(kind=8), dimension(:,:,:), intent(in)    :: fn
    real(kind=8), dimension(:,:,:), intent(inout) :: fb
!
!-------------------------------------------------------------------------------
!
! update fluxes for each direction separately
!
    select case(nc)

! X direction
!
    case(1)

#if NDIMS == 2
! average fluxes from higher level neighbor
!
      fb(1:nv,:,:) = (fn(1:nv,1:jn:2,1:kn) + fn(1:nv,2:jn:2,1:kn)) / 2.0d+00
#endif /* NDIMS == 2 */
#if NDIMS == 3
! average fluxes from higher level neighbor
!
      fb(1:nv,:,:) = ((fn(1:nv,1:in:2,1:kn:2) + fn(1:nv,2:in:2,2:kn:2))        &
                    + (fn(1:nv,1:in:2,2:kn:2) + fn(1:nv,2:in:2,1:kn:2)))       &
                    / 4.0d+00
#endif /* NDIMS == 3 */

! Y direction
!
    case(2)

#if NDIMS == 2
! average fluxes from higher level neighbor
!
      fb(1:nv,:,:) = (fn(1:nv,1:in:2,1:kn) + fn(1:nv,2:in:2,1:kn)) / 2.0d+00
#endif /* NDIMS == 2 */
#if NDIMS == 3
! average fluxes from higher level neighbor
!
      fb(1:nv,:,:) = ((fn(1:nv,1:in:2,1:kn:2) + fn(1:nv,2:in:2,2:kn:2))        &
                    + (fn(1:nv,1:in:2,2:kn:2) + fn(1:nv,2:in:2,1:kn:2)))       &
                    / 4.0d+00
#endif /* NDIMS == 3 */

#if NDIMS == 3
! Z direction
!
    case(3)

! average fluxes from higher level neighbor
!
      fb(1:nv,:,:) = ((fn(1:nv,1:in:2,1:jn:2) + fn(1:nv,2:in:2,2:jn:2))        &
                    + (fn(1:nv,1:in:2,2:jn:2) + fn(1:nv,2:in:2,1:jn:2)))       &
                    / 4.0d+00
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_update_flux
!
!===============================================================================
!
!  OTHER BOUNDARY SUBROUTINES
!
!===============================================================================
!
!===============================================================================
!
! subroutine UPDATE_GHOST_CELLS:
! -----------------------------
!
!   Subroutine updates conservative variables in all ghost cells from
!   already updated primitive variables.
!
!
!===============================================================================
!
  subroutine update_ghost_cells()

! include external variables
!
    use blocks        , only : block_data, list_data
    use coordinates   , only : im , jm , km , in , jn , kn
    use coordinates   , only : ib , jb , kb , ie , je , ke
    use coordinates   , only : ibl, jbl, kbl, ieu, jeu, keu
    use equations     , only : nv
    use equations     , only : prim2cons

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer                   :: i, j, k

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time subroutine
!
    call start_timer(imu)
#endif /* PROFILE */

! assign the pointer to the first block on the list
!
    pdata => list_data

! scan all data blocks until the last is reached
!
    do while(associated(pdata))

! update the X and Y boundary ghost cells
!
      do k = 1, km

! update lower layers of the Y boundary
!
        do j = 1, jbl
          call prim2cons(im, pdata%q(1:nv,1:im,j,k), pdata%u(1:nv,1:im,j,k))
        end do ! j = 1, jbl

! update upper layers of the Y boundary
!
        do j = jeu, jm
          call prim2cons(im, pdata%q(1:nv,1:im,j,k), pdata%u(1:nv,1:im,j,k))
        end do ! j = jeu, jm

! update remaining left layers of the X boundary
!
        do i = 1, ibl
          call prim2cons(jn, pdata%q(1:nv,i,jb:je,k), pdata%u(1:nv,i,jb:je,k))
        end do ! i = 1, ibl

! update remaining right layers of the X boundary
!
        do i = ieu, im
          call prim2cons(jn, pdata%q(1:nv,i,jb:je,k), pdata%u(1:nv,i,jb:je,k))
        end do ! i = 1, ibl

      end do ! k = 1, km

#if NDIMS == 3
! update the Z boundary ghost cells
!
      do j = jb, je

! update the remaining front layers of the Z boundary
!
        do k = 1, kbl
          call prim2cons(in, pdata%q(1:nv,ib:ie,j,k), pdata%u(1:nv,ib:ie,j,k))
        end do ! k = 1, kbl

! update the remaining back layers of the Z boundary
!
        do k = keu, km
          call prim2cons(in, pdata%q(1:nv,ib:ie,j,k), pdata%u(1:nv,ib:ie,j,k))
        end do ! k = keu, km

      end do ! j = jb, je
#endif /* NDIMS == 3 */

! assign the pointer to the next block on the list
!
      pdata => pdata%next

    end do ! data blocks

#ifdef PROFILE
! stop accounting time subroutine
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_ghost_cells
#ifdef MPI
!
!===============================================================================
!
! subroutine PREPARE_EXCHANGE_ARRAY:
! ---------------------------------
!
!   Subroutine prepares the arrays for block exchange lists and their counters.
!
!
!===============================================================================
!
  subroutine prepare_exchange_array()

! include external variables
!
    use mpitools      , only : npmax

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: icol, irow
!
!-------------------------------------------------------------------------------
!
! iterate over all elements of the block exchange array
!
    do irow = 0, npmax
      do icol = 0, npmax

! nullify the array element pointer
!
        nullify(barray(irow,icol)%ptr)

! reset the corresponding counter
!
        bcount(irow,icol) = 0

      end do ! icol = 0, npmax
    end do ! irow = 0, npmax

!-------------------------------------------------------------------------------
!
  end subroutine prepare_exchange_array
!
!===============================================================================
!
! subroutine RELEASE_EXCHANGE_ARRAY:
! ---------------------------------
!
!   Subroutine releases objects on the array of block exchange lists.
!
!
!===============================================================================
!
  subroutine release_exchange_array()

! include external variables
!
    use blocks        , only : block_info, pointer_info
    use mpitools      , only : npmax

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer :: icol, irow

! local pointers
!
    type(block_info), pointer :: pinfo
!
!-------------------------------------------------------------------------------
!
! iterate over all elements of the block exchange array
!
    do irow = 0, npmax
      do icol = 0, npmax

! associate pinfo with the first block in the exchange list
!
        pinfo => barray(irow,icol)%ptr

! scan all elements on the exchange list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          barray(irow,icol)%ptr => pinfo%prev

! nullify pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate info block
!
          deallocate(pinfo)

! associate pinfo with the next block
!
          pinfo => barray(irow,icol)%ptr

        end do ! %ptr blocks

      end do ! icol = 0, npmax
    end do ! irow = 0, npmax

!-------------------------------------------------------------------------------
!
  end subroutine release_exchange_array
!
!===============================================================================
!
! subroutine APPEND_EXCHANGE_BLOCK:
! ---------------------------------
!
!   Subroutine appends an info block to the element of array of block
!   exchange lists. The element is determined by the processes of the meta
!   and neighbor blocks.
!
!   Arguments:
!
!     pmeta      - the pointer to meta block;
!     pneigh     - the pointer to the neighbor of pmeta;
!     n, i, j, k - the location of the neighbor;
!
!===============================================================================
!
  subroutine append_exchange_block(pmeta, pneigh, n, i, j, k)

! include external variables
!
    use blocks        , only : block_info, block_meta

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_meta), pointer, intent(inout) :: pmeta, pneigh
    integer                  , intent(in)    :: n, i, j, k

! local variables
!
    integer :: icol, irow

! local pointers
!
    type(block_info), pointer :: pinfo
!
!-------------------------------------------------------------------------------
!
! get the column and row indices
!
    irow = pneigh%process
    icol = pmeta%process

! increase the counter for the number of blocks to exchange
!
    bcount(irow,icol) = bcount(irow,icol) + 1

! allocate a new info object
!
    allocate(pinfo)

! fill out its fields
!
    pinfo%block            => pmeta
    pinfo%neigh            => pneigh
    pinfo%direction        =  n
    pinfo%corner(1)        =  i
    pinfo%corner(2)        =  j
#if NDIMS == 3
    pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */
    pinfo%level_difference =  pmeta%level - pneigh%level

! nullify pointer fields
!
    nullify(pinfo%prev)
    nullify(pinfo%next)

! check if the list is empty
!
    if (associated(barray(irow,icol)%ptr)) then

! if it is, associate the newly created block with it
!
      pinfo%prev => barray(irow,icol)%ptr

    end if ! %ptr associated

! point the list to the newly created block
!
    barray(irow,icol)%ptr => pinfo

!-------------------------------------------------------------------------------
!
  end subroutine append_exchange_block
#endif /* MPI */

!===============================================================================
!
end module
