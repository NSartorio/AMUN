!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2014 Grzegorz Kowal <grzegorz@amuncode.org>
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
  integer                , save :: imi, imv, imf, ims, imc, imr, imp
#endif /* PROFILE */

! parameters corresponding to the boundary type
!
  integer, parameter            :: bnd_periodic   = 0
  integer, parameter            :: bnd_open       = 1
  integer, parameter            :: bnd_outflow    = 2
  integer, parameter            :: bnd_reflective = 3

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
#ifdef MPI
    use mpitools       , only : pdims, pcoords, periodic, npmax
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
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('boundaries:: initialization', imi)
    call set_timer('boundaries:: variables'     , imv)
    call set_timer('boundaries:: fluxes'        , imf)
    call set_timer('boundaries:: specific'      , ims)
    call set_timer('boundaries:: copy'          , imc)
    call set_timer('boundaries:: restrict'      , imr)
    call set_timer('boundaries:: prolong'       , imp)

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
    case default
      bnd_type(3,2) = bnd_periodic
    end select

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
!
!===============================================================================
!
  subroutine boundary_variables()

! import external procedures and variables
!
    use blocks         , only : ndims

! local variables are not implicit by default
!
    implicit none

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
    call boundaries_specific()

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

! update specific boundaries
!
    call boundaries_specific()

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
!   Subroutine updates the numerical fluxes for blocks which have neighbors
!   at higher level. The fluxes of neighbors at higher level are calulated
!   with smaller error, therefore they are restricted down and the flux
!   of lower level meta block is updated.
!
!
!===============================================================================
!
  subroutine boundary_fluxes()

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_meta
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
    use mpitools       , only : nprocs, nproc, npmax
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
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
    integer :: irecv, isend, nblocks, itag, iret
    integer :: l, p

! local arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable  :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
! do not correct fluxes if we do not use adaptive mesh
!
    if (minlev == maxlev) return

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

#ifdef PROFILE
! start accounting time for flux boundary update
!
    call start_timer(imf)
#endif /* PROFILE */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! prepare the array of exchange block lists and its counters
!
    call prepare_exchange_array()
#endif /* MPI */

!! 2. UPDATE THE FLUX BOUNDARIES BETWEEN LOCAL BLOCKS
!!
! associate pmeta with the first block on the meta list
!
    pmeta => list_meta

! scan all meta blocks in the list
!
    do while(associated(pmeta))

! check if the meta block is leaf
!
      if (pmeta%leaf) then

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

! associate pneigh with the current neighbor
!
#if NDIMS == 2
                pneigh => pmeta%edges(i,j,m)%ptr
#endif /* NDIMS == 2 */
#if NDIMS == 3
                pneigh => pmeta%faces(i,j,k,n)%ptr
#endif /* NDIMS == 3 */

! check if the neighbor is associated
!
                if (associated(pneigh)) then

! check if the neighbor is at highed level than the current block
!
                  if (pneigh%level > pmeta%level) then

#ifdef MPI
! check if the current block and its neighbor belong to the same process, if so,
! update fluxes directly
!
                    if (pmeta%process == nproc .and.                           &
                                                 pneigh%process == nproc) then
#endif /* MPI */

! update directional flux from the neighbor
!
                      select case(n)
                      case(1)

! prepare the boundary layer indices depending on the corner position
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

! update the flux edge from the neighbor at higher level
!
                        call block_update_flux(i, j, k, n                      &
                                 , pneigh%data%f(n,1:nv,is,jb:je,kb:ke)        &
                                 ,  pmeta%data%f(n,1:nv,it,jl:ju,kl:ku))

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

! update the flux edge from the neighbor at higher level
!
                        call block_update_flux(i, j, k, n                      &
                                 , pneigh%data%f(n,1:nv,ib:ie,js,kb:ke)        &
                                 ,  pmeta%data%f(n,1:nv,il:iu,jt,kl:ku))

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
                          ks = ke
                          kt = kbl
                        else
                          ks = kbl
                          kt = ke
                        end if

! update the flux edge from the neighbor at higher level
!
                        call block_update_flux(i, j, k, n                      &
                                 , pneigh%data%f(n,1:nv,ib:ie,jb:je,ks)        &
                                 ,  pmeta%data%f(n,1:nv,il:iu,jl:ju,kt))
#endif /* NDIMS == 3 */

                      end select

#ifdef MPI
! blocks belong to different processes, therefore prepare the block exchange
! object
!
                    else

! append the block to the exchange list
!
                      call append_exchange_block(pmeta, pneigh, n, i, j, k)

                    end if ! pmeta and pneigh on local process
#endif /* MPI */

                  end if ! pmeta level < pneigh level

                end if ! pneigh associated

              end do ! i = 1, nsides
            end do ! j = 1, nsides
#if NDIMS == 3
          end do ! k = 1, nsides
#endif /* NDIMS == 3 */
        end do ! n = 1, ndims

      end if ! leaf

! associate pmeta with the next block
!
      pmeta => pmeta%next ! assign pointer to the next meta block in the list

    end do ! meta blocks

#ifdef MPI
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have anything to exchange
!
      if (bcount(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = bcount(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 1

! allocate the buffer for variable exchange
!
        allocate(rbuf(nblocks,nv,ih,kh))

! if isend == nproc we are sending data
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

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
                                 ,          rbuf(l,1:nv,1:jh,1:kh))

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
                                 ,          rbuf(l,1:nv,1:ih,1:kh))

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
                                 ,          rbuf(l,1:nv,1:ih,1:jh))
#endif /* NDIMS == 3 */

            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr blocks

! send the data buffer to another process
!
          call send_real_array(size(rbuf(:,:,:,:)), irecv, itag                &
                                                        , rbuf(:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:)), isend, itag             &
                                                        , rbuf(:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

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

        end if ! irecv = nproc

! deallocate data buffer
!
        deallocate(rbuf)

      end if ! if bcount > 0

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for flux boundary update
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
!   neighbors, then updates its boundaries for selected type.
!
!
!===============================================================================
!
  subroutine boundaries_specific()

! import external procedures and variables
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : ndims, nsides
    use coordinates    , only : im, jm, km
    use equations      , only : nv
#ifdef MPI
    use mpitools       , only : nproc
#endif /* MPI */
    use mpitools       , only : periodic

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh

! local variables
!
    integer                   :: i, j, k, n, m
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for specific boundary update
!
    call start_timer(ims)
#endif /* PROFILE */

! associate pmeta with the first block on the meta list
!
    pmeta => list_meta

! scan all blocks on meta block list
!
    do while(associated(pmeta))

! check if the current meta block is a leaf
!
      if (pmeta%leaf) then

! process only if this block is marked for update
!
        if (pmeta%update) then

#ifdef MPI
! check if the current block belongs to the local process
!
          if (pmeta%process == nproc) then
#endif /* MPI */

#if NDIMS == 2
! iterate over all directions
!
            do n = 1, ndims

! process boundaries only if they are not periodic in a given direction
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

! process boundaries only if they are not periodic in a given direction
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
                                          , pmeta%data%q(1:nv,1:im,1:jm,1:km))

                    end do ! i = 1, sides
                  end do ! j = 1, sides
                end do ! k = 1, sides

              end if ! not periodic

            end do ! n = 1, ndims
#endif /* NDIMS == 3 */

#ifdef MPI
          end if ! block belong to the local process
#endif /* MPI */

        end if ! pmeta is marked for update

      end if ! leaf

! associate pmeta with the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef PROFILE
! stop accounting time for specific boundary update
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
!   Subroutine scans over all leaf blocks in order to find face neighbors which
!   are the same level, and perform the update of the face boundaries between
!   them.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for copy boundary update
!
    call start_timer(imc)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
    kh = kn / 2

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! prepare the array of exchange block lists and its counters
!
    call prepare_exchange_array()
#endif /* MPI */

!! 2. UPDATE VARIABLE FACE BOUNDARIES BETWEEN BLOCKS BELONGING TO THE SAME
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pneigh%process == nproc) then
#endif /* MPI */

! prepare region indices for the face boundary update
!
                        select case(idir)
                        case(1)
                          if (i == 1) then
                            il = 1
                            iu = ibl
                          else
                            il = ieu
                            iu = im
                          end if
                          if (j == 1) then
                            jl = jb
                            ju = jb + jh - 1
                          else
                            jl = je - jh + 1
                            ju = je
                          end if
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
                        case(2)
                          if (i == 1) then
                            il = ib
                            iu = ib + ih - 1
                          else
                            il = ie - ih + 1
                            iu = ie
                          end if
                          if (j == 1) then
                            jl = 1
                            ju = jbl
                          else
                            jl = jeu
                            ju = jm
                          end if
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
                        case(3)
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
                            kl = 1
                            ku = kbl
                          else
                            kl = keu
                            ku = km
                          end if
                        end select

! extract the corresponding face region from the neighbor and insert it in
! the current data block
!
                        call block_face_copy(idir, i, j, k                     &
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

      end if ! leaf

! associate pmeta with the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have anything to exchange
!
      if (bcount(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = bcount(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 2

! allocate data buffer for variables to exchange
!
        select case(idir)
        case(1)
          allocate(rbuf(nblocks,nv,ng,jh,kh))
        case(2)
          allocate(rbuf(nblocks,nv,ih,ng,kh))
        case(3)
          allocate(rbuf(nblocks,nv,ih,jh,ng))
        end select

! if isend == nproc we are sending data
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! associate pneigh with pinfo%neigh
!
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
              call block_face_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:kh))
            case(2)
              call block_face_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:kh))
            case(3)
              call block_face_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:jh,1:ng))
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! associate pmeta with pinfo%block
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
              if (i == 1) then
                il = 1
                iu = ibl
              else
                il = ieu
                iu = im
              end if
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:kh)
            case(2)
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
              if (j == 1) then
                jl = 1
                ju = jbl
              else
                jl = jeu
                ju = jm
              end if
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:kh)
            case(3)
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
                kl = 1
                ku = kbl
              else
                kl = keu
                ku = km
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                 rbuf(l,1:nv,1:ih,1:jh,1:ng)
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

      end if ! if bcount > 0

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for copy boundary update
!
    call stop_timer(imc)
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
!   Subroutine scans over all leaf blocks in order to find face neighbors which
!   are on different levels, and perform the update of face boundaries of
!   lower blocks by restricting them from higher level neighbors.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for restrict boundary update
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
    kh = kn / 2

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! prepare the array of exchange block lists and its counters
!
    call prepare_exchange_array()
#endif /* MPI */

!! 2. UPDATE VARIABLE FACE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for face boundary update
!
                        select case(idir)
                        case(1)
                          if (i == 1) then
                            il = 1
                            iu = ibl
                          else
                            il = ieu
                            iu = im
                          end if
                          if (j == 1) then
                            jl = jb
                            ju = jb + jh - 1
                          else
                            jl = je - jh + 1
                            ju = je
                          end if
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
                        case(2)
                          if (i == 1) then
                            il = ib
                            iu = ib + ih - 1
                          else
                            il = ie - ih + 1
                            iu = ie
                          end if
                          if (j == 1) then
                            jl = 1
                            ju = jbl
                          else
                            jl = jeu
                            ju = jm
                          end if
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
                        case(3)
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
                            kl = 1
                            ku = kbl
                          else
                            kl = keu
                            ku = km
                          end if
                        end select

! extract the corresponding face region from the neighbor and insert it in
! the current data block
!
                        call block_face_restrict(idir, i, j, k                 &
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

      end if ! leaf

! associate pmeta with the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (bcount(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = bcount(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 3

! allocate data buffer for variables to exchange
!
        select case(idir)
        case(1)
          allocate(rbuf(nblocks,nv,ng,jh,kh))
        case(2)
          allocate(rbuf(nblocks,nv,ih,ng,kh))
        case(3)
          allocate(rbuf(nblocks,nv,ih,jh,ng))
        end select

! if isend == nproc we are sending data
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! associate pneigh with pinfo%neigh
!
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
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:kh))
            case(2)
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:kh))
            case(3)
              call block_face_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:jh,1:ng))
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! associate pmeta with pinfo%block
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
              if (i == 1) then
                il = 1
                iu = ibl
              else
                il = ieu
                iu = im
              end if
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:kh)
            case(2)
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
              if (j == 1) then
                jl = 1
                ju = jbl
              else
                jl = jeu
                ju = jm
              end if
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:kh)
            case(3)
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
                kl = 1
                ku = kbl
              else
                kl = keu
                ku = km
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:jh,1:ng)
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

      end if ! if bcount > 0

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for restrict boundary update
!
    call stop_timer(imr)
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
!   Subroutine scans over all leaf blocks in order to find face neighbors which
!   are on different levels, and perform the update of face boundaries of
!   higher blocks by prolongating them from lower level neighbors.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local arrays
!
    real(kind=8), dimension(:,:,:,:,:), allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for prolong boundary update
!
    call start_timer(imp)
#endif /* PROFILE */

! calculate the sizes
!
    ih = in + ng
    jh = jn + ng
    kh = kn + ng

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! prepare the array of exchange block lists and its counters
!
    call prepare_exchange_array()
#endif /* MPI */

!! 2. UPDATE VARIABLE FACE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! extract the corresponding face region from the neighbor and insert it in
! the current data block
!
                        select case(idir)
                        case(1)
                          jc = pmeta%pos(2)
                          kc = pmeta%pos(3)
                          if (i == 1) then
                            il = 1
                            iu = ibl
                          else
                            il = ieu
                            iu = im
                          end if
                          if (jc == 0) then
                            jl = jb
                            ju = jm
                          else
                            jl = 1
                            ju = je
                          end if
                          if (kc == 0) then
                            kl = kb
                            ku = km
                          else
                            kl = 1
                            ku = ke
                          end if

                        case(2)
                          ic = pmeta%pos(1)
                          kc = pmeta%pos(3)
                          if (ic == 0) then
                            il = ib
                            iu = im
                          else
                            il = 1
                            iu = ie
                          end if
                          if (j == 1) then
                            jl = 1
                            ju = jbl
                          else
                            jl = jeu
                            ju = jm
                          end if
                          if (kc == 0) then
                            kl = kb
                            ku = km
                          else
                            kl = 1
                            ku = ke
                          end if

                        case(3)
                          ic = pmeta%pos(1)
                          jc = pmeta%pos(2)
                          if (ic == 0) then
                            il = ib
                            iu = im
                          else
                            il = 1
                            iu = ie
                          end if
                          if (jc == 0) then
                            jl = jb
                            ju = jm
                          else
                            jl = 1
                            ju = je
                          end if
                          if (k == 1) then
                            kl = 1
                            ku = kbl
                          else
                            kl = keu
                            ku = km
                          end if
                        end select

! extract the corresponding face region from the neighbor and insert it in
! the current data block
!
                        call block_face_prolong(idir, ic, jc, kc               &
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

      end if ! leaf

! associate pmeta with the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (bcount(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = bcount(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 4

! allocate data buffer for variables to exchange
!
        select case(idir)
        case(1)
          allocate(rbuf(nblocks,nv,ng,jh,kh))
        case(2)
          allocate(rbuf(nblocks,nv,ih,ng,kh))
        case(3)
          allocate(rbuf(nblocks,nv,ih,jh,ng))
        end select

! if isend == nproc we are sending data
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! scan over all blocks on the block exchange list
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! prepare pointer for updated meta block and its neighbor
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
              j = pmeta%pos(2)
              k = pmeta%pos(3)
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:kh))
            case(2)
              i = pmeta%pos(1)
              k = pmeta%pos(3)
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:kh))
            case(3)
              i = pmeta%pos(1)
              j = pmeta%pos(2)
              call block_face_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:jh,1:ng))
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate pinfo with the first block in the exchange list
!
          pinfo => barray(isend,irecv)%ptr

! iterate over all received blocks and update boundaries of the corresponding
! data blocks
!
          do while(associated(pinfo))

! increase the block counter
!
            l = l + 1

! prepare the pointer to updated block
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
              if (i == 1) then
                il = 1
                iu = ibl
              else
                il = ieu
                iu = im
              end if
              if (pmeta%pos(2) == 0) then
                jl = jb
                ju = jm
              else
                jl =  1
                ju = je
              end if
              if (pmeta%pos(3) == 0) then
                kl = kb
                ku = km
              else
                kl =  1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:jh,1:kh)
            case(2)
              if (j == 1) then
                jl = 1
                ju = jbl
              else
                jl = jeu
                ju = jm
              end if
              if (pmeta%pos(1) == 0) then
                il = ib
                iu = im
              else
                il =  1
                iu = ie
              end if
              if (pmeta%pos(3) == 0) then
                kl = kb
                ku = km
              else
                kl =  1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:kh)
            case(3)
              if (k == 1) then
                kl = 1
                ku = kbl
              else
                kl = keu
                ku = km
              end if
              if (pmeta%pos(1) == 0) then
                il = ib
                iu = im
              else
                il =  1
                iu = ie
              end if
              if (pmeta%pos(2) == 0) then
                jl = jb
                ju = jm
              else
                jl =  1
                ju = je
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:jh,1:ng)
            end select

! associate pinfo with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

      end if ! if bcount > 0

    end do ! p = 1, npairs

! release the memory used by the array of exchange block lists
!
    call release_exchange_array()
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for prolong boundary update
!
    call stop_timer(imp)
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
!   Subroutine scans over all leaf blocks in order to find edge neighbors which
!   are the same level, and perform the update of the edge boundaries between
!   them.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for copy boundary update
!
    call start_timer(imc)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
#if NDIMS == 3
    kh = kn / 2
#endif /* NDIMS == 3 */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE EDGE BOUNDARIES BETWEEN BLOCKS BELONGING TO THE SAME
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for edge boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor and insert it in
! the current data block
!
                        select case(idir)
                        case(1)
                          if (i == 1) then
                            il = ib
                            iu = ib + ih - 1
                          else
                            il = ie - ih + 1
                            iu = ie
                          end if
                        case(2)
                          if (j == 1) then
                            jl = jb
                            ju = jb + jh - 1
                          else
                            jl = je - jh + 1
                            ju = je
                          end if
#if NDIMS == 3
                        case(3)
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
#endif /* NDIMS == 3 */
                        end select
#if NDIMS == 2
                        call block_edge_copy(idir, i, j, k                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_edge_copy(idir, i, j, k                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! pneigh on the current process

                    else ! block and neighbor belong to different processes

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%direction        =  idir
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! associate the pointer to the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 5

! allocate data buffer for variables to exchange
!
        select case(idir)
#if NDIMS == 2
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,km))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,ng))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,ng))
        case(3)
          allocate(rbuf(nblocks,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
        end select

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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
              call block_edge_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:ng))
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              call block_edge_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:ng))
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              call block_edge_copy(idir, i, j, k                               &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:kh))
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
            select case(idir)
            case(1)
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */
            case(2)
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
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
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for copy boundary update
!
    call stop_timer(imc)
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
!   Subroutine scans over all leaf blocks in order to find edge neighbors which
!   are on different levels, and perform the update of edge boundaries of
!   lower blocks by restricting them from higher level neighbors.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for restrict boundary update
!
    call start_timer(imr)
#endif /* PROFILE */

! calculate half sizes
!
    ih = in / 2
    jh = jn / 2
#if NDIMS == 3
    kh = kn / 2
#endif /* NDIMS == 3 */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE EDGE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for edge boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor and insert it in
! the current data block
!
                        select case(idir)
                        case(1)
                          if (i == 1) then
                            il = ib
                            iu = ib + ih - 1
                          else
                            il = ie - ih + 1
                            iu = ie
                          end if
                        case(2)
                          if (j == 1) then
                            jl = jb
                            ju = jb + jh - 1
                          else
                            jl = je - jh + 1
                            ju = je
                          end if
#if NDIMS == 3
                        case(3)
                          if (k == 1) then
                            kl = kb
                            ku = kb + kh - 1
                          else
                            kl = ke - kh + 1
                            ku = ke
                          end if
#endif /* NDIMS == 3 */
                        end select
#if NDIMS == 2
                        call block_edge_restrict(idir, i, j, k                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_edge_restrict(idir, i, j, k                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! pneigh on the current process

                    else ! block and neighbor belong to different processes

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%direction        =  idir
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! associate the pointer to the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 6

! allocate data buffer for variables to exchange
!
        select case(idir)
#if NDIMS == 2
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,km))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,ng))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,ng))
        case(3)
          allocate(rbuf(nblocks,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
        end select

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:ng))
#endif /* NDIMS == 3 */
            case(2)
#if NDIMS == 2
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:ng))
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              call block_edge_restrict(idir, i, j, k                           &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:kh))
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
            select case(idir)
            case(1)
              if (i == 1) then
                il = ib
                iu = ib + ih - 1
              else
                il = ie - ih + 1
                iu = ie
              end if
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */
            case(2)
              if (j == 1) then
                jl = jb
                ju = jb + jh - 1
              else
                jl = je - jh + 1
                ju = je
              end if
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
              if (k == 1) then
                kl = kb
                ku = kb + kh - 1
              else
                kl = ke - kh + 1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for restrict boundary update
!
    call stop_timer(imr)
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
!   Subroutine scans over all leaf blocks in order to find edge neighbors which
!   are on different levels, and perform the update of edge boundaries of
!   higher blocks by prolongating them from lower level neighbors.
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
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
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
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for prolong boundary update
!
    call start_timer(imp)
#endif /* PROFILE */

! calculate the sizes
!
    ih = in + ng
    jh = jn + ng
#if NDIMS == 3
    kh = kn + ng
#endif /* NDIMS == 3 */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE EDGE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for edge boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! extract the corresponding edge region from the neighbor and insert it in
! the current data block
!
                        select case(idir)
                        case(1)
                          ic = pmeta%pos(1)
                          if (ic == 0) then
                            il = ib
                            iu = im
                          else
                            il = 1
                            iu = ie
                          end if
                        case(2)
                          jc = pmeta%pos(2)
                          if (jc == 0) then
                            jl = jb
                            ju = jm
                          else
                            jl = 1
                            ju = je
                          end if
#if NDIMS == 3
                        case(3)
                          kc = pmeta%pos(3)
                          if (kc == 0) then
                            kl = kb
                            ku = km
                          else
                            kl = 1
                            ku = ke
                          end if
#endif /* NDIMS == 3 */
                        end select
#if NDIMS == 2
                        call block_edge_prolong(idir, ic, jc, kc               &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_edge_prolong(idir, ic, jc, kc               &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! pneigh on the current process

                    else ! block and neighbor belong to different processes

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%direction        =  idir
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! associate the pointer to the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 7

! allocate data buffer for variables to exchange
!
        select case(idir)
#if NDIMS == 2
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,km))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        case(1)
          allocate(rbuf(nblocks,nv,ih,ng,ng))
        case(2)
          allocate(rbuf(nblocks,nv,ng,jh,ng))
        case(3)
          allocate(rbuf(nblocks,nv,ng,ng,kh))
#endif /* NDIMS == 3 */
        end select

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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
              i = pmeta%pos(1)
#if NDIMS == 2
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ih,1:ng,1:ng))
#endif /* NDIMS == 3 */
            case(2)
              j = pmeta%pos(2)
#if NDIMS == 2
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:jh,1:ng))
#endif /* NDIMS == 3 */
#if NDIMS == 3
            case(3)
              k = pmeta%pos(3)
              call block_edge_prolong(idir, i, j, k                            &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:kh))
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
#endif /* NDIMS == 3 */

! update the corresponding corner region of the current block
!
            select case(idir)
            case(1)
              if (pmeta%pos(1) == 0) then
                il = ib
                iu = im
              else
                il =  1
                iu = ie
              end if
#if NDIMS == 2
              pmeta%data%q(1:nv,il:iu,jl:ju, 1:km) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ih,1:ng,1:ng)
#endif /* NDIMS == 3 */
            case(2)
              if (pmeta%pos(2) == 0) then
                jl = jb
                ju = jm
              else
                jl =  1
                ju = je
              end if
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
              if (pmeta%pos(3) == 0) then
                kl = kb
                ku = km
              else
                kl =  1
                ku = ke
              end if
              pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku) =                           &
                                                   rbuf(l,1:nv,1:ng,1:ng,1:kh)
#endif /* NDIMS == 3 */
            end select

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for prolong boundary update
!
    call stop_timer(imp)
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
!   Subroutine scans over all leaf blocks in order to find corner neighbors at
!   the same level, and perform the update of the corner boundaries between
!   them.
!
!
!===============================================================================
!
  subroutine boundaries_corner_copy()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for copy boundary update
!
    call start_timer(imc)
#endif /* PROFILE */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE CORNER BOUNDARIES BETWEEN BLOCKS BELONGING TO THE SAME
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for corner boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! extract the corresponding corner region from the neighbor and insert it in
! the current data block
!
#if NDIMS == 2
                        call block_corner_copy(i, j, k                         &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_corner_copy(i, j, k                         &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! pneigh on the current process

                    else ! block and neighbor belong to different processes

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! associate the pointer to the next meta block
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 8

! allocate data buffer for variables to exchange
!
#if NDIMS == 2
        allocate(rbuf(nblocks,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        allocate(rbuf(nblocks,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! extract the corresponding corner region from the neighbor and insert it
! to the buffer
!
#if NDIMS == 2
            call block_corner_copy(i, j, k                                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
            call block_corner_copy(i, j, k                                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:ng))
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
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

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for copy boundary update
!
    call stop_timer(imc)
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
!   Subroutine scans over all leaf blocks in order to find corner neighbors
!   which are on different levels, and perform the update of corner boundaries
!   of lower blocks by restricting them from higher level neighbors.
!
!
!===============================================================================
!
  subroutine boundaries_corner_restrict()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for restrict boundary update
!
    call start_timer(imr)
#endif /* PROFILE */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE CORNER BOUNDARIES BETWEEN BLOCKS BELONGING TO THE SAME
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for corner boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! restrict and extract the corresponding corner region from the neighbor and
! insert it in the current data block
!
#if NDIMS == 2
                        call block_corner_restrict(i, j, k                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_corner_restrict(i, j, k                     &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! block on the current processor

                    else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! assign the pointer to the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 9

! allocate data buffer for variables to exchange
!
#if NDIMS == 2
        allocate(rbuf(nblocks,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        allocate(rbuf(nblocks,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
            call block_corner_restrict(i, j, k                                 &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:ng))
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
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

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for restrict boundary update
!
    call stop_timer(imr)
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
!   Subroutine scans over all leaf blocks in order to find corner neighbors at
!   different levels, and update the corner boundaries of blocks at higher
!   levels by prolongating variables from lower level blocks.
!
!
!===============================================================================
!
  subroutine boundaries_corner_prolong()

! import external procedures and variables
!
    use blocks         , only : nsides
    use blocks         , only : block_meta, block_data
    use blocks         , only : list_meta
    use blocks         , only : block_info, pointer_info
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : ibl, jbl, kbl
    use coordinates    , only : ieu, jeu, keu
    use equations      , only : nv
    use mpitools       , only : nproc, nprocs, npmax
#ifdef MPI
    use mpitools       , only : npairs, pairs
    use mpitools       , only : send_real_array, receive_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local pointers
!
    type(block_meta), pointer :: pmeta, pneigh
#ifdef MPI
    type(block_info), pointer :: pinfo
#endif /* MPI */

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: iret
#ifdef MPI
    integer :: isend, irecv, nblocks, itag
    integer :: l, p

! local pointer arrays
!
    type(pointer_info), dimension(0:npmax,0:npmax) :: block_array

! local arrays
!
    integer     , dimension(0:npmax,0:npmax)              :: block_counter
    real(kind=8), dimension(:,:,:,:,:)      , allocatable :: rbuf
#endif /* MPI */
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for prolong boundary update
!
    call start_timer(imp)
#endif /* PROFILE */

#ifdef MPI
!! 1. PREPARE THE BLOCK EXCHANGE ARRAYS FOR MPI
!!
! reset the exchange block counters
!
    block_counter(:,:) = 0

! nullify the info pointers
!
    do irecv = 0, npmax
      do isend = 0, npmax
        nullify(block_array(isend,irecv)%ptr)
      end do
    end do
#endif /* MPI */

!! 2. UPDATE VARIABLE CORNER BOUNDARIES BETWEEN BLOCKS BELONGING TO THE SAME
!!    PROCESS AND PREPARE THE EXCHANGE BLOCK LIST OF BLOCKS WHICH BELONG TO
!!    DIFFERENT PROCESSES
!!
! associate pmeta with the first block on the meta block list
!
    pmeta => list_meta

! scan all meta blocks
!
    do while(associated(pmeta))

! check if the block is leaf
!
      if (pmeta%leaf) then

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
                  if (pmeta%update .and. pneigh%update) then

#ifdef MPI
! check if the block and its neighbor belong to the same process
!
                    if (pmeta%process == pneigh%process) then

! check if the neighbor belongs to the current process
!
                      if (pmeta%process == nproc) then
#endif /* MPI */

! prepare the region indices for corner boundary update
!
                        if (i == 1) then
                          il = 1
                          iu = ibl
                        else
                          il = ieu
                          iu = im
                        end if
                        if (j == 1) then
                          jl = 1
                          ju = jbl
                        else
                          jl = jeu
                          ju = jm
                        end if
#if NDIMS == 3
                        if (k == 1) then
                          kl = 1
                          ku = kbl
                        else
                          kl = keu
                          ku = km
                        end if
#endif /* NDIMS == 3 */

! restrict and extract the corresponding corner region from the neighbor and
! insert it in the current data block
!
#if NDIMS == 2
                        call block_corner_prolong(i, j, k                      &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju, 1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
                        call block_corner_prolong(i, j, k                      &
                                   , pneigh%data%q(1:nv, 1:im, 1:jm, 1:km)     &
                                   ,  pmeta%data%q(1:nv,il:iu,jl:ju,kl:ku))
#endif /* NDIMS == 3 */

#ifdef MPI
                      end if ! block on the current processor

                    else ! block and neighbor on different processors

! increase the counter for number of blocks to exchange
!
                      block_counter(pneigh%process,pmeta%process) =            &
                               block_counter(pneigh%process,pmeta%process) + 1

! allocate a new info object
!
                      allocate(pinfo)

! fill out only fields which are used
!
                      pinfo%block            => pmeta
                      pinfo%neigh            => pneigh
                      pinfo%corner(1)        =  i
                      pinfo%corner(2)        =  j
#if NDIMS == 3
                      pinfo%corner(3)        =  k
#endif /* NDIMS == 3 */

! nullify pointer fields of the object
!
                      nullify(pinfo%prev)
                      nullify(pinfo%next)

! if the list is not empty append the newly created block to it
!
                      if (associated(block_array(pneigh%process                &
                                                        ,pmeta%process)%ptr))  &
                        pinfo%prev => block_array(pneigh%process               &
                                                        ,pmeta%process)%ptr

! point the list to the newly created block
!
                      block_array(pneigh%process,pmeta%process)%ptr => pinfo

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

      end if ! leaf

! assign the pointer to the next block on the list
!
      pmeta => pmeta%next

    end do ! meta blocks

#ifdef MPI
!! 3. UPDATE VARIABLE BOUNDARIES BETWEEN BLOCKS BELONGING TO DIFFERENT PROCESSES
!!
! iterate over all process pairs
!
    do p = 1, npairs

! get sending and receiving process identifiers
!
      isend = pairs(p,1)
      irecv = pairs(p,2)

! process only pairs which have something to exchange
!
      if (block_counter(isend,irecv) > 0) then

! obtain the number of blocks to exchange
!
        nblocks = block_counter(isend,irecv)

! prepare the tag for communication
!
        itag = 16 * (irecv * nprocs + isend) + 10

! allocate data buffer for variables to exchange
!
#if NDIMS == 2
        allocate(rbuf(nblocks,nv,ng,ng,km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        allocate(rbuf(nblocks,nv,ng,ng,ng))
#endif /* NDIMS == 3 */

! if isend == nproc we are sending data from the neighbor block
!
        if (isend == nproc) then

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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
            call block_corner_prolong(i, j, k                                  &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:km))
#endif /* NDIMS == 2 */
#if NDIMS == 3
            call block_corner_prolong(i, j, k                                  &
                                   , pneigh%data%q(1:nv,1:im,1:jm,1:km)        &
                                   ,        rbuf(l,1:nv,1:ng,1:ng,1:ng))
#endif /* NDIMS == 3 */

! associate the pointer with the next block
!
            pinfo => pinfo%prev

          end do ! %ptr block list

! send the data buffer to another process
!
          call send_real_array(size(rbuf), irecv, itag, rbuf(:,:,:,:,:), iret)

        end if ! isend = nproc

! if irecv == nproc we are receiving data from the neighbor block
!
        if (irecv == nproc) then

! receive the data buffer
!
          call receive_real_array(size(rbuf(:,:,:,:,:)), isend, itag           &
                                                      , rbuf(:,:,:,:,:), iret)

! reset the block counter
!
          l = 0

! associate the pointer with the first block in the exchange list
!
          pinfo => block_array(isend,irecv)%ptr

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

! calculate the insertion indices
!
            if (i == 1) then
              il = 1
              iu = ibl
            else
              il = ieu
              iu = im
            end if
            if (j == 1) then
              jl = 1
              ju = jbl
            else
              jl = jeu
              ju = jm
            end if
#if NDIMS == 3
            if (k == 1) then
              kl = 1
              ku = kbl
            else
              kl = keu
              ku = km
            end if
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

        end if ! irecv = nproc

! deallocate data buffer
!
        if (allocated(rbuf)) deallocate(rbuf)

! associate the pointer with the first block in the exchange list
!
        pinfo => block_array(isend,irecv)%ptr

! scan over all blocks on the exchange block list
!
        do while(associated(pinfo))

! associate the exchange list pointer
!
          block_array(isend,irecv)%ptr => pinfo%prev

! nullify the pointer fields
!
          nullify(pinfo%prev)
          nullify(pinfo%next)
          nullify(pinfo%block)
          nullify(pinfo%neigh)

! deallocate the object
!
          deallocate(pinfo)

! associate the pointer with the next block
!
          pinfo => block_array(isend,irecv)%ptr

        end do ! %ptr block list

      end if ! if block_count > 0

    end do ! p = 1, npairs
#endif /* MPI */

#ifdef PROFILE
! stop accounting time for prolong boundary update
!
    call stop_timer(imp)
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
!     qn         - the variable array;
!
!===============================================================================
!
  subroutine block_boundary_specific(ic, jc, kc, nc, qn)

! import external procedures and variables
!
    use coordinates    , only : im , jm , km , ng
    use coordinates    , only : ib , jb , kb , ie , je , ke
    use coordinates    , only : ibl, jbl, kbl, ieu, jeu, keu
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ibx, iby, ibz, ibp
    use error          , only : print_error, print_warning

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: ic, jc, kc
    integer                                     , intent(in)    :: nc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn

! local variables
!
    integer :: i , j , k
    integer :: il, jl, kl
    integer :: iu, ju, ku
    integer :: is, js, ks
    integer :: it, jt, kt
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
          end do
        else
          do i = 1, ng
            it = ie  + i
            is = ieu - i

            qn(1:nv,it,jl:ju,kl:ku) =   qn(1:nv,is,jl:ju,kl:ku)
            qn(ivx ,it,jl:ju,kl:ku) = - qn(ivx ,is,jl:ju,kl:ku)
          end do
        end if

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
          end do
        else
          do j = 1, ng
            jt = je  + j
            js = jeu - j

            qn(1:nv,il:iu,jt,kl:ku) =   qn(1:nv,il:iu,js,kl:ku)
            qn(ivy ,il:iu,jt,kl:ku) = - qn(ivy ,il:iu,js,kl:ku)
          end do
        end if

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
          end do
        else
          do k = 1, ng
            kt = ke  + k
            ks = keu - k

            qn(1:nv,il:iu,jl:ju,kt) =   qn(1:nv,il:iu,jl:ju,ks)
            qn(ivz ,il:iu,jl:ju,kt) = - qn(ivz ,il:iu,jl:ju,ks)
          end do
        end if

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
! subroutine BLOCK_FACE_COPY:
! --------------------------
!
!   Subroutine returns the face boundary region copied from the provided
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
  subroutine block_face_copy(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibu, jbu, kbu
    use coordinates    , only : iel, jel, kel
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local indices
!
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half sizes
!
      jh = jn / 2
      kh = kn / 2

! prepare indices for the face region
!
      if (ic == 1) then
        il = iel
        iu = ie
      else
        il = ib
        iu = ibu
      end if
      if (jc == 1) then
        jl = jb
        ju = jb + jh - 1
      else
        jl = je - jh + 1
        ju = je
      end if
      if (kc == 1) then
        kl = kb
        ku = kb + kh - 1
      else
        kl = ke - kh + 1
        ku = ke
      end if

! copy the face region to the output array
!
      qb(1:nv,1:ng,1:jh,1:kh) = qn(1:nv,il:iu,jl:ju,kl:ku)

    case(2)

! calculate half sizes
!
      ih = in / 2
      kh = kn / 2

! prepare indices for the face region
!
      if (ic == 1) then
        il = ib
        iu = ib + ih - 1
      else
        il = ie - ih + 1
        iu = ie
      end if
      if (jc == 1) then
        jl = jel
        ju = je
      else
        jl = jb
        ju = jbu
      end if
      if (kc == 1) then
        kl = kb
        ku = kb + kh - 1
      else
        kl = ke - kh + 1
        ku = ke
      end if

! copy the face region to the output array
!
      qb(1:nv,1:ih,1:ng,1:kh) = qn(1:nv,il:iu,jl:ju,kl:ku)

    case(3)

! calculate half sizes
!
      ih = in / 2
      jh = jn / 2

! prepare indices for the face region
!
      if (ic == 1) then
        il = ib
        iu = ib + ih - 1
      else
        il = ie - ih + 1
        iu = ie
      end if
      if (jc == 1) then
        jl = jb
        ju = jb + jh - 1
      else
        jl = je - jh + 1
        ju = je
      end if
      if (kc == 1) then
        kl = kel
        ku = ke
      else
        kl = kb
        ku = kbu
      end if

! copy the face region to the output array
!
      qb(1:nv,1:ih,1:jh,1:ng) = qn(1:nv,il:iu,jl:ju,kl:ku)

    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_face_copy
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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
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
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half sizes
!
      jh = jn / 2
      kh = kn / 2

! prepare indices for the face region
!
      if (ic == 1) then
        il = ie - nd + 1
        ip = il + 1
        iu = ie
      else
        il = ib
        ip = il + 1
        iu = ib + nd - 1
      end if
      jl = jb
      jp = jl + 1
      ju = je
      kl = kb
      kp = kl + 1
      ku = ke

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

! prepare indices for the face region
!
      il = ib
      ip = il + 1
      iu = ie
      if (jc == 1) then
        jl = je - nd + 1
        jp = jl + 1
        ju = je
      else
        jl = jb
        jp = jl + 1
        ju = jb + nd - 1
      end if
      kl = kb
      kp = kl + 1
      ku = ke

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

! prepare indices for the face region
!
      il = ib
      ip = il + 1
      iu = ie
      jl = jb
      jp = jl + 1
      ju = je
      if (kc == 1) then
        kl = ke - nd + 1
        kp = kl + 1
        ku = ke
      else
        kl = kb
        kp = kl + 1
        ku = kb + nd - 1
      end if

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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use equations      , only : nv
    use interpolations , only : limiter

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
    real(kind=8) :: dqx, dqy, dqz
    real(kind=8) :: dq1, dq2, dq3, dq4
!
!-------------------------------------------------------------------------------
!
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half sizes
!
      jh = jn / 2
      kh = kn / 2

! prepare indices for the face region
!
      if (ic == 1) then
        il = ie - nh + 1
        iu = ie
      else
        il = ib
        iu = ib + nh - 1
      end if
      if (jc == 0) then
        jl = jb
        ju = jb + jh + nh - 1
      else
        jl = je - jh - nh + 1
        ju = je
      end if
      if (kc == 0) then
        kl = kb
        ku = kb + kh + nh - 1
      else
        kl = ke - kh - nh + 1
        ku = ke
      end if

    case(2)

! calculate half sizes
!
      ih = in / 2
      kh = kn / 2

! prepare indices for the face region
!
      if (ic == 0) then
        il = ib
        iu = ib + ih + nh - 1
      else
        il = ie - ih - nh + 1
        iu = ie
      end if
      if (jc == 1) then
        jl = je - nh + 1
        ju = je
      else
        jl = jb
        ju = jb + nh - 1
      end if
      if (kc == 0) then
        kl = kb
        ku = kb + kh + nh - 1
      else
        kl = ke - kh - nh + 1
        ku = ke
      end if

    case(3)

! calculate half sizes
!
      ih = in / 2
      jh = jn / 2

! prepare indices for the face region
!
      if (ic == 0) then
        il = ib
        iu = ib + ih + nh - 1
      else
        il = ie - ih - nh + 1
        iu = ie
      end if
      if (jc == 0) then
        jl = jb
        ju = jb + jh + nh - 1
      else
        jl = je - jh - nh + 1
        ju = je
      end if
      if (kc == 1) then
        kl = ke - nh + 1
        ku = ke
      else
        kl = kb
        ku = kb + nh - 1
      end if

    end select

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
            dql = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dqx = limiter(0.25d+00, dql, dqr)

            dql = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dqy = limiter(0.25d+00, dql, dqr)

            dql = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dqz = limiter(0.25d+00, dql, dqr)

! calculate the derivative terms
!
            dq1 = dqx + dqy + dqz
            dq2 = dqx - dqy - dqz
            dq3 = dqx - dqy + dqz
            dq4 = dqx + dqy - dqz

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
! subroutine BLOCK_EDGE_COPY:
! --------------------------
!
!   Subroutine returns the edge boundary region by copying the corresponding
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
  subroutine block_edge_copy(nc, ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng
    use coordinates    , only : in , jn , kn
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibu, jbu, kbu
    use coordinates    , only : iel, jel, kel
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)  :: nc, ic, jc, kc
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(in)  :: qn
    real(kind=8), dimension( :  , :  , :  , :  ), intent(out) :: qb

! local indices
!
    integer :: ih, jh, kh
    integer :: il, jl, kl
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half size
!
      ih = in / 2

! prepare indices for the edge region
!
      if (ic == 1) then
        il = ib
        iu = ib + ih - 1
      else
        il = ie - ih + 1
        iu = ie
      end if
      if (jc == 1) then
        jl = jel
        ju = je
      else
        jl = jb
        ju = jbu
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = kel
        ku = ke
      else
        kl = kb
        ku = kbu
      end if
#endif /* NDIMS == 3 */

! copy the edge region to the output array
!
#if NDIMS == 2
      qb(1:nv,1:ih,1:ng,1:km) = qn(1:nv,il:iu,jl:ju, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      qb(1:nv,1:ih,1:ng,1:ng) = qn(1:nv,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

    case(2)

! calculate half size
!
      jh = jn / 2

! prepare indices for the edge region
!
      if (ic == 1) then
        il = iel
        iu = ie
      else
        il = ib
        iu = ibu
      end if
      if (jc == 1) then
        jl = jb
        ju = jb + jh - 1
      else
        jl = je - jh + 1
        ju = je
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = kb
        ku = kbu
      else
        kl = kel
        ku = ke
      end if
#endif /* NDIMS == 3 */

! copy the edge region to the output array
!
#if NDIMS == 2
      qb(1:nv,1:ng,1:jh,1:km) = qn(1:nv,il:iu,jl:ju, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      qb(1:nv,1:ng,1:jh,1:ng) = qn(1:nv,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

! calculate half size
!
      kh = kn / 2

! prepare source edge region indices
!
      if (ic == 1) then
        il = iel
        iu = ie
      else
        il = ib
        iu = ibu
      end if
      if (jc == 1) then
        jl = jel
        ju = je
      else
        jl = jb
        ju = jbu
      end if
      if (kc == 1) then
        kl = kb
        ku = kb + kh - 1
      else
        kl = ke - kh + 1
        ku = ke
      end if

! copy the edge region to the output array
!
      qb(1:nv,1:ng,1:ng,1:kh) = qn(1:nv,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

    end select

!-------------------------------------------------------------------------------
!
  end subroutine block_edge_copy
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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
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
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half size
!
      ih = in / 2

! prepare indices for the edge region
!
      il = ib
      ip = il + 1
      iu = ie
      if (jc == 1) then
        jl = je - nd + 1
        jp = jl + 1
        ju = je
      else
        jl = jb
        jp = jl + 1
        ju = jb + nd - 1
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = ke - nd + 1
        kp = kl + 1
        ku = ke
      else
        kl = kb
        kp = kl + 1
        ku = kb + nd - 1
      end if
#endif /* NDIMS == 3 */

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

! prepare indices for the edge region
!
      if (ic == 1) then
        il = ie - nd + 1
        ip = il + 1
        iu = ie
      else
        il = ib
        ip = il + 1
        iu = ib + nd - 1
      end if
      jl = jb
      jp = jl + 1
      ju = je
#if NDIMS == 3
      if (kc == 1) then
        kl = ke - nd + 1
        kp = kl + 1
        ku = ke
      else
        kl = kb
        kp = kl + 1
        ku = kb + nd - 1
      end if
#endif /* NDIMS == 3 */

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

! prepare indices for the edge region
!
      if (ic == 1) then
        il = ie - nd + 1
        ip = il + 1
        iu = ie
      else
        il = ib
        ip = il + 1
        iu = ib + nd - 1
      end if
      if (jc == 1) then
        jl = je - nd + 1
        jp = jl + 1
        ju = je
      else
        jl = jb
        jp = jl + 1
        ju = jb + nd - 1
      end if
      kl = kb
      kp = kl + 1
      ku = ke

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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use equations      , only : nv
    use interpolations , only : limiter

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
    real(kind=8) :: dqx, dqy, dqz
    real(kind=8) :: dq1, dq2, dq3, dq4
!
!-------------------------------------------------------------------------------
!
! process depending on the direction
!
    select case(nc)
    case(1)

! calculate half size
!
      ih = in / 2

! prepare indices for the edge region
!
      if (ic == 0) then
        il = ib
        iu = ib + ih + nh - 1
      else
        il = ie - ih - nh + 1
        iu = ie
      end if
      if (jc == 1) then
        jl = je - nh + 1
        ju = je
      else
        jl = jb
        ju = jb + nh - 1
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = ke - nh + 1
        ku = ke
      else
        kl = kb
        ku = kb + nh - 1
      end if
#endif /* NDIMS == 3 */

    case(2)

! calculate half size
!
      jh = jn / 2

! prepare indices for the edge region
!
      if (ic == 1) then
        il = ie - nh + 1
        iu = ie
      else
        il = ib
        iu = ib + nh - 1
      end if
      if (jc == 0) then
        jl = jb
        ju = jb + jh + nh - 1
      else
        jl = je - jh - nh + 1
        ju = je
      end if
#if NDIMS == 3
      if (kc == 1) then
        kl = ke - nh + 1
        ku = ke
      else
        kl = kb
        ku = kb + nh - 1
      end if
#endif /* NDIMS == 3 */

#if NDIMS == 3
    case(3)

! calculate half size
!
      kh = kn / 2

! prepare indices for the edge region
!
      if (ic == 1) then
        il = ie - nh + 1
        iu = ie
      else
        il = ib
        iu = ib + nh - 1
      end if
      if (jc == 1) then
        jl = je - nh + 1
        ju = je
      else
        jl = jb
        ju = jb + nh - 1
      end if
      if (kc == 0) then
        kl = kb
        ku = kb + kh + nh - 1
      else
        kl = ke - kh - nh + 1
        ku = ke
      end if
#endif /* NDIMS == 3 */

    end select

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
            dql = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dqx = limiter(0.25d+00, dql, dqr)

            dql = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dqy = limiter(0.25d+00, dql, dqr)

#if NDIMS == 3
            dql = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dqz = limiter(0.25d+00, dql, dqr)
#endif /* NDIMS == 3 */

#if NDIMS == 2
! calculate the derivative terms
!
            dq1 = dqx + dqy
            dq2 = dqx - dqy

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
            dq1 = dqx + dqy + dqz
            dq2 = dqx - dqy - dqz
            dq3 = dqx - dqy + dqz
            dq4 = dqx + dqy - dqz

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
! subroutine BLOCK_CORNER_COPY:
! ----------------------------
!
!   Subroutine returns the corner boundary region by copying the corresponding
!   region from the provided input variable array.
!
!   Arguments:
!
!     ic, jc, kc - the corner position;
!     qn         - the input neighbor variable array;
!     qb         - the output corner boundary array;
!
!===============================================================================
!
  subroutine block_corner_copy(ic, jc, kc, qn, qb)

! import external procedures and variables
!
    use coordinates    , only : ng
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use coordinates    , only : ibu, jbu, kbu
    use coordinates    , only : iel, jel, kel
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

! local indices
!
    integer :: il, jl, kl
    integer :: iu, ju, ku
!
!-------------------------------------------------------------------------------
!
! prepare indices for the corner region
!
    if (ic == 1) then
      il = iel
      iu = ie
    else
      il = ib
      iu = ibu
    end if
    if (jc == 1) then
      jl = jel
      ju = je
    else
      jl = jb
      ju = jbu
    end if
#if NDIMS == 3
    if (kc == 1) then
      kl = kel
      ku = ke
    else
      kl = kb
      ku = kbu
    end if
#endif /* NDIMS == 3 */

! copy the corner region to the output array
!
#if NDIMS == 2
    qb(1:nv,1:ng,1:ng,1:km) = qn(1:nv,il:iu,jl:ju, 1:km)
#endif /* NDIMS == 2 */
#if NDIMS == 3
    qb(1:nv,1:ng,1:ng,1:ng) = qn(1:nv,il:iu,jl:ju,kl:ku)
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
  end subroutine block_corner_copy
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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
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
    if (ic == 1) then
      il = ie - nd + 1
      ip = il + 1
      iu = ie
    else
      il = ib
      ip = il + 1
      iu = ib + nd - 1
    end if
    if (jc == 1) then
      jl = je - nd + 1
      jp = jl + 1
      ju = je
    else
      jl = jb
      jp = jl + 1
      ju = jb + nd - 1
    end if
#if NDIMS == 3
    if (kc == 1) then
      kl = ke - nd + 1
      kp = kl + 1
      ku = ke
    else
      kl = kb
      kp = kl + 1
      ku = kb + nd - 1
    end if
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
    use coordinates    , only : ib , jb , kb
    use coordinates    , only : ie , je , ke
    use equations      , only : nv
    use interpolations , only : limiter

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
    real(kind=8) :: dqx, dqy, dqz
    real(kind=8) :: dq1, dq2, dq3, dq4
!
!-------------------------------------------------------------------------------
!
! prepare indices for the corner region
!
    if (ic == 1) then
      il = ie - nh + 1
      iu = ie
    else
      il = ib
      iu = ib + nh - 1
    end if
    if (jc == 1) then
      jl = je - nh + 1
      ju = je
    else
      jl = jb
      ju = jb + nh - 1
    end if
#if NDIMS == 3
    if (kc == 1) then
      kl = ke - nh + 1
      ku = ke
    else
      kl = kb
      ku = kb + nh - 1
    end if
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
            dql = qn(p,i  ,j,k) - qn(p,im1,j,k)
            dqr = qn(p,ip1,j,k) - qn(p,i  ,j,k)
            dqx = limiter(0.25d+00, dql, dqr)

            dql = qn(p,i,j  ,k) - qn(p,i,jm1,k)
            dqr = qn(p,i,jp1,k) - qn(p,i,j  ,k)
            dqy = limiter(0.25d+00, dql, dqr)

#if NDIMS == 3
            dql = qn(p,i,j,k  ) - qn(p,i,j,km1)
            dqr = qn(p,i,j,kp1) - qn(p,i,j,k  )
            dqz = limiter(0.25d+00, dql, dqr)
#endif /* NDIMS == 3 */

#if NDIMS == 2
! calculate the derivative terms
!
            dq1 = dqx + dqy
            dq2 = dqx - dqy

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
            dq1 = dqx + dqy + dqz
            dq2 = dqx - dqy - dqz
            dq3 = dqx - dqy + dqz
            dq4 = dqx + dqy - dqz

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
