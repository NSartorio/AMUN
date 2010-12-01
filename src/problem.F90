!!******************************************************************************
!!
!! module: problem - handling the initial problem definition
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
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
!!******************************************************************************
!!
!
module problem

  implicit none

  contains
!
!===============================================================================
!
! init_domain: subroutine initializes the domain for a given problem
!
!===============================================================================
!
  subroutine init_domain

    use config, only : problem
!
!-------------------------------------------------------------------------------
!
    select case(trim(problem))
#if NDIMS == 2
    case("blast")
      call domain_blast()
#endif /* NDIMS == 2 */
    case default
      call domain_default()
    end select

!-------------------------------------------------------------------------------
!
  end subroutine init_domain
!
!===============================================================================
!
! init_problem: subroutine initializes the variables according to
!               the studied problem
!
!===============================================================================
!
  subroutine init_problem(pblock)

    use blocks, only : block_data
    use config, only : problem

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local arguments
!
    type(block_data), pointer :: pb
!
!-------------------------------------------------------------------------------
!
    pb => pblock

    select case(trim(problem))
    case("blast")
      call init_blast(pb)
    case("implosion")
      call init_implosion(pb)
    case("binaries")
      call init_binaries(pb)
    end select

    nullify(pb)

!-------------------------------------------------------------------------------
!
  end subroutine init_problem
#ifdef SHAPE
!
!===============================================================================
!
! shapes: subroutine resets the update for a give shape
!
!===============================================================================
!
  subroutine update_shapes(pblock, du)

    use blocks, only : block_data
    use config, only : problem

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock
    real, dimension(:,:,:,:) , intent(inout) :: du

! local arguments
!
    type(block_data), pointer :: pb
!
!-------------------------------------------------------------------------------
!
    pb => pblock

    select case(trim(problem))
    case("binaries")
      call shape_binaries(pb, du)
    case default
    end select

    nullify(pb)

!-------------------------------------------------------------------------------
!
  end subroutine update_shapes
#endif /* SHAPE */
!
!===============================================================================
!
! domain_default: subroutine initializes the default domain of 2x2 blocks in
!                'N' configuration
!
!===============================================================================
!
  subroutine domain_default

    use blocks, only : block_meta, block_data, append_metablock                &
                     , append_datablock, associate_blocks, metablock_setleaf   &
                     , metablock_setconfig, metablock_setlevel                 &
                     , metablock_set_coord, metablock_setbounds                &
                     , nsides, nfaces
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , xmin, xmax, ymin, ymax, zmin, zmax

    implicit none

! local variables
!
    integer :: i, j

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create root meta blocks
!
    call append_metablock(pmeta)

! mark block as a leaf
!
    call metablock_setleaf(pmeta)

! set block config flag
!
    call metablock_setconfig(pmeta, 12)

! set block level
!
    call metablock_setlevel(pmeta, 1)

! set block coordinates
!
    call metablock_set_coord(pmeta, 0, 0, 0)

! set block bounds
!
    call metablock_setbounds(pmeta, xmin, xmax, ymin, ymax, zmin, zmax)

! if periodic boundary conditions set all neighbors to itself
!
    if (xlbndry .eq. 'periodic' .and. xubndry .eq. 'periodic') then
      do j = 1, nfaces
        do i = 1, nsides
          pmeta%neigh(1,i,j)%ptr => pmeta
        end do
      end do
    end if
    if (ylbndry .eq. 'periodic' .and. yubndry .eq. 'periodic') then
      do j = 1, nfaces
        do i = 1, nsides
          pmeta%neigh(2,i,j)%ptr => pmeta
        end do
      end do
    end if
#if NDIMS == 3
    if (zlbndry .eq. 'periodic' .and. zubndry .eq. 'periodic') then
      do j = 1, nfaces
        do i = 1, nsides
          pmeta%neigh(3,i,j)%ptr => pmeta
        end do
      end do
    end if
#endif /* NDIMS == 3 */

! create root data block
!
    call append_datablock(pdata)

! associate root blocks
!
    call associate_blocks(pmeta, pdata)
!
!-------------------------------------------------------------------------------
!
  end subroutine domain_default
!
!===============================================================================
!
! domain_default: subroutine initializes the default domain of 2x2 blocks in
!                'N' configuration
!
!===============================================================================
!
  subroutine domain_blast

    use blocks, only : block_meta, block_data, pointer_meta, append_metablock  &
                     , append_datablock, associate_blocks, metablock_setleaf   &
                     , metablock_setconfig, metablock_setlevel                 &
                     , metablock_setbounds, metablock_set_coord                &
                     , nsides, nfaces, res
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , xmin, xmax, ymin, ymax, zmin, zmax, rdims, maxlev, ncells

    implicit none

! local variables
!
    integer :: i, j, ih, jl, ju
    real    :: xm, yl, yu

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata

! local pointer array
!
    type(pointer_meta)        :: block_array(2,3)
!
!-------------------------------------------------------------------------------
!
!! TO DO:
!!
!! create a general way to initiate NxM blocks at level 1
!!
!!
! create root meta blocks
!
    call append_metablock(block_array(1,1)%ptr)
    call append_metablock(block_array(2,1)%ptr)
    call append_metablock(block_array(2,2)%ptr)
    call append_metablock(block_array(1,2)%ptr)
    call append_metablock(block_array(1,3)%ptr)
    call append_metablock(block_array(2,3)%ptr)

! mark block as a leaf
!
    call metablock_setleaf(block_array(1,1)%ptr)
    call metablock_setleaf(block_array(2,1)%ptr)
    call metablock_setleaf(block_array(2,2)%ptr)
    call metablock_setleaf(block_array(1,2)%ptr)
    call metablock_setleaf(block_array(1,3)%ptr)
    call metablock_setleaf(block_array(2,3)%ptr)

! set block config flag
!
    call metablock_setconfig(block_array(1,1)%ptr, 12)
    call metablock_setconfig(block_array(2,1)%ptr, 13)
    call metablock_setconfig(block_array(2,2)%ptr, 13)
    call metablock_setconfig(block_array(1,2)%ptr, 43)
    call metablock_setconfig(block_array(1,3)%ptr, 12)
    call metablock_setconfig(block_array(2,3)%ptr, 12)

! set block level
!
    call metablock_setlevel(block_array(1,1)%ptr, 1)
    call metablock_setlevel(block_array(2,1)%ptr, 1)
    call metablock_setlevel(block_array(2,2)%ptr, 1)
    call metablock_setlevel(block_array(1,2)%ptr, 1)
    call metablock_setlevel(block_array(1,3)%ptr, 1)
    call metablock_setlevel(block_array(2,3)%ptr, 1)

! set boundary conditions
!
    block_array(1,1)%ptr%neigh(1,1,1)%ptr => block_array(2,1)%ptr
    block_array(1,1)%ptr%neigh(1,1,2)%ptr => block_array(2,1)%ptr
    block_array(1,1)%ptr%neigh(1,2,1)%ptr => block_array(2,1)%ptr
    block_array(1,1)%ptr%neigh(1,2,2)%ptr => block_array(2,1)%ptr
    block_array(1,1)%ptr%neigh(2,1,1)%ptr => block_array(1,3)%ptr
    block_array(1,1)%ptr%neigh(2,1,2)%ptr => block_array(1,3)%ptr
    block_array(1,1)%ptr%neigh(2,2,1)%ptr => block_array(1,2)%ptr
    block_array(1,1)%ptr%neigh(2,2,2)%ptr => block_array(1,2)%ptr

    block_array(1,2)%ptr%neigh(1,1,1)%ptr => block_array(2,2)%ptr
    block_array(1,2)%ptr%neigh(1,1,2)%ptr => block_array(2,2)%ptr
    block_array(1,2)%ptr%neigh(1,2,1)%ptr => block_array(2,2)%ptr
    block_array(1,2)%ptr%neigh(1,2,2)%ptr => block_array(2,2)%ptr
    block_array(1,2)%ptr%neigh(2,1,1)%ptr => block_array(1,1)%ptr
    block_array(1,2)%ptr%neigh(2,1,2)%ptr => block_array(1,1)%ptr
    block_array(1,2)%ptr%neigh(2,2,1)%ptr => block_array(1,3)%ptr
    block_array(1,2)%ptr%neigh(2,2,2)%ptr => block_array(1,3)%ptr

    block_array(1,3)%ptr%neigh(1,1,1)%ptr => block_array(2,3)%ptr
    block_array(1,3)%ptr%neigh(1,1,2)%ptr => block_array(2,3)%ptr
    block_array(1,3)%ptr%neigh(1,2,1)%ptr => block_array(2,3)%ptr
    block_array(1,3)%ptr%neigh(1,2,2)%ptr => block_array(2,3)%ptr
    block_array(1,3)%ptr%neigh(2,1,1)%ptr => block_array(1,2)%ptr
    block_array(1,3)%ptr%neigh(2,1,2)%ptr => block_array(1,2)%ptr
    block_array(1,3)%ptr%neigh(2,2,1)%ptr => block_array(1,1)%ptr
    block_array(1,3)%ptr%neigh(2,2,2)%ptr => block_array(1,1)%ptr

    block_array(2,1)%ptr%neigh(1,1,1)%ptr => block_array(1,1)%ptr
    block_array(2,1)%ptr%neigh(1,1,2)%ptr => block_array(1,1)%ptr
    block_array(2,1)%ptr%neigh(1,2,1)%ptr => block_array(1,1)%ptr
    block_array(2,1)%ptr%neigh(1,2,2)%ptr => block_array(1,1)%ptr
    block_array(2,1)%ptr%neigh(2,1,1)%ptr => block_array(2,3)%ptr
    block_array(2,1)%ptr%neigh(2,1,2)%ptr => block_array(2,3)%ptr
    block_array(2,1)%ptr%neigh(2,2,1)%ptr => block_array(2,2)%ptr
    block_array(2,1)%ptr%neigh(2,2,2)%ptr => block_array(2,2)%ptr

    block_array(2,2)%ptr%neigh(1,1,1)%ptr => block_array(1,2)%ptr
    block_array(2,2)%ptr%neigh(1,1,2)%ptr => block_array(1,2)%ptr
    block_array(2,2)%ptr%neigh(1,2,1)%ptr => block_array(1,2)%ptr
    block_array(2,2)%ptr%neigh(1,2,2)%ptr => block_array(1,2)%ptr
    block_array(2,2)%ptr%neigh(2,1,1)%ptr => block_array(2,1)%ptr
    block_array(2,2)%ptr%neigh(2,1,2)%ptr => block_array(2,1)%ptr
    block_array(2,2)%ptr%neigh(2,2,1)%ptr => block_array(2,3)%ptr
    block_array(2,2)%ptr%neigh(2,2,2)%ptr => block_array(2,3)%ptr

    block_array(2,3)%ptr%neigh(1,1,1)%ptr => block_array(1,3)%ptr
    block_array(2,3)%ptr%neigh(1,1,2)%ptr => block_array(1,3)%ptr
    block_array(2,3)%ptr%neigh(1,2,1)%ptr => block_array(1,3)%ptr
    block_array(2,3)%ptr%neigh(1,2,2)%ptr => block_array(1,3)%ptr
    block_array(2,3)%ptr%neigh(2,1,1)%ptr => block_array(2,2)%ptr
    block_array(2,3)%ptr%neigh(2,1,2)%ptr => block_array(2,2)%ptr
    block_array(2,3)%ptr%neigh(2,2,1)%ptr => block_array(2,1)%ptr
    block_array(2,3)%ptr%neigh(2,2,2)%ptr => block_array(2,1)%ptr

! prepare the coordinates
!
    ih = res(1)
    jl = res(1)
    ju = res(1) + jl

! set the coordinates
!
    call metablock_set_coord(block_array(1,1)%ptr,  0,  0, 0)
    call metablock_set_coord(block_array(2,1)%ptr, ih,  0, 0)
    call metablock_set_coord(block_array(1,2)%ptr,  0, jl, 0)
    call metablock_set_coord(block_array(2,2)%ptr, ih, jl, 0)
    call metablock_set_coord(block_array(1,3)%ptr,  0, ju, 0)
    call metablock_set_coord(block_array(2,3)%ptr, ih, ju, 0)

! calculate bounds
!
    xm = 0.5 * (xmin + xmax)
    yl = (ymax - ymin) / 3.0 + ymin
    yu = (ymax - ymin) / 3.0 + yl

! set block bounds
!
    call metablock_setbounds(block_array(1,1)%ptr, xmin, xm  , ymin, yl  , zmin, zmax)
    call metablock_setbounds(block_array(2,1)%ptr, xm  , xmax, ymin, yl  , zmin, zmax)
    call metablock_setbounds(block_array(2,2)%ptr, xm  , xmax, yl  , yu  , zmin, zmax)
    call metablock_setbounds(block_array(1,2)%ptr, xmin, xm  , yl  , yu  , zmin, zmax)
    call metablock_setbounds(block_array(1,3)%ptr, xmin, xm  , yu  , ymax, zmin, zmax)
    call metablock_setbounds(block_array(2,3)%ptr, xm  , xmax, yu  , ymax, zmin, zmax)

! create root data block
!
    call append_datablock(pdata)
    call associate_blocks(block_array(1,1)%ptr, pdata)
    call append_datablock(pdata)
    call associate_blocks(block_array(2,1)%ptr, pdata)
    call append_datablock(pdata)
    call associate_blocks(block_array(2,2)%ptr, pdata)
    call append_datablock(pdata)
    call associate_blocks(block_array(1,2)%ptr, pdata)
    call append_datablock(pdata)
    call associate_blocks(block_array(1,3)%ptr, pdata)
    call append_datablock(pdata)
    call associate_blocks(block_array(2,3)%ptr, pdata)

! set the block dimensions for the lowest level
!
    rdims(1) = 2
    rdims(2) = 3
!
!-------------------------------------------------------------------------------
!
  end subroutine domain_blast
!
!===============================================================================
!
! init_blast: subroutine initializes the variables for the blast problem
!
!===============================================================================
!
  subroutine init_blast(pblock)

    use blocks   , only : block_data
    use config   , only : in, jn, kn, im, jm, km, ng                           &
                        , gamma, csnd2, rcut, dens, dnrat
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt, idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: r, dx, dy, dz, pr

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(km)     :: z
    real, dimension(nvr,im) :: q, u
!
!-------------------------------------------------------------------------------
!
! calculate parameters
!
#ifdef ADI
    pr = dens * csnd2 / gamma
#endif /* ADI */

! calculate cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%meta%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%meta%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0
#endif /* NDIMS == 3 */

! set variables
!
    q(idn,:) = dens
    q(ivx,:) = 0.0d0
    q(ivy,:) = 0.0d0
    q(ivz,:) = 0.0d0
#ifdef MHD
    q(ibx,:) = 0.70710678118654752440
    q(iby,:) = 0.70710678118654752440
    q(ibz,:) = 0.0d0
#ifdef FLUXCT
    q(icx,:) = 0.70710678118654752440
    q(icy,:) = 0.70710678118654752440
    q(icz,:) = 0.0d0
#endif /* FLUXCT */
#endif /* MHD */

! set initial pressure
!
    do k = 1, km
      do j = 1, jm

! fill out the pressure profile
!
        do i = 1, im

          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

#ifdef ISO
          if (r .lt. rcut) then
            q(idn,i) = dens * dnrat
          else
            q(idn,i) = dens
          end if
#endif /* ISO */
#ifdef ADI
          if (r .lt. rcut) then
            q(ipr,i) = pr * dnrat
          else
            q(ipr,i) = pr
          end if
#endif /* ADI */

        end do

! convert primitive variables to conserved
!
        call prim2cons(im, q, u)

! copy conservative variables to the current block
!
        pblock%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_blast
!
!===============================================================================
!
! init_implosion: subroutine initializes variables for the implosion problem
!
!===============================================================================
!
  subroutine init_implosion(pblock)

    use blocks   , only : block_data
    use config   , only : in, jn, kn, im, jm, km, ng, dens, pres, rmid, gammam1i
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: dx, dy, dxh, dyh, rc, rl, ru, ds, dl, dr

! local arrays
!
    real, dimension(:), allocatable :: x, y
!
!-------------------------------------------------------------------------------
!
! allocate coordinates
!
    allocate(x(im))
    allocate(y(jm))

! calculate cell sizes
!
    dx  = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy  = (pblock%meta%ymax - pblock%meta%ymin) / jn
    dxh = 0.5d0 * dx
    dyh = 0.5d0 * dy
    ds  = dx * dy

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng) * dx - dxh + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng) * dy - dyh + pblock%meta%ymin

! set variables
!
    pblock%u(idn,:,:,:) = dens
    pblock%u(imx,:,:,:) = 0.0d0
    pblock%u(imy,:,:,:) = 0.0d0
    pblock%u(imz,:,:,:) = 0.0d0

! set initial pressure
!
    do j = 1, jm
      do i = 1, im
        rc = x(i) + y(j)
        rl = rc - dxh - dyh
        ru = rc + dxh + dyh

        if (ru .le. rmid) then
          pblock%u(idn,i,j,:) = 0.125d0
#ifdef ADI
          pblock%u(ien,i,j,:) = gammam1i * 0.140d0
#endif /* ADI */
        else if (rl .ge. rmid) then
          pblock%u(idn,i,j,:) = 1.0d0
#ifdef ADI
          pblock%u(ien,i,j,:) = gammam1i
#endif /* ADI */
        else
          if (rc .ge. rmid) then
            dl = 0.125d0 * (rmid - rl)**2 / ds
            dr = 1.0d0 - dl
          else
            dr = 0.125d0 * (ru - rmid)**2 / ds
            dl = 1.0d0 - dr
          endif
          pblock%u(idn,i,j,:) = 0.125d0 * dl + dr
#ifdef ADI
          pblock%u(ien,i,j,:) = gammam1i * (0.140d0 * dl + dr)
#endif /* ADI */
        endif
      end do
    end do

! deallocate coordinates
!
    deallocate(x)
    deallocate(y)

!-------------------------------------------------------------------------------
!
  end subroutine init_implosion
!
!===============================================================================
!
! init_binaries: subroutine initializes the variables for the binary star
!                problem
!
!===============================================================================
!
  subroutine init_binaries(pblock)

    use blocks   , only : block_data
    use config   , only : ng, in, jn, kn, im, jm, km, dens, pres, dnfac, dnrat &
                        , x1c, y1c, z1c, r1c, x2c, y2c, z2c, r2c, v1ini, v2ini &
                        , csnd2, gamma, gammam1i
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k
    real    :: dx, dy, dz, dnamb, enamb, ekin
    real    :: dnstar1, enstar1, x1l, y1l, z1l, r1
    real    :: dnstar2, enstar2, x2l, y2l, z2l, r2

! local arrays
!
    real, dimension(:), allocatable :: x, y, z
!
!-------------------------------------------------------------------------------
!
! calculate pressure from sound speed
!
#ifdef ISO
    pres = csnd2 * dens
#endif /* ISO */
#ifdef ADI
    pres = csnd2 * dens / gamma
#endif /* ADI */

! calculate parameters
!
    dnamb   = dens
    dnstar2 = dnamb*dnfac
    dnstar1 = dnstar2*dnrat
    enamb   = gammam1i*pres
    enstar1 = enamb*dnfac
    enstar2 = enstar1/dnrat

! allocate coordinates
!
    allocate(x(im))
    allocate(y(jm))
    allocate(z(km))

! calculate cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%meta%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%meta%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0

    z1l  = 0.0
    z2l  = 0.0
#endif /* NDIMS == 3 */

! set variables
!
    pblock%u(idn,:,:,:) = dnamb
    pblock%u(imx,:,:,:) = 0.0d0
    pblock%u(imy,:,:,:) = 0.0d0
    pblock%u(imz,:,:,:) = 0.0d0
#ifdef ADI
    pblock%u(ien,:,:,:) = enamb
#endif /* ADI */

! set initial pressure
!
    do k = 1, km
#if NDIMS == 3
      z1l = z(k) - z1c
      z2l = z(k) - z2c
#endif /* NDIMS == 3 */

      do j = 1, jm
        y1l = y(j) - y1c
        y2l = y(j) - y2c

        do i = 1, im
          x1l = x(i) - x1c
          x2l = x(i) - x2c

          r1 = sqrt(x1l**2 + y1l**2 + z1l**2)
          r2 = sqrt(x2l**2 + y2l**2 + z2l**2)

          if (r1 .le. r1c) then
            pblock%u(idn,i,j,k) = dnstar1
            pblock%u(imx,i,j,k) = dnstar1*v1ini*x1l
            pblock%u(imy,i,j,k) = dnstar1*v1ini*y1l
#if NDIMS == 3
            pblock%u(imz,i,j,k) = dnstar1*v1ini*z1l
#endif /* NDIMS == 3 */
            ekin = 0.5 * dnstar1 * v1ini**2 * (x1l**2 + y1l**2 + z1l**2)
#ifdef ADI
            pblock%u(ien,i,j,k) = enstar1 + ekin
#endif /* ADI */
          endif

          if (r2 .le. r2c) then
            pblock%u(idn,i,j,k) = dnstar2
            pblock%u(imx,i,j,k) = dnstar2*v2ini*x2l
            pblock%u(imy,i,j,k) = dnstar2*v2ini*y2l
#if NDIMS == 3
            pblock%u(imz,i,j,k) = dnstar2*v2ini*z2l
#endif /* NDIMS == 3 */
            ekin = 0.5 * dnstar2 * v2ini**2 * (x2l**2 + y2l**2 + z2l**2)
#ifdef ADI
            pblock%u(ien,i,j,k) = enstar2 + ekin
#endif /* ADI */
          endif
        end do
      end do
    end do

! deallocate coordinates
!
    deallocate(x)
    deallocate(y)
    deallocate(z)

!-------------------------------------------------------------------------------
!
  end subroutine init_binaries
#ifdef SHAPE
!
!===============================================================================
!
! shape_binaries: subroutine initializes the variables for the binary star
!                problem
!
!===============================================================================
!
  subroutine shape_binaries(pblock, du)

    use blocks   , only : block_data
    use config   , only : ng, in, jn, kn, im, jm, km, dens, pres, dnfac, dnrat &
                        , x1c, y1c, z1c, r1c, x2c, y2c, z2c, r2c               &
                        , csnd2, gamma, gammam1i
    use variables, only : idn, imx, imy, imz, ien

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock
    real, dimension(:,:,:,:) , intent(inout) :: du

! local variables
!
    integer :: i, j, k
    real    :: dx, dy, dz
    real    :: x1l, y1l, z1l, r1
    real    :: x2l, y2l, z2l, r2

! local arrays
!
    real, dimension(:), allocatable :: x, y, z
!
!-------------------------------------------------------------------------------
!
! allocate coordinates
!
    allocate(x(im))
    allocate(y(jm))
    allocate(z(km))

! calculate cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%meta%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%meta%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0

    z1l  = 0.0
    z2l  = 0.0
#endif /* NDIMS == 3 */

! reset update
!
    do k = 1, km
#if NDIMS == 3
      z1l = z(k) - z1c
      z2l = z(k) - z2c
#endif /* NDIMS == 3 */

      do j = 1, jm
        y1l = y(j) - y1c
        y2l = y(j) - y2c

        do i = 1, im
          x1l = x(i) - x1c
          x2l = x(i) - x2c

          r1 = sqrt(x1l**2 + y1l**2 + z1l**2)
          r2 = sqrt(x2l**2 + y2l**2 + z2l**2)

          if (r1 .le. r1c) then
            du(:,i,j,k) = 0.0
          endif

          if (r2 .le. r2c) then
            du(:,i,j,k) = 0.0
          endif
        end do
      end do
    end do

! deallocate coordinates
!
    deallocate(x)
    deallocate(y)
    deallocate(z)

!-------------------------------------------------------------------------------
!
  end subroutine shape_binaries
#endif /* SHAPE */
!
!===============================================================================
!
! check_ref: function checks refinement criterium and returns +1 if
!            the criterium fullfiled and block is selected for
!            refinement, 0 there is no need for refinement, and -1 if
!            block is selected for refinement
!
!===============================================================================
!
  function check_ref(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km, ibl, ieu, jbl, jeu, kbl, keu, gammam1i   &
                        , crefmin, crefmax
    use variables, only : idn, imx, imy, imz, nvr
#ifdef ADI
    use variables, only : ien
#endif /* ADI */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! return variable
!
    integer(kind=1) :: check_ref

! local variables
!
    integer :: i, j, k
    real    :: dpmax, vx, vy, vz, en, ek, ei
    real    :: dnl, dnr, prl, prr, ddndx, ddndy, ddndz, ddn, dprdx, dprdy      &
             , dprdz, dpr

! local arrays
!
    real, dimension(im,jm,km) :: dn, pr
!
!-------------------------------------------------------------------------------
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          dn(i,j,k) = pblock%u(idn,i,j,k)
#ifdef ADI
          vx = pblock%u(imx,i,j,k) / dn(i,j,k)
          vy = pblock%u(imy,i,j,k) / dn(i,j,k)
          vz = pblock%u(imz,i,j,k) / dn(i,j,k)
          en = pblock%u(ien,i,j,k)

          ek = 0.5 * dn(i,j,k) * (vx*vx + vy*vy + vz*vz)
          ei = en - ek
          pr(i,j,k) = gammam1i * ei
#endif /* ADI */
        end do
      end do
    end do

! check gradient of pressure
!
    dpmax = 0.0d0

    do k = kbl, keu
      do j = jbl, jeu
        do i = ibl, ieu

          dnl = dn(i-1,j,k)
          dnr = dn(i+1,j,k)
          ddndx = abs(dnr-dnl)/(dnr+dnl)
          dnl = dn(i,j-1,k)
          dnr = dn(i,j+1,k)
          ddndy = abs(dnr-dnl)/(dnr+dnl)
#if NDIMS == 3
          dnl = dn(i,j,k-1)
          dnr = dn(i,j,k+1)
          ddndz = abs(dnr-dnl)/(dnr+dnl)

          ddn = sqrt(ddndx**2 + ddndy**2 + ddndz**2)
#else /* NDIMS == 3 */

          ddn = sqrt(ddndx**2 + ddndy**2)
#endif /* NDIMS == 3 */

          dpmax = max(dpmax, ddn)

#ifdef ADI
          prl = pr(i-1,j,k)
          prr = pr(i+1,j,k)
          dprdx = abs(prr-prl)/(prr+prl)
          prl = pr(i,j-1,k)
          prr = pr(i,j+1,k)
          dprdy = abs(prr-prl)/(prr+prl)
#if NDIMS == 3
          prl = pr(i,j,k-1)
          prr = pr(i,j,k+1)
          dprdz = abs(prr-prl)/(prr+prl)

          dpr = sqrt(dprdx**2 + dprdy**2 + dprdz**2)
#else /* NDIMS == 3 */

          dpr = sqrt(dprdx**2 + dprdy**2)
#endif /* NDIMS == 3 */

          dpmax = max(dpmax, dpr)
#endif /* ADI */

        end do
      end do
    end do

! check condition
!
    check_ref = 0

    if (dpmax .ge. crefmax) then
      check_ref =  1
    end if
    if (dpmax .le. crefmin) then
      check_ref = -1
    end if

    return

!-------------------------------------------------------------------------------
!
  end function check_ref

!===============================================================================
!
end module
