!!******************************************************************************
!!
!! module: problem - handling the initial problem definition
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
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
  subroutine init_domain()

    use config, only : problem
!
!-------------------------------------------------------------------------------
!
    select case(trim(problem))
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
    case("reconnection")
      call init_reconnection(pb)
    case("multi_current_sheet")
      call init_multi_current_sheet(pb)
    case("turbulence")
      call init_turbulence(pb)
    case("orszag_tang")
      call init_orszag_tang(pb)
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
! domain_default: subroutine initializes the default domain of N1xN2xN3 blocks
!                 in the proper configuration; the dimensions N1xN2xN3 are set
!                 in variable rdims
!
!===============================================================================
!
  subroutine domain_default()

    use blocks, only : pointer_meta, block_meta, block_data, append_metablock  &
                     , append_datablock, associate_blocks, metablock_setleaf   &
                     , metablock_setconfig, metablock_setlevel                 &
                     , metablock_set_coord, metablock_setbounds                &
                     , nsides, nfaces, res
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , xmin, xmax, ymin, ymax, zmin, zmax, rdims

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
    n = product(rdims(:))

!! PREPARE BLOCK CONFIGURATION ARRAY
!!
! allocate the configuration array
!
    allocate(cfg(rdims(1),rdims(2),rdims(3)))

! set the block configurations
!
    cfg(1:rdims(1),1:rdims(2):2,1:rdims(3):2) = 12

    if (rdims(2) .gt. 1) then
      cfg(1:rdims(1),2:rdims(2):2,1:rdims(3):2) = 43
      cfg(  rdims(1),1:rdims(2)  ,1:rdims(3):2) = 13
    end if

    if (rdims(3) .gt. 1) then
      cfg(1:rdims(1),1:rdims(2):2,2:rdims(3):2) = 65
      if (rdims(2) .gt. 1) then
        cfg(1:rdims(1),2:rdims(2):2,2:rdims(3):2) = 78
        cfg(  rdims(1),1:rdims(2)  ,2:rdims(3):2) = 75
      end if
      if (rdims(1) .eq. 1 .or. mod(rdims(2),2) .eq. 1) then
        cfg(  rdims(1),  rdims(2)  ,1:rdims(3)  ) = 15
      else
        cfg(  1       ,  rdims(2)  ,1:rdims(3)  ) = 48
      end if
    end if

!! ALLOCATE AND GENERATE META BLOCK CHAIN AND SET BLOCK CONFIGURATIONS
!!
! allocate the block pointer array
!
    allocate(block_array(rdims(1),rdims(2),rdims(3)))

! generate the gray code for a given configuration and link the block in
! the proper order
!
    loc(:) = (/ 0, 0, 0 /)
    del(:) = (/ 1, 1, 1 /)

    p = 1
    do k = 1, rdims(3)
      if (del(3) .eq. 1) loc(3) = loc(3) + del(3)
      do j = 1, rdims(2)
        if (del(2) .eq. 1) loc(2) = loc(2) + del(2)
        do i = 1, rdims(1)
          if (del(1) .eq. 1) loc(1) = loc(1) + del(1)

! append a new metablock
!
          call append_metablock(block_array(loc(1),loc(2),loc(3))%ptr)

! set the configuration type
!
          call metablock_setconfig(block_array(loc(1),loc(2),loc(3))%ptr       &
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
    xl = (xmax - xmin) / rdims(1)
    yl = (ymax - ymin) / rdims(2)
    zl = (zmax - zmin) / rdims(3)

! fill out block structure fields
!
    do k = 1, rdims(3)

! claculate the block position along Z
!
      kl  = (k - 1) * res(1)

! calculate the Z bounds
!
      zmn = zmin + (k - 1) * zl
      zmx = zmin +  k      * zl

      do j = 1, rdims(2)

! claculate the block position along Y
!
        jl  = (j - 1) * res(1)

! calculate the Y bounds
!
        ymn = ymin + (j - 1) * yl
        ymx = ymin +  j      * yl

        do i = 1, rdims(1)

! claculate the block position along Y
!
          il  = (i - 1) * res(1)

! calculate the Z bounds
!
          xmn = xmin + (i - 1) * xl
          xmx = xmin +  i      * xl

! assign a pointer
!
          pmeta => block_array(i,j,k)%ptr

! mark it as the leaf
!
          call metablock_setleaf(pmeta)

! set the level
!
          call metablock_setlevel(pmeta, 1)

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
          call metablock_setbounds(pmeta, xmn, xmx, ymn, ymx, zmn, zmx)
        end do
      end do
    end do

!! ASSIGN THE BLOCK NEIGHBORS
!!
! assign boundaries along the X direction
!
    do k = 1, rdims(3)
      do j = 1, rdims(2)
        do i = 1, rdims(1) - 1

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
      do k = 1, rdims(3)
        do j = 1, rdims(2)

! assign pointers
!
          pmeta => block_array(      1 ,j,k)%ptr
          pnext => block_array(rdims(1),j,k)%ptr

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
    do k = 1, rdims(3)
      do j = 1, rdims(2) - 1
        do i = 1, rdims(1)

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
      do k = 1, rdims(3)
        do i = 1, rdims(1)

! assign pointers
!
          pmeta => block_array(i,      1 ,k)%ptr
          pnext => block_array(i,rdims(2),k)%ptr

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
    do k = 1, rdims(3) - 1
      do j = 1, rdims(2)
        do i = 1, rdims(1)

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
      do j = 1, rdims(2)
        do i = 1, rdims(1)

! assign pointers
!
          pmeta => block_array(i,j,      1 )%ptr
          pnext => block_array(i,j,rdims(3))%ptr

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
!
!-------------------------------------------------------------------------------
!
  end subroutine domain_default
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
                        , gamma, dens, pres, ratio, rcut
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt, idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: dx, dy, dz, r

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(km)     :: z
    real, dimension(nvr,im) :: q, u
!
!-------------------------------------------------------------------------------
!
! calculate the cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0d0
#endif /* NDIMS == 3 */

! generate the coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5d0) * dx + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5d0) * dy + pblock%meta%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5d0) * dz + pblock%meta%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0d0
#endif /* NDIMS == 3 */

! set the uniform variables
!
    q(idn,:) = dens
    q(ivx,:) = 0.0d0
    q(ivy,:) = 0.0d0
    q(ivz,:) = 0.0d0
#ifdef ADI
    q(ipr,:) = pres
#endif /* ADI */
#ifdef MHD
    q(ibx,:) = 1.0d0 / sqrt(2.0d0)
    q(iby,:) = 1.0d0 / sqrt(2.0d0)
    q(ibz,:) = 0.0d0
#ifdef GLM
    q(iph,:) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial star profile (density for ISO or pressure for ADI)
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im

          r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)

#ifdef ISO
          if (r .lt. rcut) then
            q(idn,i) = dens * ratio
          else
            q(idn,i) = dens
          end if
#endif /* ISO */
#ifdef ADI
          if (r .lt. rcut) then
            q(ipr,i) = pres * ratio
          else
            q(ipr,i) = pres
          end if
#endif /* ADI */
        end do

! convert the primitive variables to conserved ones
!
        call prim2cons(im, q(:,:), u(:,:))

! copy the conserved variables to the current block
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
!
!===============================================================================
!
! init_reconnection: subroutine initializes the variables for the reconnection
!                    problem
!
!===============================================================================
!
  subroutine init_reconnection(pblock)

    use blocks   , only : block_data
    use config   , only : in, jn, kn, im, jm, km, ng                           &
                        , xmin, xmax, dens, pres, bamp, bper, ydel, ycut
#ifdef ISO
    use config   , only : csnd2
#endif /* ISO */
    use constants, only : dpi
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt
    use variables, only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: xlen, dx, dy, dz, yexp
#ifdef MHD
    real                          :: ptot, pmag
#endif /* MHD */

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(km)     :: z
    real, dimension(nvr,im) :: q, u
!
!-------------------------------------------------------------------------------
!
! calculate the length of the box
!
    xlen = xmax - xmin

#ifdef MHD
! calculate the total pressure
!
    ptot = 0.5d0 * bamp * bamp
#endif /* MHD */

! calculate the cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate the coordinates
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
#ifdef ADI
    q(ipr,:) = pres
#endif /* ADI */
#ifdef MHD
    q(ibx,:) = 0.0d0
    q(iby,:) = 0.0d0
    q(ibz,:) = 0.0d0
#ifdef GLM
    q(iph,:) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial profiles
!
    do k = 1, km
      do j = 1, jm

! calculate the exponent factor
!
        yexp = exp(- 0.5d0 * (y(j) / ycut)**2)

#ifdef MHD
        do i = 1, im

! initialize the magnetic field components
!
          q(ibx,i) = bamp * tanh(y(j) / ydel)

! calculate magnetic pressure
!
          pmag = 0.5d0 * q(ibx,i) * q(ibx,i)

! add perturbation
!
          q(ibx,i) = q(ibx,i) &
                   - bper * yexp * y(j) * cos(dpi * x(i) / xlen) / ycut**2
          q(iby,i) = bper * yexp * dpi  * sin(dpi * x(i) / xlen)

! initialize density or pressure depending on EOS, so the total pressure is
! uniform
!

#ifdef ADI
          q(ipr,i) = pres + (ptot - pmag)
#endif /* ADI */
#ifdef ISO
          q(idn,i) = dens + (ptot - pmag) / csnd2
#endif /* ISO */
        end do
#endif /* MHD */

! convert primitive variables to conserved
!
        call prim2cons(im, q(:,:), u(:,:))

! copy conservative variables to the current block
!
        pblock%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_reconnection
!
!===============================================================================
!
! init_multi_current_sheet: subroutine initializes the set up for the multiple
!                           current sheet problem
!
!===============================================================================
!
  subroutine init_multi_current_sheet(pblock)

    use blocks   , only : block_data
    use config   , only : in, jn, kn, im, jm, km, ng                           &
                        , xmin, xmax, dens, pres, bamp, vper, ydel
#ifdef ISO
    use config   , only : csnd2
#endif /* ISO */
    use constants, only : dpi
    use mpitools , only : ncpu
    use random   , only : randomn
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt
    use variables, only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: xlen, dx, dy, dz, yi, yt
#ifdef MHD
    real                          :: ptot, pmag
#endif /* MHD */

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(km)     :: z
    real, dimension(nvr,im) :: q, u
!
!-------------------------------------------------------------------------------
!
! calculate the length of the box
!
    xlen = xmax - xmin

#ifdef MHD
! calculate the total pressure
!
    ptot = 0.5d0 * bamp * bamp
#endif /* MHD */

! calculate the cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn
#if NDIMS == 3
    dz = (pblock%meta%zmax - pblock%meta%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate the coordinates
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
#ifdef ADI
    q(ipr,:) = pres
#endif /* ADI */
#ifdef MHD
    q(ibx,:) = 0.0d0
    q(iby,:) = 0.0d0
    q(ibz,:) = 0.0d0
#ifdef GLM
    q(iph,:) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial profiles
!
    do k = 1, km
      do j = 1, jm

! calculate the exponent factor
!
        yi = y(j) - floor(y(j)) - 0.5d0
        yt = yi - sign(0.25d0, yi)

#ifdef MHD
        do i = 1, im

! initialize random velocity field
!
          q(ivx,i) = vper * randomn(ncpu)
          q(ivy,i) = vper * randomn(ncpu)
#if NDIMS == 3
          q(ivz,i) = vper * randomn(ncpu)
#endif /* NDIMS == 3 */

! initialize the magnetic field components
!
          q(ibx,i) = - sign(1.0d0, yi) * bamp * tanh(yt / ydel)

! calculate magnetic pressure
!
          pmag = 0.5d0 * q(ibx,i) * q(ibx,i)

! initialize density or pressure depending on EOS, so the total pressure is
! uniform
!

#ifdef ADI
          q(ipr,i) = pres + (ptot - pmag)
#endif /* ADI */
#ifdef ISO
          q(idn,i) = dens + (ptot - pmag) / csnd2
#endif /* ISO */
        end do
#endif /* MHD */

! convert primitive variables to conserved
!
        call prim2cons(im, q(:,:), u(:,:))

! copy conservative variables to the current block
!
        pblock%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_multi_current_sheet
!
!===============================================================================
!
! init_turbulence: subroutine initializes the variables for the turbulence
!                  problem
!
!===============================================================================
!
  subroutine init_turbulence(pblock)

    use blocks   , only : block_data
    use config   , only : im, jm, km, dens, pres, bamp
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt
    use variables, only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k

! local arrays
!
    real, dimension(nvr,im) :: q
    real, dimension(nqt,im) :: u
!
!-------------------------------------------------------------------------------
!
! set variables
!
    q(idn,:) = dens
    q(ivx,:) = 0.0d0
    q(ivy,:) = 0.0d0
    q(ivz,:) = 0.0d0
#ifdef ADI
    q(ipr,:) = pres
#endif /* ADI */
#ifdef MHD
    q(ibx,:) = bamp
    q(iby,:) = 0.0d0
    q(ibz,:) = 0.0d0
#ifdef GLM
    q(iph,:) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial profiles
!
    do k = 1, km
      do j = 1, jm

! convert primitive variables to conserved
!
        call prim2cons(im, q(1:nvr,1:im), u(1:nqt,1:im))

! copy conservative variables to the current block
!
        pblock%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_turbulence
!
!===============================================================================
!
! init_orszag_tang: subroutine initializes the setup for Orszag-Tang problems
!
!===============================================================================
!
  subroutine init_orszag_tang(pblock)

    use constants, only : dpi, qpi
    use blocks   , only : block_data
    use config   , only : im, jm, km, in, jn, kn, ng
    use config   , only : dens, bamp
#ifdef ADI
    use config   , only : gamma, pres
#endif /* ADI */
#ifdef ISO
    use config   , only : csnd
#endif /* ISO */
    use scheme   , only : prim2cons
    use variables, only : nvr, nqt
    use variables, only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#ifdef GLM
    use variables, only : iph
#endif /* GLM */
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: dx, dy

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(nvr,im) :: q
    real, dimension(nqt,im) :: u
!
!-------------------------------------------------------------------------------
!
! calculate constants
!
#ifdef ADI
    bamp = 1.0d0 / sqrt(qpi)
    pres = gamma / qpi
    dens = gamma**2 / qpi
#endif /* ADI */
#ifdef ISO
    bamp = csnd * (5.0d0 / 3.0d0) / qpi
#endif /* ISO */

! calculate the cell sizes
!
    dx = (pblock%meta%xmax - pblock%meta%xmin) / in
    dy = (pblock%meta%ymax - pblock%meta%ymin) / jn

! generate the coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5d0) * dx + pblock%meta%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5d0) * dy + pblock%meta%ymin

! set variables
!
    q(idn,:) = dens
#ifdef ADI
    q(ipr,:) = pres
#endif /* ADI */
#ifdef MHD
#ifdef GLM
    q(iph,:) = 0.0d0
#endif /* GLM */
#endif /* MHD */

! set the initial profiles
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          q(ivx,i) = - sin(dpi * y(j))
          q(ivy,i) =   sin(dpi * x(i))
          q(ivz,i) = 0.0d0
#ifdef MHD
          q(ibx,i) = - bamp * sin(dpi * y(j))
          q(iby,i) =   bamp * sin(qpi * x(i))
          q(ibz,i) = 0.0d0
#endif /* MHD */
        end do

! convert primitive variables to conserved
!
        call prim2cons(im, q(1:nvr,1:im), u(1:nqt,1:im))

! copy conservative variables to the current block
!
        pblock%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine init_orszag_tang
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
    use config   , only : im, jm, km, ib, ie, jb, je, kb, ke                   &
                        , crefmin, crefmax, epsref
    use scheme   , only : cons2prim
    use variables, only : nvr, nqt
    use variables, only : idn
#ifdef ADI
    use variables, only : ipr
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#endif /* MHD */

! input arguments
!
    type(block_data), pointer, intent(inout) :: pblock

! return variable
!
    integer(kind=1) :: check_ref

! local variables
!
    integer :: i, j, k, im2, jm2, km2, ip2, jp2, kp2
    real    :: cref, fl, fr, fc, fx, fy, fz

! local arrays
!
    real, dimension(nvr,im)   :: u, q
    real, dimension(im,jm,km) :: dn
#ifdef ADI
    real, dimension(im,jm,km) :: pr
#endif /* ADI */
#ifdef MHD
    real, dimension(im,jm,km) :: bx, by, bz
#endif /* MHD */
    real, parameter           :: cf  = 1.0d0 / NDIMS
!
!-------------------------------------------------------------------------------
!
! convert conserved variables to primitive ones
!
    do k = 1, km
      do j = 1, jm
        dn(1:im,j,k) = pblock%u(idn,1:im,j,k)
#ifdef ADI
        u(1:nqt,1:im) = pblock%u(1:nqt,1:im,j,k)

        call cons2prim(im, u(:,:), q(:,:))

        pr(1:im,j,k) = q(ipr,1:im)
#endif /* ADI */
#ifdef MHD
        bx(1:im,j,k) = pblock%u(ibx,1:im,j,k)
        by(1:im,j,k) = pblock%u(iby,1:im,j,k)
        bz(1:im,j,k) = pblock%u(ibz,1:im,j,k)
#endif /* MHD */
      end do
    end do

! reset indicators
!
    cref = 0.0d0

! find the maximum smoothness indicator over all cells
!
    do k = kb, ke
#if NDIMS == 3
      km2 = k - 2
      kp2 = k + 2
#endif /* NDIMS == 3 */
      do j = jb, je
        jm2 = j - 2
        jp2 = j + 2
        do i = ib, ie
          im2 = i - 2
          ip2 = i + 2

! density
!
          fr   = dn(ip2,j,k) - dn(i,j,k)
          fl   = dn(im2,j,k) - dn(i,j,k)
          fc   = dn(ip2,j,k) + dn(im2,j,k) + 2.0d0 * dn(i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = dn(i,jp2,k) - dn(i,j,k)
          fl   = dn(i,jm2,k) - dn(i,j,k)
          fc   = dn(i,jp2,k) + dn(i,jm2,k) + 2.0d0 * dn(i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = dn(i,j,kp2) - dn(i,j,k)
          fl   = dn(i,j,km2) - dn(i,j,k)
          fc   = dn(i,j,kp2) + dn(i,j,km2) + 2.0d0 * dn(i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#ifdef ADI
! pressure
!
          fr   = pr(ip2,j,k) - pr(i,j,k)
          fl   = pr(im2,j,k) - pr(i,j,k)
          fc   = pr(ip2,j,k) + pr(im2,j,k) + 2.0d0 * pr(i,j,k)
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          fr   = pr(i,jp2,k) - pr(i,j,k)
          fl   = pr(i,jm2,k) - pr(i,j,k)
          fc   = pr(i,jp2,k) + pr(i,jm2,k) + 2.0d0 * pr(i,j,k)
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = pr(i,j,kp2) - pr(i,j,k)
          fl   = pr(i,j,km2) - pr(i,j,k)
          fc   = pr(i,j,kp2) + pr(i,j,km2) + 2.0d0 * pr(i,j,k)
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * fc)
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#endif /* ADI */
#ifdef MHD
! X magnetic component
!
          fr   = bx(ip2,j,k) - bx(i,j,k)
          fl   = bx(im2,j,k) - bx(i,j,k)
          fc   = abs(bx(ip2,j,k)) + abs(bx(im2,j,k)) + 2.0d0 * abs(bx(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = bx(i,jp2,k) - bx(i,j,k)
          fl   = bx(i,jm2,k) - bx(i,j,k)
          fc   = abs(bx(i,jp2,k)) + abs(bx(i,jm2,k)) + 2.0d0 * abs(bx(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = bx(i,j,kp2) - bx(i,j,k)
          fl   = bx(i,j,km2) - bx(i,j,k)
          fc   = abs(bx(i,j,kp2)) + abs(bx(i,j,km2)) + 2.0d0 * abs(bx(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Y magnetic component
!
          fr   = by(ip2,j,k) - by(i,j,k)
          fl   = by(im2,j,k) - by(i,j,k)
          fc   = abs(by(ip2,j,k)) + abs(by(im2,j,k)) + 2.0d0 * abs(by(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = by(i,jp2,k) - by(i,j,k)
          fl   = by(i,jm2,k) - by(i,j,k)
          fc   = abs(by(i,jp2,k)) + abs(by(i,jm2,k)) + 2.0d0 * abs(by(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = by(i,j,kp2) - by(i,j,k)
          fl   = by(i,j,km2) - by(i,j,k)
          fc   = abs(by(i,j,kp2)) + abs(by(i,j,km2)) + 2.0d0 * abs(by(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

! Z magnetic component
!
          fr   = bz(ip2,j,k) - bz(i,j,k)
          fl   = bz(im2,j,k) - bz(i,j,k)
          fc   = abs(bz(ip2,j,k)) + abs(bz(im2,j,k)) + 2.0d0 * abs(bz(i,j,k))
          fx   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          fr   = bz(i,jp2,k) - bz(i,j,k)
          fl   = bz(i,jm2,k) - bz(i,j,k)
          fc   = abs(bz(i,jp2,k)) + abs(bz(i,jm2,k)) + 2.0d0 * abs(bz(i,j,k))
          fy   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
#if NDIMS == 2
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
          fr   = bz(i,j,kp2) - bz(i,j,k)
          fl   = bz(i,j,km2) - bz(i,j,k)
          fc   = abs(bz(i,j,kp2)) + abs(bz(i,j,km2)) + 2.0d0 * abs(bz(i,j,k))
          fz   = abs(fr + fl) / (abs(fr) + abs(fl) + epsref * (fc + 1.0e-8))
          cref = max(cref, sqrt(cf * (fx * fx + fy * fy + fz * fz)))
#endif /* NDIMS == 3 */

#endif /* MHD */
        end do
      end do
    end do

! return the refinement flag depending on the condition value
!
    check_ref = 0

    if (cref .ge. crefmax) then
      check_ref =  1
    end if
    if (cref .le. crefmin) then
      check_ref = -1
    end if

    return

!-------------------------------------------------------------------------------
!
  end function check_ref

!===============================================================================
!
end module
