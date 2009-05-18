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

    use blocks, only : block
    use config, only : problem

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

! local arguments
!
    type(block), pointer :: pb
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

    use blocks, only : block
    use config, only : problem

! input arguments
!
    type(block), pointer    , intent(inout) :: pblock
    real, dimension(:,:,:,:), intent(inout) :: du

! local arguments
!
    type(block), pointer :: pb
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

    use blocks, only : block, append_block
    use config, only : xlbndry, xubndry, ylbndry, yubndry                      &
                     , xmin, xmax, ymin, ymax

    implicit none

! local variables
!
    real :: xl, xc, xr, yl, yc, yr

! local pointers
!
    type(block), pointer :: pbl, pbr, ptl, ptr
!
!-------------------------------------------------------------------------------
!
! create root blocks
!
    call append_block(pbl)
    call append_block(ptl)
    call append_block(ptr)
    call append_block(pbr)

! set configurations
!
    pbl%config = 'D'
    ptl%config = 'N'
    ptr%config = 'N'
    pbr%config = 'C'

! set positions
!
    pbl%pos(1) = 1
    pbl%pos(2) = 1
    pbr%pos(1) = 2
    pbr%pos(2) = 1
    ptl%pos(1) = 1
    ptl%pos(2) = 2
    ptr%pos(1) = 2
    ptr%pos(2) = 2

! set leaf flags
!
    pbl%leaf   = .true.
    pbr%leaf   = .true.
    ptl%leaf   = .true.
    ptr%leaf   = .true.

! set levels
!
    pbl%level  = 1
    pbr%level  = 1
    ptl%level  = 1
    ptr%level  = 1

! set neighbors
!
    if (xlbndry .eq. 'periodic') &
      pbl%neigh(1,1,:)%cpu = pbr%cpu
    pbl%neigh(1,2,:)%cpu = pbr%cpu
    if (ylbndry .eq. 'periodic') &
      pbl%neigh(2,1,:)%cpu = ptl%cpu
    pbl%neigh(2,2,:)%cpu = ptl%cpu

    pbr%neigh(1,1,:)%cpu = pbl%cpu
    if (xubndry .eq. 'periodic') &
      pbr%neigh(1,2,:)%cpu = pbl%cpu
    if (ylbndry .eq. 'periodic') &
      pbr%neigh(2,1,:)%cpu = ptr%cpu
    pbr%neigh(2,2,:)%cpu = ptr%cpu

    if (xlbndry .eq. 'periodic') &
      ptl%neigh(1,1,:)%cpu = ptr%cpu
    ptl%neigh(1,2,:)%cpu = ptr%cpu
    ptl%neigh(2,1,:)%cpu = pbl%cpu
    if (yubndry .eq. 'periodic') &
      ptl%neigh(2,2,:)%cpu = pbl%cpu

    ptr%neigh(1,1,:)%cpu = ptl%cpu
    if (xubndry .eq. 'periodic') &
      ptr%neigh(1,2,:)%cpu = ptl%cpu
    ptr%neigh(2,1,:)%cpu = pbr%cpu
    if (yubndry .eq. 'periodic') &
      ptr%neigh(2,2,:)%cpu = pbr%cpu

    if (xlbndry .eq. 'periodic') &
      pbl%neigh(1,1,:)%id = pbr%id
    pbl%neigh(1,2,:)%id = pbr%id
    if (ylbndry .eq. 'periodic') &
      pbl%neigh(2,1,:)%id = ptl%id
    pbl%neigh(2,2,:)%id = ptl%id

    pbr%neigh(1,1,:)%id = pbl%id
    if (xubndry .eq. 'periodic') &
      pbr%neigh(1,2,:)%id = pbl%id
    if (ylbndry .eq. 'periodic') &
      pbr%neigh(2,1,:)%id = ptr%id
    pbr%neigh(2,2,:)%id = ptr%id

    if (xlbndry .eq. 'periodic') &
      ptl%neigh(1,1,:)%id = ptr%id
    ptl%neigh(1,2,:)%id = ptr%id
    ptl%neigh(2,1,:)%id = pbl%id
    if (yubndry .eq. 'periodic') &
      ptl%neigh(2,2,:)%id = pbl%id

    ptr%neigh(1,1,:)%id = ptl%id
    if (xubndry .eq. 'periodic') &
      ptr%neigh(1,2,:)%id = ptl%id
    ptr%neigh(2,1,:)%id = pbr%id
    if (yubndry .eq. 'periodic') &
      ptr%neigh(2,2,:)%id = pbr%id

! set the bounds of the blocks
!
    xl = xmin
    xc = 0.5 * (xmax + xmin)
    xr = xmax
    yl = ymin
    yc = 0.5 * (ymax + ymin)
    yr = ymax

    pbl%xmin = xl
    pbl%xmax = xc
    pbl%ymin = yl
    pbl%ymax = yc

    ptl%xmin = xl
    ptl%xmax = xc
    ptl%ymin = yc
    ptl%ymax = yr

    ptr%xmin = xc
    ptr%xmax = xr
    ptr%ymin = yc
    ptr%ymax = yr

    pbr%xmin = xc
    pbr%xmax = xr
    pbr%ymin = yl
    pbr%ymax = yc

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

    use blocks, only : block, idn, imx, imy, imz, ien
    use config, only : in, jn, kn, im, jm, km, ng, dens, pres, gammam1i

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
!-------------------------------------------------------------------------------
!
! calculate parameters
!
    enamb = gammam1i * pres
    en    = 100.0 * enamb

! allocate coordinates
!
    allocate(x(im))
    allocate(y(jm))
    allocate(z(km))

! calculate cell sizes
!
    dx = (pblock%xmax - pblock%xmin) / in
    dy = (pblock%ymax - pblock%ymin) / jn
#if NDIMS == 3
    dz = (pblock%zmax - pblock%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%zmin
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
    do k = 1, km
      do j = 1, jm
        do i = 1, im

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

    use blocks, only : block, idn, imx, imy, imz, ien
    use config, only : in, jn, kn, im, jm, km, ng, dens, pres, rmid, gammam1i

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

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
    dx  = (pblock%xmax - pblock%xmin) / in
    dy  = (pblock%ymax - pblock%ymin) / jn
    dxh = 0.5d0 * dx
    dyh = 0.5d0 * dy
    ds  = dx * dy

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng) * dx - dxh + pblock%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng) * dy - dyh + pblock%ymin

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

    use blocks, only : block, idn, imx, imy, imz, ien
    use config, only : ng, in, jn, kn, im, jm, km, dens, pres, dnfac, dnrat    &
                     , x1c, y1c, z1c, r1c, x2c, y2c, z2c, r2c                  &
                     , csnd2, gamma, gammam1i

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

! local variables
!
    integer :: i, j, k
    real    :: dx, dy, dz, dnamb, enamb
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
    dx = (pblock%xmax - pblock%xmin) / in
    dy = (pblock%ymax - pblock%ymin) / jn
#if NDIMS == 3
    dz = (pblock%zmax - pblock%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0
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
      z1l = z(k) - z1c
      z2l = z(k) - z2c

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
#ifdef ADI
            pblock%u(ien,i,j,k) = enstar1
#endif /* ADI */
          endif

          if (r2 .le. r2c) then
            pblock%u(idn,i,j,k) = dnstar2
#ifdef ADI
            pblock%u(ien,i,j,k) = enstar2
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

    use blocks, only : block, idn, imx, imy, imz, ien
    use config, only : ng, in, jn, kn, im, jm, km, dens, pres, dnfac, dnrat    &
                     , x1c, y1c, z1c, r1c, x2c, y2c, z2c, r2c                  &
                     , csnd2, gamma, gammam1i

! input arguments
!
    type(block), pointer, intent(inout) :: pblock
    real, dimension(:,:,:,:), intent(inout) :: du

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
    dx = (pblock%xmax - pblock%xmin) / in
    dy = (pblock%ymax - pblock%ymin) / jn
#if NDIMS == 3
    dz = (pblock%zmax - pblock%zmin) / kn
#else /* NDIMS == 3 */
    dz = 1.0
#endif /* NDIMS == 3 */

! generate coordinates
!
    x(:) = ((/(i, i = 1, im)/) - ng - 0.5) * dx + pblock%xmin
    y(:) = ((/(j, j = 1, jm)/) - ng - 0.5) * dy + pblock%ymin
#if NDIMS == 3
    z(:) = ((/(k, k = 1, km)/) - ng - 0.5) * dz + pblock%zmin
#else /* NDIMS == 3 */
    z(1) = 0.0
#endif /* NDIMS == 3 */

! reset update
!
    do k = 1, km
      z1l = z(k) - z1c
      z2l = z(k) - z2c

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

    use blocks, only : block, idn, imx, imy, imz, ien, nv => nvars
    use config, only : im, jm, km, ib, ie, jb, je, kb, ke, gammam1i            &
                     , crefmin, crefmax

! input arguments
!
    type(block), pointer, intent(inout) :: pblock

! return variable
!
    integer(kind=1) :: check_ref

! local variables
!
    integer                       :: i, j, k
    real                          :: dpmax, vx, vy, vz, en, ek, ei
    real                          :: dnl, dnr, prl, prr, ddndx, ddndy, ddn     &
                                   , dprdx, dprdy, dpr
    real, dimension(im,jm,km)     :: dn, pr
!
!-------------------------------------------------------------------------------
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          dn(i,j,k) = pblock%u(idn,i,j,k)
          vx = pblock%u(imx,i,j,k) / dn(i,j,k)
          vy = pblock%u(imy,i,j,k) / dn(i,j,k)
          vz = pblock%u(imz,i,j,k) / dn(i,j,k)
          en = pblock%u(ien,i,j,k)

          ek = 0.5 * dn(i,j,k) * (vx*vx + vy*vy + vz*vz)
          ei = en - ek
          pr(i,j,k) = gammam1i * ei
        end do
      end do
    end do

! check gradient of pressure
!
    dpmax = 0.0d0

    do k = kb, ke
      do j = jb-2, je+2
        do i = ib-2, ie+2
          dnl = dn(i-1,j,k)
          dnr = dn(i+1,j,k)
          ddndx = abs(dnr-dnl)/(dnr+dnl)
          dnl = dn(i,j-1,k)
          dnr = dn(i,j+1,k)
          ddndy = abs(dnr-dnl)/(dnr+dnl)

          ddn = sqrt(ddndx**2 + ddndy**2)

          prl = pr(i-1,j,k)
          prr = pr(i+1,j,k)
          dprdx = abs(prr-prl)/(prr+prl)
          prl = pr(i,j-1,k)
          prr = pr(i,j+1,k)
          dprdy = abs(prr-prl)/(prr+prl)

          dpr = sqrt(dprdx**2 + dprdy**2)

          pblock%c(i,j,k) = max(ddn,dpr)

          dpmax = max(dpmax, ddn, dpr)

        end do
      end do
    end do

! check condition
!
    check_ref = 0

    if (dpmax .ge. crefmax) then
      check_ref =  1
    endif
    if (dpmax .le. crefmin) then
      check_ref = -1
    endif

    return

!-------------------------------------------------------------------------------
!
  end function check_ref

!===============================================================================
!
end module
