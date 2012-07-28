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
!! module: PROBLEMS
!!
!!  This module handles the initialization of various test and research
!!  problems.
!!
!!
!!******************************************************************************
!
module problems

! include external procedures and variables
!
  use parameters, only : get_parameter_integer, get_parameter_real             &
                       , get_parameter_string

! module variables are not implicit by default
!
  implicit none

! module variable to store the problem name
!
  character(len=32), save :: problem = "blast"

! refinement criterion parameters
!
  real             , save :: crefmin = 2.0d-1
  real             , save :: crefmax = 5.0d-1
  real             , save :: epsref  = 1.0d-2

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_problems, setup_domain, setup_problem
  public :: check_refinement_criterion

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
! subroutine INITIALIZE_PROBLEMS:
! ------------------------------
!
!   Subroutine prepares module PROBLEMS.
!
!
!===============================================================================
!
  subroutine initialize_problems()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! get the problem name
!
    call get_parameter_string("problem", problem)

! get the refinement parameters
!
    call get_parameter_real  ("crefmin", crefmin)
    call get_parameter_real  ("crefmax", crefmax)
    call get_parameter_real  ("epsref" , epsref )

!-------------------------------------------------------------------------------
!
  end subroutine initialize_problems
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
!
! subroutine SETUP_PROBLEM:
! ------------------------
!
!   Subroutine sets the initial conditions for selected problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!
!===============================================================================
!
  subroutine setup_problem(pdata)

! include external procedures and variables
!
    use blocks, only : block_data
    use error , only : print_error

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata
!
!-------------------------------------------------------------------------------
!
! select the setup subroutine depending on the problem name
!
    select case(problem)

! general test problems
!
    case("blast")
      call setup_problem_blast(pdata)

    case default
      call print_error("problems::init_problem()"                              &
                     , "Setup subroutime is not implemented for this problem!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem
!
!===============================================================================
!
! function CHECK_REFINEMENT_CRITERION:
! -----------------------------------
!
!   Function scans the given data block and checks for the refinement
!   criterion.  It returns +1 if the criterion is met, which indicates that
!   the block needs to be refined, 0 if there is no need for the refinement, and
!   -1 if the block can be derefined.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  function check_refinement_criterion(pdata) result(criterion)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km, ib, ie, jb, je, kb, ke
    use scheme     , only : cons2prim
    use variables  , only : nvr, nqt
    use variables  , only : idn
#ifdef ADI
    use variables  , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables  , only : ibx, iby, ibz
#endif /* MHD */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! return variable
!
    integer(kind=4) :: criterion

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
        dn(1:im,j,k) = pdata%u(idn,1:im,j,k)
#ifdef ADI
        u(1:nqt,1:im) = pdata%u(1:nqt,1:im,j,k)

        call cons2prim(im, u(:,:), q(:,:))

        pr(1:im,j,k) = q(ipr,1:im)
#endif /* ADI */
#ifdef MHD
        bx(1:im,j,k) = pdata%u(ibx,1:im,j,k)
        by(1:im,j,k) = pdata%u(iby,1:im,j,k)
        bz(1:im,j,k) = pdata%u(ibz,1:im,j,k)
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
    criterion = 0

    if (cref .ge. crefmax) then
      criterion =  1
    end if
    if (cref .le. crefmin) then
      criterion = -1
    end if

    return

!-------------------------------------------------------------------------------
!
  end function check_refinement_criterion
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
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_BLAST:
! ------------------------------
!
!   Subroutine sets the initial conditions for the spherical blast problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!
!===============================================================================
!
  subroutine setup_problem_blast(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az
    use equations  , only : gamma
    use scheme     , only : prim2cons
    use variables  , only : nvr, nqt, idn, ivx, ivy, ivz
#ifdef ADI
    use variables  , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables  , only : ibx, iby, ibz
#ifdef GLM
    use variables  , only : iph
#endif /* GLM */
#endif /* MHD */

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! default parameter values
!
    real   , save :: dens = 1.0d0, ratio = 1.0e2, radius = 0.1d0
#ifdef ADI
    real   , save :: csnd = 0.40824829046386301635d0
    real   , save :: pres = 1.0d0
#endif /* ADI */
    logical, save :: first = .true.

! local variables
!
    integer(kind=4), dimension(3) :: dm
    integer                       :: i, j, k
    real                          :: r

! local arrays
!
    real, dimension(im)     :: x
    real, dimension(jm)     :: y
    real, dimension(km)     :: z
    real, dimension(nvr,im) :: q, u
!
!-------------------------------------------------------------------------------
!
! prepare problem constants during the first subroutine call
!
    if (first) then

! get problem parameters
!
      call get_parameter_real("dens"  , dens  )
      call get_parameter_real("ratio" , ratio )
      call get_parameter_real("radius", radius)
#ifdef ADI
      call get_parameter_real("csnd"  , csnd  )

! calculate pressure
!
      pres = dens * csnd * csnd / gamma
#endif /* ADI */

! reset the first execution flag
!
      first = .false.

    end if

! obtain block coordinates
!
    x(:) = pdata%meta%xmin + ax(pdata%meta%level,:)
    y(:) = pdata%meta%ymin + ay(pdata%meta%level,:)
#if NDIMS == 3
    z(:) = pdata%meta%zmin + az(pdata%meta%level,:)
#else /* NDIMS == 3 */
    z(:) = 0.0d0
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
          if (r .lt. radius) then
            q(idn,i) = dens * ratio
          else
            q(idn,i) = dens
          end if
#endif /* ISO */
#ifdef ADI
          if (r .lt. radius) then
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
        pdata%u(1:nqt,1:im,j,k) = u(1:nqt,1:im)

      end do
    end do

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_blast

!===============================================================================
!
end module problems
