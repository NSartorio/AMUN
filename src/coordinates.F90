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
!! module: COORDINATES
!!
!!  This module provides variables and subroutines handling the coordinates
!!  for all refinement levels.
!!
!!******************************************************************************
!
module coordinates

! module variables are not implicit by default
!
  implicit none

! the domain block dimensions
!
  integer, save :: nn =  8, in =  8, jn =  8, kn =  1

! the number of ghost zones
!
  integer, save :: ng =  2

! the domain block dimensions including the ghost zones
!
  integer, save :: im = 12, jm = 12, km =  1

! the domain division
!
  integer, save :: ir =  1, jr =  1, kr =  1

! the limits of refinement level
!
  integer, save :: minlev = 1, maxlev = 1, toplev = 1

! the domain bounds
!
  real, save :: xmin = 0.0d0
  real, save :: xmax = 1.0d0
  real, save :: xlen = 1.0d0
  real, save :: ymin = 0.0d0
  real, save :: ymax = 1.0d0
  real, save :: ylen = 1.0d0
  real, save :: zmin = 0.0d0
  real, save :: zmax = 1.0d0
  real, save :: zlen = 1.0d0

! the effective resolution of the full domain
!
  integer,dimension(NDIMS), save :: effres  = 1

! the block coordinates for all levels of refinement
!
  real, dimension(:,:), allocatable, save :: ax  , ay  , az
  real, dimension(:  ), allocatable, save :: adx , ady , adz, adr
  real, dimension(:  ), allocatable, save :: adxi, adyi, adzi
  real, dimension(:  ), allocatable, save :: advol

! the block resolution in the units of effective resolution for all levels
!
  integer(kind=4), dimension(:,:), allocatable, save :: res

! by default everything is private
!
  public

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_COORDINATES:
! ---------------------------------
!
!   Subroutine initializes mesh coordinates and other coordinate parameters.
!
!===============================================================================
!
  subroutine initialize_coordinates(flag)

! include external procedures and variables
!
    use parameters, only : get_parameter_integer, get_parameter_real

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    logical, intent(in) :: flag

! local variables
!
    integer :: i, j, k, l
    integer :: ni, nj, nk
    logical :: info

! local arrays
!
    integer(kind=4), dimension(3) :: cm, rm, dm
!
!-------------------------------------------------------------------------------
!
! obtain the number of cells along each block dimension in the case when all
! dimensions are the same
!
    call get_parameter_integer("ncells", nn    )

! set the block dimensions
!
    in = nn
    jn = nn
#if NDIMS == 3
    kn = nn
#endif /* NDIMS == 3 */

! change individual block dimension if set
!
    call get_parameter_integer("icells", in    )
    call get_parameter_integer("jcells", jn    )
#if NDIMS == 3
    call get_parameter_integer("kcells", kn    )
#endif /* NDIMS == 3 */

! obtain the number of ghost zones
!
    call get_parameter_integer("nghost", ng    )

! calculate the block dimensions including ghost cells
!
    im = in + 2 * ng
    jm = jn + 2 * ng
#if NDIMS == 3
    km = kn + 2 * ng
#endif /* NDIMS == 3 */

! obtain the refinement level bounds
!
    call get_parameter_integer("minlev", minlev)
    call get_parameter_integer("maxlev", maxlev)

! set the top level
!
    toplev = maxlev

! obtain the domain base division
!
    call get_parameter_integer("rdimx" , ir    )
    call get_parameter_integer("rdimy" , jr    )
#if NDIMS == 3
    call get_parameter_integer("rdimz" , kr    )
#endif /* NDIMS == 3 */

! obtain the domain bounds
!
    call get_parameter_real   ("xmin"  , xmin  )
    call get_parameter_real   ("xmax"  , xmax  )
    call get_parameter_real   ("ymin"  , ymin  )
    call get_parameter_real   ("ymax"  , ymax  )
    call get_parameter_real   ("zmin"  , zmin  )
    call get_parameter_real   ("zmax"  , zmax  )

! allocate space for coordinate variables and resolutions
!
    allocate(ax   (toplev, im))
    allocate(ay   (toplev, jm))
    allocate(az   (toplev, km))
    allocate(adx  (toplev))
    allocate(ady  (toplev))
    allocate(adz  (toplev))
    allocate(adr  (toplev))
    allocate(adxi (toplev))
    allocate(adyi (toplev))
    allocate(adzi (toplev))
    allocate(advol(toplev))
    allocate(res  (toplev, 3))

! reset all coordinate variables to initial values
!
    ax   (:,:) = 0.0d0
    ay   (:,:) = 0.0d0
    az   (:,:) = 0.0d0
    adx  (:)   = 1.0d0
    ady  (:)   = 1.0d0
    adz  (:)   = 1.0d0
    adr  (:)   = 1.0d0
    adxi (:)   = 1.0d0
    adyi (:)   = 1.0d0
    adzi (:)   = 1.0d0
    advol(:)   = 1.0d0
    res  (:,:) = 1

! generate the coordinate variables for each level
!
    do l = 1, toplev

! calculate the block resolution at each level
!
      j  = 2**(l - 1)
      ni = in * j
      nj = jn * j
      nk = kn * j

! calculate the cell sizes for each level
!
      adx (l) = (xmax - xmin) / (ir * ni)
      ady (l) = (ymax - ymin) / (jr * nj)
#if NDIMS == 3
      adz (l) = (zmax - zmin) / (kr * nk)
#endif /* NDIMS == 3 */
#if NDIMS == 2
      adr (l) = sqrt(adx(l)**2 + ady(l)**2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      adr (l) = sqrt(adx(l)**2 + ady(l)**2 + adz(l)**2)
#endif /* NDIMS == 3 */

! calculate the inverse of cell size
!
      adxi(l) = 1.0d0 / adx(l)
      adyi(l) = 1.0d0 / ady(l)
#if NDIMS == 3
      adzi(l) = 1.0d0 / adz(l)
#endif /* NDIMS == 3 */

! calculate the block coordinates for each level
!
      ax(l,:) = ((/(i, i = 1, im)/) - ng - 0.5d0) * adx(l)
      ay(l,:) = ((/(j, j = 1, jm)/) - ng - 0.5d0) * ady(l)
#if NDIMS == 3
      az(l,:) = ((/(k, k = 1, km)/) - ng - 0.5d0) * adz(l)
#endif /* NDIMS == 3 */

! calculate the cell volume at each level
!
      advol(l) = adx(l) * ady(l) * adz(l)

! calculate the effective block resolution at each level
!
      j         = 2**(maxlev - l)
      res(l,1)  = max(1, in * j)
      res(l,2)  = max(1, jn * j)
#if NDIMS == 3
      res(l,3)  = max(1, kn * j)
#endif /* NDIMS == 3 */

    end do

! print general information about the level resolutions
!
    if (flag) then

! the base resolution
!
      cm(1) = ir * in
      cm(2) = jr * jn
      cm(3) = kr * kn

! the effective resolution
!
      rm(1) = ir * res(1,1)
      rm(2) = jr * res(1,2)
      rm(3) = kr * res(1,3)

! the effective resolution
!
      effres(1:NDIMS) = rm(1:NDIMS)

! the top level block division
!
      dm(1) = rm(1) / in
      dm(2) = rm(2) / jn
      dm(3) = rm(3) / kn

! print info
!
      write(*,*)
      write(*,"(1x,a)"         ) "Geometry:"
      write(*,"(4x,a,  1x,i6)" ) "refinement to level    =", toplev
      write(*,"(4x,a,3(1x,i6))") "base configuration     =", ir, jr, kr
      write(*,"(4x,a,3(1x,i6))") "top level blocks       =", dm(:)
      write(*,"(4x,a,  1x,i6)" ) "maxium cover blocks    =", product(dm(:))
      write(*,"(4x,a,3(1x,i6))") "base resolution        =", cm(1:NDIMS)
      write(*,"(4x,a,3(1x,i6))") "effective resolution   =", rm(1:NDIMS)

    end if ! master

!-------------------------------------------------------------------------------
!
  end subroutine initialize_coordinates
!
!===============================================================================
!
! subroutine FINALIZE_COORDINATES:
! -------------------------------
!
!   Subroutine deallocates mesh coordinates.
!
!===============================================================================
!
  subroutine finalize_coordinates()

! local variables are not implicit by default
!
    implicit none

!-------------------------------------------------------------------------------
!
! deallocating coordinate variables
!
    if (allocated(ax)   ) deallocate(ax)
    if (allocated(ay)   ) deallocate(ay)
    if (allocated(az)   ) deallocate(az)
    if (allocated(adx)  ) deallocate(adx)
    if (allocated(ady)  ) deallocate(ady)
    if (allocated(adz)  ) deallocate(adz)
    if (allocated(adr)  ) deallocate(adr)
    if (allocated(adxi) ) deallocate(adxi)
    if (allocated(adyi) ) deallocate(adyi)
    if (allocated(adzi) ) deallocate(adzi)
    if (allocated(advol)) deallocate(advol)
    if (allocated(res)  ) deallocate(res)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_coordinates

!===============================================================================
!
end module
