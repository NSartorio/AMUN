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

! MODULE PARAMETERS:
! =================
!
! the domain block dimensions
!
  integer     , save :: nc =  8, in =  8, jn =  8, kn =  1

! the number of ghost zones
!
  integer     , save :: ng =  2, nh =  1, nd =  4

! the domain block dimensions including the ghost zones
!
  integer     , save :: im = 12, jm = 12, km =  1

! the domain block dimensions including the ghost zones
!
  integer     , save :: it = 11, jt = 11, kt =  1

! the domain division
!
  integer     , save :: ir =  1, jr =  1, kr =  1

! the limits of refinement level
!
  integer     , save :: minlev = 1, maxlev = 1, toplev = 1

! block indices
!
  integer     , save :: ih  =  6, jh  =  6, kh  =  1
  integer     , save :: ib  =  3, jb  =  3, kb  =  1
  integer     , save :: ie  = 10, je  = 10, ke  =  1
  integer     , save :: ibl =  2, jbl =  2, kbl =  1
  integer     , save :: ibu =  4, jbu =  4, kbu =  1
  integer     , save :: iel =  9, jel =  9, kel =  1
  integer     , save :: ieu = 11, jeu = 11, keu =  1

! the domain bounds
!
  real(kind=8), save :: xmin = 0.0d+00
  real(kind=8), save :: xmax = 1.0d+00
  real(kind=8), save :: xlen = 1.0d+00
  real(kind=8), save :: ymin = 0.0d+00
  real(kind=8), save :: ymax = 1.0d+00
  real(kind=8), save :: ylen = 1.0d+00
  real(kind=8), save :: zmin = 0.0d+00
  real(kind=8), save :: zmax = 1.0d+00
  real(kind=8), save :: zlen = 1.0d+00

! the domain volume and its inversion
!
  real(kind=8), save ::  vol = 1.0d+00
  real(kind=8), save :: voli = 1.0d+00

! the domain boundary areas
!
  real(kind=8), save :: xarea = 1.0d+00
  real(kind=8), save :: yarea = 1.0d+00
  real(kind=8), save :: zarea = 1.0d+00

! the block coordinates for all levels of refinement
!
  real(kind=8), dimension(:,:), allocatable, save :: ax  , ay  , az
  real(kind=8), dimension(:  ), allocatable, save :: adx , ady , adz, adr
  real(kind=8), dimension(:  ), allocatable, save :: adxi, adyi, adzi
  real(kind=8), dimension(:  ), allocatable, save :: advol

! define type for rectangular subarray description
!
  type rectangular
    integer, dimension(NDIMS) :: l ! indices of the lower corner
    integer, dimension(NDIMS) :: u ! indices of the upper corner
  end type rectangular

! the subarray indices to ghost and domain areas used for boundary exchange
! ('c' for copy, 'p' for prolongation, 'r' for restriction)
!
#if NDIMS == 2
  type(rectangular), dimension(2,2,NDIMS)  , save :: edges_dc  , edges_gc
  type(rectangular), dimension(2,2,NDIMS)  , save :: edges_dp  , edges_gp
  type(rectangular), dimension(2,2,NDIMS)  , save :: edges_dr  , edges_gr
  type(rectangular), dimension(2,2)        , save :: corners_dc, corners_gc
  type(rectangular), dimension(2,2)        , save :: corners_dp, corners_gp
  type(rectangular), dimension(2,2)        , save :: corners_dr, corners_gr
#endif /* NDIMS == 2 */
#if NDIMS == 3
  type(rectangular), dimension(2,2,2,NDIMS), save :: faces_dc  , faces_gc
  type(rectangular), dimension(2,2,2,NDIMS), save :: faces_dp  , faces_gp
  type(rectangular), dimension(2,2,2,NDIMS), save :: faces_dr  , faces_gr
  type(rectangular), dimension(2,2,2,NDIMS), save :: edges_dc  , edges_gc
  type(rectangular), dimension(2,2,2,NDIMS), save :: edges_dp  , edges_gp
  type(rectangular), dimension(2,2,2,NDIMS), save :: edges_dr  , edges_gr
  type(rectangular), dimension(2,2,2)      , save :: corners_dc, corners_gc
  type(rectangular), dimension(2,2,2)      , save :: corners_dp, corners_gp
  type(rectangular), dimension(2,2,2)      , save :: corners_dr, corners_gr
#endif /* NDIMS == 3 */

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
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_coordinates(verbose, iret)

! include external procedures
!
    use parameters, only : get_parameter_integer, get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    integer :: i, j, k, l, p, q, r, ff
    integer :: fi, fj, fk
    integer :: ni, nj, nk, nm, np, nr, ns
    logical :: info

! local arrays
!
    integer(kind=4), dimension(3) :: cm, rm, dm
!
!-------------------------------------------------------------------------------
!
! obtain the number of cells along each block dimension
!
    call get_parameter_integer("ncells" , nc )

! obtain the number of ghost zones
!
    call get_parameter_integer("nghosts", ng )

! set the block dimensions
!
    in = nc
    jn = nc
#if NDIMS == 3
    kn = nc
#endif /* NDIMS == 3 */

! calculate half and double of the number of ghose zones
!
    nh = ng / 2
    nd = ng * 2

! calculate the block dimensions including ghost cells
!
    im = in + 2 * ng
    jm = jn + 2 * ng
#if NDIMS == 3
    km = kn + 2 * ng
#endif /* NDIMS == 3 */

! prepare indices
!
    it = im - nh + 1
    jt = jm - nh + 1
#if NDIMS == 3
    kt = km - nh + 1
#endif /* NDIMS == 3 */

! calculate block indices
!
    ih  = im / 2
    ib  = ng +  1
    ie  = ng + in
    ibl = ib - 1
    ibu = ib + ng - 1
    iel = ie - ng + 1
    ieu = ie + 1
    jh  = jm / 2
    jb  = ng +  1
    je  = ng + jn
    jbl = jb - 1
    jbu = jb + ng - 1
    jel = je - ng + 1
    jeu = je + 1
#if NDIMS == 3
    kh  = km / 2
    kb  = ng +  1
    ke  = ng + kn
    kbl = kb - 1
    kbu = kb + ng - 1
    kel = ke - ng + 1
    keu = ke + 1
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
#if NDIMS == 3
    call get_parameter_real   ("zmin"  , zmin  )
    call get_parameter_real   ("zmax"  , zmax  )
#endif /* NDIMS == 3 */

! calculate the domain sizes
!
    xlen = xmax - xmin
    ylen = ymax - ymin
#if NDIMS == 3
    zlen = zmax - zmin
#endif /* NDIMS == 3 */

! calculate the domain volume
!
    vol  = xlen * ylen * zlen
    voli = 1.0d+00 / vol

! calculate the boundary areas
!
    xarea = ylen * zlen
    yarea = xlen * zlen
    zarea = xlen * ylen

! allocate space for coordinate variables
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

! reset all coordinate variables to initial values
!
    ax   (:,:) = 0.0d+00
    ay   (:,:) = 0.0d+00
    az   (:,:) = 0.0d+00
    adx  (:)   = 1.0d+00
    ady  (:)   = 1.0d+00
    adz  (:)   = 1.0d+00
    adr  (:)   = 1.0d+00
    adxi (:)   = 1.0d+00
    adyi (:)   = 1.0d+00
    adzi (:)   = 1.0d+00
    advol(:)   = 1.0d+00

! generate the coordinate variables for each level
!
    do l = 1, toplev

! calculate the block resolution at each level
!
      ff      = 2**(l - 1)
      ni      = in * ff
      nj      = jn * ff
      nk      = kn * ff

! calculate the cell sizes for each level
!
      adx (l) = xlen / (ir * ni)
      ady (l) = ylen / (jr * nj)
#if NDIMS == 3
      adz (l) = zlen / (kr * nk)
#endif /* NDIMS == 3 */
#if NDIMS == 2
      adr (l) = sqrt(adx(l)**2 + ady(l)**2)
#endif /* NDIMS == 2 */
#if NDIMS == 3
      adr (l) = sqrt(adx(l)**2 + ady(l)**2 + adz(l)**2)
#endif /* NDIMS == 3 */

! calculate the inverse of cell size
!
      adxi(l) = 1.0d+00 / adx(l)
      adyi(l) = 1.0d+00 / ady(l)
#if NDIMS == 3
      adzi(l) = 1.0d+00 / adz(l)
#endif /* NDIMS == 3 */

! calculate the block coordinates for each level
!
      ax(l,:) = ((/(i, i = 1, im)/) - ng - 5.0d-01) * adx(l)
      ay(l,:) = ((/(j, j = 1, jm)/) - ng - 5.0d-01) * ady(l)
#if NDIMS == 3
      az(l,:) = ((/(k, k = 1, km)/) - ng - 5.0d-01) * adz(l)
#endif /* NDIMS == 3 */

! calculate the cell volume at each level
!
      advol(l) = adx(l) * ady(l) * adz(l)

    end do ! l = 1, toplev

! initialize ghost subarray indices
!
    np = nc + ng
    nm = nc - ng
    nr = nc - nd
    ns = nc / 2
#if NDIMS == 2
    do j = 1, 2
      fj = j - 1
      q  = 3 - j
      do i = 1, 2
        fi = i - 1
        p  = 3 - i

! for edges copy
!
        edges_gc(i,j,1)%l(1) = ib + fi * ns
        edges_gc(i,j,1)%l(2) =  1 + fj * np
        edges_gc(i,j,2)%l(1) =  1 + fi * np
        edges_gc(i,j,2)%l(2) = jb + fj * ns

        edges_dc(i,q,1)%l(1) = ib + fi * ns
        edges_dc(i,q,1)%l(2) = jb + fj * nm
        edges_dc(p,j,2)%l(1) = ib + fi * nm
        edges_dc(p,j,2)%l(2) = jb + fj * ns

        edges_gc(i,j,1)%u(:) = edges_gc(i,j,1)%l(:) + (/ ns, ng /) - 1
        edges_gc(i,j,2)%u(:) = edges_gc(i,j,2)%l(:) + (/ ng, ns /) - 1

        edges_dc(i,q,1)%u(:) = edges_dc(i,q,1)%l(:) + (/ ns, ng /) - 1
        edges_dc(p,j,2)%u(:) = edges_dc(p,j,2)%l(:) + (/ ng, ns /) - 1

! for edges restrict
!
        edges_gr(i,j,1)%l(1) = ib + fi * ns
        edges_gr(i,j,1)%l(2) =  1 + fj * np
        edges_gr(i,j,2)%l(1) =  1 + fi * np
        edges_gr(i,j,2)%l(2) = jb + fj * ns

        edges_dr(i,q,1)%l(1) = ib
        edges_dr(i,q,1)%l(2) = jb + fj * nr
        edges_dr(p,j,2)%l(1) = ib + fi * nr
        edges_dr(p,j,2)%l(2) = jb

        edges_gr(i,j,1)%u(:) = edges_gr(i,j,1)%l(:) + (/ ns, ng /) - 1
        edges_gr(i,j,2)%u(:) = edges_gr(i,j,2)%l(:) + (/ ng, ns /) - 1

        edges_dr(i,q,1)%u(:) = edges_dr(i,q,1)%l(:) + (/ nc, nd /) - 1
        edges_dr(p,j,2)%u(:) = edges_dr(p,j,2)%l(:) + (/ nd, nc /) - 1

! for corners copy
!
        corners_gc(i,j)%l(1) =  1 + fi * np
        corners_gc(i,j)%l(2) =  1 + fj * np

        corners_dc(p,q)%l(1) = ib + fi * nm
        corners_dc(p,q)%l(2) = jb + fj * nm

        corners_gc(i,j)%u(:) = corners_gc(i,j)%l(:) + ng - 1

        corners_dc(p,q)%u(:) = corners_dc(p,q)%l(:) + ng - 1

! for corners restrict
!
        corners_gr(i,j)%l(1) =  1 + fi * np
        corners_gr(i,j)%l(2) =  1 + fj * np

        corners_dr(p,q)%l(1) = ib + fi * nr
        corners_dr(p,q)%l(2) = jb + fj * nr

        corners_gr(i,j)%u(:) = corners_gr(i,j)%l(:) + ng - 1

        corners_dr(p,q)%u(:) = corners_dr(p,q)%l(:) + nd - 1

      end do ! i = 1, 2
    end do ! j = 1, 2
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = 1, 2
      fk = k - 1
      r  = 3 - k
      do j = 1, 2
        fj = j - 1
        q  = 3 - j
        do i = 1, 2
          fi = i - 1
          p  = 3 - i

! for faces copy
!
          faces_gc(i,j,k,1)%l(1) =  1 + fi * np
          faces_gc(i,j,k,1)%l(2) = jb + fj * ns
          faces_gc(i,j,k,1)%l(3) = kb + fk * ns
          faces_gc(i,j,k,2)%l(1) = ib + fi * ns
          faces_gc(i,j,k,2)%l(2) =  1 + fj * np
          faces_gc(i,j,k,2)%l(3) = kb + fk * ns
          faces_gc(i,j,k,3)%l(1) = ib + fi * ns
          faces_gc(i,j,k,3)%l(2) = jb + fj * ns
          faces_gc(i,j,k,3)%l(3) =  1 + fk * np

          faces_dc(p,j,k,1)%l(1) = ib + fi * nm
          faces_dc(p,j,k,1)%l(2) = jb + fj * ns
          faces_dc(p,j,k,1)%l(3) = kb + fk * ns
          faces_dc(i,q,k,2)%l(1) = ib + fi * ns
          faces_dc(i,q,k,2)%l(2) = jb + fj * nm
          faces_dc(i,q,k,2)%l(3) = kb + fk * ns
          faces_dc(i,j,r,3)%l(1) = ib + fi * ns
          faces_dc(i,j,r,3)%l(2) = jb + fj * ns
          faces_dc(i,j,r,3)%l(3) = kb + fk * nm

          faces_gc(i,j,k,1)%u(:) = faces_gc(i,j,k,1)%l(:) + (/ ng, ns, ns /) - 1
          faces_gc(i,j,k,2)%u(:) = faces_gc(i,j,k,2)%l(:) + (/ ns, ng, ns /) - 1
          faces_gc(i,j,k,3)%u(:) = faces_gc(i,j,k,3)%l(:) + (/ ns, ns, ng /) - 1

          faces_dc(p,j,k,1)%u(:) = faces_dc(p,j,k,1)%l(:) + (/ ng, ns, ns /) - 1
          faces_dc(i,q,k,2)%u(:) = faces_dc(i,q,k,2)%l(:) + (/ ns, ng, ns /) - 1
          faces_dc(i,j,r,3)%u(:) = faces_dc(i,j,r,3)%l(:) + (/ ns, ns, ng /) - 1

! for faces restrict
!
          faces_gr(i,j,k,1)%l(1) =  1 + fi * np
          faces_gr(i,j,k,1)%l(2) = jb + fj * ns
          faces_gr(i,j,k,1)%l(3) = kb + fk * ns
          faces_gr(i,j,k,2)%l(1) = ib + fi * ns
          faces_gr(i,j,k,2)%l(2) =  1 + fj * np
          faces_gr(i,j,k,2)%l(3) = kb + fk * ns
          faces_gr(i,j,k,3)%l(1) = ib + fi * ns
          faces_gr(i,j,k,3)%l(2) = jb + fj * ns
          faces_gr(i,j,k,3)%l(3) =  1 + fk * np

          faces_dr(p,j,k,1)%l(1) = ib + fi * nr
          faces_dr(p,j,k,1)%l(2) = jb
          faces_dr(p,j,k,1)%l(3) = kb
          faces_dr(i,q,k,2)%l(1) = ib
          faces_dr(i,q,k,2)%l(2) = jb + fj * nr
          faces_dr(i,q,k,2)%l(3) = kb
          faces_dr(i,j,r,3)%l(1) = ib
          faces_dr(i,j,r,3)%l(2) = jb
          faces_dr(i,j,r,3)%l(3) = kb + fk * nr

          faces_gr(i,j,k,1)%u(:) = faces_gr(i,j,k,1)%l(:) + (/ ng, ns, ns /) - 1
          faces_gr(i,j,k,2)%u(:) = faces_gr(i,j,k,2)%l(:) + (/ ns, ng, ns /) - 1
          faces_gr(i,j,k,3)%u(:) = faces_gr(i,j,k,3)%l(:) + (/ ns, ns, ng /) - 1

          faces_dr(p,j,k,1)%u(:) = faces_dr(p,j,k,1)%l(:) + (/ nd, nc, nc /) - 1
          faces_dr(i,q,k,2)%u(:) = faces_dr(i,q,k,2)%l(:) + (/ nc, nd, nc /) - 1
          faces_dr(i,j,r,3)%u(:) = faces_dr(i,j,r,3)%l(:) + (/ nc, nc, nd /) - 1

! for edges copy
!
          edges_gc(i,j,k,1)%l(1) = ib + fi * ns
          edges_gc(i,j,k,1)%l(2) =  1 + fj * np
          edges_gc(i,j,k,1)%l(3) =  1 + fk * np
          edges_gc(i,j,k,2)%l(1) =  1 + fi * np
          edges_gc(i,j,k,2)%l(2) = jb + fj * ns
          edges_gc(i,j,k,2)%l(3) =  1 + fk * np
          edges_gc(i,j,k,3)%l(1) =  1 + fi * np
          edges_gc(i,j,k,3)%l(2) =  1 + fj * np
          edges_gc(i,j,k,3)%l(3) = kb + fk * ns

          edges_dc(i,q,r,1)%l(1) = ib + fi * ns
          edges_dc(i,q,r,1)%l(2) = jb + fj * nm
          edges_dc(i,q,r,1)%l(3) = kb + fk * nm
          edges_dc(p,j,r,2)%l(1) = jb + fi * nm
          edges_dc(p,j,r,2)%l(2) = jb + fj * ns
          edges_dc(p,j,r,2)%l(3) = kb + fk * nm
          edges_dc(p,q,k,3)%l(1) = ib + fi * nm
          edges_dc(p,q,k,3)%l(2) = jb + fj * nm
          edges_dc(p,q,k,3)%l(3) = kb + fk * ns

          edges_gc(i,j,k,1)%u(:) = edges_gc(i,j,k,1)%l(:) + (/ ns, ng, ng /) - 1
          edges_gc(i,j,k,2)%u(:) = edges_gc(i,j,k,2)%l(:) + (/ ng, ns, ng /) - 1
          edges_gc(i,j,k,3)%u(:) = edges_gc(i,j,k,3)%l(:) + (/ ng, ng, ns /) - 1

          edges_dc(i,q,r,1)%u(:) = edges_dc(i,q,r,1)%l(:) + (/ ns, ng, ng /) - 1
          edges_dc(p,j,r,2)%u(:) = edges_dc(p,j,r,2)%l(:) + (/ ng, ns, ng /) - 1
          edges_dc(p,q,k,3)%u(:) = edges_dc(p,q,k,3)%l(:) + (/ ng, ng, ns /) - 1

! for edges restrict
!
          edges_gr(i,j,k,1)%l(1) = ib + fi * ns
          edges_gr(i,j,k,1)%l(2) =  1 + fj * np
          edges_gr(i,j,k,1)%l(3) =  1 + fk * np
          edges_gr(i,j,k,2)%l(1) =  1 + fi * np
          edges_gr(i,j,k,2)%l(2) = jb + fj * ns
          edges_gr(i,j,k,2)%l(3) =  1 + fk * np
          edges_gr(i,j,k,3)%l(1) =  1 + fi * np
          edges_gr(i,j,k,3)%l(2) =  1 + fj * np
          edges_gr(i,j,k,3)%l(3) = kb + fk * ns

          edges_dr(i,q,r,1)%l(1) = ib
          edges_dr(i,q,r,1)%l(2) = jb + fj * nr
          edges_dr(i,q,r,1)%l(3) = kb + fk * nr
          edges_dr(p,j,r,2)%l(1) = ib + fi * nr
          edges_dr(p,j,r,2)%l(2) = jb
          edges_dr(p,j,r,2)%l(3) = kb + fk * nr
          edges_dr(p,q,k,3)%l(1) = ib + fi * nr
          edges_dr(p,q,k,3)%l(2) = jb + fj * nr
          edges_dr(p,q,k,3)%l(3) = kb

          edges_gr(i,j,k,1)%u(:) = edges_gr(i,j,k,1)%l(:) + (/ ns, ng, ng /) - 1
          edges_gr(i,j,k,2)%u(:) = edges_gr(i,j,k,2)%l(:) + (/ ng, ns, ng /) - 1
          edges_gr(i,j,k,3)%u(:) = edges_gr(i,j,k,3)%l(:) + (/ ng, ng, ns /) - 1

          edges_dr(i,q,r,1)%u(:) = edges_dr(i,q,r,1)%l(:) + (/ nc, nd, nd /) - 1
          edges_dr(p,j,r,2)%u(:) = edges_dr(p,j,r,2)%l(:) + (/ nd, nc, nd /) - 1
          edges_dr(p,q,k,3)%u(:) = edges_dr(p,q,k,3)%l(:) + (/ nd, nd, nc /) - 1

! for corners copy
!
          corners_gc(i,j,k)%l(1) =  1 + fi * np
          corners_gc(i,j,k)%l(2) =  1 + fj * np
          corners_gc(i,j,k)%l(3) =  1 + fk * np

          corners_dc(p,q,r)%l(1) = ib + fi * nm
          corners_dc(p,q,r)%l(2) = jb + fj * nm
          corners_dc(p,q,r)%l(3) = kb + fk * nm

          corners_gc(i,j,k)%u(:) = corners_gc(i,j,k)%l(:) + ng - 1

          corners_dc(p,q,r)%u(:) = corners_dc(p,q,r)%l(:) + ng - 1

! for corners restrict
!
          corners_gr(i,j,k)%l(1) =  1 + fi * np
          corners_gr(i,j,k)%l(2) =  1 + fj * np
          corners_gr(i,j,k)%l(3) =  1 + fk * np

          corners_dr(p,q,r)%l(1) = ib + fi * nr
          corners_dr(p,q,r)%l(2) = jb + fj * nr
          corners_dr(p,q,r)%l(3) = kb + fk * nr

          corners_gr(i,j,k)%u(:) = corners_gr(i,j,k)%l(:) + ng - 1

          corners_dr(p,q,r)%u(:) = corners_dr(p,q,r)%l(:) + nd - 1

        end do ! i = 1, 2
      end do ! j = 1, 2
    end do ! k = 1, 2
#endif /* NDIMS == 3 */

! print general information about the level resolutions
!
    if (verbose) then

! the base resolution
!
      cm(1) = ir * in
      cm(2) = jr * jn
      cm(3) = kr * kn

! the effective resolution
!
      ff    = 2**(maxlev - 1)
      rm(1) = cm(1) * ff
      rm(2) = cm(2) * ff
      rm(3) = cm(3) * ff

! the top level block division
!
      dm(1) = rm(1) / in
      dm(2) = rm(2) / jn
      dm(3) = rm(3) / kn

! obtain the maximum number of block
!
      ff    = product(dm(1:NDIMS))

! print info
!
      write(*,"(4x,a,  1x,i6 )" ) "refinement to level    =", toplev
      write(*,"(4x,a,3(1x,i6 ))") "base configuration     =", ir, jr, kr
      write(*,"(4x,a,3(1x,i6 ))") "top level blocks       =", dm(1:NDIMS)
      write(*,"(4x,a,  3x,i18)" ) "maximum cover blocks   =", ff
      write(*,"(4x,a,3(1x,i6 ))") "base resolution        =", cm(1:NDIMS)
      write(*,"(4x,a,3(1x,i6 ))") "effective resolution   =", rm(1:NDIMS)

    end if ! verbose

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
!   Arguments:
!
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine finalize_coordinates(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
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

!-------------------------------------------------------------------------------
!
  end subroutine finalize_coordinates

!===============================================================================
!
end module
