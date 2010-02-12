!!*****************************************************************************
!!
!! module: interpolation - subroutines for different kinds of interpolation
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
module interpolation

  implicit none

  contains
!
!===============================================================================
!
! reconstruct: subroutine for the reconstruction of the values at the right and
!              left interfaces of cells from their cell centered representation
!
!===============================================================================
!
  subroutine reconstruct(n, vx, vl, vr)

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real, dimension(n), intent(in)  :: vx
    real, dimension(n), intent(out) :: vl, vr

! local variables
!
    integer            :: i
    real               :: dv
    real, dimension(n) :: ds, dvl, dvr
!
!------------------------------------------------------------------------------
!
! second order interpolation
!
    dvl(1) = 0.0
    dvr(n) = 0.0

    do i = 1, n-1
      dvr(i  ) = vx(i+1) - vx(i)
      dvl(i+1) = dvr(i)
    enddo

    do i = 1, n
      ds (i) = dvr(i) * dvl(i)

      if (ds(i) .gt. 0.0) then
        dv  = ds(i) / (dvr(i) + dvl(i))

        vl(i) = vx(i) + dv
        vr(i) = vx(i) - dv
      else
        vl(i) = vx(i)
        vr(i) = vx(i)
      endif
    enddo

    do i = 1, n-1
      vr(i) = vr(i+1)
    enddo
    vr(n) = vx(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
!===============================================================================
!
! expand: Expands multi-dimentional array similar to EXPAND for 2D images
!         but using second-order TVD interpolation. EXPAND_TVD expands
!         array only in one selected direction at time.
!
!===============================================================================
!
  subroutine expand(dm, fm, ng, u, v, xflag, yflag, zflag)

    implicit none

! input parameters
!
    integer, dimension(3)         , intent(in)  :: dm, fm
    integer                       , intent(in)  :: ng
    real        , dimension(:,:,:), intent(in)  :: u
    real        , dimension(:,:,:), intent(out) :: v
    character                     , intent(in)  :: xflag, yflag, zflag

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real(kind=8), dimension(:)    , allocatable :: x, y
    real(kind=8), dimension(:,:,:), allocatable :: w, z
!
!-------------------------------------------------------------------------------
!
    v(:,:,:) = 0.0

    allocate(x(maxval(dm)))
    allocate(y(maxval(fm)))
    allocate(w(fm(1),dm(2),dm(3)))
    allocate(z(fm(1),fm(2),dm(3)))

! expand in X direction
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        x(1:dm(1)) = u(1:dm(1),j,k)

        call expand_1d(dm(1),fm(1),ng,x(1:dm(1)),y(1:fm(1)),xflag)

        w(1:fm(1),j,k) = y(1:fm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        x(1:dm(2)) = w(i,1:dm(2),k)

        call expand_1d(dm(2),fm(2),ng,x(1:dm(2)),y(1:fm(2)),yflag)

        z(i,1:fm(2),k) = y(1:fm(2))
      end do
    end do

! expand in Z-direction
!
    if (dm(3) .gt. 1) then
      do j = 1, fm(2)
        do i = 1, fm(1)
          x(1:dm(3)) = z(i,j,1:dm(3))

          call expand_1d(dm(3),fm(3),ng,x(1:dm(3)),y(1:fm(3)),zflag)

          v(i,j,1:fm(3)) = y(1:fm(3))
        end do
      end do
    else
      do k = 1, fm(3)
        do j = 1, fm(2)
          do i = 1, fm(1)
            v(i,j,k) = z(i,j,k)
          end do
        end do
      end do
    endif

    deallocate(w)
    deallocate(z)
    deallocate(x)
    deallocate(y)

!-------------------------------------------------------------------------------
!
  end subroutine expand
!
!===============================================================================
!
! shrink: shrinks multi-dimentional array using different interpolations
!
!===============================================================================
!
  subroutine shrink(dm, fm, ng, u, v, xflag, yflag, zflag)

    implicit none

! input parameters
!
    integer, dimension(3)         , intent(in)  :: dm, fm
    integer                       , intent(in)  :: ng
    real        , dimension(:,:,:), intent(in)  :: u
    real        , dimension(:,:,:), intent(out) :: v
    character                     , intent(in)  :: xflag, yflag, zflag

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real(kind=8), dimension(:)    , allocatable :: x, y
    real(kind=8), dimension(:,:,:), allocatable :: w, z
!
!-------------------------------------------------------------------------------
!
    v(:,:,:) = 0.0

    allocate(x(maxval(dm)))
    allocate(y(maxval(fm)))
    allocate(w(fm(1),dm(2),dm(3)))
    allocate(z(fm(1),fm(2),dm(3)))

! expand in X direction
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        x(1:dm(1)) = u(1:dm(1),j,k)

        call shrink_1d(dm(1),fm(1),ng,x(1:dm(1)),y(1:fm(1)),xflag)

        w(1:fm(1),j,k) = y(1:fm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        x(1:dm(2)) = w(i,1:dm(2),k)

        call shrink_1d(dm(2),fm(2),ng,x(1:dm(2)),y(1:fm(2)),yflag)

        z(i,1:fm(2),k) = y(1:fm(2))
      end do
    end do

! expand in Z-direction
!
    if (dm(3) .gt. 1) then
      do j = 1, fm(2)
        do i = 1, fm(1)
          x(1:dm(3)) = z(i,j,1:dm(3))

          call shrink_1d(dm(3),fm(3),ng,x(1:dm(3)),y(1:fm(3)),zflag)

          v(i,j,1:fm(3)) = y(1:fm(3))
        end do
      end do
    else
      do k = 1, fm(3)
        do j = 1, fm(2)
          do i = 1, fm(1)
            v(i,j,k) = z(i,j,k)
          end do
        end do
      end do
    endif

    deallocate(w)
    deallocate(z)
    deallocate(x)
    deallocate(y)

!-------------------------------------------------------------------------------
!
  end subroutine shrink
!
!===============================================================================
!
! expand_1d: one dimensional expansion using different interpolations for
!            different location of the new points
!
!===============================================================================
!
  subroutine expand_1d(n, m, ng, u, v, flag)

    implicit none

! input parameters
!
    integer                   , intent(in)  :: n, m, ng
    real        , dimension(n), intent(in)  :: u
    real        , dimension(m), intent(out) :: v
    character                 , intent(in)  :: flag

! local variables
!
    integer      :: i, ib, ie, j0, j1
    real(kind=8) :: du, ddu, dum, dup, du0, ds
!
!-------------------------------------------------------------------------------
!
    v(:) = 0.0

    ib = ng / 2 + 1
    ie = n - ng / 2

    select case(flag)
    case('m')
      do i = ib, ie
        j1 = 2 * i - ng
        j0 = j1 - 1

        v(j0) = u(i)
        v(j1) = u(i)
      end do
    case('l')
      do i = ib, ie
        du = 0.5 * (u(i+1) - u(i-1))

        j1 = 2 * i - ng
        j0 = j1 - 1

        v(j0) = u(i) - 0.25 * du
        v(j1) = u(i) + 0.25 * du
      end do
    case('c')
      do i = ib, ie
        du  = u(i) - u(i-1)
        ddu = - 0.25 * (u(i+1) - u(i-1) - u(i) + u(i-2))

        j1 = 2 * i - ng
        j0 = j1 - 1

        v(j0) = u(i) - 0.5 * du + 0.25 * ddu
        v(j1) = u(i)
      end do
    case('u')
      do i = ib, ie
        dup = u(i+1) - u(i)
        dum = u(i) - u(i-1)

        j1 = 2 * i - ng
        j0 = j1 - 1

        v(j0) = u(i) - 0.25d0*dum
        v(j1) = u(i) + 0.25d0*dup
      end do
    case default
      do i = ib, ie
        dup = u(i+1) - u(i)
        dum = u(i) - u(i-1)
        du0 = dup + dum
        ds  = dup * dum

        j1 = 2 * i - ng
        j0 = j1 - 1

        if (ds .gt. 0.0) then
          du = ds / du0

          v(j0) = u(i) - 0.25 * du
          v(j1) = u(i) + 0.25 * du
        else
          v(j0) = u(i)
          v(j1) = u(i)
        endif
      end do
    end select

!-------------------------------------------------------------------------------
!
  end subroutine expand_1d
!
!===============================================================================
!
! shrink_1d: one dimensional shrinking using different interpolations for
!            different location of the new points
!
!===============================================================================
!
  subroutine shrink_1d(n, m, ng, u, v, flag)

    implicit none

! input parameters
!
    integer                   , intent(in)  :: n, m, ng
    real        , dimension(n), intent(in)  :: u
    real        , dimension(m), intent(out) :: v
    character                 , intent(in)  :: flag

! local variables
!
    integer      :: i, ib, ie, j0, j1
    real(kind=8) :: du, ddu, dum, dup, du0, ds
!
!-------------------------------------------------------------------------------
!
    v(:) = 0.0

    ib = ng / 2 + 1
    ie = m - ng / 2

    select case(flag)
    case('l', 'u')
      do i = ib, ie
        j1    = 2 * i - ng
        j0    = j1 - 1

        v(i)  = 0.5 * (u(j0) + u(j1))
      end do
    case('c')
      do i = ib, ie
        j1 = 2 * i - ng
        v(i)  = u(j1)
      end do
    case default
      do i = ib, ie
        j1    = 2 * i - ng
        j0    = j1 - 1

        v(i)  = 0.5 * (u(j0) + u(j1))
      end do
    end select

!-------------------------------------------------------------------------------
!
  end subroutine shrink_1d
!
!===============================================================================
!
! magtocen: subroutine for interpolation of magnetic field from cell interface
!           to cell center
!             [Bx(i+1/2,j,k), By(i,j+1/2,k), Bz(i,j,k+1/2)] ->
!             [Bx(i,j,k)    , By(i,j,k)    , Bz(i,j,k)]
!
!===============================================================================
!
  subroutine magtocen(im, jm, km, bi, bc)

    implicit none

! input and output variables
!
    integer                    , intent(in)    :: im, jm, km
    real, dimension(3,im,jm,km), intent(in)    :: bi
    real, dimension(3,im,jm,km), intent(inout) :: bc

! local  variables
!
    integer :: i, j, k, n

    real, dimension(:), allocatable, save :: al
!
!------------------------------------------------------------------------------
!
! allocate a vector for interpolation
!
    allocate(al(max(im,jm,km)))

! perform interpolation
!
    do k = 1, km
      do j = 1, jm
        al(1:im) = bi(1,1:im,j,k)

        bc(1,2:im,j,k) = 0.5 * (al(1:im-1) + al(2:im))
        bc(1,1   ,j,k) =        al(1     )
      enddo

      do i = 1, im
        al(1:jm) = bi(2,i,1:jm,k)

        bc(2,i,2:jm,k) = 0.5 * (al(1:jm-1) + al(2:jm))
        bc(2,i,1   ,k) =        al(1     )
      enddo
    enddo

#if NDIMS == 3
    do j = 1, jm
      do i = 1, im
        al(1:km) = bi(3,i,j,1:km)


        bc(3,i,j,2:km) = 0.5 * (al(1:km-1) + ar(1:km-1))
        bc(3,i,j,1   ) =        al(1     )
      enddo
    enddo
#endif /* NDIMS == 3 */

! deallocate temporary vector
!
    deallocate(al)

  end subroutine magtocen

!===============================================================================
!
end module
