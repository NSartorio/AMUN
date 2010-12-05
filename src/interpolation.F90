!!******************************************************************************
!!
!! module: interpolation - subroutines for different kinds of interpolation
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
!-------------------------------------------------------------------------------
!
! second order interpolation
!
    do i = 1, n-1
      dvr(i  ) = vx(i+1) - vx(i)
      dvl(i+1) = dvr(i)
    enddo

    dvl(1) = dvr(1)
    dvr(n) = dvl(n)

    do i = 1, n
      ds (i) = dvr(i) * dvl(i)

      if (ds(i) .gt. 0.0) then
#ifdef MINMOD
        dv  = sign(0.5, dvr(i)) * min(abs(dvr(i)), abs(dvl(i)))
#endif /* MINMOD */
#ifdef LF
        dv  = ds(i) / (dvr(i) + dvl(i))
#endif /* LF */

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
    real, dimension(:)    , allocatable :: x, y
    real, dimension(:,:,:), allocatable :: w, z
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
! expand_tvd: expands a multi-dimentional array similar to EXPAND but using
!             only the second-order TVD interpolation
!
!===============================================================================
!
  subroutine expand_tvd(cm, dm, ng, u, v)

    implicit none

! input parameters
!
    integer, dimension(3)                     , intent(in)  :: cm, dm
    integer                                   , intent(in)  :: ng
    real        , dimension(cm(1),cm(2),cm(3)), intent(in)  :: u
    real        , dimension(dm(1),dm(2),dm(3)), intent(out) :: v

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real, dimension(:)    , allocatable :: x, y
    real, dimension(:,:,:), allocatable :: w, z
!
!-------------------------------------------------------------------------------
!
! allocate temporary arrays
!
    allocate(x(maxval(cm)))
    allocate(y(maxval(dm)))
    allocate(w(dm(1),cm(2),cm(3)))
    allocate(z(dm(1),dm(2),cm(3)))

! expand in X direction
!
    do k = 1, cm(3)
      do j = 1, cm(2)
        x(1:cm(1)) = u(1:cm(1),j,k)

        call expand_1d_tvd(cm(1), dm(1), ng, x(1:cm(1)), y(1:dm(1)))

        w(1:dm(1),j,k) = y(1:dm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, cm(3)
      do i = 1, dm(1)
        x(1:cm(2)) = w(i,1:cm(2),k)

        call expand_1d_tvd(cm(2), dm(2), ng, x(1:cm(2)), y(1:dm(2)))

        z(i,1:dm(2),k) = y(1:dm(2))
      end do
    end do

! expand in Z-direction
!
    if (cm(3) .gt. 1) then
      do j = 1, dm(2)
        do i = 1, dm(1)
          x(1:cm(3)) = z(i,j,1:cm(3))

          call expand_1d_tvd(cm(3), dm(3), ng, x(1:cm(3)), y(1:dm(3)))

          v(i,j,1:dm(3)) = y(1:dm(3))
        end do
      end do
    else
      v(:,:,:) = z(:,:,:)
    endif

! deallocate temporary arrays
!
    deallocate(w)
    deallocate(z)
    deallocate(x)
    deallocate(y)
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_tvd
!
!===============================================================================
!
! shrink: shrinks multi-dimentional array using different interpolations
!
!===============================================================================
!
  subroutine shrink(dm, fm, u, v, xflag, yflag, zflag)

    implicit none

! input parameters
!
    integer, dimension(3)         , intent(in)  :: dm, fm
    real        , dimension(:,:,:), intent(in)  :: u
    real        , dimension(:,:,:), intent(out) :: v
    character                     , intent(in)  :: xflag, yflag, zflag

! local variables
!
    integer :: i, j, k

! allocatable variables
!
    real, dimension(:)    , allocatable :: x, y
    real, dimension(:,:,:), allocatable :: w, z
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

        call shrink_1d(dm(1), fm(1), x(1:dm(1)), y(1:fm(1)), xflag)

        w(1:fm(1),j,k) = y(1:fm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        x(1:dm(2)) = w(i,1:dm(2),k)

        call shrink_1d(dm(2), fm(2), x(1:dm(2)), y(1:fm(2)), yflag)

        z(i,1:fm(2),k) = y(1:fm(2))
      end do
    end do

! expand in Z-direction
!
    if (dm(3) .gt. 1) then
      do j = 1, fm(2)
        do i = 1, fm(1)
          x(1:dm(3)) = z(i,j,1:dm(3))

          call shrink_1d(dm(3), fm(3), x(1:dm(3)), y(1:fm(3)), zflag)

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
    integer  :: i, ib, ie, il, ir
    real     :: du, dum, dup, du0, ds
!
!-------------------------------------------------------------------------------
!
    v(:) = 0.0

    ib = ng / 2 + 1
    ie = n - ng / 2

    select case(flag)
    case('c')
      do i = ib, ie
        du  = u(i) - u(i-1)

        ir = 2 * i - ng
        il = ir - 1

        v(il) = u(i) - 0.5 * du
        v(ir) = u(i)
      end do
    case('t')
      do i = ib, ie

        ir = 2 * i - ng
        il = ir - 1

        dup = u(i+1) - u(i)
        dum = u(i-1) - u(i)
        ds  = - dup * dum

        if (ds .gt. 0.0) then
          du    = 0.5 * ds / (dup - dum)

          v(il) = u(i) - du
          v(ir) = u(i) + du
        else
          v(il) = u(i)
          v(ir) = u(i)
        end if
      end do
    case default
      do i = ib, ie
        ir = 2 * i - ng
        il = ir - 1

        v(il) = u(i)
        v(ir) = u(i)
      end do
    end select

!-------------------------------------------------------------------------------
!
  end subroutine expand_1d
!
!===============================================================================
!
! expand_1d_tvd: one dimensional expansion using the TVD interpolation
!
!===============================================================================
!
  subroutine expand_1d_tvd(nu, nv, ng, u, v)

    implicit none

! input parameters
!
    integer                    , intent(in)    :: nu, nv, ng
    real        , dimension(nu), intent(in)    :: u
    real        , dimension(nv), intent(inout) :: v

! local variables
!
    integer  :: i, ib, ie, il, ir
    real     :: dum, dup, ds, du
!
!-------------------------------------------------------------------------------
!
    ib =  1 + ng
    ie = nu - ng

    do i = ib, ie

      ir = 2 * (i - ng)
      il = ir - 1

      dup = u(i+1) - u(i)
      dum = u(i) - u(i-1)
      ds  = dup * dum

      if (ds .gt. 0.0) then
#ifdef MINMOD
        du  = sign(0.25, dup) * min(abs(dup), abs(dum))
#endif /* MINMOD */
#ifdef LF
        du    = 0.5 * ds / (dup + dum)
#endif /* LF */
        v(il) = u(i) - du
        v(ir) = u(i) + du
      else
        v(il) = u(i)
        v(ir) = u(i)
      end if

    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_1d_tvd
!
!===============================================================================
!
! shrink_1d: one dimensional shrinking using different interpolations for
!            different location of the new points
!
!===============================================================================
!
  subroutine shrink_1d(n, m, u, v, flag)

    implicit none

! input parameters
!
    integer                   , intent(in)  :: n, m
    real        , dimension(n), intent(in)  :: u
    real        , dimension(m), intent(out) :: v
    character                 , intent(in)  :: flag

! local variables
!
    integer      :: i, il, ir
!
!-------------------------------------------------------------------------------
!
    v(:) = 0.0

    select case(flag)
    case('c')
      do i = 1, m
        ir   = 2 * i
        v(i) = u(ir)
      end do
    case default
      do i = 1, m
        ir   = 2 * i
        il   = ir - 1

        v(i) = 0.5 * (u(il) + u(ir))
      end do
    end select
!
!-------------------------------------------------------------------------------
!
  end subroutine shrink_1d
#ifdef MHD
!
!===============================================================================
!
! divergence: function calculates divergence of the input staggered field
!             [Bx(i+1/2,j,k), By(i,j+1/2,k), Bz(i,j,k+1/2)] -> div B
!
!===============================================================================
!
  subroutine divergence(b, db, dx, dy, dz)

    use config   , only : im, jm, km
    use variables, only : nvr

    implicit none

! input and output variables
!
    real, dimension(3,im,jm,km), intent(in)  :: b
    real, dimension(  im,jm,km), intent(out) :: db
    real, optional             , intent(in)  :: dx, dy, dz

! local  variables
!
    integer :: i, j, k
    real    :: dxi, dyi, dzi
!
!------------------------------------------------------------------------------
!
! check optional arguments
!
    dxi = 1.0
    dyi = 1.0
    dzi = 1.0

    if (present(dx)) dxi = 1.0 / dx
    if (present(dy)) dyi = 1.0 / dy
    if (present(dz)) dzi = 1.0 / dz

! reset output array
!
    db(:,:,:) = 0.0

! iterate over all points
!
#if NDIMS == 3
    do k = 2, km
      do j = 2, jm
        do i = 2, im
          db(i,j,k) = dxi * (b(1,i,j,k) - b(1,i-1,j,k))                        &
                    + dyi * (b(2,i,j,k) - b(2,i,j-1,k))                        &
                    + dzi * (b(3,i,j,k) - b(3,i,j,k-1))
#else /* NDIMS == 3 */
    do k = 1, km
      do j = 2, jm
        do i = 2, im
          db(i,j,k) = dxi * (b(1,i,j,k) - b(1,i-1,j,k))                        &
                    + dyi * (b(2,i,j,k) - b(2,i,j-1,k))
#endif /* NDIMS == 3 */
        end do
      end do
    end do
!
  end subroutine divergence
#endif /* MHD */
!
!===============================================================================
!
end module
