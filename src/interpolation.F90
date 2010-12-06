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
  subroutine reconstruct(n, h, f, fl, fr)

    implicit none

! input/output arguments
!
    integer           , intent(in)  :: n
    real              , intent(in)  :: h
    real, dimension(n), intent(in)  :: f
    real, dimension(n), intent(out) :: fl, fr

! local variables
!
    integer            :: i
    real               :: df, ds
    real, dimension(n) :: dfl, dfr
!
!-------------------------------------------------------------------------------
!
#ifdef TVD
!! second order TVD interpolation
!!
! calculate left and right derivatives
!
    do i = 1, n - 1
      dfr(i  ) = f(i+1) - f(i)
      dfl(i+1) = dfr(i)
    end do
    dfl(1) = dfr(1)
    dfr(n) = dfl(n)

! interpolate the values at i-1/2 and i+1/2
!
    do i = 1, n
      ds = dfr(i) * dfl(i)

      if (ds .gt. 0.0d0) then
#ifdef MINMOD
        df = sign(0.5d0, dfr(i)) * min(abs(dfr(i)), abs(dfl(i)))
#endif /* MINMOD */
#ifdef LF
        df  = ds / (dfr(i) + dfl(i))
#endif /* LF */

        fl(i) = f(i) + df
        fr(i) = f(i) - df
      else
        fl(i) = f(i)
        fr(i) = f(i)
      end if
    end do

! shift i-1/2 to the left
!
    do i = 1, n - 1
      fr(i) = fr(i+1)
    end do
    fr(n) = f(n)
#endif /* TVD */
!
!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
!===============================================================================
!
! expand: expands a multi-dimentional array similar using only TVD or high order
!         interpolation
!
!===============================================================================
!
  subroutine expand(cm, dm, ng, u, v)

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

        call expand_1d(cm(1), dm(1), ng, x(1:cm(1)), y(1:dm(1)))

        w(1:dm(1),j,k) = y(1:dm(1))
      end do
    end do

! expand in Y-direction
!
    do k = 1, cm(3)
      do i = 1, dm(1)
        x(1:cm(2)) = w(i,1:cm(2),k)

        call expand_1d(cm(2), dm(2), ng, x(1:cm(2)), y(1:dm(2)))

        z(i,1:dm(2),k) = y(1:dm(2))
      end do
    end do

! expand in Z-direction
!
    if (cm(3) .gt. 1) then
      do j = 1, dm(2)
        do i = 1, dm(1)
          x(1:cm(3)) = z(i,j,1:cm(3))

          call expand_1d(cm(3), dm(3), ng, x(1:cm(3)), y(1:dm(3)))

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
  end subroutine expand
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
! expand_1d: one dimensional expansion using the TVD or high ordred
!            interpolation
!
!===============================================================================
!
  subroutine expand_1d(nu, nv, ng, u, v)

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

#ifdef TVD
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
#endif /* TVD */

    end do
!
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
