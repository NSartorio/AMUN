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
    do i = 1, n-1
      dvr(i  ) = vx(i+1) - vx(i)
      dvl(i+1) = dvr(i)
    enddo

    dvl(1) = dvr(1)
    dvr(n) = dvl(n)

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
        du    = 0.5 * ds / (dup + dum)
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
#ifdef FLUXCT
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
  subroutine magtocen(u)

    use blocks, only : nvr
    use blocks, only : ibx, iby, ibz, icx, icy, icz
    use config, only : im, jm, km

    implicit none

! input and output variables
!
    real, dimension(nvr,im,jm,km), intent(inout) :: u

! local  variables
!
    integer :: i, j, k, im1, ip1, jm1, jp1, km1, kp1
    real    :: dur, dul, du, ds

! local arrays
!
    real, dimension(im) :: ax
    real, dimension(jm) :: ay
#if NDIMS == 3
    real, dimension(km) :: az
#endif /* NDIMS == 3 */
!
!------------------------------------------------------------------------------
!
! perform interpolation
!
    do k = 1, km
      do j = 1, jm
        ax(1:im) = u(ibx,1:im,j,k)

        u(icx,2:im,j,k) = 0.5 * (ax(1:im-1) + ax(2:im))
        u(icx,1   ,j,k) =        ax(1     )
      enddo

      do i = 1, im
        ay(1:jm) = u(iby,i,1:jm,k)

        u(icy,i,2:jm,k) = 0.5 * (ay(1:jm-1) + ay(2:jm))
        u(icy,i,1   ,k) =        ay(1     )
      enddo
    enddo

#if NDIMS == 3
    do j = 1, jm
      do i = 1, im
        az(1:km) = u(ibz,i,j,1:km)

        u(icz,i,j,2:km) = 0.5 * (az(1:km-1) + az(2:km))
        u(icz,i,j,1   ) =        az(1     )
      enddo
    enddo
#else /* NDIMS == 3 */
    u(icz,:,:,:) = u(ibz,:,:,:)
#endif /* NDIMS == 3 */

! correct components to get the divergence free interpolation
!
#if NDIMS == 3
    do k = 2, km - 1
      km1 = k - 1
      kp1 = k + 1
      do j = 2, jm - 1
        jm1 = j - 1
        jp1 = j + 1
        do i = 2, im - 1
          im1 = i - 1
          ip1 = i + 1

! Bx correction from By
!
          dur = u(iby,ip1,j,k) - u(iby,i,j,k)
          dul = u(iby,im1,j,k) - u(iby,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icx,i,j  ,k) = u(icx,i,j  ,k) + 0.25 * du
            u(icx,i,jp1,k) = u(icx,i,jp1,k) - 0.25 * du
          end if

! Bx correction from Bz
!
          dur = u(ibz,ip1,j,k) - u(ibz,i,j,k)
          dul = u(ibz,im1,j,k) - u(ibz,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icx,i,j,k  ) = u(icx,i,j,k  ) + 0.25 * du
            u(icx,i,j,kp1) = u(icx,i,j,kp1) - 0.25 * du
          end if

! By correction from Bx
!
          dur = u(ibx,i,jp1,k) - u(ibx,i,j,k)
          dul = u(ibx,i,jm1,k) - u(ibx,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icy,i  ,j,k) = u(icy,i  ,j,k) + 0.25 * du
            u(icy,ip1,j,k) = u(icy,ip1,j,k) - 0.25 * du
          end if

! By correction from Bz
!
          dur = u(ibz,i,jp1,k) - u(ibz,i,j,k)
          dul = u(ibz,i,jm1,k) - u(ibz,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icy,i,j,k  ) = u(icy,i,j,k  ) + 0.25 * du
            u(icy,i,j,kp1) = u(icy,i,j,kp1) - 0.25 * du
          end if

! Bz correction from Bx
!
          dur = u(ibx,i,j,kp1) - u(ibx,i,j,k)
          dul = u(ibx,i,j,km1) - u(ibx,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icz,i  ,j,k) = u(icz,i  ,j,k) + 0.25 * du
            u(icz,ip1,j,k) = u(icz,ip1,j,k) - 0.25 * du
          end if

! Bz correction from By
!
          dur = u(iby,i,j,kp1) - u(iby,i,j,k)
          dul = u(iby,i,j,km1) - u(iby,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icz,i,j  ,k) = u(icz,i,j  ,k) + 0.25 * du
            u(icz,i,jp1,k) = u(icz,i,jp1,k) - 0.25 * du
          end if

        end do
      end do
    end do
#else /* NDIMS == 3 */
    do k = 1, km
      do j = 2, jm - 1
        jm1 = j - 1
        jp1 = j + 1
        do i = 2, im - 1
          im1 = i - 1
          ip1 = i + 1

! Bx correction from By
!
          dur = u(iby,ip1,j,k) - u(iby,i,j,k)
          dul = u(iby,im1,j,k) - u(iby,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icx,i,j  ,k) = u(icx,i,j  ,k) + 0.25 * du
            u(icx,i,jp1,k) = u(icx,i,jp1,k) - 0.25 * du
          end if

! By correction from Bx
!
          dur = u(ibx,i,jp1,k) - u(ibx,i,j,k)
          dul = u(ibx,i,jm1,k) - u(ibx,i,j,k)

          ds  = - dur * dul

          if (ds .gt. 0.0) then
            du = ds / (dur - dul)

            u(icy,i  ,j,k) = u(icy,i  ,j,k) + 0.25 * du
            u(icy,ip1,j,k) = u(icy,ip1,j,k) - 0.25 * du
          end if
        end do
      end do
    end do
#endif /* NDIMS == 3 */
!
  end subroutine magtocen
!
!===============================================================================
!
! divergence: function calculates divergence of the input staggered field
!             [Bx(i+1/2,j,k), By(i,j+1/2,k), Bz(i,j,k+1/2)] -> div B
!
!===============================================================================
!
  subroutine divergence(b, db, dx, dy, dz)

    use blocks, only : nvr
    use config, only : im, jm, km

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
!
!===============================================================================
!
! expand_mag: subroutine for prolongation of the magnetic field with preserving
!             the divergence free condition (Li & Li, 2004, JCoPh, 199, 1)
!
!===============================================================================
!
  subroutine expand_mag(dm, fm, ng, ux, uy, uz, wx, wy, wz)

    implicit none

! input and output variables
!
    integer, dimension(3)                , intent(in)    :: dm, fm
    integer                              , intent(in)    :: ng
    real   , dimension(dm(1),dm(2),dm(3)), intent(in)    :: ux, uy, uz
    real   , dimension(fm(1),fm(2),fm(3)), intent(inout) :: wx, wy, wz

! local variables
!
    integer :: i, j, k, im1, jm1, km1

! temporary arrays
!
    real, dimension(:)      , allocatable :: sx, sy, sz, rx, ry, rz
    real, dimension(:,:,:,:), allocatable :: tx, ty
!
!------------------------------------------------------------------------------
!
    wx(:,:,:) = 0.0
    wy(:,:,:) = 0.0
    wz(:,:,:) = 0.0

#if NDIMS == 2
! allocate dimensional vectors for expansion
!
    allocate(sx(maxval(dm(:))))
    allocate(sy(maxval(dm(:))))
    allocate(sz(maxval(dm(:))))
    allocate(rx(maxval(fm(:))))
    allocate(ry(maxval(fm(:))))
    allocate(rz(maxval(fm(:))))

! prepare temporary arrays for expansion along the X direction
!
    allocate(tx(3,fm(1),dm(2),dm(3)))

! expand By and Bz along the X direction using TVD reconstruction
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        sx(1:dm(1)) = ux(1:dm(1),j,k)
        sy(1:dm(1)) = uy(1:dm(1),j,k)
        sz(1:dm(1)) = uz(1:dm(1),j,k)

        call expand_1d(dm(1), fm(1), ng, sx(1:dm(1)), rx(1:fm(1)), 'c')
        call expand_1d(dm(1), fm(1), ng, sy(1:dm(1)), ry(1:fm(1)), 't')
        call expand_1d(dm(1), fm(1), ng, sz(1:dm(1)), rz(1:fm(1)), 't')

        tx(1,1:fm(1),j,k) = rx(1:fm(1))
        tx(2,1:fm(1),j,k) = ry(1:fm(1))
        tx(3,1:fm(1),j,k) = rz(1:fm(1))
      end do
    end do

! calculate Bx from expanded By and Bz using the divergence condition
!
    do k = 1, dm(3)
      do j = 2, dm(2)
        jm1 = j - 1
        do i = 2, fm(1), 2
          im1 = i - 1

          tx(1,im1,j,k) = tx(1,im1,j,k)                                        &
                        + 0.25 * ((tx(2,i,j  ,k) - tx(2,im1,j  ,k))            &
                                - (tx(2,i,jm1,k) - tx(2,im1,jm1,k)))
        end do
      end do
    end do

! expand Bx and Bz along the Y direction using TVD reconstruction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        sx(1:dm(2)) = tx(1,i,1:dm(2),k)
        sy(1:dm(2)) = tx(2,i,1:dm(2),k)
        sz(1:dm(2)) = tx(3,i,1:dm(2),k)

        call expand_1d(dm(2), fm(2), ng, sx(1:dm(2)), rx(1:fm(2)), 't')
        call expand_1d(dm(2), fm(2), ng, sy(1:dm(2)), ry(1:fm(2)), 'c')
        call expand_1d(dm(2), fm(2), ng, sz(1:dm(2)), rz(1:fm(2)), 't')

        wx(i,1:fm(2),k) = rx(1:fm(2))
        wy(i,1:fm(2),k) = ry(1:fm(2))
        wz(i,1:fm(2),k) = rz(1:fm(2))
      end do
    end do

! calculate By from expanded Bx and Bz using the divergence condition
!
    do k = 1, dm(3)
      do j = 2, fm(2), 2
        jm1 = j - 1
        do i = 2, fm(1)
          im1 = i - 1

          wy(i,jm1,k) = wy(i,jm1,k)                                            &
                      + 0.5 * ((wx(i  ,j,k) - wx(i  ,jm1,k))                   &
                             - (wx(im1,j,k) - wx(im1,jm1,k)))
        end do
      end do
    end do

! deallocate direction one dimensional vector
!
    deallocate(sx)
    deallocate(sy)
    deallocate(sz)
    deallocate(rx)
    deallocate(ry)
    deallocate(rz)

! deallocate temporary arrays
!
    deallocate(tx)
#endif /* NDIMS == 2 */
#if NDIMS == 3
! allocate dimensional vectors for expansion
!
    allocate(sx(maxval(dm(:))))
    allocate(sy(maxval(dm(:))))
    allocate(sz(maxval(dm(:))))
    allocate(rx(maxval(fm(:))))
    allocate(ry(maxval(fm(:))))
    allocate(rz(maxval(fm(:))))

! prepare temporary arrays for expansion along the X direction
!
    allocate(tx(3,fm(1),dm(2),dm(3)))
    allocate(ty(3,fm(1),fm(2),dm(3)))

! expand By and Bz along the X direction using TVD reconstruction
!
    do k = 1, dm(3)
      do j = 1, dm(2)
        sx(1:dm(1)) = ux(1:dm(1),j,k)
        sy(1:dm(1)) = uy(1:dm(1),j,k)
        sz(1:dm(1)) = uz(1:dm(1),j,k)

        call expand_1d(dm(1), fm(1), ng, sx(1:dm(1)), rx(1:fm(1)), 'c')
        call expand_1d(dm(1), fm(1), ng, sy(1:dm(1)), ry(1:fm(1)), 't')
        call expand_1d(dm(1), fm(1), ng, sz(1:dm(1)), rz(1:fm(1)), 't')

        tx(1,1:fm(1),j,k) = rx(1:fm(1))
        tx(2,1:fm(1),j,k) = ry(1:fm(1))
        tx(3,1:fm(1),j,k) = rz(1:fm(1))
      end do
    end do

! calculate Bx from expanded By and Bz using the divergence condition
!
    do k = 2, dm(3)
      km1 = k - 1
      do j = 2, dm(2)
        jm1 = j - 1
        do i = 2, fm(1), 2
          im1 = i - 1

          tx(1,im1,j,k) = tx(1,im1,j,k)                                        &
                        + 0.25 * ((tx(2,i,j  ,k) - tx(2,im1,j  ,k))            &
                                - (tx(2,i,jm1,k) - tx(2,im1,jm1,k))            &
                                + (tx(3,i,j,k  ) - tx(3,im1,j,k  ))            &
                                - (tx(3,i,j,km1) - tx(3,im1,j,km1)))
        end do
      end do
    end do

! expand Bx and Bz along the Y direction using TVD reconstruction
!
    do k = 1, dm(3)
      do i = 1, fm(1)
        sx(1:dm(2)) = tx(1,i,1:dm(2),k)
        sy(1:dm(2)) = tx(2,i,1:dm(2),k)
        sz(1:dm(2)) = tx(3,i,1:dm(2),k)

        call expand_1d(dm(2), fm(2), ng, sx(1:dm(2)), rx(1:fm(2)), 't')
        call expand_1d(dm(2), fm(2), ng, sy(1:dm(2)), ry(1:fm(2)), 'c')
        call expand_1d(dm(2), fm(2), ng, sz(1:dm(2)), rz(1:fm(2)), 't')

        ty(1,i,1:fm(2),k) = rx(1:fm(2))
        ty(2,i,1:fm(2),k) = ry(1:fm(2))
        ty(3,i,1:fm(2),k) = rz(1:fm(2))
      end do
    end do

! calculate Bx from expanded By and Bz using the divergence condition
!
    do k = 2, dm(3)
      km1 = k - 1
      do j = 2, fm(2), 2
        jm1 = j - 1
        do i = 2, dm(1)
          im1 = i - 1

          ty(2,i,jm1,k) = ty(2,i,jm1,k)                                        &
                        + 0.5  * ((ty(1,i  ,j,k) - ty(1,i  ,jm1,k))            &
                                - (ty(1,im1,j,k) - ty(1,im1,jm1,k)))           &
                        + 0.25 * ((ty(3,i,j,k  ) - ty(3,i,jm1,k  ))            &
                                - (ty(3,i,j,km1) - ty(3,i,jm1,km1)))
        end do
      end do
    end do

! expand Bx and Bz along the Y direction using TVD reconstruction
!
    do j = 1, fm(2)
      do i = 1, fm(1)
        sx(1:dm(3)) = tx(1,i,j,1:dm(3))
        sy(1:dm(3)) = tx(2,i,j,1:dm(3))
        sz(1:dm(3)) = tx(3,i,j,1:dm(3))

        call expand_1d(dm(3), fm(3), ng, sx(1:dm(3)), rx(1:fm(3)), 't')
        call expand_1d(dm(3), fm(3), ng, sy(1:dm(3)), ry(1:fm(3)), 't')
        call expand_1d(dm(3), fm(3), ng, sz(1:dm(3)), rz(1:fm(3)), 'c')

        wx(i,j,1:fm(2)) = rx(1:fm(3))
        wy(i,j,1:fm(2)) = ry(1:fm(3))
        wz(i,j,1:fm(2)) = rz(1:fm(3))
      end do
    end do

! calculate Bx from expanded By and Bz using the divergence condition
!
    do k = 2, fm(3), 2
      km1 = k - 1
      do j = 2, fm(2)
        jm1 = j - 1
        do i = 2, fm(1)
          im1 = i - 1

          wz(i,j,km1) = wz(i,j,km1)                                            &
                      + 0.5 * ((wx(i  ,j,k) - wx(i  ,j,km1))                   &
                             - (wx(im1,j,k) - wx(im1,j,km1))                   &
                             + (wy(i  ,j,k) - wy(i  ,j,km1))                   &
                             - (wy(i,jm1,k) - wy(i,jm1,km1)))
        end do
      end do
    end do

! deallocate direction one dimensional vector
!
    deallocate(sx)
    deallocate(sy)
    deallocate(sz)
    deallocate(rx)
    deallocate(ry)
    deallocate(rz)

! deallocate temporary arrays
!
    deallocate(tx)
    deallocate(ty)
#endif /* NDIMS == 3 */
!
  end subroutine expand_mag
!
!===============================================================================
!
! expand_mag_bnd: subroutine for prolongation of the magnetic field with
!                 preserving the divergence free condition and applied
!                 boundaries (Li & Li, 2004, JCoPh, 199, 1)
!
! arguments:
!
!   id        - the direction of the expansion
!   if        - the side of the boundary values
!   dl        - number of safety zones
!   dm(3)     - the dimensions of the expanded region without boundaries
!   bn(:,:)   - the finer values of the boundary
!   ux(:,:,:) - X component of the coarse field
!   uy(:,:,:) - Y component of the coarse field
!   uz(:,:,:) - Z component of the coarse field
!   wx(:,:,:) - X component of the fine field
!   wy(:,:,:) - Y component of the fine field
!   wz(:,:,:) - Z component of the fine field
!
!===============================================================================
!
  subroutine expand_mag_bnd(id, if, dl, dm, bn, ux, uy, uz, wx, wy, wz)

    implicit none

! input and output variables
!
    integer                  , intent(in)    :: id, if, dl
    integer, dimension(3)    , intent(in)    :: dm
    real   , dimension(:,:)  , intent(in)    :: bn
    real   , dimension(:,:,:), intent(in)    :: ux, uy, uz
    real   , dimension(:,:,:), intent(inout) :: wx, wy, wz

! local variables
!
    integer :: i, j, k, im1, jm1, km1, ip1, jp1, kp1
    integer :: ib, ie, jb, je, kb, ke
    integer :: il, jl, kl

! local arrays
!
    integer, dimension(3) :: cm, fm

! temporary arrays
!
    real, dimension(:,:,:), allocatable :: tx, ty
!
!------------------------------------------------------------------------------
!
! notes:
!
!  - the dimensions of parallel and perpendicular components are different;
!    the parallel components is [ dm(1) + 1, dm(2) + 2 * ng, dm(3) + 2 * ng];
!    the perpendicular is [ dm(1) + 2 * ng, dm(2) + 2 * ng + 1, dm(3) + 2 * ng]
!    or [ dm(1) + 2 * ng, dm(2) + 2 * ng, dm(3) + 2 * ng + 1]
!
!  - the boundary is an array of the dimensions
!    [ 4 * dm(2) + 2 * ng, 4 * dm(3) + 2 * ng ]
!
!  - the parallel output vector is of the dimensions
!    [ 2 * dm(1) + 1, 2 * dm(2) + 2 * ng, 2 * dm(2) + 2 * ng ];
!
!  - the perpendicular output vector is of the dimensions
!    [ 2 * dm(1), 2 * dm(2) + 2 * ng + 1, 2 * dm(2) + 2 * ng ] or
!    [ 2 * dm(1), 2 * dm(2) + 2 * ng, 2 * dm(2) + 2 * ng + 1 ];
!
!  - depending on the direction and side of the boundary vector choose the right
!    interpolation;
!
!  - perpendicular components to the direction of boundary should be larger by
!    at least 2 ghost zones in order to handle properly TVD interpolation;
!
!  - boundary array should be two dimensional and have perpendicular dimensions
!    of the fine output array, e.g. if idir = 1, bn(fm(2),fm(3));

! prepare dimensions
!
    cm(:) = dm(:) + 2 * dl
#if NDIMS == 2
    cm(3) = 1
#endif /* NDIMS == 2 */

    fm(:) = 2 * dm(:)
#if NDIMS == 2
    fm(3) = 1
#endif /* NDIMS == 2 */

    wx(:,:,:) = 0.0
    wy(:,:,:) = 0.0
    wz(:,:,:) = 0.0

    select case(id)

!! boundary perpendicular to the X direction
!!
    case(1)

! allocate temporary arrays
!
      allocate(tx(fm(1)+1,cm(2)  ,cm(3)))
      allocate(ty(fm(1)  ,cm(2)+1,cm(3)))

      tx(:,:,:) = 0.0
      ty(:,:,:) = 0.0

! prepare indices
!
      il = 1 + (1 - if) * fm(1)

      jb = 1     + dl
      je = cm(2) - dl

! steps:
!
! 1. check if the input arrays have the right dimensions
!

! 2. interpolate components perpendicular to the direction of boundary using the
!    TVD interpolation; place the calculated values in the output array;
!
      do k = 1, cm(3)
        do j = 1, cm(2) + 1
          call expand_1d_tvd(cm(1), fm(1), dl, uy(1:cm(1),j,k), ty(1:fm(1),j,k))
        end do
      end do

! 3. fill the already known values of the parallel component
!
      tx(1:fm(1)+1:2,1:cm(2),1:cm(3)) = ux(1:dm(1)+1,1:cm(2),1:cm(3))

! 4. average the fine boundary in order to get the coarse values and substitute
!    these values in the coarse parallel component;
!
      tx(il,jb:je,1:cm(3)) = 0.50 * (bn(1:fm(2):2,1:cm(3))                     &
                                   + bn(2:fm(2):2,1:cm(3)))

! 5. calculate the even values of the parallel component using the divergence
!    free condition;
!
      do k = 1, cm(3)
        do j = 1, cm(2)
          jp1 = j + 1
          do i = 2, fm(1) + 1, 2
            im1 = i - 1
            ip1 = i + 1

            tx(i,j,k) = 0.50 * (tx(im1,j  ,k) + tx(ip1,j  ,k))                 &
                      + 0.25 * (ty(i  ,jp1,k) - ty(im1,jp1,k))                 &
                      - 0.25 * (ty(i  ,j  ,k) - ty(im1,j  ,k))
          end do
        end do
      end do

! 6. interpolate components perpendicular to the direction of boundary using the
!    TVD interpolation; place the calculated values in the output array;
!
      do k = 1, cm(3)
        do i = 1, fm(1) + 1
          call expand_1d_tvd(cm(2), fm(2), dl, tx(i,1:cm(2),k), wx(i,1:fm(2),k))
        end do
      end do

! 7. substitute the original values of the fine boundary in the parallel
!    component;
!
      wx(il,1:fm(2),1:fm(3)) = bn(1:fm(2),1:fm(3))

! 8. substitute the odd values of the perpendicular component from the previous
!    array;
!
      wy(1:fm(1),1:fm(2)+1:2,1:fm(3)) = ty(1:fm(1),jb:je+1,1:fm(3))

! 9. calculate the even values of the perpendicular component using the
!    divergence free condition;
!
      do k = 1, fm(3)
        do j = 2, fm(2), 2
          jm1 = j - 1
          jp1 = j + 1
          do i = 1, fm(1)
            ip1 = i + 1

            wy(i,j,k) = 0.50 * (wy(i  ,jm1,k) + wy(i  ,jp1,k))                 &
                      + 0.50 * (wx(ip1,j  ,k) - wx(ip1,jm1,k))                 &
                      - 0.50 * (wx(i  ,j  ,k) - wx(i  ,jm1,k))
          end do
        end do
      end do

! deallocate temporary arrays
!
      deallocate(tx)
      deallocate(ty)

!! boundary perpendicular to the Y direction
!!
    case(2)

! allocate temporary arrays
!
      allocate(tx(cm(1)+1,fm(2)  ,cm(3)))
      allocate(ty(cm(1)  ,fm(2)+1,cm(3)))

      tx(:,:,:) = 0.0
      ty(:,:,:) = 0.0

! prepare indices
!
      jl = 1 + (1 - if) * fm(2)

      ib = 1     + dl
      ie = cm(1) - dl

! steps:
!
! 1. check if the input arrays have the right dimensions
!

! 2. interpolate components perpendicular to the direction of boundary using the
!    TVD interpolation; place the calculated values in the output array;
!
      do k = 1, cm(3)
        do i = 1, cm(1) + 1
          call expand_1d_tvd(cm(2), fm(2), dl, ux(i,1:cm(2),k), tx(i,1:fm(2),k))
        end do
      end do

! 3. fill the already known values of the parallel component
!
      ty(1:cm(1),1:fm(2)+1:2,1:cm(3)) = uy(1:cm(1),1:dm(2)+1,1:cm(3))

! 4. average the fine boundary in order to get the coarse values and substitute
!    these values in the coarse parallel component;
!
      ty(ib:ie,jl,1:cm(3)) = 0.50 * (bn(1:fm(1):2,1:cm(3))                     &
                                   + bn(2:fm(1):2,1:cm(3)))

! 5. calculate the even values of the parallel component using the divergence
!    free condition;
!
      do k = 1, cm(3)
        do j = 2, fm(2) + 1, 2
          jm1 = j - 1
          jp1 = j + 1
          do i = 1, cm(1)
            ip1 = i + 1

            ty(i,j,k) = 0.50 * (ty(i  ,jm1,k) + ty(i  ,jp1,k))                 &
                      + 0.25 * (tx(ip1,j  ,k) - tx(ip1,jm1,k))                 &
                      - 0.25 * (tx(i  ,j  ,k) - tx(i  ,jm1,k))
          end do
        end do
      end do

! 6. interpolate components perpendicular to the direction of boundary using the
!    TVD interpolation; place the calculated values in the output array;
!
      do k = 1, cm(3)
        do j = 1, fm(2) + 1
          call expand_1d_tvd(cm(1), fm(1), dl, ty(1:cm(1),j,k), wy(1:fm(1),j,k))
        end do
      end do

! 7. substitute the original values of the fine boundary in the parallel
!    component;
!
      wy(1:fm(1),jl,1:fm(3)) = bn(1:fm(1),1:fm(3))

! 8. substitute the odd values of the perpendicular component from the previous
!    array;
!
      wx(1:fm(1)+1:2,1:fm(2),1:fm(3)) = tx(ib:ie+1,1:fm(2),1:fm(3))

! 9. calculate the even values of the perpendicular component using the
!    divergence free condition;
!
      do k = 1, fm(3)
        do j = 1, fm(2)
          jp1 = j + 1
          do i = 2, fm(1), 2
            im1 = i - 1
            ip1 = i + 1

            wx(i,j,k) = 0.50 * (wx(im1,j  ,k) + wx(ip1,j  ,k))                 &
                      + 0.50 * (wy(i  ,jp1,k) - wy(im1,jp1,k))                 &
                      - 0.50 * (wy(i  ,j  ,k) - wy(im1,j  ,k))
          end do
        end do
      end do


! deallocate temporary arrays
!
      deallocate(tx)
      deallocate(ty)

    case default
    end select
!
  end subroutine expand_mag_bnd
#endif /* FLUXCT */
!
!===============================================================================
!
end module
