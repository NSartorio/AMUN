!!******************************************************************************
!!
!! module: forcing - handles turbulence driving
!!
!! Copyright (C) 2007-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module forcing

  implicit none

#ifdef FORCE
! number of driven components
!
  integer                             , save :: nf

! array of k vectors, mode amplitudes, and unit vectors
!
  integer, dimension(:,:), allocatable, save :: ktab
  complex, dimension(:,:), allocatable, save :: ftab
  real   , dimension(:,:), allocatable, save :: e1, e2
#endif /* FORCE */

  contains
!
!===============================================================================
!
! init_forcing: subroutine allocates and initializes the forcing variables
!               required e.g. for driving turbulence
!
!===============================================================================
!
  subroutine init_forcing()

#ifdef FORCE
    use config, only : fpow, fani, fdt, kf, kl, ku, kc, kd
#endif /* FORCE */

    implicit none

#ifdef FORCE
! local variables
!
    integer :: kmx, i, j, k, l
    real    :: rk, kr, fnor, famp, fa, kx, ky, kz, kxy, kyz
#endif /* FORCE */

!-------------------------------------------------------------------------------
!
#ifdef FORCE
! initialize the number of drived components, normalization factor, and
! maximum wave number
!
    nf   = 0
    fnor = 0.0d0
    kmx  = int(ku + 1)

! calculate the number of drived components and normalization factor
!
#if NDIMS == 2
    do j = - kmx, kmx, kd
      do i = - kmx, kmx, kd
        rk = sqrt(real(i * i + j * j))
        if (rk .ge. kl .and. rk .le. ku) then
          nf   = nf + 1
          kr   = (rk - kf) / kc
          fnor = fnor + exp(-0.5d0 * kr * kr)
        end if
      end do
    end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = - kmx, kmx, kd
      do j = - kmx, kmx, kd
        do i = - kmx, kmx, kd
          rk = sqrt(real(i * i + j * j + k * k))
          if (rk .ge. kl .and. rk .le. ku) then
            nf   = nf + 1
            kr   = (rk - kf) / kc
            fnor = fnor + exp(-0.5d0 * kr * kr)
          end if
        end do
      end do
    end do
#endif /* NDIMS == 3 */

! calculate the maximum driving amplitude
!
    famp = sqrt(fpow * fdt / fnor)

! allocate arrays for k vectors, mode amplitudes and unit vectors
!
    allocate(ktab(nf,3))
    allocate(ftab(nf,3))
    allocate(e1  (nf,3))
    allocate(e2  (nf,3))

! prepare k vector, amplitude and unit vectors for each node
!
    l = 0
#if NDIMS == 2
    do j = - kmx, kmx, kd
      do i = - kmx, kmx, kd
        rk = sqrt(real(i * i + j * j))
        if (rk .ge. kl .and. rk .le. ku) then
          l    = l + 1

! prepare k vector
!
          ktab(l,1) = i
          ktab(l,2) = j
          ktab(l,3) = 0

! compute its amplitude
!
          kr  = (rk - kf) / kc
          fa  = famp * exp(-0.5d0 * kr * kr)

! prepare the unit vectors
!
          kx  = real(i)
          ky  = real(j)
          kxy = sqrt(kx * kx + ky * ky)

          e1(l,1) =   ky / kxy
          e1(l,2) = - kx / kxy
          e1(l,3) = 0.0d0

          e2(l,1) = 0.0d0
          e2(l,2) = 0.0d0
          e2(l,3) = 1.0d0

          e1(l,:) = fa * e1(l,:)
          e2(l,:) = fa * e2(l,:)
        end if
      end do
    end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do k = - kmx, kmx, kd
      do j = - kmx, kmx, kd
        do i = - kmx, kmx, kd
          rk = sqrt(real(i * i + j * j + k * k))
          if (rk .ge. kl .and. rk .le. ku) then
            l    = l + 1

! prepare k vector
!
            ktab(l,1) = i
            ktab(l,2) = j
            ktab(l,3) = k

! compute its amplitude
!
            kr = (rk - kf) / kc
            fa = famp * exp(-0.5d0 * kr * kr)

! prepare the unit vectors
!
            kx  = real(i)
            ky  = real(j)
            kz  = real(k)
            kxy = sqrt(kx * kx + ky * ky)

            if (kxy .gt. 0.0d0) then
              e1(l,1) =   ky / kxy
              e1(l,2) = - kx / kxy
              e1(l,3) = 0.0d0

              e2(l,1) =   kx * kz / (rk * kxy)
              e2(l,2) =   ky * kz / (rk * kxy)
              e2(l,3) = - kxy / rk
            else
              kyz = sqrt(ky * ky + kz * kz)
              e1(l,1) = 0.0d0
              e1(l,2) =   kz / kyz
              e1(l,3) = - ky / kyz

              e2(l,1) = - kyz / rk
              e2(l,2) =   ky * kx / (rk * kyz)
              e2(l,3) =   kz * kx / (rk * kyz)
            end if

            e1(l,:) = fa * e1(l,:)
            e2(l,:) = fa * e2(l,:)
          end if
        end do
      end do
    end do
#endif /* NDIMS == 3 */
#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
  end subroutine init_forcing
!
!===============================================================================
!
! clear_forcing: subroutine deallocates the forcing module allocatable arrays
!
!===============================================================================
!
  subroutine clear_forcing()

    implicit none

!-------------------------------------------------------------------------------
!
#ifdef FORCE
    if (allocated(ktab)) deallocate(ktab)
    if (allocated(ftab)) deallocate(ftab)
    if (allocated(e1))   deallocate(e1)
    if (allocated(e2))   deallocate(e2)
#endif /* FORCE */
!
!-------------------------------------------------------------------------------
!
  end subroutine clear_forcing
#ifdef FORCE
!
!===============================================================================
!
! evolve_forcing: subroutine evolves the forcing terms in Fourier space
!
!===============================================================================
!
  subroutine evolve_forcing(dt)

    use config   , only : fdt
    use constants, only : dpi
    use random   , only : randomu

    implicit none

! input arguments
!
    real, intent(in) :: dt

! local variables
!
    integer :: l, n, ni
    complex :: aran, bran
    real    :: phi, th1, th2

!-------------------------------------------------------------------------------
!
! calculate the number of forcing integration iteration for the current timestep
!
    ni = int(dt / fdt)

! iterate over all perturbed Fourier components
!
    do l = 1, nf

! reset complex wave coefficient
!
      aran = cmplx(0.0d0, 0.0d0)
      bran = cmplx(0.0d0, 0.0d0)

! integrate over all forcing substeps
!
      do n = 1, ni

! obtain random phases phi, th1, and th2
!
        phi = dpi * randomu(0)
        th1 = dpi * randomu(0)
        th2 = dpi * randomu(0)

! update coefficients for the current k-vector
!
        aran = aran + sin(phi) * cmplx(cos(th1), sin(th1))
        bran = bran + cos(phi) * cmplx(cos(th2), sin(th2))

      end do

! update the Fourier coefficients for the current k-vector
!
      ftab(l,:) = aran * e1(l,:) + bran * e2(l,:)

    end do

! normalize the forcing terms by the time step
!
    ftab(:,:) = ftab(:,:) / dt
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_forcing
!
!===============================================================================
!
! real_forcing: subroutine returns the forcing terms in transformed to the real
!               space for a given position and level
!
!===============================================================================
!
  subroutine real_forcing(l, xmn, ymn, zmn, f)

    use config   , only : im, jm, km
    use constants, only : dpi
    use mesh     , only : ax, ay, az

    implicit none

! input/output arguments
!
    integer                    , intent(in)    :: l
    real                       , intent(in)    :: xmn, ymn, zmn
    real, dimension(3,im,jm,km), intent(inout) :: f

! local variables
!
    integer :: i, j, k, p
    real    :: fx, fy, fz, kr, cs, sn

! local arrays
!
    real, dimension(im) :: x
    real, dimension(jm) :: y
#if NDIMS == 3
    real, dimension(km) :: z
#endif /* NDIMS == 3 */

!-------------------------------------------------------------------------------
!
! prepare local block coordinates
!
    x(:) = dpi * (xmn + ax(l,:))
    y(:) = dpi * (ymn + ay(l,:))
#if NDIMS == 3
    z(:) = dpi * (zmn + az(l,:))
#endif /* NDIMS == 3 */

! perform inverse Fourier transform
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          fx = 0.0d0
          fy = 0.0d0
          fz = 0.0d0

          do p = 1, nf
#if NDIMS == 2
            kr = ktab(p,1) * x(i) + ktab(p,2) * y(j)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            kr = ktab(p,1) * x(i) + ktab(p,2) * y(j) + ktab(p,3) * z(k)
#endif /* NDIMS == 3 */

            cs = cos(kr)
            sn = sin(kr)

            fx = fx + real(ftab(p,1)) * cs + aimag(ftab(p,1)) * sn
            fy = fy + real(ftab(p,2)) * cs + aimag(ftab(p,2)) * sn
            fz = fz + real(ftab(p,3)) * cs + aimag(ftab(p,3)) * sn
          end do

! update coefficients
!
          f(1,i,j,k) = fx
          f(2,i,j,k) = fy
          f(3,i,j,k) = fz
        end do
      end do
    end do
!
!-------------------------------------------------------------------------------
!
  end subroutine real_forcing
#endif /* FORCE */
!===============================================================================
!
end module forcing
