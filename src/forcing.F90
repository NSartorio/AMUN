!!******************************************************************************
!!
!! module: forcing - handles turbulence driving
!!
!! Copyright (C) 2007-2011 Grzegorz Kowal <grzegorz@gkowal.info>
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

! force-force correlation (energy injected) in real and Fourier spaces
!
  real(kind=8)                        , save :: fcor, finp

! array of k vectors, mode amplitudes, and unit vectors
!
  integer, dimension(:,:), allocatable, save :: ktab
  complex, dimension(:,:), allocatable, save :: vtab
  complex, dimension(:,:), allocatable, save :: ftab
  real   , dimension(:  ), allocatable, save :: famp
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
    use config   , only : fpow, fani, fdt, kf, kl, ku, kc, kd
    use constants, only : pi
    use mpitools , only : is_master
    use timer    , only : start_timer, stop_timer
#endif /* FORCE */

    implicit none

#ifdef FORCE
! local variables
!
    integer :: kmx, jbg, kbg, i, j, k, l
    real    :: kl2, ku2, kc2
    real    :: rk2, rk, kr2
    real    :: fnor, fa
    real    :: kx, ky, kz, kxy, kyz
#endif /* FORCE */

!-------------------------------------------------------------------------------
!
#ifdef FORCE
! start the timer for forcing
!
    call start_timer(6)

! initialize the number of drived components, normalization factor, and
! maximum wave number
!
    nf   = 0
    fnor = 0.0d0
    fcor = 0.0d0
    finp = 0.0d0
    kmx  = int(ku + 1)
    kl2  = kl * kl
    ku2  = ku * ku
    kc2  = kc * kc

! calculate the number of drived components and normalization factor
!
#if NDIMS == 2
    do i = 0, kmx, kd
      if (i .eq. 0) then
        jbg = 0
      else
        jbg = - kmx
      end if
      do j = jbg, kmx, kd

        rk2 = real(i * i + j * j)

        if (rk2 .ge. kl2 .and. rk2 .le. ku2) then

! increase the number of forcing components
!
          nf   = nf + 1

! update the normalization coefficient
!
          rk   = sqrt(rk2)
          kr2  = (rk - kf)**2 / kc2
          fnor = fnor + exp(- kr2)

        end if

      end do
    end do

#endif /* NDIMS == 2 */
#if NDIMS == 3
    do i = 0, kmx, kd
      if (i .eq. 0) then
        jbg = 0
      else
        jbg = - kmx
      end if
      do j = jbg, kmx, kd
        if (i .eq. 0 .and. j .eq. 0) then
          kbg = 0
        else
          kbg = - kmx
        end if
        do k = kbg, kmx, kd

          rk2 = real(i * i + j * j + k * k)

          if (rk2 .ge. kl2 .and. rk2 .le. ku2) then

! increase the number of forcing components
!
            nf   = nf + 1

! update the normalization coefficient
!
            rk   = sqrt(rk2)
            kr2  = (rk - kf)**2 / kc2
            fnor = fnor + exp(- kr2)

          end if

        end do
      end do
    end do

#endif /* NDIMS == 3 */
! calculate the maximum driving amplitude
!
    fa = sqrt((4.0d0 / 3.0d0) * pi * fpow / fnor)

! print information about forcing
!
    if (is_master()) then
      write (*,'(a)')       ''
      write (*,'(a,a,a)'     ) 'INFO    : ', ' Forcing:'
      write (*,'(a,a,1i7)'   ) 'INFO    : ', '    no. of components   : ', nf
      write (*,'(a,a,1f12.6)') 'INFO    : ', '    normalization coeff.: ', fnor
      write (*,'(a,a,1f12.6)') 'INFO    : ', '    power               : ', fpow
      write (*,'(a,a,1f12.6)') 'INFO    : ', '    amplitude           : ', fa
    endif

! multiply aplitude by the square of forcing time step so we don't have to
! multiply this at each step
!
    fa = fa * sqrt(fdt)

! allocate arrays for k vectors, mode amplitudes and unit vectors
!
    allocate(ktab(nf,3))
    allocate(vtab(nf,3))
    allocate(ftab(nf,3))
    allocate(famp(nf)  )
    allocate(e1  (nf,3))
    allocate(e2  (nf,3))

! initialize the velocity fourier components
!
    vtab(:,:) = cmplx(0.0d0, 0.0d0)

! prepare k vector, amplitude and unit vectors for each node
!
    l = 0
#if NDIMS == 2
    do i = 0, kmx, kd
      if (i .eq. 0) then
        jbg = 0
      else
        jbg = - kmx
      end if
      do j = jbg, kmx, kd

        rk2 = real(i * i + j * j)

        if (rk2 .ge. kl2 .and. rk2 .le. ku2) then

! increase the number of forcing component
!
          l = l + 1

! prepare k vector
!
          ktab(l,1) = i
          ktab(l,2) = j
          ktab(l,3) = 0

! compute its amplitude
!
          rk      = sqrt(rk2)
          kr2     = (rk - kf)**2 / kc2
          famp(l) = fa * exp(- 0.5d0 * kr2)

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
        end if

      end do
    end do
#endif /* NDIMS == 2 */
#if NDIMS == 3
    do i = 0, kmx, kd
      if (i .eq. 0) then
        jbg = 0
      else
        jbg = - kmx
      end if
      do j = jbg, kmx, kd
        if (i .eq. 0 .and. j .eq. 0) then
          kbg = 0
        else
          kbg = - kmx
        end if
        do k = kbg, kmx, kd

          rk2 = real(i * i + j * j + k * k)

          if (rk2 .ge. kl2 .and. rk2 .le. ku2) then

! increase the number of forcing component
!
            l = l + 1

! prepare k vector
!
            ktab(l,1) = i
            ktab(l,2) = j
            ktab(l,3) = k

! compute its amplitude
!
            rk      = sqrt(rk2)
            kr2     = (rk - kf)**2 / kc2
            famp(l) = fa * exp(- 0.5d0 * kr2)

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
          end if

        end do
      end do
    end do
#endif /* NDIMS == 3 */

! stop the timer
!
    call stop_timer(6)
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

#ifdef FORCE
    use timer , only : start_timer, stop_timer

#endif /* FORCE */
    implicit none

!-------------------------------------------------------------------------------
!
#ifdef FORCE
! start the timer for forcing
!
    call start_timer(6)

! deallocate all module arrays
!
    if (allocated(ktab)) deallocate(ktab)
    if (allocated(vtab)) deallocate(vtab)
    if (allocated(ftab)) deallocate(ftab)
    if (allocated(famp)) deallocate(famp)
    if (allocated(e1))   deallocate(e1)
    if (allocated(e2))   deallocate(e2)

! stop the timer
!
    call stop_timer(6)
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

    use config   , only : fpow, fdt
    use constants, only : dpi
#ifdef MPI
    use mpitools , only : mallreducesumc
#endif /* MPI */
    use random   , only : randomu
    use timer    , only : start_timer, stop_timer

    implicit none

! input arguments
!
    real, intent(in) :: dt

! local variables
!
    integer :: l, n, ni
    complex :: aran, bran, xi1, xi2
    real    :: phi, psi, th1, th2, div, tanth1, fpw, fac

!-------------------------------------------------------------------------------
!
! start the timer for forcing
!
    call start_timer(6)

#ifdef MPI
! reduce velocity fourier components from all processors
!
    call mallreducesumc(nf, 3, vtab(:,:))

#endif /* MPI */
! calculate the number of forcing integration iteration for the current timestep
!
    ni = int(dt / fdt)

! iterate over all perturbed Fourier components
!
    do l = 1, nf

! project velocity vector onto the vectors e1 and e2
!
      xi1 = dot_product(vtab(l,:), e1(l,:))
      xi2 = dot_product(vtab(l,:), e2(l,:))

! reset complex wave coefficient
!
      aran = cmplx(0.0d0, 0.0d0)
      bran = cmplx(0.0d0, 0.0d0)

! integrate over all forcing substeps
!
      do n = 1, ni

! obtain random phases phi and psi
!
        phi = dpi * randomu(0)
        psi = dpi * randomu(0)

! obtain phases th1 and th2 from the condition for minimizing
! the velocity-force correlation
!
        div    = - sin(phi) * aimag(xi1)                                       &
                 + cos(phi) * (sin(psi) * real(xi2) - cos(psi) * aimag(xi2))

        if (div .ne. 0.0) then
          tanth1 = (sin(phi) * real(xi1)                                       &
                  + cos(phi) * (sin(psi) * aimag(xi2) + cos(psi) * real(xi2))) &
                 / div
        else
          tanth1 = 0.0
        end if

        th1 = atan(tanth1)
        th2 = psi + th1

! update coefficients for the current k-vector
!
        aran = aran + sin(phi) * cmplx(cos(th1), sin(th1))
        bran = bran + cos(phi) * cmplx(cos(th2), sin(th2))

      end do

! update the Fourier coefficients for the current k-vector
!
      ftab(l,:) = famp(l) * (aran * e1(l,:) + bran * e2(l,:))

    end do

! reset the velocity fourier components
!
    vtab(:,:) = cmplx(0.0, 0.0)

! normalize amplitudes to the correct input power
!
    fpw       = abs(sum(ftab(:,:) * conjg(ftab(:,:))))
    fac       = sqrt(4.0d0 * dt * fpow / fpw)
    ftab(:,:) = fac * ftab(:,:)

! calculate the force-force correlation in Fourier space
!
    fpw  = 0.25d0 * abs(sum(ftab(:,:) * conjg(ftab(:,:))))
    finp = finp + fpw

! stop the timer
!
    call stop_timer(6)
!
!-------------------------------------------------------------------------------
!
  end subroutine evolve_forcing
!
!===============================================================================
!
! fourier_transform: subroutine transforms the velocity to the Fourier space;
!                    only components corresponding to forcing are calculated;
!
!===============================================================================
!
  subroutine fourier_transform(l, xmn, ymn, zmn, u)

    use config   , only : im, jm, km, ib, ie, jb, je, kb, ke
    use constants, only : dpi
    use mesh     , only : ax, ay, az, advol
    use timer    , only : start_timer, stop_timer
    use variables, only : idn, imx, imy, imz

    implicit none

! input/output arguments
!
    integer                    , intent(in)    :: l
    real                       , intent(in)    :: xmn, ymn, zmn
    real, dimension(4,im,jm,km), intent(inout) :: u

! local variables
!
    integer :: i, j, k, p, kmn, kmx
    real    :: kx, ky, kz
    real    :: vx, vy, vz
    real    :: snx, sny, snz, snp, sn
    real    :: csx, csy, csz, csp, cs

! local arrays
!
    real, dimension(im) :: x
    real, dimension(jm) :: y
#if NDIMS == 3
    real, dimension(km) :: z
#endif /* NDIMS == 3 */
    real, dimension(:,:), allocatable :: asnx, asny, asnz
    real, dimension(:,:), allocatable :: acsx, acsy, acsz

!-------------------------------------------------------------------------------
!
! start the timer for forcing
!
    call start_timer(6)

! prepare local block coordinates
!
    x(:) = dpi * (xmn + ax(l,:))
    y(:) = dpi * (ymn + ay(l,:))
#if NDIMS == 3
    z(:) = dpi * (zmn + az(l,:))
#endif /* NDIMS == 3 */

! allocate arrays for directional sinuses and cosinuses
!
    kmn = minval(ktab(:,:))
    kmx = maxval(ktab(:,:))

    allocate(asnx(kmn:kmx,im))
    allocate(acsx(kmn:kmx,im))
    allocate(asny(kmn:kmx,jm))
    allocate(acsy(kmn:kmx,jm))
#if NDIMS == 3
    allocate(asnz(kmn:kmx,km))
    allocate(acsz(kmn:kmx,km))
#endif /* NDIMS == 3 */

! calculate directional sinuses and cosinuses for each mode
!
    do p = kmn, kmx
      do i = 1, im
        kx        = p * x(i)
        asnx(p,i) = sin(kx)
        acsx(p,i) = cos(kx)
      end do
      do j = 1, jm
        ky        = p * y(j)
        asny(p,j) = sin(ky)
        acsy(p,j) = cos(ky)
      end do
#if NDIMS == 3
      do k = 1, km
        kz        = p * z(k)
        asnz(p,k) = sin(kz)
        acsz(p,k) = cos(kz)
      end do
#endif /* NDIMS == 3 */
    end do

! perform the inverse Fourier transform
!
    do k = kb, ke
      do j = jb, je
        do i = ib, ie

! prepare velocity components at the current position
!
          vx = u(imx,i,j,k) / u(idn,i,j,k)
          vy = u(imy,i,j,k) / u(idn,i,j,k)
          vz = u(imz,i,j,k) / u(idn,i,j,k)

! iterate over all forcing components
!
          do p = 1, nf

! obtain directional sinuses and cosinuses for each mode
!
            snx = asnx(ktab(p,1),i)
            csx = acsx(ktab(p,1),i)
            sny = asny(ktab(p,2),j)
            csy = acsy(ktab(p,2),j)
#if NDIMS == 3
            snz = asnz(ktab(p,3),k)
            csz = acsz(ktab(p,3),k)
#endif /* NDIMS == 3 */

! calculate total sinus and cosinus
!
#if NDIMS == 2
            sn  = (snx * csy + csx * sny) * advol(l)
            cs  = (csx * csy - snx * sny) * advol(l)
#endif /* NDIMS == 2 */
#if NDIMS == 3
            snp = snx * csy + csx * sny
            csp = csx * csy - snx * sny

            sn  = (snp * csz + csp * snz) * advol(l)
            cs  = (csp * csz - snp * snz) * advol(l)
#endif /* NDIMS == 3 */

! update the Fourier coefficient
!
            vtab(p,1) = vtab(p,1) + cmplx(vx * cs, vx * sn)
            vtab(p,2) = vtab(p,2) + cmplx(vy * cs, vy * sn)
            vtab(p,3) = vtab(p,3) + cmplx(vz * cs, vz * sn)

          end do

        end do
      end do
    end do

! deallocate local arrays
!
    deallocate(asnx)
    deallocate(acsx)
    deallocate(asny)
    deallocate(acsy)
#if NDIMS == 3
    deallocate(asnz)
    deallocate(acsz)
#endif /* NDIMS == 3 */

! stop the timer
!
    call stop_timer(6)
!
!-------------------------------------------------------------------------------
!
  end subroutine fourier_transform
!
!===============================================================================
!
! real_forcing: subroutine returns the forcing terms in transformed to the real
!               space for a given position and level
!
!===============================================================================
!
  subroutine real_forcing(l, xmn, ymn, zmn, f)

    use config   , only : im, jm, km, ib, ie, jb, je, kb, ke
    use constants, only : dpi
    use mesh     , only : ax, ay, az, advol
    use timer    , only : start_timer, stop_timer

    implicit none

! input/output arguments
!
    integer                    , intent(in)    :: l
    real                       , intent(in)    :: xmn, ymn, zmn
    real, dimension(3,im,jm,km), intent(inout) :: f

! local variables
!
    integer :: i, j, k, p, kmn, kmx
    real    :: fx, fy, fz
    real    :: kx, ky, kz
    real    :: snx, sny, snz, snp, sn
    real    :: csx, csy, csz, csp, cs

! local arrays
!
    real, dimension(im) :: x
    real, dimension(jm) :: y
#if NDIMS == 3
    real, dimension(km) :: z
#endif /* NDIMS == 3 */
    real, dimension(:,:), allocatable :: asnx, asny, asnz
    real, dimension(:,:), allocatable :: acsx, acsy, acsz

!-------------------------------------------------------------------------------
!
! start the timer for forcing
!
    call start_timer(6)

! prepare local block coordinates
!
    x(:) = dpi * (xmn + ax(l,:))
    y(:) = dpi * (ymn + ay(l,:))
#if NDIMS == 3
    z(:) = dpi * (zmn + az(l,:))
#endif /* NDIMS == 3 */

! allocate arrays for directional sinuses and cosinuses
!
    kmn = minval(ktab(:,:))
    kmx = maxval(ktab(:,:))

    allocate(asnx(kmn:kmx,im))
    allocate(acsx(kmn:kmx,im))
    allocate(asny(kmn:kmx,jm))
    allocate(acsy(kmn:kmx,jm))
#if NDIMS == 3
    allocate(asnz(kmn:kmx,km))
    allocate(acsz(kmn:kmx,km))
#endif /* NDIMS == 3 */

! calculate directional sinuses and cosinuses for each mode
!
    do p = kmn, kmx
      do i = 1, im
        kx        = p * x(i)
        asnx(p,i) = sin(kx)
        acsx(p,i) = cos(kx)
      end do
      do j = 1, jm
        ky        = p * y(j)
        asny(p,j) = sin(ky)
        acsy(p,j) = cos(ky)
      end do
#if NDIMS == 3
      do k = 1, km
        kz        = p * z(k)
        asnz(p,k) = sin(kz)
        acsz(p,k) = cos(kz)
      end do
#endif /* NDIMS == 3 */
    end do

! perform the inverse Fourier transform
!
    do k = 1, km
      do j = 1, jm
        do i = 1, im

          fx = 0.0d0
          fy = 0.0d0
          fz = 0.0d0

          do p = 1, nf

! obtain directional sinuses and cosinuses for each mode
!
            snx = asnx(ktab(p,1),i)
            csx = acsx(ktab(p,1),i)
            sny = asny(ktab(p,2),j)
            csy = acsy(ktab(p,2),j)
#if NDIMS == 3
            snz = asnz(ktab(p,3),k)
            csz = acsz(ktab(p,3),k)
#endif /* NDIMS == 3 */

! calculate total sinus and cosinus
!
#if NDIMS == 2
            sn  = snx * csy + csx * sny
            cs  = csx * csy - snx * sny
#endif /* NDIMS == 2 */
#if NDIMS == 3
            snp = snx * csy + csx * sny
            csp = csx * csy - snx * sny

            sn  = snp * csz + csp * snz
            cs  = csp * csz - snp * snz
#endif /* NDIMS == 3 */

! update the real value
!
            fx = fx + real(ftab(p,1)) * cs - aimag(ftab(p,1)) * sn
            fy = fy + real(ftab(p,2)) * cs - aimag(ftab(p,2)) * sn
            fz = fz + real(ftab(p,3)) * cs - aimag(ftab(p,3)) * sn

          end do

! update the local value
!
          f(1,i,j,k) = fx
          f(2,i,j,k) = fy
          f(3,i,j,k) = fz

        end do
      end do
    end do

! calculate the force-force correlation
!
    fcor = fcor + 0.5d0 * sum(f(:,ib:ie,jb:je,kb:ke)**2) * advol(l)

! deallocate local arrays
!
    deallocate(asnx)
    deallocate(acsx)
    deallocate(asny)
    deallocate(acsy)
#if NDIMS == 3
    deallocate(asnz)
    deallocate(acsz)
#endif /* NDIMS == 3 */

! stop the timer
!
    call stop_timer(6)
!
!-------------------------------------------------------------------------------
!
  end subroutine real_forcing
#endif /* FORCE */
!===============================================================================
!
end module forcing
