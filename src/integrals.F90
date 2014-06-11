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
!! module: INTEGRALS
!!
!!  This module provides subroutines to calculate and store integrals of
!!  the conserved variables, and other useful statistics in the integrals file.
!!
!!******************************************************************************
!
module integrals

#ifdef PROFILE
! import external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer            , save :: imi, ims
#endif /* PROFILE */

! MODULE PARAMETERS:
! =================
!
!   funit - a file handler to the integrals file;
!   sunit - a file handler to the statistics file;
!   iintd - the number of steps between subsequent intervals storing;
!
  integer(kind=4), save :: funit = 7
  integer(kind=4), save :: sunit = 8
  integer(kind=4), save :: iintd = 1

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_integrals, finalize_integrals
  public :: store_integrals

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
!===============================================================================
!
! subroutine INITIALIZE_INTEGRALS:
! -------------------------------
!
!   Subroutine initializes module INTEGRALS by preparing the integrals file.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     irun    - job execution counter;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_integrals(verbose, irun, iret)

! import external variables and subroutines
!
    use mpitools       , only : master
    use parameters     , only : get_parameter_integer

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: irun, iret

! local variables
!
    character(len=32)      :: fname
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('integrals:: initialization'   , imi)
    call set_timer('integrals:: integrals storing', ims)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the integrals storing interval
!
    call get_parameter_integer("integrals_interval", iintd)

! make sure storage interval is larger than zero
!
    iintd = max(1, iintd)

! only master process handles the integral file
!
    if (master) then

! generate the integrals file name
!
      write(fname, "('integrals_',i2.2,'.dat')") irun

! create a new integrals file
!
#ifdef INTEL
      open (newunit = funit, file = fname, form = 'formatted'                  &
                                       , status = 'replace', buffered = 'yes')
#else /* INTEL */
      open (newunit = funit, file = fname, form = 'formatted'                  &
                                       , status = 'replace')
#endif /* INTEL */

! write the integral file header
!
      write(funit,"('#',a8,10(1x,a18))") 'step', 'time', 'dt'                  &
                                       , 'mass', 'momx', 'momy', 'momz'        &
                                       , 'ener', 'ekin', 'emag', 'eint'
      write(funit,"('#')")

! generate the statistics file name
!
      write(fname, "('statistics_',i2.2,'.dat')") irun

! create a new statistics file
!
#ifdef INTEL
      open (newunit = sunit, file = fname, form = 'formatted'                  &
                                       , status = 'replace', buffered = 'yes')
#else /* INTEL */
      open (newunit = sunit, file = fname, form = 'formatted'                  &
                                       , status = 'replace')
#endif /* INTEL */

! write the integral file header
!
      write(sunit,"('#',a8,23(1x,a18))") 'step', 'time'                        &
                                 , 'mean(dens)', 'min(dens)', 'max(dens)'      &
                                 , 'mean(pres)', 'min(pres)', 'max(pres)'      &
                                 , 'mean(velo)', 'min(velo)', 'max(velo)'      &
                                 , 'mean(magn)', 'min(magn)', 'max(magn)'      &
                                 , 'mean(bpsi)', 'min(bpsi)', 'max(bpsi)'      &
                                 , 'mean(Mson)', 'min(Mson)', 'max(Mson)'      &
                                 , 'mean(Malf)', 'min(Malf)', 'max(Malf)'
      write(sunit,"('#')")

    end if ! master

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_integrals
!
!===============================================================================
!
! subroutine FINALIZE_INTEGRALS:
! -----------------------------
!
!   Subroutine finalizes module INTEGRALS by closing the integrals file.
!
!
!===============================================================================
!
  subroutine finalize_integrals()

! import external variables
!
    use mpitools, only : master

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! close the integrals and statistics files
!
    if (master) then
      close(funit)
      close(sunit)
    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_integrals
!
!===============================================================================
!
! subroutine STORE_INTEGRALS:
! --------------------------
!
!   Subroutine calculates the integrals, collects from all processes and
!   stores them in the integrals file.
!
!
!===============================================================================
!
  subroutine store_integrals()

! import external variables and subroutines
!
    use blocks         , only : block_meta, block_data, list_data
    use coordinates    , only : in, jn, kn, ib, jb, kb, ie, je, ke
    use coordinates    , only : advol, voli
    use equations      , only : idn, ipr, ivx, ivy, ivz, ibx, iby, ibz, ibp
    use equations      , only : ien, imx, imy, imz
    use equations      , only : gamma, csnd
    use evolution      , only : step, time, dt
    use mpitools       , only : master
#ifdef MPI
    use mpitools       , only : reduce_sum_real_array
    use mpitools       , only : reduce_minimum_real_array
    use mpitools       , only : reduce_maximum_real_array
#endif /* MPI */

! local variables are not implicit by default
!
    implicit none

! local variables
!
    integer                       :: iret
    real(kind=8)                  :: dvol, dvolh

! local pointers
!
    type(block_data), pointer     :: pdata

! local parameters
!
    integer, parameter            :: narr = 16

! local arrays
!
    real(kind=8), dimension(narr)     :: inarr, avarr, mnarr, mxarr
    real(kind=8), dimension(in,jn,kn) :: vel, mag, sqd, tmp

! parameters
!
    real(kind=8), parameter :: eps = epsilon(1.0d+00)
    real(kind=8), parameter :: big = huge(1.0d+00)
!
!-------------------------------------------------------------------------------
!
! return if the storage interval was not reached
!
    if (mod(step, iintd) > 0) return

#ifdef PROFILE
! start accounting time for the integrals storing
!
    call start_timer(ims)
#endif /* PROFILE */

! reset the integrals array
!
    inarr(:) =   0.0d+00
    avarr(:) =   0.0d+00
    mnarr(:) =   big
    mxarr(:) = - big

! reset some statistics if they are not used
!
    if (ipr < 1) then
      mnarr(2) = 0.0d+00
      mxarr(2) = 0.0d+00
    end if
    if (ibx < 1) then
      mnarr(4) = 0.0d+00
      mxarr(4) = 0.0d+00
      mnarr(5) = 0.0d+00
      mxarr(5) = 0.0d+00
      mnarr(7) = 0.0d+00
      mxarr(7) = 0.0d+00
    end if

! associate the pointer with the first block on the data block list
!
    pdata => list_data

! iterate over all data blocks on the list
!
    do while(associated(pdata))

! obtain the volume elements for the current block
!
      dvol  = advol(pdata%meta%level)
      dvolh = 0.5d+00 * dvol

! sum up density and momenta components
!
      inarr(1) = inarr(1) + sum(pdata%u(idn,ib:ie,jb:je,kb:ke)) * dvol
      inarr(2) = inarr(2) + sum(pdata%u(imx,ib:ie,jb:je,kb:ke)) * dvol
      inarr(3) = inarr(3) + sum(pdata%u(imy,ib:ie,jb:je,kb:ke)) * dvol
      inarr(4) = inarr(4) + sum(pdata%u(imz,ib:ie,jb:je,kb:ke)) * dvol

! sum up total energy
!
      if (ien > 0) then
        inarr(5) = inarr(5) + sum(pdata%u(ien,ib:ie,jb:je,kb:ke)) * dvol
      end if

! sum up kinetic energy
!
      inarr(6) = inarr(6) + sum((pdata%u(imx,ib:ie,jb:je,kb:ke)**2             &
                               + pdata%u(imy,ib:ie,jb:je,kb:ke)**2             &
                               + pdata%u(imz,ib:ie,jb:je,kb:ke)**2)            &
                               / pdata%u(idn,ib:ie,jb:je,kb:ke)) * dvolh

! sum up magnetic energy
!
      if (ibx > 0) then
        inarr(7) = inarr(7) + sum(pdata%u(ibx,ib:ie,jb:je,kb:ke)**2            &
                                + pdata%u(iby,ib:ie,jb:je,kb:ke)**2            &
                                + pdata%u(ibz,ib:ie,jb:je,kb:ke)**2) * dvolh
      end if

! get average, minimum and maximum values of density
!
      tmp(:,:,:) = pdata%q(idn,ib:ie,jb:je,kb:ke)
      avarr(1) = avarr(1) + sum(tmp(:,:,:)) * dvol
      mnarr(1) = min(mnarr(1), minval(tmp(:,:,:)))
      mxarr(1) = max(mxarr(1), maxval(tmp(:,:,:)))

! get average, minimum and maximum values of pressure
!
      if (ipr > 0) then
        tmp(:,:,:) = pdata%q(ipr,ib:ie,jb:je,kb:ke)
        avarr(2) = avarr(2) + sum(tmp(:,:,:)) * dvol
        mnarr(2) = min(mnarr(2), minval(tmp(:,:,:)))
        mxarr(2) = max(mxarr(2), maxval(tmp(:,:,:)))
      end if

! get average, minimum and maximum values of velocity amplitude
!
      vel(:,:,:) = sqrt(sum(pdata%q(ivx:ivz,ib:ie,jb:je,kb:ke)**2, 1))
      avarr(3) = avarr(3) + sum(vel(:,:,:)) * dvol
      mnarr(3) = min(mnarr(3), minval(vel(:,:,:)))
      mxarr(3) = max(mxarr(3), maxval(vel(:,:,:)))

! get average, minimum and maximum values of magnetic field amplitude, and
! divergence potential
!
      if (ibx > 0) then
        mag(:,:,:) = sqrt(sum(pdata%q(ibx:ibz,ib:ie,jb:je,kb:ke)**2, 1))
        avarr(4) = avarr(4) + sum(mag(:,:,:)) * dvol
        mnarr(4) = min(mnarr(4), minval(mag(:,:,:)))
        mxarr(4) = max(mxarr(4), maxval(mag(:,:,:)))

        tmp(:,:,:) = pdata%q(ibp,ib:ie,jb:je,kb:ke)
        avarr(5) = avarr(5) + sum(tmp(:,:,:)) * dvol
        mnarr(5) = min(mnarr(5), minval(tmp(:,:,:)))
        mxarr(5) = max(mxarr(5), maxval(tmp(:,:,:)))
      end if

! calculate square root of density
!
      sqd(:,:,:) = sqrt(pdata%q(idn,ib:ie,jb:je,kb:ke))

! get average, minimum and maximum values of sonic Mach number
!
      if (ipr > 0) then
        tmp(:,:,:) = sqd(:,:,:) * vel(:,:,:)                                   &
                                / sqrt(gamma * pdata%q(ipr,ib:ie,jb:je,kb:ke))
      else
        tmp(:,:,:) = vel(:,:,:) / csnd
      end if
      avarr(6) = avarr(6) + sum(tmp(:,:,:)) * dvol
      mnarr(6) = min(mnarr(6), minval(tmp(:,:,:)))
      mxarr(6) = max(mxarr(6), maxval(tmp(:,:,:)))

! get average, minimum and maximum values of Alfvénic Mach number
!
      if (ibx > 0) then
        tmp(:,:,:) = sqd(:,:,:) * vel(:,:,:) / max(eps, mag(:,:,:))
        avarr(7) = avarr(7) + sum(tmp(:,:,:)) * dvol
        mnarr(7) = min(mnarr(7), minval(tmp(:,:,:)))
        mxarr(7) = max(mxarr(7), maxval(tmp(:,:,:)))
      end if

! associate the pointer with the next block on the list
!
      pdata => pdata%next

    end do ! data blocks

#ifdef MPI
! sum the integral array from all processes
!
    call reduce_sum_real_array(narr, inarr(:), iret)

! reduce average, minimum and maximum values
!
    call reduce_sum_real_array(narr, avarr(:), iret)
    call reduce_minimum_real_array(narr, mnarr(:), iret)
    call reduce_maximum_real_array(narr, mxarr(:), iret)
#endif /* MPI */

! calculate the internal energy
!
    if (ien > 0) inarr(8) = inarr(5) - inarr(6) - inarr(7)

! normalize the averages by the volume of domain
!
    avarr(:) = avarr(:) * voli

! write down the integrals and statistics to appropriate files
!
    if (master) then
      write(funit,"(i9,10(1x,1e18.8))") step, time, dt, inarr(1:8)
      write(sunit,"(i9,23(1x,1e18.8))") step, time                             &
                                            , avarr(1), mnarr(1), mxarr(1)     &
                                            , avarr(2), mnarr(2), mxarr(2)     &
                                            , avarr(3), mnarr(3), mxarr(3)     &
                                            , avarr(4), mnarr(4), mxarr(4)     &
                                            , avarr(5), mnarr(5), mxarr(5)     &
                                            , avarr(6), mnarr(6), mxarr(6)     &
                                            , avarr(7), mnarr(7), mxarr(7)
    end if

#ifdef PROFILE
! stop accounting time for the integrals storing
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine store_integrals

!===============================================================================
!
end module
