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
!   iintd - the number of steps between subsequent intervals storing;
!
  integer(kind=4), save :: funit = 7
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

! close the integrals file
!
    if (master) close(funit)

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
    use coordinates    , only : ib, ie, jb, je, kb, ke
    use coordinates    , only : advol
    use equations      , only : idn, imx, imy, imz, ien, ibx, iby, ibz
    use evolution      , only : step, time, dt
    use mpitools       , only : master
#ifdef MPI
    use mpitools       , only : reduce_sum_real_array
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
    real(kind=8), dimension(narr) :: arr
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
    arr(:) = 0.0d+00

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
      arr(1) = arr(1) + sum(pdata%u(idn,ib:ie,jb:je,kb:ke)) * dvol
      arr(2) = arr(2) + sum(pdata%u(imx,ib:ie,jb:je,kb:ke)) * dvol
      arr(3) = arr(3) + sum(pdata%u(imy,ib:ie,jb:je,kb:ke)) * dvol
      arr(4) = arr(4) + sum(pdata%u(imz,ib:ie,jb:je,kb:ke)) * dvol

! sum up total energy
!
      if (ien > 0) then
        arr(5) = arr(5) + sum(pdata%u(ien,ib:ie,jb:je,kb:ke)) * dvol
      end if

! sum up kinetic energy
!
      arr(6) = arr(6) + sum((pdata%u(imx,ib:ie,jb:je,kb:ke)**2                 &
                           + pdata%u(imy,ib:ie,jb:je,kb:ke)**2                 &
                           + pdata%u(imz,ib:ie,jb:je,kb:ke)**2)                &
                           / pdata%u(idn,ib:ie,jb:je,kb:ke)) * dvolh

! sum up magnetic energy
!
      if (ibx > 0) then
        arr(7) = arr(7) + sum(pdata%u(ibx,ib:ie,jb:je,kb:ke)**2                &
                            + pdata%u(iby,ib:ie,jb:je,kb:ke)**2                &
                            + pdata%u(ibz,ib:ie,jb:je,kb:ke)**2) * dvolh
      end if

! associate the pointer with the next block on the list
!
      pdata => pdata%next

    end do ! data blocks

#ifdef MPI
! sum the integral array from all processes
!
    call reduce_sum_real_array(narr, arr(:), iret)
#endif /* MPI */

! calculate the internal energy
!
    if (ien > 0) arr(8) = arr(5) - arr(6) - arr(7)

! write down the integrals to the integrals file
!
    if (master) write(funit,"(i9,10(1x,1e18.8))") step, time, dt, arr(1:8)

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
