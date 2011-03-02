!!******************************************************************************
!!
!! module: integrals - handles calculation of the integrals such as total mass,
!!                     momenta, energies, etc., and stores them in a file
!!
!! Copyright (C) 2011 Grzegorz Kowal <grzegorz@gkowal.info>
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
module integrals

  implicit none

  character(len=32) :: fname = "integrals.dat"
  integer           :: funit = 10

  contains
!
!===============================================================================
!
! init_integrals: subroutine initializes the integral variables to zero and
!                 prepares the file where the integrals are stored
!
!===============================================================================
!
  subroutine init_integrals()

    use mpitools , only : is_master

    implicit none
!
!-------------------------------------------------------------------------------
!
! only master process creates the file
!
    if (is_master()) then

! create a new integrals.dat file
!
      open(unit=funit, file=fname, form='formatted', status='replace')

! write the integral file header
!
      write(funit,"('#')")
      write(funit,"('#',a6,10(1x,a15))") 'step', 'time', 'dt', 'mass'          &
                                       , 'momx', 'momy', 'momz'                &
                                       , 'ener', 'ekin', 'emag', 'eint'
      write(funit,"('#')")

    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine init_integrals
!
!===============================================================================
!
! clear_integrals: subroutine closes the integrals file and frees all allocated
!                  module variables
!
!===============================================================================
!
  subroutine clear_integrals()

    use mpitools, only : is_master

    implicit none
!
!-------------------------------------------------------------------------------
!
! close integrals.dat
!
    if (is_master()) close(funit)
!
!-------------------------------------------------------------------------------
!
  end subroutine clear_integrals
!
!===============================================================================
!
! store_integrals: subroutine calculates the integrals and stores them in the
!                  integrals file
!
!===============================================================================
!
  subroutine store_integrals()

    use blocks   , only : block_meta, block_data, list_data
    use blocks   , only : dblocks, ndims
    use config   , only : ib, ie, jb, je, kb, ke
    use evolution, only : n, t, dt
    use mesh     , only : advol
    use mpitools , only : is_master, mallreducesumd
    use variables, only : idn, imx, imy, imz
#ifdef ADI
    use variables, only : ien
#endif /* ADI */
#ifdef MHD
    use variables, only : ibx, iby, ibz
#endif /* MHD */

    implicit none

! local parameters
!
    integer, parameter :: narr = 18

! local variables
!
    real(kind=PREC) :: dvol

! local arrays
!
    real(kind=8), dimension(narr) :: arr

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! reset the total mass
!
    arr(:) = 0.0d0

! iterate over all data blocks and sum up the density
!
    pdata => list_data
    do while(associated(pdata))

      dvol = advol(pdata%meta%level)

      arr(1) = arr(1) + sum(pdata%u(idn,ib:ie,jb:je,kb:ke)) * dvol
      arr(2) = arr(2) + sum(pdata%u(imx,ib:ie,jb:je,kb:ke)) * dvol
      arr(3) = arr(3) + sum(pdata%u(imy,ib:ie,jb:je,kb:ke)) * dvol
      arr(4) = arr(4) + sum(pdata%u(imz,ib:ie,jb:je,kb:ke)) * dvol
#ifdef ADI
      arr(5) = arr(5) + sum(pdata%u(ien,ib:ie,jb:je,kb:ke)) * dvol
#endif /* ADI */
      arr(6) = arr(6) + sum((pdata%u(imx,ib:ie,jb:je,kb:ke)**2                 &
                           + pdata%u(imy,ib:ie,jb:je,kb:ke)**2                 &
                           + pdata%u(imz,ib:ie,jb:je,kb:ke)**2)                &
                           / pdata%u(idn,ib:ie,jb:je,kb:ke)) * 0.5d0 * dvol
#ifdef MHD
      arr(7) = arr(7) + sum(pdata%u(ibx,ib:ie,jb:je,kb:ke)**2                  &
                          + pdata%u(iby,ib:ie,jb:je,kb:ke)**2                  &
                          + pdata%u(ibz,ib:ie,jb:je,kb:ke)**2) * 0.5d0 * dvol
#endif /* MHD */

      pdata => pdata%next
    end do

#ifdef MPI
! sum the integrals from all processors
!
    call mallreducesumd(narr, arr(:))

#endif /* MPI */
#ifdef ADI
! calculate the internal energy
!
    arr(8) = arr(5) - arr(6) - arr(7)

#endif /* ADI */
! close integrals.dat
!
    if (is_master()) then
      write(funit,"(i8,10(1x,1pe15.8))") n, t, dt, arr(1:8)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine store_integrals

!===============================================================================
!
end module
