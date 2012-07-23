!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2012 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: INTEGRALS - handles calculation of the integrals such as total mass,
!!                     momenta, energies, etc., and stores them in a file
!!
!!******************************************************************************
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
  subroutine init_integrals(flag)

    use mpitools , only : master

    implicit none

! input arguments
!
    logical, intent(in) :: flag

! local arguments
!
    logical :: info
!
!-------------------------------------------------------------------------------
!
! only master process creates the file
!
    if (master) then

! check if the integrals file exists
!
      inquire(file = fname, exist = info)

! if flag is set to .true. or the original integrals file does not exist, create
! a new one and write the header in it, otherwise open the original file and
! move the writing position to the end to allow for appending
!
      if (flag .or. .not. info) then

! create a new integrals file
!
        open(unit=funit, file=fname, form='formatted', status='replace')

! write the integral file header
!
        write(funit,"('#')")
        write(funit,"('#',a6,12(1x,a15))") 'step', 'time', 'dt', 'mass'        &
                                         , 'momx', 'momy', 'momz'              &
                                         , 'ener', 'ekin', 'emag', 'eint'      &
                                         , 'fcor', 'finp'
        write(funit,"('#')")

      else

! open the integrals file and set the position at the end of file
!
        open(unit=funit, file=fname, form='formatted', position='append')

! write a marker that the job has been restarted from here
!
        write(funit,"('#',1x,a)") "job restarted from this point"

      end if

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

    use mpitools, only : master

    implicit none
!
!-------------------------------------------------------------------------------
!
! close integrals.dat
!
    if (master) close(funit)
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
    use config   , only : ib, ie, jb, je, kb, ke
    use evolution, only : n, t, dt
#ifdef FORCE
    use forcing  , only : fcor, finp
#endif /* FORCE */
    use coordinates, only : advol
    use mpitools , only : master
#ifdef MPI
    use mpitools , only : reduce_sum_real_array
#endif /* MPI */
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
    integer, parameter :: narr = 16

! local variables
!
    integer         :: iret
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

#ifdef FORCE
! update the force-force correlation
!
    arr(9) = fcor

#endif /* FORCE */
#ifdef MPI
! sum the integrals from all processors
!
    call reduce_sum_real_array(narr, arr(:), iret)

#endif /* MPI */
#ifdef ADI
! calculate the internal energy
!
    arr(8) = arr(5) - arr(6) - arr(7)

#endif /* ADI */
#ifdef FORCE
! update the force-force correlation
!
    arr(10) = finp

#endif /* FORCE */
! close integrals.dat
!
    if (master) then
      write(funit,"(i8,12(1x,1pe15.8))") n, t, dt, arr(1:10)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine store_integrals

!===============================================================================
!
end module
