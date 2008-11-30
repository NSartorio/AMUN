!!*****************************************************************************
!!
!! module: io - handling data input/output
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
module io

  implicit none

  contains
#ifdef HDF5
!
!======================================================================
!
! write_data: subroutine writes all data to file for a given time step
!
!======================================================================
!
  subroutine write_data(ftype, nfile, nproc)

    use error , only : print_error
    use hdf5  , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f  &
                     , hid_t, H5F_ACC_TRUNC_F

    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc

! HDF5 variables
!
    integer(hid_t)    :: fid


! local variables
!
    character(len=64) :: fl
    integer           :: err
!
!----------------------------------------------------------------------
!
! initialize FORTRAN interface
!
    call h5open_f(err)

    if (err .ge. 0) then

! prepare filename
!
      write (fl,'(a1,i6.6,"_",i5.5,a3)') ftype, nfile, nproc, '.h5'

! create file
!
      call h5fcreate_f(fl, H5F_ACC_TRUNC_F, fid, err)

      if (err .ge. 0) then

! TODO: write all attributes
! TODO: write coordinates for all refinement levels
! TODO: iterate over all blocks and write complete structure of each of them
!

! terminate access to the file
!
        call h5fclose_f(fid, err)

        if (err .gt. 0) then

! print error of closing the file
!
          call print_error("io::write_data", "Could not close HDF5 file!")

        endif

      else

! print error of creating HDF5 file
!
        call print_error("io::write_data", "Could not create HDF5 file!")

      endif

! close FORTRAN interface
!
      call h5close_f(err)

      if (err .gt. 0) then

! print error of closing HDF5 Fortran interface
!
        call print_error("io::write_data", "Could not close HDF5 Fortran interface!")

      endif

    else

! print error of initializing HDF5 Fortran interface
!
      call print_error("io::write_data", "Could not initialize HDF5 Fortran interface!")

    endif

!----------------------------------------------------------------------
!
  end subroutine write_data
#endif /* HDF5 */
!======================================================================
!
end module
