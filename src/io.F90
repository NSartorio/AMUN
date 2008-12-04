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

    use blocks, only : block, pfirst
    use error , only : print_error
    use hdf5  , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f         &
                     , h5gcreate_f, h5gclose_f, h5acreate_f, h5aclose_f     &
                     , h5awrite_f, h5screate_simple_f, h5sclose_f           &
                     , hid_t, hsize_t, H5F_ACC_TRUNC_F                      &
                     , H5T_NATIVE_INTEGER, H5T_NATIVE_DOUBLE

    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc

! HDF5 variables
!
    integer(hid_t)    :: fid, gid, sid, aid
    integer(hsize_t)  :: am(1)

! local variables
!
    character(len=64) :: fl, gnm
    integer           :: err

! pointers
!
    type(block), pointer :: pcurr, ptemp
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
        pcurr => pfirst

        do while(associated(pcurr))

          write(gnm,"('blk',i8.8)") pcurr%id

          call h5gcreate_f(fid, gnm, gid, err)

! TODO: create block attributes (xmin, xmax, ..., level, parent, all entries of structure)
!
          am(1) = 1
          call h5screate_simple_f(1, am, sid, err)

          call h5acreate_f(gid, 'id', H5T_NATIVE_INTEGER, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, pcurr%id, am, err)
          call h5aclose_f(aid, err)

          ptemp => pcurr%prev
          call h5acreate_f(gid, 'prev', H5T_NATIVE_INTEGER, sid, aid, err)
          if (associated(ptemp)) then
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, ptemp%id, am, err)
          else
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, -1, am, err)
          endif
          call h5aclose_f(aid, err)

          ptemp => pcurr%next
          call h5acreate_f(gid, 'next', H5T_NATIVE_INTEGER, sid, aid, err)
          if (associated(ptemp)) then
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, ptemp%id, am, err)
          else
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, -1, am, err)
          endif
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'parent', H5T_NATIVE_INTEGER, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, pcurr%parent, am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'level', H5T_NATIVE_INTEGER, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_INTEGER, pcurr%level, am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'xmin', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%xmin,8), am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'xmax', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%xmax,8), am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'ymin', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%ymin,8), am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'ymax', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%ymax,8), am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'zmin', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%zmin,8), am, err)
          call h5aclose_f(aid, err)

          call h5acreate_f(gid, 'zmax', H5T_NATIVE_DOUBLE, sid, aid, err)
          call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pcurr%zmax,8), am, err)
          call h5aclose_f(aid, err)

          call h5sclose_f(sid, err)

          call h5gclose_f(gid, err)

          pcurr => pcurr.next
        end do

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
