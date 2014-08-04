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
!! module: IO
!!
!!  This module handles data storage and job restart from restart files.
!!
!!
!!******************************************************************************
!
module io

! import external subroutines
!
  use blocks, only : pointer_meta
#ifdef PROFILE
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

! subroutine interfaces
!
  interface write_attribute
#ifdef HDF5
    module procedure write_scalar_attribute_integer_h5
    module procedure write_scalar_attribute_double_h5
    module procedure write_vector_attribute_integer_h5
    module procedure write_vector_attribute_double_h5
#endif /* HDF5 */
  end interface
  interface read_attribute
#ifdef HDF5
    module procedure read_scalar_attribute_integer_h5
    module procedure read_scalar_attribute_double_h5
    module procedure read_vector_attribute_integer_h5
    module procedure read_vector_attribute_double_h5
#endif /* HDF5 */
  end interface
  interface write_array
#ifdef HDF5
    module procedure write_1d_array_integer_h5
    module procedure write_2d_array_integer_h5
    module procedure write_3d_array_integer_h5
    module procedure write_4d_array_integer_h5
    module procedure write_5d_array_integer_h5
    module procedure write_1d_array_double_h5
    module procedure write_2d_array_double_h5
    module procedure write_3d_array_double_h5
    module procedure write_4d_array_double_h5
    module procedure write_5d_array_double_h5
#endif /* HDF5 */
  end interface
  interface read_array
#ifdef HDF5
    module procedure read_1d_array_integer_h5
    module procedure read_2d_array_integer_h5
    module procedure read_3d_array_integer_h5
    module procedure read_4d_array_integer_h5
    module procedure read_5d_array_integer_h5
    module procedure read_1d_array_double_h5
    module procedure read_2d_array_double_h5
    module procedure read_3d_array_double_h5
    module procedure read_4d_array_double_h5
    module procedure read_5d_array_double_h5
#endif /* HDF5 */
  end interface

#ifdef PROFILE
! timer indices
!
  integer            , save :: ioi, iow, ios
#endif /* PROFILE */

! MODULE PARAMETERS:
! =================
!
!   respath - the directory from which the restart snapshots should be read;
!   ftype   - the type of snapshots to write:
!        'p' -> all primitive variables (default);
!        'c' -> all conserved variables;
!   nrest   - for job restarting, this is the number of restart snapshot;
!   irest   - the local counter for the restart snapshots;
!   isnap   - the local counter for the regular snapshots;
!   ishift  - the shift of the snapshot counter for restarting job with
!             different snapshot interval;
!   hrest   - the execution time interval for restart snapshot storing
!             (in hours); the minimum allowed value is 3 minutes;
!   hsnap   - the problem time interval for regular snapshot storing;
!   tsnap   - the next snapshot time;
!
  character(len=255), save :: respath     = "./"
  character         , save :: ftype       = "p"
  integer           , save :: nrest       = -1
  integer(kind=4)   , save :: irest       =  1
  integer(kind=4)   , save :: isnap       =  0
  integer(kind=4)   , save :: ishift      =  0
  real(kind=8)      , save :: hrest       =  6.0e+00
  real(kind=8)      , save :: hsnap       =  1.0e+00
  real(kind=8)      , save :: tsnap       =  0.0e+00

! flags to determine the way of data writing
!
  logical           , save :: with_ghosts = .true.

! local variables to store the number of processors and maximum level read from
! the restart file
!
  integer(kind=4)   , save :: rtoplev = 1, nfiles = 1

! the coefficient related to the difference between the maximum level stored in
! the restart file and set through the configuration file
!
  integer(kind=4)   , save :: ucor = 1, dcor = 1

! array of pointer used during job restart
!
  type(pointer_meta), dimension(:), allocatable, save :: block_array

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_io
  public :: read_restart_snapshot, write_restart_snapshot, write_snapshot
  public :: restart_from_snapshot
  public :: next_tout

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
! subroutine INITIALIZE_IO:
! ------------------------
!
!   Subroutine initializes module IO by setting its parameters.
!
!   Arguments:
!
!     verbose - flag determining if the subroutine should be verbose;
!     irun    - job execution counter;
!     iret    - return flag of the procedure execution status;
!
!===============================================================================
!
  subroutine initialize_io(verbose, irun, iret)

! import external procedures
!
    use parameters     , only : get_parameter_integer, get_parameter_real      &
                              , get_parameter_string

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: irun, iret

! local variables
!
    character(len=255) :: ghosts = "on"
    integer            :: dd, hh, mm, ss
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('io:: initialization'  , ioi)
    call set_timer('io:: snapshot writing', iow)
    call set_timer('io:: snapshot reading', ios)

! start accounting time for module initialization/finalization
!
    call start_timer(ioi)
#endif /* PROFILE */

! get restart parameters
!
    call get_parameter_string ("restart_path"     , respath)
    call get_parameter_integer("restart_number"   , nrest  )
    call get_parameter_real   ("restart_interval" , hrest  )

! get the interval between snapshots
!
    call get_parameter_string ("snapshot_type"    , ftype  )
    call get_parameter_real   ("snapshot_interval", hsnap  )

! get the flag determining if the ghost cells are stored
!
    call get_parameter_string ("include_ghosts"   , ghosts )

! check ghost cell storing flag
!
    select case(trim(ghosts))
    case ("off", "OFF", "n", "N", "false", "FALSE", "no", "NO")
      with_ghosts = .false.
    case default
      with_ghosts = .true.
    end select

! return the run number
!
    irun = max(1, nrest)

! print info about snapshot parameters
!
    if (verbose) then
      if (ftype == 'p') write (*,"(4x,a13,10x,'=',1x,a)")                      &
                                     "snapshot type", "primitive variables"
      if (ftype == 'c') write (*,"(4x,a13,10x,'=',1x,a)")                      &
                                     "snapshot type", "conservative variables"
      if (with_ghosts) then
        write (*,"(4x,a21,2x,'=',1x,a)") "with ghosts cells    ", "on"
      else
        write (*,"(4x,a21,2x,'=',1x,a)") "with ghosts cells    ", "off"
      end if
      write (*,"(4x,a21,2x,'=',1x,e9.2)") "snapshot interval    ", hsnap
      if (hrest > 0.0d+00) then
        dd = int(hrest / 2.4d+01)
        hh = int(mod(hrest, 2.4d+01))
        mm = int(mod(6.0d+01 * hrest, 6.0d+01))
        ss = int(mod(3.6d+03 * hrest, 6.0d+01))
        write (*,"(4x,a16,7x,'=',1x,i2.2,'d',i2.2,'h',i2.2,'m',i2.2,'s')")     &
                                            "restart interval", dd, hh, mm, ss
      end if
      if (restart_from_snapshot()) then
        write (*,"(4x,a18,5x,'=',1x,'[',a,']')") "restart from path ", trim(respath)
        write (*,"(4x,a21,2x,'=',1x,i4)") "restart from snapshot", nrest
      end if
    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(ioi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_io
!
!===============================================================================
!
! function RESTART_FROM_SNAPSHOT:
! ------------------------------
!
!   Subroutine returns true if the job was selected to be restarted from
!   a snapshot.
!
!
!===============================================================================
!
  logical function restart_from_snapshot()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
    restart_from_snapshot = (nrest > 0)

!-------------------------------------------------------------------------------
!
  end function restart_from_snapshot
!
!===============================================================================
!
! subroutine READ_RESTART_SNAPSHOT:
! --------------------------------
!
!   Subroutine reads restart snapshot files in order to resume the job.
!   This is a wrapper calling specific format subroutine.
!
!
!===============================================================================
!
  subroutine read_restart_snapshot()

! import external variables
!
    use evolution      , only : time

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the data reading
!
    call start_timer(ios)
#endif /* PROFILE */

#ifdef HDF5
! read HDF5 restart file and rebuild the meta and data block structures
!
    call read_restart_snapshot_h5()
#endif /* HDF5 */

! calculate the shift of the snapshot counter, and the next snapshot time
!
    ishift = int(time / hsnap) - isnap + 1
    tsnap  = (ishift + isnap) * hsnap

#ifdef PROFILE
! stop accounting time for the data reading
!
    call stop_timer(ios)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine read_restart_snapshot
!
!===============================================================================
!
! subroutine WRITE_RESTART_SNAPSHOTS:
! ----------------------------------
!
!   Subroutine stores current restart snapshot files.  This is a wrapper
!   calling specific format subroutine.
!
!   Arguments:
!
!     thrs - the current execution time in hours;
!     nrun - the run number;
!     iret - the return flag;
!
!===============================================================================
!
  subroutine write_restart_snapshot(thrs, nrun, iret)

! local variables are not implicit by default
!
    implicit none

! input and output arguments
!
    real(kind=8), intent(in)  :: thrs
    integer     , intent(in)  :: nrun
    integer     , intent(out) :: iret
!
!-------------------------------------------------------------------------------
!
! check if conditions for storing the restart snapshot have been met
!
    if (hrest < 5.0e-02 .or. thrs < irest * hrest) return

#ifdef PROFILE
! start accounting time for the data writing
!
    call start_timer(iow)
#endif /* PROFILE */

#ifdef HDF5
! store restart file
!
    call write_restart_snapshot_h5(nrun, iret)
#endif /* HDF5 */

! increase the restart snapshot counter
!
    irest = irest + 1

#ifdef PROFILE
! stop accounting time for the data writing
!
    call stop_timer(iow)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine write_restart_snapshot
!
!===============================================================================
!
! subroutine WRITE_SNAPSHOT:
! -------------------------
!
!   Subroutine stores block data in snapshots.  Block variables are grouped
!   together and stored in big 4D arrays separately.  This is a wrapper for
!   specific format storing.
!
!
!===============================================================================
!
  subroutine write_snapshot()

! import external variables
!
    use evolution      , only : time

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
! check if conditions for storing the regular snapshot have been met
!
    if (hsnap <= 0.0e+00 .or. time < tsnap) return

#ifdef PROFILE
! start accounting time for the data writing
!
    call start_timer(iow)
#endif /* PROFILE */

#ifdef HDF5
! store variable snapshot file
!
    call write_snapshot_h5()
#endif /* HDF5 */

! increase the snapshot counter and calculate the next snapshot time
!
    isnap = isnap + 1
    tsnap = (ishift + isnap) * hsnap

#ifdef PROFILE
! stop accounting time for the data writing
!
    call stop_timer(iow)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine write_snapshot
!
!===============================================================================
!
! function NEXT_TOUT:
! ------------------
!
!   Function returns the next data snapshot time.
!
!
!===============================================================================
!
  real(kind=8) function next_tout()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
    if (hsnap > 0.0d+00) then
      next_tout = tsnap
    else
      next_tout = huge(hsnap)
    end if

!-------------------------------------------------------------------------------
!
  end function next_tout
!
!===============================================================================
!!
!!***  PRIVATE SUBROUTINES  ****************************************************
!!
!===============================================================================
!
#ifdef HDF5
!===============================================================================
!
! subroutine READ_RESTART_SNAPSHOT_H5:
! -----------------------------------
!
!   Subroutine reads restart snapshot, i.e. parameters, meta and data blocks
!   stored in the HDF5 format restart files and reconstructs the data structure
!   in order to resume a terminated job.
!
!
!===============================================================================
!
  subroutine read_restart_snapshot_h5()

! import external procedures and variables
!
    use blocks         , only : change_blocks_process
    use error          , only : print_error
    use hdf5           , only : hid_t
    use hdf5           , only : H5F_ACC_RDONLY_F
    use hdf5           , only : h5open_f, h5close_f
    use hdf5           , only : h5fis_hdf5_f, h5fopen_f, h5fclose_f
#ifdef MPI
    use mesh           , only : redistribute_blocks
#endif /* MPI */
    use mpitools       , only : nprocs, nproc

! local variables are not implicit by default
!
    implicit none

! local variables
!
    character(len=255) :: fl
    integer(hid_t)     :: fid
    integer            :: err, lfile
    logical            :: info
!
!-------------------------------------------------------------------------------
!
!! 1. RESTORE PARAMETERS AND META BLOCKS FROM THE FIRST FILE
!!
! prepare the filename
!
    write (fl, "(a,'r',i6.6,'_',i5.5,'.h5')") trim(respath), nrest, 0

! check if the HDF5 file exists
!
    inquire(file = fl, exist = info)

! if file does not exist, print error and quit the subroutine
!
    if (.not. info) then
      call print_error("io::read_restart_snapshot_h5"                          &
                                  , "File " // trim(fl) // " does not exist!")
      return
    end if ! file does not exist

! initialize the FORTRAN interface
!
    call h5open_f(err)

! in the case of error, print a message and quit the subroutine
!
    if (err < 0) then
      call print_error("io::read_restart_snapshot_h5"                          &
                            , "Cannot initialize the HDF5 Fortran interface!")
      return
    end if

! check if this file is in the HDF5 format
!
    call h5fis_hdf5_f(fl, info, err)

! if the format verification failed, close the interface, print error and exit
!
    if (err < 0) then
      call print_error("io::read_restart_snapshot_h5"                          &
                                            , "Cannot check the file format!")
      call h5close_f(err)
      return
    end if

! the file is not in the HDF5 format, print message and quit
!
    if (.not. info) then
      call print_error("io::read_restart_snapshot_h5"                          &
                             , "File " // trim(fl) // " is not an HDF5 file!")
      call h5close_f(err)
      return
    end if

! open the HDF5 file
!
    call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! if the file could not be opened, print message and quit
!
    if (err < 0) then
      call print_error("io::read_restart_snapshot_h5"                          &
                                           , "Cannot open file: " // trim(fl))
      call h5close_f(err)
      return
    end if

! read global attributes
!
    call read_attributes_h5(fid)

! read meta blocks and recreate the meta block hierarchy
!
    call read_metablocks_h5(fid)

! close the file
!
    call h5fclose_f(fid, err)

! if the file could not be closed print message and quit
!
    if (err > 0) then
      call print_error("io::read_restart_snapshot_h5"                          &
                                          , "Cannot close file: " // trim(fl))
      call h5close_f(err)
      return
    end if

!! 1. RESTORE DATA BLOCKS
!!
! separate data blocks reading into two cases, when the number of processors is
! larger or equal to the number of files, and when we have less processors than
! files
!
! first, read data blocks to processes which have corresponding restart file
!
    if (nproc < nfiles) then

! prepare the filename
!
      write (fl, "(a,'r',i6.6,'_',i5.5,'.h5')") trim(respath), nrest, nproc

! check if the HDF5 file exists
!
      inquire(file = fl, exist = info)

! if file does not exist, print error and quit the subroutine
!
      if (.not. info) then
        call print_error("io::read_restart_snapshot_h5"                        &
                                  , "File " // trim(fl) // " does not exist!")
        call h5close_f(err)
        return
      end if ! file does not exist

! check if this file is in the HDF5 format
!
      call h5fis_hdf5_f(fl, info, err)

! if the format verification failed, close the interface, print error and exit
!
      if (err < 0) then
        call print_error("io::read_restart_snapshot_h5"                        &
                                            , "Cannot check the file format!")
        call h5close_f(err)
        return
      end if

! the file is not in the HDF5 format, print message and quit
!
      if (.not. info) then
        call print_error("io::read_restart_snapshot_h5"                        &
                             , "File " // trim(fl) // " is not an HDF5 file!")
        call h5close_f(err)
        return
      end if

! open the HDF5 file
!
      call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! if the file could not be opened, print message and quit
!
      if (err < 0) then
        call print_error("io::read_restart_snapshot_h5"                        &
                                           , "Cannot open file: " // trim(fl))
        call h5close_f(err)
        return
      end if

! read data blocks
!
      call read_datablocks_h5(fid)

! close the file
!
      call h5fclose_f(fid, err)

! if the file could not be closed print message and quit
!
      if (err > 0) then
        call print_error("io::read_restart_snapshot_h5"                        &
                                          , "Cannot close file: " // trim(fl))
        call h5close_f(err)
        return
      end if

    end if ! nproc < nfiles

! if there are more files than processes, read the remaining files by
! the last process and redistribute blocks after each processed file,
! otherwise only redistribute blocks
!
    if (nprocs < nfiles) then

! iterate over remaining files and read one by one to the last block
!
      do lfile = nprocs, nfiles - 1

! switch meta blocks from the read file to belong to the reading process
!
        call change_blocks_process(lfile, nprocs - 1)

! read the remaining files by the last process only
!
        if (nproc == nprocs - 1) then

! prepare the filename
!
          write (fl, "(a,'r',i6.6,'_',i5.5,'.h5')") trim(respath), nrest, lfile

! check if the HDF5 file exists
!
          inquire(file = fl, exist = info)

! if file does not exist, print error and quit the subroutine
!
          if (.not. info) then
            call print_error("io::read_restart_snapshot_h5"                    &
                                  , "File " // trim(fl) // " does not exist!")
            call h5close_f(err)
            return
          end if ! file does not exist

! check if this file is in the HDF5 format
!
          call h5fis_hdf5_f(fl, info, err)

! if the format verification failed, close the interface, print error and exit
!
          if (err < 0) then
            call print_error("io::read_restart_snapshot_h5"                    &
                                            , "Cannot check the file format!")
            call h5close_f(err)
            return
          end if

! the file is not in the HDF5 format, print message and quit
!
          if (.not. info) then
            call print_error("io::read_restart_snapshot_h5"                    &
                             , "File " // trim(fl) // " is not an HDF5 file!")
            call h5close_f(err)
            return
          end if

! open the HDF5 file
!
          call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! if the file could not be opened, print message and quit
!
          if (err < 0) then
            call print_error("io::read_restart_snapshot_h5"                    &
                                           , "Cannot open file: " // trim(fl))
            call h5close_f(err)
            return
          end if

! read data blocks
!
          call read_datablocks_h5(fid)

! close the file
!
          call h5fclose_f(fid, err)

! if the file could not be closed print message and quit
!
          if (err > 0) then
            call print_error("io::read_restart_snapshot_h5"                    &
                                          , "Cannot close file: " // trim(fl))
            call h5close_f(err)
            return
          end if

        end if ! nproc == nprocs - 1

#ifdef MPI
! redistribute blocks between processors
!
        call redistribute_blocks()
#endif /* MPI */

      end do ! lfile = nprocs, nfiles - 1

    else ! nprocs < nfiles

#ifdef MPI
! redistribute blocks between processors
!
      call redistribute_blocks()
#endif /* MPI */

    end if ! nprocs < nfiles

! deallocate the array of block pointers
!
    if (allocated(block_array)) deallocate(block_array)

! close the FORTRAN interface
!
    call h5close_f(err)

! check if the interface has been closed successfuly
!
    if (err > 0) then
      call print_error("io::read_restart_snapshot_h5"                          &
                                 , "Cannot close the HDF5 Fortran interface!")
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_restart_snapshot_h5
!
!===============================================================================
!
! subroutine WRITE_RESTART_SNAPSHOT_H5:
! ------------------------------------
!
!   Subroutine writes restart snapshot, i.e. parameters, meta and data blocks
!   to the HDF5 format restart files in order to resume a terminated job later.
!
!   Arguments:
!
!     nrun - the snapshot number;
!     iret - the return flag to inform if subroutine succeeded or failed;
!
!===============================================================================
!
  subroutine write_restart_snapshot_h5(nrun, iret)

! import external procedures and variables
!
    use error          , only : print_error
    use hdf5           , only : hid_t
    use hdf5           , only : H5F_ACC_TRUNC_F
    use hdf5           , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f
    use mpitools       , only : nproc

! local variables are not implicit by default
!
    implicit none

! input and output arguments
!
    integer, intent(in)  :: nrun
    integer, intent(out) :: iret

! local variables
!
    character(len=64) :: fl
    integer(hid_t)    :: fid
    integer           :: err
!
!-------------------------------------------------------------------------------
!
! initialize the FORTRAN interface
!
    call h5open_f(err)

! in the case of error, print a message and quit the subroutine
!
    if (err < 0) then
      call print_error("io::write_restart_snapshot_h5"                         &
                            , "Cannot initialize the HDF5 Fortran interface!")
      iret = 200
      return
    end if

! prepare the restart snapshot filename
!
    write (fl, "('r',i6.6,'_',i5.5,'.h5')") nrun, nproc

! create the new HDF5 file to store the snapshot
!
    call h5fcreate_f(fl, H5F_ACC_TRUNC_F, fid, err)

! if the file could not be created, print message and quit
!
    if (err < 0) then
      call print_error("io::write_restart_snapshot_h5"                         &
                                         , "Cannot create file: " // trim(fl))
      call h5close_f(err)
      iret = 201
      return
    end if

! write the global attributes
!
    call write_attributes_h5(fid)

! write all metablocks which represent the internal structure of domain
!
    call write_metablocks_h5(fid)

! write all datablocks which represent the all variables
!
    call write_datablocks_h5(fid)

! close the file
!
    call h5fclose_f(fid, err)

! if the file could not be closed print message and quit
!
    if (err > 0) then
      call print_error("io::write_restart_snapshot_h5"                         &
                                          , "Cannot close file: " // trim(fl))
      call h5close_f(err)
      iret = 203
      return
    end if

! close the FORTRAN interface
!
    call h5close_f(err)

! check if the interface has been closed successfuly
!
    if (err > 0) then
      call print_error("io::write_restart_snapshot_h5"                         &
                                 , "Cannot close the HDF5 Fortran interface!")
      iret = 204
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_restart_snapshot_h5
!
!===============================================================================
!
! subroutine WRITE_SNAPSHOT_H5:
! ----------------------------
!
!   Subroutine writes the current simulation snapshot, i.e. parameters,
!   coordinates and variables to the HDF5 format files for further processing.
!
!
!===============================================================================
!
  subroutine write_snapshot_h5()

! import external procedures and variables
!
    use error          , only : print_error
    use hdf5           , only : hid_t
    use hdf5           , only : H5F_ACC_TRUNC_F
    use hdf5           , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f
    use mpitools       , only : nproc

! local variables are not implicit by default
!
    implicit none

! local variables
!
    character(len=64) :: fl
    integer(hid_t)    :: fid
    integer           :: err
!
!-------------------------------------------------------------------------------
!
! initialize the FORTRAN interface
!
    call h5open_f(err)

! in the case of error, print a message and quit the subroutine
!
    if (err < 0) then
      call print_error("io::write_snapshot_h5"                                 &
                            , "Cannot initialize the HDF5 Fortran interface!")
      return
    end if

! prepare the restart snapshot filename
!
    write (fl, "(a1,i6.6,'_',i5.5,'.h5')") ftype, isnap, nproc

! create the new HDF5 file to store the snapshot
!
    call h5fcreate_f(fl, H5F_ACC_TRUNC_F, fid, err)

! if the file could not be created, print message and quit
!
    if (err < 0) then
      call print_error("io::write_snapshot_h5"                                 &
                                         , "Cannot create file: " // trim(fl))
      call h5close_f(err)
      return
    end if

! write the global attributes
!
    call write_attributes_h5(fid)

! depending on the selected type of output file write the right groups
!
    select case(ftype)

    case('c')

! write the coordinates (data block bounds, refinement levels, etc.)
!
      call write_coordinates_h5(fid)

! write the variables stored in data blocks (leafs)
!
      call write_conservative_variables_h5(fid)

    case('p')

! write the coordinates (data block bounds, refinement levels, etc.)
!
      call write_coordinates_h5(fid)

! write the variables stored in data blocks (leafs)
!
      call write_primitive_variables_h5(fid)

    case default

! print information about unsupported file format and quit
!
      call print_error("io::write_snapshot_h5", "File type is not suppoerted!")
      call h5fclose_f(fid, err)
      call h5close_f(err)
      return

    end select

! close the file
!
    call h5fclose_f(fid, err)

! if the file could not be closed print message and quit
!
    if (err > 0) then
      call print_error("io::write_snapshot_h5"                                 &
                                          , "Cannot close file: " // trim(fl))
      call h5close_f(err)
      return
    end if

! close the FORTRAN interface
!
    call h5close_f(err)

! check if the interface has been closed successfuly
!
    if (err > 0) then
      call print_error("io::write_snapshot_h5"                                 &
                                 , "Cannot close the HDF5 Fortran interface!")
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_snapshot_h5
!
!===============================================================================
!
! subroutine WRITE_ATTRIBUTES_H5:
! ------------------------------
!
!   Subroutine stores global attributes in the HDF5 file provided by an
!   identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_attributes_h5(fid)

! import external procedures and variables
!
    use blocks         , only : get_mblocks, get_dblocks, get_nleafs
    use blocks         , only : get_last_id
    use coordinates    , only : minlev, maxlev, toplev
    use coordinates    , only : nn, ng, in, jn, kn, ir, jr, kr
    use coordinates    , only : xmin, xmax, ymin, ymax, zmin, zmax
    use error          , only : print_error
    use evolution      , only : step, time, dt, dtn
    use hdf5           , only : hid_t
    use hdf5           , only : h5gcreate_f, h5gclose_f
    use mpitools       , only : nprocs, nproc
    use random         , only : nseeds, get_seeds

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                :: gid
    integer                       :: err

! local allocatable arrays
!
    integer(kind=4), dimension(:), allocatable :: seeds
!
!-------------------------------------------------------------------------------
!
! create a group to store the global attributes
!
    call h5gcreate_f(fid, 'attributes', gid, err)

! check if the group has been created successfuly
!
    if (err < 0) then

! print error about the problem with creating the group
!
      call print_error("io::write_attributes_h5", "Cannot create the group!")

! return from the subroutine
!
      return

    end if

! store the integer attributes
!
    call write_attribute(gid, 'ndims'  , NDIMS        )
    call write_attribute(gid, 'last_id', get_last_id())
    call write_attribute(gid, 'mblocks', get_mblocks())
    call write_attribute(gid, 'dblocks', get_dblocks())
    call write_attribute(gid, 'nleafs' , get_nleafs() )
    call write_attribute(gid, 'ncells' , nn           )
    call write_attribute(gid, 'nghost' , ng           )
    call write_attribute(gid, 'minlev' , minlev       )
    call write_attribute(gid, 'maxlev' , maxlev       )
    call write_attribute(gid, 'toplev' , toplev       )
    call write_attribute(gid, 'nprocs' , nprocs       )
    call write_attribute(gid, 'nproc'  , nproc        )
    call write_attribute(gid, 'nseeds' , nseeds       )
    call write_attribute(gid, 'step'   , step         )
    call write_attribute(gid, 'isnap'  , isnap        )

! store the real attributes
!
    call write_attribute(gid, 'xmin', xmin)
    call write_attribute(gid, 'xmax', xmax)
    call write_attribute(gid, 'ymin', ymin)
    call write_attribute(gid, 'ymax', ymax)
    call write_attribute(gid, 'zmin', zmin)
    call write_attribute(gid, 'zmax', zmax)
    call write_attribute(gid, 'time', time)
    call write_attribute(gid, 'dt'  , dt  )
    call write_attribute(gid, 'dtn' , dtn )

! store the vector attributes
!
    call write_attribute(gid, 'dims' , (/ in, jn, kn /))
    call write_attribute(gid, 'rdims', (/ ir, jr, kr /))

! store random number generator seed values
!
    if (nseeds > 0) then

! allocate space for seeds
!
      allocate(seeds(nseeds))

! get the seed values
!
      call get_seeds(seeds)

! store them in the current group
!
      call write_attribute(gid, 'seeds', seeds(:))

! deallocate seed array
!
      deallocate(seeds)

    end if ! nseeds > 0

! close the group
!
    call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
    if (err < 0) then

! print error about the problem with closing the group
!
      call print_error("io::write_attributes_h5", "Cannot close the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_attributes_h5
!
!===============================================================================
!
! subroutine READ_ATTRIBUTES_H5:
! -----------------------------
!
!   Subroutine restores global attributes from an HDF5 file provided by its
!   identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine read_attributes_h5(fid)

! import external procedures and variables
!
    use blocks         , only : block_meta
    use blocks         , only : append_metablock
    use blocks         , only : set_last_id, get_last_id
    use blocks         , only : get_mblocks, get_dblocks, get_nleafs
    use coordinates    , only : nn, ng, in, jn, kn, ir, jr, kr
    use coordinates    , only : maxlev, toplev
    use coordinates    , only : xmin, xmax, ymin, ymax, zmin, zmax
    use coordinates    , only : initialize_coordinates, finalize_coordinates
    use error          , only : print_error
    use evolution      , only : step, time, dt, dtn
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5gopen_f, h5gclose_f
    use mpitools       , only : nprocs, nproc
    use random         , only : nseeds, set_seeds

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t) :: gid
    integer        :: ierr, l
    integer        :: lndims, lmaxlev, lmblocks, lnleafs, llast_id
    integer        :: lncells, lnghost, lnproc, lnseeds

! local pointers
!
    type(block_meta), pointer :: pmeta

! allocatable arrays
!
    integer(kind=4), dimension(:), allocatable :: seeds
!
!-------------------------------------------------------------------------------
!
! open the global attributes group
!
    call h5gopen_f(fid, 'attributes', gid, ierr)

! check if the group has been opened successfuly
!
    if (ierr < 0) then
      call print_error("io::read_attributes_h5", "Cannot open the group!")
      return
    end if

! restore integer attributes
!
    call read_attribute(gid, 'ndims'  , lndims  )
    call read_attribute(gid, 'maxlev' , lmaxlev )
    call read_attribute(gid, 'nprocs' , nfiles  )
    call read_attribute(gid, 'nproc'  , lnproc  )
    call read_attribute(gid, 'mblocks', lmblocks)
    call read_attribute(gid, 'nleafs' , lnleafs )
    call read_attribute(gid, 'last_id', llast_id)
    call read_attribute(gid, 'ncells' , lncells )
    call read_attribute(gid, 'nghost' , lnghost )
    call read_attribute(gid, 'nseeds' , lnseeds )
    call read_attribute(gid, 'step'   , step    )
    call read_attribute(gid, 'isnap'  , isnap   )

! restore double precision attributes
!
    call read_attribute(gid, 'xmin', xmin)
    call read_attribute(gid, 'xmax', xmax)
    call read_attribute(gid, 'ymin', ymin)
    call read_attribute(gid, 'ymax', ymax)
    call read_attribute(gid, 'zmin', zmin)
    call read_attribute(gid, 'zmax', zmax)
    call read_attribute(gid, 'time', time)
    call read_attribute(gid, 'dt'  , dt  )
    call read_attribute(gid, 'dtn' , dtn )

! check the number of dimensions
!
    if (lndims /= NDIMS) then
      call print_error("io::read_attributes_h5"                                &
                                 , "The number of dimensions does not match!")
      return
    end if

! check the block dimensions
!
    if (lncells /= nn) then
      call print_error("io::read_attributes_h5"                                &
                                       , "The block dimensions do not match!")
    end if

! check the number of ghost layers
!
    if (lnghost /= ng) then
      call print_error("io::read_attributes_h5"                                &
                               , "The number of ghost layers does not match!")
    end if

! prepare coordinates and rescaling factors if the maximum level has changed
!
    if (lmaxlev > toplev) then

! subtitute the new value of toplev
!
      toplev = lmaxlev

! regenerate coordinates
!
      call finalize_coordinates(ierr)
      call initialize_coordinates(.false., ierr)

! calculate a factor to rescale the block coordinates
!
      dcor = 2**(toplev - maxlev)

    else

! calculate a factor to rescale the block coordinates
!
      ucor = 2**(maxlev - lmaxlev)

    end if

! allocate space for seeds
!
    allocate(seeds(lnseeds))

! store them in the current group
!
    call read_attribute(gid, 'seeds', seeds(:))

! set the seed values
!
    call set_seeds(lnseeds, seeds(:))

! deallocate seed array
!
    deallocate(seeds)

! allocate all metablocks
!
    do l = 1, lmblocks
      call append_metablock(pmeta)
    end do

! check if the number of created metablocks is equal to lbmcloks
!
    if (lmblocks /= get_mblocks()) then
      call print_error("io::read_attributes_h5"                                &
                                     , "Number of metablocks does not match!")
    end if

! allocate an array of pointers with the size llast_id
!
    allocate(block_array(llast_id))

! set the last_id
!
    call set_last_id(llast_id)

! close the group
!
    call h5gclose_f(gid, ierr)

! check if the group has been closed successfuly
!
    if (ierr /= 0) then
      call print_error("io::read_attributes_h5", "Cannot close the group!")
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_attributes_h5
!
!===============================================================================
!
! subroutine WRITE_METABLOCKS_H5:
! ------------------------------
!
!   Subroutine stores all meta blocks with their complete fields in 'metablock'
!   group in a provided file identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_metablocks_h5(fid)

! import procedures and variables from other modules
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : ndims, nchildren, nsides
    use blocks         , only : get_last_id, get_mblocks
    use error          , only : print_error
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5gcreate_f, h5gclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: i, j, k, l, n, p
    integer                        :: iret
    integer(hsize_t), dimension(1) :: am, cm
    integer(hsize_t), dimension(2) :: dm, pm
#if NDIMS == 2
    integer(hsize_t), dimension(4) :: nm
#endif /* NDIMS == 2 */
#if NDIMS == 3
    integer(hsize_t), dimension(5) :: nm
#endif /* NDIMS == 3 */

! local allocatable arrays
!
    integer(kind=4), dimension(:)  , allocatable :: idx
    integer(kind=4), dimension(:)  , allocatable :: par, dat
    integer(kind=4), dimension(:)  , allocatable ::  id, cpu, lev, cfg, ref, lea
    real   (kind=8), dimension(:)  , allocatable :: xmn, xmx, ymn, ymx, zmn, zmx
    integer(kind=4), dimension(:,:), allocatable :: chl, pos, cor
#if NDIMS == 2
    integer(kind=4), dimension(:,:,:,:)  , allocatable :: edges
    integer(kind=4), dimension(:,:,:)    , allocatable :: corners
#endif /* NDIMS == 2 */
#if NDIMS == 3
    integer(kind=4), dimension(:,:,:,:,:), allocatable :: faces
    integer(kind=4), dimension(:,:,:,:,:), allocatable :: edges
    integer(kind=4), dimension(:,:,:,:)  , allocatable :: corners
#endif /* NDIMS == 3 */

! local pointers
!
    type(block_meta), pointer :: pmeta

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_metablocks_h5"
!
!-------------------------------------------------------------------------------
!
! create the group for metadata
!
    call h5gcreate_f(fid, 'metablocks', gid, iret)

! check if the group has been created successfuly
!
    if (iret >= 0) then

! prepate dimensions
!
      am(1) = get_mblocks()
      cm(1) = get_last_id()
      dm(1) = get_mblocks()
      dm(2) = nchildren
      pm(1) = get_mblocks()
      pm(2) = NDIMS
      nm(1) = get_mblocks()
      nm(2) = nsides
      nm(3) = nsides
#if NDIMS == 2
      nm(4) = ndims
#endif /* NDIMS == 2 */
#if NDIMS == 3
      nm(4) = nsides
      nm(5) = ndims
#endif /* NDIMS == 3 */

! only store data from processes that have any meta blocks
!
      if (am(1) > 0) then

! allocate arrays to store meta block fields
!
        allocate(idx(cm(1)))
        allocate(par(am(1)))
        allocate(dat(am(1)))
        allocate(id (am(1)))
        allocate(cpu(am(1)))
        allocate(lev(am(1)))
        allocate(cfg(am(1)))
        allocate(ref(am(1)))
        allocate(lea(am(1)))
        allocate(xmn(am(1)))
        allocate(xmx(am(1)))
        allocate(ymn(am(1)))
        allocate(ymx(am(1)))
        allocate(zmn(am(1)))
        allocate(zmx(am(1)))
        allocate(chl(dm(1),dm(2)))
        allocate(pos(pm(1),pm(2)))
        allocate(cor(pm(1),pm(2)))
#if NDIMS == 2
        allocate(edges  (nm(1),nm(2),nm(3),nm(4)))
        allocate(corners(nm(1),nm(2),nm(3)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        allocate(faces  (nm(1),nm(2),nm(3),nm(4),nm(5)))
        allocate(edges  (nm(1),nm(2),nm(3),nm(4),nm(5)))
        allocate(corners(nm(1),nm(2),nm(3),nm(4)))
#endif /* NDIMS == 3 */

! reset stored arrays
!
        idx(:)           = -1
        par(:)           = -1
        dat(:)           = -1
        lea(:)           = -1
        chl(:,:)         = -1
#if NDIMS == 2
        edges(:,:,:,:)   = -1
        corners(:,:,:)   = -1
#endif /* NDIMS == 2 */
#if NDIMS == 3
        faces(:,:,:,:,:) = -1
        edges(:,:,:,:,:) = -1
        corners(:,:,:,:) = -1
#endif /* NDIMS == 3 */

! reset the block counter
!
        l = 0

! associate pmeta with the first block on the meta block list
!
        pmeta => list_meta

! iterate over all meta blocks and fill in the arrays for storage
!
        do while(associated(pmeta))

! increase the block counter
!
          l = l + 1

! store meta block fields
!
          idx(pmeta%id) = l

          if (associated(pmeta%parent)) par(l) = pmeta%parent%id
          if (associated(pmeta%data)  ) dat(l) = 1

          id (l)   = pmeta%id
          cpu(l)   = pmeta%process
          lev(l)   = pmeta%level
          cfg(l)   = pmeta%conf
          ref(l)   = pmeta%refine
          pos(l,:) = pmeta%pos(:)
          cor(l,:) = pmeta%coords(:)

          if (pmeta%leaf) lea(l) = 1

          xmn(l)   = pmeta%xmin
          xmx(l)   = pmeta%xmax
          ymn(l)   = pmeta%ymin
          ymx(l)   = pmeta%ymax
          zmn(l)   = pmeta%zmin
          zmx(l)   = pmeta%zmax

          do p = 1, nchildren
            if (associated(pmeta%child(p)%ptr)) chl(l,p) = pmeta%child(p)%ptr%id
          end do

! store face, edge and corner neighbor pointers
!
#if NDIMS == 2
          do i = 1, nsides
            do j = 1, nsides
              do n = 1, ndims
                if (associated(pmeta%edges(i,j,n)%ptr))                        &
                                    edges(l,i,j,n) = pmeta%edges(i,j,n)%ptr%id
              end do ! ndims
              if (associated(pmeta%corners(i,j)%ptr))                          &
                                    corners(l,i,j) = pmeta%corners(i,j)%ptr%id
            end do ! i = 1, nsides
          end do ! j = 1, nsides
#endif /* NDIMS == 2 */
#if NDIMS == 3
          do i = 1, nsides
            do j = 1, nsides
              do k = 1, nsides
                do n = 1, ndims
                  if (associated(pmeta%faces(i,j,k,n)%ptr))                    &
                                faces(l,i,j,k,n) = pmeta%faces(i,j,k,n)%ptr%id
                  if (associated(pmeta%edges(i,j,k,n)%ptr))                    &
                                edges(l,i,j,k,n) = pmeta%edges(i,j,k,n)%ptr%id
                end do ! ndims
                if (associated(pmeta%corners(i,j,k)%ptr))                      &
                                corners(l,i,j,k) = pmeta%corners(i,j,k)%ptr%id
              end do ! i = 1, nsides
            end do ! j = 1, nsides
          end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pmeta with the next block on the list
!
          pmeta => pmeta%next

        end do ! over all meta blocks

! store meta block data in the HDF5 file
!
        call write_array(gid, 'indices', cm(1)  , idx)
        call write_array(gid, 'parent' , am(1)  , par)
        call write_array(gid, 'data'   , am(1)  , dat)
        call write_array(gid, 'id'     , am(1)  , id )
        call write_array(gid, 'cpu'    , am(1)  , cpu)
        call write_array(gid, 'level'  , am(1)  , lev)
        call write_array(gid, 'config' , am(1)  , cfg)
        call write_array(gid, 'refine' , am(1)  , ref)
        call write_array(gid, 'leaf'   , am(1)  , lea)
        call write_array(gid, 'xmin'   , am(1)  , xmn)
        call write_array(gid, 'xmax'   , am(1)  , xmx)
        call write_array(gid, 'ymin'   , am(1)  , ymn)
        call write_array(gid, 'ymax'   , am(1)  , ymx)
        call write_array(gid, 'zmin'   , am(1)  , zmn)
        call write_array(gid, 'zmax'   , am(1)  , zmx)
        call write_array(gid, 'child'  , dm(:)  , chl(:,:))
        call write_array(gid, 'pos'    , pm(:)  , pos(:,:))
        call write_array(gid, 'coord'  , pm(:)  , cor(:,:))
#if NDIMS == 2
        call write_array(gid, 'edges'  , nm(1:4), edges(:,:,:,:))
        call write_array(gid, 'corners', nm(1:3), corners(:,:,:))
#endif /* NDIMS == 2 */
#if NDIMS == 3
        call write_array(gid, 'faces'  , nm(1:5), faces(:,:,:,:,:))
        call write_array(gid, 'edges'  , nm(1:5), edges(:,:,:,:,:))
        call write_array(gid, 'corners', nm(1:4), corners(:,:,:,:))
#endif /* NDIMS == 3 */

! deallocate allocated arrays
!
        if (allocated(idx))     deallocate(idx)
        if (allocated(par))     deallocate(par)
        if (allocated(dat))     deallocate(dat)
        if (allocated(id) )     deallocate(id)
        if (allocated(cpu))     deallocate(cpu)
        if (allocated(lev))     deallocate(lev)
        if (allocated(cfg))     deallocate(cfg)
        if (allocated(ref))     deallocate(ref)
        if (allocated(lea))     deallocate(lea)
        if (allocated(xmn))     deallocate(xmn)
        if (allocated(xmx))     deallocate(xmx)
        if (allocated(ymn))     deallocate(ymn)
        if (allocated(ymx))     deallocate(ymx)
        if (allocated(zmn))     deallocate(zmn)
        if (allocated(zmx))     deallocate(zmx)
        if (allocated(chl))     deallocate(chl)
        if (allocated(cor))     deallocate(cor)
#if NDIMS == 3
        if (allocated(faces))   deallocate(faces)
#endif /* NDIMS == 3 */
        if (allocated(edges))   deallocate(edges)
        if (allocated(corners)) deallocate(corners)

      end if ! meta blocks > 0

! close the group
!
      call h5gclose_f(gid, iret)

! check if the group has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the group
!
        call print_error(fname, "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error(fname, "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_metablocks_h5
!
!===============================================================================
!
! subroutine READ_METABLOCKS_H5:
! -----------------------------
!
!   Subroutine restores all meta blocks with their complete fields from
!   'metablock' group in a provided restart file identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine read_metablocks_h5(fid)

! import procedures and variables from other modules
!
    use blocks         , only : block_meta, list_meta
    use blocks         , only : ndims, nchildren, nsides
    use blocks         , only : get_mblocks
    use blocks         , only : metablock_set_id, metablock_set_process
    use blocks         , only : metablock_set_refinement
    use blocks         , only : metablock_set_configuration
    use blocks         , only : metablock_set_level, metablock_set_position
    use blocks         , only : metablock_set_coordinates, metablock_set_bounds
    use blocks         , only : metablock_set_leaf
    use error          , only : print_error
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5gopen_f, h5gclose_f
    use mpitools       , only : nprocs

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: i, j, k, l, p, n, ip, lcpu
    integer                        :: err
    integer(hsize_t), dimension(1) :: am
    integer(hsize_t), dimension(2) :: dm, pm
#if NDIMS == 2
    integer(hsize_t), dimension(4) :: nm
#endif /* NDIMS == 2 */
#if NDIMS == 3
    integer(hsize_t), dimension(5) :: nm
#endif /* NDIMS == 3 */

! local allocatable arrays
!
    integer(kind=4), dimension(:)  , allocatable :: idx
    integer(kind=4), dimension(:)  , allocatable :: par, dat
    integer(kind=4), dimension(:)  , allocatable ::  id, cpu, lev, cfg, ref, lea
    real   (kind=8), dimension(:)  , allocatable :: xmn, xmx, ymn, ymx, zmn, zmx
    integer(kind=4), dimension(:,:), allocatable :: chl, pos, cor
#if NDIMS == 2
    integer(kind=4), dimension(:,:,:,:)  , allocatable :: edges
    integer(kind=4), dimension(:,:,:)    , allocatable :: corners
#endif /* NDIMS == 2 */
#if NDIMS == 3
    integer(kind=4), dimension(:,:,:,:,:), allocatable :: faces
    integer(kind=4), dimension(:,:,:,:,:), allocatable :: edges
    integer(kind=4), dimension(:,:,:,:)  , allocatable :: corners
#endif /* NDIMS == 3 */

! local pointers
!
    type(block_meta), pointer :: pmeta

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_metablocks_h5"
!
!-------------------------------------------------------------------------------
!
! prepare last cpu index
!
    lcpu = nprocs - 1

! open metablock group
!
    call h5gopen_f(fid, 'metablocks', gid, err)

! check if the group has been opened successfuly
!
    if (err >= 0) then

! prepate dimensions
!
      am(1) = get_mblocks()
      dm(1) = get_mblocks()
      dm(2) = nchildren
      pm(1) = get_mblocks()
      pm(2) = NDIMS
      nm(1) = get_mblocks()
      nm(2) = nsides
      nm(3) = nsides
#if NDIMS == 2
      nm(4) = ndims
#endif /* NDIMS == 2 */
#if NDIMS == 3
      nm(4) = nsides
      nm(5) = ndims
#endif /* NDIMS == 3 */

! allocate arrays to restore metablocks data
!
      allocate(id (am(1)))
      allocate(cpu(am(1)))
      allocate(lev(am(1)))
      allocate(par(am(1)))
      allocate(dat(am(1)))
      allocate(cfg(am(1)))
      allocate(ref(am(1)))
      allocate(lea(am(1)))
      allocate(xmn(am(1)))
      allocate(xmx(am(1)))
      allocate(ymn(am(1)))
      allocate(ymx(am(1)))
      allocate(zmn(am(1)))
      allocate(zmx(am(1)))
      allocate(chl(dm(1),dm(2)))
      allocate(pos(pm(1),pm(2)))
      allocate(cor(pm(1),pm(2)))
#if NDIMS == 2
      allocate(edges  (nm(1),nm(2),nm(3),nm(4)))
      allocate(corners(nm(1),nm(2),nm(3)))
#endif /* NDIMS == 2 */
#if NDIMS == 3
      allocate(faces  (nm(1),nm(2),nm(3),nm(4),nm(5)))
      allocate(edges  (nm(1),nm(2),nm(3),nm(4),nm(5)))
      allocate(corners(nm(1),nm(2),nm(3),nm(4)))
#endif /* NDIMS == 3 */

! reset vectors
!
      par(:)           = -1
      dat(:)           = -1
      lea(:)           = -1
      chl(:,:)         = -1
#if NDIMS == 2
      edges(:,:,:,:)   = -1
      corners(:,:,:)   = -1
#endif /* NDIMS == 2 */
#if NDIMS == 3
      faces(:,:,:,:,:) = -1
      edges(:,:,:,:,:) = -1
      corners(:,:,:,:) = -1
#endif /* NDIMS == 3 */

! read metablock fields from the HDF5 file
!
      call read_array(gid, 'id'     , am(:), id (:))
      call read_array(gid, 'cpu'    , am(:), cpu(:))
      call read_array(gid, 'level'  , am(:), lev(:))
      call read_array(gid, 'config' , am(:), cfg(:))
      call read_array(gid, 'refine' , am(:), ref(:))
      call read_array(gid, 'leaf'   , am(:), lea(:))
      call read_array(gid, 'parent' , am(:), par(:))
      call read_array(gid, 'xmin'   , am(:), xmn(:))
      call read_array(gid, 'xmax'   , am(:), xmx(:))
      call read_array(gid, 'ymin'   , am(:), ymn(:))
      call read_array(gid, 'ymax'   , am(:), ymx(:))
      call read_array(gid, 'zmin'   , am(:), zmn(:))
      call read_array(gid, 'zmax'   , am(:), zmx(:))
      call read_array(gid, 'pos'    , pm(:), pos(:,:))
      call read_array(gid, 'coord'  , pm(:), cor(:,:))
      call read_array(gid, 'child'  , dm(:), chl(:,:))
#if NDIMS == 2
      call read_array(gid, 'edges'  , nm(1:4), edges(:,:,:,:))
      call read_array(gid, 'corners', nm(1:3), corners(:,:,:))
#endif /* NDIMS == 2 */
#if NDIMS == 3
      call read_array(gid, 'faces'  , nm(1:5), faces(:,:,:,:,:))
      call read_array(gid, 'edges'  , nm(1:5), edges(:,:,:,:,:))
      call read_array(gid, 'corners', nm(1:4), corners(:,:,:,:))
#endif /* NDIMS == 3 */

! check if the maximum level has been changed, is so, rescale block coordinates
!
      if (dcor .gt. 1) then
        cor(:,:) = cor(:,:) / dcor
      end if
      if (ucor .gt. 1) then
        cor(:,:) = cor(:,:) * ucor
      end if

! reset the block counter
!
      l = 0

! associate pmeta with the first block on the meta block list
!
      pmeta => list_meta

! iterate over all meta blocks and restore their fields
!
      do while(associated(pmeta))

! increase the block counter
!
        l = l + 1

! restore meta block fields
!
        block_array(id(l))%ptr => pmeta

        call metablock_set_id           (pmeta, id (l))
        call metablock_set_process      (pmeta, cpu(l))
        call metablock_set_refinement   (pmeta, ref(l))
        call metablock_set_configuration(pmeta, cfg(l))
        call metablock_set_level        (pmeta, lev(l))
        call metablock_set_position     (pmeta, pos(l,1), pos(l,2), pos(l,3))
        call metablock_set_coordinates  (pmeta, cor(l,1), cor(l,2), cor(l,3))
        call metablock_set_bounds       (pmeta, xmn(l), xmx(l), ymn(l), ymx(l) &
                                                              , zmn(l), zmx(l))

        if (lea(l) == 1) call metablock_set_leaf(pmeta)

! associate pmeta with the next block on the list
!
        pmeta => pmeta%next

      end do ! over all meta blocks

! reset the block counter
!
      l = 0

! associate pmeta with the first block on the meta block list
!
      pmeta => list_meta

! iterate over all meta blocks and restore their pointers
!
      do while(associated(pmeta))

! increase the block counter
!
        l = l + 1

! restore %parent pointer
!
        if (par(l) > 0) pmeta%parent => block_array(par(l))%ptr

! restore %child pointers
!
        do p = 1, nchildren
          if (chl(l,p) > 0) then
            pmeta%child(p)%ptr => block_array(chl(l,p))%ptr
          end if
        end do ! p = 1, nchildren

! restore %faces, %edges and %corners neighbor pointers
!
#if NDIMS == 2
        do i = 1, nsides
          do j = 1, nsides
            do n = 1, ndims
              ip = edges(l,i,j,n)
              if (ip > 0) pmeta%edges(i,j,n)%ptr => block_array(ip)%ptr
            end do ! n = 1, ndims
            ip = corners(l,i,j)
            if (ip > 0) pmeta%corners(i,j)%ptr => block_array(ip)%ptr
          end do ! i = 1, nsides
        end do ! j = 1, nsides
#endif /* NDIMS == 2 */
#if NDIMS == 3
        do i = 1, nsides
          do j = 1, nsides
            do k = 1, nsides
              do n = 1, ndims
                ip = faces(l,i,j,k,n)
                if (ip > 0) pmeta%faces(i,j,k,n)%ptr => block_array(ip)%ptr
                ip = edges(l,i,j,k,n)
                if (ip > 0) pmeta%edges(i,j,k,n)%ptr => block_array(ip)%ptr
              end do ! n = 1, ndims
              ip = corners(l,i,j,k)
              if (ip > 0) pmeta%corners(i,j,k)%ptr => block_array(ip)%ptr
            end do ! i = 1, nsides
          end do ! j = 1, nsides
        end do ! k = 1, nsides
#endif /* NDIMS == 3 */

! associate pmeta with the next block on the list
!
        pmeta => pmeta%next

      end do ! over all meta blocks

! deallocate allocatable arrays
!
      if (allocated(id) )     deallocate(id )
      if (allocated(par))     deallocate(par)
      if (allocated(dat))     deallocate(dat)
      if (allocated(cpu))     deallocate(cpu)
      if (allocated(lev))     deallocate(lev)
      if (allocated(cfg))     deallocate(cfg)
      if (allocated(ref))     deallocate(ref)
      if (allocated(lea))     deallocate(lea)
      if (allocated(xmn))     deallocate(xmn)
      if (allocated(xmx))     deallocate(xmx)
      if (allocated(ymn))     deallocate(ymn)
      if (allocated(ymx))     deallocate(ymx)
      if (allocated(zmn))     deallocate(zmn)
      if (allocated(zmx))     deallocate(zmx)
      if (allocated(chl))     deallocate(chl)
      if (allocated(cor))     deallocate(cor)
#if NDIMS == 3
      if (allocated(faces))   deallocate(faces)
#endif /* NDIMS == 3 */
      if (allocated(edges))   deallocate(edges)
      if (allocated(corners)) deallocate(corners)

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err > 0) then

! print error about the problem with closing the group
!
        call print_error(fname, "Cannot close metablock group!")

      end if

    else

! print error about the problem with opening the group
!
      call print_error(fname, "Cannot open metablock group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_metablocks_h5
!
!===============================================================================
!
! subroutine WRITE_DATABLOCKS_H5:
! ------------------------------
!
!   Subroutine writes all data block fields in the new group 'datablocks'
!   in the provided handler to the HDF5 file.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_datablocks_h5(fid)

! import external procedures and variables
!
    use blocks         , only : ndims
    use blocks         , only : block_meta, block_data, list_data
    use blocks         , only : get_dblocks
    use coordinates    , only : im, jm, km
    use equations      , only : nv
    use error          , only : print_error
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5gcreate_f, h5gclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine variables
!
    integer(hid_t), intent(in) :: fid

! local pointers
!
    type(block_meta), pointer  :: pmeta
    type(block_data), pointer  :: pdata

! local variables
!
    integer(hid_t)             :: gid
    integer(kind=4)            :: l
    integer                    :: err

! local arrays
!
    integer(hsize_t), dimension(1) :: am
    integer(hsize_t), dimension(5) :: dm

! local allocatable arrays
!
    integer(kind=4), dimension(:)          , allocatable :: id
    real(kind=8)   , dimension(:,:,:,:,:)  , allocatable :: uv, qv
!
!-------------------------------------------------------------------------------
!
! create a new group for storing data blocks
!
    call h5gcreate_f(fid, 'datablocks', gid, err)

! check if the group has been created successfuly
!
    if (err >= 0) then

! store data blocks only if there is at least one belonging to
! the current process
!
      if (get_dblocks() > 0) then

! prepate dimensions
!
        am(1) = get_dblocks()
        dm(1) = get_dblocks()
        dm(2) = nv
        dm(3) = im
        dm(4) = jm
        dm(5) = km

! allocate arrays to store associated meta block identifiers, conserved and
! primitive variables
!
        allocate(id(am(1)))
        allocate(uv(dm(1),dm(2),dm(3),dm(4),dm(5)))
        allocate(qv(dm(1),dm(2),dm(3),dm(4),dm(5)))

! reset the block counter
!
        l = 0

! associate the pointer with the first block in the data block list
!
        pdata => list_data

! iterate over all data blocks and fill in the arrays id, u, and q
!
        do while(associated(pdata))

#ifdef DEBUG
! store only data from data blocks associated with meta blocks
!
          if (associated(pdata%meta)) then

! increase the block counter
!
            l             = l + 1

! fill in the meta block ID array
!
            id(l)         = pdata%meta%id

! fill in the conservative and primitive variable arrays
!
            uv(l,:,:,:,:) = pdata%u(:,:,:,:)
            qv(l,:,:,:,:) = pdata%q(:,:,:,:)

          else ! meta block not associated

! print error about the lack of associated meta block
!
            call print_error("io::write_datablocks_h5"                         &
                                                , "Meta block no associated!")

          end if ! meta block not associated
#else /* DEBUG */
! increase the block counter
!
          l             = l + 1

! fill in the meta block ID array
!
          id(l)         = pdata%meta%id

! fill in the conservative and primitive variable arrays
!
          uv(l,:,:,:,:) = pdata%u(:,:,:,:)
          qv(l,:,:,:,:) = pdata%q(:,:,:,:)
#endif /* DEBUG */

! associate the pointer with the next data block on the list
!
          pdata => pdata%next

        end do ! data blocks

! store data arrays in the current group
!
        call write_array(gid, 'meta', am(1), id)
        call write_array(gid, 'uvar', dm(:), uv)
        call write_array(gid, 'qvar', dm(:), qv)

! deallocate allocatable arrays
!
        if (allocated(id)) deallocate(id)
        if (allocated(uv)) deallocate(uv)
        if (allocated(qv)) deallocate(qv)

      end if ! dblocks > 0

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err > 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_datablocks_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_datablocks_h5", "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_datablocks_h5
!
!===============================================================================
!
! subroutine READ_DATABLOCKS_H5:
! -----------------------------
!
!   Subroutine reads all data blocks stored in the group 'datablocks' of
!   the provided handler to the HDF5 restart file.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine read_datablocks_h5(fid)

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_data
    use blocks         , only : append_datablock, link_blocks
    use coordinates    , only : im, jm, km
    use equations      , only : nv
    use error          , only : print_error
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5gopen_f, h5gclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine variables
!
    integer(hid_t), intent(in)     :: fid

! local pointers
!
    type(block_data), pointer      :: pdata

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: l
    integer                        :: dblocks, ierr

! local arrays
!
    integer(hsize_t), dimension(5) :: dm

! local allocatable arrays
!
    integer(kind=4), dimension(:)        , allocatable :: id
    real(kind=8)   , dimension(:,:,:,:,:), allocatable :: uv, qv
!
!-------------------------------------------------------------------------------
!
! read the number of data blocks
!
    call h5gopen_f(fid, 'attributes', gid, ierr)
    if (ierr /= 0) then
      call print_error("io::read_datablocks_h5"                                &
                                         , "Cannot open the attribute group!")
      return
    end if
    call read_attribute(gid, 'dblocks', dblocks)
    call h5gclose_f(gid, ierr)
    if (ierr /= 0) then
      call print_error("io::read_datablocks_h5"                                &
                                        , "Cannot close the attribute group!")
      return
    end if

! restore data blocks only if there are any
!
    if (dblocks > 0) then

! fill out dimensions dm(:)
!
      dm(1) = dblocks
      dm(2) = nv
      dm(3) = im
      dm(4) = jm
      dm(5) = km

! allocate arrays to read data
!
      allocate(id(dm(1)))
      allocate(uv(dm(1),dm(2),dm(3),dm(4),dm(5)))
      allocate(qv(dm(1),dm(2),dm(3),dm(4),dm(5)))

! open the group 'datablocks'
!
      call h5gopen_f(fid, 'datablocks', gid, ierr)
      if (ierr /= 0) then
        call print_error("io::read_datablocks_h5"                              &
                                        , "Cannot open the data block group!")
        return
      end if

! read array data from the HDF5 file
!
      call read_array(gid, 'meta', dm(1:1), id(:)        )
      call read_array(gid, 'uvar', dm(1:5), uv(:,:,:,:,:))
      call read_array(gid, 'qvar', dm(1:5), qv(:,:,:,:,:))

! close the data block group
!
      call h5gclose_f(gid, ierr)
      if (ierr /= 0) then
        call print_error("io::read_datablocks_h5"                              &
                                       , "Cannot close the data block group!")
        return
      end if

! iterate over data blocks, allocate them and fill out their fields
!
      do l = 1, dm(1)

! allocate and append to the end of the list a new datablock
!
        call append_datablock(pdata)

! associate the corresponding meta block with the current data block
!
        call link_blocks(block_array(id(l))%ptr, pdata)

! fill out the array of conservative and primitive variables
!
        pdata%u(:,:,:,:) = uv(l,:,:,:,:)
        pdata%q(:,:,:,:) = qv(l,:,:,:,:)

      end do ! l = 1, dm(1)

! deallocate allocatable arrays
!
      if (allocated(id)) deallocate(id)
      if (allocated(uv)) deallocate(uv)
      if (allocated(qv)) deallocate(qv)

    end if ! dblocks > 0

!-------------------------------------------------------------------------------
!
  end subroutine read_datablocks_h5
!
!===============================================================================
!
! write_coordinates_h5: subroutine writes data block coordinates and other
!                       variables which determine geometrical position of
!                       the blocks
!
! info: this subroutine stores coordinates
!
! arguments:
!   fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_coordinates_h5(fid)

! references to other modules
!
    use blocks, only : block_meta, block_data, list_data
    use blocks, only : nsides
    use blocks, only : get_dblocks
    use coordinates, only : maxlev, toplev
    use error , only : print_error
    use hdf5  , only : hid_t, hsize_t
    use hdf5  , only : h5gcreate_f, h5gclose_f
    use coordinates, only : adx, ady, adz

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! HDF5 variables
!
    integer(hid_t)    :: gid

! local variables
!
    integer           :: err
    integer(kind=4)   :: l
    integer(hsize_t)  :: am(1), cm(2), rm(2), dm(3)

! local allocatable arrays
!
    integer(kind=4), dimension(:)    , allocatable :: lev, ref
    integer(kind=4), dimension(:,:)  , allocatable :: cor
    real   (kind=8), dimension(:,:,:), allocatable :: bnd

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create a group to store global attributes
!
    call h5gcreate_f(fid, 'coordinates', gid, err)

! check if the group has been created successfuly
!
    if (err .ge. 0) then

! store coordinates only if there are some data blocks on the current processor
!
      if (get_dblocks() .gt. 0) then

! prepare dimensions
!
        am(1) = maxlev
        cm(1) = get_dblocks()
        cm(2) = NDIMS
        rm(1) = maxlev
        rm(2) = NDIMS
        dm(1) = get_dblocks()
        dm(2) = NDIMS
        dm(3) = nsides

! allocate arrays to store coordinates
!
        allocate(lev(cm(1)))
        allocate(ref(cm(1)))
        allocate(cor(cm(1),cm(2)))
        allocate(bnd(dm(1),dm(2),dm(3)))

! iterate over all data blocks and fill in the arrays
!
        l = 1
        pdata => list_data
        do while(associated(pdata))

! fill in the level array
!
          lev(l)     = pdata%meta%level

! fill in the refinement flag
!
          ref(l)     = pdata%meta%refine

! fill in the coordinate array
!
          cor(l,:)   = pdata%meta%coords(:)

! fill in the bounds array
!
          bnd(l,1,1) = pdata%meta%xmin
          bnd(l,1,2) = pdata%meta%xmax
          bnd(l,2,1) = pdata%meta%ymin
          bnd(l,2,2) = pdata%meta%ymax
#if NDIMS == 3
          bnd(l,3,1) = pdata%meta%zmin
          bnd(l,3,2) = pdata%meta%zmax
#endif /* NDIMS == 3 */

          l = l + 1
          pdata => pdata%next
        end do

! write the arrays to the HDF5 file
!
        call write_array(gid, 'levels', cm(1), lev)
        call write_array(gid, 'refine', cm(1), ref)
        call write_array(gid, 'coords', cm(:), cor)
        call write_array(gid, 'bounds', dm(:), bnd)
        call write_array(gid, 'dx'    , am(1), adx(1:maxlev))
        call write_array(gid, 'dy'    , am(1), ady(1:maxlev))
        call write_array(gid, 'dz'    , am(1), adz(1:maxlev))

! deallocate temporary arrays
!
        if (allocated(lev)) deallocate(lev)
        if (allocated(ref)) deallocate(ref)
        if (allocated(cor)) deallocate(cor)
        if (allocated(bnd)) deallocate(bnd)

      end if ! dblocks > 0

! close the attribute group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_coordinates_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_coordinates_h5", "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_coordinates_h5
!
!===============================================================================
!
! subroutine WRITE_PRIMITIVE_VARIABLES_H5:
! ---------------------------------------
!
!   Subroutine groups each primitive variable from all data blocks and writes
!   it as an array in the HDF5 dataset connected to the input HDF file
!   identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_primitive_variables_h5(fid)

! references to other modules
!
    use blocks       , only : block_data, list_data
    use blocks       , only : get_dblocks
    use coordinates  , only : im, jm, km, in, jn, kn
    use coordinates  , only : ib, jb, kb, ie, je, ke
    use equations    , only : nv, pvars
    use error        , only : print_error
    use hdf5         , only : hid_t, hsize_t
    use hdf5         , only : h5gcreate_f, h5gclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! HDF5 variables
!
    integer(hid_t)             :: gid
    integer(hsize_t)           :: dm(4)

! local variables
!
    integer         :: err
    integer(kind=4) :: i, j, k, l, n
    integer(kind=4) :: il, jl, kl, iu, ju, ku

! local allocatable arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable :: qarr

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create a group to store variables
!
    call h5gcreate_f(fid, 'variables', gid, err)

! check if the group was created successfuly
!
    if (err >= 0) then

! store variables only if there is at least one data block associated with
! the current process
!
      if (get_dblocks() > 0) then

! prepare dimensions and index limits
!
        dm(1) = get_dblocks()
        if (with_ghosts) then
          dm(2) = im
          dm(3) = jm
          dm(4) = km

          il    =  1
          jl    =  1
          kl    =  1
          iu    = im
          ju    = jm
          ku    = km
        else
          dm(2) = in
          dm(3) = jn
          dm(4) = kn

          il    = ib
          jl    = jb
          kl    = kb
          iu    = ie
          ju    = je
          ku    = ke
        end if

! allocate array to group a variable from all data blocks
!
        allocate(qarr(dm(1),dm(2),dm(3),dm(4)))

! iterate over all variables
!
        do n = 1, nv

! reset the block counter
!
          l = 0

! assosiate the block pointer with the first data block on the list
!
          pdata => list_data

! iterate over all data blocks and copy the variable from each of them to
! the allocate array
!
          do while(associated(pdata))

! increase the data block counter
!
            l = l + 1

! copy the variable from the current data block
!
            qarr(l,1:dm(2),1:dm(3),1:dm(4)) = pdata%q(n,il:iu,jl:ju,kl:ku)

! assign the pointer with the next data block on the list
!
            pdata => pdata%next

          end do ! pdata=>list_data

! write the variable array to the HDF5 file
!
          call write_array(gid, trim(pvars(n)), dm, qarr)

        end do ! n = 1, nv

! deallocate allocatable array
!
        if (allocated(qarr)) deallocate(qarr)

      end if ! dblocks > 0

! close the variable group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err > 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_primitive_variables_h5"                    &
                                                  , "Cannot close the group!")

      end if

    else ! error with creating a group

! print error about the problem with creating the group
!
      call print_error("io::write_primitive_variables_h5"                      &
                                                 , "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_primitive_variables_h5
!
!===============================================================================
!
! subroutine WRITE_CONSERVATIVE_VARIABLES_H5:
! ------------------------------------------
!
!   Subroutine groups each conservative variable from all data blocks and writes
!   it as an array in the HDF5 dataset connected to the input HDF file
!   identifier.
!
!   Arguments:
!
!     fid - the HDF5 file identifier;
!
!===============================================================================
!
  subroutine write_conservative_variables_h5(fid)

! references to other modules
!
    use blocks       , only : block_data, list_data
    use blocks       , only : get_dblocks
    use coordinates  , only : im, jm, km, in, jn, kn
    use coordinates  , only : ib, jb, kb, ie, je, ke
    use equations    , only : nv, cvars
    use error        , only : print_error
    use hdf5         , only : hid_t, hsize_t
    use hdf5         , only : h5gcreate_f, h5gclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t), intent(in) :: fid

! HDF5 variables
!
    integer(hid_t)             :: gid
    integer(hsize_t)           :: dm(4)

! local variables
!
    integer         :: err
    integer(kind=4) :: i, j, k, l, n
    integer(kind=4) :: il, jl, kl, iu, ju, ku

! local allocatable arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable :: qarr

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create a group to store variables
!
    call h5gcreate_f(fid, 'variables', gid, err)

! check if the group was created successfuly
!
    if (err >= 0) then

! store variables only if there is at least one data block associated with
! the current process
!
      if (get_dblocks() > 0) then

! prepare dimensions and index limits
!
        dm(1) = get_dblocks()
        if (with_ghosts) then
          dm(2) = im
          dm(3) = jm
          dm(4) = km

          il    =  1
          jl    =  1
          kl    =  1
          iu    = im
          ju    = jm
          ku    = km
        else
          dm(2) = in
          dm(3) = jn
          dm(4) = kn

          il    = ib
          jl    = jb
          kl    = kb
          iu    = ie
          ju    = je
          ku    = ke
        end if

! allocate array to group a variable from all data blocks
!
        allocate(qarr(dm(1),dm(2),dm(3),dm(4)))

! iterate over all variables
!
        do n = 1, nv

! reset the block counter
!
          l = 0

! assosiate the block pointer with the first data block on the list
!
          pdata => list_data

! iterate over all data blocks and copy the variable from each of them to
! the allocate array
!
          do while(associated(pdata))

! increase the data block counter
!
            l = l + 1

! copy the variable from the current data block
!
            qarr(l,1:dm(2),1:dm(3),1:dm(4)) = pdata%u(n,il:iu,jl:ju,kl:ku)

! assign the pointer with the next data block on the list
!
            pdata => pdata%next

          end do ! pdata=>list_data

! write the variable array to the HDF5 file
!
          call write_array(gid, trim(cvars(n)), dm, qarr)

        end do ! n = 1, nv

! deallocate allocatable array
!
        if (allocated(qarr)) deallocate(qarr)

      end if ! dblocks > 0

! close the variable group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err > 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_conservative_variables_h5"                 &
                                                  , "Cannot close the group!")

      end if

    else ! error with creating a group

! print error about the problem with creating the group
!
      call print_error("io::write_conservative_variables_h5"                   &
                                                 , "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_conservative_variables_h5
!
!===============================================================================
!
! WRITE_ATTRIBUTE SUBROUTINES
!
!===============================================================================
!
! subroutine WRITE_SCALAR_ATTRIBUTE_INTEGER_H5:
! --------------------------------------------
!
!   Subroutine stores a value of the integer attribute in the group provided
!   by an identifier and the attribute name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute should be linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine write_scalar_attribute_integer_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5acreate_f, h5awrite_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)  , intent(in) :: gid
    character(len=*), intent(in) :: aname
    integer(kind=4) , intent(in) :: avalue

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                       "io::write_scalar_attribute_integer_h5"
!
!-------------------------------------------------------------------------------
!
! create space for the attribute value
!
    call h5screate_simple_f(1, am, sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                       , "Cannot create space for attribute :" // trim(aname))
      return
    end if

! create the attribute in the given group
!
    call h5acreate_f(gid, aname, H5T_NATIVE_INTEGER, sid, aid, ierr)
    if (ierr == 0) then

! write the attribute data
!
      call h5awrite_f(aid, H5T_NATIVE_INTEGER, avalue, am, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                      , "Cannot write the attribute data in :" // trim(aname))
      end if

! close the attribute
!
      call h5aclose_f(aid, ierr)
      if (ierr /= 0) then
        call print_error(fname, "Cannot close attribute :" // trim(aname))
      end if

    else
      call print_error(fname, "Cannot create attribute :" // trim(aname))
    end if

! release the space
!
    call h5sclose_f(sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot close space for attribute :" // trim(aname))
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_scalar_attribute_integer_h5
!
!===============================================================================
!
! subroutine WRITE_SCALAR_ATTRIBUTE_DOUBLE_H5:
! -------------------------------------------
!
!   Subroutine stores a value of the double precision attribute in the group
!   provided by an identifier and the attribute name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute should be linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine write_scalar_attribute_double_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5acreate_f, h5awrite_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)  , intent(in) :: gid
    character(len=*), intent(in) :: aname
    real(kind=8)    , intent(in) :: avalue

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                        "io::write_scalar_attribute_double_h5"
!
!-------------------------------------------------------------------------------
!
! create space for the attribute value
!
    call h5screate_simple_f(1, am, sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                       , "Cannot create space for attribute :" // trim(aname))
      return
    end if

! create the attribute in the given group
!
    call h5acreate_f(gid, aname, H5T_NATIVE_DOUBLE, sid, aid, ierr)
    if (ierr == 0) then

! write the attribute data
!
      call h5awrite_f(aid, H5T_NATIVE_DOUBLE, avalue, am, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                      , "Cannot write the attribute data in :" // trim(aname))
      end if

! close the attribute
!
      call h5aclose_f(aid, ierr)
      if (ierr /= 0) then
        call print_error(fname, "Cannot close attribute :" // trim(aname))
      end if

    else
      call print_error(fname, "Cannot create attribute :" // trim(aname))
    end if

! release the space
!
    call h5sclose_f(sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot close space for attribute :" // trim(aname))
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_scalar_attribute_double_h5
!
!===============================================================================
!
! subroutine WRITE_VECTOR_ATTRIBUTE_INTEGER_H5:
! --------------------------------------------
!
!   Subroutine stores a vector of the integer attribute in the group provided
!   by an identifier and the attribute name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute should be linked;
!     aname  - the attribute name;
!     avalue - the attribute values;
!
!===============================================================================
!
  subroutine write_vector_attribute_integer_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5acreate_f, h5awrite_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: aname
    integer(kind=4) , dimension(:), intent(in) :: avalue

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                       "io::write_vector_attribute_integer_h5"
!
!-------------------------------------------------------------------------------
!
! set the proper attribute length
!
    am(1) = size(avalue)

! create space for the attribute value
!
    call h5screate_simple_f(1, am, sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                       , "Cannot create space for attribute :" // trim(aname))
      return
    end if

! create the attribute in the given group
!
    call h5acreate_f(gid, aname, H5T_NATIVE_INTEGER, sid, aid, ierr)
    if (ierr == 0) then

! write the attribute data
!
      call h5awrite_f(aid, H5T_NATIVE_INTEGER, avalue, am, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                      , "Cannot write the attribute data in :" // trim(aname))
      end if

! close the attribute
!
      call h5aclose_f(aid, ierr)
      if (ierr /= 0) then
        call print_error(fname, "Cannot close attribute :" // trim(aname))
      end if

    else
      call print_error(fname, "Cannot create attribute :" // trim(aname))
    end if

! release the space
!
    call h5sclose_f(sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot close space for attribute :" // trim(aname))
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_vector_attribute_integer_h5
!
!===============================================================================
!
! subroutine WRITE_VECTOR_ATTRIBUTE_DOUBLE_H5:
! -------------------------------------------
!
!   Subroutine stores a vector of the double precision attribute in the group
!   provided by an identifier and the attribute name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute should be linked;
!     aname  - the attribute name;
!     avalue - the attribute values;
!
!===============================================================================
!
  subroutine write_vector_attribute_double_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5acreate_f, h5awrite_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: aname
    real(kind=8)    , dimension(:), intent(in) :: avalue

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                        "io::write_vector_attribute_double_h5"
!
!-------------------------------------------------------------------------------
!
! set the proper attribute length
!
    am(1) = size(avalue)

! create space for the attribute value
!
    call h5screate_simple_f(1, am, sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                       , "Cannot create space for attribute :" // trim(aname))
      return
    end if

! create the attribute in the given group
!
    call h5acreate_f(gid, aname, H5T_NATIVE_DOUBLE, sid, aid, ierr)
    if (ierr == 0) then

! write the attribute data
!
      call h5awrite_f(aid, H5T_NATIVE_DOUBLE, avalue, am, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                      , "Cannot write the attribute data in :" // trim(aname))
      end if

! close the attribute
!
      call h5aclose_f(aid, ierr)
      if (ierr /= 0) then
        call print_error(fname, "Cannot close attribute :" // trim(aname))
      end if

    else
      call print_error(fname, "Cannot create attribute :" // trim(aname))
    end if

! release the space
!
    call h5sclose_f(sid, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot close space for attribute :" // trim(aname))
    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_vector_attribute_double_h5

!===============================================================================
!
! READ_ATTRIBUTE SUBROUTINES
!
!===============================================================================
!
! subroutine READ_SCALAR_ATTRIBUTE_INTEGER_H5:
! -------------------------------------------
!
!   Subroutine reads a value of the integer attribute provided by the group
!   identifier to which it is linked and its name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute is linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine read_scalar_attribute_integer_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5aexists_by_name_f
    use hdf5           , only : h5aopen_by_name_f, h5aread_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)  , intent(in)    :: gid
    character(len=*), intent(in)    :: aname
    integer(kind=4) , intent(inout) :: avalue

! local variables
!
    logical                         :: exists = .false.
    integer(hid_t)                  :: aid
    integer(hsize_t), dimension(1)  :: am     = (/ 1 /)
    integer                         :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                        "io::read_scalar_attribute_integer_h5"
!
!-------------------------------------------------------------------------------
!
! check if the attribute exists in the group provided by gid
!
    call h5aexists_by_name_f(gid, '.', aname, exists, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot check if attribute exists :" // trim(aname))
      return
    end if
    if (.not. exists) then
      call print_error(fname, "Attribute does not exist :" // trim(aname))
      return
    end if

! open the attribute
!
    call h5aopen_by_name_f(gid, '.', aname, aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot open attribute :" // trim(aname))
      return
    end if

! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_INTEGER, avalue, am(:), ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot read attribute :" // trim(aname))
    end if

! close the attribute
!
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot close attribute :" // trim(aname))
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_scalar_attribute_integer_h5
!
!===============================================================================
!
! subroutine READ_SCALAR_ATTRIBUTE_DOUBLE_H5:
! ------------------------------------------
!
!   Subroutine reads a value of the double precision attribute provided by
!   the group identifier to which it is linked and its name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute is linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine read_scalar_attribute_double_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5aexists_by_name_f
    use hdf5           , only : h5aopen_by_name_f, h5aread_f, h5aclose_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)  , intent(in)    :: gid
    character(len=*), intent(in)    :: aname
    real(kind=8)    , intent(inout) :: avalue

! local variables
!
    logical                         :: exists = .false.
    integer(hid_t)                  :: aid
    integer(hsize_t), dimension(1)  :: am     = (/ 1 /)
    integer                         :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                         "io::read_scalar_attribute_double_h5"
!
!-------------------------------------------------------------------------------
!
! check if the attribute exists in the group provided by gid
!
    call h5aexists_by_name_f(gid, '.', aname, exists, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot check if attribute exists :" // trim(aname))
      return
    end if
    if (.not. exists) then
      call print_error(fname, "Attribute does not exist :" // trim(aname))
      return
    end if

! open the attribute
!
    call h5aopen_by_name_f(gid, '.', aname, aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot open attribute :" // trim(aname))
      return
    end if

! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_DOUBLE, avalue, am(:), ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot read attribute :" // trim(aname))
    end if

! close the attribute
!
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot close attribute :" // trim(aname))
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_scalar_attribute_double_h5
!
!===============================================================================
!
! subroutine READ_VECTOR_ATTRIBUTE_INTEGER_H5:
! -------------------------------------------
!
!   Subroutine reads a vector of the integer attribute provided by the group
!   identifier to which it is linked and its name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute is linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine read_vector_attribute_integer_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5aexists_by_name_f, h5aget_space_f
    use hdf5           , only : h5aopen_by_name_f, h5aread_f, h5aclose_f
    use hdf5           , only : h5sclose_f, h5sget_simple_extent_dims_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: aname
    integer(kind=4) , dimension(:), intent(inout) :: avalue

! local variables
!
    logical                         :: exists = .false.
    integer(hid_t)                  :: aid, sid
    integer(hsize_t), dimension(1)  :: am, bm
    integer(hsize_t)                :: alen
    integer                         :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                        "io::read_vector_attribute_integer_h5"
!
!-------------------------------------------------------------------------------
!
! check if the attribute exists in the group provided by gid
!
    call h5aexists_by_name_f(gid, '.', aname, exists, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot check if attribute exists :" // trim(aname))
      return
    end if
    if (.not. exists) then
      call print_error(fname, "Attribute does not exist :" // trim(aname))
      return
    end if

! open the attribute
!
    call h5aopen_by_name_f(gid, '.', aname, aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot open attribute :" // trim(aname))
      return
    end if

! get the attribute space
!
    call h5aget_space_f(aid, sid, ierr)
    if (ierr == 0) then
      call h5sget_simple_extent_dims_f(sid, am, bm, ierr)
      if (ierr /= 1) then
        call print_error(fname                                                 &
                         , "Cannot get attribute dimensions :" // trim(aname))
      end if
      call h5sclose_f(sid, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                        , "Cannot close the attribute space :" // trim(aname))
      end if
    else
      call print_error(fname                                                   &
                          , "Cannot get the attribute space :" // trim(aname))
      return
    end if

! check if the output array is large enough
!
    if (am(1) > size(avalue)) then
      call print_error(fname                                                   &
                 , "Attribute too large for output argument :" // trim(aname))
      return
    end if

! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_INTEGER, avalue, am(:), ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot read attribute :" // trim(aname))
    end if

! close the attribute
!
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot close attribute :" // trim(aname))
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_vector_attribute_integer_h5
!
!===============================================================================
!
! subroutine READ_VECTOR_ATTRIBUTE_DOUBLE_H5:
! ------------------------------------------
!
!   Subroutine reads a vector of the double precision attribute provided by
!   the group identifier to which it is linked and its name.
!
!   Arguments:
!
!     gid    - the group identifier to which the attribute is linked;
!     aname  - the attribute name;
!     avalue - the attribute value;
!
!===============================================================================
!
  subroutine read_vector_attribute_double_h5(gid, aname, avalue)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5aexists_by_name_f, h5aget_space_f
    use hdf5           , only : h5aopen_by_name_f, h5aread_f, h5aclose_f
    use hdf5           , only : h5sclose_f, h5sget_simple_extent_dims_f

! local variables are not implicit by default
!
    implicit none

! attribute arguments
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: aname
    real(kind=8)    , dimension(:), intent(inout) :: avalue

! local variables
!
    logical                         :: exists = .false.
    integer(hid_t)                  :: aid, sid
    integer(hsize_t), dimension(1)  :: am, bm
    integer(hsize_t)                :: alen
    integer                         :: ierr

! subroutine name string
!
    character(len=*), parameter :: fname =                                     &
                                         "io::read_vector_attribute_double_h5"
!
!-------------------------------------------------------------------------------
!
! check if the attribute exists in the group provided by gid
!
    call h5aexists_by_name_f(gid, '.', aname, exists, ierr)
    if (ierr /= 0) then
      call print_error(fname                                                   &
                        , "Cannot check if attribute exists :" // trim(aname))
      return
    end if
    if (.not. exists) then
      call print_error(fname, "Attribute does not exist :" // trim(aname))
      return
    end if

! open the attribute
!
    call h5aopen_by_name_f(gid, '.', aname, aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot open attribute :" // trim(aname))
      return
    end if

! get the attribute space
!
    call h5aget_space_f(aid, sid, ierr)
    if (ierr == 0) then
      call h5sget_simple_extent_dims_f(sid, am, bm, ierr)
      if (ierr /= 1) then
        call print_error(fname                                                 &
                         , "Cannot get attribute dimensions :" // trim(aname))
      end if
      call h5sclose_f(sid, ierr)
      if (ierr /= 0) then
        call print_error(fname                                                 &
                        , "Cannot close the attribute space :" // trim(aname))
      end if
    else
      call print_error(fname                                                   &
                          , "Cannot get the attribute space :" // trim(aname))
      return
    end if

! check if the output array is large enough
!
    if (am(1) > size(avalue)) then
      call print_error(fname                                                   &
                 , "Attribute too large for output argument :" // trim(aname))
      return
    end if

! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_DOUBLE, avalue, am(:), ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot read attribute :" // trim(aname))
    end if

! close the attribute
!
    call h5aclose_f(aid, ierr)
    if (ierr /= 0) then
      call print_error(fname, "Cannot close attribute :" // trim(aname))
      return
    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_vector_attribute_double_h5
!
!===============================================================================
!
! WRITE_ARRAY SUBROUTINES
!
!===============================================================================
!
! subroutine WRITE_1D_ARRAY_INTEGER_H5:
! ------------------------------------
!
!   Subroutine stores a one-dimensional integer array in a group specified by
!   identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_1d_array_integer_h5(gid, name, ln, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_INTEGER
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: name
    integer(hsize_t)              , intent(in) :: ln
    integer(kind=4) , dimension(:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! array dimensions
!
    integer(hsize_t), dimension(1) :: dm

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_1d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! substitute array dimensions
!
    dm(1) = ln

! create a space for the array
!
    call h5screate_simple_f(1, dm(1:1), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 1, dm(1:1), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (dm >= 32) call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:), dm(1:1), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_1d_array_integer_h5
!
!===============================================================================
!
! subroutine WRITE_2D_ARRAY_INTEGER_H5:
! ------------------------------------
!
!   Subroutine stores a two-dimensional integer array in a group specified by
!   identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_2d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_INTEGER
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                  , intent(in) :: gid
    character(len=*)                , intent(in) :: name
    integer(hsize_t), dimension(2)  , intent(in) :: dm
    integer(kind=4) , dimension(:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_2d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(2, dm(1:2), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 2, dm(1:2), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:2)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:), dm(1:2), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_2d_array_integer_h5
!
!===============================================================================
!
! subroutine WRITE_3D_ARRAY_INTEGER_H5:
! ------------------------------------
!
!   Subroutine stores a three-dimensional integer array in a group specified by
!   identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_3d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_INTEGER
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                    , intent(in) :: gid
    character(len=*)                  , intent(in) :: name
    integer(hsize_t), dimension(3)    , intent(in) :: dm
    integer(kind=4) , dimension(:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_3d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(3, dm(1:3), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 3, dm(1:3), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:3)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:,:), dm(1:3), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_3d_array_integer_h5
!
!===============================================================================
!
! subroutine WRITE_4D_ARRAY_INTEGER_H5:
! ------------------------------------
!
!   Subroutine stores a four-dimensional integer array in a group specified by
!   identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_4d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_INTEGER
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                      , intent(in) :: gid
    character(len=*)                    , intent(in) :: name
    integer(hsize_t), dimension(4)      , intent(in) :: dm
    integer(kind=4) , dimension(:,:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_4d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(4, dm(1:4), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 4, dm(1:4), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:4)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:), dm(1:4), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_4d_array_integer_h5
!
!===============================================================================
!
! subroutine WRITE_5D_ARRAY_INTEGER_H5:
! ------------------------------------
!
!   Subroutine stores a five-dimensional integer array in a group specified by
!   identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_5d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_INTEGER
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                        , intent(in) :: gid
    character(len=*)                      , intent(in) :: name
    integer(hsize_t), dimension(5)        , intent(in) :: dm
    integer(kind=4) , dimension(:,:,:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_5d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(5, dm(1:5), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 5, dm(1:5), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:5)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:,:), dm(1:5)         &
                                                                  , iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_5d_array_integer_h5
!
!===============================================================================
!
! subroutine WRITE_1D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine stores a one-dimensional double precision array in a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_1d_array_double_h5(gid, name, ln, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_DOUBLE
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: name
    integer(hsize_t)              , intent(in) :: ln
    real(kind=8)    , dimension(:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! array dimensions
!
    integer(hsize_t), dimension(1) :: dm

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_1d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! substitute array dimensions
!
    dm(1) = ln

! create a space for the array
!
    call h5screate_simple_f(1, dm(1:1), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 1, dm(1:1), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (dm >= 32) call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:), dm(1:1), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_1d_array_double_h5
!
!===============================================================================
!
! subroutine WRITE_2D_ARRAY_DOUBLE_H5:
! ------------------------------------
!
!   Subroutine stores a two-dimensional double precision array in a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_2d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_DOUBLE
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                  , intent(in) :: gid
    character(len=*)                , intent(in) :: name
    integer(hsize_t), dimension(2)  , intent(in) :: dm
    real(kind=8)    , dimension(:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_2d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(2, dm(1:2), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 2, dm(1:2), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:2)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:), dm(1:2), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_2d_array_double_h5
!
!===============================================================================
!
! subroutine WRITE_3D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine stores a three-dimensional double precision array in a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_3d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_DOUBLE
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                    , intent(in) :: gid
    character(len=*)                  , intent(in) :: name
    integer(hsize_t), dimension(3)    , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_3d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(3, dm(1:3), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 3, dm(1:3), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:3)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:), dm(1:3), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_3d_array_double_h5
!
!===============================================================================
!
! subroutine WRITE_4D_ARRAY_DOUBLE_H5:
! ------------------------------------
!
!   Subroutine stores a four-dimensional double precision array in a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_4d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_DOUBLE
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                      , intent(in) :: gid
    character(len=*)                    , intent(in) :: name
    integer(hsize_t), dimension(4)      , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_4d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(4, dm(1:4), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 4, dm(1:4), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:4)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:), dm(1:4), iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_4d_array_double_h5
!
!===============================================================================
!
! subroutine WRITE_5D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine stores a five-dimensional double precision array in a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine write_5d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error, print_warning
    use hdf5           , only : H5T_NATIVE_DOUBLE
#ifdef COMPRESS
    use hdf5           , only : H5P_DATASET_CREATE_F
#ifdef SZIP
    use hdf5           , only : H5_SZIP_NN_OM_F
#endif /* SZIP */
#endif /* COMPRESS */
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5screate_simple_f, h5sclose_f
    use hdf5           , only : h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5           , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5           , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5           , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                        , intent(in) :: gid
    character(len=*)                      , intent(in) :: name
    integer(hsize_t), dimension(5)        , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:,:,:), intent(in) :: var

#ifdef COMPRESS
! test for compression
!
    logical        :: compress = .false.
#endif /* COMPRESS */

! HDF5 object identifiers
!
    integer(hid_t) :: sid, pid, did

! procedure return value
!
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::write_5d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! create a space for the array
!
    call h5screate_simple_f(5, dm(1:5), sid, iret)

! check if the space has been created successfuly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the space
!
      call print_error(fname, "Cannot create space for dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

#ifdef COMPRESS
! prepare the property object for compression
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, iret)

! check if the object has been created properly, if not quit
!
    if (iret < 0) then

! print error about the problem with creating the compression property
!
      call print_error(fname, "Cannot create property for dataset: "           &
                                                                // trim(name))

! quit the subroutine
!
      return

    end if

! so far ok, so turn on the compression
!
    compress = .true.

! set the chunk size
!
    call h5pset_chunk_f(pid, 5, dm(1:5), iret)

! check if the chunk size has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the chunk size
!
      call print_warning(fname, "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
      compress = .false.

    end if

! set the compression algorithm
!
#ifdef DEFLATE
    call h5pset_deflate_f(pid, 9, iret)
#endif /* DEFLATE */
#ifdef SZIP
    if (product(dm(1:5)) >= 32)                                                &
                            call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, iret)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
    if (iret > 0) then

! print error about the problem with setting the compression method
!
      call print_warning(fname, "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
      compress = .false.

    end if

! check if it is safe to use compression
!
    if (compress) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret, pid)

    else

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)

    end if
#else /* COMPRESS */
! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, iret)
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
    if (iret >= 0) then

! write the dataset data
!
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:,:), dm(1:5)          &
                                                                  , iret, sid)

! check if the dataset has been written successfuly
!
      if (iret > 0) then

! print error about the problem with writing down the dataset
!
        call print_error(fname, "Cannot write dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
      if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the dataset
!
      call print_error(fname, "Cannot create dataset: " // trim(name))

    end if

! release the space
!
    call h5sclose_f(sid, iret)

! check if the space has been released successfuly
!
    if (iret > 0) then

! print error about the problem with closing the space
!
      call print_error(fname, "Cannot close space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_5d_array_double_h5

!===============================================================================
!
! READ_ARRAY SUBROUTINES
!
!===============================================================================
!
! subroutine READ_1D_ARRAY_INTEGER_H5:
! -----------------------------------
!
!   Subroutine restores a one-dimensional integer array from a group specified
!   by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_1d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: name
    integer(hsize_t), dimension(1), intent(inout) :: dm
    integer(kind=4) , dimension(:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_1d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_INTEGER, var(:), dm(1:1), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_1d_array_integer_h5
!
!===============================================================================
!
! subroutine READ_2D_ARRAY_INTEGER_H5:
! -----------------------------------
!
!   Subroutine restores a two-dimensional integer array from a group specified
!   by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_2d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                  , intent(in)    :: gid
    character(len=*)                , intent(in)    :: name
    integer(hsize_t), dimension(2)  , intent(inout) :: dm
    integer(kind=4) , dimension(:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_2d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:), dm(1:2), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_2d_array_integer_h5
!
!===============================================================================
!
! subroutine READ_3D_ARRAY_INTEGER_H5:
! -----------------------------------
!
!   Subroutine restores a three-dimensional integer array from a group specified
!   by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_3d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                    , intent(in)    :: gid
    character(len=*)                  , intent(in)    :: name
    integer(hsize_t), dimension(3)    , intent(inout) :: dm
    integer(kind=4) , dimension(:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_3d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:,:), dm(1:3), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_3d_array_integer_h5
!
!===============================================================================
!
! subroutine READ_4D_ARRAY_INTEGER_H5:
! -----------------------------------
!
!   Subroutine restores a four-dimensional integer array from a group specified
!   by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_4d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                      , intent(in)    :: gid
    character(len=*)                    , intent(in)    :: name
    integer(hsize_t), dimension(4)      , intent(inout) :: dm
    integer(kind=4) , dimension(:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_4d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:), dm(1:4), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_4d_array_integer_h5
!
!===============================================================================
!
! subroutine READ_5D_ARRAY_INTEGER_H5:
! -----------------------------------
!
!   Subroutine restores a five-dimensional integer array from a group specified
!   by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_5d_array_integer_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_INTEGER
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                        , intent(in)    :: gid
    character(len=*)                      , intent(in)    :: name
    integer(hsize_t), dimension(5)        , intent(inout) :: dm
    integer(kind=4) , dimension(:,:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_5d_array_integer_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:,:), dm(1:5), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_5d_array_integer_h5
!
!===============================================================================
!
! subroutine READ_1D_ARRAY_DOUBLE_H5:
! ----------------------------------
!
!   Subroutine restores a one-dimensional double precision array from a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_1d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: name
    integer(hsize_t), dimension(1), intent(inout) :: dm
    real(kind=8)    , dimension(:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_1d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:), dm(1:1), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_1d_array_double_h5
!
!===============================================================================
!
! subroutine READ_2D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine restores a two-dimensional double precision array from a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_2d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                  , intent(in)    :: gid
    character(len=*)                , intent(in)    :: name
    integer(hsize_t), dimension(2)  , intent(inout) :: dm
    real(kind=8)    , dimension(:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_2d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:,:), dm(1:2), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_2d_array_double_h5
!
!===============================================================================
!
! subroutine READ_3D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine restores a three-dimensional double precision array from a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_3d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                    , intent(in)    :: gid
    character(len=*)                  , intent(in)    :: name
    integer(hsize_t), dimension(3)    , intent(inout) :: dm
    real(kind=8)    , dimension(:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_3d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:,:,:), dm(1:3), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_3d_array_double_h5
!
!===============================================================================
!
! subroutine READ_4D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine restores a four-dimensional double precision array from a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_4d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                      , intent(in)    :: gid
    character(len=*)                    , intent(in)    :: name
    integer(hsize_t), dimension(4)      , intent(inout) :: dm
    real(kind=8)    , dimension(:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_4d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:), dm(1:4), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_4d_array_double_h5
!
!===============================================================================
!
! subroutine READ_5D_ARRAY_DOUBLE_H5:
! -----------------------------------
!
!   Subroutine restores a five-dimensional double precision array from a group
!   specified by identifier.
!
!   Arguments:
!
!     gid   - the HDF5 group identifier
!     name  - the string name describing the array
!     dm    - the array dimensions
!     value - the array values
!
!===============================================================================
!
  subroutine read_5d_array_double_h5(gid, name, dm, var)

! import procedures and variables from other modules
!
    use error          , only : print_error
    use hdf5           , only : H5T_NATIVE_DOUBLE
    use hdf5           , only : hid_t, hsize_t
    use hdf5           , only : h5dopen_f, h5dread_f, h5dclose_f

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer(hid_t)                        , intent(in)    :: gid
    character(len=*)                      , intent(in)    :: name
    integer(hsize_t), dimension(5)        , intent(inout) :: dm
    real(kind=8)    , dimension(:,:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: iret

! subroutine name string
!
    character(len=*), parameter :: fname = "io::read_5d_array_double_h5"
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, iret)

! check if the dataset has been opened successfuly
!
    if (iret < 0) then

! print error about the problem with opening the data space
!
      call print_error(fname, "Cannot open dataset: " // trim(name))

! quit the subroutine
!
      return

    end if

! read dataset data
!
    call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:,:), dm(1:5), iret)

! check if the dataset has been read successfuly
!
    if (iret > 0) then

! print error about the problem with reading the dataset
!
      call print_error(fname, "Cannot read dataset: " // trim(name))

    end if

! close the dataset
!
    call h5dclose_f(did, iret)

! check if the dataset has been closed successfuly
!
    if (iret > 0) then

! print error about the problem with closing the dataset
!
        call print_error(fname, "Cannot close dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_5d_array_double_h5
#endif /* HDF5 */

!===============================================================================
!
end module
