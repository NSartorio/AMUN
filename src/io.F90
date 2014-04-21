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
  real              , save :: hrest       =  6.0e+00
  real              , save :: hsnap       =  1.0e+00
  real              , save :: tsnap       =  0.0e+00

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
    real   , intent(in)  :: thrs
    integer, intent(in)  :: nrun
    integer, intent(out) :: iret
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
!   todether and stored in big 4D arrays separately.  This is a wrapper for
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
  real function next_tout()

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

! opent the HDF5 file
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
! write_attributes_h5: subroutine writes attributes in the HDF5 format
!                      connected to the provided identificator
!
! info: this subroutine stores only the global attributes
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine write_attributes_h5(fid)

! references to other modules
!
    use blocks   , only : get_mblocks, get_dblocks, get_nleafs
    use blocks   , only : get_last_id
    use coordinates, only : nn, ng, in, jn, kn, minlev, maxlev, toplev, ir, jr, kr
    use coordinates, only : xmin, xmax, ymin, ymax, zmin, zmax
    use error    , only : print_error
    use evolution, only : step, time, dt, dtn
    use hdf5     , only : hid_t
    use hdf5     , only : h5gcreate_f, h5gclose_f
    use mpitools , only : nprocs, nproc
    use random   , only : nseeds, get_seeds

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)  :: gid
    integer(kind=4) :: dm(3)
    integer         :: err

! allocatable arrays
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
    if (err .ge. 0) then

! store the integer attributes
!
      call write_attribute_integer_h5(gid, 'ndims'  , NDIMS)
      call write_attribute_integer_h5(gid, 'last_id', get_last_id())
      call write_attribute_integer_h5(gid, 'mblocks', get_mblocks())
      call write_attribute_integer_h5(gid, 'dblocks', get_dblocks())
      call write_attribute_integer_h5(gid, 'nleafs' , get_nleafs())
      call write_attribute_integer_h5(gid, 'ncells' , nn)
      call write_attribute_integer_h5(gid, 'nghost' , ng)
      call write_attribute_integer_h5(gid, 'minlev' , minlev)
      call write_attribute_integer_h5(gid, 'maxlev' , maxlev)
      call write_attribute_integer_h5(gid, 'toplev' , toplev)
      call write_attribute_integer_h5(gid, 'nprocs' , nprocs)
      call write_attribute_integer_h5(gid, 'nproc'  , nproc)
      call write_attribute_integer_h5(gid, 'nseeds' , nseeds)
      call write_attribute_integer_h5(gid, 'step'   , step  )
      call write_attribute_integer_h5(gid, 'isnap'  , isnap )

! store the real attributes
!
      call write_attribute_double_h5(gid, 'xmin', xmin)
      call write_attribute_double_h5(gid, 'xmax', xmax)
      call write_attribute_double_h5(gid, 'ymin', ymin)
      call write_attribute_double_h5(gid, 'ymax', ymax)
      call write_attribute_double_h5(gid, 'zmin', zmin)
      call write_attribute_double_h5(gid, 'zmax', zmax)
      call write_attribute_double_h5(gid, 'time', time)
      call write_attribute_double_h5(gid, 'dt'  , dt  )
      call write_attribute_double_h5(gid, 'dtn' , dtn )

! store the vector attributes
!
      dm(:) = (/ in, jn, kn /)
      call write_attribute_vector_integer_h5(gid, 'dims' , 3, dm(:))
      call write_attribute_vector_integer_h5(gid, 'rdims', 3, (/ ir, jr, kr /))

! store random number generator seed values
!
      if (nseeds .gt. 0) then

! allocate space for seeds
!
        allocate(seeds(nseeds))

! get the seed values
!
        call get_seeds(seeds)

! store them in the current group
!
        call write_attribute_vector_integer_h5(gid, 'seeds', nseeds, seeds(:))

! deallocate seed array
!
        deallocate(seeds)

      end if ! nseeds > 0

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_attributes_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_attributes_h5", "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_attributes_h5
!
!===============================================================================
!
! read_attributes_h5: subroutine restores attributes from an HDF5 file linked
!                     to the HDF5 file identificator
!
! info: this subroutine restores only the global attributes
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine read_attributes_h5(fid)

! references to other modules
!
    use blocks   , only : block_meta, block_data
    use blocks   , only : append_metablock
    use blocks   , only : set_last_id, get_last_id, get_mblocks, get_dblocks   &
                        , get_nleafs
    use coordinates, only : nn, ng, in, jn, kn, maxlev, toplev, ir, jr, kr
    use coordinates, only : initialize_coordinates, finalize_coordinates
    use coordinates, only : xmin, xmax, ymin, ymax, zmin, zmax
    use error    , only : print_error
    use evolution, only : step, time, dt, dtn
    use hdf5     , only : hid_t, hsize_t
    use hdf5     , only : h5gopen_f, h5gclose_f, h5aget_num_attrs_f            &
                        , h5aopen_idx_f, h5aclose_f, h5aget_name_f
    use mpitools , only : nprocs, nproc
    use random   , only : nseeds, set_seeds

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    character(len=16) :: aname
    integer(hid_t)    :: gid, aid
    integer(hsize_t)  :: alen = 16
    integer(kind=4)   :: dm(3)
    integer           :: err, i, l
    integer           :: nattrs, lndims, llast_id, lmblocks, lnleafs           &
                       , lncells, lnghost, lnseeds, lmaxlev, lnproc

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata

! allocatable arrays
!
    integer(kind=4), dimension(:), allocatable :: seeds
!
!-------------------------------------------------------------------------------
!
! open the global attributes group
!
    call h5gopen_f(fid, 'attributes', gid, err)

! check if the group has been opened successfuly
!
    if (err .ge. 0) then

! read the number of global attributes
!
      call h5aget_num_attrs_f(gid, nattrs, err)

! check if the number of attributes has been read successfuly
!
      if (err .ge. 0) then

! iterate over all attributes
!
        do i = 0, nattrs - 1

! open the current attribute
!
          call h5aopen_idx_f(gid, i, aid, err)

! check if the attribute has been opened successfuly
!
          if (err .ge. 0) then

! obtain the attribute name
!
            call h5aget_name_f(aid, alen, aname, err)

! depending on the attribute name use proper subroutine to read its value
!
            select case(trim(aname))
            case('ndims')
              call read_attribute_integer_h5(aid, aname, lndims)

! check if the restart file and compiled program have the same number of
! dimensions
!
              if (lndims .ne. NDIMS) then
                call print_error("io::read_attributes_h5"                      &
                              , "File and program dimensions are incompatible!")
              end if
            case('maxlev')
              call read_attribute_integer_h5(aid, aname, lmaxlev)
              if (lmaxlev .gt. toplev) then

! subtitute the new value of toplev
!
                toplev = lmaxlev

! regenerate coordinates
!
                call finalize_coordinates(err)
                call initialize_coordinates(.false., err)

! calculate a factor to rescale the block coordinates
!
                dcor = 2**(toplev - maxlev)

              else

! calculate a factor to rescale the block coordinates
!
                ucor = 2**(maxlev - lmaxlev)
              end if
            case('nprocs')
              call read_attribute_integer_h5(aid, aname, nfiles)
            case('nproc')
              call read_attribute_integer_h5(aid, aname, lnproc)
            case('last_id')
              call read_attribute_integer_h5(aid, aname, llast_id)
            case('mblocks')
              call read_attribute_integer_h5(aid, aname, lmblocks)
            case('nleafs')
              call read_attribute_integer_h5(aid, aname, lnleafs)
            case('ncells')
              call read_attribute_integer_h5(aid, aname, lncells)

! check if the block dimensions are compatible
!
              if (lncells .ne. nn) then
                call print_error("io::read_attributes_h5"                      &
                        , "File and program block dimensions are incompatible!")
              end if
            case('nghost')
              call read_attribute_integer_h5(aid, aname, lnghost)

! check if the ghost layers are compatible
!
              if (lnghost .ne. ng) then
                call print_error("io::read_attributes_h5"                      &
                      , "File and program block ghost layers are incompatible!")
              end if
            case('step')
              call read_attribute_integer_h5(aid, aname, step)
            case('isnap')
              call read_attribute_integer_h5(aid, aname, isnap)
            case('time')
              call read_attribute_double_h5(aid, aname, time)
            case('dt')
              call read_attribute_double_h5(aid, aname, dt)
            case('dtn')
              call read_attribute_double_h5(aid, aname, dtn)
            case('xmin')
              call read_attribute_double_h5(aid, aname, xmin)
            case('xmax')
              call read_attribute_double_h5(aid, aname, xmax)
            case('ymin')
              call read_attribute_double_h5(aid, aname, ymin)
            case('ymax')
              call read_attribute_double_h5(aid, aname, ymax)
            case('zmin')
              call read_attribute_double_h5(aid, aname, zmin)
            case('zmax')
              call read_attribute_double_h5(aid, aname, zmax)
            case('nseeds')
              call read_attribute_integer_h5(aid, aname, lnseeds)

! ! check if the numbers of seeds are compatible
! !
!               if (lnseeds .ne. nseeds) then
!                 call print_error("io::read_attributes_h5"                      &
!                 , "The number of seeds from file and program are incompatible!")
!               end if
            case('seeds')

! check if the numbers of seeds are compatible
!
              if (lnseeds .eq. nseeds) then

! allocate space for seeds
!
                allocate(seeds(nseeds))

! store them in the current group
!
                call read_attribute_vector_integer_h5(aid, aname, nseeds       &
                                                                   , seeds(:))

! set the seed values
!
                call set_seeds(nseeds, seeds(:))

! deallocate seed array
!
                deallocate(seeds)

              end if
            case default
            end select

! close the current attribute
!
            call h5aclose_f(aid, err)

          else

! print error about the problem with obtaining the number of attributes
!
            call print_error("io::read_attributes_h5",                         &
                                           "Cannot open the current attribute!")

          end if

        end do

! allocate all metablocks
!
        do l = 1, lmblocks
          call append_metablock(pmeta)
        end do

! check if the number of created metablocks is equal to lbmcloks
!
        if (lmblocks .ne. get_mblocks()) then
          call print_error("io::read_attributes_h5"                            &
                                        , "Number of metablocks doesn't match!")
        end if

! allocate an array of pointers with the size llast_id
!
        allocate(block_array(llast_id))
        call set_last_id(llast_id)

      else

! print error about the problem with obtaining the number of attributes
!
        call print_error("io::read_attributes_h5",                             &
                                  "Cannot get the number of global attributes!")

      end if

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::read_attributes_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::read_attributes_h5", "Cannot open the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_attributes_h5
!
!===============================================================================
!
! read_datablock_dims_h5: subroutine reads the data block dimensions from the
!                         attributes group of the file given by the file
!                         identificator
!
! arguments:
!   fid - the HDF5 file identificator
!   dm  - the data block dimensions
!
!===============================================================================
!
  subroutine read_datablock_dims_h5(fid, dm)

! references to other modules
!
    use error    , only : print_error
    use hdf5     , only : hid_t, hsize_t
    use hdf5     , only : h5gopen_f, h5gclose_f, h5aget_num_attrs_f            &
                        , h5aopen_idx_f, h5aclose_f, h5aget_name_f
    use equations, only : nv

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in) :: fid
    integer(hsize_t), dimension(5), intent(out) :: dm

! local variables
!
    character(len=16) :: aname
    integer(hid_t)    :: gid, aid
    integer(hsize_t)  :: alen = 16
    integer           :: err, i
    integer           :: nattrs, ldblocks, lnghost

! local arrays
!
    integer(kind=4), dimension(3) :: lm
!
!-------------------------------------------------------------------------------
!
! initiate the output vector
!
    dm(:) = 0
    dm(2) = nv

! open the global attributes group
!
    call h5gopen_f(fid, 'attributes', gid, err)

! check if the group has been opened successfuly
!
    if (err .ge. 0) then

! read the number of global attributes
!
      call h5aget_num_attrs_f(gid, nattrs, err)

! check if the number of attributes has been read successfuly
!
      if (err .ge. 0) then

! iterate over all attributes
!
        do i = 0, nattrs - 1

! open the current attribute
!
          call h5aopen_idx_f(gid, i, aid, err)

! check if the attribute has been opened successfuly
!
          if (err .ge. 0) then

! obtain the attribute name
!
            call h5aget_name_f(aid, alen, aname, err)

! check if the attribute name has been read successfuly
!
              if (err .ge. 0) then

! depending on the attribute name use proper subroutine to read its value
!
              select case(trim(aname))
              case('dblocks')

! obtain the number of data blocks
!
                call read_attribute_integer_h5(aid, aname, ldblocks)

              case('dims')

! obtain the block dimensions
!
                call read_attribute_vector_integer_h5(aid, aname, 3, lm(:))

              case('nghost')

! obtain the number of data blocks
!
                call read_attribute_integer_h5(aid, aname, lnghost)

              case default
              end select

            else

! print error about the problem with reading the attribute name
!
              call print_error("io::read_datablock_dims_h5",                   &
                                    "Cannot read the current attribute name!")

            end if

! close the current attribute
!
            call h5aclose_f(aid, err)

          else

! print error about the problem with opening the current attribute
!
            call print_error("io::read_datablock_dims_h5",                     &
                                         "Cannot open the current attribute!")

          end if

        end do ! i = 0, nattrs - 1

! prepare the output array
!
        dm(1) = ldblocks
        do i = 1, 3
          if (lm(i) .gt. 1) lm(i) = lm(i) + 2 * lnghost
        end do
        dm(3:5) = lm(1:3)

      else

! print error about the problem with obtaining the number of attributes
!
        call print_error("io::read_datablock_dims_h5",                         &
                                "Cannot get the number of global attributes!")

      end if

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::read_datablock_dims_h5"                          &
                                       , "Cannot close the attributes group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::read_datablock_dims_h5"                            &
                                        , "Cannot open the attributes group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_datablock_dims_h5
!
!===============================================================================
!
! write_attribute_integer_h5: subroutine writes an attribute with the integer
!                             value in a group given by its indentificator
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name of the attribute
!   value - the attribute value
!
!===============================================================================
!
  subroutine write_attribute_integer_h5(gid, name, value)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5acreate_f, h5awrite_f, h5aclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)  , intent(in) :: gid
    character(len=*), intent(in) :: name
    integer         , intent(in) :: value

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! create space for the attribute
!
    call h5screate_simple_f(1, am, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

! create the attribute
!
      call h5acreate_f(gid, name, H5T_NATIVE_INTEGER, sid, aid, err)

! check if the attribute has been created successfuly
!
      if (err .ge. 0) then

! write the attribute data
!
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, value, am, err)

! check if the attribute data has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with storing the attribute data
!
          call print_error("io::write_attribute_integer_h5"  &
                         , "Cannot write the attribute data in :" // trim(name))

        end if

! close the attribute
!
        call h5aclose_f(aid, err)

! check if the attribute has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the attribute
!
          call print_error("io::write_attribute_integer_h5"  &
                                     , "Cannot close attribute :" // trim(name))

        end if

      else

! print error about the problem with creating the attribute
!
        call print_error("io::write_attribute_integer_h5"  &
                                    , "Cannot create attribute :" // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_attribute_integer_h5"  &
                           , "Cannot close space for attribute :" // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_attribute_integer_h5"  &
                          , "Cannot create space for attribute :" // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_attribute_integer_h5
!
!===============================================================================
!
! read_attribute_integer_h5: subroutine read an integer value from an attribute
!                            given by the identificator
!
! arguments:
!   aid   - the HDF5 attribute identificator
!   value - the attribute value
!
!===============================================================================
!
  subroutine read_attribute_integer_h5(aid, name, value)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5aread_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)  , intent(in)    :: aid
    character(len=*), intent(in)    :: name
    integer         , intent(inout) :: value

! local variables
!
    integer(hsize_t), dimension(1)  :: am = (/ 1 /)
    integer                         :: err
!
!-------------------------------------------------------------------------------
!
! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_INTEGER, value, am(:), err)

! check if the attribute has been read successfuly
!
    if (err .gt. 0) then

! print error about the problem with reading the attribute
!
      call print_error("io::read_attribute_integer_h5"                         &
                          , "Cannot read attribute :" // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_attribute_integer_h5
!
!===============================================================================
!
! write_attribute_vector_integer_h5: subroutine writes an vector intiger
!                   attribute in a group given by its indentificator
!
! arguments:
!   gid    - the HDF5 group identificator
!   name   - the string name of the attribute
!   length - the vector length
!   data   - the attribute value
!
!===============================================================================
!
  subroutine write_attribute_vector_integer_h5(gid, name, length, data)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5acreate_f, h5awrite_f, h5aclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: name
    integer                       , intent(in) :: length
    integer(kind=4) , dimension(:), intent(in) :: data

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! prepare the space dimensions
!
    am(1) = length

! create space for the attribute
!
    call h5screate_simple_f(1, am, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

! create the attribute
!
      call h5acreate_f(gid, name, H5T_NATIVE_INTEGER, sid, aid, err)

! check if the attribute has been created successfuly
!
      if (err .ge. 0) then

! write the attribute data
!
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, data, am, err)

! check if the attribute data has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with storing the attribute data
!
          call print_error("io::write_attribute_vector_integer_h5"  &
                         , "Cannot write the attribute data in :" // trim(name))

        end if

! close the attribute
!
        call h5aclose_f(aid, err)

! check if the attribute has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the attribute
!
          call print_error("io::write_attribute_vector_integer_h5"  &
                                     , "Cannot close attribute :" // trim(name))

        end if

      else

! print error about the problem with creating the attribute
!
        call print_error("io::write_attribute_vector_integer_h5"  &
                                    , "Cannot create attribute :" // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_attribute_vector_integer_h5"  &
                           , "Cannot close space for attribute :" // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_attribute_vector_integer_h5"  &
                          , "Cannot create space for attribute :" // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_attribute_vector_integer_h5
!
!===============================================================================
!
! read_attribute_vector_integer_h5: subroutine reads a vector of integer values
!                                   from an attribute given by the identificator
!
! arguments:
!   aid    - the HDF5 attribute identificator
!   name   - the attribute name
!   length - the vector length
!   value  - the attribute value
!
!===============================================================================
!
  subroutine read_attribute_vector_integer_h5(aid, name, length, value)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5aread_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)  , intent(in)         :: aid
    character(len=*), intent(in)         :: name
    integer         , intent(in)         :: length
    integer, dimension(:), intent(inout) :: value

! local variables
!
    integer(hsize_t), dimension(1)  :: am
    integer                         :: err
!
!-------------------------------------------------------------------------------
!
! check if the length is larger than 0
!
    if (length .gt. 0) then

! prepare dimension array
!
      am(1) = length

! read attribute value
!
      call h5aread_f(aid, H5T_NATIVE_INTEGER, value(:), am(:), err)

! check if the attribute has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the attribute
!
        call print_error("io::read_attribute_vector_integer_h5"                &
                                    , "Cannot read attribute :" // trim(name))

      end if

    else ! length > 0

! print error about the wrong vector size
!
      call print_error("io::read_attribute_vector_integer_h5"                  &
                                                   , "Wrong length of vector")

    end if ! length > 0

!-------------------------------------------------------------------------------
!
  end subroutine read_attribute_vector_integer_h5
!
!===============================================================================
!
! write_attribute_double_h5: subroutine writes a double precision attribute in
!                          a group given by its indentificator
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name of the attribute
!   value - the attribute value
!
!===============================================================================
!
  subroutine write_attribute_double_h5(gid, name, value)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f, h5acreate_f, h5awrite_f  &
                    , h5aclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)  , intent(in) :: gid
    character(len=*), intent(in) :: name
    real(kind=8)    , intent(in) :: value

! local variables
!
    integer(hid_t)                 :: sid, aid
    integer(hsize_t), dimension(1) :: am = (/ 1 /)
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! create space for the attribute
!
    call h5screate_simple_f(1, am, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

! create the attribute
!
      call h5acreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, aid, err)

! check if the attribute has been created successfuly
!
      if (err .ge. 0) then

! write the attribute data
!
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, value, am, err)

! check if the attribute data has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with storing the attribute data
!
          call print_error("io::write_attribute_double_h5"  &
                         , "Cannot write the attribute data in :" // trim(name))

        end if

! close the attribute
!
        call h5aclose_f(aid, err)

! check if the attribute has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the attribute
!
          call print_error("io::write_attribute_double_h5"  &
                                     , "Cannot close attribute :" // trim(name))

        end if

      else

! print error about the problem with creating the attribute
!
        call print_error("io::write_attribute_double_h5"  &
                                    , "Cannot create attribute :" // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_attribute_double_h5"  &
                           , "Cannot close space for attribute :" // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_attribute_double_h5"  &
                          , "Cannot create space for attribute :" // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_attribute_double_h5
!
!===============================================================================
!
! read_attribute_double_h5: subroutine reads a double precision attribute
!
! arguments:
!   aid   - the HDF5 attribute identificator
!   value - the attribute value
!
!===============================================================================
!
  subroutine read_attribute_double_h5(aid, name, value)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5aread_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)  , intent(in)    :: aid
    character(len=*), intent(in)    :: name
    real(kind=8)    , intent(inout) :: value

! local variables
!
    integer(hsize_t), dimension(1)  :: am = (/ 1 /)
    integer                         :: err
!
!-------------------------------------------------------------------------------
!
! read attribute value
!
    call h5aread_f(aid, H5T_NATIVE_DOUBLE, value, am(:), err)

! check if the attribute has been read successfuly
!
    if (err .gt. 0) then

! print error about the problem with reading the attribute
!
      call print_error("io::read_attribute_double_h5"                          &
                                      , "Cannot read attribute :" // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_attribute_double_h5
!
!===============================================================================
!
! write_metablocks_h5: subroutine writes metablocks in the HDF5 format connected
!                      to the provided identificator
!
! info: this subroutine stores only the metablocks
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine write_metablocks_h5(fid)

! references to other modules
!
    use blocks  , only : block_meta, list_meta
    use blocks  , only : get_last_id, get_mblocks, nchildren, nsides, nfaces
    use error   , only : print_error
    use hdf5    , only : hid_t, hsize_t
    use hdf5    , only : h5gcreate_f, h5gclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: l, p, i, j, k
    integer                        :: err
    integer(hsize_t), dimension(1) :: am, cm
    integer(hsize_t), dimension(2) :: dm, pm
    integer(hsize_t), dimension(4) :: qm

! local allocatable arrays
!
    integer(kind=4), dimension(:)  , allocatable :: idx
    integer(kind=4), dimension(:)  , allocatable :: par, dat
    integer(kind=4), dimension(:)  , allocatable ::  id, cpu, lev, cfg, ref, lea
    real   (kind=8), dimension(:)  , allocatable :: xmn, xmx, ymn, ymx, zmn, zmx
    integer(kind=4), dimension(:,:), allocatable :: chl, pos, cor
    integer(kind=4), dimension(:,:,:,:), allocatable :: ngh

! local pointers
!
    type(block_meta), pointer :: pmeta
!
!-------------------------------------------------------------------------------
!
! create the group for metadata
!
    call h5gcreate_f(fid, 'metablocks', gid, err)

! check if the group has been created successfuly
!
    if (err .ge. 0) then

! prepate dimensions
!
      am(1) = get_mblocks()
      cm(1) = get_last_id()
      dm(1) = get_mblocks()
      dm(2) = nchildren
      pm(1) = get_mblocks()
      pm(2) = NDIMS
      qm(1) = get_mblocks()
      qm(2) = NDIMS
      qm(3) = nsides
      qm(4) = nfaces

! only process if there are some metablocks
!
      if (am(1) > 0) then

! allocate arrays to store metablocks data
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
        allocate(ngh(qm(1),qm(2),qm(3),qm(4)))

! reset vectors
!
        idx(:)       = -1
        par(:)       = -1
        dat(:)       = -1
        lea(:)       = -1
        chl(:,:)     = -1
        ngh(:,:,:,:) = -1

! iterate over all metablocks and fill in the arrays for storage
!
        l = 1
        pmeta => list_meta
        do while(associated(pmeta))

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

          do i = 1, NDIMS
            do j = 1, nsides
              do k = 1, nfaces
                if (associated(pmeta%neigh(i,j,k)%ptr))  &
                                          ngh(l,i,j,k) = pmeta%neigh(i,j,k)%ptr%id
              end do
            end do
          end do

          l = l + 1
          pmeta => pmeta%next
        end do

! store metadata in the HDF5 file
!
        call write_vector_integer_h5(gid, 'indices', cm(1), idx)
        call write_vector_integer_h5(gid, 'parent' , am(1), par)
        call write_vector_integer_h5(gid, 'data'   , am(1), dat)
        call write_vector_integer_h5(gid, 'id'     , am(1), id)
        call write_vector_integer_h5(gid, 'cpu'    , am(1), cpu)
        call write_vector_integer_h5(gid, 'level'  , am(1), lev)
        call write_vector_integer_h5(gid, 'config' , am(1), cfg)
        call write_vector_integer_h5(gid, 'refine' , am(1), ref)
        call write_vector_integer_h5(gid, 'leaf'   , am(1), lea)
        call write_vector_double_h5 (gid, 'xmin'   , am(1), xmn)
        call write_vector_double_h5 (gid, 'xmax'   , am(1), xmx)
        call write_vector_double_h5 (gid, 'ymin'   , am(1), ymn)
        call write_vector_double_h5 (gid, 'ymax'   , am(1), ymx)
        call write_vector_double_h5 (gid, 'zmin'   , am(1), zmn)
        call write_vector_double_h5 (gid, 'zmax'   , am(1), zmx)
        call write_array2_integer_h5(gid, 'child'  , dm(:), chl)
        call write_array2_integer_h5(gid, 'pos'    , pm(:), pos)
        call write_array2_integer_h5(gid, 'coord'  , pm(:), cor)
        call write_array4_integer_h5(gid, 'neigh'  , qm(:), ngh)

! deallocate allocatable arrays
!
        if (allocated(idx)) deallocate(idx)
        if (allocated(par)) deallocate(par)
        if (allocated(dat)) deallocate(dat)
        if (allocated(id) ) deallocate(id)
        if (allocated(cpu)) deallocate(cpu)
        if (allocated(lev)) deallocate(lev)
        if (allocated(cfg)) deallocate(cfg)
        if (allocated(ref)) deallocate(ref)
        if (allocated(lea)) deallocate(lea)
        if (allocated(xmn)) deallocate(xmn)
        if (allocated(xmx)) deallocate(xmx)
        if (allocated(ymn)) deallocate(ymn)
        if (allocated(ymx)) deallocate(ymx)
        if (allocated(zmn)) deallocate(zmn)
        if (allocated(zmx)) deallocate(zmx)
        if (allocated(chl)) deallocate(chl)
        if (allocated(cor)) deallocate(cor)
        if (allocated(ngh)) deallocate(ngh)

      end if ! meta blocks > 0

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_metablocks_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_metablocks_h5", "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_metablocks_h5
!
!===============================================================================
!
! read_metablocks_h5: subroutine reads metablocks from the restart HDF5 file
!                     and restores all their structure fields
!
! info: this subroutine restores metablocks only
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine read_metablocks_h5(fid)

! references to other modules
!
    use blocks  , only : block_meta, list_meta
    use blocks  , only : nchildren, nsides, nfaces
    use blocks  , only : get_mblocks
    use blocks  , only : metablock_set_id, metablock_set_process               &
                       , metablock_set_refinement, metablock_set_configuration &
                       , metablock_set_level, metablock_set_position           &
                       , metablock_set_coordinates, metablock_set_bounds       &
                       , metablock_set_leaf
    use error   , only : print_error
    use hdf5    , only : hid_t, hsize_t
    use hdf5    , only : h5gopen_f, h5gclose_f
    use mpitools, only : nprocs

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: l, p, i, j, k, lcpu
    integer                        :: err
    integer(hsize_t), dimension(1) :: am
    integer(hsize_t), dimension(2) :: dm, pm
    integer(hsize_t), dimension(4) :: qm

! local allocatable arrays
!
    integer(kind=4), dimension(:)  , allocatable :: idx
    integer(kind=4), dimension(:)  , allocatable :: par, dat
    integer(kind=4), dimension(:)  , allocatable ::  id, cpu, lev, cfg, ref, lea
    real   (kind=8), dimension(:)  , allocatable :: xmn, xmx, ymn, ymx, zmn, zmx
    integer(kind=4), dimension(:,:), allocatable :: chl, pos, cor
    integer(kind=4), dimension(:,:,:,:), allocatable :: ngh

! local pointers
!
    type(block_meta), pointer :: pmeta
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
    if (err .ge. 0) then

! prepate dimensions
!
      am(1) = get_mblocks()
      dm(1) = get_mblocks()
      dm(2) = nchildren
      pm(1) = get_mblocks()
      pm(2) = NDIMS
      qm(1) = get_mblocks()
      qm(2) = NDIMS
      qm(3) = nsides
      qm(4) = nfaces

! allocate arrays to store metablocks data
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
      allocate(ngh(qm(1),qm(2),qm(3),qm(4)))

! reset vectors
!
      par(:)       = -1
      dat(:)       = -1
      lea(:)       = -1
      chl(:,:)     = -1
      ngh(:,:,:,:) = -1

! read metablock fields from the HDF5 file
!
      call read_vector_integer_h5(gid, 'id'     , am(:), id (:))
      call read_vector_integer_h5(gid, 'cpu'    , am(:), cpu(:))
      call read_vector_integer_h5(gid, 'level'  , am(:), lev(:))
      call read_vector_integer_h5(gid, 'config' , am(:), cfg(:))
      call read_vector_integer_h5(gid, 'refine' , am(:), ref(:))
      call read_vector_integer_h5(gid, 'leaf'   , am(:), lea(:))
      call read_vector_integer_h5(gid, 'parent' , am(:), par(:))
      call read_vector_double_h5 (gid, 'xmin'   , am(:), xmn(:))
      call read_vector_double_h5 (gid, 'xmax'   , am(:), xmx(:))
      call read_vector_double_h5 (gid, 'ymin'   , am(:), ymn(:))
      call read_vector_double_h5 (gid, 'ymax'   , am(:), ymx(:))
      call read_vector_double_h5 (gid, 'zmin'   , am(:), zmn(:))
      call read_vector_double_h5 (gid, 'zmax'   , am(:), zmx(:))
      call read_array2_integer_h5(gid, 'pos'    , pm(:), pos(:,:))
      call read_array2_integer_h5(gid, 'coord'  , pm(:), cor(:,:))
      call read_array2_integer_h5(gid, 'child'  , dm(:), chl(:,:))
      call read_array4_integer_h5(gid, 'neigh'  , qm(:), ngh(:,:,:,:))

! check if the maximum level has been changed, is so, rescale block coordinates
!
      if (dcor .gt. 1) then
        cor(:,:) = cor(:,:) / dcor
      end if
      if (ucor .gt. 1) then
        cor(:,:) = cor(:,:) * ucor
      end if

! prepare the array of pointers to metablocks
!
      l = 1
      pmeta => list_meta
      do while(associated(pmeta))

        block_array(id(l))%ptr => pmeta

        call metablock_set_id           (pmeta, id (l))
        call metablock_set_process      (pmeta, min(lcpu, cpu(l)))
        call metablock_set_refinement   (pmeta, ref(l))
        call metablock_set_configuration(pmeta, cfg(l))
        call metablock_set_level        (pmeta, lev(l))
        call metablock_set_position     (pmeta, pos(l,1), pos(l,2), pos(l,3))
        call metablock_set_coordinates  (pmeta, cor(l,1), cor(l,2), cor(l,3))
        call metablock_set_bounds       (pmeta, xmn(l), xmx(l), ymn(l), ymx(l) &
                                                             , zmn(l), zmx(l))

        if (lea(l) .eq. 1) call metablock_set_leaf(pmeta)

        l = l + 1
        pmeta => pmeta%next
      end do

! iterate over all metablocks and restore pointers
!
      l = 1
      pmeta => list_meta
      do while(associated(pmeta))

        if (par(l) .gt. 0) pmeta%parent => block_array(par(l))%ptr

        do p = 1, nchildren
          if (chl(l,p) .gt. 0) then
            pmeta%child(p)%ptr => block_array(chl(l,p))%ptr
          end if
        end do

        do i = 1, NDIMS
          do j = 1, nsides
            do k = 1, nfaces
              if (ngh(l,i,j,k) .gt. 0) then
                pmeta%neigh(i,j,k)%ptr => block_array(ngh(l,i,j,k))%ptr
              end if
            end do
          end do
        end do

        l = l + 1
        pmeta => pmeta%next
      end do

! deallocate allocatable arrays
!
      if (allocated(id) ) deallocate(id )
      if (allocated(par)) deallocate(par)
      if (allocated(dat)) deallocate(dat)
      if (allocated(cpu)) deallocate(cpu)
      if (allocated(lev)) deallocate(lev)
      if (allocated(cfg)) deallocate(cfg)
      if (allocated(ref)) deallocate(ref)
      if (allocated(lea)) deallocate(lea)
      if (allocated(xmn)) deallocate(xmn)
      if (allocated(xmx)) deallocate(xmx)
      if (allocated(ymn)) deallocate(ymn)
      if (allocated(ymx)) deallocate(ymx)
      if (allocated(zmn)) deallocate(zmn)
      if (allocated(zmx)) deallocate(zmx)
      if (allocated(chl)) deallocate(chl)
      if (allocated(cor)) deallocate(cor)
      if (allocated(ngh)) deallocate(ngh)

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::read_metablocks_h5"                              &
                                              , "Cannot close metablock group!")

      end if

    else

! print error about the problem with opening the group
!
      call print_error("io::read_metablocks_h5", "Cannot open metablock group!")

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
!     fid - the HDF5 file identificator;
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

! allocate arrays to store associated meta block identificators, conserved and
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
        call write_vector_integer_h5(gid, 'meta', am(1), id)
        call write_array5_double_h5 (gid, 'uvar', dm(:), uv)
        call write_array5_double_h5 (gid, 'qvar', dm(:), qv)

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
!     fid - the HDF5 file identificator;
!
!===============================================================================
!
  subroutine read_datablocks_h5(fid)

! import external procedures and variables
!
    use blocks         , only : block_meta, block_data, list_data
    use blocks         , only : append_datablock, link_blocks
    use coordinates    , only : im, jm, km
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
    integer                        :: err

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
! get datablock array dimensions
!
    call read_datablock_dims_h5(fid, dm(:))

! open the group 'datablocks'
!
    call h5gopen_f(fid, 'datablocks', gid, err)

! check if the datablock group has been opened successfuly
!
    if (err >= 0) then

! restore data blocks only if there are any
!
      if (dm(1) > 0) then

! allocate arrays to read data
!
        allocate(id(dm(1)))
        allocate(uv(dm(1),dm(2),dm(3),dm(4),dm(5)))
        allocate(qv(dm(1),dm(2),dm(3),dm(4),dm(5)))

! read array data from the HDF5 file
!
        call read_vector_integer_h5(gid, 'meta', dm(1:1), id(:)        )
        call read_array5_double_h5 (gid, 'uvar', dm(1:5), uv(:,:,:,:,:))
        call read_array5_double_h5 (gid, 'qvar', dm(1:5), qv(:,:,:,:,:))

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

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err > 0) then

! print error about the problem with closing the group
!
        call print_error("io::read_datablocks_h5"                              &
                                           , "Cannot close data block group!")

      end if

    else

! print error about the problem with opening the group
!
      call print_error("io::read_datablocks_h5"                                &
                                            , "Cannot open data block group!")

    end if

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
!   fid - the HDF5 file identificator
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
        call write_vector_integer_h5(gid, 'levels', cm(1), lev)
        call write_vector_integer_h5(gid, 'refine', cm(1), ref)
        call write_array2_integer_h5(gid, 'coords', cm(:), cor)
        call write_array3_double_h5 (gid, 'bounds', dm(:), bnd)
        call write_vector_double_h5 (gid, 'dx'    , am(1), adx(1:maxlev))
        call write_vector_double_h5 (gid, 'dy'    , am(1), ady(1:maxlev))
        call write_vector_double_h5 (gid, 'dz'    , am(1), adz(1:maxlev))

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
!   identificator.
!
!   Arguments:
!
!     fid - the HDF5 file identificator;
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
          call write_array4_double_h5(gid, trim(pvars(n)), dm, qarr)

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
!   identificator.
!
!   Arguments:
!
!     fid - the HDF5 file identificator;
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
          call write_array4_double_h5(gid, trim(cvars(n)), dm, qarr)

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
! write_vector_integer_h5: subroutine stores a 1D integer vector in a group
!
! arguments:
!   gid    - the HDF5 group identificator
!   name   - the string name representing the dataset
!   length - the vector length
!   value  - the data
!
!===============================================================================
!
  subroutine write_vector_integer_h5(gid, name, length, data)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: name
    integer(hsize_t)              , intent(in) :: length
    integer(kind=4) , dimension(:), intent(in) :: data

! local variables
!
    integer(hid_t)                 :: sid, did
    integer(hsize_t), dimension(1) :: am
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! prepare the vector dimensions
!
    am(1) = length

! create space for the vector
!
    call h5screate_simple_f(1, am, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, err)

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, data(:), am, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_vector_integer_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_vector_integer_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_vector_integer_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_vector_integer_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_vector_integer_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_vector_integer_h5
!
!===============================================================================
!
! read_vector_integer_h5: subroutine reads a 1D integer vector
!
! arguments:
!   gid    - the HDF5 group identificator
!   name   - the string name representing the dataset
!   length - the vector length
!   value  - the data
!
!===============================================================================
!
  subroutine read_vector_integer_h5(gid, name, dm, data)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5dopen_f, h5dread_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: name
    integer(hsize_t), dimension(1), intent(inout) :: dm
    integer(kind=4) , dimension(:), intent(inout) :: data

! local variables
!
    integer(hid_t)                 :: did
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, err)

! check if the dataset has been opened successfuly
!
    if (err .ge. 0) then

! read the dataset data
!
      call h5dread_f(did, H5T_NATIVE_INTEGER, data(:), dm(:), err)

! check if the dataset has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the dataset
!
        call print_error("io::read_vector_integer_h5"                          &
                                        , "Cannot read dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the dataset
!
        call print_error("io::read_vector_integer_h5"                          &
                                       , "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with opening the dataset
!
      call print_error("io::read_vector_integer_h5"                            &
                                        , "Cannot open dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_vector_integer_h5
!
!===============================================================================
!
! write_array2_integer_h5: subroutine stores a 2D integer array in a group
!
! arguments:
!   gid - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array2_integer_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                  , intent(in) :: gid
    character(len=*)                , intent(in) :: name
    integer(hsize_t), dimension(2)  , intent(in) :: dm
    integer(kind=4) , dimension(:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(2, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 2, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array4_integer_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array4_integer_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array2_integer_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array2_integer_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array2_integer_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array2_integer_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array2_integer_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array2_integer_h5
!
!===============================================================================
!
! read_array2_integer_h5: subroutine reads a 2D integer array
!
! arguments:
!   gid - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine read_array2_integer_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5dopen_f, h5dread_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                  , intent(in)    :: gid
    character(len=*)                , intent(in)    :: name
    integer(hsize_t), dimension(2)  , intent(inout) :: dm
    integer(kind=4) , dimension(:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: err
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, err)

! check if the dataset has been opened successfuly
!
    if (err .ge. 0) then

! read dataset data
!
      call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:), dm(:), err)

! check if the dataset has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the dataset
!
        call print_error("io::read_array2_integer_h5"                          &
                                       , "Cannot read dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the dataset
!
        call print_error("io::read_array2_integer_h5"                          &
                                       , "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with opening the dataset
!
      call print_error("io::read_array2_integer_h5"                            &
                                        , "Cannot open dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_array2_integer_h5
!
!===============================================================================
!
! write_array4_integer_h5: subroutine stores a 4D integer array in a group
!
! arguments:
!   gid - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array4_integer_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                      , intent(in) :: gid
    character(len=*)                    , intent(in) :: name
    integer(hsize_t), dimension(4)      , intent(in) :: dm
    integer(kind=4) , dimension(:,:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(4, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 4, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array4_integer_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array4_integer_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_INTEGER, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array4_integer_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array4_integer_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array4_integer_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array4_integer_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array4_integer_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array4_integer_h5
!
!===============================================================================
!
! read_array4_integer_h5: subroutine reads a 4D integer array
!
! arguments:
!   gid - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine read_array4_integer_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_INTEGER
    use hdf5 , only : h5dopen_f, h5dread_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                      , intent(in)    :: gid
    character(len=*)                    , intent(in)    :: name
    integer(hsize_t), dimension(4)      , intent(inout) :: dm
    integer(kind=4) , dimension(:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: err
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, err)

! check if the dataset has been opened successfuly
!
    if (err .ge. 0) then

! read dataset data
!
      call h5dread_f(did, H5T_NATIVE_INTEGER, var(:,:,:,:), dm(:), err)

! check if the dataset has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the dataset
!
        call print_error("io::read_array4_integer_h5"                          &
                                       , "Cannot read dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the dataset
!
        call print_error("io::read_array4_integer_h5"                          &
                                       , "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with opening the dataset
!
      call print_error("io::read_array4_integer_h5"                            &
                                        , "Cannot open dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_array4_integer_h5
!
!===============================================================================
!
! write_vector_double_h5: subroutine stores a 1D double precision vector in
!                         a group
!
! arguments:
!   gid - the HDF5 group identificator
!
!===============================================================================
!
  subroutine write_vector_double_h5(gid, name, length, data)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in) :: gid
    character(len=*)              , intent(in) :: name
    integer(hsize_t)              , intent(in) :: length
    real(kind=8)    , dimension(:), intent(in) :: data

! local variables
!
    integer(hid_t)                 :: sid, did
    integer(hsize_t), dimension(1) :: am
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! prepare the vector dimensions
!
    am(1) = length

! create space for the vector
!
    call h5screate_simple_f(1, am, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

! create the dataset
!
      call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err)

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, data(:), am, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_vector_double_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_vector_double_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_vector_double_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_vector_double_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_vector_double_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_vector_double_h5
!
!===============================================================================
!
! read_vector_double_h5: subroutine reads a 1D double precision vector
!
! arguments:
!   gid    - the HDF5 group identificator
!   name   - the string name representing the dataset
!   dm     - the vector dimensions
!   value  - the data
!
!===============================================================================
!
  subroutine read_vector_double_h5(gid, name, dm, data)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5dopen_f, h5dread_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                , intent(in)    :: gid
    character(len=*)              , intent(in)    :: name
    integer(hsize_t), dimension(1), intent(inout) :: dm
    real(kind=8)    , dimension(:), intent(inout) :: data

! local variables
!
    integer(hid_t)                 :: did
    integer                        :: err
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, err)

! check if the dataset has been opened successfuly
!
    if (err .ge. 0) then

! read the dataset data
!
      call h5dread_f(did, H5T_NATIVE_DOUBLE, data(:), dm(:), err)

! check if the dataset has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the dataset
!
        call print_error("io::read_vector_double_h5"                           &
                                        , "Cannot read dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the dataset
!
        call print_error("io::read_vector_double_h5"                           &
                                       , "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with opening the dataset
!
      call print_error("io::read_vector_double_h5"                             &
                                        , "Cannot open dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_vector_double_h5
!
!===============================================================================
!
! write_array4_float_h5: subroutine stores a 4D single precision array
!
! arguments:
!   gid   - the HDF5 group identificator where the dataset should be located
!   name  - the string name representing the dataset
!   dm    - the dataset dimensions
!   value - the dataset values
!
!===============================================================================
!
  subroutine write_array4_float_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error, print_warning
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_REAL
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */

! define default variables
!
    implicit none

! input variables
!
    integer(hid_t)                      , intent(in) :: gid
    character(len=*)                    , intent(in) :: name
    integer(hsize_t), dimension(4)      , intent(in) :: dm
    real(kind=4)    , dimension(:,:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
!
!-------------------------------------------------------------------------------
!
! create a space for the dataset dimensions
!
    call h5screate_simple_f(4, dm(:), sid, err)

! print an error, if the space for dimensions couldn't be created
!
    if (err .eq. -1) call print_error("io::write_array4_float_h5"              &
                    , "Cannot create a space for the dataset: " // trim(name))

! prepare the compression properties
!
    call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! if the compression properties could be created properly, set the compression
! algorithm and strength
!
    if (err .eq. 0) then

! set the chunk size
!
      call h5pset_chunk_f(pid, 4, dm(:), err)

! print a warning, if the chunk size couldn't be set properly
!
      if (err .eq. -1) call print_warning("io::write_array4_float_h5"          &
                                            , "Cannot set the size of chunk!")

! set the compression algorithm
!
#ifdef DEFLATE
      call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
      if (product(dm) .ge. 32)                                                 &
        call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! print a warning, if the compression algorithm couldn't be set
!
      if (err .eq. -1) call print_warning("io::write_array4_float_h5"          &
                                       , "Cannot set the compression method!")

    else

! print a warning, if the property list couldn't be created
!
      call print_warning("io::write_array4_float_h5"                           &
                                           , "Cannot create a property list!")

    end if

! create the dataset
!
    call h5dcreate_f(gid, name, H5T_NATIVE_REAL, sid, did, err, pid)

! print an error, if the dataset couldn't be created
!
    if (err .eq. -1) call print_error("io::write_array4_float_h5"              &
                                , "Cannot create the dataset: " // trim(name))

! write the dataset values
!
    call h5dwrite_f(did, H5T_NATIVE_REAL, var(:,:,:,:), dm, err, sid)

! print an error, if the dataset couldn't be written successfuly
!
    if (err .eq. -1) call print_error("io::write_array4_float_h5"              &
                                 , "Cannot write the dataset: " // trim(name))

! close the dataset
!
    call h5dclose_f(did, err)

! print an error, if the dataset couldn't be closed
!
    if (err .eq. -1) call print_error("io::write_array4_float_h5"              &
                                 , "Cannot close the dataset: " // trim(name))

! if the property list is created
!
    if (pid .ne. -1) then

! terminate access to the property list
!
      call h5pclose_f(pid, err)

! print a warning, if the property list couldn't be closed
!
      if (err .eq. -1)  call print_warning("io::write_array4_float_h5"         &
                                          , "Cannot close the property list!")

    end if

! release the dataspace of the current dataset
!
    call h5sclose_f(sid, err)

! print an error, if the space couldn't be released successfuly
!
    if (err .eq. -1) call print_error("io::write_array4_float_h5"              &
                   , "Cannot close the space for the dataset: " // trim(name))

!-------------------------------------------------------------------------------
!
  end subroutine write_array4_float_h5
!
!===============================================================================
!
! write_array3_double_h5: subroutine stores a 3D double precision array
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array3_double_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                    , intent(in) :: gid
    character(len=*)                  , intent(in) :: name
    integer(hsize_t), dimension(3)    , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(3, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 3, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array3_double_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array3_double_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array3_double_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array3_double_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array3_double_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array3_double_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array3_double_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array3_double_h5
!
!===============================================================================
!
! write_array4_double_h5: subroutine stores a 4D double precision array
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array4_double_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                      , intent(in) :: gid
    character(len=*)                    , intent(in) :: name
    integer(hsize_t), dimension(4)      , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(4, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 4, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array4_double_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array4_double_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array4_double_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array4_double_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array4_double_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array4_double_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array4_double_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array4_double_h5
!
!===============================================================================
!
! write_array5_double_h5: subroutine stores a 5D double precision array
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array5_double_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                        , intent(in) :: gid
    character(len=*)                      , intent(in) :: name
    integer(hsize_t), dimension(5)        , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(5, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 5, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array5_double_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array5_double_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array5_double_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array5_double_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array5_double_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! release the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array5_double_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array5_double_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array5_double_h5
!
!===============================================================================
!
! read_array5_double_h5: subroutine reads a 5D double precision array
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine read_array5_double_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5dopen_f, h5dread_f, h5dclose_f

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                        , intent(in)    :: gid
    character(len=*)                      , intent(in)    :: name
    integer(hsize_t), dimension(5)        , intent(inout) :: dm
    real(kind=8)    , dimension(:,:,:,:,:), intent(inout) :: var

! local variables
!
    integer(hid_t) :: did
    integer        :: err
!
!-------------------------------------------------------------------------------
!
! open the dataset
!
    call h5dopen_f(gid, name, did, err)

! check if the dataset has been opened successfuly
!
    if (err .ge. 0) then

! read dataset
!
      call h5dread_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:,:), dm(:), err)

! check if the dataset has been read successfuly
!
      if (err .gt. 0) then

! print error about the problem with reading the dataset
!
        call print_error("io::read_array5_double_h5"                           &
                                        , "Cannot read dataset: " // trim(name))

      end if

! close the dataset
!
      call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the dataset
!
        call print_error("io::read_array5_double_h5"                           &
                                       , "Cannot close dataset: " // trim(name))

      end if

    else

! print error about the problem with opening the dataset
!
      call print_error("io::read_array5_double_h5"                             &
                                        , "Cannot open dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_array5_double_h5
!
!===============================================================================
!
! write_array6_double_h5: subroutine stores a 6D double precision array
!
! arguments:
!   gid   - the HDF5 group identificator
!   name  - the string name representing the dataset
!   dm    - the data dimensions
!   value - the data
!
!===============================================================================
!
  subroutine write_array6_double_h5(gid, name, dm, var)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, hsize_t, H5T_NATIVE_DOUBLE
    use hdf5 , only : h5screate_simple_f, h5sclose_f                           &
                    , h5dcreate_f, h5dwrite_f, h5dclose_f
#ifdef COMPRESS
    use hdf5 , only : H5P_DATASET_CREATE_F
    use hdf5 , only : h5pcreate_f, h5pset_chunk_f, h5pclose_f
#ifdef DEFLATE
    use hdf5 , only : h5pset_deflate_f
#endif /* DEFLATE */
#ifdef SZIP
    use hdf5 , only : H5_SZIP_NN_OM_F
    use hdf5 , only : h5pset_szip_f
#endif /* SZIP */
#endif /* COMPRESS */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t)                          , intent(in) :: gid
    character(len=*)                        , intent(in) :: name
    integer(hsize_t), dimension(6)          , intent(in) :: dm
    real(kind=8)    , dimension(:,:,:,:,:,:), intent(in) :: var

! local variables
!
    integer(hid_t) :: sid, pid, did
    integer        :: err
#ifdef COMPRESS
    logical        :: compress = .false.
#endif /* COMPRESS */
!
!-------------------------------------------------------------------------------
!
! create space for the vector
!
    call h5screate_simple_f(6, dm, sid, err)

! check if the space has been created successfuly
!
    if (err .ge. 0) then

#ifdef COMPRESS
! prepare compression
!
      call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)

! check if the properties have been created properly
!
      if (err .ge. 0) then

! so far ok, so turn on the compression
!
        compress = .true.

! set the chunk size
!
        call h5pset_chunk_f(pid, 6, dm, err)

! check if the chunk size has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the chunk size
!
          call print_error("io::write_array6_double_h5"  &
                                          , "Cannot set the size of the chunk!")

! setting the size of the chunk failed, so turn off the compression
!
          compress = .false.

        end if

! set the compression algorithm
!
#ifdef DEFLATE
        call h5pset_deflate_f(pid, 9, err)
#endif /* DEFLATE */
#ifdef SZIP
        if (product(dm) .ge. 32)                                               &
          call h5pset_szip_f(pid, H5_SZIP_NN_OM_F, 32, err)
#endif /* SZIP */

! check if the compression algorithm has been set properly
!
        if (err .gt. 0) then

! print error about the problem with setting the compression method
!
          call print_error("io::write_array6_double_h5"  &
                                         , "Cannot set the compression method!")

! setting compression method failed, so turn off the compression
!
          compress = .false.

        end if

      end if

! check if it is safe to use compression
!
      if (compress) then

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err, pid)

      else
#endif /* COMPRESS */

! create the dataset
!
        call h5dcreate_f(gid, name, H5T_NATIVE_DOUBLE, sid, did, err)

#ifdef COMPRESS
      end if
#endif /* COMPRESS */

! check if the dataset has been created successfuly
!
      if (err .ge. 0) then

! write the dataset data
!
        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, var(:,:,:,:,:,:), dm, err, sid)

! check if the dataset has been written successfuly
!
        if (err .gt. 0) then

! print error about the problem with writing down the dataset
!
          call print_error("io::write_array6_double_h5"  &
                                       , "Cannot write dataset: " // trim(name))

        end if

! close the dataset
!
        call h5dclose_f(did, err)

! check if the dataset has been closed successfuly
!
        if (err .gt. 0) then

! print error about the problem with closing the dataset
!
          call print_error("io::write_array6_double_h5"  &
                                       , "Cannot close dataset: " // trim(name))

        end if

      else

! print error about the problem with creating the dataset
!
        call print_error("io::write_array6_double_h5"  &
                                      , "Cannot create dataset: " // trim(name))

      end if

! close the space
!
      call h5sclose_f(sid, err)

! check if the space has been released successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the space
!
        call print_error("io::write_array6_double_h5"  &
                             , "Cannot close space for dataset: " // trim(name))

      end if

    else

! print error about the problem with creating the space for the attribute
!
      call print_error("io::write_array6_double_h5"  &
                            , "Cannot create space for dataset: " // trim(name))

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_array6_double_h5
#endif /* HDF5 */

!===============================================================================
!
end module
