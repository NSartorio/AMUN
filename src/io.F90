!!******************************************************************************
!!
!! module: io - handling the data input and output for storing the snapshots
!!              and restarting the jobs
!!
!! Copyright (C) 2008-2011 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of AMUN.
!!
!!  AMUN is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  AMUN is distributed in the hope that it will be useful,
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
module io

  use blocks, only : pointer_meta

  implicit none

! the maximum level stored in the restart file
!
  integer(kind=4), save :: fcor = 1

! array of pointer used during job restart
!
  type(pointer_meta), dimension(:), allocatable, save :: block_array

  contains
!
!===============================================================================
!
! write_data: wrapper subroutine for storing the data
!
! info: subroutine selects the writing subroutine from supported output formats
!       depending on the compilation time options, e.g. at this moment only the
!       HDF5 format is supported;
!
! arguments:
!   ftype - character variable specifying the type of output file (i.e. restart
!           file or only the primitive variables);
!   nfile - integer variable counting the number of file for a given snapshots;
!   nproc - integer variable specifying the current processor number; the output
!           is divided into seperate files, one file per processor;
!
!===============================================================================
!
  subroutine write_data(ftype, nfile, nproc)

! declare variables
!
    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc
!
!-------------------------------------------------------------------------------
!
#ifdef HDF5
    call write_data_h5(ftype, nfile, nproc)
#endif /* HDF5 */

!-------------------------------------------------------------------------------
!
  end subroutine write_data
!
!===============================================================================
!
! restart_job: wrapper subroutine for the job restart from a data file
!
! info: subroutine selects the restoring subroutine from supported output
!       formats depending on the compilation time options, e.g. at this moment
!       only the HDF5 format is supported;
!
! arguments:
!   nfile - integer variable counting the number of file for a given snapshots;
!   nproc - integer variable specifying the current processor number; the output
!           is divided into seperate files, one file per processor;
!
!===============================================================================
!
  subroutine restart_job(nfile, nproc)

! declare variables
!
    implicit none

! input variables
!
    integer  , intent(in) :: nfile, nproc
!
!-------------------------------------------------------------------------------
!
#ifdef HDF5
    call read_data_h5(nfile, nproc)
#endif /* HDF5 */

!-------------------------------------------------------------------------------
!
  end subroutine restart_job
#ifdef HDF5
!
!===============================================================================
!
! write_data_h5: wrapper subroutine for the HDF5 format
!
! info: subroutine performs the initialization and finalization of the HDF5
!       interface, creates and closes the file, and also stores the parameters
!        andvariables in the chosen format (i.e. for the restart, visualization,
!       etc.)
!
! arguments: same as in write_data()
!
!===============================================================================
!
  subroutine write_data_h5(ftype, nfile, nproc)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t, H5F_ACC_TRUNC_F
    use hdf5 , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f

! declare variables
!
    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc

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

! check if the interface has been initialized successfuly
!
    if (err .ge. 0) then

! prepare the filename
!
      write (fl,'(a1,i6.6,"_",i5.5,a3)') ftype, nfile, nproc, '.h5'

! create the new HDF5 file
!
      call h5fcreate_f(fl, H5F_ACC_TRUNC_F, fid, err)

! check if the file has been created successfuly
!
      if (err .ge. 0) then

! write the global attributes
!
        call write_attributes_h5(fid)

! depending on the selected type of output file write the right groups
!
        select case(ftype)
        case('f')

! write the coordinates (data block bounds, refinement levels, etc.)
!
          call write_coordinates_h5(fid)

! write the variables stored in data blocks (leafs)
!
          call write_variables_full_h5(fid)

        case('p')

! write the coordinates (data block bounds, refinement levels, etc.)
!
          call write_coordinates_h5(fid)

! write the variables stored in data blocks (leafs)
!
          call write_variables_h5(fid)

        case('r')

! write all metablocks which represent the internal structure of domain
!
          call write_metablocks_h5(fid)

! write all datablocks which represent the all variables
!
          call write_datablocks_h5(fid)

        case default
        end select

! terminate access to the current file
!
        call h5fclose_f(fid, err)

! check if the file has been closed successfully
!
        if (err .gt. 0) then

! print error about the problem with closing the current file
!
          call print_error("io::write_data_h5"  &
                                            , "Cannot close file: " // trim(fl))

        end if

      else

! print error about the problem with creating the HDF5 file
!
        call print_error("io::write_data_h5"  &
                                           , "Cannot create file: " // trim(fl))

      end if

! close the FORTRAN interface
!
      call h5close_f(err)

! check if the interface has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the HDF5 Fortran interface
!
        call print_error("io::write_data_h5"  &
                                   , "Cannot close the HDF5 Fortran interface!")

      end if

    else

! print the error about the problem with initialization of the HDF5 Fortran
! interface
!
      call print_error("io::write_data_h5"  &
                              , "Cannot initialize the HDF5 Fortran interface!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_data_h5
!
!===============================================================================
!
! read_data_h5: wrapper subroutine for job restart from a HDF5 format
!
! info: subroutine reads and restores the meta and data block structures from
!       an HDF5 file
!
! arguments: same as in restart_job()
!
!===============================================================================
!
  subroutine read_data_h5(nfile, nproc)

! references to other modules
!
    use error, only : print_error
    use hdf5 , only : hid_t
    use hdf5 , only : H5F_ACC_RDONLY_F
    use hdf5 , only : h5open_f, h5close_f, h5fis_hdf5_f, h5fopen_f, h5fclose_f

! declare variables
!
    implicit none

! input variables
!
    integer  , intent(in) :: nfile, nproc

! local variables
!
    character(len=64) :: fl
    integer(hid_t)    :: fid
    integer           :: err
    logical           :: info
!
!-------------------------------------------------------------------------------
!
! initialize the FORTRAN interface
!
    call h5open_f(err)

! check if the interface has been initialized successfuly
!
    if (err .ge. 0) then

! prepare the filename
!
      write (fl,'("r",i6.6,"_",i5.5,a3)') nfile, nproc, '.h5'

! check if the HDF5 file exists
!
      inquire(file = fl, exist = info)

      if (info) then

! check if this is an HDF5 file
!
        call h5fis_hdf5_f(fl, info, err)

        if (err .ge. 0) then

          if (info) then

! opent the current HDF5 file
!
            call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! check if the file has been opened successfuly
!
            if (err .ge. 0) then

! read global attributes
!
              call read_attributes_h5(fid)

! read meta blocks
!
              call read_metablocks_h5(fid)

! read data blocks
!
              call read_datablocks_h5(fid)

! deallocate the array of block pointers
!
              if (allocated(block_array)) deallocate(block_array)

! terminate access to the current file
!
              call h5fclose_f(fid, err)

! check if the file has been closed successfully
!
              if (err .gt. 0) then

! print error about the problem with closing the current file
!
                call print_error("io::read_data_h5"                            &
                                            , "Cannot close file: " // trim(fl))

              end if

            else

! print error about the problem with opening the HDF5 file
!
              call print_error("io::read_data_h5"  &
                                             , "Cannot open file: " // trim(fl))

            end if

          else

! print error about the wrong file format
!
            call print_error("io::read_data_h5", "File " // trim(fl)           &
                                                     // " is not an HDF5 file!")
          end if

        else

! print error about the problem with checking the file format
!
          call print_error("io::read_data_h5", "Cannot check the file format!")

        end if

      else

! print error since files does not exist
!
        call print_error("io::read_data_h5", "File " // trim(fl)               &
                                                     // " does not exist!")
      end if

! close the FORTRAN interface
!
      call h5close_f(err)

! check if the interface has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the HDF5 Fortran interface
!
        call print_error("io::read_data_h5"  &
                                   , "Cannot close the HDF5 Fortran interface!")

      end if

    else

! print the error about the problem with initialization of the HDF5 Fortran
! interface
!
      call print_error("io::read_data_h5"  &
                              , "Cannot initialize the HDF5 Fortran interface!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_data_h5
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
    use blocks   , only : ndims, last_id, mblocks, dblocks, nleafs
    use config   , only : ncells, nghost
    use config   , only : xmin, xmax, ymin, ymax, zmin, zmax
    use config   , only : in, jn, kn, rdims, maxlev
    use error    , only : print_error
    use evolution, only : n, t, dt, dtn
    use hdf5     , only : hid_t
    use hdf5     , only : h5gcreate_f, h5gclose_f
    use mpitools , only : ncpus, ncpu
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
      call write_attribute_integer_h5(gid, 'ndims'  , ndims)
      call write_attribute_integer_h5(gid, 'last_id', last_id)
      call write_attribute_integer_h5(gid, 'mblocks', mblocks)
      call write_attribute_integer_h5(gid, 'dblocks', dblocks)
      call write_attribute_integer_h5(gid, 'nleafs' , nleafs)
      call write_attribute_integer_h5(gid, 'ncells' , ncells)
      call write_attribute_integer_h5(gid, 'nghost' , nghost)
      call write_attribute_integer_h5(gid, 'maxlev' , maxlev)
      call write_attribute_integer_h5(gid, 'ncpus'  , ncpus)
      call write_attribute_integer_h5(gid, 'ncpu'   , ncpu)
      call write_attribute_integer_h5(gid, 'nseeds' , nseeds)
      call write_attribute_integer_h5(gid, 'iter'   , n)

! store the real attributes
!
      call write_attribute_double_h5(gid, 'xmin', xmin)
      call write_attribute_double_h5(gid, 'xmax', xmax)
      call write_attribute_double_h5(gid, 'ymin', ymin)
      call write_attribute_double_h5(gid, 'ymax', ymax)
      call write_attribute_double_h5(gid, 'zmin', zmin)
      call write_attribute_double_h5(gid, 'zmax', zmax)
      call write_attribute_double_h5(gid, 'time', t   )
      call write_attribute_double_h5(gid, 'dt'  , dt  )
      call write_attribute_double_h5(gid, 'dtn' , dtn )

! store the vector attributes
!
      dm(:) = (/ in, jn, kn /)
      call write_attribute_vector_integer_h5(gid, 'dims' , 3, dm(:))
      call write_attribute_vector_integer_h5(gid, 'rdims', 3, rdims(:))

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
    use blocks   , only : append_metablock, append_datablock
    use blocks   , only : ndims, last_id, mblocks, dblocks, nleafs
    use config   , only : ncells, nghost
    use config   , only : in, jn, kn, rdims, maxlev
    use config   , only : xmin, xmax, ymin, ymax, zmin, zmax
    use error    , only : print_error, print_warning
    use evolution, only : n, t, dt, dtn
    use hdf5     , only : hid_t, hsize_t
    use hdf5     , only : h5gopen_f, h5gclose_f, h5aget_num_attrs_f            &
                        , h5aopen_idx_f, h5aclose_f, h5aget_name_f
    use mpitools , only : ncpus, ncpu
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
    integer           :: nattrs, lndims, llast_id, lmblocks, ldblocks          &
                       , lncells, lnghost, lnseeds, lmaxlev

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
              if (lmaxlev .gt. maxlev) then
                call print_warning("io::read_attributes_h5"                    &
                         , "The maximum refinement level has been decreased!")
              else
                fcor = 2**(maxlev - lmaxlev)
              end if
            case('last_id')
              call read_attribute_integer_h5(aid, aname, llast_id)
            case('mblocks')
              call read_attribute_integer_h5(aid, aname, lmblocks)
            case('dblocks')
              call read_attribute_integer_h5(aid, aname, ldblocks)
            case('nleafs')
              call read_attribute_integer_h5(aid, aname, nleafs)
            case('ncells')
              call read_attribute_integer_h5(aid, aname, lncells)

! check if the block dimensions are compatible
!
              if (lncells .ne. ncells) then
                call print_error("io::read_attributes_h5"                      &
                        , "File and program block dimensions are incompatible!")
              end if
            case('nghost')
              call read_attribute_integer_h5(aid, aname, lnghost)

! check if the ghost layers are compatible
!
              if (lnghost .ne. nghost) then
                call print_error("io::read_attributes_h5"                      &
                      , "File and program block ghost layers are incompatible!")
              end if
            case('iter')
              call read_attribute_integer_h5(aid, aname, n)
            case('time')
              call read_attribute_double_h5(aid, aname, t)
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

! check if the numbers of seeds are compatible
!
              if (lnseeds .ne. nseeds) then
                call print_error("io::read_attributes_h5"                      &
                , "The number of seeds from file and program are incompatible!")
              end if
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
                call set_seeds(seeds(:))

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
        if (lmblocks .ne. mblocks) then
          call print_error("io::read_attributes_h5"                            &
                                        , "Number of metablocks doesn't match!")
        end if

! allocate all datablocks
!
        do l = 1, ldblocks
          call append_datablock(pdata)
        end do

! check if the number of created datablocks is equal to the ldblocks
!
        if (ldblocks .ne. dblocks) then
          call print_error("io::read_attributes_h5"                            &
                                        , "Number of datablocks doesn't match!")
        end if

! allocate an array of pointers with the size llast_id
!
        allocate(block_array(llast_id))
        last_id = llast_id

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
    use blocks  , only : last_id, mblocks, nchild, ndims, nsides, nfaces
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
      am(1) = mblocks
      cm(1) = last_id
      dm(1) = mblocks
      dm(2) = nchild
      pm(1) = mblocks
      pm(2) = ndims
      qm(1) = mblocks
      qm(2) = ndims
      qm(3) = nsides
      qm(4) = nfaces

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
        cpu(l)   = pmeta%cpu
        lev(l)   = pmeta%level
        cfg(l)   = pmeta%config
        ref(l)   = pmeta%refine
        pos(l,:) = pmeta%pos(:)
        cor(l,:) = pmeta%coord(:)

        if (pmeta%leaf) lea(l) = 1

        xmn(l)   = pmeta%xmin
        xmx(l)   = pmeta%xmax
        ymn(l)   = pmeta%ymin
        ymx(l)   = pmeta%ymax
        zmn(l)   = pmeta%zmin
        zmx(l)   = pmeta%zmax

        do p = 1, nchild
          if (associated(pmeta%child(p)%ptr)) chl(l,p) = pmeta%child(p)%ptr%id
        end do

        do i = 1, ndims
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
    use blocks  , only : last_id, mblocks, nchild, ndims, nsides, nfaces
    use error   , only : print_error
    use hdf5    , only : hid_t, hsize_t
    use hdf5    , only : h5gopen_f, h5gclose_f

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
! open metablock group
!
    call h5gopen_f(fid, 'metablocks', gid, err)

! check if the group has been opened successfuly
!
    if (err .ge. 0) then

! prepate dimensions
!
      am(1) = mblocks
      dm(1) = mblocks
      dm(2) = nchild
      pm(1) = mblocks
      pm(2) = ndims
      qm(1) = mblocks
      qm(2) = ndims
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
      if (fcor .gt. 1) then
        cor(:,:) = cor(:,:) * fcor
      end if

! prepare the array of pointers to metablocks
!
      l = 1
      pmeta => list_meta
      do while(associated(pmeta))

        block_array(id(l))%ptr => pmeta

        pmeta%id       = id(l)
        pmeta%cpu      = cpu(l)
        pmeta%level    = lev(l)
        pmeta%config   = cfg(l)
        pmeta%refine   = ref(l)
        pmeta%xmin     = xmn(l)
        pmeta%xmax     = xmx(l)
        pmeta%ymin     = ymn(l)
        pmeta%ymax     = ymx(l)
        pmeta%zmin     = zmn(l)
        pmeta%zmax     = zmx(l)
        pmeta%pos(:)   = pos(l,:)
        pmeta%coord(:) = cor(l,:)

        if (lea(l) .eq. 1) then
          pmeta%leaf = .true.
        else
          pmeta%leaf = .false.
        end if

        l = l + 1
        pmeta => pmeta%next
      end do

! iterate over all metablocks and restore pointers
!
      l = 1
      pmeta => list_meta
      do while(associated(pmeta))

        pmeta%parent => block_array(par(l))%ptr

        do p = 1, nchild
          if (chl(l,p) .gt. 0) then
            pmeta%child(p)%ptr => block_array(chl(l,p))%ptr
          end if
        end do

        do i = 1, ndims
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
! write_datablocks_h5: subroutine writes datablocks in the HDF5 format connected
!                      to the provided identificator
!
! info: this subroutine stores only the datablocks
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine write_datablocks_h5(fid)

! references to other modules
!
    use blocks   , only : block_meta, block_data, list_data
    use blocks   , only : dblocks, ndims
    use config   , only : im, jm, km
    use error    , only : print_error
    use hdf5     , only : hid_t, hsize_t
    use hdf5     , only : h5gcreate_f, h5gclose_f
    use variables, only : nqt

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: l
    integer                        :: err
    integer(hsize_t), dimension(1) :: am
    integer(hsize_t), dimension(5) :: cm, dm
    integer(hsize_t), dimension(6) :: qm

! local allocatable arrays
!
    integer(kind=4), dimension(:)          , allocatable :: met
    real(kind=8)   , dimension(:,:,:,:,:)  , allocatable :: u

! local pointers
!
    type(block_meta), pointer :: pmeta
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create the group for datablocks
!
    call h5gcreate_f(fid, 'datablocks', gid, err)

! check if the group has been created successfuly
!
    if (err .ge. 0) then

! store data blocks only if there are some on the current processor
!
      if (dblocks .gt. 0) then

! prepate dimensions
!
        am(1) = dblocks
        cm(1) = dblocks
        dm(1) = dblocks
        dm(2) = nqt
        dm(3) = im
        dm(4) = jm
        dm(5) = km
        qm(1) = dblocks
        qm(2) = ndims
        qm(3) = nqt
        qm(4) = im
        qm(5) = jm
        qm(6) = km
        cm(2) = 3
        cm(3) = im
        cm(4) = jm
        cm(5) = km

! allocate arrays to store datablocks data
!
        allocate(met(am(1)))
        allocate(u  (dm(1),dm(2),dm(3),dm(4),dm(5)))

! iterate over all metablocks and fill in the arrays for storage
!
        l = 1
        pdata => list_data
        do while(associated(pdata))

          if (associated(pdata%meta)) met(l) = pdata%meta%id

          u(l,:,:,:,:)   = pdata%u(:,:,:,:)

          l = l + 1
          pdata => pdata%next
        end do

! store datablocks in the HDF5 file
!
        call write_vector_integer_h5(gid, 'meta', am(1), met)
        call write_array5_double_h5 (gid, 'u'   , dm(:), u)

! deallocate allocatable arrays
!
        if (allocated(met)) deallocate(met)
        if (allocated(u)  ) deallocate(u)

      end if ! dblocks > 0

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

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
! read_datablocks_h5: subroutine restored datablocks from the HDF5 restart file
!
! info: this subroutine restores only the datablocks
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine read_datablocks_h5(fid)

! references to other modules
!
    use blocks   , only : block_meta, block_data, list_data
    use blocks   , only : associate_blocks
    use blocks   , only : dblocks, ndims
    use config   , only : im, jm, km
    use error    , only : print_error
    use hdf5     , only : hid_t, hsize_t
    use hdf5     , only : h5gopen_f, h5gclose_f
    use variables, only : nqt

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! local variables
!
    integer(hid_t)                 :: gid
    integer(kind=4)                :: l
    integer                        :: err
    integer(hsize_t), dimension(1) :: am
    integer(hsize_t), dimension(5) :: dm

! local allocatable arrays
!
    integer(kind=4), dimension(:)        , allocatable :: m
    real(kind=8)   , dimension(:,:,:,:,:), allocatable :: u

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! open the datablock group
!
    call h5gopen_f(fid, 'datablocks', gid, err)

! check if the datablock group has been opened successfuly
!
    if (err .ge. 0) then

! restore all data blocks
!
      if (dblocks .gt. 0) then

! prepate dimensions
!
        am(1) = dblocks
        dm(1) = dblocks
        dm(2) = nqt
        dm(3) = im
        dm(4) = jm
        dm(5) = km

! allocate array to restore datablocks data
!
        allocate(m(am(1)))
        allocate(u(dm(1),dm(2),dm(3),dm(4),dm(5)))

! read datablocks from the HDF5 file
!
        call read_vector_integer_h5(gid, 'meta', am(:), m(:))
        call read_array5_double_h5 (gid, 'u'   , dm(:), u(:,:,:,:,:))

! iterate over all data blocks and fill their U arrays
!
        l = 1
        pdata => list_data
        do while(associated(pdata))

          call associate_blocks(block_array(m(l))%ptr, pdata)

          pdata%u(:,:,:,:) = u(l,:,:,:,:)

          l = l + 1
          pdata => pdata%next
        end do

! deallocate allocatable arrays
!
        if (allocated(m)) deallocate(m)
        if (allocated(u)) deallocate(u)

      end if ! dblocks > 0

! close the group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::read_datablocks_h5"                              &
                                              , "Cannot close datablock group!")

      end if

    else

! print error about the problem with opening the group
!
      call print_error("io::read_datablocks_h5", "Cannot open datablock group!")

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
    use blocks, only : dblocks, ndims, nsides, res
    use config, only : maxlev
    use error , only : print_error
    use hdf5  , only : hid_t, hsize_t
    use hdf5  , only : h5gcreate_f, h5gclose_f
    use mesh  , only : adx, ady, adz

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
    integer(hsize_t)  :: am(1), cm(2), dm(3)

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
      if (dblocks .gt. 0) then

! prepare dimensions
!
        am(1) = maxlev
        cm(1) = dblocks
        cm(2) = ndims
        dm(1) = dblocks
        dm(2) = ndims
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
          cor(l,:)   = pdata%meta%coord(:)

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
        call write_vector_integer_h5(gid, 'blkres', am(1), res)
        call write_vector_integer_h5(gid, 'levels', cm(1), lev)
        call write_vector_integer_h5(gid, 'refine', cm(1), ref)
        call write_array2_integer_h5(gid, 'coords', cm(:), cor)
        call write_array3_double_h5 (gid, 'bounds', dm(:), bnd)
        call write_vector_double_h5 (gid, 'dx'    , am(1), adx)
        call write_vector_double_h5 (gid, 'dy'    , am(1), ady)
        call write_vector_double_h5 (gid, 'dz'    , am(1), adz)

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
! write_variables_full_h5: subroutine writes each variable from datablocks in a
!                          separate array in the HDF5 format connected to the
!                          provided identificator
!
! info: this subroutine stores variables with ghost cells
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine write_variables_full_h5(fid)

! references to other modules
!
    use blocks       , only : block_data, list_data
    use blocks       , only : dblocks
    use config       , only : im, jm, km
    use error        , only : print_error
    use hdf5         , only : hid_t, hsize_t
    use hdf5         , only : h5gcreate_f, h5gclose_f
    use scheme       , only : cons2prim
    use variables    , only : nvr, nqt
    use variables    , only : idn, imx, imy, imz, ivx, ivy, ivz
#ifdef ADI
    use variables    , only : ien, ipr
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#ifdef GLM
    use variables    , only : iph
#endif /* GLM */
#endif /* MHD */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! HDF5 variables
!
    integer(hid_t)    :: gid
    integer(hsize_t)  :: dm(4)

! local variables
!
    integer           :: err
    integer(kind=4)   :: i, j, k, l

! local allocatable arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable :: u, q
    real(kind=8), dimension(:,:,:,:), allocatable :: dens
    real(kind=8), dimension(:,:,:,:), allocatable :: momx, momy, momz
    real(kind=8), dimension(:,:,:,:), allocatable :: velx, vely, velz
#ifdef ADI
    real(kind=8), dimension(:,:,:,:), allocatable :: ener, pres
#endif /* ADI */
#ifdef MHD
    real(kind=8), dimension(:,:,:,:), allocatable :: magx, magy, magz
#ifdef GLM
    real(kind=8), dimension(:,:,:,:), allocatable :: bpsi
#endif /* GLM */
#endif /* MHD */

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create a group to store global attributes
!
    call h5gcreate_f(fid, 'variables', gid, err)

! check if the group has been created successfuly
!
    if (err .ge. 0) then

! store variables only if there are some data blocks on the current processor
!
      if (dblocks .gt. 0) then

! prepare dimensions
!
        dm(1) = dblocks
        dm(2) = im
        dm(3) = jm
        dm(4) = km

! allocate arrays to store variables from all datablocks
!
        allocate(u(nvr,im,jm,km))
        allocate(q(nvr,im,jm,km))

        allocate(dens(dm(1),dm(2),dm(3),dm(4)))
        allocate(momx(dm(1),dm(2),dm(3),dm(4)))
        allocate(momy(dm(1),dm(2),dm(3),dm(4)))
        allocate(momz(dm(1),dm(2),dm(3),dm(4)))
        allocate(velx(dm(1),dm(2),dm(3),dm(4)))
        allocate(vely(dm(1),dm(2),dm(3),dm(4)))
        allocate(velz(dm(1),dm(2),dm(3),dm(4)))
#ifdef ADI
        allocate(ener(dm(1),dm(2),dm(3),dm(4)))
        allocate(pres(dm(1),dm(2),dm(3),dm(4)))
#endif /* ADI */
#ifdef MHD
        allocate(magx(dm(1),dm(2),dm(3),dm(4)))
        allocate(magy(dm(1),dm(2),dm(3),dm(4)))
        allocate(magz(dm(1),dm(2),dm(3),dm(4)))
#ifdef GLM
        allocate(bpsi(dm(1),dm(2),dm(3),dm(4)))
#endif /* GLM */
#endif /* MHD */

! iterate over all data blocks and fill in the arrays
!
        l = 1
        pdata => list_data
        do while(associated(pdata))

! copy conserved variables from the current block to the temporary array
!
          u(1:nqt,1:im,1:jm,1:km) = pdata%u(1:nqt,1:im,1:jm,1:km)

! obtain the primitive variables from the conserved ones
!
          do k = 1, km
            do j = 1, jm
              call cons2prim(im, u(:,:,j,k), q(:,:,j,k))
            end do
          end do

          dens(l,1:im,1:jm,1:km) = u(idn,1:im,1:jm,1:km)
          momx(l,1:im,1:jm,1:km) = u(imx,1:im,1:jm,1:km)
          momy(l,1:im,1:jm,1:km) = u(imy,1:im,1:jm,1:km)
          momz(l,1:im,1:jm,1:km) = u(imz,1:im,1:jm,1:km)
          velx(l,1:im,1:jm,1:km) = q(ivx,1:im,1:jm,1:km)
          vely(l,1:im,1:jm,1:km) = q(ivy,1:im,1:jm,1:km)
          velz(l,1:im,1:jm,1:km) = q(ivz,1:im,1:jm,1:km)
#ifdef ADI
          ener(l,1:im,1:jm,1:km) = u(ien,1:im,1:jm,1:km)
          pres(l,1:im,1:jm,1:km) = q(ipr,1:im,1:jm,1:km)
#endif /* ADI */
#ifdef MHD
          magx(l,1:im,1:jm,1:km) = u(ibx,1:im,1:jm,1:km)
          magy(l,1:im,1:jm,1:km) = u(iby,1:im,1:jm,1:km)
          magz(l,1:im,1:jm,1:km) = u(ibz,1:im,1:jm,1:km)
#ifdef GLM
          bpsi(l,1:im,1:jm,1:km) = q(iph,1:im,1:jm,1:km)
#endif /* GLM */
#endif /* MHD */

          l = l + 1
          pdata => pdata%next
        end do

! write the variables to the HDF5 file
!
        call write_array4_double_h5(gid, 'dens', dm, dens)
        call write_array4_double_h5(gid, 'momx', dm, momx)
        call write_array4_double_h5(gid, 'momy', dm, momy)
        call write_array4_double_h5(gid, 'momz', dm, momz)
        call write_array4_double_h5(gid, 'velx', dm, velx)
        call write_array4_double_h5(gid, 'vely', dm, vely)
        call write_array4_double_h5(gid, 'velz', dm, velz)
#ifdef ADI
        call write_array4_double_h5(gid, 'ener', dm, ener)
        call write_array4_double_h5(gid, 'pres', dm, pres)
#endif /* ADI */
#ifdef MHD
        call write_array4_double_h5(gid, 'magx', dm, magx)
        call write_array4_double_h5(gid, 'magy', dm, magy)
        call write_array4_double_h5(gid, 'magz', dm, magz)
#ifdef GLM
        call write_array4_double_h5(gid, 'bpsi', dm, bpsi)
#endif /* GLM */
#endif /* MHD */

! deallocate allocatable arrays
!
        if (allocated(dens)) deallocate(dens)
        if (allocated(momx)) deallocate(momx)
        if (allocated(momy)) deallocate(momy)
        if (allocated(momz)) deallocate(momz)
        if (allocated(velx)) deallocate(velx)
        if (allocated(vely)) deallocate(vely)
        if (allocated(velz)) deallocate(velz)
#ifdef ADI
        if (allocated(ener)) deallocate(ener)
        if (allocated(pres)) deallocate(pres)
#endif /* ADI */
#ifdef MHD
        if (allocated(magx)) deallocate(magx)
        if (allocated(magy)) deallocate(magy)
        if (allocated(magz)) deallocate(magz)
#ifdef GLM
        if (allocated(bpsi)) deallocate(bpsi)
#endif /* GLM */
#endif /* MHD */
        if (allocated(u))    deallocate(u)
        if (allocated(q))    deallocate(q)

      end if ! dblocks > 0

! close the attribute group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_variables_full_h5"  &
                                                    , "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_variables_full_h5"  &
                                                   , "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_variables_full_h5
!
!===============================================================================
!
! write_variables_h5: subroutine writes each variable from datablocks in a
!                     separate array in the HDF5 format connected to the
!                     provided identificator
!
! info: this subroutine stores variables
!
! arguments:
!   fid - the HDF5 file identificator
!
!===============================================================================
!
  subroutine write_variables_h5(fid)

! references to other modules
!
    use blocks       , only : block_data, list_data
    use blocks       , only : dblocks
    use config       , only : im, jm, km, in, jn, kn, ib, ie, jb, je, kb, ke
    use error        , only : print_error
    use hdf5         , only : hid_t, hsize_t
    use hdf5         , only : h5gcreate_f, h5gclose_f
    use scheme       , only : cons2prim
    use variables    , only : nvr, nqt
    use variables    , only : idn, ivx, ivy, ivz
#ifdef ADI
    use variables    , only : ipr
#endif /* ADI */
#ifdef MHD
    use variables    , only : ibx, iby, ibz
#endif /* MHD */

! declare variables
!
    implicit none

! input variables
!
    integer(hid_t), intent(in) :: fid

! HDF5 variables
!
    integer(hid_t)    :: gid
    integer(hsize_t)  :: dm(4)

! local variables
!
    integer           :: err
    integer(kind=4)   :: i, j, k, l

! local allocatable arrays
!
    real(kind=8), dimension(:,:,:,:), allocatable :: u, q
    real(kind=8), dimension(:,:,:,:), allocatable :: dens, velx, vely, velz
#ifdef ADI
    real(kind=8), dimension(:,:,:,:), allocatable :: pres
#endif /* ADI */
#ifdef MHD
    real(kind=8), dimension(:,:,:,:), allocatable :: magx, magy, magz
#endif /* MHD */

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
!
! create a group to store global attributes
!
    call h5gcreate_f(fid, 'variables', gid, err)

! check if the group has been created successfuly
!
    if (err .ge. 0) then

! store variables only if there are some data blocks on the current processor
!
      if (dblocks .gt. 0) then

! prepare dimensions
!
        dm(1) = dblocks
        dm(2) = in
        dm(3) = jn
        dm(4) = kn

! allocate arrays to store variables from all datablocks
!
        allocate(u(nvr       ,im,jm,km))
        allocate(q(nvr       ,im,jm,km))

        allocate(dens(dm(1),dm(2),dm(3),dm(4)))
        allocate(velx(dm(1),dm(2),dm(3),dm(4)))
        allocate(vely(dm(1),dm(2),dm(3),dm(4)))
        allocate(velz(dm(1),dm(2),dm(3),dm(4)))
#ifdef ADI
        allocate(pres(dm(1),dm(2),dm(3),dm(4)))
#endif /* ADI */
#ifdef MHD
        allocate(magx(dm(1),dm(2),dm(3),dm(4)))
        allocate(magy(dm(1),dm(2),dm(3),dm(4)))
        allocate(magz(dm(1),dm(2),dm(3),dm(4)))
#endif /* MHD */

! iterate over all data blocks and fill in the arrays
!
        l = 1
        pdata => list_data
        do while(associated(pdata))

! copy conserved variables from the current block to the temporary array
!
          u(1:nqt,1:im,1:jm,1:km) = pdata%u(1:nqt,1:im,1:jm,1:km)

! obtain the primitive variables from the conserved ones
!
          do k = 1, km
            do j = 1, jm
              call cons2prim(im, u(:,:,j,k), q(:,:,j,k))
            end do
          end do

          dens(l,1:in,1:jn,1:kn) = q(idn,ib:ie,jb:je,kb:ke)
          velx(l,1:in,1:jn,1:kn) = q(ivx,ib:ie,jb:je,kb:ke)
          vely(l,1:in,1:jn,1:kn) = q(ivy,ib:ie,jb:je,kb:ke)
          velz(l,1:in,1:jn,1:kn) = q(ivz,ib:ie,jb:je,kb:ke)
#ifdef ADI
          pres(l,1:in,1:jn,1:kn) = q(ipr,ib:ie,jb:je,kb:ke)
#endif /* ADI */
#ifdef MHD
          magx(l,1:in,1:jn,1:kn) = q(ibx,ib:ie,jb:je,kb:ke)
          magy(l,1:in,1:jn,1:kn) = q(iby,ib:ie,jb:je,kb:ke)
          magz(l,1:in,1:jn,1:kn) = q(ibz,ib:ie,jb:je,kb:ke)
#endif /* MHD */

          l = l + 1
          pdata => pdata%next
        end do

! write the variables to the HDF5 file
!
        call write_array4_double_h5(gid, 'dens', dm, dens)
        call write_array4_double_h5(gid, 'velx', dm, velx)
        call write_array4_double_h5(gid, 'vely', dm, vely)
        call write_array4_double_h5(gid, 'velz', dm, velz)
#ifdef ADI
        call write_array4_double_h5(gid, 'pres', dm, pres)
#endif /* ADI */
#ifdef MHD
        call write_array4_double_h5(gid, 'magx', dm, magx)
        call write_array4_double_h5(gid, 'magy', dm, magy)
        call write_array4_double_h5(gid, 'magz', dm, magz)
#endif /* MHD */

! deallocate allocatable arrays
!
        if (allocated(dens)) deallocate(dens)
        if (allocated(velx)) deallocate(velx)
        if (allocated(vely)) deallocate(vely)
        if (allocated(velz)) deallocate(velz)
#ifdef ADI
        if (allocated(pres)) deallocate(pres)
#endif /* ADI */
#ifdef MHD
        if (allocated(magx)) deallocate(magx)
        if (allocated(magy)) deallocate(magy)
        if (allocated(magz)) deallocate(magz)
#endif /* MHD */
        if (allocated(u))    deallocate(u)
        if (allocated(q))    deallocate(q)

      end if ! dblocks > 0

! close the attribute group
!
      call h5gclose_f(gid, err)

! check if the group has been closed successfuly
!
      if (err .gt. 0) then

! print error about the problem with closing the group
!
        call print_error("io::write_variables_h5", "Cannot close the group!")

      end if

    else

! print error about the problem with creating the group
!
      call print_error("io::write_variables_h5", "Cannot create the group!")

    end if

!-------------------------------------------------------------------------------
!
  end subroutine write_variables_h5
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
