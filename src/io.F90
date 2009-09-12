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
!===============================================================================
!
! write_data: subroutine writes all data to file for a given time step
!
!===============================================================================
!
  subroutine write_data(ftype, nfile, nproc)

    use blocks  , only : block_data, list_data, nv => nvars
    use config  , only : ncells, nghost, ngrids, igrids, jgrids, kgrids        &
                       , im, jm, km, maxlev, xmin, xmax, ymin, ymax, zmin, zmax
    use error   , only : print_error
    use hdf5    , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f          &
                       , h5gcreate_f, h5gclose_f, h5acreate_f, h5aclose_f      &
                       , h5awrite_f, h5screate_simple_f, h5sclose_f            &
                       , h5dcreate_f, h5dwrite_f, h5dclose_f                   &
                       , hid_t, hsize_t, H5F_ACC_TRUNC_F                       &
                       , H5T_NATIVE_CHARACTER, H5T_NATIVE_INTEGER              &
                       , H5T_NATIVE_DOUBLE
    use mesh    , only : ax, ay, az, adx, ady, adz
    use mpitools, only : ncpus, ncpu
    use scheme  , only : cons2prim
    use problem , only : check_ref

    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc

! HDF5 variables
!
    integer(hid_t)    :: fid, gid, sid, aid, did, bid
    integer(hsize_t)  :: am(1), cm(2), dm(3), pm(3)

! local variables
!
    character(len=64) :: fl, gnm
    integer           :: err, i, j, k, r

! pointers
!
    type(block_data), pointer :: pblock

! local arrays
!
    real, dimension(im,jm,km)    :: tmp, c
    real, dimension(nv,im,jm,km) :: u, v
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

! prepare dimensions
!
      dm(1) = im
      dm(2) = jm
      dm(3) = km

      pm(1) = 2
      pm(2) = 2
      pm(3) = 2

! create file
!
      call h5fcreate_f(fl, H5F_ACC_TRUNC_F, fid, err)

      if (err .ge. 0) then

! create a group for the global attributes
!
        call h5gcreate_f(fid, 'attributes', gid, err)
        am(1) = 1
        call h5screate_simple_f(1, am, sid, err)

        call h5acreate_f(gid, 'ncells', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ncells, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ngrids', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ngrids, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'igrids', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, igrids, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'jgrids', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, jgrids, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'kgrids', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, kgrids, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'nghost', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, nghost, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'maxlev', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, maxlev, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'xmin', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(xmin,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'xmax', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(xmax,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ymin', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(ymin,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ymax', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(ymax,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'zmin', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(zmin,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'zmax', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(zmax,8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ncpus', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ncpus, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ncpu', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ncpu, am, err)
        call h5aclose_f(aid, err)

        call h5sclose_f(sid, err)
        call h5gclose_f(gid, err)

! create a group for the coordinates
!
        call h5gcreate_f(fid, 'coordinates', gid, err)
        cm(1) = maxlev
        cm(2) = ngrids
        call h5screate_simple_f(1, cm, sid, err)

        call h5acreate_f(gid, 'x', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(ax(:,:),8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'y', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(ay(:,:),8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'z', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(az(:,:),8), am, err)
        call h5aclose_f(aid, err)

        call h5sclose_f(sid, err)

        am(1) = maxlev
        call h5screate_simple_f(1, am, sid, err)

        call h5acreate_f(gid, 'dx', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(adx(:),8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'dy', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(ady(:),8), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'dz', H5T_NATIVE_DOUBLE, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(adz(:),8), am, err)
        call h5aclose_f(aid, err)

        call h5sclose_f(sid, err)

        call h5gclose_f(gid, err)

! create a group for the block storage
!
        call h5gcreate_f(fid, 'blocks', gid, err)

! TODO: iterate over all blocks and write complete structure of each of them
!
        pblock => list_data

        do while(associated(pblock))

          print *, pblock%meta%id, pblock%meta%leaf, pblock%meta%neigh(1,2,1)%ptr%id

          if (pblock%meta%leaf) then

            write(gnm,"('blk',i8.8)") pblock%meta%id

            call h5gcreate_f(gid, gnm, bid, err)

! TODO: create block attributes (xmin, xmax, ..., level, parent, all entries of structure)
!
            am(1) = 1
            call h5screate_simple_f(1, am, sid, err)

            call h5acreate_f(bid, 'id', H5T_NATIVE_INTEGER, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, pblock%meta%id, am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'refine', H5T_NATIVE_INTEGER, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, pblock%meta%refine, am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'config', H5T_NATIVE_INTEGER, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, pblock%meta%config, am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'level', H5T_NATIVE_INTEGER, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_INTEGER, pblock%meta%level, am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'xmin', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%xmin,8), am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'xmax', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%xmax,8), am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'ymin', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%ymin,8), am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'ymax', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%ymax,8), am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'zmin', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%zmin,8), am, err)
            call h5aclose_f(aid, err)

            call h5acreate_f(bid, 'zmax', H5T_NATIVE_DOUBLE, sid, aid, err)
            call h5awrite_f(aid, H5T_NATIVE_DOUBLE, real(pblock%zmax,8), am, err)
            call h5aclose_f(aid, err)

            call h5sclose_f(sid, err)

! prepare field variables for writing
!
            do k = 1, km
              do j = 1, jm
                u(:,:,j,k) = pblock%u(:,:,j,k)
                call cons2prim(nv,im,u(:,:,j,k),v(:,:,j,k))
                c(:,j,k) = pblock%c(:,j,k)
              end do
            end do

! write field variables
!
            call h5screate_simple_f(3, dm(1:3), sid, err)

            call h5dcreate_f(bid, 'dens', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, v(1,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

            call h5dcreate_f(bid, 'velx', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, v(2,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

            call h5dcreate_f(bid, 'vely', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, v(3,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

            call h5dcreate_f(bid, 'velz', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, v(4,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

#ifdef ADI
            call h5dcreate_f(bid, 'pres', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, v(5,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

            call h5dcreate_f(bid, 'ener', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, u(5,:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)
#endif /* ADI */

            call h5dcreate_f(bid, 'crit', H5T_NATIVE_DOUBLE, sid, did, err)
            call h5dwrite_f(did, H5T_NATIVE_DOUBLE, c(:,:,:), dm(1:3), err, sid)
            call h5dclose_f(did, err)

            call h5sclose_f(sid, err)

            call h5gclose_f(bid, err)
          endif
          pblock => pblock%next
        end do

! close group 'blocks'
!
        call h5gclose_f(gid, err)

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
