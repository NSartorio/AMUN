!!******************************************************************************
!!
!! module: io - handling data input/output
!!
!! Copyright (C) 2008-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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

    use blocks  , only : block_data, list_data, ndims, nsides, nvr, nblocks    &
                       , nleafs, dblocks                                       &
                       , nqt, idn, ivx, ivy, ivz
#ifdef ADI
    use blocks  , only : ipr, ien
#endif /* ADI */
#ifdef MHD
    use blocks  , only : ibx, iby, ibz, icx, icy, icz
#endif /* MHD */
    use config  , only : nghost, in, jn, kn, im, jm, km, maxlev                &
                       , ib, ie, jb, je, kb, ke                                &
                       , xmin, xmax, ymin, ymax, zmin, zmax, rdims
    use error   , only : print_error
    use hdf5    , only : h5open_f, h5close_f, h5fcreate_f, h5fclose_f          &
                       , h5gcreate_f, h5gclose_f, h5acreate_f, h5aclose_f      &
                       , h5awrite_f, h5screate_simple_f, h5sclose_f            &
                       , h5dcreate_f, h5dwrite_f, h5dclose_f                   &
                       , h5pcreate_f, h5pset_chunk_f, h5pset_deflate_f         &
                       , h5pclose_f, hid_t, hsize_t, H5F_ACC_TRUNC_F           &
                       , H5T_NATIVE_CHARACTER, H5T_NATIVE_INTEGER              &
                       , H5T_NATIVE_DOUBLE, H5P_DATASET_CREATE_F
#if defined MHD && defined FLUXCT
    use interpolation, only : magtocen, divergence
#endif /* MHD & FLUXCT */
    use mesh    , only : ax, ay, az, adx, ady, adz
    use mpitools, only : ncpus, ncpu
    use scheme  , only : cons2prim

    implicit none

! input variables
!
    character, intent(in) :: ftype
    integer  , intent(in) :: nfile, nproc

! HDF5 variables
!
    integer(hid_t)    :: fid, gid, sid, aid, did, pid
    integer(hsize_t)  :: am(1), pm(3), qm(4)

! local variables
!
    character(len=64) :: fl
    integer           :: err
    integer(kind=4)   :: i, j, k, l
    integer           :: it, jt, kt

! local arrays
!
    integer(kind=4), dimension(3) :: dm

! local allocatable arrays
!
    integer(kind=4), dimension(:)      , allocatable :: indices
    integer(kind=4), dimension(:)      , allocatable :: levels
    real(kind=8)   , dimension(:,:,:)  , allocatable :: bounds
    real(kind=8)   , dimension(:,:,:,:), allocatable :: dens, pres, ener
    real(kind=8)   , dimension(:,:,:,:), allocatable :: velx, vely, velz
#ifdef MHD
    real(kind=8)   , dimension(:,:,:,:), allocatable :: magx, magy, magz
#ifdef FLUXCT
    real           , dimension(  :,:,:), allocatable :: db
    real(kind=8)   , dimension(:,:,:,:), allocatable :: divb
#endif /* FLUXCT */
#endif /* MHD */
    real           , dimension(:,:,:,:), allocatable :: u, q

! local pointers
!
    type(block_data), pointer :: pdata
!
!-------------------------------------------------------------------------------
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

! prepare dimensions
!
        dm(1) = in
        dm(2) = jn
        dm(3) = kn

        qm(1) = 1
        qm(2) = in
        qm(3) = jn
        qm(4) = kn

! allocate arrays for storing indices, levels, and variables
!
        allocate(indices(dblocks))
        allocate(levels (dblocks))
        allocate(bounds (dblocks,ndims,nsides))
        allocate(dens   (dblocks,in,jn,kn))
        allocate(velx   (dblocks,in,jn,kn))
        allocate(vely   (dblocks,in,jn,kn))
        allocate(velz   (dblocks,in,jn,kn))
#ifdef MHD
        allocate(magx   (dblocks,in,jn,kn))
        allocate(magy   (dblocks,in,jn,kn))
        allocate(magz   (dblocks,in,jn,kn))
#ifdef FLUXCT
        allocate(divb   (dblocks,in,jn,kn))
        allocate(db     (        im,jm,km))
#endif /* FLUXCT */
#endif /* MHD */
        allocate(pres   (dblocks,in,jn,kn))
        allocate(ener   (dblocks,in,jn,kn))
        allocate(u      (nvr    ,im,jm,km))
        allocate(q      (nvr    ,im,jm,km))

! iterate over all data blocks and fill the stored arrays
!
        l = 1
        pdata => list_data
        do while(associated(pdata))
          indices(l) = pdata%meta%id
          levels (l) = pdata%meta%level
          bounds (l,1,1) = pdata%meta%xmin
          bounds (l,1,2) = pdata%meta%xmax
          bounds (l,2,1) = pdata%meta%ymin
          bounds (l,2,2) = pdata%meta%ymax
#if NDIMS == 3
          bounds (l,3,1) = pdata%meta%zmin
          bounds (l,3,2) = pdata%meta%zmax
#endif /* NDIMS == 3 */
          u(1:nqt,1:im,1:jm,1:km) = pdata%u(1:nqt,1:im,1:jm,1:km)
#if defined MHD && defined FLUXCT
! interpolate cell centered components of the magnetic field
!
          call magtocen(u)
#endif /* MHD & FLUXCT */
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
          ener(l,1:in,1:jn,1:kn) = u(ien,ib:ie,jb:je,kb:ke)
#endif /* ADI */
#ifdef MHD
          magx(l,1:in,1:jn,1:kn) = q(icx,ib:ie,jb:je,kb:ke)
          magy(l,1:in,1:jn,1:kn) = q(icy,ib:ie,jb:je,kb:ke)
          magz(l,1:in,1:jn,1:kn) = q(icz,ib:ie,jb:je,kb:ke)
#ifdef FLUXCT
          call divergence(u(ibx:ibz,:,:,:), db)
          divb(l,1:in,1:jn,1:kn) = db(   ib:ie,jb:je,kb:ke)
#endif /* FLUXCT */
#endif /* MHD */
          l = l + 1
          pdata => pdata%next
        end do

! deallocate unused arrays
!
        if (allocated(u)) deallocate(u)
        if (allocated(q)) deallocate(q)

! prepare properties / compression
!
        call h5pcreate_f(H5P_DATASET_CREATE_F, pid, err)
        call h5pset_chunk_f(pid, 4, qm(:), err)
        call h5pset_deflate_f(pid, 9, err)

! create a group for the global attributes
!
        call h5gcreate_f(fid, 'attributes', gid, err)

        am(1) = 1
        call h5screate_simple_f(1, am, sid, err)

        call h5acreate_f(gid, 'dblocks', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, dblocks, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'nleafs', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, nleafs, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'nblocks', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, nblocks, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'nghost', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, nghost, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'maxlev', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, maxlev, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ncpus', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ncpus, am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'ncpu', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, ncpu, am, err)
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

        call h5sclose_f(sid, err)

        am(1) = 3
        call h5screate_simple_f(1, am, sid, err)

        call h5acreate_f(gid, 'dims', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, dm(:), am, err)
        call h5aclose_f(aid, err)

        call h5acreate_f(gid, 'rdims', H5T_NATIVE_INTEGER, sid, aid, err)
        call h5awrite_f(aid, H5T_NATIVE_INTEGER, rdims(:), am, err)
        call h5aclose_f(aid, err)

        call h5sclose_f(sid, err)

        am(1) = dblocks
        call h5screate_simple_f(1, am, sid, err)

        call h5dcreate_f(gid, 'indices', H5T_NATIVE_INTEGER, sid, did, err)
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, indices(:), am, err)
        call h5dclose_f(did, err)

        call h5dcreate_f(gid, 'levels', H5T_NATIVE_INTEGER, sid, did, err)
        call h5dwrite_f(did, H5T_NATIVE_INTEGER, levels(:), am, err)
        call h5dclose_f(did, err)

        call h5sclose_f(sid, err)

        call h5gclose_f(gid, err)

! create a group for the coordinates
!
        call h5gcreate_f(fid, 'coordinates', gid, err)

        if (ftype .eq. 'p') then
          pm(1) = dblocks
          pm(2) = ndims
          pm(3) = nsides
          call h5screate_simple_f(3, pm, sid, err)

          call h5dcreate_f(gid, 'bounds', H5T_NATIVE_DOUBLE, sid, did, err)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, bounds(:,:,:), pm(:), err, sid)
          call h5dclose_f(did, err)

          call h5sclose_f(sid, err)
        end if

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

! if file type is 'p' write data in a new group 'variables'
!
        if (ftype .eq. 'p') then

! create group for storing variables
!
          call h5gcreate_f(fid, 'variables', gid, err)

          qm(1) = dblocks
          qm(2) = in
          qm(3) = jn
          qm(4) = kn
          call h5screate_simple_f(4, qm, sid, err)

          call h5dcreate_f(gid, 'dens', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, dens(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'velx', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, velx(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'vely', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, vely(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'velz', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, velz(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
#ifdef ADI
          call h5dcreate_f(gid, 'pres', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, pres(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'ener', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, ener(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
#endif /* ADI */
#ifdef MHD
          call h5dcreate_f(gid, 'magx', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, magx(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'magy', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, magy(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
          call h5dcreate_f(gid, 'magz', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, magz(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
#ifdef FLUXCT
          call h5dcreate_f(gid, 'divb', H5T_NATIVE_DOUBLE, sid, did, err, pid)
          call h5dwrite_f(did, H5T_NATIVE_DOUBLE, divb(:,:,:,:), qm(:), err, sid)
          call h5dclose_f(did, err)
#endif /* FLUXCT */
#endif /* MHD */

          call h5sclose_f(sid, err)

! close group 'variables'
!
          call h5gclose_f(gid, err)

        end if

! deallocate all arrays
!
        if (allocated(indices)) deallocate(indices)
        if (allocated(levels )) deallocate(levels)
        if (allocated(bounds )) deallocate(bounds)
        if (allocated(dens   )) deallocate(dens)
        if (allocated(velx   )) deallocate(velx)
        if (allocated(vely   )) deallocate(vely)
        if (allocated(velz   )) deallocate(velz)
        if (allocated(pres   )) deallocate(pres)
#ifdef MHD
        if (allocated(magx   )) deallocate(magx)
        if (allocated(magy   )) deallocate(magy)
        if (allocated(magz   )) deallocate(magz)
#ifdef FLUXCT
        if (allocated(divb   )) deallocate(divb)
        if (allocated(db     )) deallocate(db)
#endif /* FLUXCT */
#endif /* MHD */

! release properties
!
        call h5pclose_f(pid, err)

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

!-------------------------------------------------------------------------------
!
  end subroutine write_data
#endif /* HDF5 */
!===============================================================================
!
end module
