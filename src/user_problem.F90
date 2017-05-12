!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2017 Grzegorz Kowal <grzegorz@amuncode.org>
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
!!*****************************************************************************
!!
!! module: USER_PROBLEM
!!
!!  This module provides subroutines to setup custom problem.
!!
!!*****************************************************************************
!
module user_problem

#ifdef PROFILE
! include external procedures
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer, save :: imi, imp, ims, imu, img, imb
#endif /* PROFILE */

! default parameter values
!
  real(kind=8), save :: rjet   = 1.00d+00
  real(kind=8), save :: ljet   = 0.00d+00
  real(kind=8), save :: dn_jet = 1.00d-01
  real(kind=8), save :: pr_jet = 1.00d-02
  real(kind=8), save :: lf_jet = 1.00d+01
  real(kind=8), save :: bt_jet = 1.00d+99
  real(kind=8), save :: an_jet = 0.00d+00
  real(kind=8), save :: tm_jet = 1.00d+99
  real(kind=8), save :: dn_amb = 1.00d+01
  real(kind=8), save :: pr_amb = 1.00d-02
  real(kind=8), save :: bt_amb = 1.00d+99
  real(kind=8), save :: tm_on  = 1.00d+99
  real(kind=8), save :: tm_off = 0.00d+00

! module parameters
!
  logical     , save :: state  = .true.
  real(kind=8), save :: rjet2  = 1.00d+00
  real(kind=8), save :: vjet   = 0.00d+00
  real(kind=8), save :: bjet   = 0.00d+00
  real(kind=8), save :: ajet   = 0.00d+00
  real(kind=8), save :: bamb   = 0.00d+00
  real(kind=8), save :: sina   = 0.00d+00
  real(kind=8), save :: cosa   = 1.00d+00
  real(kind=8), save :: tnext  = 1.00d+99

! flag indicating if the gravitational source term is enabled
!
  logical     , save :: gravity_enabled_user = .false.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_USER_PROBLEM:
! ----------------------------------
!
!   Subroutine initializes user problem. It could read problem parameters which
!   are used in all subroutines defining this specific problem.
!
!   Arguments:
!
!     verbose - a logical flag turning the information printing;
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine initialize_user_problem(verbose, iret)

! include external procedures and variables
!
    use constants , only : d2r
    use equations , only : ipr, ibx
    use parameters, only : get_parameter_string, get_parameter_real
    use random    , only : randomu

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=64) :: problem_name   = "none"
    character(len=64) :: sfmts, sfmtf, sfmti
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('user_problem:: initialize'   , imi)
    call set_timer('user_problem:: problem setup', imp)
    call set_timer('user_problem:: shape'        , ims)
    call set_timer('user_problem:: sources'      , imu)
    call set_timer('user_problem:: gravity'      , img)
    call set_timer('user_problem:: boundaries'   , imb)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! get the problem name
!
    call get_parameter_string("problem", problem_name)

! get problem parameters
!
    call get_parameter_real("rjet"   , rjet  )
    call get_parameter_real("ljet"   , ljet  )
    call get_parameter_real("dn_jet" , dn_jet)
    call get_parameter_real("lf_jet" , lf_jet)
    call get_parameter_real("pr_jet" , pr_jet)
    call get_parameter_real("bt_jet" , bt_jet)
    call get_parameter_real("an_jet" , an_jet)
    call get_parameter_real("tm_jet" , tm_jet)
    call get_parameter_real("dn_amb" , dn_amb)
    call get_parameter_real("pr_amb" , pr_amb)
    call get_parameter_real("bt_amb" , bt_amb)
    call get_parameter_real("tactive", tm_on )
    call get_parameter_real("tquiet" , tm_off)

! calculate Rjet²
!
    rjet2 = rjet * rjet

! calculate plasma parameters
!
    vjet  = sqrt(1.0d+00 - 1.0d+00 / lf_jet**2)
    if (bt_jet > 0.0d+00 .and. bt_jet < 1.0d+30) then
      bjet  = sqrt(2.0d+00 * pr_jet / bt_jet)
    else
      bjet  = 0.0d+00
    end if
    if (bt_amb > 0.0d+00 .and. bt_amb < 1.0d+30) then
      bamb  = sqrt(2.0d+00 * pr_amb / bt_amb)
    else
      bamb  = 0.0d+00
    end if

! convert jet precession angle to radians
!
    ajet  = d2r * an_jet
    sina  = sin(ajet)
    cosa  = cos(ajet)

! update the next jet engine switch moment
!
    state = .true.
    if (tm_on > 0.0d+00 .and. tm_off > 0.0d+00) tnext = tm_on  * randomu()

! print information about the user problem such as problem name, its
! parameters, etc.
!
    if (verbose) then

      sfmts = "(4x,a14,9x,'=',2x,a)"
      write (*,sfmts) "problem name  ", trim(problem_name)
      sfmts = "(6x,a24)"
      write (*,sfmts)          "jet parameters:         "
      sfmtf = "(8x,a10,9x,'=',1pe12.4,2x,a)"
      write (*,sfmtf) "rjet      ", rjet  , "(jet radius)"
      write (*,sfmtf) "ljet      ", ljet  , "(jet length)"
      write (*,sfmtf) "dn_jet    ", dn_jet, "(jet density)"
      write (*,sfmtf) "lf_jet    ", lf_jet, "(jet Lorentz factor)"
      if (ipr > 0) &
        write (*,sfmtf) "pr_jet      ", pr_jet, "(jet pressure)"
      if (ibx > 0) then
        if (bjet == 0.0d+00) then
          sfmti = "(8x,a10,9x,'=',2x,'infinity',4x,a)"
          write (*,sfmti) "bt_jet      ", "(jet plasma-β)"
        else
          write (*,sfmtf) "bt_jet      ", bt_jet, "(jet plasma-β)"
        end if
        write (*,sfmtf) "Bjet      ", bjet  , "(jet magnetic field amplitude)"
      end if
      write (*,sfmts)          "precession parameters:  "
      write (*,sfmtf) "an_jet    ", an_jet, "(inclination angle in deg.)"
      write (*,sfmtf) "tm_jet    ", tm_jet, "(precession period)"
      write (*,sfmts)          "ambient parameters:     "
      write (*,sfmtf) "dn_amb    ", dn_amb, "(ambient density)"
      if (ipr > 0) &
        write (*,sfmtf) "pr_amb      ", pr_amb, "(ambient pressure)"
      if (ibx > 0) then
        if (bamb == 0.0d+00) then
          sfmti = "(8x,a10,9x,'=',2x,'infinity',4x,a)"
          write (*,sfmti) "bt_amb      ", "(ambient plasma-β)"
        else
          write (*,sfmtf) "bt_amb      ", bt_amb, "(ambient plasma-β)"
        end if
        write (*,sfmtf) "Bamb      ", bamb  , "(ambient magnetic field amplitude)"
      end if
      write (*,sfmts)          "engine parameters:      "
      write (*,sfmtf) "tactive   ", tm_on , "(maximum active period length)"
      write (*,sfmtf) "tquiet    ", tm_off, "(maximum quiet period length)"

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_user_problem
!
!===============================================================================
!
! subroutine FINALIZE_USER_PROBLEM:
! --------------------------------
!
!   Subroutine releases memory used by the module.
!
!   Arguments:
!
!     iret    - an integer flag for error return value;
!
!===============================================================================
!
  subroutine finalize_user_problem(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_user_problem
!
!===============================================================================
!
! subroutine SETUP_PROBLEM_USER:
! -----------------------------
!
!   Subroutine sets the initial conditions for the user specific problem.
!
!   Arguments:
!
!     pdata - pointer to the datablock structure of the currently initialized
!             block;
!
!===============================================================================
!
  subroutine setup_problem_user(pdata)

! include external procedures and variables
!
    use blocks     , only : block_data
    use coordinates, only : im, jm, km
    use coordinates, only : ax, ay, az
    use equations  , only : prim2cons
    use equations  , only : nv
    use equations  , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    type(block_data), pointer, intent(inout) :: pdata

! local variables
!
    integer      :: i, j, k
    real(kind=8) :: dx, dy, dz, rm, rr

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(nv)    :: qj
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the problem setup
!
    call start_timer(imp)
#endif /* PROFILE */

! set the conditions inside the jet radius
!
    qj(idn) = dn_jet
    if (ipr > 0) qj(ipr) = pr_jet
    qj(ivx) = vjet * cosa
    qj(ivy) = vjet * sina
    qj(ivz) = 0.0d+00
    if (ibx > 0) then
      qj(ibx) = 0.0d+00
      qj(iby) = 0.0d+00
      qj(ibz) = bjet
      qj(ibp) = 0.0d+00
    end if ! ibx > 0

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    dx = x(2) - x(1)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
    dy = y(2) - y(1)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
    dz = z(2) - z(1)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
    dz      = 0.0d+00
#endif /* NDIMS == 3 */
    rm = dy * dy + dz * dz

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do j = 1, jm

! calculate radius
!
        rr = y(j) * y(j) + z(k) * z(k)

! set the ambient density, pressure, and velocity
!
        q(idn,1:im) = dn_amb
        if (ipr > 0) q(ipr,1:im) = pr_amb
        q(ivx,1:im) = 0.0d+00
        q(ivy,1:im) = 0.0d+00
        q(ivz,1:im) = 0.0d+00

! if magnetic field is present, set it to be uniform with the desired strength
! and orientation
!
        if (ibx > 0) then
          q(ibx,1:im) = 0.0d+00
          q(iby,1:im) = 0.0d+00
          q(ibz,1:im) = bamb
          q(ibp,1:im) = 0.0d+00
        end if ! ibx > 0

! set the jet injection
!
        if (rr <= max(rm, rjet2)) then
          do i = 1, im
            if (x(i) <= max(dx, ljet)) then
              q(1:nv,i) = qj(1:nv)
            end if
          end do ! i = 1, im
        end if ! R < Rjet

! convert the primitive variables to conservative ones
!
        call prim2cons(im, q(1:nv,1:im), u(1:nv,1:im))

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the problems setup
!
    call stop_timer(imp)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine setup_problem_user
!
!===============================================================================
!
! subroutine UPDATE_SHAPES_USER:
! -----------------------------
!
!   Subroutine defines the regions updated by user.
!
!   Arguments:
!
!     pdata - pointer to the data block structure of the currently initialized
!             block;
!     time  - time at the moment of update;
!     dt    - time step since the last update;
!
!===============================================================================
!
  subroutine update_shapes_user(pdata, time, dt)

! include external procedures and variables
!
    use blocks         , only : block_data
    use constants      , only : pi2
    use coordinates    , only : im, jm, km
    use coordinates    , only : ax, ay, az
    use equations      , only : prim2cons
    use equations      , only : nv
    use equations      , only : idn, ivx, ivy, ivz, ipr, ibx, iby, ibz, ibp
    use random         , only : randomu

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer, intent(inout) :: pdata
    real(kind=8)             , intent(in)    :: time, dt

! local variables
!
    integer       :: i, j, k
    real(kind=8)  :: dx, dy, dz, rm, rr
    real(kind=8)  :: tph, sint, cost

! local arrays
!
    real(kind=8), dimension(nv,im) :: q, u
    real(kind=8), dimension(nv)    :: qj, uj
    real(kind=8), dimension(im)    :: x
    real(kind=8), dimension(jm)    :: y
    real(kind=8), dimension(km)    :: z
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the shape update
!
    call start_timer(ims)
#endif /* PROFILE */

! calculate the precession angle in the perpendicular direction
!
    tph  = pi2 * time / tm_jet
    sint = sin(tph)
    cost = cos(tph)

! get the next engine switch moment
!
    if (time >= tnext) then
      if (state) then
        tnext = tnext + tm_off * randomu()
      else
        tnext = tnext + tm_on  * randomu()
      end if
      state = .not. state
    end if

! set the conditions inside the jet radius
!
    qj(idn) = dn_jet
    if (ipr > 0) qj(ipr) = pr_jet
    if (state) then
      qj(ivx) = vjet * cosa
      qj(ivy) = vjet * sina * cost
      qj(ivz) = vjet * sina * sint
    else
      qj(ivx) = 0.0d+00
      qj(ivy) = 0.0d+00
      qj(ivz) = 0.0d+00
    end if
    if (ibx > 0) then
      qj(ibx) = 0.0d+00
      qj(iby) = 0.0d+00
      qj(ibz) = bjet
      qj(ibp) = 0.0d+00
    end if ! ibx > 0
    call prim2cons(1, qj(1:nv), uj(1:nv))

! prepare block coordinates
!
    x(1:im) = pdata%meta%xmin + ax(pdata%meta%level,1:im)
    dx = x(2) - x(1)
    y(1:jm) = pdata%meta%ymin + ay(pdata%meta%level,1:jm)
    dy = y(2) - y(1)
#if NDIMS == 3
    z(1:km) = pdata%meta%zmin + az(pdata%meta%level,1:km)
    dz = z(2) - z(1)
#else /* NDIMS == 3 */
    z(1:km) = 0.0d+00
    dz      = 0.0d+00
#endif /* NDIMS == 3 */
    rm = dy * dy + dz * dz

! iterate over all positions in the YZ plane
!
    do k = 1, km
      do j = 1, jm

! copy the primitive variable vector
!
        q(1:nv,1:im) = pdata%q(1:nv,1:im,j,k)
        u(1:nv,1:im) = pdata%u(1:nv,1:im,j,k)

! calculate radius
!
        rr = y(j) * y(j) + z(k) * z(k)

        if (rr <= max(rm, rjet2)) then
          do i = 1, im
            if (x(i) <= max(dx, ljet)) then
              q(1:nv,i) = qj(1:nv)
              u(1:nv,i) = uj(1:nv)
            end if
          end do ! i = 1, im
        end if ! R < Rjet

! copy the primitive variables to the current block
!
        pdata%q(1:nv,1:im,j,k) = q(1:nv,1:im)

! copy the conserved variables to the current block
!
        pdata%u(1:nv,1:im,j,k) = u(1:nv,1:im)

      end do ! j = 1, jm
    end do ! k = 1, km

#ifdef PROFILE
! stop accounting time for the shape update
!
    call stop_timer(ims)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_shapes_user
!
!===============================================================================
!
! subroutine UPDATE_SOURCES_USER:
! ------------------------------
!
!   Subroutine adds the user defined source terms.
!
!   Arguments:
!
!     pdata - the pointer to a data block;
!     t, dt - the time and time increment;
!     du    - the array of variable increment;
!
!===============================================================================
!
  subroutine update_sources_user(pdata, t, dt, du)

! include external variables
!
    use blocks         , only : block_data
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    type(block_data), pointer           , intent(inout) :: pdata
    real(kind=8)                        , intent(in)    :: t, dt
    real(kind=8), dimension(nv,im,jm,km), intent(inout) :: du
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for source terms
!
    call start_timer(imu)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for source terms
!
    call stop_timer(imu)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine update_sources_user
!
!===============================================================================
!
! subroutine GRAVITATIONAL_ACCELERATION_USER:
! ------------------------------------------
!
!   Subroutine returns the user defined gravitational acceleration.
!
!   Arguments:
!
!     t, dt   - time and the time increment;
!     x, y, z - rectangular coordinates;
!     acc     - vector of the gravitational acceleration;
!
!===============================================================================
!
  subroutine gravitational_acceleration_user(t, dt, x, y, z, acc)

! include external procedures and variables
!
    use parameters , only : get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8)              , intent(in)  :: t, dt
    real(kind=8)              , intent(in)  :: x, y, z
    real(kind=8), dimension(3), intent(out) :: acc
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the gravitational acceleration calculation
!
    call start_timer(img)
#endif /* PROFILE */

! reset gravitational acceleration
!
    acc(:) = 0.0d+00

#ifdef PROFILE
! stop accounting time for the gravitational acceleration calculation
!
    call stop_timer(img)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine gravitational_acceleration_user
!
!===============================================================================
!
! subroutine BOUNDARY_USER_X:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the X direction.
!
!   Arguments:
!
!     ic      - the block side along the X direction for the ghost zone update;
!     jl, ju  - the cell index limits for the Y direction;
!     kl, ku  - the cell index limits for the Z direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_x(ic, jl, ju, kl, ku, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: ic
    integer                                     , intent(in)    :: jl, ju
    integer                                     , intent(in)    :: kl, ku
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_x
!
!===============================================================================
!
! subroutine BOUNDARY_USER_Y:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the Y direction.
!
!   Arguments:
!
!     jc      - the block side along the Y direction for the ghost zone update;
!     il, iu  - the cell index limits for the X direction;
!     kl, ku  - the cell index limits for the Z direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_y(jc, il, iu, kl, ku, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: jc
    integer                                     , intent(in)    :: il, iu
    integer                                     , intent(in)    :: kl, ku
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_y
!
!===============================================================================
!
! subroutine BOUNDARY_USER_Z:
! --------------------------
!
!   Subroutine updates ghost zones within the specific region along
!   the Z direction.
!
!   Arguments:
!
!     kc      - the block side along the Z direction for the ghost zone update;
!     il, iu  - the cell index limits for the X direction;
!     jl, ju  - the cell index limits for the Y direction;
!     t, dt   - time and time increment;
!     x, y, z - the block coordinates;
!     qn      - the array of variables to update;
!
!===============================================================================
!
  subroutine boundary_user_z(kc, il, iu, jl, ju, t, dt, x, y, z, qn)

! import external procedures and variables
!
    use coordinates    , only : im, jm, km
    use equations      , only : nv

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                                     , intent(in)    :: kc
    integer                                     , intent(in)    :: il, iu
    integer                                     , intent(in)    :: jl, ju
    real(kind=8)                                , intent(in)    :: t, dt
    real(kind=8), dimension(1:im)               , intent(in)    :: x
    real(kind=8), dimension(1:jm)               , intent(in)    :: y
    real(kind=8), dimension(1:km)               , intent(in)    :: z
    real(kind=8), dimension(1:nv,1:im,1:jm,1:km), intent(inout) :: qn
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for the boundary update
!
    call start_timer(imb)
#endif /* PROFILE */

#ifdef PROFILE
! stop accounting time for the boundary update
!
    call stop_timer(imb)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine boundary_user_z

!===============================================================================
!
end module user_problem
