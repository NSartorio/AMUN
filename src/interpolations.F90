!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2018 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: INTERPOLATIONS
!!
!!  This module provides subroutines to interpolate variables and reconstruct
!!  the Riemann states.
!!
!!
!!******************************************************************************
!
module interpolations

#ifdef PROFILE
! import external subroutines
!
  use timers, only : set_timer, start_timer, stop_timer
#endif /* PROFILE */

! module variables are not implicit by default
!
  implicit none

#ifdef PROFILE
! timer indices
!
  integer            , save :: imi, imr, imf, imc
#endif /* PROFILE */

! pointers to the reconstruction and limiter procedures
!
  procedure(interfaces_tvd)    , pointer, save :: interfaces         => null()
  procedure(reconstruct)       , pointer, save :: reconstruct_states => null()
  procedure(limiter_zero)      , pointer, save :: limiter_tvd        => null()
  procedure(limiter_zero)      , pointer, save :: limiter_prol       => null()
  procedure(limiter_zero)      , pointer, save :: limiter_clip       => null()

! module parameters
!
  real(kind=8), save :: eps        = epsilon(1.0d+00)
  real(kind=8), save :: rad        = 0.5d+00

! monotonicity preserving reconstruction coefficients
!
  real(kind=8), save :: kappa      = 1.0d+00
  real(kind=8), save :: kbeta      = 1.0d+00

! number of ghost zones (required for compact schemes)
!
  integer     , save :: ng         = 2

! number of cells used in the Gaussian process reconstruction
!
  integer     , save :: ngp        = 5
  integer     , save :: mgp        = 1
  integer     , save :: dgp        = 9

! normal distribution width in the Gaussian process reconstruction
!
  real(kind=8), save :: sgp        = 9.0d+00

! PPM constant
!
  real(kind=8), save :: ppm_const  = 1.25d+00

! Gaussian process reconstruction coefficients vector
!
  real(kind=8), dimension(:)    , allocatable, save :: cgp
  real(kind=8), dimension(:,:,:), allocatable, save :: ugp

! flags for reconstruction corrections
!
  logical     , save :: mlp        = .false.
  logical     , save :: positivity = .false.
  logical     , save :: clip       = .false.

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_interpolations, finalize_interpolations
  public :: interfaces, reconstruct, limiter_prol
  public :: fix_positivity

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_INTERPOLATIONS:
! ------------------------------------
!
!   Subroutine initializes the interpolation module by reading the module
!   parameters.
!
!
!===============================================================================
!
  subroutine initialize_interpolations(verbose, iret)

! include external procedures
!
    use iso_fortran_env, only : error_unit
    use parameters     , only : get_parameter_string, get_parameter_integer    &
                              , get_parameter_real

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)    :: verbose
    integer, intent(inout) :: iret

! local variables
!
    character(len=255) :: sreconstruction = "tvd"
    character(len=255) :: tlimiter        = "mm"
    character(len=255) :: plimiter        = "mm"
    character(len=255) :: climiter        = "mm"
    character(len=255) :: mlp_limiting    = "off"
    character(len=255) :: positivity_fix  = "off"
    character(len=255) :: clip_extrema    = "off"
    character(len=255) :: name_rec        = ""
    character(len=255) :: name_tlim       = ""
    character(len=255) :: name_plim       = ""
    character(len=255) :: name_clim       = ""
    character(len= 16) :: stmp
    real(kind=8)       :: cfl             = 0.5d+00

! local parameters
!
    character(len=*), parameter :: loc = 'INTERPOLATIONS:initialize_interpolation()'
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! set timer descriptions
!
    call set_timer('interpolations:: initialization', imi)
    call set_timer('interpolations:: reconstruction', imr)
    call set_timer('interpolations:: fix positivity', imf)
    call set_timer('interpolations:: clip extrema'  , imc)

! start accounting time for module initialization/finalization
!
    call start_timer(imi)
#endif /* PROFILE */

! obtain the user defined interpolation methods and coefficients
!
    call get_parameter_string ("reconstruction"      , sreconstruction)
    call get_parameter_string ("limiter"             , tlimiter       )
    call get_parameter_string ("fix_positivity"      , positivity_fix )
    call get_parameter_string ("clip_extrema"        , clip_extrema   )
    call get_parameter_string ("extrema_limiter"     , climiter       )
    call get_parameter_string ("prolongation_limiter", plimiter       )
    call get_parameter_string ("mlp_limiting"        , mlp_limiting   )
    call get_parameter_integer("nghosts"             , ng             )
    call get_parameter_integer("ngp"                 , ngp            )
    call get_parameter_real   ("sgp"                 , sgp            )
    call get_parameter_real   ("eps"                 , eps            )
    call get_parameter_real   ("limo3_rad"           , rad            )
    call get_parameter_real   ("kappa"               , kappa          )
    call get_parameter_real   ("kbeta"               , kbeta          )
    call get_parameter_real   ("ppm_const"           , ppm_const      )
    call get_parameter_real   ("cfl"                 , cfl            )

! calculate κ = (1 - ν) / ν
!
    kappa = min(kappa, (1.0d+00 - cfl) / cfl)

! calculate mgp
!
    mgp   = (ngp - 1) / 2
    dgp   = ngp**NDIMS

! select the reconstruction method
!
    select case(trim(sreconstruction))
    case ("tvd", "TVD")
      name_rec           =  "2nd order TVD"
      interfaces         => interfaces_tvd
      reconstruct_states => reconstruct_tvd
      if (verbose .and. ng < 2)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 2)."
    case ("weno3", "WENO3")
      name_rec           =  "3rd order WENO"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_weno3
      if (verbose .and. ng < 2)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 2)."
    case ("limo3", "LIMO3", "LimO3")
      name_rec           =  "3rd order logarithmic limited"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_limo3
      if (verbose .and. ng < 2)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 2)."
      eps = max(1.0d-12, eps)
    case ("ppm", "PPM")
      name_rec           =  "3rd order PPM"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_ppm
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("weno5z", "weno5-z", "WENO5Z", "WENO5-Z")
      name_rec           =  "5th order WENO-Z (Borges et al. 2008)"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_weno5z
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("weno5yc", "weno5-yc", "WENO5YC", "WENO5-YC")
      name_rec           =  "5th order WENO-YC (Yamaleev & Carpenter 2009)"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_weno5yc
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("weno5ns", "weno5-ns", "WENO5NS", "WENO5-NS")
      name_rec           =  "5th order WENO-NS (Ha et al. 2013)"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_weno5ns
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crweno5z", "crweno5-z", "CRWENO5Z", "CRWENO5-Z")
      name_rec           =  "5th order Compact WENO-Z"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crweno5z
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crweno5yc", "crweno5-yc", "CRWENO5YC", "CRWENO5-YC")
      name_rec           =  "5th order Compact WENO-YC"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crweno5yc
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crweno5ns", "crweno5-ns", "CRWENO5NS", "CRWENO5-NS")
      name_rec           =  "5th order Compact WENO-NS"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crweno5ns
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("mp5", "MP5")
      name_rec           =  "5th order Monotonicity Preserving"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_mp5
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("mp7", "MP7")
      name_rec           =  "7th order Monotonicity Preserving"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_mp7
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crmp5", "CRMP5")
      name_rec           =  "5th order Compact Monotonicity Preserving"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crmp5
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crmp5l", "crmp5ld", "CRMP5L", "CRMP5LD")
      name_rec           =  "5th order Low-Dissipation Compact Monotonicity Preserving"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crmp5ld
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("crmp7", "CRMP7")
      name_rec           =  "7th order Compact Monotonicity Preserving"
      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_crmp7
      if (verbose .and. ng < 4)                                                &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "Increase the number of ghost cells (at least 4)."
    case ("gp", "GP")
      write(stmp, '(f16.1)') sgp
      write(name_rec, '("Gaussian Process (",i1,"-point, δ=",a,")")') ngp      &
                                                         , trim(adjustl(stmp))

! allocate the Gaussian process reconstruction matrix and position vector
!
      allocate(cgp(ngp))

! prepare matrix coefficients
!
      call prepare_gp(iret)
      if (iret > 0) return

      interfaces         => interfaces_dir
      reconstruct_states => reconstruct_gp
      if (verbose .and. 2 * ng <= ngp - 1)                                     &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                 , "Increase the number of ghost cells (at least (ngp+1)/2)."
      if (verbose .and. mod(ngp,2) == 0)                                       &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                      , "The parameter ngp has to be integer with odd value."
    case ("mgp", "MGP")
      write(stmp, '(f16.1)') sgp
      write(name_rec, &
          '("Multidimensional Gaussian Process (",i1,"-point, δ=",a,")")')     &
                                                      ngp, trim(adjustl(stmp))

! allocate the Gaussian process reconstruction matrix and position vector
!
      allocate(ugp(dgp,2,NDIMS))

! prepare matrix coefficients
!
      call prepare_mgp(iret)
      if (iret > 0) return

      interfaces         => interfaces_mgp
      if (verbose .and. 2 * ng <= ngp - 1)                                     &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                 , "Increase the number of ghost cells (at least (ngp+1)/2)."
      if (verbose .and. mod(ngp,2) == 0)                                       &
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                      , "The parameter ngp has to be integer with odd value."
    case default
      if (verbose) then
        write(error_unit,"('[',a,']: ',a)") trim(loc)                          &
                        , "The selected reconstruction method is not " //      &
                          "implemented: " // trim(sreconstruction)
        iret = 200
        return
      end if
    end select

! select the TVD limiter
!
    select case(trim(tlimiter))
    case ("mm", "minmod")
      name_tlim          =  "minmod"
      limiter_tvd        => limiter_minmod
    case ("mc", "monotonized_central")
      name_tlim          =  "monotonized central"
      limiter_tvd        => limiter_monotonized_central
    case ("sb", "superbee")
      name_tlim          =  "superbee"
      limiter_tvd        => limiter_superbee
    case ("vl", "vanleer")
      name_tlim          =  "van Leer"
      limiter_tvd        => limiter_vanleer
    case ("va", "vanalbada")
      name_tlim          =  "van Albada"
      limiter_tvd        => limiter_vanalbada
    case default
      name_tlim          =  "zero derivative"
      limiter_tvd        => limiter_zero
    end select

! select the prolongation limiter
!
    select case(trim(plimiter))
    case ("mm", "minmod")
      name_plim          =  "minmod"
      limiter_prol       => limiter_minmod
    case ("mc", "monotonized_central")
      name_plim          =  "monotonized central"
      limiter_prol       => limiter_monotonized_central
    case ("sb", "superbee")
      name_plim          =  "superbee"
      limiter_prol       => limiter_superbee
    case ("vl", "vanleer")
      name_plim          =  "van Leer"
      limiter_prol       => limiter_vanleer
    case default
      name_plim          =  "zero derivative"
      limiter_prol       => limiter_zero
    end select

! select the clipping limiter
!
    select case(trim(climiter))
    case ("mm", "minmod")
      name_clim          =  "minmod"
      limiter_clip       => limiter_minmod
    case ("mc", "monotonized_central")
      name_clim          =  "monotonized central"
      limiter_clip       => limiter_monotonized_central
    case ("sb", "superbee")
      name_clim          =  "superbee"
      limiter_clip       => limiter_superbee
    case ("vl", "vanleer")
      name_clim          =  "van Leer"
      limiter_clip       => limiter_vanleer
    case default
      name_clim          =  "zero derivative"
      limiter_clip       => limiter_zero
    end select

! check additional reconstruction limiting
!
    select case(trim(mlp_limiting))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      mlp = .true.
    case default
      mlp = .false.
    end select
    select case(trim(positivity_fix))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      positivity = .true.
    case default
      positivity = .false.
    end select
    select case(trim(clip_extrema))
    case ("on", "ON", "t", "T", "y", "Y", "true", "TRUE", "yes", "YES")
      clip = .true.
    case default
      clip = .false.
    end select

! print informations about the reconstruction methods and parameters
!
    if (verbose) then

      write (*,"(4x,a14, 9x,'=',1x,a)") "reconstruction"      , trim(name_rec)
      write (*,"(4x,a11,12x,'=',1x,a)") "TVD limiter"         , trim(name_tlim)
      write (*,"(4x,a20, 3x,'=',1x,a)") "prolongation limiter", trim(name_plim)
      write (*,"(4x,a12,11x,'=',1x,a)") "MLP limiting"        , trim(mlp_limiting)
      write (*,"(4x,a14, 9x,'=',1x,a)") "fix positivity"      , trim(positivity_fix)
      write (*,"(4x,a12,11x,'=',1x,a)") "clip extrema"        , trim(clip_extrema)
      if (clip) then
        write (*,"(4x,a15,8x,'=',1x,a)") "extrema limiter", trim(name_clim)
      end if

    end if

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine initialize_interpolations
!
!===============================================================================
!
! subroutine FINALIZE_INTERPOLATIONS:
! ----------------------------------
!
!   Subroutine finalizes the interpolation module by releasing all memory used
!   by its module variables.
!
!
!===============================================================================
!
  subroutine finalize_interpolations(iret)

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

! deallocate Gaussian process reconstruction coefficient vector if used
!
    if (allocated(cgp)) deallocate(cgp)
    if (allocated(ugp)) deallocate(ugp)

! release the procedure pointers
!
    nullify(reconstruct_states)
    nullify(limiter_tvd)
    nullify(limiter_prol)
    nullify(limiter_clip)

#ifdef PROFILE
! stop accounting time for module initialization/finalization
!
    call stop_timer(imi)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine finalize_interpolations
!
!===============================================================================
!
! subroutine PREPARE_MGP:
! ---------------------
!
!   Subroutine prepares matrixes for the multidimensional
!   Gaussian Process (GP) method.
!
!===============================================================================
!
  subroutine prepare_mgp(iret)

! include external procedures
!
    use algebra        , only : max_real_kind, invert
    use constants      , only : pi
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret

! local variables
!
    logical                  :: flag
    integer                  :: i, j, i1, j1, k1, i2, j2, k2
    real(kind=max_real_kind) :: sig, fc, fx, fy, fz, xl, xr, yl, yr, zl, zr

! local arrays for derivatives
!
    real(kind=max_real_kind), dimension(:,:)  , allocatable :: cov, inv
    real(kind=max_real_kind), dimension(:,:,:), allocatable :: xgp

! local parameters
!
    character(len=*), parameter :: loc = 'INTERPOLATIONS::prepare_mgp()'
!
!-------------------------------------------------------------------------------
!
! calculate normal distribution sigma
!
    sig = sqrt(2.0d+00) * sgp

! allocate the convariance matrix and interpolation position vector
!
    allocate(cov(dgp,dgp))
    allocate(inv(dgp,dgp))
    allocate(xgp(dgp,2,NDIMS))

! prepare the covariance matrix
!
    fc = 0.5d+00 * sqrt(pi) * sig
    i = 0
#if NDIMS == 3
    do k1 = - mgp, mgp
#endif /* NDIMS == 3 */
      do j1 = - mgp, mgp
        do i1 = - mgp, mgp
          i = i + 1
          j = 0
#if NDIMS == 3
          do k2 = - mgp, mgp
#endif /* NDIMS == 3 */
            do j2 = - mgp, mgp
              do i2 = - mgp, mgp
                j = j + 1

                xl = (1.0d+00 * (i1 - i2) - 0.5d+00) / sig
                xr = (1.0d+00 * (i1 - i2) + 0.5d+00) / sig
                yl = (1.0d+00 * (j1 - j2) - 0.5d+00) / sig
                yr = (1.0d+00 * (j1 - j2) + 0.5d+00) / sig
#if NDIMS == 3
                zl = (1.0d+00 * (k1 - k2) - 0.5d+00) / sig
                zr = (1.0d+00 * (k1 - k2) + 0.5d+00) / sig
#endif /* NDIMS == 3 */

                cov(i,j) = fc * (erf(xr) - erf(xl)) * (erf(yr) - erf(yl))
#if NDIMS == 3
                cov(i,j) = cov(i,j) * (erf(zr) - erf(zl))
#endif /* NDIMS == 3 */
              end do
            end do
          end do
        end do
#if NDIMS == 3
      end do
    end do
#endif /* NDIMS == 3 */

! prepare the interpolation position vector
!
    i = 0
#if NDIMS == 3
    do k1 = - mgp, mgp
#endif /* NDIMS == 3 */
      do j1 = - mgp, mgp
        do i1 = - mgp, mgp
          i = i + 1

          xl = (1.0d+00 * i1 - 0.5d+00) / sig
          xr = (1.0d+00 * i1 + 0.5d+00) / sig
          yl = (1.0d+00 * j1 - 0.5d+00) / sig
          yr = (1.0d+00 * j1 + 0.5d+00) / sig
#if NDIMS == 3
          zl = (1.0d+00 * k1 - 0.5d+00) / sig
          zr = (1.0d+00 * k1 + 0.5d+00) / sig
#endif /* NDIMS == 3 */

          fx = erf(xr) - erf(xl)
          fy = erf(yr) - erf(yl)
#if NDIMS == 3
          fz = erf(zr) - erf(zl)

          xgp(i,1,1) = exp(- xl**2) * fy * fz
          xgp(i,2,1) = exp(- xr**2) * fy * fz
          xgp(i,1,2) = exp(- yl**2) * fx * fz
          xgp(i,2,2) = exp(- yr**2) * fx * fz
          xgp(i,1,3) = exp(- zl**2) * fx * fy
          xgp(i,2,3) = exp(- zr**2) * fx * fy
#else /* NDIMS == 3 */
          xgp(i,1,1) = exp(- xl**2) * fy
          xgp(i,2,1) = exp(- xr**2) * fy
          xgp(i,1,2) = exp(- yl**2) * fx
          xgp(i,2,2) = exp(- yr**2) * fx
#endif /* NDIMS == 3 */
        end do
      end do
#if NDIMS == 3
    end do
#endif /* NDIMS == 3 */

! invert the matrix
!
    call invert(dgp, cov(1:dgp,1:dgp), inv(1:dgp,1:dgp), flag)

! prepare the interpolation coefficients vector
!
    do j = 1, NDIMS
      do i = 1, 2
        ugp(1:dgp,i,j) = matmul(xgp(1:dgp,i,j), inv(1:dgp,1:dgp))
      end do
    end do

! deallocate the convariance matrix and interpolation position vector
!
    deallocate(cov)
    deallocate(inv)
    deallocate(xgp)

! check if the matrix was inverted successfully
!
    if (.not. flag) then
      write(error_unit,"('[',a,']: ',a)") trim(loc)                            &
                      , "Could not invert the covariance matrix!"
      iret = 201
    end if

!-------------------------------------------------------------------------------
!
  end subroutine prepare_mgp
!
!===============================================================================
!
! subroutine INTERFACES_TVD:
! -------------------------
!
!   Subroutine reconstructs both side interfaces of variable using TVD methods.
!
!   Arguments:
!
!     positive - the variable positivity flag;
!     h        - the spatial step;
!     q        - the variable array;
!     qi       - the array of reconstructed interfaces (2 in each direction);
!
!===============================================================================
!
  subroutine interfaces_tvd(positive, h, q, qi)

! include external procedures
!
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb , ie , je , ke
    use coordinates    , only : ibl, jbl, kbl, ieu, jeu, keu
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical                                  , intent(in)  :: positive
    real(kind=8), dimension(NDIMS)           , intent(in)  :: h
    real(kind=8), dimension(im,jm,km)        , intent(in)  :: q
    real(kind=8), dimension(im,jm,km,2,NDIMS), intent(out) :: qi

! local variables
!
    integer :: i, im1, ip1
    integer :: j, jm1, jp1
    integer :: k, km1, kp1

! local vectors
!
    real(kind=8), dimension(NDIMS) :: dql, dqr, dq

! local parameters
!
    character(len=*), parameter :: loc = 'INTERPOLATIONS::interfaces_tvd()'
!
!-------------------------------------------------------------------------------
!
! copy ghost zones
!
    do k = 1, NDIMS
      do j = 1, 2
        qi( 1:ib, 1:jm, 1:km,j,k) = q( 1:ib, 1:jm, 1:km)
        qi(ie:im, 1:jm, 1:km,j,k) = q(ie:im, 1:jm, 1:km)
        qi(ib:ie, 1:jb, 1:km,j,k) = q(ib:ie, 1:jb, 1:km)
        qi(ib:ie,je:jm, 1:km,j,k) = q(ib:ie,je:jm, 1:km)
#if NDIMS == 3
        qi(ib:ie,jb:je, 1:kb,j,k) = q(ib:ie,jb:je, 1:kb)
        qi(ib:ie,jb:je,ke:km,j,k) = q(ib:ie,jb:je,ke:km)
#endif /* NDIMS == 3 */
      end do
    end do

! interpolate interfaces
!
    do k = kbl, keu
#if NDIMS == 3
      km1 = k - 1
      kp1 = k + 1
#endif /* NDIMS == 3 */
      do j = jbl, jeu
        jm1 = j - 1
        jp1 = j + 1
        do i = ibl, ieu
          im1 = i - 1
          ip1 = i + 1

! calculate the TVD derivatives
!
          dql(1) = q(i  ,j,k) - q(im1,j,k)
          dqr(1) = q(ip1,j,k) - q(i  ,j,k)
          dq (1) = limiter_tvd(0.5d+00, dql(1), dqr(1))

          dql(2) = q(i,j  ,k) - q(i,jm1,k)
          dqr(2) = q(i,jp1,k) - q(i,j  ,k)
          dq (2) = limiter_tvd(0.5d+00, dql(2), dqr(2))

#if NDIMS == 3
          dql(3) = q(i,j,k  ) - q(i,j,km1)
          dqr(3) = q(i,j,kp1) - q(i,j,k  )
          dq (3) = limiter_tvd(0.5d+00, dql(3), dqr(3))
#endif /* NDIMS == 3 */

! limit the derivatives if they produce negative interpolation for positive
! variables
!
          if (positive) then
            if (q(i,j,k) > 0.0d+00) then
              do while (q(i,j,k) <= sum(abs(dq(1:NDIMS))))
                dq(:) = 0.5d+00 * dq(:)
              end do
            else
              write(error_unit,"('[',a,']: ',a,3i4,a)") trim(loc)              &
                    , "Positive variable is not positive at ( ", i, j, k, " )"
              dq(:) = 0.0d+00
            end if
          end if

! interpolate states
!
          qi(i  ,j,k,1,1) = q(i,j,k) + dq(1)
          qi(im1,j,k,2,1) = q(i,j,k) - dq(1)

          qi(i,j  ,k,1,2) = q(i,j,k) + dq(2)
          qi(i,jm1,k,2,2) = q(i,j,k) - dq(2)

#if NDIMS == 3
          qi(i,j,k  ,1,3) = q(i,j,k) + dq(3)
          qi(i,j,km1,2,3) = q(i,j,k) - dq(3)
#endif /* NDIMS == 3 */

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
    end do ! k = kbl, keu

! apply Multi-dimensional Limiting Process
!
    if (mlp) call mlp_limiting(q(:,:,:), qi(:,:,:,:,:))

!-------------------------------------------------------------------------------
!
  end subroutine interfaces_tvd
!
!===============================================================================
!
! subroutine INTERFACES_DIR:
! -------------------------
!
!   Subroutine reconstructs both side interfaces of variable separately
!   along each direction.
!
!   Arguments:
!
!     positive - the variable positivity flag;
!     h        - the spatial step;
!     q        - the variable array;
!     qi       - the array of reconstructed interfaces (2 in each direction);
!
!===============================================================================
!
  subroutine interfaces_dir(positive, h, q, qi)

! include external procedures
!
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb , ie , je , ke
    use coordinates    , only : ibl, jbl, kbl, ieu, jeu, keu

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical                                  , intent(in)  :: positive
    real(kind=8), dimension(NDIMS)           , intent(in)  :: h
    real(kind=8), dimension(im,jm,km)        , intent(in)  :: q
    real(kind=8), dimension(im,jm,km,2,NDIMS), intent(out) :: qi

! local variables
!
    integer :: i, j, k
!
!-------------------------------------------------------------------------------
!
! copy ghost zones
!
    do k = 1, NDIMS
      do j = 1, 2
        qi( 1:ib, 1:jm, 1:km,j,k) = q( 1:ib, 1:jm, 1:km)
        qi(ie:im, 1:jm, 1:km,j,k) = q(ie:im, 1:jm, 1:km)
        qi(ib:ie, 1:jb, 1:km,j,k) = q(ib:ie, 1:jb, 1:km)
        qi(ib:ie,je:jm, 1:km,j,k) = q(ib:ie,je:jm, 1:km)
#if NDIMS == 3
        qi(ib:ie,jb:je, 1:kb,j,k) = q(ib:ie,jb:je, 1:kb)
        qi(ib:ie,jb:je,ke:km,j,k) = q(ib:ie,jb:je,ke:km)
#endif /* NDIMS == 3 */
      end do
    end do

! interpolate interfaces
!
    do k = kbl, keu
      do j = jbl, jeu
        call reconstruct(im, h(1), q(1:im,j,k)                                 &
                                         , qi(1:im,j,k,1,1), qi(1:im,j,k,2,1))
      end do ! j = jbl, jeu
      do i = ibl, ieu
        call reconstruct(jm, h(2), q(i,1:jm,k)                                 &
                                         , qi(i,1:jm,k,1,2), qi(i,1:jm,k,2,2))
      end do ! i = ibl, ieu
    end do ! k = kbl, keu
#if NDIMS == 3
    do j = jbl, jeu
      do i = ibl, ieu
        call reconstruct(km, h(3), q(i,j,1:km)                                 &
                                         , qi(i,j,1:km,1,3), qi(i,j,1:km,2,3))
      end do ! i = ibl, ieu
    end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

! make sure the interface states are positive for positive variables
!
    if (positive) then

      do k = kbl, keu
        do j = jbl, jeu
          call fix_positivity(im, q(1:im,j,k)                                  &
                                         , qi(1:im,j,k,1,1), qi(1:im,j,k,2,1))
        end do ! j = jbl, jeu
        do i = ibl, ieu
          call fix_positivity(jm, q(i,1:jm,k)                                  &
                                         , qi(i,1:jm,k,1,2), qi(i,1:jm,k,2,2))
        end do ! i = ibl, ieu
      end do ! k = kbl, keu
#if NDIMS == 3
      do j = jbl, jeu
        do i = ibl, ieu
          call fix_positivity(km, q(i,j,1:km)                                  &
                                         , qi(i,j,1:km,1,3), qi(i,j,1:km,2,3))
        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
#endif /* NDIMS == 3 */

    end if

! apply Multi-dimensional Limiting Process
!
    if (mlp) call mlp_limiting(q(:,:,:), qi(:,:,:,:,:))

!-------------------------------------------------------------------------------
!
  end subroutine interfaces_dir
!
!===============================================================================
!
! subroutine INTERFACES_MGP:
! -------------------------
!
!   Subroutine reconstructs both side interfaces of variable using
!   multidimensional Gaussian Process method.
!
!   Arguments:
!
!     positive - the variable positivity flag;
!     h        - the spatial step;
!     q        - the variable array;
!     qi       - the array of reconstructed interfaces (2 in each direction);
!
!===============================================================================
!
  subroutine interfaces_mgp(positive, h, q, qi)

! include external procedures
!
    use coordinates    , only : im , jm , km
    use coordinates    , only : ib , jb , kb , ie , je , ke
    use coordinates    , only : ibl, jbl, kbl, ieu, jeu, keu
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical                                  , intent(in)  :: positive
    real(kind=8), dimension(NDIMS)           , intent(in)  :: h
    real(kind=8), dimension(im,jm,km)        , intent(in)  :: q
    real(kind=8), dimension(im,jm,km,2,NDIMS), intent(out) :: qi

! local variables
!
    logical       :: flag
    integer       :: i, il, iu, im1, ip1
    integer       :: j, jl, ju, jm1, jp1
    integer       :: k, kl, ku, km1, kp1

! local arrays for derivatives
!
    real(kind=8), dimension(NDIMS) :: dql, dqr, dq
    real(kind=8), dimension(dgp)   :: u

! local parameters
!
    character(len=*), parameter :: loc = 'INTERPOLATIONS::interfaces_mgp()'
!
!-------------------------------------------------------------------------------
!
! copy ghost zones
!
    do k = 1, NDIMS
      do j = 1, 2
        qi( 1:ib, 1:jm, 1:km,j,k) = q( 1:ib, 1:jm, 1:km)
        qi(ie:im, 1:jm, 1:km,j,k) = q(ie:im, 1:jm, 1:km)
        qi(ib:ie, 1:jb, 1:km,j,k) = q(ib:ie, 1:jb, 1:km)
        qi(ib:ie,je:jm, 1:km,j,k) = q(ib:ie,je:jm, 1:km)
#if NDIMS == 3
        qi(ib:ie,jb:je, 1:kb,j,k) = q(ib:ie,jb:je, 1:kb)
        qi(ib:ie,jb:je,ke:km,j,k) = q(ib:ie,jb:je,ke:km)
#endif /* NDIMS == 3 */
      end do
    end do

! interpolate interfaces using precomputed interpolation vectors
!
    do k = kbl, keu
#if NDIMS == 3
      kl  = k - mgp
      ku  = k + mgp
      km1 = k - 1
      kp1 = k + 1
#endif /* NDIMS == 3 */
      do j = jbl, jeu
        jl  = j - mgp
        ju  = j + mgp
        jm1 = j - 1
        jp1 = j + 1
        do i = ibl, ieu
          il  = i - mgp
          iu  = i + mgp
          im1 = i - 1
          ip1 = i + 1

#if NDIMS == 3
          u(:) = reshape(q(il:iu,jl:ju,kl:ku), (/ dgp /)) - q(i,j,k)

          qi(i  ,j,k,1,1) = sum(ugp(1:dgp,1,1) * u(1:dgp)) + q(i,j,k)
          qi(im1,j,k,2,1) = sum(ugp(1:dgp,2,1) * u(1:dgp)) + q(i,j,k)
          qi(i,j  ,k,1,2) = sum(ugp(1:dgp,1,2) * u(1:dgp)) + q(i,j,k)
          qi(i,jm1,k,2,2) = sum(ugp(1:dgp,2,2) * u(1:dgp)) + q(i,j,k)
          qi(i,j,k  ,1,3) = sum(ugp(1:dgp,1,3) * u(1:dgp)) + q(i,j,k)
          qi(i,j,km1,2,3) = sum(ugp(1:dgp,2,3) * u(1:dgp)) + q(i,j,k)
#else /* NDIMS == 3 */
          u(:) = reshape(q(il:iu,jl:ju,k    ), (/ dgp /)) - q(i,j,k)

          qi(i  ,j,k,1,1) = sum(ugp(1:dgp,1,1) * u(1:dgp)) + q(i,j,k)
          qi(im1,j,k,2,1) = sum(ugp(1:dgp,2,1) * u(1:dgp)) + q(i,j,k)
          qi(i,j  ,k,1,2) = sum(ugp(1:dgp,1,2) * u(1:dgp)) + q(i,j,k)
          qi(i,jm1,k,2,2) = sum(ugp(1:dgp,2,2) * u(1:dgp)) + q(i,j,k)
#endif /* NDIMS == 3 */

! if the interpolation is not monotonic, apply a TVD slope
!
          flag =           ((qi(i  ,j,k,1,1) - q(ip1,j,k))                     &
                          * (qi(i  ,j,k,1,1) - q(i  ,j,k)) > eps)
          flag = flag .or. ((qi(im1,j,k,2,1) - q(im1,j,k))                     &
                          * (qi(im1,j,k,2,1) - q(i  ,j,k)) > eps)
          flag = flag .or. ((qi(i,j  ,k,1,2) - q(i,jp1,k))                     &
                          * (qi(i,j  ,k,1,2) - q(i,j  ,k)) > eps)
          flag = flag .or. ((qi(i,jm1,k,2,2) - q(i,jm1,k))                     &
                          * (qi(i,jm1,k,2,2) - q(i,j  ,k)) > eps)
#if NDIMS == 3
          flag = flag .or. ((qi(i,j,k  ,1,3) - q(i,j,kp1))                     &
                          * (qi(i,j,k  ,1,3) - q(i,j,k  )) > eps)
          flag = flag .or. ((qi(i,j,km1,2,3) - q(i,j,km1))                     &
                          * (qi(i,j,km1,2,3) - q(i,j,k  )) > eps)
#endif /* NDIMS == 3 */

          if (flag) then

! calculate the TVD derivatives
!
            dql(1) = q(i  ,j,k) - q(im1,j,k)
            dqr(1) = q(ip1,j,k) - q(i  ,j,k)
            dq (1) = limiter_tvd(0.5d+00, dql(1), dqr(1))

            dql(2) = q(i,j  ,k) - q(i,jm1,k)
            dqr(2) = q(i,jp1,k) - q(i,j  ,k)
            dq (2) = limiter_tvd(0.5d+00, dql(2), dqr(2))

#if NDIMS == 3
            dql(3) = q(i,j,k  ) - q(i,j,km1)
            dqr(3) = q(i,j,kp1) - q(i,j,k  )
            dq (3) = limiter_tvd(0.5d+00, dql(3), dqr(3))
#endif /* NDIMS == 3 */

! limit the derivatives if they produce negative interpolation for positive
! variables
!
            if (positive) then
              if (q(i,j,k) > 0.0d+00) then
                do while (q(i,j,k) <= sum(abs(dq(1:NDIMS))))
                  dq(:) = 0.5d+00 * dq(:)
                end do
              else
                write(error_unit,"('[',a,']: ',a,3i4,a)") trim(loc)            &
                    , "Positive variable is not positive at ( ", i, j, k, " )"
                dq(:) = 0.0d+00
              end if
            end if

! interpolate states
!
            qi(i  ,j,k,1,1) = q(i,j,k) + dq(1)
            qi(im1,j,k,2,1) = q(i,j,k) - dq(1)

            qi(i,j  ,k,1,2) = q(i,j,k) + dq(2)
            qi(i,jm1,k,2,2) = q(i,j,k) - dq(2)

#if NDIMS == 3
            qi(i,j,k  ,1,3) = q(i,j,k) + dq(3)
            qi(i,j,km1,2,3) = q(i,j,k) - dq(3)
#endif /* NDIMS == 3 */

          end if

        end do ! i = ibl, ieu
      end do ! j = jbl, jeu
    end do ! k = kbl, keu

! apply Multi-dimensional Limiting Process
!
    if (mlp) call mlp_limiting(q(:,:,:), qi(:,:,:,:,:))

!-------------------------------------------------------------------------------
!
  end subroutine interfaces_mgp
!
!===============================================================================
!
! subroutine MLP_LIMITING:
! -----------------------
!
!   Subroutine applies the multi-dimensional limiting process to
!   the reconstructed states in order to control oscillations due to
!   the Gibbs phenomena near discontinuities.
!
!   Arguments:
!
!     q        - the variable array;
!     qi       - the array of reconstructed interfaces (2 in each direction);
!
!   References:
!
!     [1] Gerlinger, P.,
!         "Multi-dimensional limiting for high-order schemes including
!          turbulence and combustion",
!         Journal of Computational Physics,
!         2012, vol. 231, pp. 2199-2228,
!         http://dx.doi.org/10.1016/j.jcp.2011.10.024
!
!===============================================================================
!
  subroutine mlp_limiting(q, qi)

! include external procedures
!
    use coordinates    , only : im , jm , km

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    real(kind=8), dimension(im,jm,km)        , intent(in)    :: q
    real(kind=8), dimension(im,jm,km,2,NDIMS), intent(inout) :: qi

! local variables
!
    integer      :: i, im1, ip1
    integer      :: j, jm1, jp1
    integer      :: k, km1, kp1
    integer      :: m
#if NDIMS == 3
    integer      :: n, np1, np2
#endif /* NDIMS == 3 */
    real(kind=8) :: qmn, qmx, dqc
    real(kind=8) :: fc, fl

! local vectors
!
    real(kind=8), dimension(NDIMS) :: dql, dqr, dq
    real(kind=8), dimension(NDIMS) :: dqm, ap
#if NDIMS == 3
    real(kind=8), dimension(NDIMS) :: hh, uu
#endif /* NDIMS == 3 */
!
!-------------------------------------------------------------------------------
!
    do k = 1, km
#if NDIMS == 3
      km1 = max( 1, k - 1)
      kp1 = min(km, k + 1)
#endif /* NDIMS == 3 */
      do j = 1, jm
        jm1 = max( 1, j - 1)
        jp1 = min(jm, j + 1)
        do i = 1, im
          im1 = max( 1, i - 1)
          ip1 = min(im, i + 1)

! calculate the minmod TVD derivatives
!
          dql(1) = q(i  ,j,k) - q(im1,j,k)
          dqr(1) = q(ip1,j,k) - q(i  ,j,k)
          dq (1) = limiter_minmod(0.5d+00, dql(1), dqr(1))

          dql(2) = q(i,j  ,k) - q(i,jm1,k)
          dqr(2) = q(i,jp1,k) - q(i,j  ,k)
          dq (2) = limiter_minmod(0.5d+00, dql(2), dqr(2))

#if NDIMS == 3
          dql(3) = q(i,j,k  ) - q(i,j,km1)
          dqr(3) = q(i,j,kp1) - q(i,j,k  )
          dq (3) = limiter_minmod(0.5d+00, dql(3), dqr(3))
#endif /* NDIMS == 3 */

! calculate dqc
!
#if NDIMS == 3
          qmn = minval(q(im1:ip1,jm1:jp1,km1:kp1))
          qmx = maxval(q(im1:ip1,jm1:jp1,km1:kp1))
#else /* NDIMS == 3 */
          qmn = minval(q(im1:ip1,jm1:jp1,k))
          qmx = maxval(q(im1:ip1,jm1:jp1,k))
#endif /* NDIMS == 3 */
          dqc = min(qmx - q(i,j,k), q(i,j,k) - qmn)

! if needed, apply the multi-dimensional limiting process
!
          if (sum(abs(dq(1:NDIMS))) > dqc) then

            dqm(1) = abs(q(ip1,j,k) - q(im1,j,k))
            dqm(2) = abs(q(i,jp1,k) - q(i,jm1,k))
#if NDIMS == 3
            dqm(3) = abs(q(i,j,kp1) - q(i,j,km1))
#endif /* NDIMS == 3 */

            fc = dqc / max(1.0d-16, sum(abs(dqm(1:NDIMS))))
            do m = 1, NDIMS
              ap(m) = fc * abs(dqm(m))
            end do

#if NDIMS == 3
            do n = 1, NDIMS
              hh(n)   = max(ap(n) - abs(dq(n)), 0.0d+00)
              np1     = mod(n    , NDIMS) + 1
              np2     = mod(n + 1, NDIMS) + 1
              uu(n)   = ap(n)   -           hh(n)
              uu(np1) = ap(np1) + 0.5d+00 * hh(n)
              uu(np2) = ap(np2) + 0.5d+00 * hh(n)
              fc      = hh(n) / (hh(n) + 1.0d-16)
              fl      = fc * (max(ap(np1) - abs(dq(np1)), 0.0d+00)             &
                            - max(ap(np2) - abs(dq(np2)), 0.0d+00))
              ap(n  ) = uu(n)
              ap(np1) = uu(np1) - fl
              ap(np2) = uu(np2) + fl
            end do
#else /* NDIMS == 3 */
            fl     = max(ap(1) - abs(dq(1)), 0.0d+00)                          &
                   - max(ap(2) - abs(dq(2)), 0.0d+00)
            ap(1)  = ap(1) - fl
            ap(2)  = ap(2) + fl
#endif /* NDIMS == 3 */

            do m = 1, NDIMS
              dq(m) = sign(ap(m), dq(m))
            end do

          end if

! update the interpolated states
!
          dql(1) = qi(i  ,j,k,1,1) - q(i,j,k)
          dqr(1) = qi(im1,j,k,2,1) - q(i,j,k)
          if (max(abs(dql(1)), abs(dqr(1))) > abs(dq(1))) then
            qi(i  ,j,k,1,1) = q(i,j,k) + dq(1)
            qi(im1,j,k,2,1) = q(i,j,k) - dq(1)
          end if

          dql(2) = qi(i,j  ,k,1,2) - q(i,j,k)
          dqr(2) = qi(i,jm1,k,2,2) - q(i,j,k)
          if (max(abs(dql(2)), abs(dqr(2))) > abs(dq(2))) then
            qi(i,j  ,k,1,2) = q(i,j,k) + dq(2)
            qi(i,jm1,k,2,2) = q(i,j,k) - dq(2)
          end if

#if NDIMS == 3
          dql(3) = qi(i,j,k  ,1,3) - q(i,j,k)
          dqr(3) = qi(i,j,km1,2,3) - q(i,j,k)
          if (max(abs(dql(3)), abs(dqr(3))) > abs(dq(3))) then
            qi(i,j,k  ,1,3) = q(i,j,k) + dq(3)
            qi(i,j,km1,2,3) = q(i,j,k) - dq(3)
          end if
#endif /* NDIMS == 3 */

        end do ! i = 1, im
      end do ! j = 1, jm
    end do ! k = 1, km

!-------------------------------------------------------------------------------
!
  end subroutine mlp_limiting
!
!===============================================================================
!
! subroutine RECONSTRUCT:
! ----------------------
!
!   Subroutine calls a reconstruction procedure, depending on the compilation
!   flag SPACE, in order to interpolate the left and right states from their
!   cell integrals.  These states are required by any approximate Riemann
!   solver.
!
!   Arguments:
!
!     n  - the length of the input vector;
!     h  - the spatial step; this is required for some reconstruction methods;
!     f  - the input vector of cell averaged values;
!     fl - the left side state reconstructed for location (i+1/2);
!     fr - the right side state reconstructed for location (i+1/2);
!
!===============================================================================
!
  subroutine reconstruct(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr
!
!-------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for reconstruction
!
    call start_timer(imr)
#endif /* PROFILE */

! reconstruct the states using the selected subroutine
!
    call reconstruct_states(n, h, f(:), fl(:), fr(:))

! correct the reconstruction near extrema by clipping them in order to improve
! the stability of scheme
!
    if (clip) call clip_extrema(n, f(:), fl(:), fr(:))

#ifdef PROFILE
! stop accounting time for reconstruction
!
    call stop_timer(imr)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct
!
!===============================================================================
!
! subroutine RECONSTRUCT_TVD:
! --------------------------
!
!   Subroutine reconstructs the interface states using the second order TVD
!   method with a selected limiter.
!
!   Arguments are described in subroutine reconstruct().
!
!
!===============================================================================
!
  subroutine reconstruct_tvd(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      ::  i, im1, ip1
    real(kind=8) :: df, dfl, dfr
!
!-------------------------------------------------------------------------------
!
! calculate the left- and right-side interface interpolations
!
    do i = 2, n - 1

! calculate left and right indices
!
      im1     = i - 1
      ip1     = i + 1

! calculate left and right side derivatives
!
      dfl     = f(i  ) - f(im1)
      dfr     = f(ip1) - f(i  )

! obtain the TVD limited derivative
!
      df      = limiter_tvd(0.5d+00, dfl, dfr)

! update the left and right-side interpolation states
!
      fl(i  ) = f(i) + df
      fr(im1) = f(i) - df

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_tvd
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO3:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the third order
!   Weighted Essentially Non-Oscillatory (WENO) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Yamaleev & Carpenter, 2009, J. Comput. Phys., 228, 3025
!
!===============================================================================
!
  subroutine reconstruct_weno3(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: bp, bm, ap, am, wp, wm, ww
    real(kind=8) :: dfl, dfr, df, fp, fm, fc, h2

! selection weights
!
    real(kind=8), parameter :: dp = 2.0d+00 / 3.0d+00, dm = 1.0d+00 / 3.0d+00
!
!-------------------------------------------------------------------------------
!
! prepare common parameters
!
    h2 = h * h

! iterate along the vector
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1     = i - 1
      ip1     = i + 1

! calculate the left and right derivatives
!
      dfl     = f(i  ) - f(im1)
      dfr     = f(ip1) - f(i  )

! calculate coefficient omega
!
      ww      = (dfr - dfl)**2

! calculate corresponding betas
!
      bp      = dfr * dfr
      bm      = dfl * dfl

! calculate improved alphas
!
      ap      = 1.0d+00 + ww / (bp + h2)
      am      = 1.0d+00 + ww / (bm + h2)

! calculate weights
!
      wp      = dp * ap
      wm      = dm * am
      ww      = 2.0d+00 * (wp + wm)

! calculate central interpolation
!
      fp      =   f(i  ) +         f(ip1)

! calculate left side interpolation
!
      fm      = - f(im1) + 3.0d+00 * f(i  )

! calculate the left state
!
      fl( i ) = (wp * fp + wm * fm) / ww

! calculate weights
!
      wp      = dp * am
      wm      = dm * ap
      ww      = 2.0d+00 * (wp + wm)

! calculate central interpolation
!
      fp      =   f(i  ) +         f(im1)

! calculate right side interpolation
!
      fm      = - f(ip1) + 3.0d+00 * f(i  )

! calculate the right state
!
      fr(im1) = (wp * fp + wm * fm) / ww

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno3
!
!===============================================================================
!
! subroutine RECONSTRUCT_LIMO3:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the third order method
!   with a limiter function LimO3.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Cada, M. & Torrilhon, M.,
!         "Compact third-order limiter functions for finite volume methods",
!         Journal of Computational Physics, 2009, 228, 4118-4145
!     [2] Mignone, A., Tzeferacos, P., & Bodo, G.,
!         "High-order conservative finite divergence GLM-MHD schemes for
!          cell-centered MHD",
!         Journal of Computational Physics, 2010, 229, 5896-5920
!
!===============================================================================
!
  subroutine reconstruct_limo3(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: dfl, dfr
    real(kind=8) :: th, et, f1, f2, xl, xi, rdx, rdx2
!
!-------------------------------------------------------------------------------
!
! prepare parameters
!
    rdx   = rad * h
    rdx2  = rdx * rdx

! iterate over positions and interpolate states
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! prepare left and right differences
!
      dfl = f(i  ) - f(im1)
      dfr = f(ip1) - f(i  )

! calculate the indicator function (eq. 3.17 in [1])
!
      et = (dfl * dfl + dfr * dfr) / rdx2

! the switching function (embedded in eq. 3.22 in [1], eq. 32 in [2])
!
      xi = max(0.0d+00, 0.5d+00 * min(2.0d+00, 1.0d+00 + (et - 1.0d+00) / eps))
      xl = 1.0d+00 - xi

! calculate values at i + ½
!
      if (abs(dfr) > eps) then

! calculate the slope ratio (eq. 2.8 in [1])
!
        th = dfl / dfr

! calculate the quadratic reconstruction (eq. 3.8 in [1], divided by 2)
!
        f1 = (2.0d+00 + th) / 6.0d+00

! calculate the third order limiter (eq. 3.13 in [1], cofficients divided by 2)
!
        if (th >= 0.0d+00) then
          f2 = max(0.0d+00, min(f1, th, 0.8d+00))
        else
          f2 = max(0.0d+00, min(f1, - 0.25d+00 * th))
        end if

! interpolate the left state (eq. 3.5 in [1], eq. 30 in [2])
!
        fl(i) = f(i) + dfr * (xl * f1 + xi * f2)

      else

        fl(i) = f(i)

      end if

! calculate values at i - ½
!
      if (abs(dfl) > eps) then

! calculate the slope ratio (eq. 2.8 in [1])
!
        th = dfr / dfl

! calculate the quadratic reconstruction (eq. 3.8 in [1], divided by 2)
!
        f1 = (2.0d+00 + th) / 6.0d+00

! calculate the third order limiter (eq. 3.13 in [1], cofficients divided by 2)
!
        if (th >= 0.0d+00) then
          f2 = max(0.0d+00, min(f1, th, 0.8d+00))
        else
          f2 = max(0.0d+00, min(f1, - 0.25d+00 * th))
        end if

! interpolate the right state (eq. 3.5 in [1], eq. 30 in [2])
!
        fr(im1) = f(i) - dfl * (xl * f1 + xi * f2)

      else

        fr(im1) = f(i)

      end if

    end do ! i = 2, n - 1

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_limo3
!
!===============================================================================
!
! subroutine RECONSTRUCT_PPM:
! --------------------------
!
!   Subroutine reconstructs the interface states using the third order
!   Piecewise Parabolic Method (PPM) with new limiters. This version lacks
!   the support for flattening important for keeping the spurious oscillations
!   near strong shocks/discontinuties under control.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Colella, P., Woodward, P. R.,
!         "The Piecewise Parabolic Method (PPM) for Gas-Dynamical Simulations",
!         Journal of Computational Physics, 1984, vol. 54, pp. 174-201,
!         https://dx.doi.org/10.1016/0021-9991(84)90143-8
!     [2] Colella, P., Sekora, M. D.,
!         "A limiter for PPM that preserves accuracy at smooth extrema",
!         Journal of Computational Physics, 2008, vol. 227, pp. 7069-7076,
!         https://doi.org/10.1016/j.jcp.2008.03.034
!
!===============================================================================
!
  subroutine reconstruct_ppm(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    logical      :: cm, cp, ext
    integer      :: i, im1, ip1, im2
    real(kind=8) :: df2c, df2m, df2p, df2s, df2l
    real(kind=8) :: dfm, dfp, dcm, dcp
    real(kind=8) :: alm, alp, amx, sgn, del

! local arrays
!
    real(kind=8), dimension(n) :: fi, df2
!
!-------------------------------------------------------------------------------
!
! calculate the high-order interface interpolation; eq. (16) in [2]
!
    fi(1  ) = 5.0d-01 * (f(1  ) + f(2  ))
    do i = 2, n - 2
      fi(i) = (7.0d+00 * (f(i) + f(i+1)) - (f(i-1) + f(i+2))) / 1.2d+01
    end do
    fi(n-1) = 5.0d-01 * (f(n-1) + f(n  ))
    fi(n  ) = f(n)

! calculate second-order central derivative
!
    df2(1) = 0.0d+00
    do i = 2, n - 1
      df2(i) = (f(i+1) + f(i-1)) - 2.0d+00 * f(i)
    end do
    df2(n) = 0.0d+00

! limit the interpolation; eqs. (18) and (19) in [2]
!
    do i = 2, n - 2
      im1 = i - 1
      ip1 = i + 1

      if ((f(ip1) - fi(i)) * (fi(i) - f(i)) < 0.0d+00) then
        df2c = 3.0d+00 * ((f(ip1) + f(i  )) - 2.0d+00 * fi(i))
        df2m = df2(i  )
        df2p = df2(ip1)
        if (min(df2c, df2m, df2p) * max(df2c, df2m, df2p) > 0.0d+00) then
          df2l = sign(min(ppm_const * min(abs(df2m), abs(df2p)), abs(df2c)), df2c)
        else
          df2l = 0.0d+00
        end if
        fi(i) = 5.0d-01 * (f(i) + f(ip1)) - df2l / 6.0d+00
      end if
    end do

! iterate along the vector
!
    do i = 2, n - 1
      im1 = i - 1
      ip1 = i + 1

! limit states if local extremum is detected or the interpolation is not
! monotone
!
      alm = fi(im1) - f(i)
      alp = fi(i  ) - f(i)
      cm  = abs(alm) >= 2.0d+00 * abs(alp)
      cp  = abs(alp) >= 2.0d+00 * abs(alm)
      ext = .false.

! check if we have local extremum here
!
      if (alm * alp > 0.0d+00) then
        ext = .true.
      else if (cm .or. cp) then
        im2 = max(1, i - 1)

        dfm = fi(im1) - fi(im2)
        dfp = fi(ip1) - fi(i  )
        dcm = f (i  ) - f (im1)
        dcp = f (ip1) - f (i  )
        if (min(abs(dfm),abs(dfp)) >= min(abs(dcm),abs(dcp))) then
          ext = dfm * dfp < 0.0d+00
        else
          ext = dcm * dcp < 0.0d+00
        end if
      end if

! limit states if the local extremum is detected
!
      if (ext) then
        df2s = 6.0d+00 * (alm + alp)
        df2m = df2(im1)
        df2c = df2(i  )
        df2p = df2(ip1)
        if (min(df2s, df2c, df2m, df2p)                                        &
                      * max(df2s, df2c, df2m, df2p) > 0.0d+00) then
          df2l = sign(min(ppm_const * min(abs(df2m), abs(df2c), abs(df2p))     &
                                                   , abs(df2s)), df2s)
          if (abs(df2l) > 0.0d+00) then
            alm = alm * df2l / df2s
            alp = alp * df2l / df2s
          else
            alm = 0.0d+00
            alp = 0.0d+00
          end if
        else
          alm = 0.0d+00
          alp = 0.0d+00
        end if
      else

! the interpolation is not monotonic so apply additional limiter
!
        if (cp) then
          sgn = sign(1.0d+00, alp)
          amx = - 2.5d-01 * alp**2 / (alp + alm)
          dfp = f(ip1) - f(i)
          if (sgn * amx >= sgn * dfp) then
            del = dfp * (dfp - alm)
            if (del >= 0.0d+00) then
              alp = - 2.0d+00 * (dfp + sgn * sqrt(del))
            else
              alp = - 2.0d+00 * alm
            end if
          end if
        end if

        if (cm) then
          sgn = sign(1.0d+00, alm)
          amx = - 2.5d-01 * alm**2 / (alp + alm)
          dfm = f(im1) - f(i)
          if (sgn * amx >= sgn * dfm) then
            del = dfm * (dfm - alp)
            if (del >= 0.0d+00) then
              alm = - 2.0d+00 * (dfm + sgn * sqrt(del))
            else
              alm = - 2.0d+00 * alp
            end if
          end if
        end if
      end if

! update the states
!
      fr(im1) = f(i) + alm
      fl(i  ) = f(i) + alp

    end do ! i = 2, n

! update the interpolation of the first and last points
!
    fl(1  ) = fi(1)
    fr(n-1) = fi(n)
    fl(n  ) = fi(n)
    fr(n  ) = f (n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_ppm
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5Z:
! -----------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with
!   stencil weights by Borges et al. (2008).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Borges, R., Carmona, M., Costa, B., & Don, W.-S.,
!         "An improved weighted essentially non-oscillatory scheme for
!          hyperbolic conservation laws"
!         Journal of Computational Physics,
!         2008, vol. 227, pp. 3191-3211,
!         http://dx.doi.org/10.1016/j.jcp.2007.11.038
!
!===============================================================================
!
  subroutine reconstruct_weno5z(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt, df
    real(kind=8) :: al, ac, ar
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! iterate along the vector
!
    do i = 3, n - 2

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate βₖ (eqs. 9-11 in [1])
!
      bl  = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc  = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br  = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 25 in [1])
!
      tt  = abs(br - bl)

! calculate αₖ (eq. 28 in [1])
!
      al  = 1.0d+00 + tt / (bl + eps)
      ac  = 1.0d+00 + tt / (bc + eps)
      ar  = 1.0d+00 + tt / (br + eps)

! calculate weights
!
      wl  = dl * al
      wc  = dc * ac
      wr  = dr * ar
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! normalize weights
!
      wl  = dl * ar
      wc  = dc * ac
      wr  = dr * al
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 3, n - 2

! update the interpolation of the first and last two points
!
    fl(1)   = 0.5d+00 * (f(1) + f(2))
    df      = limiter_tvd(0.5d+00, dfm(2), dfp(2))
    fr(1)   = f(2) - df
    fl(2)   = f(2) + df
    i       = n - 1
    df      = limiter_tvd(0.5d+00, dfm(i), dfp(i))
    fr(i-1) = f(i) - df
    fl(i)   = f(i) + df
    fr(i)   = 0.5d+00 * (f(i) + f(n))
    fl(n)   = f(n)
    fr(n)   = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5z
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5YC:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with
!   stencil weights by Yamaleev & Carpenter (2009).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Yamaleev, N. K. & Carpenter, H. C.,
!         "A Systematic Methodology for Constructing High-Order Energy Stable
!          WENO Schemes"
!         Journal of Computational Physics,
!         2009, vol. 228, pp. 4248-4272,
!         http://dx.doi.org/10.1016/j.jcp.2009.03.002
!
!===============================================================================
!
  subroutine reconstruct_weno5yc(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br, tt, df
    real(kind=8) :: al, ac, ar
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! iterate along the vector
!
    do i = 3, n - 2

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate βₖ (eq. 19 in [1])
!
      bl  = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc  = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br  = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 20 in [1])
!
      tt  = (6.0d+00 * f(i) - 4.0d+00 * (f(im1) + f(ip1))                      &
                                                       + (f(im2) + f(ip2)))**2

! calculate αₖ (eqs. 18 or 58 in [1])
!
      al  = 1.0d+00 + tt / (bl + eps)
      ac  = 1.0d+00 + tt / (bc + eps)
      ar  = 1.0d+00 + tt / (br + eps)

! calculate weights
!
      wl  = dl * al
      wc  = dc * ac
      wr  = dr * ar
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! normalize weights
!
      wl  = dl * ar
      wc  = dc * ac
      wr  = dr * al
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 3, n - 2

! update the interpolation of the first and last two points
!
    fl(1)   = 0.5d+00 * (f(1) + f(2))
    df      = limiter_tvd(0.5d+00, dfm(2), dfp(2))
    fr(1)   = f(2) - df
    fl(2)   = f(2) + df
    i       = n - 1
    df      = limiter_tvd(0.5d+00, dfm(i), dfp(i))
    fr(i-1) = f(i) - df
    fl(i)   = f(i) + df
    fr(i)   = 0.5d+00 * (f(i) + f(n))
    fl(n)   = f(n)
    fr(n)   = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5yc
!
!===============================================================================
!
! subroutine RECONSTRUCT_WENO5NS:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Explicit Weighted Essentially Non-Oscillatory (WENO5) method with new
!   smoothness indicators and stencil weights by Ha et al. (2013).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ha, Y., Kim, C. H., Lee, Y. J., & Yoon, J.,
!         "An improved weighted essentially non-oscillatory scheme with a new
!          smoothness indicator",
!         Journal of Computational Physics,
!         2013, vol. 232, pp. 68-86
!         http://dx.doi.org/10.1016/j.jcp.2012.06.016
!
!===============================================================================
!
  subroutine reconstruct_weno5ns(n, h, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    real(kind=8) :: bl, bc, br
    real(kind=8) :: al, ac, ar, aa
    real(kind=8) :: wl, wc, wr
    real(kind=8) :: df, lq, l3, zt
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: dfm, dfp, df2

! weight coefficients
!
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00

! the free parameter for smoothness indicators (see Eq. 3.6 in [1])
!
    real(kind=8), parameter :: xi  =   4.0d-01
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = 0.5d+00 * abs(dfp(:) - dfm(:))

! iterate along the vector
!
    do i = 3, n - 2

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate βₖ (eq. 3.6 in [1])
!
      df  = abs(dfp(i))
      lq  = xi * df
      bl  = df2(im1) + xi * abs(2.0d+00 * dfm(i) - dfm(im1))
      bc  = df2(i  ) + lq
      br  = df2(ip1) + lq

! calculate ζ (below eq. 3.6 in [1])
!
      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ (eq. 3.9 in [4])
!
      al  = dl * (1.0d+00 + zt / (bl + eps)**2)
      ac  = dc * (1.0d+00 + zt / (bc + eps)**2)
      ar  = dr * (1.0d+00 + zt / (br + eps)**2)

! calculate weights
!
      aa  = (al + ar) + ac
      wl  = al / aa
      wr  = ar / aa
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i  ) = (wl * ql + wr * qr) + wc * qc

! calculate βₖ (eq. 3.6 in [1])
!
      df  = abs(dfm(i))
      lq  = xi * df
      bl  = df2(ip1) + xi * abs(2.0d+00 * dfp(i) - dfp(ip1))
      bc  = df2(i  ) + lq
      br  = df2(im1) + lq

! calculate ζ (below eq. 3.6 in [1])

      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ (eq. 3.9 in [4])
!
      al  = dl * (1.0d+00 + zt / (bl + eps)**2)
      ac  = dc * (1.0d+00 + zt / (bc + eps)**2)
      ar  = dr * (1.0d+00 + zt / (br + eps)**2)

! normalize weights
!
      aa  = (al + ar) + ac
      wl  = al / aa
      wr  = ar / aa
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(im1) = (wl * ql + wr * qr) + wc * qc

    end do ! i = 3, n - 2

! update the interpolation of the first and last two points
!
    fl(1)   = 0.5d+00 * (f(1) + f(2))
    df      = limiter_tvd(0.5d+00, dfm(2), dfp(2))
    fr(1)   = f(2) - df
    fl(2)   = f(2) + df
    i       = n - 1
    df      = limiter_tvd(0.5d+00, dfm(i), dfp(i))
    fr(i-1) = f(i) - df
    fl(i)   = f(i) + df
    fr(i)   = 0.5d+00 * (f(i) + f(n))
    fl(n)   = f(n)
    fr(n)   = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_weno5ns
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5Z:
! -------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method and smoothness indicators by Borges et al. (2008).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Borges, R., Carmona, M., Costa, B., & Don, W.-S.,
!         "An improved weighted essentially non-oscillatory scheme for
!          hyperbolic conservation laws"
!         Journal of Computational Physics,
!         2008, vol. 227, pp. 3191-3211,
!         http://dx.doi.org/10.1016/j.jcp.2007.11.038
!
!===============================================================================
!
  subroutine reconstruct_crweno5z(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    integer      :: iret
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n)   :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! prepare smoothness indicators
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1   = i - 1
      ip1   = i + 1

! calculate βₖ (eqs. 9-11 in [1])
!
      bl    = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc    = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br    = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 25 in [1])
!
      tt    = abs(br - bl)

! calculate αₖ (eq. 28 in [1])
!
      al(i) = 1.0d+00 + tt / (bl + eps)
      ac(i) = 1.0d+00 + tt / (bc + eps)
      ar(i) = 1.0d+00 + tt / (br + eps)

    end do ! i = 2, n - 1

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng + 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! calculate weights
!
      wl  = cl * al(i)
      wc  = cc * ac(i)
      wr  = cr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i)
      wc  = cc * ac(i)
      wr  = cr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 2, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 2, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n - 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n), iret)

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n), iret)

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5z
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5YC:
! --------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method and smoothness indicators by Yamaleev & Carpenter (2009).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Yamaleev, N. K. & Carpenter, H. C.,
!         "A Systematic Methodology for Constructing High-Order Energy Stable
!          WENO Schemes"
!         Journal of Computational Physics,
!         2009, vol. 228, pp. 4248-4272,
!         http://dx.doi.org/10.1016/j.jcp.2009.03.002
!
!===============================================================================
!
  subroutine reconstruct_crweno5yc(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    integer      :: iret
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n)   :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! smoothness indicator coefficients
!
    real(kind=8), parameter :: c1 = 1.3d+01 / 1.2d+01, c2 = 2.5d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! interpolation coefficients
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = c1 * (dfp(:) - dfm(:))**2

! prepare smoothness indicators
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im2   = max(1, i - 2)
      im1   = i - 1
      ip1   = i + 1
      ip2   = min(n, i + 2)

! calculate βₖ (eqs. 9-11 in [1])
!
      bl    = df2(im1) + c2 * (3.0d+00 * dfm(i  ) - dfm(im1))**2
      bc    = df2(i  ) + c2 * (          dfp(i  ) + dfm(i  ))**2
      br    = df2(ip1) + c2 * (3.0d+00 * dfp(i  ) - dfp(ip1))**2

! calculate τ (below eq. 64 in [3])
!
      tt  = (6.0d+00 * f(i) + (f(im2) + f(ip2))                                &
                                             - 4.0d+00 * (f(im1) + f(ip1)))**2

! calculate αₖ (eq. 28 in [1])
!
      al(i) = 1.0d+00 + tt / (bl + eps)
      ac(i) = 1.0d+00 + tt / (bc + eps)
      ar(i) = 1.0d+00 + tt / (br + eps)

    end do ! i = 2, n - 1

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng + 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! calculate weights
!
      wl  = cl * al(i)
      wc  = cc * ac(i)
      wr  = cr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i)
      wc  = cc * ac(i)
      wr  = cr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 2, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i)
      wc  = dc * ac(i)
      wr  = dr * ar(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 2, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i)
      wc  = dc * ac(i)
      wr  = dr * al(i)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n), iret)

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n), iret)

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5yc
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRWENO5NS:
! --------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact-Reconstruction Weighted Essentially Non-Oscillatory (CRWENO)
!   method combined with the smoothness indicators by Ha et al. (2013).
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Ghosh, D. & Baeder, J. D.,
!         "Compact Reconstruction Schemes with Weighted ENO Limiting for
!          Hyperbolic Conservation Laws"
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1706,
!         http://dx.doi.org/10.1137/110857659
!     [2] Ghosh, D. & Baeder, J. D.,
!         "Weighted Non-linear Compact Schemes for the Direct Numerical
!          Simulation of Compressible, Turbulent Flows"
!         Journal on Scientific Computing,
!         2014,
!         http://dx.doi.org/10.1007/s10915-014-9818-0
!     [3] Ha, Y., Kim, C. H., Lee, Y. J., & Yoon, J.,
!         "An improved weighted essentially non-oscillatory scheme with a new
!          smoothness indicator",
!         Journal of Computational Physics,
!         2013, vol. 232, pp. 68-86
!         http://dx.doi.org/10.1016/j.jcp.2012.06.016
!
!===============================================================================
!
  subroutine reconstruct_crweno5ns(n, h, f, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: f
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, im2, ip2
    integer      :: iret
    real(kind=8) :: bl, bc, br, tt
    real(kind=8) :: wl, wc, wr, ww
    real(kind=8) :: df, lq, l3, zt
    real(kind=8) :: ql, qc, qr

! local arrays for derivatives
!
    real(kind=8), dimension(n)   :: dfm, dfp, df2
    real(kind=8), dimension(n,2) :: al, ac, ar
    real(kind=8), dimension(n)   :: u
    real(kind=8), dimension(n,2) :: a, b, c, r

! the free parameter for smoothness indicators (see eq. 3.6 in [3])
!
    real(kind=8), parameter :: xi  =   4.0d-01

! weight coefficients for implicit (c) and explicit (d) interpolations
!
    real(kind=8), parameter :: cl = 1.0d+00 / 9.0d+00
    real(kind=8), parameter :: cc = 5.0d+00 / 9.0d+00
    real(kind=8), parameter :: cr = 1.0d+00 / 3.0d+00
    real(kind=8), parameter :: dl = 1.0d-01, dc = 6.0d-01, dr = 3.0d-01

! implicit method coefficients
!
    real(kind=8), parameter :: dq = 5.0d-01

! 3rd order interpolation coefficients for three stencils
!
    real(kind=8), parameter :: a11 =   2.0d+00 / 6.0d+00                       &
                             , a12 = - 7.0d+00 / 6.0d+00                       &
                             , a13 =   1.1d+01 / 6.0d+00
    real(kind=8), parameter :: a21 = - 1.0d+00 / 6.0d+00                       &
                             , a22 =   5.0d+00 / 6.0d+00                       &
                             , a23 =   2.0d+00 / 6.0d+00
    real(kind=8), parameter :: a31 =   2.0d+00 / 6.0d+00                       &
                             , a32 =   5.0d+00 / 6.0d+00                       &
                             , a33 = - 1.0d+00 / 6.0d+00
!
!-------------------------------------------------------------------------------
!
! calculate the left and right derivatives
!
    do i = 1, n - 1
      ip1      = i + 1
      dfp(i  ) = f(ip1) - f(i)
      dfm(ip1) = dfp(i)
    end do
    dfm(1) = dfp(1)
    dfp(n) = dfm(n)

! calculate the absolute value of the second derivative
!
    df2(:) = 0.5d+00 * abs(dfp(:) - dfm(:))

! prepare smoothness indicators
!
    do i = 2, n - 1

! prepare neighbour indices
!
      im1  = i - 1
      ip1  = i + 1

! calculate βₖ
!
      df  = abs(dfp(i))
      lq  = xi * df
      bl  = df2(im1) + xi * abs(2.0d+00 * dfm(i) - dfm(im1))
      bc  = df2(i  ) + lq
      br  = df2(ip1) + lq

! calculate ζ
!
      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ
!
      al(i,1) = 1.0d+00 + zt / (bl + eps)**2
      ac(i,1) = 1.0d+00 + zt / (bc + eps)**2
      ar(i,1) = 1.0d+00 + zt / (br + eps)**2

! calculate βₖ
!
      df  = abs(dfm(i))
      lq  = xi * df
      bl  = df2(im1) + lq
      bc  = df2(i  ) + lq
      br  = df2(ip1) + xi * abs(2.0d+00 * dfp(i) - dfp(ip1))

! calculate ζ

      l3  = df**3
      zt  = 0.5d+00 * ((bl - br)**2 + (l3 / (1.0d+00 + l3))**2)

! calculate αₖ
!
      al(i,2) = 1.0d+00 + zt / (bl + eps)**2
      ac(i,2) = 1.0d+00 + zt / (bc + eps)**2
      ar(i,2) = 1.0d+00 + zt / (br + eps)**2

    end do ! i = 2, n - 1

! prepare tridiagonal system coefficients
!
    do i = ng, n - ng + 1

! prepare neighbour indices
!
      im1 = i - 1
      ip1 = i + 1

! calculate weights
!
      wl  = cl * al(i,1)
      wc  = cc * ac(i,1)
      wr  = cr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,1) = 2.0d+00 * wl +            wc
      b(i,1) =           wl + 2.0d+00 * (wc + wr)
      c(i,1) =           wr

! prepare right hand side of tridiagonal equation
!
      r(i,1) = (wl * f(im1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(ip1)) * dq

! calculate weights
!
      wl  = cl * ar(i,2)
      wc  = cc * ac(i,2)
      wr  = cr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate tridiagonal matrix coefficients
!
      a(i,2) =           wr
      b(i,2) =           wl + 2.0d+00 * (wc + wr)
      c(i,2) = 2.0d+00 * wl +            wc

! prepare right hand side of tridiagonal equation
!
      r(i,2) = (wl * f(ip1) + (5.0d+00 * (wl + wc) + wr) * f(i  )              &
                                   + (wc + 5.0d+00 * wr) * f(im1)) * dq

    end do ! i = ng, n - ng + 1

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = 2, ng

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! calculate weights
!
      wl  = dl * al(i,1)
      wc  = dc * ac(i,1)
      wr  = dr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = 2, ng
    a(1,1) = 0.0d+00
    b(1,1) = 1.0d+00
    c(1,1) = 0.0d+00
    r(1,1) = 0.5d+00 * (f(1) + f(2))

! interpolate ghost zones using explicit solver (left-side reconstruction)
!
    do i = n - ng, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! calculate weights
!
      wl  = dl * al(i,1)
      wc  = dc * ac(i,1)
      wr  = dr * ar(i,1)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the left state
!
      ql = a11 * f(im2) + a12 * f(im1) + a13 * f(i  )
      qc = a21 * f(im1) + a22 * f(i  ) + a23 * f(ip1)
      qr = a31 * f(i  ) + a32 * f(ip1) + a33 * f(ip2)

! calculate the left state
!
      fl(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,1) = 0.0d+00
      b(i,1) = 1.0d+00
      c(i,1) = 0.0d+00
      r(i,1) = fl(i)

    end do ! i = n - ng, n - 1
    a(n,1) = 0.0d+00
    b(n,1) = 1.0d+00
    c(n,1) = 0.0d+00
    r(n,1) = f(n)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = 2, ng + 1

! prepare neighbour indices
!
      im2 = max(1, i - 2)
      im1 = i - 1
      ip1 = i + 1
      ip2 = i + 2

! normalize weights
!
      wl  = dl * ar(i,2)
      wc  = dc * ac(i,2)
      wr  = dr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = 2, ng + 1
    a(1,2) = 0.0d+00
    b(1,2) = 1.0d+00
    c(1,2) = 0.0d+00
    r(1,2) = f(1)

! interpolate ghost zones using explicit solver (right-side reconstruction)
!
    do i = n - ng + 1, n - 1

! prepare neighbour indices
!
      im2 = i - 2
      im1 = i - 1
      ip1 = i + 1
      ip2 = min(n, i + 2)

! normalize weights
!
      wl  = dl * ar(i,2)
      wc  = dc * ac(i,2)
      wr  = dr * al(i,2)
      ww  = (wl + wr) + wc
      wl  = wl / ww
      wr  = wr / ww
      wc  = 1.0d+00 - (wl + wr)

! calculate the interpolations of the right state
!
      ql = a11 * f(ip2) + a12 * f(ip1) + a13 * f(i  )
      qc = a21 * f(ip1) + a22 * f(i  ) + a23 * f(im1)
      qr = a31 * f(i  ) + a32 * f(im1) + a33 * f(im2)

! calculate the right state
!
      fr(i) = (wl * ql + wr * qr) + wc * qc

! prepare coefficients of the tridiagonal system
!
      a(i,2) = 0.0d+00
      b(i,2) = 1.0d+00
      c(i,2) = 0.0d+00
      r(i,2) = fr(i)

    end do ! i = n - ng + 1, n - 1
    a(n,2) = 0.0d+00
    b(n,2) = 1.0d+00
    c(n,2) = 0.0d+00
    r(n,2) = 0.5d+00 * (f(n-1) + f(n))

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,1), b(1:n,1), c(1:n,1), r(1:n,1), u(1:n), iret)

! substitute the left-side values
!
    fl(1:n  ) = u(1:n)

! solve the tridiagonal system of equations for the left-side interpolation
!
    call tridiag(n, a(1:n,2), b(1:n,2), c(1:n,2), r(1:n,2), u(1:n), iret)

! substitute the right-side values
!
    fr(1:n-1) = u(2:n)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (f(1) + f(2))
    fr(i) = 0.5d+00 * (f(i) + f(n))
    fl(n) = f(n)
    fr(n) = f(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crweno5ns
!
!===============================================================================
!
! subroutine RECONSTRUCT_MP5:
! --------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Monotonicity Preserving (MP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_mp5(n, h, fc, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.7d+01             &
                                                , 2.7d+01,-3.0d+00 /) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00 /) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
!! === left-side interpolation ===
!!
! reconstruct the interface state using the 5th order interpolation
!
    do i = 3, n - 2
      u(i) = sum(ce5(:) * fc(i-2:i+2))
    end do

! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fc(  1:  2))
    u(  2) = sum(ce3(:) * fc(  1:  3))
    u(n-1) = sum(ce3(:) * fc(n-2:  n))
    u(n  ) =              fc(n      )

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! reconstruct the interface state using the 5th order interpolation
!
    do i = 3, n - 2
      u(i) = sum(ce5(:) * fi(i-2:i+2))
    end do

! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fi(  1:  2))
    u(  2) = sum(ce3(:) * fi(  1:  3))
    u(n-1) = sum(ce3(:) * fi(n-2:  n))
    u(n  ) =              fi(n      )

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_mp5
!
!===============================================================================
!
! subroutine RECONSTRUCT_MP7:
! --------------------------
!
!   Subroutine reconstructs the interface states using the seventh order
!   Monotonicity Preserving (MP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_mp7(n, h, fc, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(7), parameter ::                                   &
                                ce7 = (/-3.0d+00, 2.5d+01,-1.01d+02, 3.19d+02  &
                                      , 2.14d+02,-3.8d+01, 4.0d+00 /) / 4.2d+02
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.7d+01             &
                                                , 2.7d+01,-3.0d+00 /) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00 /) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
!! === left-side interpolation ===
!!
! reconstruct the interface state using the 5th order interpolation
!
    do i = 4, n - 3
      u(i) = sum(ce7(:) * fc(i-3:i+3))
    end do

! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fc(  1:  2))
    u(  2) = sum(ce3(:) * fc(  1:  3))
    u(  3) = sum(ce5(:) * fc(  1:  5))
    u(n-2) = sum(ce5(:) * fc(n-4:  n))
    u(n-1) = sum(ce3(:) * fc(n-2:  n))
    u(n  ) =              fc(n      )

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! reconstruct the interface state using the 5th order interpolation
!
    do i = 4, n - 3
      u(i) = sum(ce7(:) * fi(i-3:i+3))
    end do

! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fi(  1:  2))
    u(  2) = sum(ce3(:) * fi(  1:  3))
    u(  3) = sum(ce5(:) * fi(  1:  5))
    u(n-2) = sum(ce5(:) * fi(n-4:  n))
    u(n-1) = sum(ce3(:) * fi(n-2:  n))
    u(n  ) =              fi(n      )

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_mp7
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRMP5:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Compact Reconstruction Monotonicity Preserving (CRMP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_crmp5(n, h, fc, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i, iret

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: a, b, c
    real(kind=8), dimension(n) :: r
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(3), parameter ::                                   &
                                di  = (/ 3.0d-01, 6.0d-01, 1.0d-01 /)
    real(kind=8), dimension(3), parameter ::                                   &
                                ci5 = (/ 1.0d+00, 1.9d+01, 1.0d+01 /) / 3.0d+01
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.7d+01             &
                                                , 2.7d+01,-3.0d+00 /) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00 /) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
! prepare the diagonals of the tridiagonal matrix
!
    do i = 1, ng
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do
    do i = ng + 1, n - ng - 1
      a(i) = di(1)
      b(i) = di(2)
      c(i) = di(3)
    end do
    do i = n - ng, n
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do

!! === left-side interpolation ===
!!
! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum(ci5(:) * fc(i-1:i+1))
    end do

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fc(  1:  2))
    r(  2) = sum(ce3(:) * fc(  1:  3))
    do i = 3, ng
      r(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    do i = n - ng, n - 2
      r(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    r(n-1) = sum(ce3(:) * fc(n-2:  n))
    r(n  ) =              fc(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum(ci5(:) * fi(i-1:i+1))
    end do ! i = ng, n - ng + 1

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fi(  1:  2))
    r(  2) = sum(ce3(:) * fi(  1:  3))
    do i = 3, ng
      r(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    do i = n - ng, n - 2
      r(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    r(n-1) = sum(ce3(:) * fi(n-2:  n))
    r(n  ) =              fi(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crmp5
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRMP5LD:
! ------------------------------
!
!   Subroutine reconstructs the interface states using the fifth order
!   Low-Dissipation Compact Reconstruction Monotonicity Preserving (CRMP)
!   method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!     [3] Ghosh, D. & Baeder, J.,
!         "Compact Reconstruction Schemes With Weighted ENO Limiting For
!          Hyperbolic Conservation Laws",
!         SIAM Journal on Scientific Computing,
!         2012, vol. 34, no. 3, pp. A1678-A1705,
!         http://dx.doi.org/10.1137/110857659
!
!===============================================================================
!
  subroutine reconstruct_crmp5ld(n, h, fc, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i, iret

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: a, b, c
    real(kind=8), dimension(n) :: r
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(3), parameter ::                                   &
                                di  = (/ 2.5d-01, 6.0d-01, 1.5d-01 /)
    real(kind=8), dimension(4), parameter ::                                   &
                                ci5 = (/ 3.0d+00, 6.7d+01, 4.9d+01             &
                                                         , 1.0d+00 /) / 1.2d+02
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.7d+01             &
                                                , 2.7d+01,-3.0d+00 /) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00 /) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
! prepare the diagonals of the tridiagonal matrix
!
    do i = 1, ng
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do
    do i = ng + 1, n - ng - 1
      a(i) = di(1)
      b(i) = di(2)
      c(i) = di(3)
    end do
    do i = n - ng, n
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do

!! === left-side interpolation ===
!!
! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum(ci5(:) * fc(i-1:i+2))
    end do

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fc(  1:  2))
    r(  2) = sum(ce3(:) * fc(  1:  3))
    do i = 3, ng
      r(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    do i = n - ng, n - 2
      r(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    r(n-1) = sum(ce3(:) * fc(n-2:  n))
    r(n  ) =              fc(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum(ci5(:) * fi(i-1:i+2))
    end do ! i = ng, n - ng + 1

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fi(  1:  2))
    r(  2) = sum(ce3(:) * fi(  1:  3))
    do i = 3, ng
      r(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    do i = n - ng, n - 2
      r(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    r(n-1) = sum(ce3(:) * fi(n-2:  n))
    r(n  ) =              fi(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crmp5ld
!
!===============================================================================
!
! subroutine RECONSTRUCT_CRMP7:
! ----------------------------
!
!   Subroutine reconstructs the interface states using the seventh order
!   Compact Reconstruction Monotonicity Preserving (CRMP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine reconstruct_crmp7(n, h, fc, fl, fr)

! include external procedures
!
    use algebra   , only : tridiag

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i, iret

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: a, b, c
    real(kind=8), dimension(n) :: r
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(3), parameter ::                                   &
                                di  = (/ 2.0d+00, 4.0d+00, 1.0d+00/) / 7.0d+00
    real(kind=8), dimension(5), parameter ::                                   &
                                ci  = (/-1.0d+00, 1.9d+01, 2.39d+02            &
                                              , 1.59d+02, 4.0d+00 /) / 4.2d+02
    real(kind=8), dimension(7), parameter ::                                   &
                                ce7 = (/-3.0d+00, 2.5d+01,-1.01d+02, 3.19d+02  &
                                      , 2.14d+02,-3.8d+01, 4.0d+00/) / 4.2d+02
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.70d+01            &
                                               , 2.70d+01,-3.0d+00/) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00/) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
! prepare the diagonals of the tridiagonal matrix
!
    do i = 1, ng
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do
    do i = ng + 1, n - ng - 1
      a(i) = di(1)
      b(i) = di(2)
      c(i) = di(3)
    end do
    do i = n - ng, n
      a(i) = 0.0d+00
      b(i) = 1.0d+00
      c(i) = 0.0d+00
    end do

!! === left-side interpolation ===
!!
! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum( ci(:) * fc(i-2:i+2))
    end do

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fc(  1:  2))
    r(  2) = sum(ce3(:) * fc(  1:  3))
    r(  3) = sum(ce5(:) * fc(  1:  5))
    do i = 4, ng
      r(i) = sum(ce7(:) * fc(i-3:i+3))
    end do
    do i = n - ng, n - 3
      r(i) = sum(ce7(:) * fc(i-3:i+3))
    end do
    r(n-2) = sum(ce5(:) * fc(n-4:  n))
    r(n-1) = sum(ce3(:) * fc(n-2:  n))
    r(n  ) =              fc(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! prepare the tridiagonal system coefficients for the interior part
!
    do i = ng, n - ng + 1
      r(i) = sum( ci(:) * fi(i-2:i+2))
    end do ! i = ng, n - ng + 1

! interpolate ghost zones using the explicit method
!
    r(  1) = sum(ce2(:) * fi(  1:  2))
    r(  2) = sum(ce3(:) * fi(  1:  3))
    r(  3) = sum(ce5(:) * fi(  1:  5))
    do i = 4, ng
      r(i) = sum(ce7(:) * fi(i-3:i+3))
    end do
    do i = n - ng, n - 3
      r(i) = sum(ce7(:) * fi(i-3:i+3))
    end do
    r(n-2) = sum(ce5(:) * fi(n-4:  n))
    r(n-1) = sum(ce3(:) * fi(n-2:  n))
    r(n  ) =              fi(n      )

! solve the tridiagonal system of equations
!
    call tridiag(n, a(1:n), b(1:n), c(1:n), r(1:n), u(1:n), iret)

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_crmp7
!
!===============================================================================
!
! subroutine PREPARE_GP:
! ---------------------
!
!   Subroutine prepares matrixes for the Gaussian Process (GP) method.
!
!===============================================================================
!
  subroutine prepare_gp(iret)

! include external procedures
!
    use algebra        , only : max_real_kind, invert
    use constants      , only : pi
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret

! local variables
!
    logical                  :: flag
    integer                  :: i, j, ip, jp
    real(kind=max_real_kind) :: sig, zl, zr, fc

! local arrays for derivatives
!
    real(kind=max_real_kind), dimension(:,:), allocatable :: cov, agp
    real(kind=max_real_kind), dimension(:)  , allocatable :: xgp

! local parameters
!
    character(len=*), parameter :: loc = 'INTERPOLATIONS::prepare_gp()'
!
!-------------------------------------------------------------------------------
!
! calculate normal distribution sigma
!
    sig = sqrt(2.0d+00) * sgp

! allocate the convariance matrix and interpolation position vector
!
    allocate(cov(ngp,ngp))
    allocate(agp(ngp,ngp))
    allocate(xgp(ngp))

! prepare the covariance matrix
!
    fc = 0.5d+00 * sqrt(pi) * sig
    i  = 0
    do ip = - mgp, mgp
      i = i + 1
      j = 0
      do jp = - mgp, mgp
        j = j + 1
        zl = (1.0d+00 * (ip - jp) - 0.5d+00) / sig
        zr = (1.0d+00 * (ip - jp) + 0.5d+00) / sig
        cov(i,j) = fc * (erf(zr) - erf(zl))
      end do
    end do

! invert the matrix
!
    call invert(ngp, cov(:,:), agp(:,:), flag)

! continue if the matrix inversion was successful
!
    if (flag) then

! make the inversed matrix symmetric
!
      do j = 1, ngp
        do i = 1, ngp
          if (i /= j) then
            fc       = 0.5d+00 * (agp(i,j) + agp(j,i))
            agp(i,j) = fc
            agp(j,i) = fc
          end if
        end do
      end do

! prepare the interpolation position vector
!
      j = 0
      do jp = - mgp, mgp
        j      = j + 1
        zr     = (0.5d+00 - 1.0d+00 * jp) / sig
        xgp(j) = exp(- zr**2)
      end do

! prepare the interpolation coefficients vector
!
      cgp(1:ngp) = matmul(xgp(1:ngp), agp(1:ngp,1:ngp))

    end if

! deallocate the convariance matrix and interpolation position vector
!
    deallocate(cov)
    deallocate(agp)
    deallocate(xgp)

! check if the matrix was inverted successfully
!
    if (.not. flag) then
      write(error_unit,"('[',a,']: ',a)") trim(loc)                            &
                      , "Could not invert the covariance matrix!"
      iret = 202
    end if

!-------------------------------------------------------------------------------
!
  end subroutine prepare_gp
!
!===============================================================================
!
! subroutine RECONSTRUCT_GP:
! -------------------------
!
!   Subroutine reconstructs the interface states using one dimensional
!   high order Gaussian Process (GP) method.
!
!   Arguments are described in subroutine reconstruct().
!
!===============================================================================
!
  subroutine reconstruct_gp(n, h, fc, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)  :: n
    real(kind=8)              , intent(in)  :: h
    real(kind=8), dimension(n), intent(in)  :: fc
    real(kind=8), dimension(n), intent(out) :: fl, fr

! local variables
!
    integer :: i

! local arrays for derivatives
!
    real(kind=8), dimension(n) :: fi
    real(kind=8), dimension(n) :: u

! local parameters
!
    real(kind=8), dimension(5), parameter ::                                   &
                                ce5 = (/ 2.0d+00,-1.3d+01, 4.70d+01            &
                                               , 2.70d+01,-3.0d+00/) / 6.0d+01
    real(kind=8), dimension(3), parameter :: &
                                ce3 = (/-1.0d+00, 5.0d+00, 2.0d+00/) / 6.0d+00
    real(kind=8), dimension(2), parameter :: ce2 = (/ 5.0d-01, 5.0d-01 /)
!
!-------------------------------------------------------------------------------
!
!! === left-side interpolation ===
!!
! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fc(  1:  2))
    u(  2) = sum(ce3(:) * fc(  1:  3))
    do i = 3, mgp
      u(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    do i = n + 1 - mgp, n - 2
      u(i) = sum(ce5(:) * fc(i-2:i+2))
    end do
    u(n-1) = sum(ce3(:) * fc(n-2:  n))
    u(n  ) =              fc(n      )

! reconstruct the interface state using the Gauss process interpolation
!
    do i = 1 + mgp, n - mgp
      u(i) = sum(cgp(1:ngp) * (fc(i-mgp:i+mgp) - fc(i))) + fc(i)
    end do

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fc(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fl(1:n) = u(1:n)

!! === right-side interpolation ===
!!
! invert the cell-centered value vector
!
    fi(1:n) = fc(n:1:-1)

! interpolate the interface state of the ghost zones using the interpolations
! of lower orders
!
    u(  1) = sum(ce2(:) * fi(  1:  2))
    u(  2) = sum(ce3(:) * fi(  1:  3))
    do i = 3, mgp
      u(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    do i = n + 1 - mgp, n - 2
      u(i) = sum(ce5(:) * fi(i-2:i+2))
    end do
    u(n-1) = sum(ce3(:) * fi(n-2:  n))
    u(n  ) =              fi(n      )

! reconstruct the interface state using the Gauss process interpolation
!
    do i = 1 + mgp, n - mgp
      u(i) = sum(cgp(1:ngp) * (fi(i-mgp:i+mgp) - fi(i))) + fi(i)
    end do

! apply the monotonicity preserving limiting
!
    call mp_limiting(n, fi(1:n), u(1:n))

! copy the interpolation to the respective vector
!
    fr(1:n-1) = u(n-1:1:-1)

! update the interpolation of the first and last points
!
    i     = n - 1
    fl(1) = 0.5d+00 * (fc(1) + fc(2))
    fr(i) = 0.5d+00 * (fc(i) + fc(n))
    fl(n) = fc(n)
    fr(n) = fc(n)

!-------------------------------------------------------------------------------
!
  end subroutine reconstruct_gp
!
!===============================================================================
!
! function LIMITER_ZERO:
! ---------------------
!
!   Function returns zero.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_zero(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.0d+00

!-------------------------------------------------------------------------------
!
  end function limiter_zero
!
!===============================================================================
!
! function LIMITER_MINMOD:
! -----------------------
!
!   Function returns the minimum module value among two arguments using
!   minmod limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_minmod(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.5d+00 * (sign(x, a) + sign(x, b)) * min(abs(a), abs(b))

!-------------------------------------------------------------------------------
!
  end function limiter_minmod
!
!===============================================================================
!
! function LIMITER_MONOTONIZED_CENTRAL:
! ------------------------------------
!
!   Function returns the minimum module value among two arguments using
!   the monotonized central TVD limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_monotonized_central(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = (sign(x, a) + sign(x, b)) * min(abs(a), abs(b), 2.5d-01 * abs(a + b))

!-------------------------------------------------------------------------------
!
  end function limiter_monotonized_central
!
!===============================================================================
!
! function LIMITER_SUPERBEE:
! -------------------------
!
!   Function returns the minimum module value among two arguments using
!   superbee limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_superbee(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = 0.5d+00 * (sign(x, a) + sign(x, b))                                    &
           * max(min(2.0d+00 * abs(a), abs(b)), min(abs(a), 2.0d+00 * abs(b)))

!-------------------------------------------------------------------------------
!
  end function limiter_superbee
!
!===============================================================================
!
! function LIMITER_VANLEER:
! ------------------------
!
!   Function returns the minimum module value among two arguments using
!   van Leer's limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_vanleer(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = a * b
    if (c > 0.0d+00) then
      c = 2.0d+00 * x * c / (a + b)
    else
      c = 0.0d+00
    end if

!-------------------------------------------------------------------------------
!
  end function limiter_vanleer
!
!===============================================================================
!
! function LIMITER_VANALBADA:
! --------------------------
!
!   Function returns the minimum module value among two arguments using
!   van Albada's limiter.
!
!   Arguments:
!
!     x    - scaling factor;
!     a, b - the input values;
!
!===============================================================================
!
  function limiter_vanalbada(x, a, b) result(c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: x, a, b
    real(kind=8)             :: c
!
!-------------------------------------------------------------------------------
!
    c = x * a * b * (a + b) / max(eps, a * a + b * b)

!-------------------------------------------------------------------------------
!
  end function limiter_vanalbada
!
!===============================================================================
!
! function MINMOD:
! ===============
!
!   Function returns the minimum module value among two arguments.
!
!   Arguments:
!
!     a, b - the input values;
!
!===============================================================================
!
  real(kind=8) function minmod(a, b)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b
!
!-------------------------------------------------------------------------------
!
    minmod = (sign(0.5d+00, a) + sign(0.5d+00, b)) * min(abs(a), abs(b))

    return

!-------------------------------------------------------------------------------
!
  end function minmod
!
!===============================================================================
!
! function MINMOD4:
! ================
!
!   Function returns the minimum module value among four arguments.
!
!   Arguments:
!
!     a, b, c, d - the input values;
!
!===============================================================================
!
  real(kind=8) function minmod4(a, b, c, d)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b, c, d
!
!-------------------------------------------------------------------------------
!
    minmod4 = minmod(minmod(a, b), minmod(c, d))

    return

!-------------------------------------------------------------------------------
!
  end function minmod4
!
!===============================================================================
!
! function MEDIAN:
! ===============
!
!   Function returns the median of three argument values.
!
!   Arguments:
!
!     a, b, c - the input values;
!
!===============================================================================
!
  real(kind=8) function median(a, b, c)

! local variables are not implicit by default
!
    implicit none

! input arguments
!
    real(kind=8), intent(in) :: a, b, c
!
!-------------------------------------------------------------------------------
!
    median = a + minmod(b - a, c - a)

    return

  end function median
!
!===============================================================================
!
! subroutine FIX_POSITIVITY:
! -------------------------
!
!   Subroutine scans the input arrays of the left and right states fl(:) and
!   fr(:) for negative values.  If it finds a negative value, it repeates the
!   state reconstruction from f(:) using the zeroth order interpolation.
!
!===============================================================================
!
  subroutine fix_positivity(n, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! input/output arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(in)    :: f
    real(kind=8), dimension(n), intent(inout) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1
    real(kind=8) :: fmn, fmx
!
!------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for positivity fix
!
    call start_timer(imf)
#endif /* PROFILE */

! check positivity only if desired
!
    if (positivity) then

! look for negative values in the states along the vector
!
      do i = 1, n

! check if the left state has a negative value
!
        if (fl(i) <= 0.0d+00) then

! calculate the left neighbour index
!
          im1 = max(1, i - 1)

! limit the states using the zeroth-order reconstruction
!
          fl(i  ) = f(i)
          fr(im1) = f(i)

        end if ! fl ≤ 0

! check if the right state has a negative value
!
        if (fr(i) <= 0.0d+00) then

! calculate the right neighbour index
!
          ip1 = min(n, i + 1)

! limit the states using the zeroth-order reconstruction
!
          fl(ip1) = f(ip1)
          fr(i  ) = f(ip1)

        end if ! fr ≤ 0

      end do ! i = 1, n

    end if ! positivity == .true.

#ifdef PROFILE
! stop accounting time for positivity fix
!
    call stop_timer(imf)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine fix_positivity
!
!===============================================================================
!
! subroutine CLIP_EXTREMA:
! -----------------------
!
!   Subroutine scans the reconstructed states and check if they didn't leave
!   the allowed limits.  In the case where the limits where exceeded,
!   the states are limited using constant reconstruction.
!
!   Arguments:
!
!     n      - the length of input vectors;
!     f      - the cell centered integrals of variable;
!     fl, fr - the left and right states of variable;
!
!===============================================================================
!
  subroutine clip_extrema(n, f, fl, fr)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(in)    :: f
    real(kind=8), dimension(n), intent(inout) :: fl, fr

! local variables
!
    integer      :: i, im1, ip1, ip2
    real(kind=8) :: fmn, fmx
    real(kind=8) :: dfl, dfr, df
!
!------------------------------------------------------------------------------
!
#ifdef PROFILE
! start accounting time for extrema clipping
!
    call start_timer(imc)
#endif /* PROFILE */

! iterate over all points
!
    do i = 1, n

! calculate indices
!
      im1 = max(1, i - 1)
      ip1 = min(n, i + 1)

! estimate the bounds of the allowed interval for reconstructed states
!
      fmn = min(f(i), f(ip1))
      fmx = max(f(i), f(ip1))

! check if the left state lays in the allowed range
!
      if (fl(i) < fmn .or. fl(i) > fmx) then

! calculate the left and right derivatives
!
        dfl = f(i  ) - f(im1)
        dfr = f(ip1) - f(i  )

! get the limited slope
!
        df  = limiter_clip(0.5d+00, dfl, dfr)

! calculate new states
!
        fl(i  ) = f(i  ) + df
        fr(im1) = f(i  ) - df

      end if

! check if the right state lays in the allowed range
!
      if (fr(i) < fmn .or. fr(i) > fmx) then

! calculate the missing index
!
        ip2 = min(n, i + 2)

! calculate the left and right derivatives
!
        dfl = f(ip1) - f(i  )
        dfr = f(ip2) - f(ip1)

! get the limited slope
!
        df  = limiter_clip(0.5d+00, dfl, dfr)

! calculate new states
!
        fl(ip1) = f(ip1) + df
        fr(i  ) = f(ip1) - df

      end if

    end do ! i = 1, n

#ifdef PROFILE
! stop accounting time for extrema clipping
!
    call stop_timer(imc)
#endif /* PROFILE */

!-------------------------------------------------------------------------------
!
  end subroutine clip_extrema
!
!===============================================================================
!
! subroutine MP_LIMITING:
! ----------------------
!
!   Subroutine applies the monotonicity preserving (MP) limiter to a vector of
!   high order reconstructed interface values.
!
!   Arguments:
!
!     n - the length of vectors;
!     f - the vector of cell-centered values;
!     v - the vector of interface values obtained from the high order
!         interpolation as input and its monotonicity limited values as output;
!
!   References:
!
!     [1] Suresh, A. & Huynh, H. T.,
!         "Accurate Monotonicity-Preserving Schemes with Runge-Kutta
!          Time Stepping"
!         Journal on Computational Physics,
!         1997, vol. 136, pp. 83-99,
!         http://dx.doi.org/10.1006/jcph.1997.5745
!     [2] He, ZhiWei, Li, XinLiang, Fu, DeXun, & Ma, YanWen,
!         "A 5th order monotonicity-preserving upwind compact difference
!          scheme",
!         Science China Physics, Mechanics and Astronomy,
!         Volume 54, Issue 3, pp. 511-522,
!         http://dx.doi.org/10.1007/s11433-010-4220-x
!
!===============================================================================
!
  subroutine mp_limiting(n, f, v)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer                   , intent(in)    :: n
    real(kind=8), dimension(n), intent(in)    :: f
    real(kind=8), dimension(n), intent(inout) :: v

! local variables
!
    integer      :: i, im1, ip1, ip2
    real(kind=8) :: df, ds, dc0, dc4, dm1, dp1, dml, dmr
    real(kind=8) :: flc, fmd, fmp, fmn, fmx, ful

! local vectors
!
    real(kind=8), dimension(0:n+2) :: dm
!
!-------------------------------------------------------------------------------
!
! calculate derivatives
!
    dm(0  ) = 0.0d+00
    dm(1  ) = 0.0d+00
    dm(2:n) = f(2:n) - f(1:n-1)
    dm(n+1) = 0.0d+00
    dm(n+2) = 0.0d+00

! check monotonicity condition for all elements and apply limiting if required
!
    do i = 1, n - 1

      ip1 = i + 1

      if (dm(i) * dm(ip1) >= 0.0d+00) then
        df = kappa * dm(i)
      else
        df = kbeta * dm(i)
      end if

      fmp  = f(i) + minmod(dm(ip1), df)
      ds   = (v(i) - f(i)) * (v(i) - fmp)

      if (ds > eps) then

        im1 = i - 1
        ip2 = i + 2

        dm1   = dm(i  ) - dm(im1)
        dc0   = dm(ip1) - dm(i  )
        dp1   = dm(ip2) - dm(ip1)
        dc4   = 4.0d+00 * dc0

        dml   = 0.5d+00 * minmod4(dc4 - dm1, 4.0d+00 * dm1 - dc0, dc0, dm1)
        dmr   = 0.5d+00 * minmod4(dc4 - dp1, 4.0d+00 * dp1 - dc0, dc0, dp1)

        fmd   = f(i) + 0.5d+00 * dm(ip1) - dmr
        ful   = f(i) +           df
        flc   = f(i) + 0.5d+00 * df      + dml

        fmx   = max(min(f(i), f(ip1), fmd), min(f(i), ful, flc))
        fmn   = min(max(f(i), f(ip1), fmd), max(f(i), ful, flc))

        v(i) = median(v(i), fmn, fmx)

      end if

    end do ! i = 1, n - 1

!-------------------------------------------------------------------------------
!
  end subroutine mp_limiting

!===============================================================================
!
end module interpolations
