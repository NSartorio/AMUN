!!******************************************************************************
!!
!! module: config - handling configuration file
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
module config

  implicit none

! dimensional variables
!
  integer(kind=4), save :: iblocks = 1, jblocks = 1, kblocks = 1

! block dimensions, number of ghost zones
!
  integer(kind=4), save :: ncells =  8, nghost =  4, ngrids = 16  &
                         , igrids = 16, jgrids = 16, kgrids = 16

! derived dimensional variables
!
  integer(kind=4), save :: in, im, ib, ie, ibl, ibu, iel, ieu, ih, ng  &
                         , jn, jm, jb, je, jbl, jbu, jel, jeu, jh, nd  &
                         , kn, km, kb, ke, kbl, kbu, kel, keu, kh, nh

! the dimensions of the lowest level
!
  integer(kind=4), dimension(3), save :: rdims = (/ 1, 1, 1 /)

! mesh refinement control
!
  integer(kind=4), save :: maxlev  = 2
  real           , save :: crefmin = 0.1
  real           , save :: crefmax = 0.2

! the maximum number of iterations
!
  integer(kind=4), save :: nmax   = 1

! domain bounds
!
  real           , save :: xmin = 0.0, xmax = 1.0  &
                         , ymin = 0.0, ymax = 1.0  &
                         , zmin = 0.0, zmax = 1.0

! the maximum time and initial time step
!
  real           , save :: tmax = 1.0, dtini = 1.0e-8

! data file type
!
  character      , save :: ftype = 'r'

! data storing time increment
!
  real           , save :: dtout = 0.0

! equation of state parameters (gamma, sound speed)
!
  real           , save :: gamma = 5./3.  , csnd = 1.0, csnd2 = 1.0  &
                         , gammam1 = 2./3., gammam1i = 1.5

! CFL time step condition
!
  real           , save :: cfl = 0.6

! problem selection
!
  character(len = 64), save :: problem = "blast"

! initial values
!
  real               , save :: dens = 1.0, pres = 1.0, rmid = 0.5, rcut = 0.1
  real               , save :: x1c = 0.0, y1c = 0.0, z1c = 0.0, r1c = 0.1
  real               , save :: x2c = 0.0, y2c = 0.0, z2c = 0.0, r2c = 0.1
  real               , save :: dnfac = 1.0, dnrat = 1.0
  real               , save :: v1ini = 0.0, v2ini = 2.0
  real               , save :: ratio = 1.0
  real               , save :: bamp = 1.0, bper = 0.0
  real               , save :: ydel = 0.1, ycut = 0.1

! boundary conditions
!
  character(len = 32), save :: xlbndry = "periodic"
  character(len = 32), save :: xubndry = "periodic"
  character(len = 32), save :: ylbndry = "periodic"
  character(len = 32), save :: yubndry = "periodic"
  character(len = 32), save :: zlbndry = "periodic"
  character(len = 32), save :: zubndry = "periodic"

  logical, dimension(NDIMS), save :: periodic = .true.

#if defined MHD && defined GLM
! coefficient controlling the decay of scalar potential Psi
!
  real               , save :: alpha_p = 0.5d0
#endif /* MHD & GLM */

  contains
!
!===============================================================================
!
! read_config: subroutine reads configuration file 'config.ini'
!
!===============================================================================
!
  subroutine read_config

    use error, only : print_error, print_warning
!
!-------------------------------------------------------------------------------
!
    implicit none

! local parameters
!
    character(len=*), parameter :: fl = "./config.in"

! local variables
!
    logical            :: info
    character(len=255) :: line, name, value
    integer            :: l
!
!-------------------------------------------------------------------------------
!
! check if file exists
!
    inquire(file = fl, exist = info)
    if (.not. info) then
      call print_error("config::read_config", "Could not find config.in file!")
    endif

! open file
!
    open(unit=1, file=fl, err=100)

10  read(1, fmt="(a)", end=20, err=200) line

! skip if line is empty or a comment (begins with '#')
!
    if (len_trim(line) .eq. 0) go to 10
    if (index(trim(adjustl(line)), '#') .eq. 1) go to 10

! parse line
!
    call parse_line(line, name, value)

! print error if the line is incorrect
!
    if (len_trim(name) .eq. 0 .or. len_trim(value) .eq. 0) go to 300

! substitute the parameter value
!
    select case(name)
    case("iblocks")
      read(value, "(i9.9)") iblocks
    case("jblocks")
      read(value, "(i9.9)") jblocks
#if NDIMS == 3
    case("kblocks")
      read(value, "(i9.9)") kblocks
#endif /* NDIMS == 3 */
    case("rdimsx")
      read(value, "(i9.9)") rdims(1)
    case("rdimsy")
      read(value, "(i9.9)") rdims(2)
#if NDIMS == 3
    case("rdimsz")
      read(value, "(i9.9)") rdims(3)
#endif /* NDIMS == 3 */
    case("ncells")
      read(value, "(i9.9)") ncells
    case("nghost")
      read(value, "(i9.9)") nghost
    case("maxlev")
      read(value, "(i9.9)") maxlev
    case("xmin")
      read(value,        *) xmin
    case("xmax")
      read(value,        *) xmax
    case("ymin")
      read(value,        *) ymin
    case("ymax")
      read(value,        *) ymax
    case("zmin")
      read(value,        *) zmin
    case("zmax")
      read(value,        *) zmax
    case("ftype")
      l = len_trim(value)
      write(ftype, "(a1)") value(2:l-1)
    case("nmax")
      read(value, "(i9.9)") nmax
    case("tmax")
      read(value,        *) tmax
    case("dtini")
      read(value,        *) dtini
    case("gamma")
      read(value,        *) gamma
    case("csnd")
      read(value,        *) csnd
    case("dtout")
      read(value,        *) dtout
    case("cfl")
      read(value,        *) cfl
    case("dens")
      read(value,        *) dens
    case("pres")
      read(value,        *) pres
    case("rmid")
      read(value,        *) rmid
    case("rcut")
      read(value,        *) rcut
    case("x1c")
      read(value,        *) x1c
    case("y1c")
      read(value,        *) y1c
    case("z1c")
      read(value,        *) z1c
    case("r1c")
      read(value,        *) r1c
    case("x2c")
      read(value,        *) x2c
    case("y2c")
      read(value,        *) y2c
    case("z2c")
      read(value,        *) z2c
    case("r2c")
      read(value,        *) r2c
    case("dnfac")
      read(value,        *) dnfac
    case("dnrat")
      read(value,        *) dnrat
    case("v1ini")
      read(value,        *) v1ini
    case("v2ini")
      read(value,        *) v2ini
    case("ratio")
      read(value,        *) ratio
    case("bamp")
      read(value,        *) bamp
    case("bper")
      read(value,        *) bper
    case("ydel")
      read(value,        *) ydel
    case("ycut")
      read(value,        *) ycut
    case("crefmin")
      read(value,        *) crefmin
    case("crefmax")
      read(value,        *) crefmax
    case ('problem')
      l = len_trim(value)
      write(problem  , "(a)") value(2:l-1)
    case ('xlbndry')
      l = len_trim(value)
      write(xlbndry  , "(a)") value(2:l-1)
    case ('xubndry')
      l = len_trim(value)
      write(xubndry  , "(a)") value(2:l-1)
    case ('ylbndry')
      l = len_trim(value)
      write(ylbndry  , "(a)") value(2:l-1)
    case ('yubndry')
      l = len_trim(value)
      write(yubndry  , "(a)") value(2:l-1)
    case ('zlbndry')
      l = len_trim(value)
      write(zlbndry  , "(a)") value(2:l-1)
    case ('zubndry')
      l = len_trim(value)
      write(zubndry  , "(a)") value(2:l-1)
#if defined MHD && defined GLM
    case("alpha_p")
      read(value,        *) alpha_p
#endif /* MHD & GLM */
    case default
      call print_warning("config::read_config", "Parameter '" // trim(name) // "' not implemented!")
    end select

! return to read the next line
!
    go to 10

! close if end of file
!
20  continue
    close(1)

! limit the minimum number of domain and ghost cells
#ifdef MHD
    nghost = max(nghost, 4)
#else /* MHD */
    nghost = max(nghost, 6)
#endif /* MHD */
    ncells = max(ncells, 2 * nghost)

! compute additional parameters
!
    ngrids = ncells + 2 * nghost
    igrids = ngrids
    jgrids = ngrids
#if NDIMS == 2
    kgrids = 1
#endif /* */
#if NDIMS == 3
    kgrids = ngrids
#endif /* */

    ng  = nghost
    nd  = ng * 2
    nh  = ng / 2

    in  = ncells
    im  = in + nd
    ih  = im / 2
    ib  = ng + 1
    ie  = ng + in
    ibl = ib - 1
    ibu = ib + ng - 1
    iel = ie - ng + 1
    ieu = ie + 1

    jn  = ncells
    jm  = jn + nd
    jh  = jm / 2
    jb  = ng + 1
    je  = ng + jn
    jbl = jb - 1
    jbu = jb + ng - 1
    jel = je - ng + 1
    jeu = je + 1

#if NDIMS == 2
    kn  = 1
    km  = 1
    kh  = 1
    kb  = 1
    ke  = 1
    kbl = 1
    kbu = 1
    kel = 1
    keu = 1
#endif /* NDIMS == 2 */
#if NDIMS == 3
    kn  = ncells
    km  = kn + nd
    kh  = km / 2
    kb  = ng + 1
    ke  = ng + kn
    kbl = kb - 1
    kbu = kb + ng - 1
    kel = ke - ng + 1
    keu = ke + 1
#endif /* NDIMS == 3 */

    gammam1  = gamma - 1.0
    gammam1i = 1.0 / gammam1
    csnd2    = csnd * csnd

    if (xlbndry .ne. 'periodic' .or. xubndry .ne. 'periodic')  &
      periodic(1) = .false.
    if (ylbndry .ne. 'periodic' .or. yubndry .ne. 'periodic')  &
      periodic(2) = .false.
#if NDIMS == 3
    if (zlbndry .ne. 'periodic' .or. zubndry .ne. 'periodic')  &
      periodic(3) = .false.
#endif /* NDIMS == 3 */

! return before error messages
!
    return

! print error if could not open file
!
100 continue
    call print_error("config::read_config", "Could not open config.in file!")

! print error if could not read next line
!
200 continue
    call print_error("config::read_config", "Could not read from config.in file!")

! print error if the line is incorrect
!
300 continue
    call print_error("config::read_config", "Improper line in config.in: " // trim(line))

!-------------------------------------------------------------------------------
!
  end subroutine read_config
!
!===============================================================================
!
! parse_line: subroutine parses line extracting the name and value of
!             parameter
!
!===============================================================================
!
  subroutine parse_line(line, name, value)
!
!-------------------------------------------------------------------------------
!
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: name, value

! local variables
!
    integer :: l, i
!
!-------------------------------------------------------------------------------
!
! get the length of line
!
    l = len_trim(line)

! find the index of '='
!
    i = index(line, '=')

! extract the name and value of parameter
!
    name  = trim( adjustl( line(1:i-1) ) )
    value = trim( adjustl( line(i+1:l) ) )

!-------------------------------------------------------------------------------
!
  end subroutine parse_line

!===============================================================================
!
end module
