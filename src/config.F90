!!*****************************************************************************
!!
!! module: config - handling configuration file
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
module config

  implicit none

! dimensional variables
!
  integer(kind=4), save :: iblocks = 1, jblocks = 1, kblocks = 1

! block dimensions, number of ghost zones
!
  integer(kind=4), save :: ncells = 8, nghost = 2, ngrids = 10

! mesh refinement control
!
  integer(kind=4), save :: maxlev = 2

! domain bounds
!
  real           , save :: xmin = 0.0, xmax = 1.0  &
                         , ymin = 0.0, ymax = 1.0  &
                         , zmin = 0.0, zmax = 1.0

  contains
!
!======================================================================
!
! read_config: subroutine reads configuration file 'config.ini'
!
!======================================================================
!
  subroutine read_config

    use error, only : print_error, print_warning
!
!----------------------------------------------------------------------
!
    implicit none

! local parameters
!
    character(len=*), parameter :: fl = "./config.in"

! local variables
!
    logical            :: info
    character(len=255) :: line, name, value
!
!----------------------------------------------------------------------
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
#ifdef R3D
    case("kblocks")
      read(value, "(i9.9)") kblocks
#endif /* R3D */
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

! compute additional parameters
!
    ngrids = ncells + 2 * nghost

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

!----------------------------------------------------------------------
!
  end subroutine read_config
!
!======================================================================
!
! parse_line: subroutine parses line extracting the name and value of
!             parameter
!
!======================================================================
!
  subroutine parse_line(line, name, value)
!
!----------------------------------------------------------------------
!
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: name, value

! local variables
!
    integer :: l, i
!
!----------------------------------------------------------------------
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

!----------------------------------------------------------------------
!
  end subroutine parse_line

!======================================================================
!
end module
