!!******************************************************************************
!!
!!  This file is part of the AMUN source code, a program to perform
!!  Newtonian or relativistic magnetohydrodynamical simulations on uniform or
!!  adaptive mesh.
!!
!!  Copyright (C) 2008-2013 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: PARAMETERS
!!
!!  This module handles runtime parameters by reading them from a parameter
!!  file and distributing them among all processors.
!!
!!******************************************************************************
!
module parameters

! module variables are not implicit by default
!
  implicit none

! module parameters determining the name and value field lengths, and the
! maximum length
!
  integer, parameter         :: nlen = 32, vlen = 128, mlen = 256

! the name of the parameter file
!
  character(len=mlen), save  :: fname   = './params.in'

! the number of parameters stored in the parameter file
!
  integer            , save  :: nparams = 0

! allocatable arrays to store parameter names and values
!
  character(len=nlen), dimension(:), allocatable, save :: pnames
  character(len=vlen), dimension(:), allocatable, save :: pvalues

! by default everything is private
!
  private

! declare public subroutines
!
  public :: read_parameters, finalize_parameters
  public :: get_parameter_integer, get_parameter_real, get_parameter_string
#ifdef MPI
  public :: redistribute_parameters
#endif /* MPI */

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine READ_PARAMETERS:
! --------------------------
!
!   Subroutine checks if the parameter file exists, checks its length,
!   allocates structures to store all parameters provided in the parameters
!   file and read parameters into these structures.
!
!   Arguments:
!
!     iret - the return value; if it is 0 everything went successfully,
!            otherwise there was a problem;
!
!   Note:
!
!     There is a possibility to specify customized input file by adding a
!     command line option -i or --input followed by the name of the input file.
!
!===============================================================================
!
  subroutine read_parameters(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret

! local variables
!
    character(len=mlen)    :: line
    integer                :: l
    logical                :: info

! external functions required to obtain comand line parameters
!
    integer                :: iargc
#ifdef GNU
    intrinsic              :: iargc, getarg
#else /* GNU */
    external               :: iargc, getarg
#endif /* GNU */
!
!-------------------------------------------------------------------------------
!
! parse the command line to check if another configuration file name has been
! specified
!
    do l = 1, iargc()
      call getarg(l, line)
      if (trim(line) .eq. '-i' .or. trim(line) .eq. '--input') then
        call getarg(l + 1, line)
        if (trim(line) .ne. '') then
          fname = trim(line)
        end if
      end if
    end do

! print information about the file from which the parameters will be read
!
    write(*,*) 'Reading parameters from ' // trim(fname)
    write(*,*)

! check if the file exists
!
    inquire(file = fname, exist = info)

! proceed if file exists
!
    if (info) then

! obtain the number of parameters stored in the file
!
      call get_parameters_number(iret)

! check if obtaining the number of parameters was successful
!
      if (iret .gt. 0) return

! if the parameter file is empty, print an error and return from the subroutine
!
      if (nparams .le. 0) then
        write(*,*) 'The parameter file ' // trim(fname)                        &
                                                    // ' is empty! Exiting...'
        write(*,*)
        iret = 110

        return
      end if

! allocate arrays to store the parameter names and values as string variables
!
      allocate(pnames (nparams))
      allocate(pvalues(nparams))

! get the parameter names and values and copy them to the corresponding arrays
!
      call get_parameters(iret)

    else

! the parameter file does not exists, so print a warning and exit
!
      write(*,*) 'The parameter file ' // trim(fname)                          &
                                              // ' does not exist! Exiting...'
      write(*,*)
      iret = 111

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_parameters
!
!===============================================================================
!
! subroutine FINALIZE_PARAMETERS:
! ------------------------------
!
!   Subroutine releases memory used by arrays in this module.
!
!===============================================================================
!
  subroutine finalize_parameters()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(pnames) ) deallocate(pnames)
    if (allocated(pvalues)) deallocate(pvalues)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_parameters
!
!===============================================================================
!
! subroutine GET_PARAMETERS_NUMBER:
! --------------------------------
!
!   Subroutine scans the input file and accounts the number of parameters
!   stored in it.
!
!   Arguments:
!
!     iret - the return value; if it is 0 everything went successfully,
!            otherwise there was a problem;
!
!===============================================================================
!
  subroutine get_parameters_number(iret)

! local variables are not implicit by default
!
    implicit none

! input and output variables
!
    integer, intent(inout) :: iret

! local variable to store the line content
!
    character(len=mlen)    :: line
!
!-------------------------------------------------------------------------------
!
! reset the number of parameters
!
    nparams = 0

! open the parameter file
!
    open(1, file = fname, err = 30)

! read the line
!
10  read(1, fmt = "(a)", end = 20) line

! if the line is empty or it's a comment, skip the counting
!
    if ((len_trim(line) .eq. 0)                                                &
                         .or. index(trim(adjustl(line)), '#') .eq. 1) go to 10

! increase the number of parameters
!
    nparams = nparams + 1

! go to the next line
!
    go to 10

! close the file
!
20  close(1)

! quit the subroutine without printing any errors since everything went fine
!
    return

! print a massage if an error occurred
!
30  write(*,*) 'There was a problem with reading the parameters file '         &
                                                                // trim(fname)
    write(*,*)

! set the return flag
!
    iret = 112

!-------------------------------------------------------------------------------
!
  end subroutine get_parameters_number
!
!===============================================================================
!
! subroutine GET_PARAMETERS:
! -------------------------
!
!   Subroutine scans the input file, reads parameter names and values, and
!   stores them in module arrays.
!
!   Arguments:
!
!     iret - the return value; if it is 0 everything went successfully,
!            otherwise there was a problem;
!
!===============================================================================
!
  subroutine get_parameters(iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret

! the parameter counter
!
    integer                :: np, nl

! local variables to store the line content, the parameter name and value
!
    character(len=256)     :: line, name, value
!
!-------------------------------------------------------------------------------
!
! initialize the parameter counter
!
    np = 1
    nl = 0

! open the parameter file
!
    open(1, file = fname, err = 30)

! read the line
!
10  read(1, fmt = "(a)", end = 20) line

! increase the line number
!
    nl = nl + 1

! if the line is empty or it's a comment, skip the counting
!
    if ((len_trim(line) .eq. 0)                                                &
                         .or. index(trim(adjustl(line)), '#') .eq. 1) go to 10

! parse the line to get parameter name and value
!
    call parse_line(line, name, value, iret)

! check if the line was parsed successfuly
!
    if (iret > 0) then
      write (*,"(1x,a,1x,a,'.')") "Wrong parameter format in"                  &
                                                        , trim(adjustl(fname))
      write (*,"(1x,'Line',i4,' : ',a)") nl, trim(line)
      write (*,*)
      go to 30
    end if

! fill the arrays of parameter names and values
!
    pnames (np) = name (1:nlen)
    pvalues(np) = value(1:vlen)

! increase the parameter counter
!
    np = np + 1

! go to the next line
!
    go to 10

! close the file
!
20  close(1)

! quit the subroutine without printing any errors since everything went fine
!
    return

! print a massage if an error occurred
!
30  write(*,*) 'There was a problem with reading the parameters file '         &
                                                                // trim(fname)
    write(*,*)

! set the return flag
!
    iret = 140

!-------------------------------------------------------------------------------
!
  end subroutine get_parameters
!
!===============================================================================
!
! subroutine PARSE_LINE:
! ---------------------
!
!   Subroutine extracts the parameter name and value from the input line.
!
!   Arguments:
!
!     line  - the input line containing the parameter information;
!     name  - the extracted name of the parameter;
!     value - the extracted value of the parameter;
!
!===============================================================================
!
  subroutine parse_line(line, name, value, iret)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: line
    character(len=*), intent(inout) :: name, value
    integer         , intent(out)   :: iret

! local indices to store positions in the input line
!
    integer :: l, p, c, i, j
!
!-------------------------------------------------------------------------------
!
! reset the return flag
!
    iret = 0

! get the length of line
!
    l = len_trim(line)

! find the indices of '=' and '#' in the line
!
    p = index(line, '=')
    c = index(line, '#')
    i = index(line, '"')
    j = index(line, '"', back = .true.)

! remove the length of the in-line comment from the length of line
!
    if (c .gt. 0) l = c - 1
    if (i .gt. 0 .and.  j .gt. 0) then
     i = i + 1
     j = j - 1
    else
     i = p + 1
     j = l
    end if

! limit the indices, so we don't overrun the variable memory
!
    p = min(p, nlen)
    j = min(j, vlen + i)

! extract the parameter name
!
    name  = trim(adjustl(line(1:p-1)))

! extract the parameter value
!
    value = trim(adjustl(line(i:j)))

! check possible errors in formatting
!
    if (p <= 2 .or. len_trim(name) == 0 .or. len_trim(value) == 0) iret = 1

!-------------------------------------------------------------------------------
!
  end subroutine parse_line
!
!===============================================================================
!
! subroutine GET_PARAMETER_INTEGER:
! --------------------------------
!
!   Subroutine reads a given parameter name and returns its integer value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output integer value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_integer(name, value)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    integer         , intent(inout) :: value

! local parameter counter
!
    integer                         :: np
!
!-------------------------------------------------------------------------------
!
! find the selected parameter
!
    np = 1
    do while (np .le. nparams)
      if (name .eq. pnames(np)) then
        read(pvalues(np), err = 100, fmt = *) value
      end if
      np = np + 1
    end do

    return

100 write(*,"(1x,'GET_PARAMETER_INTEGER:',1x,a)")                              &
            "Wrong format of the parameter '" // trim(name) // "'"
    write(*,"(24x,a)") "or the value is too large!"
    write(*,*)
    stop

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_integer
!
!===============================================================================
!
! subroutine GET_PARAMETER_REAL:
! -----------------------------
!
!   Subroutine reads a given parameter name and returns its real value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output real value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_real(name, value)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    real            , intent(inout) :: value

! local parameter counter
!
    integer                         :: np
!
!-------------------------------------------------------------------------------
!
! find the selected parameter
!
    np = 1
    do while (np .le. nparams)
      if (name .eq. pnames(np)) then
        read(pvalues(np), err = 100, fmt = *) value
      end if
      np = np + 1
    end do

    return

100 write(*,"(1x,'GET_PARAMETER_REAL   :',1x,a)")                              &
            "Wrong format of the parameter '" // trim(name) // "'"
    write(*,"(24x,a)") "or the value is too small or too large!"
    write(*,*)
    stop

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_real
!
!===============================================================================
!
! subroutine GET_PARAMETER_STRING:
! -------------------------------
!
!   Subroutine reads a given parameter name and returns its string value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output string value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_string(name, value)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    character(len=*), intent(inout) :: value

! local parameter counters
!
    integer                         :: np, nl
!
!-------------------------------------------------------------------------------
!
! get the length of the output string
!
    nl = min(vlen, len(value))

! find the selected parameter
!
    np = 1
    do while (np .le. nparams)
      if (name .eq. pnames(np)) then
        value = pvalues(np)(1:nl)
      end if
      np = np + 1
    end do

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_string
#ifdef MPI
!
!===============================================================================
!
! subroutine REDISTRIBUTE_PARAMETERS:
! ----------------------------------
!
!   Subroutine redistributes parameters among all processors.
!
!   Arguments:
!
!     iret - the return value; if it is 0 everything went successfully,
!            otherwise there was a problem;
!
!===============================================================================
!
  subroutine redistribute_parameters(iret)

! include external procedures and variables
!
    use mpitools, only : master, bcast_integer_variable, bcast_string_variable

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(inout) :: iret

! local parameter counter
!
    integer                :: np
!
!-------------------------------------------------------------------------------
!
! broadcast the number of parameters
!
    call bcast_integer_variable(nparams, iret)

! allocate the arrays to store parameter names and values on remaining
! processors
!
    if (.not. master) then
      allocate(pnames (nparams))
      allocate(pvalues(nparams))
    end if

! send parameter names and values from master to remaining processors
!
    do np = 1, nparams
      call bcast_string_variable(pnames (np), iret)
      call bcast_string_variable(pvalues(np), iret)
    end do

!-------------------------------------------------------------------------------
!
  end subroutine redistribute_parameters
#endif /* MPI */

!===============================================================================
!
end module parameters
