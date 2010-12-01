!!******************************************************************************
!!
!! module: variables - module handling the variable indices
!!
!! Copyright (C) 2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module variables

  implicit none

! parameters
!
  integer(kind=4), parameter :: idn =  1, imx =  2, imy =  3, imz =  4 &
                                        , ivx =  2, ivy =  3, ivz =  4
#ifdef HYDRO
#ifdef ISO
  integer(kind=4), parameter :: nvr =  4, nfl =  4, nqt =  4
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien =  5, ipr =  5
  integer(kind=4), parameter :: nvr =  5, nfl =  5, nqt =  5
#endif /* ADI */
#endif /* HYDRO */
#ifdef MHD
#ifdef FIELDCD
#ifdef ISO
  integer(kind=4), parameter :: ibx =  5, iby =  6, ibz =  7
  integer(kind=4), parameter :: icx =  5, icy =  6, icz =  7
  integer(kind=4), parameter :: nvr =  7, nfl =  4, nqt =  7
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien =  5, ipr =  5
  integer(kind=4), parameter :: ibx =  6, iby =  7, ibz =  8
  integer(kind=4), parameter :: icx =  6, icy =  7, icz =  8
  integer(kind=4), parameter :: nvr =  8, nfl =  5, nqt =  8
#endif /* ADI */
#endif /* FIELDCD */
#ifdef FLUXCT
#ifdef ISO
  integer(kind=4), parameter :: ibx =  5, iby =  6, ibz =  7
  integer(kind=4), parameter :: icx =  8, icy =  9, icz = 10
  integer(kind=4), parameter :: nvr = 10, nfl =  4, nqt =  7
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien =  5, ipr =  5
  integer(kind=4), parameter :: ibx =  6, iby =  7, ibz =  8
  integer(kind=4), parameter :: icx =  9, icy = 10, icz = 11
  integer(kind=4), parameter :: nvr = 11, nfl =  5, nqt =  8
#endif /* ADI */
#endif /* FLUXCT */
#endif /* MHD */
!
!===============================================================================
!
end module variables
