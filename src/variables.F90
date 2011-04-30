!!******************************************************************************
!!
!! module: variables - module handling the variable indices
!!
!! Copyright (C) 2010,2011 Grzegorz Kowal <grzegorz@gkowal.info>
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
module variables

  implicit none

! parameters
!
  integer(kind=4), parameter :: inx =  1, iny =  2, inz =  3
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
#ifdef GLM
#ifdef ISO
  integer(kind=4), parameter :: ibx =  5, iby =  6, ibz =  7
  integer(kind=4), parameter :: iph =  8
  integer(kind=4), parameter :: nvr =  8, nfl =  4, nqt =  8
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien =  5, ipr =  5
  integer(kind=4), parameter :: ibx =  6, iby =  7, ibz =  8
  integer(kind=4), parameter :: iph =  9
  integer(kind=4), parameter :: nvr =  9, nfl =  5, nqt =  9
#endif /* ADI */
#endif /* GLM */
#endif /* MHD */
!
!===============================================================================
!
end module variables
