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
!! module: variables - module handling the variable indices
!!
!!******************************************************************************
!
module variables

  implicit none

! parameters
!
  integer(kind=4), parameter :: inx =  1, iny =  2, inz =  3
  integer(kind=4), parameter :: idn =  1, imx =  2, imy =  3, imz =  4 &
                                        , ivx =  2, ivy =  3, ivz =  4
#ifdef ISO
  integer(kind=4), parameter :: nfl =  4
#endif /* ISO */
#ifdef ADI
  integer(kind=4), parameter :: ien =  5, ipr =  5
  integer(kind=4), parameter :: nfl =  5
#endif /* ADI */
#ifdef HYDRO
  integer(kind=4), parameter :: nqt = nfl
#endif /* HYDRO */
#ifdef MHD
  integer(kind=4), parameter :: ibx = nfl + 1, iby = ibx + 1, ibz = iby + 1
#ifdef GLM
  integer(kind=4), parameter :: iph =  ibz + 1
  integer(kind=4), parameter :: nqt = iph
#else /* GLM */
  integer(kind=4), parameter :: nqt = ibz
#endif /* GLM */
#endif /* MHD */
  integer(kind=4), parameter :: nvr = nqt
  integer(kind=4), parameter :: nt  = nqt
!
!===============================================================================
!
end module variables
