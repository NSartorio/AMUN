!!*****************************************************************************
!!
!! module: boundaries - routines for handling the boundary conditions
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
module boundaries

  implicit none

  integer, save :: n
  real   , save :: t, dt, dtn

  contains
!
!===============================================================================
!
! boundary: subroutine sweeps over all leaf blocks and performs the boundary
!           update
!
!===============================================================================
!
  subroutine boundary

    use blocks  , only : block, plist, ndims, get_pointer
    use error   , only : print_error
    use mpitools, only : ncpu

    implicit none

! local variables
!
    integer :: i, j, k, dl

! local pointers
!
    type(block), pointer :: pblock, pneigh
!
!-------------------------------------------------------------------------------
!
! iterate over all blocks and perform boundary update
!
    pblock => plist
    do while (associated(pblock))

! if the current block is a leaf...
!
      if (pblock%leaf) then

! iterate over all neighbor blocks
!
        do i = 1, ndims
          do j = 1, 2
            do k = 1, 2

              if (pblock%neigh(i,j,k)%id .eq. -1) then

! neighbor is not associated, it means that we have non periodic boundary here
!
                if (k .eq. 1) &
                  call bnd_spec(pblock, i, j, k)

              else

! neighbor associated; exchange boundaries
!
                if (pblock%neigh(i,j,k)%cpu .eq. ncpu) then

! neighbor is on the same CPU, update
!
                  pneigh => get_pointer(pblock%neigh(i,j,k)%id)

! calculate the difference of current and neighbor levels
!
                  dl = pblock%level - pneigh%level

! depending on the level difference
!
                  select case(dl)
                  case(-1)  ! restriction and prolongation
                    call bnd_rest(pblock, pneigh, i, j, k)
                  case(0)   ! the same level, copying
                    if (k .eq. 1) &
                      call bnd_copy(pblock, pneigh, i, j, k)
                  case(1)   ! prolongation is handled by bnd_rest
                  case default
                    call print_error("boundaries::boundary", "Level difference unsupported!")
                  end select

                else

! neighbor is on another CPU
!
! TODO: use similar trick as in the particle exchange in Godunov, i.e.
!       1) construct an integer array info(ncpus,ncpus)
!       2) if the neighbors lay on different processes set info(source cpu, dest cpu) = 1
!       3) update info globally, so all processes know what to exchange
!       4) perform actual data exchange between processes

                endif

              endif

!               pneigh => get_pointer(pblock%neigh(i,j,k)%id)
!
! ! check if neighbor is associated
! !
!               if (associated(pneigh)) then
!
! ! neighbor is associated, which means the periodic boundary conditions
! ! or interior of the domain
! !
!                 if (pneigh%leaf) then
!
! ! calculate the difference of current and neighbor levels
! !
!                   dl = pblock%level - pneigh%level
!
! ! depending on the level difference
! !
!                   select case(dl)
!                   case(-1)  ! restriction and prolongation
!                     call bnd_rest(pblock, pneigh, i, j, k)
!                   case(0)   ! the same level, copying
!                     if (k .eq. 1) &
!                       call bnd_copy(pblock, pneigh, i, j, k)
!                   case(1)   ! prolongation is handled by bnd_rest
!                   case default
!                     call print_error("boundaries::boundary", "Level difference unsupported!")
!                   end select
!
! ! perform copying, prolongation or restriction
! !
!
!                 else
!                   print *, pneigh%id, 'is not a leaf'
!                 endif
!               else
!
!
!               endif

            end do
          end do
        end do

      endif

! assign pointer to the next block
!
      pblock => pblock%next

    end do

!-------------------------------------------------------------------------------
!
  end subroutine boundary
!
!===============================================================================
!
! bnd_copy: subroutine copies the interior of neighbor to update the boundaries
!           of current block
!
!===============================================================================
!
  subroutine bnd_copy(pb, pn, id, is, ip)

    use blocks, only : block, nv => nvars
    use config, only : im, ib, ibl, ibu, ie, iel, ieu                 &
                     , jm, jb, jbl, jbu, je, jel, jeu                 &
                     , km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning

    implicit none

! arguments
!
    type(block), pointer, intent(inout) :: pb, pn
    integer             , intent(in)    :: id, is, ip

! local variables
!
    integer :: ii, i, j, k
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * is

! perform update according to the flag
!
    select case(ii)
    case(110)
      do k = 1, km
        do j = 1, jm
          pb%u(1:nv,1:ibl,j,k) = pn%u(1:nv,iel:ie,j,k)
        end do
      end do
    case(120)
      do k = 1, km
        do j = 1, jm
          pb%u(1:nv,ieu:im,j,k) = pn%u(1:nv,ib:ibu,j,k)
        end do
      end do
    case(210)
      do k = 1, km
        do i = 1, im
          pb%u(1:nv,i,1:jbl,k) = pn%u(1:nv,i,jel:je,k)
        end do
      end do
    case(220)
      do k = 1, km
        do i = 1, im
          pb%u(1:nv,i,jeu:jm,k) = pn%u(1:nv,i,jb:jbu,k)
        end do
      end do
    case(310)
      do j = 1, jm
        do i = 1, im
          pb%u(1:nv,i,j,1:kbl) = pn%u(1:nv,i,j,kel:ke)
        end do
      end do
    case(320)
      do j = 1, jm
        do i = 1, im
          pb%u(1:nv,i,j,keu:km) = pn%u(1:nv,i,j,kb:kbu)
        end do
      end do
    case default
      call print_warning("boundaries::bnd_copy", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_copy
!
!===============================================================================
!
! bnd_rest: subroutine copies the interior of neighbor to update the boundaries
!           of current block
!
!===============================================================================
!
  subroutine bnd_rest(pb, pn, id, is, ip)

    use blocks, only : block, nv => nvars
    use config, only : ng, in, im, ib, ibl, ibu, ie, iel, ieu                 &
                         , jn, jm, jb, jbl, jbu, je, jel, jeu                 &
                         , kn, km, kb, kbl, kbu, ke, kel, keu
    use error , only : print_warning
    use interpolation, only : expand

    implicit none

! arguments
!
    type(block), pointer, intent(inout) :: pb, pn
    integer             , intent(in)    :: id, is, ip

! local variables
!
    integer :: ii, i, j, k, q, i0, i1, i2, j0, j1, j2, il, iu, jl, ju, dm(3), fm(3)
    real    :: dup, dum, du0, ds, du

! local arrays
!
    real, dimension(2*im,2*jm,km) :: ux

! local pointers
!
    type(block), pointer :: pn1, pn2
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * is + ip

! perform update according to the flag
!
    select case(ii)
    case(111)

! neighbor to current
!
      jl = ng / 2 + 1
      ju = jm / 2
      do k = 1, km
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = 1, ng
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = ng / 2 + 2
      dm(2) = jn / 2 + ng + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = ng
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + ng - 1
      j0 = ng / 2
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + jm - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,ieu:im,1:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(112)

! neighbor to current
!
      jl = jm / 2 + 1
      ju = jm - ng / 2
      do k = 1, km
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = 1, ng
            i2 = ie + 2 * (i - ng)
            i1 = i2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = ng / 2 + 2
      dm(2) = jn / 2 + ng + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = ng
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + ng - 1
      j0 = jm / 2 - ng / 2
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + jm - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,ieu:im,1:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(121)

! neighbor to current
!
      jl = ng / 2 + 1
      ju = jm / 2
      do k = 1, km
        do j = jl, ju
          j2 = 2 * j - ng
          j1 = j2 - 1
          do i = ieu, im
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = ng / 2 + 2
      dm(2) = jn / 2 + ng + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i1 = ie + 1
      i0 = i1 - dm(1) + 1
      il = 3
      iu = il + ng - 1
      j0 = ng / 2
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + jm - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:ng,1:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(122)

! neighbor to current
!
      jl = jm / 2 + 1
      ju = jm - ng / 2
      do k = 1, km
        do j = jl, ju
          j2 = 2 * j - jm + ng
          j1 = j2 - 1
          do i = ieu, im
            i2 = ng + 2 * (i - ie)
            i1 = i2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do


! current to neighbor
!
      dm(1) = ng / 2 + 2
      dm(2) = jn / 2 + ng + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i1 = ie + 1
      i0 = i1 - dm(1) + 1
      il = 3
      iu = il + ng - 1
      j0 = jm / 2 - ng / 2
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + jm - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:ng,1:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(211)

! neighbor to current
!
      il = ng / 2 + 1
      iu = im / 2
      do k = 1, km
        do i = il, iu
          i2 = 2 * i - ng
          i1 = i2 - 1
          do j = 1, ng
            j2 = je + 2 * (j - ng)
            j1 = j2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = in / 2 + ng + 2
      dm(2) = ng / 2 + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = ng / 2
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + im - 1
      j0 = ng
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + ng - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:im,jeu:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(212)

! neighbor to current
!
      il = im / 2 + 1
      iu = im - ng / 2
      do k = 1, km
        do i = il, iu
          i2 = 2 * i - im + ng
          i1 = i2 - 1
          do j = 1, ng
            j2 = je + 2 * (j - ng)
            j1 = j2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = in / 2 + ng + 2
      dm(2) = ng / 2 + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = im / 2 - ng / 2
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + im - 1
      j0 = ng
      j1 = j0 + dm(2) - 1
      jl = 3
      ju = jl + ng - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:im,jeu:jm,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(221)

! neighbor to current
!
      il = ng / 2 + 1
      iu = im / 2
      do k = 1, km
        do i = il, iu
          i2 = 2 * i - ng
          i1 = i2 - 1
          do j = jeu, jm
            j2 = ng + 2 * (j - je)
            j1 = j2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = in / 2 + ng + 2
      dm(2) = ng / 2 + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = ng / 2
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + im - 1
      j1 = je + 1
      j0 = j1 - dm(2) + 1
      jl = 3
      ju = jl + ng - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:im,1:ng,1:km) = ux(il:iu,jl:ju,1:km)

      end do

    case(222)

! neighbor to current
!
      il = im / 2 + 1
      iu = im - ng / 2
      do k = 1, km
        do i = il, iu
          i2 = 2 * i - im + ng
          i1 = i2 - 1
          do j = jeu, jm
            j2 = ng + 2 * (j - je)
            j1 = j2 - 1
            pb%u(1:nv,i,j,k) = 0.25 * (pn%u(1:nv,i1,j1,k) + pn%u(1:nv,i1,j2,k) &
                                     + pn%u(1:nv,i2,j1,k) + pn%u(1:nv,i2,j2,k))
          end do
        end do
      end do

! current to neighbor
!
      dm(1) = in / 2 + ng + 2
      dm(2) = ng / 2 + 2
      dm(3) = km
      fm(1) = 2 * dm(1)
      fm(2) = 2 * dm(2)
      fm(3) = km

      i0 = im / 2 - ng / 2
      i1 = i0 + dm(1) - 1
      il = 3
      iu = il + im - 1

      j1 = je + 1
      j0 = j1 - dm(2) + 1
      jl = 3
      ju = jl + ng - 1

! expand cube
!
      do q = 1, nv

        call expand(dm,fm,0,pb%u(q,i0:i1,j0:j1,1:km),ux,'t','t','t')

        pn%u(q,1:im,1:ng,1:km) = ux(il:iu,jl:ju,1:km)

      end do
    case default
      call print_warning("boundaries::bnd_rest", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_rest
!
!===============================================================================
!
! bnd_spec: subroutine to apply specific boundary conditions
!
!===============================================================================
!
  subroutine bnd_spec(pb, id, il, ip)

    use blocks, only : block, nv => nvars
    use config, only : xlbndry, xubndry, ylbndry, yubndry, zlbndry, zubndry    &
                     , ng, im, jm, km, ib, ibl, ie, ieu, jb, jbl, je, jeu      &
                     , kb, kbl, ke, keu
    use error , only : print_warning
    use interpolation, only : expand

    implicit none

! arguments
!
    type(block), pointer, intent(inout) :: pb
    integer             , intent(in)    :: id, il, ip

! local variables
!
    integer :: ii, i, j, k, it, jt, kt, is, js, ks, itm1, itp1, jtm1, jtp1, ktm1, ktp1
!
!-------------------------------------------------------------------------------
!
! calcuate the flag determinig the side of boundary to update
!
    ii = 100 * id + 10 * il

! perform update according to the flag
!
    select case(ii)
    case(110)
      select case(xlbndry)
      case("reflecting")
        do i = 1, ng
          it = ib  - i
          is = ibl + i
          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) =  pb%u(:,is,j,k)
              pb%u(2,it,j,k) = -pb%u(2,is,j,k)
          end do
          end do
        end do
      case default ! set "open" for default
        do i = 1, ng
          it   = ib - i
          itp1 = it + 1

          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) = pb%u(:,itp1,j,k)
            end do
          end do
        end do
      end select
    case(120)
      select case(xubndry)
      case("reflecting")
        do i = 1, ng
          it = ie  + i
          is = ieu - i
          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) =  pb%u(:,is,j,k)
              pb%u(2,it,j,k) = -pb%u(2,is,j,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do i = 1, ng
          it   = ie + i
          itm1 = it - 1

          do k = 1, km
            do j = 1, jm
              pb%u(:,it,j,k) = pb%u(:,itm1,j,k)
            end do
          end do
        end do
      end select
    case(210)
      select case(ylbndry)
      case("reflecting")
        do j = 1, ng
          jt = jb  - j
          js = jbl + j
          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) =  pb%u(:,i,js,k)
              pb%u(3,i,jt,k) = -pb%u(3,i,js,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do j = 1, ng
          jt   = jb - j
          jtp1 = jt + 1

          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) = pb%u(:,i,jtp1,k)
            end do
          end do
        end do
      end select
    case(220)
      select case(yubndry)
      case("reflecting")
        do j = 1, ng
          jt = je  + j
          js = jeu - j
          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) =  pb%u(:,i,js,k)
              pb%u(3,i,jt,k) = -pb%u(3,i,js,k)
            end do
          end do
        end do
      case default ! set "open" for default
        do j = 1, ng
          jt   = je + j
          jtm1 = jt - 1

          do k = 1, km
            do i = 1, im
              pb%u(:,i,jt,k) = pb%u(:,i,jtm1,k)
            end do
          end do
        end do
      end select
    case(310)
      select case(zlbndry)
      case("reflecting")
        do k = 1, ng
          kt = kb  - k
          ks = kbl + k
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) =  pb%u(:,i,j,ks)
              pb%u(4,i,j,kt) = -pb%u(4,i,j,ks)
            end do
          end do
        end do
      case default ! set "open" for default
        do k = 1, ng
          kt   = kb - k
          ktp1 = kt + 1

          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) = pb%u(:,i,j,ktp1)
            end do
          end do
        end do
      end select
    case(320)
      select case(zubndry)
      case("reflecting")
        do k = 1, ng
          kt = ke  + k
          ks = keu - k
          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) =  pb%u(:,i,j,ks)
              pb%u(4,i,j,kt) = -pb%u(4,i,j,ks)
            end do
          end do
        end do
      case default ! set "open" for default
        do k = 1, ng
          kt   = ke + k
          ktm1 = kt - 1

          do j = 1, jm
            do i = 1, im
              pb%u(:,i,j,kt) = pb%u(:,i,j,ktm1)
            end do
          end do
        end do
      end select
    case default
      call print_warning("boundaries::bnd_spec", "Boundary flag unsupported!")
    end select

!-------------------------------------------------------------------------------
!
  end subroutine bnd_spec

!===============================================================================
!
end module
