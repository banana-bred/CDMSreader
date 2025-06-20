! CDMSreader: reads transition data from the CDMS and determines the average lifetimes of the involved states
! Copyright (C) 2025 Josh Forer <j.forer@posteo.net>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
! ================================================================================================================================ !
module CDMSreader__system
! ================================================================================================================================ !

  implicit none

  private

  public :: die

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental subroutine die(message)
    !! Stop program execution with a message
    implicit none
    character(*), intent(in), optional :: message
    if(.not.present(message)) error stop ; error stop message
  end subroutine die

! ================================================================================================================================ !
end module CDMSreader__system
! ================================================================================================================================ !
