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
module CDMSreader__constants
  !! Contains special numbers and conversion factors
  use CDMSreader__types, only: dp

  implicit none

  private

  real(dp), parameter, public :: au2invcm  = 219474.6313710e0_dp ! -- atomic units -> inverse centimeters
  real(dp), parameter, public :: au2eV     = 27.2113834e0_dp     ! -- atomic units -> electron volts
  real(dp), parameter, public :: c_light = 299792458 ! -- speed of light in m/s
  real(dp), parameter, public :: invcm2Hz = c_light * 1e2 ! -- inverse centimeters -> Hz

! ================================================================================================================================ !
end module CDMSreader__constants
! ================================================================================================================================ !
