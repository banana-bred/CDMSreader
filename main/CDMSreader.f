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
program main
! ================================================================================================================================ !

  use iso_fortran_env,       only: input_unit, output_unit
  use CDMSreader__types,     only: asymtop_state_hfs, asymtop_state_nohfs, asymtop_transition_hfs, asymtop_transition_nohfs
  use CDMSreader__readwrite, only: CDMS_readfile

  type(asymtop_state_hfs),        allocatable :: states_hfs(:)
  type(asymtop_state_nohfs),      allocatable :: states_nohfs(:)
  type(asymtop_transition_hfs),   allocatable :: transitions_hfs(:)
  type(asymtop_transition_nohfs), allocatable :: transitions_nohfs(:)


  call CDMS_readfile( input_unit, output_unit, return_hfs = .true., return_nohfs = .true. &
                    , states_hfs = states_hfs, states_nohfs = states_nohfs                &
                    , transitions_hfs = transitions_hfs, transitions_nohfs = transitions_nohfs)

! ================================================================================================================================ !
end program main
! ================================================================================================================================ !
