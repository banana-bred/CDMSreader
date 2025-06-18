! ================================================================================================================================ !
program main
! ================================================================================================================================ !

  use iso_fortran_env,       only: input_unit, output_unit
  use CDMSreader__types,     only: asymtop_state, asymtop_transition
  use CDMSreader__readwrite, only: CDMS_readfile

  class(asymtop_state),     allocatable :: states_hfs(:)
  class(asymtop_state),     allocatable :: states_nohfs(:)
  class(asymtop_transition), allocatable :: transitions_hfs(:)
  class(asymtop_transition), allocatable :: transitions_nohfs(:)


  call CDMS_readfile( input_unit, output_unit, return_hfs = .true., return_nohfs = .true. &
                    , states_hfs = states_hfs, states_nohfs = states_nohfs                &
                    , transitions_hfs = transitions_hfs, transitions_nohfs = transitions_nohfs)

! ================================================================================================================================ !
end program main
! ================================================================================================================================ !
