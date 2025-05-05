! ================================================================================================================================ !
program CDMSreader
! ================================================================================================================================ !

  use CDMSreader__types,             only : dp, asymtop_state_hfs, asymtop_transition, add_to &
                                          , sort_last_transition, find_state_number, sort_last_state, asymtop_state_nohfs &
                                          , asymtop_state, make_asymtop_state, sort_states
  use CDMSreader__readwrite,         only : CDMS_readline, write_states, write_transitions
  use CDMSreader__system,            only : die
  use CDMSreader__constants,         only : invcm2Hz
  use, intrinsic :: iso_fortran_env, only : stdin => input_unit, stdout => output_unit, iostat_end

  implicit none

  ! -- data in the cat file
  integer      :: DR
    !! Degrees of freedom in the rotational partition function (0 for atoms, 2 for linear molecules, and 3 for nonlinear molecules)
  integer      :: GUP
    !! Upper state degeneracy, gup = gI * gF, where gF = 2F+1
  integer      :: TAG
    !! Species tag or molecular identifier. A negative value flags that the line frequency has been measured in the laboratory.
    !! The absolute value of TAG is then the species tag (as given in line 2 of file.int above) and ERR is the reported experimental
    !! error.
  integer :: qnfmt
    !! Identifies the format of the quantum numbers given in the field QN.
  real(dp)     :: freq
    !! Frequency of the transition
  real(dp)     :: EinsteinA
    !! The Einstein A coefficient for a transition
  real(dp)     :: err
    !! Estimated or experimental error (999.9999 indicates error is larger)
  real(dp)     :: EinstA
    !! The Einstein A coefficient
  real(dp)     :: sigmaA2
    !! The square of the uncertainty in the Einstein A coefficient for a transition, given by the uncertainty in the trantision
    !! frequency
  real(dp)     :: elo, eup
    !! state energy in cm–1
  integer :: QN(12)
    !! Twice the quantum numbers numbers. These are integers but are converted in the reading routine from their coding format
    !! according to QNFMT. Upper state quanta start in character 1. Lower state quanta start in character 14 (element 7).
    !! Unused quanta are blank, quanta whose magnitude is larger than 99 or smaller than –9 are shown with alphabetic characters
    !! or **. Quanta between –10 and –19 are shown as a0 through a9. Similarly, –20 is b0, etc., up to –259, which is shown as z9.
    !! Quanta between 100 and 109 are shown as A0 through A9. Similarly, 110 is B0, etc., up to 359, which is shown as Z9.
  character(:), allocatable :: mol

  integer :: degenup, degenlo, dNup, dNlo, dKaup, dKalo, dKcup, dKclo, dJup, dJlo, dItotup, dItotlo

  integer                   :: io
  logical                   :: skip
  character(:), allocatable :: fmt

  real(dp) :: dF

  type(asymtop_state_hfs) :: up_hfs
  type(asymtop_state_hfs) :: lo_hfs
  ! type(asymtop_state_hfs),  allocatable :: states_hfs(:)
  class(asymtop_state),  allocatable :: states_hfs(:)

  type(asymtop_state_nohfs) :: up_nohfs
  type(asymtop_state_nohfs) :: lo_nohfs
  ! type(asymtop_state_nohfs),  allocatable :: states_nohfs(:)
  class(asymtop_state),  allocatable :: states_nohfs(:)

  type(asymtop_transition) :: transitionul
  type(asymtop_transition), allocatable :: transitions(:)

  integer :: i
  integer :: Q
    !! the number in square brackets in the table
  integer :: H
    !! binary code to indicate which of the last three quantum numbers are half integer quanta (1 indicates that F is half integer)
    !! The least significant b it of H refers to the F quantum number and is 1 if F is half integer.
  character(3) :: Hchar
  integer :: Hbits(3)
    !! Binary code representing which of the last three quantum numbesr are integral (0) or half integral (1)
  integer :: NQN
    !! number of quanta per state
  integer :: R

  integer :: iup, ilo, iup_nohfs, ilo_nohfs

  do

    call CDMS_readline(stdin, freq, err, EinstA, dr, elo, gup, tag, qnfmt, QN(1:12), mol, io, skip, "#")

    if(io   .eq.  iostat_end) exit
    if(skip .eqv. .true.)     cycle

    ! -- we actually read log10(A)
    EinstA = 10**EinstA

    ! -- make sure freq is in inverse centimeters. If the uncertainty is < 0, then this is the case. Otherwise,
    !    it's in MHz. Energy seems to always be in inverse centimeters. The Einstein coefficients need to be converted
    !    to 1/s if the frequencies are given in inverse centimeters
    if( err .gt. 0 ) then
      err    = -err * 1e6 / invcm2Hz
      freq   = freq * 1e6 / invcm2Hz
    else
      EinstA = EinstA / (invcm2hz * 1e-6)
      err = -err
    endif

    ! -- determine the uncertainty in A
    sigmaA2 = (3*err/freq)**2

    ! -- decrypt the qnfmt message
    Q   = qnfmt / 100
    R   = mod(qnfmt, 100)
    H   = R / 10
    NQN = mod(R, 10)

    if(Q*100 + H*10 + NQN .ne. qnfmt) call die("Problem reading qnfmt")

    ! -- decrypt H
    write(Hchar, "(b3.3)") H
    read(Hchar, *) H
    R = mod(H, 10)
    Hbits = [H/100, R/10, mod(R, 10)]

    if(any(Hbits .lt. 0 .OR. Hbits .gt. 1)) call die("Somehow extracted a bit that is neither 0 nor 1")

    if(Q .ne. 23) call die("Q =/= 23 detected. No other cases have been programmed yet")

    eup = freq + elo

    dNup    = 2*QN(1)
    dKaup   = 2*QN(2)
    dKcup   = 2*QN(3)
    dNlo    = 2*QN(7)
    dKalo   = 2*QN(8)
    dKclo   = 2*QN(9)
    dJup    = 2*QN(4)   - Hbits(1)
    dJlo    = 2*QN(10)  - Hbits(1)
    dItotup = 2*QN(5)   - Hbits(2)
    dItotlo = 2*QN(11)  - Hbits(2)
    degenup = (2*QN(6)  - Hbits(3)) + 1
    degenlo = (2*QN(12) - Hbits(3)) + 1

    ! -- create the upper and lower state
    up_nohfs = make_asymtop_state( dN = dNup, dKa = dKaup, dKc = dKcup, E = 0.0_dp, EinstA = 0.0_dp, sigmaA2 = 0.0_dp, degen = 0 )
    lo_nohfs = make_asymtop_state( dN = dNlo, dKa = dKalo, dKc = dKclo, E = 0.0_dp, EinstA = 0.0_dp, sigmaA2 = 0.0_dp, degen = 0 )

    ! -- add lower non hfs state to states array if it's not already there
    !    and keep track of the lower state's index
    call find_state_number(lo_nohfs, states_nohfs, ilo_nohfs)
    if(ilo_nohfs .eq. 0) then
      call add_to(lo_nohfs, states_nohfs)
      ilo_nohfs = size(states_nohfs, 1)
    endif

    ! -- add upper non hfs state to states array if it's not already there
    !    and keep track of the upper state's index
    call find_state_number(up_nohfs, states_nohfs, iup_nohfs)
    if(iup_nohfs .eq. 0) then
      call add_to(up_nohfs, states_nohfs)
      iup_nohfs = size(states_nohfs, 1)
    endif

    ! -- create the state with hyperfine splitting
    up_hfs = make_asymtop_state( up_nohfs % dN, up_nohfs % dKa, up_nohfs % dKc, dJ = dJup, dItot = dItotup, dF = degenup - 1 &
                               , E = eup, EinstA = 0.0_dp, sigmaA2 = 0.0_dp &
         )
    lo_hfs = make_asymtop_state( lo_nohfs % dN, lo_nohfs % dKa, lo_nohfs % dKc, dJ = dJlo, dItot = dItotlo, dF = degenlo - 1 &
                               , E = elo, EinstA = 0.0_dp, sigmaA2 = 0.0_dp &
         )
    transitionul = asymtop_transition( up = up_hfs, lo = lo_hfs, freq = freq, EinstA = EinstA, err = err &
                                     , dr = dr, gup = gup                                                &
         )

    ! -- add lower hfs state to states array if it's not already there
    !    and keep track of the lower state's index. A new hfs state means
    !    that we have to add its energy and degeneracy to the non hfs state
    call find_state_number(lo_hfs, states_hfs, ilo)
    if(ilo .eq. 0) then
      call add_to(lo_hfs, states_hfs)
      call sort_last_state(states_hfs)
      ! -- new hfs state, add to corresponding non hfs state
      select type(state_nohfs => states_nohfs(ilo_nohfs))
      type is (asymtop_state_nohfs)
        state_nohfs % E     = state_nohfs % E     + Elo * degenlo
        state_nohfs % degen = state_nohfs % degen + degenlo
      class default
        call die("Wrong type in lower asymtop non hfs state assignment")
      end select
    endif

    ! -- add upper hfs state to states array if it's not already there
    !    and keep track of the upper state's index. A new hfs state means
    !    that we have to add its energy and degeneracy to the non hfs state
    call find_state_number(up_hfs, states_hfs, iup)
    if(iup .eq. 0) then
      call add_to(up_hfs, states_hfs)
      call sort_last_state(states_hfs)
      iup = size(states_hfs, 1)
      ! -- new hfs state, add to corresponding non hfs state
      select type(state_nohfs => states_nohfs(iup_nohfs))
      type is (asymtop_state_nohfs)
        state_nohfs % E      = state_nohfs % E      + Eup    * degenup
        state_nohfs % degen  = state_nohfs % degen  + degenup
      class default
        call die("Wrong type in upper asymtop non hfs state assignment")
      end select
    endif

    ! -- add the current transition's Einstein coefficient to the higher state's total Einstein coefficient
    !    HFS and non HFS
    states_nohfs(iup_nohfs) % EinstA  = states_nohfs(iup_nohfs) % EinstA  + EinstA * degenup
    states_nohfs(iup_nohfs) % sigmaA2 = states_nohfs(iup_nohfs) % sigmaA2 + sigmaA2
    states_hfs(iup)         % EinstA  = states_hfs(iup)         % EinstA  + EinstA
    states_hfs(iup)         % sigmaA2 = states_hfs(iup)         % sigmaA2 + sigmaA2

    call add_to(transitionul, transitions) ! -- assumes unique lines and does not check if we already have them in the array

    call sort_last_transition(transitions)

  enddo

  ! -- finalize the non hfs states
  do i = 1, size(states_nohfs, 1)
    select type(state_nohfs => states_nohfs(i))
    type is (asymtop_state_nohfs)
      state_nohfs % E      = state_nohfs % E      / state_nohfs % degen
      state_nohfs % EinstA = state_nohfs % EinstA / state_nohfs % degen
    class default
      call die("Wrong type in asymtop state finalization !")
    end select
  enddo

  ! -- sort the non hfs states
  call sort_states(states_nohfs)

  call write_states(stdout, states_hfs)
  write(stdout, *)
  call write_states(stdout, states_nohfs)
  write(stdout, *)
  call write_transitions(stdout, transitions)

! ================================================================================================================================ !
end program CDMSreader
! ================================================================================================================================ !
