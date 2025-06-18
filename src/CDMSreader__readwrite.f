! ================================================================================================================================ !
module CDMSreader__readwrite

  implicit none

  private

  public :: CDMS_readfile
  public :: CDMS_readfile_hfs
  public :: CDMS_readfile_nohfs

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine CDMS_readfile_hfs(funit_in, funit_out, states_hfs, transitions_hfs)
    !! Read a file containing transitions from the CDMS, get ony states and transitions with hyperfine resolution.
    !! Wrapper for CDMS_readfile.
    use CDMSreader__types, only: asymtop_state, asymtop_transition
    implicit none
    integer, intent(in) :: funit_in
      !! File unit for the CDMS data
    integer, intent(in) :: funit_out
      !! File unti for the processed data. If funit is 0, do not write output
    class(asymtop_state),     intent(inout), allocatable :: states_hfs(:)
      !! Array of states with hyperfine resolution
    class(asymtop_transition), intent(inout), allocatable :: transitions_hfs(:)
      !! Array of transitions with hyperfine resolution

    class(asymtop_state),     allocatable :: states_nohfs(:)
      !! Dummy array of states statistically averaged over hyperfine levels
    class(asymtop_transition), allocatable :: transitions_nohfs(:)
      !! Dummy array of transitions avveraged over hyperfine levels

    call CDMS_readfile(funit_in, funit_out, .true., .false., states_hfs, states_nohfs, transitions_hfs, transitions_nohfs)

  end subroutine CDMS_readfile_hfs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine CDMS_readfile_nohfs(funit_in, funit_out, states_nohfs, transitions_nohfs)
    !! Read a file containing transitions from the CDMS, get ony states and transitions without hyperfine resolution.
    !! Wrapper for CDMS_readfile.
    use CDMSreader__types, only: asymtop_state, asymtop_transition
    implicit none
    integer, intent(in) :: funit_in
      !! File unit for the CDMS data
    integer, intent(in) :: funit_out
      !! File unti for the processed data. If funit is 0, do not write output
    class(asymtop_state),     intent(inout), allocatable :: states_nohfs(:)
      !! Array of states statistically averaged over hyperfine levels
    class(asymtop_transition), intent(inout), allocatable :: transitions_nohfs(:)
      !! Array of transitions avveraged over hyperfine levels

    class(asymtop_state),  allocatable :: states_hfs(:)
      !! Dummy array of states with hyperfine resolution
    class(asymtop_transition), allocatable :: transitions_hfs(:)
      !! Dummy array of transitions with hyperfine resolution

    call CDMS_readfile(funit_in, funit_out, .false., .true., states_hfs, states_nohfs, transitions_hfs, transitions_nohfs)

  end subroutine CDMS_readfile_nohfs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine CDMS_readfile( funit_in, funit_out, return_hfs, return_nohfs &
                                        , states_hfs, states_nohfs, transitions_hfs, transitions_nohfs)
    !! Read a file containing transitions from the CDMS

    use CDMSreader__types,             only : dp, asymtop_state_hfs, asymtop_transition, add_to &
                                            , sort_last_transition, find_state_number, sort_last_state, asymtop_state_nohfs &
                                            , asymtop_state, make_asymtop_state, sort_states &
                                            , asymtop_transition, asymtop_transition_nohfs, asymtop_transition_hfs &
                                            , find_transition_number, sort_transitions
    use CDMSreader__system,            only : die
    use CDMSreader__constants,         only : invcm2Hz
    use, intrinsic :: iso_fortran_env, only : iostat_end

    implicit none

    integer, intent(in) :: funit_in
      !! File unit for the CDMS data. Must refer to a stream, like stdin or an opened file
    integer, intent(in) :: funit_out
      !! File unti for the processed data. If funit is 0, do not write output. Must refer to
      !! a stream, like stdout or an opened file
    logical, intent(in) :: return_hfs
      !! Return states and transition arrays for hyperfine data ?
    logical, intent(in) :: return_nohfs
      !! Return states and transition arrays for non-hyperfine data ?
    class(asymtop_state),     intent(inout), allocatable :: states_hfs(:)
      !! Array of states with hyperfine resolution
    class(asymtop_state),     intent(inout), allocatable :: states_nohfs(:)
      !! Array of states statistically averaged over hyperfine levels
    class(asymtop_transition), intent(inout), allocatable :: transitions_hfs(:)
      !! Array of transitions with hyperfine resolution
    class(asymtop_transition), intent(inout), allocatable :: transitions_nohfs(:)
      !! Array of transitions avveraged over hyperfine levels

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

    type(asymtop_state_nohfs) :: up_nohfs
    type(asymtop_state_nohfs) :: lo_nohfs
    ! type(asymtop_state_nohfs),  allocatable :: states_nohfs(:)

    type(asymtop_transition_hfs)   :: transition_hfs
    type(asymtop_transition_nohfs) :: transition_nohfs

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
    integer :: itran

    lines: do

      call CDMS_readline(funit_in, freq, err, EinstA, dr, elo, gup, tag, qnfmt, QN(1:12), mol, io, skip, "#")

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
      up_nohfs = make_asymtop_state( dN = dNup, dKa = dKaup, dKc = dKcup, E = 0.0_dp, EinstA = 0.0_dp, degen = 0 )
      lo_nohfs = make_asymtop_state( dN = dNlo, dKa = dKalo, dKc = dKclo, E = 0.0_dp, EinstA = 0.0_dp, degen = 0 )

      ! -- add non HFS states to array if desired
      if(return_nohfs .eqv. .true.) then

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
      endif

      ! -- create the state with hyperfine splitting
      up_hfs = make_asymtop_state( up_nohfs % dN, up_nohfs % dKa, up_nohfs % dKc, dJ = dJup, dItot = dItotup, dF = degenup - 1 &
                                 , E = eup, EinstA = 0.0_dp &
           )
      lo_hfs = make_asymtop_state( lo_nohfs % dN, lo_nohfs % dKa, lo_nohfs % dKc, dJ = dJlo, dItot = dItotlo, dF = degenlo - 1 &
                                 , E = elo, EinstA = 0.0_dp &
           )
      transition_hfs = asymtop_transition_hfs( up = up_hfs, lo = lo_hfs, freq = freq, EinstA = EinstA, err = err &
                                             , dr = dr, gup = gup                                                &
           )

      if( return_nohfs .eqv. .true. ) then
        ! -- creat non HFS transition (with appropriate weighting so that we can
        !    properly finalize these
        transition_nohfs = asymtop_transition_nohfs( up = up_nohfs, lo = lo_nohfs &
                                                   , freq = freq * degenup        &
                                                   , EinstA = EinstA * degenup    &
                                                   , err = (err*degenup)**2       &
                                                   , dr = dr, gup = gup)
      endif

      ! -- add a non HFS transition to the array of transitions if requested
      if( return_nohfs .eqv. .true. ) then
        call find_transition_number(transition_nohfs, transitions_nohfs, itran)
        if(itran .eq. 0) then
          call add_to(transition_nohfs, transitions_nohfs)
        else
          ! -- if this transition already exists, accumulate values for averaging later
          associate(t => transition_nohfs, ti => transitions_nohfs(itran))
            ti % gup    = ti % gup    + t % gup
            ti % freq   = ti % freq   + t % freq
            ti % EinstA = ti % EinstA + t % EinstA
            ti % err    = ti % err    + t % err
          end associate
        endif
      endif

      ! -- add lower hfs state to states array if it's not already there
      !    and keep track of the lower state's index. A new hfs state means
      !    that we have to add its energy and degeneracy to the non hfs state
      call find_state_number(lo_hfs, states_hfs, ilo)
      if(ilo .eq. 0) then
        call add_to(lo_hfs, states_hfs)
        call sort_last_state(states_hfs, ilo)
        ! -- new hfs state, add to corresponding non hfs state
        if(return_nohfs .eqv. .true.) then
          select type(state_nohfs => states_nohfs(ilo_nohfs))
          type is (asymtop_state_nohfs)
            state_nohfs % E     = state_nohfs % E     + Elo * degenlo
            state_nohfs % degen = state_nohfs % degen + degenlo
          class default
            call die("Wrong type in lower asymtop non hfs state assignment")
          end select
        endif
      endif

      ! -- add upper hfs state to states array if it's not already there
      !    and keep track of the upper state's index. A new hfs state means
      !    that we have to add its energy and degeneracy to the non hfs state
      call find_state_number(up_hfs, states_hfs, iup)
      if(iup .eq. 0) then
        call add_to(up_hfs, states_hfs)
        call sort_last_state(states_hfs, iup)
        ! -- new hfs state, add to corresponding non hfs state
        if(return_nohfs .eqv. .true.) then
          select type(state_nohfs => states_nohfs(iup_nohfs))
          type is (asymtop_state_nohfs)
            state_nohfs % E      = state_nohfs % E      + Eup    * degenup
            state_nohfs % degen  = state_nohfs % degen  + degenup
          class default
            call die("Wrong type in upper asymtop non hfs state assignment")
          end select
        endif
      endif

      ! -- add the current transition's Einstein coefficient to the higher state's total Einstein coefficient
      !    HFS and non HFS
      if(return_nohfs .eqv. .true.) states_nohfs(iup_nohfs) % EinstA  = states_nohfs(iup_nohfs) % EinstA  + EinstA * degenup
      states_hfs(iup) % EinstA  = states_hfs(iup)         % EinstA  + EinstA

      if(return_hfs .eqv. .true.) then
        call add_to(transition_hfs, transitions_hfs) ! -- assumes unique lines and does not check if we already have them in the array
        call sort_last_transition(transitions_hfs)
      endif

    enddo lines

    write_nohfs: if(return_nohfs .eqv. .true.) then

      ! -- finalize the non hfs states
      finalize_states_nohfs: do i = 1, size(states_nohfs, 1)
        select type(state_nohfs => states_nohfs(i))
        type is (asymtop_state_nohfs)
          state_nohfs % E      = state_nohfs % E      / state_nohfs % degen
          state_nohfs % EinstA = state_nohfs % EinstA / state_nohfs % degen
        class default
          call die("Wrong type in asymtop state finalization !")
        end select
      enddo finalize_states_nohfs

      ! -- sort these now that their energies are finalized
      call sort_states(states_nohfs)

      ! -- finalize the non hfs transitions
      finalize_transitions_nohfs: do i = 1, size(transitions_nohfs, 1)

        select type(ti => transitions_nohfs(i))
        type is (asymtop_transition_nohfs)

          ! -- divide the accumulated quantities by their summed weights
          ti % EinstA = ti % EinstA    / ti % gup
          ti % freq   = ti % freq      / ti % gup
          ti % err    = sqrt(ti % err) / ti % gup

          ! -- update the energy information of the states involved in
          !    the transition so that the transitions can be sorted
          select type(states => states_nohfs)
          type is (asymtop_state_nohfs)

            call find_state_number(ti%up, states_nohfs, iup)
            call find_state_number(ti%lo, states_nohfs, ilo)
            ti % up = states(iup)
            ti % lo = states(ilo)

            if(iup .eq. 0 .OR. ilo .eq. 0) then
              call die("Unfound state in the array of states during non-HFS transition finalization !")
            endif

          class default

            call die("Somehow, states_nohfs is NOT of type asymtop_state_nohfs !")

          end select

        class default

          call die("Wrong type in asymtop transition finalization !")

        end select

      enddo finalize_transitions_nohfs

      ! -- sort these now that their energies are finalized
      call sort_transitions(transitions_nohfs)

      call write_states(funit_out, states_nohfs)

    endif write_nohfs

    ! call sort_transitions(transitions_hfs)
    if(return_hfs   .eqv. .true.) call write_states(funit_out,      states_hfs)
    if(return_nohfs .eqv. .true.) call write_transitions(funit_out, transitions_nohfs)
    if(return_hfs   .eqv. .true.) call write_transitions(funit_out, transitions_hfs)

    ! -- deallocate arrays that don't need to be returned
    if(return_hfs .eqv. .false.) then
      deallocate(states_hfs)
      deallocate(transitions_hfs)
    endif
    if(return_nohfs .eqv. .false.) then
      deallocate(states_nohfs)
      deallocate(transitions_nohfs)
    endif

  end subroutine CDMS_readfile

! ---------------------------------------------------------------------------------------------------------------------------------!
impure module subroutine CDMS_readline(funit, freq, err, EinstA, dr, elo, gup, tag, qnfmt, QN, mol, io, skip, comment_char_in)
  !!  Read "cat" file from the CDMS contents into appropriate arrays, with the capability to skip lines that are well-commented or blank
  !!  the comment character defaults to "#", but can be set to anything not in the character NUMERIC also not whitespace

  use CDMSreader__types,             only : dp
  use CDMSreader__system,            only : die
  use, intrinsic :: iso_fortran_env, only : iostat_end

  implicit none

  integer,      intent(in)  :: funit
  integer,      intent(out) :: DR
    !! Degrees of freedom in the rotational partition function (0 for atoms, 2 for linear molecules, and 3 for nonlinear molecules)
  integer,      intent(out) :: GUP
    !! Upper state degeneracy
  integer,      intent(out) :: TAG
    !! Species tag or molecular identifier. A negative value flags that the line frequency has been measured in the laboratory.
    !! The absolute value of TAG is then the species tag (as given in line 2 of file.int above) and ERR is the reported experimental
    !! error.
  integer, intent(out) :: qnfmt
    !! Identifies the format of the quantum numbers given in the field QN.
  real(dp),     intent(out) :: freq
    !! Frequency of the line
  real(dp),     intent(out) :: err
    !! Estimated or experimental error (999.9999 indicates error is larger)
  real(dp),     intent(out) :: EinstA
    !! The Einstein A coefficient
  real(dp),     intent(out) :: ELO
    !! Lower state energy in cm–1
  integer, intent(out) :: qn(12)
    !! Twice the quantum numbers numbers. These are integers but are converted from their coding format
    !! according to QNFMT. Upper state quanta start in character 1. Lower state quanta start in character 14 (element 7).
    !! Unused quanta are blank, quanta whose magnitude is larger than 99 or smaller than –9 are shown with alphabetic characters
    !! or **. Quanta between –10 and –19 are shown as a0 through a9. Similarly, –20 is b0, etc., up to –259, which is shown as z9.
    !! Quanta between 100 and 109 are shown as A0 through A9. Similarly, 110 is B0, etc., up to 359, which is shown as Z9.
  character(:), intent(out), allocatable :: mol
  integer, intent(out) :: io
  logical, intent(out) :: skip
  character(1), intent(in),  optional :: comment_char_in
  character(2) :: qnchar(12)
    !! The quantum numbers as characters (this is what is read from the CDMS file)

  character(41), parameter :: CDMS_fmt = "(A13, 2A11, I2, A10, I3, I7, I4, 12A2, A)" ! Valid for einstein coeffs
  ! character(53), parameter :: CDMS_fmt = "(F13.6, F11.7, F11.4, I2, F10.4, I3, I7, I4, 12A2, A)" ! Valid for einstein coeffs
  character(13), parameter :: numeric = "0123456789.+-"
  character(65), parameter :: alphanumeric = "ABCDEFGHIJKLMNOPQRSTUVWXYZ&
                                             &abcdefghijklmnopqrstuvwxyz&
                                             &0123456789.+-"

  character(1) :: comment_char
  integer :: commentStart
  integer :: alphanumericStart
  integer :: alphanumericEnd
  character(100) :: line
  character(13) :: freqchar
  character(11) :: errchar, EinstAchar
  character(10) :: elochar

  skip = .false.

  ! -- set different comment character maybe
  comment_char = "#" ; if(present(comment_char_in)) comment_char = comment_char_in

  read(funit, "(A)", iostat=io) line

  if(io.eq.iostat_end) return
  if(io.ne.0) call die("Problem reading the line: " // line)

  ! -- avoid comments and prune the line
  alphanumericEnd   = scan(line,alphanumeric, .true.) ! -- find position of last  alphanumeric character
  commentStart      = scan(line,comment_char) ! -- find position of first comment character
  if(commentStart .gt. 0) then
    line = line(1:min(commentStart-1, alphanumericEnd)) ! -- prune line
  else
    line = line(1:alphanumericEnd)
  endif
  if(commentStart.gt.0 .AND. commentStart.lt.1) skip = .true. ! -- cycle reading if the line appears commented out
  if(trim(line) .eq. "") skip = .true.

  if(skip .eqv. .true.) return

  allocate(character(100) :: mol)

  ! -- read into variables
  read(line, CDMS_fmt) freqchar, errchar, EinstAchar, dr, elochar, gup, tag, qnfmt, qnchar(1:12), mol

  ! -- char -> real
  read(freqchar, *)   freq
  read(EinstAchar, *) EinstA
  read(errchar, *)    err
  read(elochar, *)    elo

  ! -- convert and trim output
  qn = charQN2int(qnchar)
  mol = trim(mol)

end subroutine CDMS_readline

! ---------------------------------------------------------------------------------------------------------------------------------!
impure elemental function charQN2int(QNchar) result(res)
  !! Converts the CDMS 2-character representation of integers to actual integers
  use CDMSreader__system, only : die

  implicit none

  character(2), intent(in) :: QNchar
  integer :: res

  character(1), parameter :: uppercase(26) = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U"&
                                             ,"V","W","X","Y","Z"]
  character(1), parameter :: lowercase(26) = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"&
                                             ,"v","w","x","y","z"]
  character(1), parameter :: integers(10) = ["0","1","2","3","4","5","6","7","8","9"]

  integer :: int1, int2

  ! -- first character
  if(any( QNchar(1:1) .eq. uppercase )) then
    int1 = ichar(QNchar(1:1)) - ichar("A") + 10
  elseif(any( QNchar(1:1) .eq. lowercase )) then
    int1 = ichar(to_uppercase(QNchar(1:1))) - ichar("A") + 10
    int1 = -int1
  elseif(QNchar(1:1) .eq. " ") then
    int1 = 0
  elseif(any( QNchar(1:1) .eq. integers )) then
    read(QNchar(1:1), "(I1)") int1
  else
    call die("Could not determine the quantum number " // QNchar)
  endif

  ! -- second character
  if(all( QNchar(2:2) .ne. integers )) call die("Could not determine the quantum number " // QNchar)

  read(QNchar(2:2), "(I1)") int2

  res = int1*10 + int2

contains

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function to_uppercase(str) result(res)
    implicit none
    character(*), intent(in) :: str
    character(len(str)) :: res
    integer :: i,n
    integer :: ic
    res = str
    n = len(res)
    do i=1,n
      ic = ichar(res(i:i))
      if(ic.ge.97 .AND. ic .le.122) res(i:i) = char(ic-32)
    enddo
  end function to_uppercase

end function charQN2int

! -------------------------------------------------------------------------------------------------------------------------------- !
subroutine write_transitions(funit, transitions)
  !! Writes the transitions in a legible format to the designated file unit
  use CDMSreader__types,     only: dp, asymtop_transition, asymtop_transition_hfs, asymtop_transition_nohfs
  use CDMSreader__system,    only: die
  use CDMSreader__constants, only: au2ev, au2invcm

  implicit none

  integer, intent(in) :: funit
    !! The file unit
  class(asymtop_transition), intent(in) :: transitions(:)
    !! The array of transitions

  real(dp), parameter :: float_lbound = 1e-2_dp
  real(dp), parameter :: float_ubound = 9.9e5_dp

  integer :: i, n
  real(dp) :: E, err, EinstA, lifetime
  character(6) :: charNup, charKaup, charKcup, charJup, charItotup, charFup
  character(6) :: charNlo, charKalo, charKclo, charJlo, charItotlo, charFlo
  character(15) :: charE, charerr, charA, charlifetime
  character(36) :: header_fmt_hfs   = '(A, 6(A6, X), 2X, A, 6(A6, X), 4A15)'
  character(36) :: header_fmt_nohfs = '(A, 3(A6, X), 2X, A, 3(A6, X), 4A15)'
  character(22) :: body_fmt_hfs     = '(2X, 6A6, A, 6A6, 4A6)'
  character(22) :: body_fmt_nohfs   = '(2X, 3A6, A, 3A6, 4A6)'
  character(7) :: float_fmt = '(F15.6)'
  character(7) :: exp_fmt   = '(E15.6)'

  if(funit .eq. 0) return

  ! -- initial space
  write(funit, *)

  n = size(transitions, 1)

  ! -- change header based on first element, assume all transitions are of same type
  select type(t => transitions(1))
  type is (asymtop_transition_hfs)
    write(funit, header_fmt_hfs) "# "                                    &
          , "N",       "Ka",        "Kc",      "J",  "Itot",  "F", "<--" &
          , "N'",      "Ka'",       "Kc'",     "J'", "Itot'", "F'"       &
          , "E (meV)", "err (meV)", "A (s⁻¹)", "τ (s)"
  type is (asymtop_transition_nohfs)
    write(funit, header_fmt_nohfs) "# "         &
          , "N",       "Ka",        "Kc", "<--" &
          , "N'",      "Ka'",       "Kc'"       &
          , "E (meV)", "err (meV)", "A (s⁻¹)", "τ (s)"
  class default
    call die("Unidentified type for the first element of the transitions array")
  end select

  do i = 1, n

    select type(t => transitions(i))
    type is (asymtop_transition_hfs)
      write(funit, '(2X)', advance = 'no')
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dN)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dKa)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dKc)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dJ)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dItot)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dF)
      write(funit, '(2X, A)', advance = 'no') "<--"
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dN)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dKa)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dKc)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dJ)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dItot)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dF)
    type is (asymtop_transition_nohfs)
      write(funit, '(2X)', advance = 'no')
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dN)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dKa)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % lo % dKc)
      write(funit, '(2X, A)', advance = 'no') "<--"
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dN)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dKa)
      write(funit, '(A6, X)', advance = 'no') doubleint2char(t % up % dKc)
    end select

    E        = transitions(i) % freq
    err      = transitions(i) % err
    EinstA   = transitions(i) % EinstA
    lifetime = 1/EinstA

    ! -- invcm -> meV
    E   = E   / au2invcm * au2ev * 1000
    err = err / au2invcm * au2ev * 1000

    ! -- write energy to character
    if(E .eq. 0 .OR. (E .ge. float_lbound .AND. E .le. float_ubound)) then
      write(funit, float_fmt, advance = 'no') E
    else
      write(funit, exp_fmt, advance = 'no') E
    endif

    ! -- write error to character
    if(err .eq. 0 .OR. (err .ge. float_lbound .AND. err .le. float_ubound)) then
      write(funit, float_fmt, advance = 'no') err
    else
      write(funit, exp_fmt, advance = 'no') err
    endif

    ! -- write Einstein A to character
    if(EinstA .eq. 0 .OR. (EinstA .ge. float_lbound .AND. EinstA .le. float_ubound)) then
      write(funit, float_fmt, advance = 'no') EinstA
    else
      write(funit, exp_fmt, advance = 'no') EinstA
    endif

    ! -- write lifetime to character
    if(lifetime .eq. 0 .OR. (lifetime .ge. float_lbound .AND. lifetime .le. float_ubound)) then
      write(funit, float_fmt) lifetime
    else
      write(funit, exp_fmt) lifetime
    endif

  enddo

end subroutine write_transitions

! -------------------------------------------------------------------------------------------------------------------------------- !
subroutine write_states(funit, states)
  !! Writes the states in a legible format to the designated file unit
  use CDMSreader__types,     only: asymtop_state, asymtop_state_hfs, asymtop_state_nohfs, dp
  use CDMSreader__system,    only: die
  use CDMSreader__constants, only: au2ev, au2invcm

  implicit none

  integer, intent(in) :: funit
    !! The file unit
  class(asymtop_state), intent(in) :: states(:)
    !! The array of states

  integer :: n, i
  integer :: degen
  real(dp) :: E, EinstA, lifetime
  real(dp), parameter :: float_lbound = 1e-2_dp
  real(dp), parameter :: float_ubound = 9.9e5_dp
  character(:), allocatable :: charN, charKa, charKc, charJ, charItot, charF
  character(19) :: header_fmt_hfs     = '(A, 6A6, A15, A15)'
  character(19) :: header_fmt_nohfs   = '(A, 4A6, A15, A15)'
  character(16) :: body_fmt_hfs       = '(2X, 6A6)'
  character(21) :: body_fmt_nohfs     = '(2X, 3A6, I6)'

  if(funit .eq. 0) return

  n = size(states, 1)

  ! -- initial empty line
  write(funit, *)

  ! -- choose appropriate header
  select type (s1 => states(1))
  type is (asymtop_state_hfs)
    write(funit, header_fmt_hfs) "# ", "N", "Ka", "Kc", "J", "Itot", "F", "energy (meV)", "lifetime (s)"
  type is (asymtop_state_nohfs)
    write(funit, header_fmt_nohfs) "# ", "N", "Ka", "Kc", "degen", "energy (meV)", "lifetime (s)"
  class default
    call die("Invalid type for element 1 of the states array in the write procedure")
  end select

  do i = 1, n
    charN    = doubleint2char(states(i) % dN)
    charKa   = doubleint2char(states(i) % dKa)
    charKc   = doubleint2char(states(i) % dKc)
    E = states(i) % E
    EinstA = states(i) % EinstA
    lifetime = 1/EinstA

    ! -- invcm -> meV
    E = E / au2invcm * au2ev * 1000

    ! -- print the corresponding infor for this type of state
    select type (si => states(i))
    type is (asymtop_state_nohfs)

      degen = si % degen
      write(funit, body_fmt_nohfs, advance = "no") charN, charKa, charKc, degen
      if(E.eq.0 .OR. (E .ge. float_lbound .AND. E .lt. float_ubound)) then
        write(funit, '(F15.6)', advance = "no") E
      else
        write(funit, '(E15.6)', advance = "no") E
      endif


    type is (asymtop_state_hfs)

      charJ    = doubleint2char(si % dJ)
      charItot = doubleint2char(si % dItot)
      charF    = doubleint2char(si % dF)
      write(funit, body_fmt_hfs, advance = "no") charN, charKa, charKc, charJ, charItot, charF

      if(E.eq.0 .OR. (E .ge. float_lbound .AND. E .lt. float_ubound)) then
        write(funit, '(F15.6)', advance = "no") E
      else
        write(funit, '(E15.6)', advance = "no") E
      endif

    class default

      call die("Invalid type detected in the states array in the write procedure")

    end select

    ! -- print the lifetimes with their calculated uncertainties
    if(EinstA .eq. 0) then
      write(funit, '(A15)') "inf"
    elseif(lifetime .lt. float_lbound .OR. lifetime .ge. float_ubound) then
      write(funit, '(E15.6)') lifetime
    else
      write(funit, '(F15.6)') lifetime
    endif

  enddo

end subroutine write_states

! -------------------------------------------------------------------------------------------------------------------------------- !
pure function doubleint2char(i) result(res)
  !! Writes the value i/2 to a character. If i is even, write i/2 as an '(I0)'.
  !! If i is odd, write i/2 as '(I0, "/", I0)'

  implicit none

  integer, intent(in) :: i
  character(:), allocatable :: res

  allocate(character(10) :: res)

  if(mod(i, 2) .eq.0) then
    write(res, '(I0)') i/2
    res = trim(res)
    return
  endif

  write(res, '(I0, "/", I0)') i, 2
  res = trim(res)

end function doubleint2char

! ================================================================================================================================ !
end module CDMSreader__readwrite
! ================================================================================================================================ !
