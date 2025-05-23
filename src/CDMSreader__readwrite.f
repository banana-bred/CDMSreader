! ================================================================================================================================ !
module CDMSreader__readwrite

  implicit none

  private

  public :: CDMS_readline
  public :: write_states
  public :: write_transitions

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

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
  use CDMSreader__types,     only: dp, asymtop_transition
  use CDMSreader__system,    only: die
  use CDMSreader__constants, only: au2ev, au2invcm

  implicit none

  integer, intent(in) :: funit
    !! The file unit
  type(asymtop_transition), intent(in) :: transitions(:)
    !! The array of transitions

  real(dp), parameter :: float_lbound = 1e-2_dp
  real(dp), parameter :: float_ubound = 9.9e5_dp

  integer :: i, n
  real(dp) :: E, err, EinstA, lifetime
  character(6) :: charNup, charKaup, charKcup, charJup, charItotup, charFup
  character(6) :: charNlo, charKalo, charKclo, charJlo, charItotlo, charFlo
  character(15) :: charE, charerr, charA, charlifetime
  character(36) :: header_fmt = '(A, 6(A6, X), 2X, A, 6(A6, X), 4A15)'
  character(22) :: body_fmt   = '(2X, 6A6, A, 6A6, 4A6)'
  character(7) :: float_fmt = '(F15.6)'
  character(7) :: exp_fmt   = '(E15.6)'

  n = size(transitions, 1)

  write(funit, header_fmt) "# "                                        &
        , "N",       "Ka",        "Kc",      "J",  "Itot",  "F", "<--" &
        , "N'",      "Ka'",       "Kc'",     "J'", "Itot'", "F'"       &
        , "E (meV)", "err (meV)", "A (s⁻¹)", "τ (s)"

  do i = 1, n

    write(funit, '(2X)', advance = 'no')
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dN)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dKa)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dKc)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dJ)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dItot)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % lo % dF)
    write(funit, '(2X, A)', advance = 'no') "<--"
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dN)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dKa)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dKc)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dJ)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dItot)
    write(funit, '(A6, X)', advance = 'no') doubleint2char(transitions(i) % up % dF)
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
      write(funit, float_fmt, advance = 'no') lifetime
    else
      write(funit, exp_fmt, advance = 'no') lifetime
    endif

    ! -- newline
    write(funit, *)

    ! write(funit, body_fmt) &
    !        charNup, charKaup, charKcup, charJup, charItotup, charFup&
    !      , charNlo, charKalo, charKclo, charJlo, charItotlo, charFlo&
    !      , charE,   charerr,  charA,    charlifetime
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

  n = size(states, 1)

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
