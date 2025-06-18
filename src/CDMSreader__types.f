! ================================================================================================================================ !
module CDMSreader__types
  !! Contains type definitions and procedures for those types

  use, intrinsic :: iso_fortran_env, only: real32, real64

  implicit none

  private

  ! -- derived types
  public :: asymtop_state
  public :: asymtop_state_hfs
  public :: asymtop_state_nohfs
  public :: asymtop_transition
  public :: asymtop_transition_hfs
  public :: asymtop_transition_nohfs

  ! -- routines
  public :: sort_last_transition
  public :: add_to
  public :: find_state_number
  public :: find_transition_number
  public :: sort_last_state
  public :: sort_states
  public :: sort_transitions
  public :: make_asymtop_state

  integer, parameter, public :: sp = real32
  integer, parameter, public :: dp = real64

  type, abstract :: asymtop_state
    !! Corresponds to an asymmetric top molecule defined by Q = 23 without the hyperfine splitting,
    !! essentially averaged over F
    integer :: dN
      !! Twice the rotational quantum number of the molecule as a rigid rotor
    integer :: dKa
      !! Twice the approximate projection of N on the A axis
    integer :: dKc
      !! Twice the approximate projection of N on the C axis
    real(dp) :: E
      !! The state energy
    real(dp) :: EinstA
      !! The total Einstein coefficient from this state
  end type asymtop_state

  type, extends(asymtop_state) :: asymtop_state_nohfs
    !! Already defined by asymtop state !
    integer :: degen
      !! The total hyperfine degeneracy of the state Σ(2F+1)
  end type asymtop_state_nohfs

  type, extends(asymtop_state) :: asymtop_state_hfs
    !! Corresponds to an asymmetric top molecule defined by Q = 23 with the hyperfine splitting
    integer :: dJ
      !! Twice the total angular momentum of the rotation (N) and the electron spin (S).
      !! \(\vec{J} = \vec{N} + \vec{S}\)
    integer :: dItot
      !! Twice the nuclear spin quantum number
    integer :: dF
      !! Twice the angular momentum from the rotation and electrons (J) and the nuclear spin (I).
      !! \(\vec{F} = \vec{J} + \vec{I} =  \vec{N} + \vec{S} + \vec{I} \)
  end type asymtop_state_hfs

  type, abstract :: asymtop_transition
    !! Corresponds to a transition
    real(dp) :: freq
      !! The frequency of the transition
    real(dp) :: EinstA
      !! The Einstein coefficient A
    real(dp) :: err
      !! The error
    integer :: dr
      !! Degrees of freedom
    integer :: gup
      !! Upper level degeneracy
  end type asymtop_transition

  type, extends(asymtop_transition) :: asymtop_transition_hfs
    type(asymtop_state_hfs) :: up
      !! Upper state
    type(asymtop_state_hfs) :: lo
      !! Lower stat
  end type asymtop_transition_hfs

  type, extends(asymtop_transition) :: asymtop_transition_nohfs
    type(asymtop_state_nohfs) :: up
      !! Upper state
    type(asymtop_state_nohfs) :: lo
      !! Lower stat
  end type asymtop_transition_nohfs

  interface operator(.precedes.)
    module procedure :: precedes_state
    module procedure :: precedes_transition
  end interface operator(.precedes.)

  interface assignment(=)
    module procedure :: state_set_eq
    module procedure :: transition_set_eq
  end interface assignment(=)

  interface operator(.eq.)
    module procedure :: state_iseq
    module procedure :: transition_iseq
  end interface operator(.eq.)

  interface operator(.ne.)
    module procedure :: state_ne
    module procedure :: transition_ne
  end interface operator(.ne.)

  interface add_to
    module procedure :: add_state_to
    module procedure :: add_transition_to
  end interface add_to

  interface make_asymtop_state
    module procedure :: make_asymtop_state_hfs
    module procedure :: make_asymtop_state_nohfs
  end interface make_asymtop_state

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function state_iseq(state1, state2) result(res)
    !! Checks if two states are equal
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_state), intent(in) :: state1
    class(asymtop_state), intent(in) :: state2
    logical :: res
    res = .false.
    if(state1 % dN    .ne. state2 % dN)    return
    if(state1 % dKa   .ne. state2 % dKa)   return
    if(state1 % dKc   .ne. state2 % dKc)   return
    select type(s1 => state1)
    type is (asymtop_state_hfs)
      select type(s2 => state2)
      type is (asymtop_state_hfs)
        if(s1 % dJ    .ne. s2 % dJ   ) return
        if(s1 % dItot .ne. s2 % dItot) return
        if(s1 % dF    .ne. s2 % dF)    return
      class default
        call die("Trying to test equality between a hfs and non hfs state !")
      end select
    end select
    res = .true.
  end function state_iseq

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function transition_iseq(transition1, transition2) result(res)
    !! Checks if two transitions are equal
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_transition), intent(in) :: transition1
    class(asymtop_transition), intent(in) :: transition2
    logical :: res
    res = .false.
    select type(transition1)
    type is (asymtop_transition_hfs)
      select type(transition2)
      type is (asymtop_transition_hfs)
        if(transition1 % lo .ne. transition2 % lo) return
        if(transition1 % up .ne. transition2 % up) return
      class default
        call die("Trying to test equality between an hfs and non hfs state in a transition !")
      end select
    type is (asymtop_transition_nohfs)
      select type(transition2)
      type is (asymtop_transition_nohfs)
        if(transition1 % lo .ne. transition2 % lo) return
        if(transition1 % up .ne. transition2 % up) return
      class default
        call die("Trying to test equality between an hfs and non hfs state in a transition !")
      end select
    end select
    res = .true.
  end function transition_iseq

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function state_ne(state1, state2) result(res)
    !! Checks if two states are not equal
    implicit none
    class(asymtop_state), intent(in) :: state1
    class(asymtop_state), intent(in) :: state2
    logical :: res
    res = .not. (state1 .eq. state2)
  end function state_ne

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function transition_ne(transition1, transition2) result(res)
    !! Checks if two transitions are not equal
    implicit none
    class(asymtop_transition), intent(in) :: transition1
    class(asymtop_transition), intent(in) :: transition2
    logical :: res
    res = .not. (transition1 .eq. transition2)
  end function transition_ne

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine find_state_number(state, states, i)
    !! Check if the state "state" is in the array of states "states"
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_state), intent(in) :: state
    class(asymtop_state), intent(in), allocatable :: states(:)
    integer, intent(out) :: i
    i = 0
    if(.not. allocated(states)) return
    do i = 1, size(states, 1)
      if( state .ne. states(i) ) cycle
      return
    enddo
    i = 0
  end subroutine find_state_number

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine find_transition_number(transition, transitions, i)
    !! Check if the transition "transition" is in the array of transitions "transitions"
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_transition), intent(in) :: transition
    class(asymtop_transition), intent(in), allocatable :: transitions(:)
    integer, intent(out) :: i
    i = 0
    if(.not. allocated(transitions)) return
    do i = 1, size(transitions, 1)
      if(transition .ne. transitions(i)) cycle
      return
    enddo
    i = 0
  end subroutine find_transition_number

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function state_is_in(state, states) result(res)
    !! Check if the state "state" is in the array of states "states"
    implicit none
    class(asymtop_state), intent(in) :: state
    class(asymtop_state), intent(in), allocatable :: states(:)
    logical :: res
    integer :: i
    res = .false.
    if(.not. allocated(states)) return
    do i = 1, size(states, 1)
      if(state .ne. states(i)) cycle
      res = .true.
      return
    enddo
  end function state_is_in

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function transition_is_in(transition, transitions) result(res)
    !! Check if the transition "transition" is in the array of transitions "transitions"
    implicit none
    class(asymtop_transition), intent(in) :: transition
    class(asymtop_transition), intent(in), allocatable :: transitions(:)
    logical :: res
    integer :: i
    res = .false.
    if(.not. allocated(transitions)) return
    do i = 1, size(transitions, 1)
      if(transition .ne. transitions(i)) cycle
      res = .true.
      return
    enddo
  end function transition_is_in

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine add_state_to(state, states)
    !! Add state "state" to the array of states "states"
    implicit none
    class(asymtop_state), intent(in) :: state
    class(asymtop_state), intent(inout), allocatable :: states(:)
    class(asymtop_state), allocatable :: tmp(:)
    integer :: n
    if(.not. allocated(states)) then
      allocate(states(1), source = state)
      states(1) = state
    else
      n = size(states, 1)
      call move_alloc(states, tmp)
      allocate(states(n+1), source = state)
      states(1:n) = tmp(1:n)
      states(n+1) = state
      deallocate(tmp)
    endif
  end subroutine add_state_to
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure module subroutine add_transition_to(transition, transitions)
    !! Add transition "transition" to the array of transitions "transitions"
    implicit none
    class(asymtop_transition), intent(in) :: transition
    class(asymtop_transition), intent(inout), allocatable :: transitions(:)
    class(asymtop_transition), allocatable :: tmp(:)
    integer :: n
    if(.not. allocated(transitions)) then
      allocate(transitions(1), source = transition)
      transitions(1) = transition
      return
    endif
    n = size(transitions, 1)
    call move_alloc(transitions, tmp)
    allocate(transitions(n+1), source = transition)
    transitions(1:n) = tmp(1:n)
    transitions(n+1) = transition
    deallocate(tmp)
  end subroutine add_transition_to

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function precedes_state(state1, state2) result(res)
    !! Test whether state 1 precedes state 2 based on their energies
    implicit none
    class(asymtop_state), intent(in) :: state1
    class(asymtop_state), intent(in) :: state2
    logical :: res
    res = .true.
    if(state1 % E .lt. state2 % E) return
    res = .false.
  end function precedes_state
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function precedes_transition(transition1, transition2) result(res)
    !! Test whether transition 1 precedes transition 2 in the sense where
    !! Elo(1) < Elo(2) and Eup(1) < Eup(2)
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_transition), intent(in) :: transition1
    class(asymtop_transition), intent(in) :: transition2
    logical :: res
    res = .false.
    select type (t1 => transition1)
    type is (asymtop_transition_hfs)
      select type (t2 => transition2)
      type is (asymtop_transition_hfs)
        if(t1 % lo .precedes. t2 % lo) then
          res = .true. ; return
        elseif(t1 % up .precedes. t2 % up) then
          res = .true. ; return
        endif
        res = .false. ; return
      class default
        call die("Trying to compare an hfs transition with a non hfs transition !")
      end select
    type is (asymtop_transition_nohfs)
      select type (t2 => transition2)
      type is (asymtop_transition_nohfs)
        if(t1 % lo .precedes. t2 % lo) then
          res = .true. ; return
        elseif(t1 % up .precedes. t2 % up) then
          res = .true. ; return
        endif
        res = .false. ; return
      class default
        call die("Trying to compare an hfs transition with a non hfs transition !")
      end select
    end select
    res = .true.
  end function precedes_transition

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine sort_last_transition(transitions, iout)
    !! Sort the last element in the array to where it should go, assuming the rest of the array is sorted.
    !! Optionally return the index of the sorted element
    use CDMSreader__system, only: die

    implicit none

    class(asymtop_transition), intent(inout), allocatable, target :: transitions(:)
      !! The array of transitions
    integer, intent(out), optional :: iout
      !! The array to which the last state was sorted

    integer :: i, k, n
    class(asymtop_transition), allocatable :: last
    class(asymtop_transition), allocatable :: tmp(:)

    n = size(transitions, 1)

    allocate(last, source = transitions(n))
    last = transitions(n)

    i = findloc( (last .precedes. transitions(1:n-1)) &
               , value = .true.             &
               , dim = 1                    &
        )

    select case(i)
    case(:-1)
      call die("Findloc returned a negative index in transition sort !")

    case(0)
      ! -- locical mask is all .false. ; last is last
      if(present(iout)) iout = n
      return

    case(1)
      ! transitions = [last, transitions(1:n-1)]
      ! -- locical mask is all .true. ; last is first
      call move_alloc(transitions, tmp)
      allocate(transitions(n), source = last)
      transitions(1)         = last
      transitions(2:n)       = tmp(1:n-1)
      if(present(iout)) iout = 1
      deallocate(tmp)

    case default

      ! transitions = [transitions(1:i-1), last, transitions(i:n-1)]

      call move_alloc(transitions, tmp)
      allocate(transitions(n), source = last)
      transitions(1:i-1)     = tmp(1:i-1)
      transitions(i)         = last
      transitions(i+1:n)     = tmp(i:n-1)
      if(present(iout)) iout = i
      deallocate(tmp)


    end select

  end subroutine sort_last_transition

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine sort_last_state(states, iout)
    !! Sort the last element in the array to where it should go, assuming the rest of the array is sorted
    use CDMSreader__system, only: die

    implicit none

    class(asymtop_state), intent(inout), allocatable, target :: states(:)
      !! The array of states
    integer, intent(out), optional :: iout
      !! The array to which the last state was sorted
    class(asymtop_state), allocatable :: last
    class(asymtop_state), allocatable :: tmp(:)
    integer :: i, k, n

    n = size(states, 1)

    allocate(last, source = states(n))
    last = states(n)

    i = findloc( (last .precedes. states(1:n-1)), value = .true., dim = 1 )

    select case(i)
    case(:-1)
      call die("Findloc returned a negative index in state sort !")

    case(0)
      ! -- locical mask is all .false. ; last is last
      if(present(iout)) iout = n
      return

    case(1)
      ! -- locical mask is all .true. ; last is first
      call move_alloc(states, tmp)
      allocate(states(n), source = last)
      states(1)   = last
      states(2:n) = tmp(1:n-1)
      deallocate(tmp)
      if(present(iout)) iout = 1

    case default

      call move_alloc(states, tmp)
      allocate(states(n), source = last)
      states(1:i-1)   = tmp(1:i-1)
      states(i)       = last
      states(i+1:n) = tmp(i:n-1)
      if(present(iout)) iout = i
      deallocate(tmp)

    end select

    if(n .ne. size(states)) call die("The size of the states array has changed during the sort !")

  end subroutine sort_last_state

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine sort_states(states)
    !! Insertion sort the array states, without assuming that it is sorted
    use CDMSreader__system, only: die

    implicit none

    class(asymtop_state), intent(inout), allocatable, target :: states(:)
      !! The array of states
    class(asymtop_state), allocatable :: state
    class(asymtop_state), allocatable :: tmp(:)
    integer :: i, k, n, pos

    n = size(states, 1)

    allocate(state,  source = states(n))
    allocate(tmp(n), source = states(n))

    tmp(1) = states(1)

    outer: do i = 2, n
      state = states(i)
      pos = i

      inner: do k = i-1, 1, -1

        if (state .precedes. tmp(k)) then
          tmp(k+1) = tmp(k)
          pos = k
          cycle inner
        end if

        exit inner

      enddo inner
      tmp(pos) = state
    enddo outer

    call move_alloc(tmp, states)

  end subroutine sort_states

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure subroutine sort_transitions(transitions)
    !! Insertion sort the array transitions, without assuming that it is sorted
    use CDMSreader__system, only: die

    implicit none

    class(asymtop_transition), intent(inout), allocatable, target :: transitions(:)
      !! The array of transitions
    class(asymtop_transition), allocatable :: transition
    class(asymtop_transition), allocatable :: tmp(:)
    integer :: i, k, n, pos

    n = size(transitions, 1)

    allocate(transition,  source = transitions(n))
    allocate(tmp(n), source = transitions(n))

    tmp(1) = transitions(1)

    outer: do i = 2, n
      transition = transitions(i)
      pos = i

      inner: do k = i-1, 1, -1

        if (transition .precedes. tmp(k)) then
          tmp(k+1) = tmp(k)
          pos = k
          cycle inner
        end if

        exit inner

      enddo inner
      tmp(pos) = transition
    enddo outer

    call move_alloc(tmp, transitions)

  end subroutine sort_transitions

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental subroutine swap_state(state1, state2)
    implicit none
    class(asymtop_state), intent(inout) :: state1, state2
    class(asymtop_state), allocatable :: tmp
    call state_set_eq(tmp, state1)
    call state_set_eq(state1, state2)
    call state_set_eq(state2, tmp)
    deallocate(tmp)
  end subroutine swap_state

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental subroutine state_set_eq(state1, state2)
    !! Define assignment for two abstract asymtop states
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_state), intent(out) :: state1
    class(asymtop_state), intent(in)  :: state2
    state1 % dN     = state2 % dN
    state1 % dKa    = state2 % dKa
    state1 % dKc    = state2 % dKc
    state1 % E      = state2 % E
    state1 % EinstA = state2 % EinstA
    select type (s1 => state1)
    ! -- two hfs states
    type is (asymtop_state_hfs)
      select type (s2 => state2)
      type is (asymtop_state_hfs)
        s1 % dJ    = s2 % dJ
        s1 % dItot = s2 % dItot
        s1 % dF    = s2 % dF
      class default
        call die("Attempting to assign an hfs state to a non hfs state !")
      end select
    ! -- two non hfs states
    type is (asymtop_state_nohfs)
      select type (s2 => state2)
      type is (asymtop_state_nohfs)
        s1 % degen = s2 % degen
      class default
        call die("Attempting to assign a non hfs state to an hfs state !")
      end select
    end select
  end subroutine state_set_eq
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental subroutine transition_set_eq(transition1, transition2)
    !! Define assignment for two abstract asymtop transitions
    use CDMSreader__system, only: die
    implicit none
    class(asymtop_transition), intent(out) :: transition1
    class(asymtop_transition), intent(in)  :: transition2
    transition1 % freq   = transition2 % freq
    transition1 % EinstA = transition2 % EinstA
    transition1 % err    = transition2 % err
    transition1 % dr     = transition2 % dr
    transition1 % gup    = transition2 % gup
    select type (t1 => transition1)
    ! -- two hfs transitions
    type is (asymtop_transition_hfs)
      select type (t2 => transition2)
      type is (asymtop_transition_hfs)
        t1 % up    = t2 % up
        t1 % lo    = t2 % lo
      class default
        call die("Attempting to assign an hfs transition to a non hfs transition !")
      end select
    ! -- two non hfs transitions
    type is (asymtop_transition_nohfs)
      select type (t2 => transition2)
      type is (asymtop_transition_nohfs)
        t1 % up    = t2 % up
        t1 % lo    = t2 % lo
      class default
        call die("Attempting to assign a non hfs transition to an hfs transition !")
      end select
    end select
  end subroutine transition_set_eq

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function make_asymtop_state_hfs(dN, dKa, dKc, dJ, dItot, dF, E, EinstA) result(state)
    implicit none
    integer,  intent(in) :: dN, dKa, dKc, dJ, dItot, dF
    real(dp), intent(in) :: E, EinstA
    type(asymtop_state_hfs) :: state
    state % dN     = dN
    state % dKa    = dKa
    state % dKc    = dKc
    state % dJ     = dJ
    state % dItot  = dItot
    state % dF     = dF
    state % E      = E
    state % EinstA = EinstA
  end function make_asymtop_state_hfs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental module function make_asymtop_state_nohfs(dN, dKa, dKc, E, EinstA, degen) result(state)
    implicit none
    integer,  intent(in) :: dN, dKa, dKc
    integer,  intent(in), optional :: degen
    real(dp), intent(in) :: E, EinstA
    type(asymtop_state_nohfs) :: state
    state % dN     = dN
    state % dKa    = dKa
    state % dKc    = dKc
    state % E      = E
    state % EinstA = EinstA
    if(.true. .eqv. present(degen)) then
      state % degen = degen
    else
      state % degen = 0
    endif
  end function make_asymtop_state_nohfs

end module CDMSreader__types
! ================================================================================================================================ !
