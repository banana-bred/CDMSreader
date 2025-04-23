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
