subroutine calc_potential(x, nwalk, mass, v)
  use iso_fortran_env, only: real64
!  use mpi
  real(real64), parameter :: au_to_wn = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)
  real(real64), parameter :: pi = acos(-1d0)
  !real(real64), parameter :: zpe_66 = 791.6373314827983d0 / au_to_wn ! J = 0
  !real(real64), parameter :: zpe_68 = 769.3702558802487d0 / au_to_wn ! J = 0
  !real(real64), parameter :: zpe_88 = 746.4311399198617d0 / au_to_wn ! J = 0
  !real(real64), parameter :: De = 9274.99560025014d0 / au_to_wn
  ! ryans stuff
  real(real64), parameter :: shift = 0d0
  parameter (natom = 3)
  parameter (ndim = 3*natom)
  integer, intent(in) :: nwalk
  double precision, dimension(nwalk,9), intent(in) :: x
  double precision, dimension (nwalk), intent(out) :: v
  ! end ryans stuff
  real(real64), intent(in) :: mass(3)
  real(real64) :: potential
  real(real64) :: r(3), r2(3), mass_copy(3)
  ! external :: INT_Cart, Cart_INT, IMLS ! Dawes' procedures
  external :: Cart_INT, IMLS ! Dawes' procedures
  mass_copy = mass
  ! call INT_Cart(r2, cart, mass_copy, 4) ! APH to Cartesian
  do k = 1, nwalk
     call Cart_INT(r2, x(k,:), mass_copy, 2) ! Cartesian to bond-angle
     r(1) = min(r2(1), r2(2))
     r(2) = max(r2(1), r2(2))
     r(3) = 180*acos(r2(3)) / pi
     call IMLS(r, potential, 1)
     ! potential = potential / au_to_wn - De + shift
     v(k) = potential
  enddo !k
end