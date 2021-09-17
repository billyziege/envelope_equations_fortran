module constants
!-------------------------------------------------------------------------------
! Module ensemble:
! A module to hold and handle parametric constants.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!type my_constants
! A data type containing helpful physical constants
!
!subroutine my_constants_init:
!  -  Creates an instance of constants
! input:
!   -  q (optional) = Charge of particles.  If not provided, assumes e.
!   -  m (opitonal) = Mass of particles.  If not provided, assumes m_e.
! output:
!   -  constants_out = A constants container with helpful constants.
!
!
implicit none

type my_constants !Contains the details of a single instance of similar data.
    sequence
	real*8 :: q
	real*8 :: m 
	real*8 :: epsilon_0
	real*8 :: c
	real*8 :: gam 
end type my_constants

contains

function new_my_constants(gam, q, m)
  real*8, intent(in), optional :: gam, q, m
  type(my_constants) :: new_my_constants

  new_my_constants%c = 299792458
  new_my_constants%epsilon_0 = 8.85418782e-12
  new_my_constants%gam = gam  

  if ( present(q) ) then
    new_my_constants%q = q
  else
    new_my_constants%q = 1.602176634e-19
  end if

  if ( present(m) ) then
    new_my_constants%m = m
  else
    new_my_constants%m =  9.109383701528e-31
  end if
end function new_my_constants

end module constants
