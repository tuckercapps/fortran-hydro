module numerical_flux
  implicit none
  
  ! This module contains (or will contrain) functions 
  ! that are used for flux calculation
contains

  ! Rusanov numerical flux
  elemental function rusanov(u_L, u_R, a) result(flux)

    real, intent(in) :: u_L, u_R, a
    real :: flux

    flux = 0.5 * (a*u_L + a*u_R) &
         - 0.5 * abs(a) * (u_R - u_L)

  end function rusanov
  
end module numerical_flux 
