program advection
  implicit none

  ! Simplest case would be to start with the advection equation
  ! u_t + a*u_x = 0
  ! Start with periodic boundary conditions [0, 1]
  ! the number of elements will be decided by delta_x

  real, parameter :: delta_x = 0.01
  integer, parameter :: length = int(1.0/delta_x) + 2

  ! Time stepping parameters
  real, parameter :: delta_t = 0.005
  integer, parameter :: num_steps = 200

  ! Advection speed
  real, parameter :: a = 1.0

  integer :: i
  !character(len=30) :: filename
  real :: uold(length) = [(0.0, i=1, length)]
  real :: unew(length) = [(0.0, i=1, length)]

  ! Initialize uold with a top hat function
  do i=1, length
    if ((i-1) * delta_x >= 0.25 .and. (i-1) * delta_x <= 0.75) then
      uold(i) = 1.0
    else
      uold(i) = 0.0
    end if
  end do

  uold(1) = uold(length-1)  ! ghost cells
  uold(length) = uold(2)  ! ghost cells

  ! Time stepping loop
  do i=1, num_steps

    ! unew(1) = uold(1) - a * delta_t / delta_x * (uold(1) - uold(length))  ! Periodic boundary condition
    unew(2:length-1) = uold(2:length-1) - a * delta_t / delta_x * (uold(2:length-1) - uold(1:length-2))

    unew(1) = unew(length-1)  ! update ghost cells
    unew(length) = unew(2)  ! update ghost cells

    ! Update uold for the next time step
    uold = unew
  end do  




  

  

end program advection

