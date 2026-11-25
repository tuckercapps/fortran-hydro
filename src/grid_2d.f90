program grid_2d
  implicit none

  ! Simplest case would be to start with the advection equation
  ! u_t + a*u_x + a*u_y= 0
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
  integer :: j
  integer :: t
  
  ! Define the 2D grid
  real, dimension(length, length) :: uold = 0.0
  real, dimension(length, length) :: unew = 0.0

  ! Initialize uold with a top hat function
  do i=1, length
    do j=1, length
        if ((i-1) * delta_x >= 0.25 .and. (i-1) * delta_x <= 0.75 .and. (j-1) * delta_x >= 0.25 .and. (j-1) * delta_x <= 0.75) then
            uold(i, j) = 1.0
        else
            uold(i, j) = 0.0
        end if
    end do 
  end do

  !print *, uold
  uold(1,:) = uold(length-1,:)  ! ghost cells
  uold(:,1) = uold(:,length-1)
  uold(length,:) = uold(2,:)  ! ghost cells
  uold(:,length) = uold(:,2)

  ! Time stepping loop
  do t=1, num_steps
    do i=2, length-1
        do j=2, length-1
            ! Looping over all physical cells. The first physical cells are down wind of the ghost first cells.
            unew(i,j) = uold(i,j) - a * delta_t / delta_x * (uold(i,j) - uold(i-1,j)) - a * delta_t / delta_x * (uold(i,j) - uold(i,j-1))
        end do
    end do

    ! Update ghost cells
    unew(1,:) = unew(length-1,:)  ! the first ghost column is equal to the last physical column
    unew(:,1) = unew(:,length-1)  ! the first ghost row is equal to the last physical row
    unew(length,:) = unew(2,:)    ! the last ghost column is equal to the first physical column
    unew(:,length) = unew(:,2)   ! the last ghost row  is equal to the first physical row

    ! Update uold for the next time step
    uold = unew
  end do  




  

  

end program grid_2d

