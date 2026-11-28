program grid_2d
  use reconstruction
  use numerical_flux
  implicit none

  ! Simplest case would be to start with the advection equation
  ! u_t + a*u_x + a*u_y= 0
  ! Start with periodic boundary conditions [0, 1]
  ! the number of elements will be decided by delta_x

  real, parameter :: delta_x = 0.1
  real, parameter :: delta_y = 0.1
  integer, parameter :: nx = int(1.0/delta_x) + 2
  integer, parameter :: ny = int(1.0/delta_y) + 2

  ! Time stepping parameters
  real, parameter :: delta_t = 0.05
  integer, parameter :: num_steps = 200

  ! Advection speed
  real, parameter :: a = 1.0
  real, parameter :: a_x = 1.0
  real, parameter :: a_y = 1.0

  ! Gaussian parameters
  real, parameter :: PI = 4*atan(1.0)
  real, parameter :: mean = 0.5
  real, parameter :: std = 0.05

  integer :: i
  integer :: j
  integer :: t

  integer :: u
  character(len=32) :: fname
  
  ! Define the 2D grid
  real, dimension(nx, ny) :: uold = 0.0
  real, dimension(nx, ny) :: unew = 0.0

  real, dimension(nx-1, ny-2) :: uL_x
  real, dimension(nx-1, ny-2) :: uR_x
  real, dimension(nx-2, ny-1) :: uL_y
  real, dimension(nx-2, ny-1) :: uR_y
  real, dimension(nx-1, ny-2) :: flux_x
  real, dimension(nx-2, ny-1) :: flux_y

  ! Initialize uold with a gaussian
  do i=1, nx
    do j=1, ny
        !if ((i-1) * delta_x >= 0.25 .and. (i-1) * delta_x <= 0.75 .and. (j-1) * delta_x >= 0.25 .and. (j-1) * delta_x <= 0.75) then
        uold(i, j) = 1/(std**2*2*PI)*exp(-(((i-1)*delta_x - mean)**2 + ((j-1)*delta_x - mean)**2)/(2*std**2)) 
    end do 
  end do

  !print *, uold
  uold(1,:) = uold(nx-1,:)  ! ghost cells
  uold(:,1) = uold(:,ny-1)
  uold(nx,:) = uold(2,:)  ! ghost cells
  uold(:,ny) = uold(:,2)

  ! Time stepping loop
  do t=1, num_steps+1
    !do i=2, length-1
    !    do j=2, length-1
    !        ! Looping over all physical cells. The first physical cells are down wind of the ghost first cells.
    !        unew(i,j) = uold(i,j) - a * delta_t / delta_x * (uold(i,j) - uold(i-1,j)) - a * delta_t / delta_x * (uold(i,j) - uold(i,j-1))
    !    end do
    !end do

    write(fname, '(a,i4.4,a)') 'u_', t, '.dat'   ! builds u_0001.dat, u_0002.dat, ...
    open(newunit=u, file=fname, status='replace')
    do i = 2, nx-1
      write(u, '(*(es12.4))') uold(i, 2:nx-1)
    end do
    close(u)

    !print *, uL_x(1,1)
    call piecewise_constant(uold, uL_x, uR_x, uL_y, uR_y, nx, ny)
    flux_x = rusanov(uL_x, uR_x, a_x)
    flux_y = rusanov(uL_y, uR_y, a_y)
    print *, shape(uR_x)!
    !printshape(flux_x)
   ! print *, uL_x(1,1)

    !unew(2:nx-1, 2:ny-1) = uold(2:nx-1, 2:ny-1) &
    !    - a_x*delta_t/delta_x*(uold(2:nx-1, 2:ny-1) - uold(2:nx-1, 1:ny-2)) &
    !    - a_y*delta_t/delta_x*(uold(2:nx-1, 2:ny-1) - uold(1:nx-2, 2:ny-1))

    unew(2:nx-1, 2:ny-1) = uold(2:nx-1, 2:ny-1) &
        - delta_t/delta_x*(flux_x(2:nx-2, 2:ny-2) - flux_x(1:nx-1, 1:ny-2)) &
        - delta_t/delta_x*(flux_y(2:nx-2, 2:ny-2) - flux_y(1:nx-2, 2:ny-1))

    

    ! Update ghost cells
    unew(1,:) = unew(nx-1,:)  ! the first ghost column is equal to the last physical column
    unew(:,1) = unew(:,ny-1)  ! the first ghost row is equal to the last physical row
    unew(nx,:) = unew(2,:)    ! the last ghost column is equal to the first physical column
    unew(:,ny) = unew(:,2)   ! the last ghost row  is equal to the first physical row

    ! Update uold for the next time step
    uold = unew
  end do  


end program grid_2d

