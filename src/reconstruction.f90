module reconstruction
    implicit none
    
contains
    subroutine piecewise_constant(u, uL_x, uR_x, uL_y, uR_y, nx, ny)
        
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        real, dimension(nx, ny), intent(in) :: u
        real, dimension(nx-1, ny-1), intent(out) :: uL_x
        real, dimension(nx-1, ny-1), intent(out) :: uR_x
        real, dimension(nx-1, ny-1), intent(out) :: uL_y
        real, dimension(nx-1, ny-1), intent(out) :: uR_y

        ! Spatial reconstruction in x
        uL_x = u(1:nx-1, :)
        uR_x = u(2:nx, :)

        ! Spatial reconstruction in y
        uL_y = u(:, 1:ny-1)
        uR_y = u(:, 2:ny)

    end subroutine piecewise_constant


    
end module reconstruction