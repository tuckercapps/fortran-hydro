module reconstruction
    implicit none
    
contains
    subroutine piecewise_constant(u, uL_x, uR_x, uL_y, uR_y)
        
        integer :: nx
        integer :: ny
        real, dimension(:,:), intent(in) :: u
        real, dimension(:,:), intent(out) :: uL_x
        real, dimension(:,:), intent(out) :: uR_x
        real, dimension(:,:), intent(out) :: uL_y
        real, dimension(:,:), intent(out) :: uR_y

        nx = size(u, 1)
        ny = size(u, 2)

        ! Spatial reconstruction in x
        uL_x = u(1:nx-1, 2:ny-1)
        uR_x = u(2:nx, 2:ny-1)

        ! Spatial reconstruction in y
        uL_y = u(2:nx-1, 1:ny-1)
        uR_y = u(2:nx-1, 2:ny)

    end subroutine piecewise_constant


    
end module reconstruction