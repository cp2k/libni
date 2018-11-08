module nao_grad_unit
USE eddi, ONLY: pi
USE lebedev, ONLY: dp, lebedev_grid, get_number_of_lebedev_grid
USE grid, ONLY: build_onecenter_grid
USE gradients, ONLY: jacobian
implicit none
contains

! The way to test the jacobian:
!  a) make sure, that grad_rad(1/r) = grad_xyz(1/r)
!     grad_rad ( 1/r ) = (-1/r^2, 0, 0)
!     grad_xyz ( 1/r ) = (-x/r^3, -y/r^3, -z/r^3)
subroutine test_jacobian()
   implicit none
   REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w, errors
   REAL(KIND=dp), DIMENSION(3) :: grad_r, grad_xyz, grad_xyz_exact
   REAL(KIND=dp), DIMENSION(3,3) :: jac
   REAL(KIND=dp) :: norm, y, xx, yy, zz, theta, phi
   INTEGER :: i, ileb, ngrid

   ileb = get_number_of_lebedev_grid(l=0)
   ngrid = lebedev_grid(ileb)%n * 25
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(errors(ngrid))

   call build_onecenter_grid(ileb=1, nshell=25, addr2=.TRUE.,&
                             quadr=1, grid_r=grid_r, grid_w=grid_w)
   do i=1,size(grid_w)
      grad_r = 0._dp; grad_xyz = 0._dp; grad_xyz_exact = 0._dp
      norm = sqrt( sum(grid_r(i, :)**2) )
      xx = grid_r(i, 1)
      yy = grid_r(i, 2)
      zz = grid_r(i, 3)

      theta = acos(zz/norm)
      phi = atan2(yy, xx)

      grad_r(1) = -norm**(-2._dp)
      grad_xyz_exact(1) = -xx/norm**3
      grad_xyz_exact(2) = -yy/norm**3
      grad_xyz_exact(3) = -zz/norm**3

      jac = jacobian(r=norm, theta=theta, phi=phi)
      grad_xyz = matmul(jac, grad_r)
      errors(i) = sqrt(sum((grad_xyz-grad_xyz_exact)**2))
   enddo
   if (sum(errors)/size(errors) .lt. 1.e-12_dp) then
      print *, '👌 test_jacobian - passed'
   else
      print *, '💣 test_jacobian - failed'
      print *, 'mean error: ', sum(errors)/size(errors)
   endif
   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(errors)
end subroutine test_jacobian

end module nao_grad_unit