module nao_grad_unit
USE eddi, ONLY: pi
USE lebedev, ONLY: dp, lebedev_grid, get_number_of_lebedev_grid
USE grid, ONLY: build_onecenter_grid, radial_grid
USE gradients, ONLY: jacobian, grad_twocenter, grad_twocenter_fd,&
                     grad_kinetic, grad_kinetic_fd
implicit none
contains

subroutine test_kinetic_grad()
   REAL(KIND=dp), DIMENSION(5000) :: r, y1, wr, y2
   REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
   REAL(KIND=dp), DIMENSION(1000,3) :: error
   INTEGER :: l1, l2, m1, m2, i, c

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-0.5_dp * r**2)

   c = 0
   error = 0._dp
   d12 = (/ .1_dp, -.5_dp, .1_dp /)
   do l1=0,1
   do l2=l1,1
   do m1=-l1,l1
   do m2=-l2,l2
   ! print *, '   ', l1, m1, l2, m2
         c = c+1
         call grad_kinetic(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                           nshell=(/100, 100/), d12=d12, grad=grad1)


         call grad_kinetic_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                              nshell=(/100, 100/), d12=d12, grad=grad2)

         error(c, :) = abs(grad1-grad2)
         if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
         if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
         if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
         ! if ( all( error(c, :)  .gt. 0.1_dp  )) then
            print *, '   ', l1, m1, l2, m2
            print *, 'e  ', grad1
            print *, 'fd ', grad2
            print *, 'ra ', error(c, :)
            print *,
         ! endif
   enddo
   enddo
   enddo
   print *,
   enddo
   print *, 'mean error', sum(error, 1)/c
   print *, 'min error', minval(error(1:c, :), 1)
   print *, 'max error', maxval(error(1:c, :), 1)
end subroutine test_kinetic_grad

subroutine test_twocenter_grad()
   REAL(KIND=dp), DIMENSION(1000) :: r, y1, wr, y2
   REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
   REAL(KIND=dp), DIMENSION(1000,3) :: error
   INTEGER :: l1, l2, m1, m2, i, c

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-0.5_dp * r**2)

   c = 0
   error = 0._dp
   d12 = (/ .1_dp, -.5_dp, .1_dp /)
   do l1=0,5
   do l2=l1,5
   do m1=-l1,l1
   do m2=-l2,l2
   ! print *, '   ', l1, m1, l2, m2
         c = c+1
         call grad_twocenter(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                             nshell=(/100, 100/), d12=d12, grad=grad1)

         call grad_twocenter_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                                nshell=(/100, 100/), d12=d12, grad=grad2)

         error(c, :) = abs(grad1-grad2)
         if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
         if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
         if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
         if ( all( error(c, :)  .gt. 0.1_dp  )) then
            print *, '   ', l1, m1, l2, m2
            print *, 'e  ', grad1
            print *, 'fd ', grad2
            print *, 'ra ', error(c, :)
            print *,
         endif
   enddo
   enddo
   enddo
   print *,
   enddo
   print *, 'error', sum(error, 1)/c
end subroutine test_twocenter_grad

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