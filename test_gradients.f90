module nao_grad_unit
USE eddi, ONLY: pi, spline
USE lebedev, ONLY: dp, lebedev_grid, get_number_of_lebedev_grid
USE grid, ONLY: build_onecenter_grid, radial_grid
USE gradients, ONLY: jacobian, grad_twocenter, grad_twocenter_fd,&
                     grad_kinetic, grad_kinetic_fd,&
                     grad_coulomb, grad_coulomb_fd
implicit none
contains

subroutine test_coulomb_grad()
   REAL(KIND=dp), DIMENSION(250) :: r, y1, wr, y2, s1, s2
   REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
   REAL(KIND=dp), DIMENSION(1000,3) :: error
   INTEGER :: l1, l2, m1, m2, c, n

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-0.5_dp * r**2)

   call spline(r, y1, size(r), s1)
   call spline(r, y2, size(r), s2)

   c = 0
   error = 0._dp
   d12 = (/ 1._dp, 3._dp, 4._dp /)
   n = 250
   do l1=0,1
   do l2=l1,l1
   do m1=-l1,l1
   do m2=-l2,l2
   ! print *, '   ', l1, m1, l2, m2
   ! open(unit=111, file='ipynb/llmm0010_nvar')
      c = c+1
      call grad_coulomb(r1=r, y1=y1, r2=r, y2=y2, s1=s1, s2=s2,&
                        l=(/l1,l2/), m=(/m1,m2/),&
                        nshell=(/n, n/), coul_n=250, d12=d12, grad=grad1)

      call grad_coulomb_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                           nshell=(/n, n/), d12=d12, grad=grad2)

      error(c, :) = abs(grad1-grad2)
      if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
      if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
      if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
      if ( any( error(c, :)  .gt. 0.1_dp  )) then
         print *, '   ', l1, m1, l2, m2
         print *, 'e  ', grad1
         print *, 'fd ', grad2
         print *, 'ra ', error(c, :)
         print *, achar(7)  ! beep, boop
      endif
      ! write (*, *) n, grad1, grad2, sum( error(c,:) )/3._dp
   ! close(111)
   enddo
   enddo
   enddo
   enddo
   print *, 'mean error', sum(error, 1)/c
   print *, 'min error', minval(error(1:c, :), 1)
   print *, 'max error', maxval(error(1:c, :), 1)
end subroutine test_coulomb_grad

subroutine test_kinetic_grad()
   REAL(KIND=dp), DIMENSION(15000) :: r, y1, wr, y2, d1y, d2y, d3y
   REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
   REAL(KIND=dp), DIMENSION(1000,3) :: error
   INTEGER :: l1, l2, m1, m2, c, n

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-r**2)
   d1y = -2._dp * r * y2
   d2y = (4._dp * r**2 - 2._dp) * y2
   d3y = (12._dp * r - 8._dp * r**3) * y2

   c = 0
   error = 0._dp
   d12 = (/ .25_dp, 0._dp, 0._dp /)
   n = 150
   l1 = 0; m1 = 0; l2 = 0; m2 = 0;
   ! do l1=0,1
   ! do l2=l1,1
   ! do m1=-l1,l1
   ! do m2=-l2,l2
   print *, '   ', l1, m1, l2, m2
   ! l1 = 0; m1 = 0; l2 = 1; m2 = 0;
   ! open(unit=111, file='ipynb/llmm0010_nvar')
   ! do n = 10, 400, 30
      ! print *, n
      c = c+1
      call grad_kinetic(r1=r, y1=y1, r2=r, y2=y2, d1y=d1y, d2y=d2y, d3y=d3y,&
                        l=(/l1,l2/), m=(/m1,m2/),&
                        nshell=(/n, n/), d12=d12, grad=grad1)

      print *, 'e  ', grad1
      call grad_kinetic_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                           nshell=(/n, n/), d12=d12, step=1.e-7_dp, grad=grad2)

      error(c, :) = abs(grad1-grad2)
      if(error(c, 1) .lt. 1.e-10) error(c, 1) = 0._dp
      if(error(c, 2) .lt. 1.e-10) error(c, 2) = 0._dp
      if(error(c, 3) .lt. 1.e-10) error(c, 3) = 0._dp
      if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
      if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
      if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
      ! if ( any( error(c, :)  .gt. 1.e-5_dp  )) then
         ! print *, '   ', l1, m1, l2, m2
         print *, 'fd ', grad2
         print *, 'ra ', error(c, :)
         print *, achar(7)  ! beep, boop
      ! endif
      ! write (111, *) n, grad1, grad2, sum( error(c,:) )/3._dp
   ! enddo
   ! close(111)
   ! enddo
   ! enddo
   ! enddo
   ! print *,
   ! enddo
   print *, 'mean error', sum(error, 1)/c
   print *, 'min error', minval(error(1:c, 1), 1, error(1:c, 1).gt.0._dp)
   print *, 'max error', maxval(error(1:c, 1), 1)
end subroutine test_kinetic_grad

subroutine test_kinetic_fd()
   REAL(KIND=dp), DIMENSION(500) :: r, y1, wr, y2
   REAL(KIND=dp), DIMENSION(24, 3) :: grad
   REAL(KIND=dp), DIMENSION(3) :: d12
   REAL(KIND=dp) :: step

   INTEGER :: i

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-0.5_dp * r**2)
   
   d12 = (/ 1._dp, 0._dp, 0._dp /)
   open(unit=111, file='ipynb/test_kinetic_fd')
   do i = 1,size(grad,1)
      step = 10._dp**(-i*0.5)
      call grad_kinetic_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/0,0/), m=(/0,0/),&
                           nshell=(/150, 150/), d12=d12, step=step, grad=grad(i, :))
      print *, step, grad(i, :)
      write (111, *) step, grad(i, :)
   enddo
   close(111)
end subroutine test_kinetic_fd

subroutine test_twocenter_grad()
   REAL(KIND=dp), DIMENSION(1000) :: r, y1, wr, y2
   REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
   REAL(KIND=dp), DIMENSION(1000,3) :: error
   INTEGER :: l1, l2, m1, m2, c, n

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-0.2_dp * r**2)
   y2 = exp(-0.5_dp * r**2)

   c = 0
   error = 0._dp
   d12 = (/ .1_dp, -.5_dp, .1_dp /)
   n = 130
   l1 = 0; m1 = 0; l2 = 1; m2 = 1;
   do l1=0,1
   do l2=l1,2
   do m1=0,l1
   do m2=0,l2
      ! print *, 'n', n
      c = c+1
      call grad_twocenter(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                          nshell=(/n, n/), d12=d12, grad=grad1)

      call grad_twocenter_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                             nshell=(/n, n/), d12=d12, grad=grad2)
      error(c, :) = abs(grad1-grad2)
      if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
      if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
      if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
      ! if ( any( error(c, :)  .gt. 0.1_dp  )) then
         print *, '   ', l1, m1, l2, m2
         print *, 'e  ', grad1
         print *, 'fd ', grad2
         print *, 'ra ', error(c, :), '/3 = ', sum(error(c, :))/3._dp
         ! print *,
      ! else
         ! print *, '👌  ', l1, m1, l2, m2
      ! endif
   enddo
   enddo
   enddo
   enddo
   if(all(sum(error, 1)/c .lt. 0.1_dp)) then
   print *, '👌 test_twocenter_grad - passed'
   print *, 'error', sum(error, 1)/c
   else
   print *, '💣 test_twocenter_grad - failed'
   print *, 'error', sum(error, 1)/c
   endif
end subroutine test_twocenter_grad

subroutine test_twocenter_fd()
   REAL(KIND=dp), DIMENSION(500) :: r, y1, wr, y2
   REAL(KIND=dp), DIMENSION(24, 3) :: grad
   REAL(KIND=dp), DIMENSION(3) :: d12
   REAL(KIND=dp) :: step

   INTEGER :: i

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

   y1 = exp(-r**2)
   y2 = exp(-0.5_dp * r**2)
   
   d12 = (/ 1._dp, 1._dp, 4._dp /)
   open(unit=111, file='ipynb/test_twocenter_fd')
   do i = 1,size(grad,1)
      step = 10._dp**(-i*0.5)
      call grad_twocenter_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/1,1/), m=(/0,-1/),&
                           nshell=(/150, 150/), d12=d12, step=step, grad=grad(i, :))
      print *, step, grad(i, :)
      write (111, *) step, grad(i, :)
   enddo
   close(111)
end subroutine test_twocenter_fd


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
   REAL(KIND=dp) :: norm, xx, yy, zz, theta, phi
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