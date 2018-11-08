module gradients
USE lebedev, ONLY: dp, get_number_of_lebedev_grid, lebedev_grid
USE eddi, ONLY: pi, spline, interpolation, integration_twocenter
USE grid, ONLY: build_onecenter_grid, build_twocenter_grid
USE spherical_harmonics, ONLY: rry_lm, dry_lm

implicit none
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ex = (/ 1._dp, 0._dp, 0._dp /)
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ey = (/ 0._dp, 1._dp, 0._dp /)
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ez = (/ 0._dp, 0._dp, 1._dp /)
contains

subroutine grad_twocenter(r1, y1, r2, y2, l, m, nshell, d12, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb, nleb
   INTEGER :: ngrid, i
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3, 4) :: findiff
   REAL(KIND=dp), DIMENSION(3) :: d12t
   REAL(KIND=dp) :: integral, h

   ! REAL(KIND=dp), DIMENSION(size(r)) :: s
   ! REAL(KIND=dp), DIMENSION(3,3) :: jac
   ! REAL(KIND=dp) :: theta, phi, integral, xx, yy, zz, dwdr, alpha, norm
   ! INTEGER :: i

   ! Get the one-center grid
   ileb(1) = get_number_of_lebedev_grid(l=l(1)+5)
   ileb(2) = get_number_of_lebedev_grid(l=l(2)+5)
   nleb(1) = lebedev_grid(ileb(1))%n
   nleb(2) = lebedev_grid(ileb(2))%n

   ngrid = sum(nleb*nshell)
   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)

   h = 1.0e-6_dp

   d12t = d12 + 2._dp*h*ex
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(1,1))

   d12t = d12 + h*ex
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(1,2))
   d12t = d12 - h*ex
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(1,3))
   d12t = d12 - 2._dp*h*ex
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(1,4))


   d12t = d12 + 2._dp*h*ey
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(2,1))

   d12t = d12 + h*ey
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(2,2))
   d12t = d12 - h*ey
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(2,3))
   d12t = d12 - 2._dp*h*ey
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(2,4))


   d12t = d12 + 2._dp*h*ez
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(3,1))

   d12t = d12 + h*ez
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(3,2))
   d12t = d12 - h*ez
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(3,3))
   d12t = d12 - 2._dp*h*ez
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(3,4))


   grad(1) = -findiff(1,1) + 8._dp*findiff(1,2) - 8._dp*findiff(1,3) + findiff(1,4)
   grad(2) = -findiff(2,1) + 8._dp*findiff(2,2) - 8._dp*findiff(2,3) + findiff(2,4)
   grad(3) = -findiff(3,1) + 8._dp*findiff(3,2) - 8._dp*findiff(3,3) + findiff(3,4)
   grad = grad/(12._dp*h)

   ! print *, findiff(1,:)
   ! print *, sum(findiff(1,:))/4._dp/findiff(1,1)

   ! print *, findiff(2,:)
   ! print *, sum(findiff(2,:))/4._dp/findiff(2,1)

   ! print *, findiff(3,:)
   ! print *, sum(findiff(3,:))/4._dp/findiff(3,1)
end subroutine grad_twocenter

subroutine grad_onecenter(r, y, l, m, nshell, grad)
   implicit none
   ! Input
   INTEGER, intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, y
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r, dylm, tmp_grad
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w, y_loc, dy, ylm
   REAL(KIND=dp), DIMENSION(size(r)) :: s
   REAL(KIND=dp), DIMENSION(3,3) :: jac
   REAL(KIND=dp) :: theta, phi, integral, xx, yy, zz, dwdr, alpha, norm
   INTEGER :: i, ileb, nleb


   ! Get the one-center grid
   ileb = get_number_of_lebedev_grid(l=l+5)
   nleb = lebedev_grid(ileb)%n

   allocate(grid_r(nleb*nshell,3))
   allocate(grid_w(nleb*nshell))

   allocate(y_loc(nleb*nshell))
   allocate(dy(nleb*nshell))

   allocate(ylm(nleb*nshell))
   allocate(dylm(nleb*nshell, 2))

   allocate(tmp_grad(size(grid_w),3))

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE.,&
                             quadr=1, grid_r=grid_r, grid_w=grid_w)
   alpha = pi/(REAL(2*nshell+2, dp))
   ! Evaluate y, dy/dr, Ylm, dYlm/dtheta, dYlm/dphi
   call spline(r=r, y=y, n=size(r), yspline=s)
   do i=1,size(grid_w)
      ! gradient (r, theta, phi)
      tmp_grad(i, :) = 0._dp

      ! Define some helpers
      norm = sqrt(sum( grid_r(i, :)**2 ))
      xx = grid_r(i, 1); yy = grid_r(i, 2); zz = grid_r(i, 3)
      theta = acos(zz/norm)
      phi = atan2(yy, xx)

      ! We will need all those
      ! f, df/dr, dw/dr
      call interpolation(gr=r, gy=y, spline=s, r=norm+5._dp*epsilon(1._dp), y=y_loc(i), yprime=dy(i))
      dwdr = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))
      ! Ylm, dYlm/dtheta, dYlm/dphi
      call rry_lm(l=l, m=m, r=grid_r(i, :), y=ylm(i))
      call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))

      !                                                          T
      ! tmp_grad = (d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi )
      !

      ! d/dr => (dw/dr * f + w * df/dr) * Ylm
      !    dw/dr = pi/(N+1) * (7*r^2.5 + 5*r^1.5)/2 
      !    df/dr <- dy
      tmp_grad(i,1) = ( dwdr * y_loc(i) + grid_w(i) * dy(i) ) * ylm(i)

      ! 1/r d/dtheta => w * f * dYlm/dtheta
      tmp_grad(i,2) = grid_w(i) * y_loc(i) * dylm(i, 1) / norm

      ! 1/(r * sin(theta)) d/dphi => w * f * dYlm/dphi / r*sin(theta)
      if (dylm(i, 2) .ne. 0._dp) then
         tmp_grad(i,3) = grid_w(i) * y_loc(i) * dylm(i, 2) / (norm * sin(theta))
      endif

      ! And transform it to xyz
      jac = jacobian(r=norm, theta=theta, phi=phi)
      tmp_grad(i,:) = matmul(jac, tmp_grad(i,:))
   enddo

   grad = 0._dp
   grad(1) = sum(tmp_grad(:,1))
   grad(2) = sum(tmp_grad(:,2))
   grad(3) = sum(tmp_grad(:,3))

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(y_loc)
   deallocate(dy)
   deallocate(ylm)
   deallocate(dylm)
   deallocate(tmp_grad)
end subroutine grad_onecenter

function jacobian(r, theta, phi)
   implicit none
   REAL(KIND=dp), intent(in) :: r, theta, phi
   REAL(KIND=dp), DIMENSION(3,3) :: jacobian
   REAL(KIND=dp) :: st, ct, sp, cp

   st = sin(theta)
   ct = cos(theta)
   sp = sin(phi)
   cp = cos(phi)
   jacobian(1,1) = st*cp
   jacobian(1,2) = r*ct*cp
   jacobian(1,3) = -r*st*sp

   jacobian(2,1) = st*sp
   jacobian(2,2) = r*ct*sp
   jacobian(2,3) = r*st*cp

   jacobian(3,1) = ct
   jacobian(3,2) = -r*st
   jacobian(3,3) = 0._dp
end function jacobian

end module gradients