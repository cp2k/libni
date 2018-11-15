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
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r, grid_dw, tmp_grad
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3) :: dr, dtheta, dphi, dylm2
   REAL(KIND=dp) :: norm, x, y, z, theta, phi, rho,&
                     f1, f2, ylm1, ylm2, df2, dwdr2

   ileb(1) = get_number_of_lebedev_grid(n=302)
   ileb(2) = get_number_of_lebedev_grid(n=302)

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(grid_r(ngrid, 3))
   allocate(tmp_grad(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(grid_dw(ngrid, 3))

   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.TRUE., grid_r=grid_r, grid_w=grid_w)
   grad = 0._dp
   tmp_grad = 0._dp
   do i=1,ngrid
      !  N                ~             ~
      !  Σ  w  * f(r) * f(r) * Y(r) * Y(r) ,
      ! i=1  i    1 i    2 i      i      i

      !       ~    ––>     ––>                  T
      ! where r = grid_r - d12 = (x-X, y-Y, z-Z)

      ! When taking the derivative wrt X, Y, Z - f1 and ylm1 won't change
      norm = sqrt( sum( grid_r(i, :)**2 ) )
      call interpolation(gr=r1, gy=y1, spline=s1, r=norm, y=f1)
      call rry_lm(l=l(1), m=m(1), r=grid_r(i, :)/norm, y=ylm1)

      ! f2, ylm2 do have non-vanishing derivatives
      norm = sqrt( sum( (grid_r(i, :) - d12)**2 ) )
      x = grid_r(i, 1)-d12(1); y = grid_r(i, 2)-d12(2); z = grid_r(i, 3)-d12(3);
      rho = x*x + y*y
      theta = acos(z/norm); phi = atan2(y, x);

      call interpolation(gr=r2, gy=y2, spline=s2, r=norm, y=f2, yprime=df2)
      call rry_lm(l=l(2), m=m(2), r=(grid_r(i, :) - d12)/norm, y=ylm2)
      call dry_lm(l=l(2), m=m(2), c=(/theta, phi/), dy=dylm2)
      ! dwdr2 = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))

      ! The partial derivatives dr/dXYZ, dtheta/XYZ, dphi/XYZ
      dr = -(grid_r(i, :) - d12)/norm

      dtheta = 0._dp
      if (rho .ne. 0._dp) then
         dtheta(1) = -x*z
         dtheta(2) = -y*z
         dtheta(3) = rho 
         dtheta = dtheta / ( norm**3 * sqrt(1._dp - z**2/norm**2) )
      endif

      dphi = 0._dp
      if (rho .ne. 0._dp) then
         dphi(1) = y
         dphi(2) = -x
         dphi = dphi/rho
      endif

      ! On to the actual calculation
      ! X: There will be 4 terms. For now we skip dw/dr, since it's nasty. TODO
      tmp_grad(i, 1) = grid_w(i) * df2 * dr(1) * ylm2 +&
                       grid_w(i) * f2 * dylm2(1) * dtheta(1) +&
                       grid_w(i) * f2 * dylm2(2) * dphi(1)

      ! Y: There will be 4 terms. For now we skip dw/dr, since it's nasty. TODO
      tmp_grad(i, 2) = grid_w(i) * df2 * dr(2) * ylm2 +&
                       grid_w(i) * f2 * dylm2(1) * dtheta(2) +&
                       grid_w(i) * f2 * dylm2(2) * dphi(2)

      ! Z: There will be 3 terms. For now we skip dw/dr, since it's nasty. TODO
      tmp_grad(i, 3) = grid_w(i) * df2 * dr(3) * ylm2 +&
                       grid_w(i) * f2 * dylm2(1) * dtheta(3)

      tmp_grad(i, :) = tmp_grad(i, :) * f1 * ylm1
   enddo

   grad = 0._dp
   grad(1) = sum(tmp_grad(:,1))
   grad(2) = sum(tmp_grad(:,2))
   grad(3) = sum(tmp_grad(:,3))

   deallocate(grid_r)
   deallocate(tmp_grad)
   deallocate(grid_w)
   deallocate(grid_dw)
end subroutine grad_twocenter

subroutine grad_twocenter_fd(r1, y1, r2, y2, l, m, nshell, d12, grad)
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

   d12t = d12 + 2._dp*ex*h
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
   d12t = d12 - 2._dp*ex*h
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(1,4))


   d12t = d12 + 2._dp*ey*h
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
   d12t = d12 - 2._dp*ey*h
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(2,4))


   d12t = d12 + 2._dp*ez*h
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
   d12t = d12 - 2._dp*ez*h
   call integration_twocenter(l=l, m=m, nshell=nshell, d12=d12t, &
                              r1=r1, y1=y1, r2=r2, y2=y2,&
                              spline1=s1, spline2=s2, integral=findiff(3,4))


   grad(1) = -findiff(1,1) + 8._dp*findiff(1,2) - 8._dp*findiff(1,3) + findiff(1,4)
   grad(2) = -findiff(2,1) + 8._dp*findiff(2,2) - 8._dp*findiff(2,3) + findiff(2,4)
   grad(3) = -findiff(3,1) + 8._dp*findiff(3,2) - 8._dp*findiff(3,3) + findiff(3,4)
   grad = grad/(12._dp*h)

   ! print *, sum(findiff(1,:))/4._dp/findiff(1,1)

   ! print *, findiff(2,:)
   ! print *, sum(findiff(2,:))/4._dp/findiff(2,1)

   ! print *, findiff(3,:)
   ! print *, sum(findiff(3,:))/4._dp/findiff(3,1)
end subroutine grad_twocenter_fd

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
   tmp_grad = 0._dp
   do i=1,size(grid_w)
      ! gradient (r, theta, phi)

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
      call rry_lm(l=l, m=m, r=grid_r(i, :)/norm, y=ylm(i))
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

subroutine grad_onecenter_cart(r, y, l, m, nshell, grad)
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
   REAL(KIND=dp) :: theta, phi, integral, xx, yy, zz, dwdr, alpha, norm, rho
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
   tmp_grad = 0._dp
   do i=1,size(grid_w)
      ! gradient (r, theta, phi)

      ! Define some helpers
      norm = sqrt(sum( grid_r(i, :)**2 ))
      xx = grid_r(i, 1); yy = grid_r(i, 2); zz = grid_r(i, 3)
      rho = xx**2 + yy**2
      theta = acos(zz/norm)
      phi = atan2(yy, xx)

      ! We will need all those
      ! f, df/dr, dw/dr
      call interpolation(gr=r, gy=y, spline=s, r=norm+5._dp*epsilon(1._dp), y=y_loc(i), yprime=dy(i))
      dwdr = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))
      ! Ylm, dYlm/dtheta, dYlm/dphi
      call rry_lm(l=l, m=m, r=grid_r(i, :)/norm, y=ylm(i))
      call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))

      !                                                          T
      ! tmp_grad = (d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi )
      !

      ! dr/dx = x/r, dr/dy = y/r, dr/dz = z/r
      tmp_grad(i,:) = ( dwdr * y_loc(i) + grid_w(i) * dy(i) ) * ylm(i)
      tmp_grad(i,1) = tmp_grad(i,1) * xx/norm
      tmp_grad(i,2) = tmp_grad(i,2) * yy/norm
      tmp_grad(i,3) = tmp_grad(i,3) * zz/norm

      ! dtheta/dx = x*z/(r^3*sqrt(1-z^2/r^2))
      ! dtheta/dy = y*z/(r^3*sqrt(1-z^2/r^2))
      ! dtheta/dz = -sqrt(rho^2/r^2)/r
      if (xx*zz .ne. 0._dp) then
      tmp_grad(i,1) = tmp_grad(i,1) + grid_w(i) * y_loc(i) * dylm(i, 1)&
         * xx*zz/( norm**3 * sqrt( 1._dp-(zz/norm)**2 ) )
      endif
      if (yy*zz .ne. 0._dp) then
      tmp_grad(i,2) = tmp_grad(i,2) + grid_w(i) * y_loc(i) * dylm(i, 1)&
         * yy*zz/( norm**3 * sqrt( 1._dp-(zz/norm)**2 ) )
      endif
      tmp_grad(i,3) = tmp_grad(i,3) + grid_w(i) * y_loc(i) * dylm(i, 1)&
         * sqrt((xx*xx+yy*yy)/norm**2)/norm

      ! dphi/dx = -y/(x^2+y^2)
      ! dphi/dy = x/(x^2+y^2)
      ! dphi/dz = 0
      if (yy .ne. 0._dp) then
         tmp_grad(i,1) = tmp_grad(i,1) - grid_w(i) * y_loc(i) * dylm(i, 2) * (yy/rho)
      endif
      if (xx .ne. 0._dp) then
         tmp_grad(i,2) = tmp_grad(i,2) + grid_w(i) * y_loc(i) * dylm(i, 2) * (xx/rho)
      endif
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
end subroutine grad_onecenter_cart

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