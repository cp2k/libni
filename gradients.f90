module gradients
USE lebedev, ONLY: dp, get_number_of_lebedev_grid, lebedev_grid
USE eddi, ONLY: pi, spline, interpolation
USE grid, ONLY: build_onecenter_grid
USE spherical_harmonics, ONLY: rry_lm, dry_lm
contains

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
   REAL(KIND=dp) :: theta, phi, integral, xx, yy, zz, dwdr, alpha, norm
   INTEGER :: i, ileb, nleb


   ! Get the one-center grid
   ileb = get_number_of_lebedev_grid(l=1)!l+5)
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
      call interpolation(gr=r, gy=y, spline=s, r=norm, y=y_loc(i), yprime=dy(i))
      dwdr = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))
      ! Ylm, dYlm/dtheta, dYlm/dphi
      call rry_lm(l=l, m=m, r=grid_r(i, :), y=ylm(i))
      call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))

      !                                                          T
      ! tmp_grad = (d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi )
      ! ============================================================
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

   ! ! Interpolate the function onto r_loc and evaluate df/dr=`dy_r`
   ! call radial_grid(r=r_loc, wr=wr_loc, n=nshell, addr2=.FALSE., quadr=1)
   ! call spline(r=r, y=y, n=size(r), yspline=s)
   ! do i=1,size(r_loc)
   !    call interpolation(gr=r, gy=y, spline=s, r=r_loc(i), y=y_loc(i), yprime=dy_r(i))
   ! enddo

   ! ! Retrieve the lebedev grid and evaluate Ylm and dYlm/d(theta,phi)
   ! allocate(ylm(nleb))
   ! allocate(dylm(nleb, 2))
   ! ylm = 0._dp; dylm = 0._dp

   ! do i=1,nleb
   !    rr => lebedev_grid(ileb)%r(:, i)
   !    call rry_lm(l=l, m=m, r=rr, y=ylm(i))

   !    theta = acos(rr(3))  ! theta = acos(z/r)
   !    phi = atan2(rr(2), rr(1))  ! phi = atan2(y, x)
   !    call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))
   ! enddo

   ! ! Build the one-center grid and integrate
   ! allocate( gg(nshell*nleb,3) )
   ! allocate( gw(nshell*nleb) )
   ! do i=1,nleb
   !    lower = 1+(i-1)*nshell
   !    upper = i*nshell

   !    do j=1,nshell
   !       ! print *, (i-1)*nshell+j, r_loc(j), dy_r(j), ylm(i), dy_r(j) * ylm(i)
   !       gg(lower+j-1,1) = dy_r(j) * ylm(i)
   !       gg(lower+j-1,2) = y_loc(j) * dylm(i, 1)
   !       gg(lower+j-1,3) = y_loc(j) * dylm(i, 2)
   !    enddo
   !    gw(lower:upper) = 4.0_dp*pi * wr_loc * lebedev_grid(ileb)%w(i)
   ! enddo

   ! !
   ! !            |  df/dr * Ylm      |
   ! ! grad Psi = |  f * dYlm/dtheta  |
   ! !            |  f * dYlm/dphi    |
   ! !

   ! grad = 0._dp
   ! grad(1) = sum( gw*gg(:, 1) )
   ! grad(2) = sum( gw*gg(:, 2) )
   ! grad(3) = sum( gw*gg(:, 3) )

   ! deallocate(ylm)
   ! deallocate(dylm)
   ! deallocate(gg)
   ! deallocate(gw)
end subroutine grad_onecenter

function jacobian(r, theta, phi)
   implicit none
   REAL(KIND=dp), intent(in) :: r, theta, phi
   REAL(KIND=dp), DIMENSION(3,3) :: jacobian
   REAL(KIND=dp) :: st, ct, sp, cp
   ! jacobian(1,:) = 1st column     1  1  1
   ! jacobian(2,:) = 2nd column     2  2  2
   ! jacobian(3,:) = 3rd column     3  3  3
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