module ni_gradients
USE ni_types, ONLY: dp, pi, type_grid, type_fun, ni_env
USE lebedev, ONLY: get_number_of_lebedev_grid, lebedev_grid
USE ni_module, ONLY: spline, interpolation, integration_twocenter,&
                kinetic_energy, coulomb_integral,&
                kah_sum
USE ni_grid, ONLY: radial_grid, build_onecenter_grid, build_twocenter_grid,&
                   type_grid, deallocate_grid
USE spherical_harmonics, ONLY: rry_lm, dry_lm
USE ni_fun, ONLY: derivatives

implicit none
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ex = (/ 1._dp, 0._dp, 0._dp /)
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ey = (/ 0._dp, 1._dp, 0._dp /)
REAL(KIND=dp), DIMENSION(3), PARAMETER :: ez = (/ 0._dp, 0._dp, 1._dp /)
contains

! **********************************************
!> \brief Calculates the gradient of the coulomb integral
!> \param nshell: number of radial points to integrate the two-center integral on
!> \param coul_n: number of points to evaluate the Coulomb potential on
!> \param d12: distance vector between both centers
!> \param l: angular quantum number
!> \param m: magnetic quantum number
!> \param r1, y1, s1: the first function and spline
!> \param r2, y2, s2: the second function and spline
!> \param grad: output vector
! **********************************************
subroutine grad_coulomb(nshell, coul_n, d12, l, m,&
                        r1, y1, r2, y2, s1, s2, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: nshell, l, m
   INTEGER, intent(in) :: coul_n
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2, s1, s2
   ! Output
   REAL(KIND=dp), DIMENSION(3), intent(out) :: grad
   ! Local variables
   INTEGER :: i, j
   ! Local variables (integration)
   TYPE(type_grid), TARGET :: grid
   TYPE(type_grid), POINTER :: pgrid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: tmp_grad
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: dy2, dy2s
   REAL(KIND=dp), DIMENSION(3) :: dylm2, dr, dtheta, dphi
   REAL(KIND=dp) :: ylm, norm, x, y, z, rho, theta, phi, dylm, df2, f1, f2
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid
   ! Local variables (potential)
   REAL(KIND=dp), DIMENSION(coul_n) :: f, gi, hi, G, H, coul_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: coul_r, pot, pots

   ! 1: Evaluate the Coulomb potential on a radial grid around A, f1
   ! ! the integral is purely radial and ≠F(R) => addr2=False
   allocate(coul_r(coul_n))
   allocate(pot(coul_n))
   allocate(pots(coul_n))
   call radial_grid(r=coul_r, &
                    wr=coul_w, &
                    n=coul_n, addr2=.FALSE., quadr=1)

   do i=1,coul_n
      call interpolation(r1, y1, s1, coul_r(i), f(i))
   enddo
   gi = coul_w * coul_r**(l(1)+2) * f
   hi = coul_w * coul_r**(1-l(1)) * f

   G(coul_n) = sum(gi)
   H(coul_n) = 0.0_dp
   do j=coul_n,1,-1
      pot(j) = coul_r(j)**(-l(1)-1) * G(j) + coul_r(j)**l(1) * H(j)
      G(j-1) = G(j) - gi(j)
      H(j-1) = H(j) + hi(j)
   enddo

   pot = pot * 4.0_dp*pi/(2.0_dp*l(1)+1.0_dp)  ! Here we have V_C
   call spline(coul_r, pot, coul_n, pots)      ! and its spline.

   ! 2: Calculate the overlap of the coulomb potential and D(y2)
   !! Build the two-center grid and integrate over it
   ileb(1) = get_number_of_lebedev_grid(n=302)
   ileb(2) = get_number_of_lebedev_grid(n=302)

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(tmp_grad(ngrid, 3))

   pgrid => grid
   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.FALSE., grid=pgrid)

   allocate(dy2(size(r2)))
   allocate(dy2s(size(r2)))
   dy2 = 0._dp; dy2s = 0._dp;
   call derivatives(r=r2, y=y2, y1=dy2)
   call spline(r2, dy2, size(r2), dy2s)

   do i=1,ngrid
      ! d CI = (dw * f2*Y2 + w * df2 * Y2 + w*f2 * dY2) * Vc*Y1 = d2 * Vc*Y1
      !! we start with Vc*Y1
      norm = sqrt(sum( grid%r(i, :)**2 ))
      call interpolation(coul_r, pot, pots, norm, f1)
      call rry_lm(l=l(1), m=m(1), r=grid%r(i, :)/norm, y=ylm)
      f1 = f1 * ylm

      !!! Helpers
      norm = sqrt( sum( (grid%r(i, :) - d12)**2 ) )
      x = grid%r(i, 1)-d12(1); y = grid%r(i, 2)-d12(2); z = grid%r(i, 3)-d12(3);
      rho = x*x + y*y
      theta = acos(z/norm); phi = atan2(y, x);
      !!! The partial derivatives dr/dXYZ, dtheta/XYZ, dphi/XYZ
      dr = -(grid%r(i, :) - d12)/norm

      dtheta = 0._dp
      if (rho .ne. 0._dp) then
         dtheta(1) = -x*z
         dtheta(2) = -y*z
         dtheta(3) = rho 
         dtheta = dtheta / ( norm**3 * sqrt(1._dp - (z/norm)**2) )
      endif

      dphi = 0._dp
      if (rho .ne. 0._dp) then
         dphi(1) = y
         dphi(2) = -x
         dphi = dphi/rho
      endif

      !!! Lets-a-go
      call interpolation(r2, y2, s2, norm, f2)
      call interpolation(r2, dy2, dy2s, norm, df2)
      call rry_lm(l=l(2), m=m(2), r=(grid%r(i, :)-d12)/norm, y=ylm)
      call dry_lm(l=l(2), m=m(2), c=(/theta, phi/), dy=dylm2)

      !!! d2 = f2 * (dw * Y + w * dY) + w * df2 * Y
      dylm = dylm2(1)*dtheta(1) + dylm2(2)*dphi(1)
      tmp_grad(i, 1) = f2 * ( grid%dw(i, 1) * ylm + grid%w(i) * dylm )&
                       + grid%w(i) * df2*dr(1) * ylm

      dylm = dylm2(1)*dtheta(2) + dylm2(2)*dphi(2)
      tmp_grad(i, 2) = f2 * ( grid%dw(i, 2) * ylm + grid%w(i) * dylm )&
                       + grid%w(i) * df2*dr(2) * ylm

      dylm = dylm2(1)*dtheta(3) + dylm2(2)*dphi(3)
      tmp_grad(i, 3) = f2 * ( grid%dw(i, 3) * ylm + grid%w(i) * dylm )&
                       + grid%w(i) * df2*dr(3) * ylm

      ! Concluding...
      tmp_grad(i, :) = tmp_grad(i, :) * f1

      if(tmp_grad(i, 1) .gt. 10.0e10_dp) then
         print *, 'happened', i
         print *, 'df2', df2
         print *, 'df2', dr(1)
         print *, 'df2*', df2*dr(1)
         print *, 'grid%r', grid%r(i, :)
         print *, 'xyz', x, y, z
         print *, 'norm', norm
      endif
   enddo

   grad = 0._dp
   grad(1) = sum(tmp_grad(:, 1))
   grad(2) = sum(tmp_grad(:, 2))
   grad(3) = sum(tmp_grad(:, 3))

   deallocate(dy2)
   deallocate(dy2s)

   deallocate(coul_r)
   deallocate(pot)
   deallocate(pots)
   call deallocate_grid(grid=pgrid)
end subroutine grad_coulomb

! **********************************************
!> \brief Calculates the gradient of the coulomb integral by finite differences
! **********************************************
subroutine grad_coulomb_fd(r1, y1, r2, y2, l, m, nshell, d12, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb, nleb
   INTEGER :: ngrid, coul_n
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3, 4) :: findiff
   REAL(KIND=dp), DIMENSION(3) :: d12t
   REAL(KIND=dp) :: h

   ileb(1) = get_number_of_lebedev_grid(n=590)
   ileb(2) = get_number_of_lebedev_grid(n=590)
   nleb(1) = lebedev_grid(ileb(1))%n
   nleb(2) = lebedev_grid(ileb(2))%n
   ngrid = sum(nleb*nshell)

   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)

   h = 1.0e-7_dp
   coul_n = 1000

   d12t = d12 + 2._dp*ex*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(1,1))

   d12t = d12 + h*ex
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(1,2))
   d12t = d12 - h*ex
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(1,3))
   d12t = d12 - 2._dp*ex*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(1,4))


   d12t = d12 + 2._dp*ey*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(2,1))

   d12t = d12 + h*ey
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(2,2))
   d12t = d12 - h*ey
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(2,3))
   d12t = d12 - 2._dp*ey*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(2,4))


   d12t = d12 + 2._dp*ez*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(3,1))

   d12t = d12 + h*ez
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(3,2))
   d12t = d12 - h*ez
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(3,3))
   d12t = d12 - 2._dp*ez*h
   call coulomb_integral(l=l, m=m, nshell=nshell, coul_n=coul_n, d12=d12t, &
                         r1=r1, y1=y1, r2=r2, y2=y2,&
                         s1=s1, s2=s2, integral=findiff(3,4))


   grad(1) = -findiff(1,1) + 8._dp*findiff(1,2) - 8._dp*findiff(1,3) + findiff(1,4)
   grad(2) = -findiff(2,1) + 8._dp*findiff(2,2) - 8._dp*findiff(2,3) + findiff(2,4)
   grad(3) = -findiff(3,1) + 8._dp*findiff(3,2) - 8._dp*findiff(3,3) + findiff(3,4)
   grad = grad/(12._dp*h)
end subroutine grad_coulomb_fd


! **********************************************
!> \brief Calculates the gradient of the kinetic energy integral
! **********************************************
subroutine grad_kinetic(r1, y1, r2, y2, d1y, d2y, d3y, l, m, nshell, d12, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2, d1y, d2y, d3y
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: i, ngrid
   TYPE(type_grid), TARGET :: grid
   TYPE(type_grid), POINTER :: pgrid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: tmp_grad
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2, ddf2_spline
   REAL(KIND=dp), DIMENSION(3) :: dr, dtheta, dphi, dylm2
   REAL(KIND=dp) :: norm, x, y, z, theta, phi, rho,&
                    f1, f2, ylm1, ylm2, df2, ddf2, dddf2, ll1

   ! Newer, better derivatives
   REAL(KIND=dp), DIMENSION(size(d1y)) :: sd1y, sd2y, sd3y
   REAL(KIND=dp) :: dylm, flap, dflap
   ! call derivatives(r=r2, y=y2, y1=d1y, y2=d2y, y3=d3y)
   call spline(r=r2, y=d1y, n=size(r2), yspline=sd1y)
   call spline(r=r2, y=d2y, n=size(r2), yspline=sd2y)
   call spline(r=r2, y=d3y, n=size(r2), yspline=sd3y)

   ileb(1) = get_number_of_lebedev_grid(n=590)
   ileb(2) = get_number_of_lebedev_grid(n=590)

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(tmp_grad(ngrid, 3))

   pgrid => grid
   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.FALSE., grid=pgrid)
   !  N                        ~      ~
   !  Σ  w  * f(r) * Y(r)  Δ f(r) * Y(r) ,
   ! i=1  i    1 i    1 i     2 i    2 i

   !       ~    ––>     ––>                  T
   ! where r = grid%r - d12 = (x-X, y-Y, z-Z)

   call spline(r=r2, y=s2, n=size(r2), yspline=ddf2_spline)

   ll1 = REAL(l(2)*(l(2)+1), dp)

   grad = 0._dp
   tmp_grad = 0._dp

   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)
   do i=1,ngrid
      if(mod(i, 10000) .eq. 0) print *, 'grad_kinetic i/ngrid', i, ngrid
      ! When taking the derivative wrt X, Y, Z - f1 and ylm1 won't change
      norm = sqrt( sum( grid%r(i, :)**2 ) )
      call interpolation(gr=r1, gy=y1, spline=s1, r=norm, y=f1)
      call rry_lm(l=l(1), m=m(1), r=grid%r(i, :)/norm, y=ylm1)

      ! f2, ylm2 do have non-vanishing derivatives + are affected by the Laplacian
      !! Helpers
      norm = sqrt( sum( (grid%r(i, :) - d12)**2 ) )
      x = grid%r(i, 1)-d12(1); y = grid%r(i, 2)-d12(2); z = grid%r(i, 3)-d12(3);
      rho = x*x + y*y; theta = acos(z/norm); phi = atan2(y, x);

      !! Derivatives and function values
      !!! f2, df2, ddf2, dddf2
      call interpolation(gr=r2, gy=y2, spline=s2, r=norm, y=f2)
      call interpolation(gr=r2, gy=d1y, spline=sd1y, r=norm, y=df2)
      call interpolation(gr=r2, gy=d2y, spline=sd2y, r=norm, y=ddf2)
      call interpolation(gr=r2, gy=d3y, spline=sd3y, r=norm, y=dddf2)
      !!! ylm2, dylm2
      call rry_lm(l=l(2), m=m(2), r=(grid%r(i, :) - d12)/norm, y=ylm2)
      call dry_lm(l=l(2), m=m(2), c=(/theta, phi/), dy=dylm2)
      if(sum(dylm2**2) .ne. 0._dp) print *, dylm2

      ! The partial derivatives dr/dXYZ, dtheta/XYZ, dphi/XYZ
      dr = (grid%r(i, :) - d12)/norm

      dtheta = 0._dp
      if (rho .ne. 0._dp) then
         dtheta(1) = -x*z
         dtheta(2) = -y*z
         dtheta(3) = rho 
         dtheta = dtheta / ( norm**3 * sqrt(1._dp - (z/norm)**2) )
      endif

      dphi = 0._dp
      if (rho .ne. 0._dp) then
         dphi(1) = y
         dphi(2) = -x
         dphi = dphi/rho
      endif

      ! On to the actual calculation
      ! X

      flap = ddf2 + 2._dp/norm * df2 - ll1/norm**2 * f2
      dflap =  -dddf2/norm - 2._dp/norm**2 * ddf2&
               + (ll1 + 2._dp)/norm**3 * df2 - (2._dp * ll1)/norm**4 * f2


      dylm = dylm2(1)*dtheta(1) + dylm2(2)*dphi(1)
      tmp_grad(i, 1) = grid%dw(i, 1) * flap     * ylm2 +&
                       grid%w(i)     * x*dflap  * ylm2 +&
                       grid%w(i)     * flap     * dylm   

      ! Y
      dylm = dylm2(1)*dtheta(2) + dylm2(2)*dphi(2)
      tmp_grad(i, 2) = grid%dw(i, 2) * flap     * ylm2 +&
                       grid%w(i)     * y*dflap  * ylm2 +&
                       grid%w(i)     * flap     * dylm   

      ! Z
      dylm = dylm2(1)*dtheta(3)
      tmp_grad(i, 3) = grid%dw(i, 3) * flap     * ylm2 +&
                       grid%w(i)     * z*dflap  * ylm2 +&
                       grid%w(i)     * flap     * dylm

      ! if (.FALSE. .and. i .eq. 60306) then
      !    print *, i, norm
      !    print *, 'd12', d12
      !    print *, 'norm1', sqrt( sum( grid%r(i, :)**2 ) )
      !    print *, 'norm2', norm
      !    print *, 'grid%r ', grid%r(i, :)
      !    print *, 'x, y, z', x, y, z
      !    print *, 'theta, phi, ll1', theta, phi, ll1
      !    print *, 
      !    print *, 'f2', f2
      !    print *, 'df2', df2
      !    print *, 'ddf2', ddf2
      !    print *, 'dddf2', dddf2
      !    print *, 'flap', flap
      !    print *, 'x*dflap', x*dflap
      !    print *,
      !    print *, 'ylm1, ylm2', ylm1, ylm2
      !    print *, '–---------- tmp_grad(i,3) –----------'
      !    print *, 'grid%w(i)', grid%w(i)
      !    print *,
      !    print *, 'dr(1),', dr(1)
      !    print *, '(dddf2 + 2._dp/norm * ddf2 - ll1/norm**2 * df2)', (dddf2 + 2._dp/norm * ddf2 - ll1/norm**2 * df2)
      !    print *, '2._dp*z*( df2/norm**3 + ll1*f2 )', 2._dp*z*( df2/norm**3 + ll1*f2 )
      !    print *, 
      !    print *, 'grid%dw(i, 3)', grid%dw(i, 3)
      !    print *,
      !    print *, 'dylm', dylm
      !    print *, 'dylm2(1)*dtheta(1)', dylm2(1)*dtheta(1)
      !    print *, 'dylm2', dylm2
      !    print *, 'dtheta(1)', dtheta(1)
      !    print *, 
      !    print *, 'tmp_grad(i, 1)', tmp_grad(i, 1)
      ! endif

      ! Concluding...
      tmp_grad(i, :) = -0.5_dp * f1 * ylm1 * tmp_grad(i, :) 
   enddo

   grad = 0._dp
   grad(1) = kah_sum(tmp_grad(:,1))
   grad(2) = kah_sum(tmp_grad(:,2))
   grad(3) = kah_sum(tmp_grad(:,3))

   deallocate(tmp_grad)
   call deallocate_grid(grid=pgrid)
end subroutine grad_kinetic

! **********************************************
!> \brief Calculates the gradient of the coulomb integral by finite differences
! **********************************************
subroutine grad_kinetic_fd(r1, y1, r2, y2, l, m, nshell, d12, step, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), intent(in), optional :: step
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb, nleb
   INTEGER :: ngrid
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3, 4) :: findiff
   REAL(KIND=dp), DIMENSION(3) :: d12t
   REAL(KIND=dp) :: h

   ileb(1) = get_number_of_lebedev_grid(n=590)
   ileb(2) = get_number_of_lebedev_grid(n=590)
   nleb(1) = lebedev_grid(ileb(1))%n
   nleb(2) = lebedev_grid(ileb(2))%n
   ngrid = sum(nleb*nshell)

   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)

   if (present(step)) then
      h = step
   else
      h = 1.0e-7_dp
   endif

   print *, '1'
   d12t = d12 + 2._dp*ex*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(1,1))

   d12t = d12 + h*ex
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(1,2))
   d12t = d12 - h*ex
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(1,3))
   d12t = d12 - 2._dp*ex*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(1,4))

   print *, '2'

   d12t = d12 + 2._dp*ey*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(2,1))

   d12t = d12 + h*ey
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(2,2))
   d12t = d12 - h*ey
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(2,3))
   d12t = d12 - 2._dp*ey*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(2,4))

   print *, '3'

   d12t = d12 + 2._dp*ez*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(3,1))

   d12t = d12 + h*ez
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(3,2))
   d12t = d12 - h*ez
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(3,3))
   d12t = d12 - 2._dp*ez*h
   call kinetic_energy(l=l, m=m, nshell=nshell, d12=d12t, &
                       r1=r1, y1=y1, r2=r2, y2=y2,&
                       spline1=s1, spline2=s2, integral=findiff(3,4))


   grad(1) = -findiff(1,1) + 8._dp*findiff(1,2) - 8._dp*findiff(1,3) + findiff(1,4)
   grad(2) = -findiff(2,1) + 8._dp*findiff(2,2) - 8._dp*findiff(2,3) + findiff(2,4)
   grad(3) = -findiff(3,1) + 8._dp*findiff(3,2) - 8._dp*findiff(3,3) + findiff(3,4)
   grad = grad/(12._dp*h)
end subroutine grad_kinetic_fd

! **********************************************
!> \brief Calculates the gradient of the two-center integral
! **********************************************
subroutine grad_twocenter(r1, y1, r2, y2, l, m, nshell, d12, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid, i
   TYPE(type_grid), TARGET :: grid
   TYPE(type_grid), POINTER :: pgrid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: tmp_grad
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3) :: dr, dtheta, dphi, dylm2
   REAL(KIND=dp) :: norm, x, y, z, theta, phi, rho,&
                     f1, f2, ylm1, ylm2, df2

   ileb(1) = get_number_of_lebedev_grid(n=590)
   ileb(2) = get_number_of_lebedev_grid(n=590)

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(tmp_grad(ngrid, 3))

   pgrid => grid
   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.TRUE., grid=pgrid)
   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)

   grad = 0._dp
   tmp_grad = 0._dp
   do i=1,ngrid
      !  N                ~             ~
      !  Σ  w  * f(r) * f(r) * Y(r) * Y(r) ,
      ! i=1  i    1 i    2 i      i      i

      !       ~    ––>     ––>                  T
      ! where r = grid%r - d12 = (x-X, y-Y, z-Z)

      ! When taking the derivative wrt X, Y, Z - f1 and ylm1 won't change
      norm = sqrt( sum( grid%r(i, :)**2 ) )
      call interpolation(gr=r1, gy=y1, spline=s1, r=norm, y=f1)
      call rry_lm(l=l(1), m=m(1), r=grid%r(i, :)/norm, y=ylm1)

      ! f2, ylm2 do have non-vanishing derivatives
      norm = sqrt( sum( (grid%r(i, :) - d12)**2 ) )
      x = grid%r(i, 1)-d12(1); y = grid%r(i, 2)-d12(2); z = grid%r(i, 3)-d12(3);
      rho = x*x + y*y
      theta = acos(z/norm); phi = atan2(y, x);

      call interpolation(gr=r2, gy=y2, spline=s2, r=norm, y=f2, yprime=df2)
      call rry_lm(l=l(2), m=m(2), r=(grid%r(i, :) - d12)/norm, y=ylm2)
      call dry_lm(l=l(2), m=m(2), c=(/theta, phi/), dy=dylm2)

      ! dwdr2 = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))

      ! The partial derivatives dr/dXYZ, dtheta/XYZ, dphi/XYZ
      dr = -(grid%r(i, :) - d12)/norm

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
      ! X  
      tmp_grad(i, 1) = grid%w(i) * df2*dr(1) * ylm2 +&
                       grid%w(i) * f2 * dylm2(1)*dtheta(1) +&
                       grid%w(i) * f2 * dylm2(2)*dphi(1)! +&
                       ! grid%dw(i, 1) * f2 * ylm2

      ! Y
      tmp_grad(i, 2) = grid%w(i) * df2 * dr(2) * ylm2 +&
                       grid%w(i) * f2 * dylm2(1) * dtheta(2) +&
                       grid%w(i) * f2 * dylm2(2) * dphi(2)! +&
                       ! grid%dw(i, 2) * f2 * ylm2

      ! Z
      tmp_grad(i, 3) = grid%w(i) * df2 * dr(3) * ylm2 +&
                       grid%w(i) * f2 * dylm2(1) * dtheta(3)! +&
                       ! grid%dw(i, 3) * f2 * ylm2

      tmp_grad(i, :) = tmp_grad(i, :) * f1 * ylm1
   enddo

   grad = 0._dp
   grad(1) = kah_sum(tmp_grad(:,1))
   grad(2) = kah_sum(tmp_grad(:,2))
   grad(3) = kah_sum(tmp_grad(:,3))

   deallocate(tmp_grad)
   call deallocate_grid(grid=pgrid)
end subroutine grad_twocenter

! **********************************************
!> \brief Calculates the gradient of the two-center integral by finite differences
! **********************************************
subroutine grad_twocenter_fd(r1, y1, r2, y2, l, m, nshell, d12, step, grad)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), intent(in), optional :: step
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   INTEGER, DIMENSION(2) :: ileb, nleb
   INTEGER :: ngrid
   REAL(KIND=dp), DIMENSION(size(r1)) :: s1
   REAL(KIND=dp), DIMENSION(size(r2)) :: s2
   REAL(KIND=dp), DIMENSION(3, 4) :: findiff
   REAL(KIND=dp), DIMENSION(3) :: d12t
   REAL(KIND=dp) :: h

   ileb(1) = get_number_of_lebedev_grid(n=302)
   ileb(2) = get_number_of_lebedev_grid(n=302)
   nleb(1) = lebedev_grid(ileb(1))%n
   nleb(2) = lebedev_grid(ileb(2))%n
   ngrid = sum(nleb*nshell)

   call spline(r=r1, y=y1, n=size(r1), yspline=s1)
   call spline(r=r2, y=y2, n=size(r2), yspline=s2)

   h = 6.5e-7_dp
   if (present(step)) h = step

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
end subroutine grad_twocenter_fd

! **********************************************
!> \brief Calculates the gradient of the one-center integral
!>        Should be zero.
! **********************************************
subroutine grad_onecenter(r, y, l, m, nshell, grad)
   implicit none
   ! Input
   INTEGER, intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, y
   ! Output
   REAL(KIND=dp), DIMENSION(3) :: grad
   ! Local variables
   TYPE(type_grid), POINTER :: grid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: dylm, tmp_grad
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y_loc, dy, ylm
   REAL(KIND=dp), DIMENSION(size(r)) :: s
   REAL(KIND=dp), DIMENSION(3,3) :: jac
   REAL(KIND=dp) :: theta, phi, xx, yy, zz, dwdr, alpha, norm
   INTEGER :: i, ileb, nleb


   ! Get the one-center grid
   ileb = get_number_of_lebedev_grid(l=l+5)
   nleb = lebedev_grid(ileb)%n


   allocate(y_loc(nleb*nshell))
   allocate(dy(nleb*nshell))

   allocate(ylm(nleb*nshell))
   allocate(dylm(nleb*nshell, 2))

   allocate(tmp_grad(size(grid%w),3))

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE.,&
                             quadr=1, grid=grid)
   alpha = pi/(REAL(2*nshell+2, dp))
   ! Evaluate y, dy/dr, Ylm, dYlm/dtheta, dYlm/dphi
   call spline(r=r, y=y, n=size(r), yspline=s)
   tmp_grad = 0._dp
   do i=1,size(grid%w)
      ! gradient (r, theta, phi)

      ! Define some helpers
      norm = sqrt(sum( grid%r(i, :)**2 ))
      xx = grid%r(i, 1); yy = grid%r(i, 2); zz = grid%r(i, 3)
      theta = acos(zz/norm)
      phi = atan2(yy, xx)

      ! We will need all those
      ! f, df/dr, dw/dr
      call interpolation(gr=r, gy=y, spline=s, r=norm+5._dp*epsilon(1._dp), y=y_loc(i), yprime=dy(i))
      dwdr = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))
      ! Ylm, dYlm/dtheta, dYlm/dphi
      call rry_lm(l=l, m=m, r=grid%r(i, :)/norm, y=ylm(i))
      call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))

      !                                                          T
      ! tmp_grad = (d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi )
      !

      ! d/dr => (dw/dr * f + w * df/dr) * Ylm
      !    dw/dr = pi/(N+1) * (7*r^2.5 + 5*r^1.5)/2 
      !    df/dr <- dy
      tmp_grad(i,1) = ( dwdr * y_loc(i) + grid%w(i) * dy(i) ) * ylm(i)

      ! 1/r d/dtheta => w * f * dYlm/dtheta
      tmp_grad(i,2) = grid%w(i) * y_loc(i) * dylm(i, 1) / norm

      ! 1/(r * sin(theta)) d/dphi => w * f * dYlm/dphi / r*sin(theta)
      if (dylm(i, 2) .ne. 0._dp) then
         tmp_grad(i,3) = grid%w(i) * y_loc(i) * dylm(i, 2) / (norm * sin(theta))
      endif

      ! And transform it to xyz
      jac = jacobian(r=norm, theta=theta, phi=phi)
      tmp_grad(i,:) = matmul(jac, tmp_grad(i,:))
   enddo

   grad = 0._dp
   grad(1) = sum(tmp_grad(:,1))
   grad(2) = sum(tmp_grad(:,2))
   grad(3) = sum(tmp_grad(:,3))

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
   TYPE(type_grid), POINTER :: grid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: dylm, tmp_grad
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y_loc, dy, ylm
   REAL(KIND=dp), DIMENSION(size(r)) :: s
   REAL(KIND=dp) :: theta, phi, xx, yy, zz, dwdr, alpha, norm, rho
   INTEGER :: i, ileb, nleb

   ! Get the one-center grid
   ileb = get_number_of_lebedev_grid(l=l+5)
   nleb = lebedev_grid(ileb)%n


   allocate(y_loc(nleb*nshell))
   allocate(dy(nleb*nshell))

   allocate(ylm(nleb*nshell))
   allocate(dylm(nleb*nshell, 2))

   allocate(tmp_grad(size(grid%w),3))

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE.,&
                             quadr=1, grid=grid)
   alpha = pi/(REAL(2*nshell+2, dp))
   ! Evaluate y, dy/dr, Ylm, dYlm/dtheta, dYlm/dphi
   call spline(r=r, y=y, n=size(r), yspline=s)
   tmp_grad = 0._dp
   do i=1,size(grid%w)
      ! gradient (r, theta, phi)

      ! Define some helpers
      norm = sqrt(sum( grid%r(i, :)**2 ))
      xx = grid%r(i, 1); yy = grid%r(i, 2); zz = grid%r(i, 3)
      rho = xx**2 + yy**2
      theta = acos(zz/norm)
      phi = atan2(yy, xx)

      ! We will need all those
      ! f, df/dr, dw/dr
      call interpolation(gr=r, gy=y, spline=s, r=norm+5._dp*epsilon(1._dp), y=y_loc(i), yprime=dy(i))
      dwdr = alpha*(7._dp*norm**(2.5_dp) + 5._dp*norm**(1.5_dp))
      ! Ylm, dYlm/dtheta, dYlm/dphi
      call rry_lm(l=l, m=m, r=grid%r(i, :)/norm, y=ylm(i))
      call dry_lm(l=l, m=m, c=(/theta, phi/), dy=dylm(i, :))

      !                                                          T
      ! tmp_grad = (d/dr, 1/r d/dtheta, 1/(r*sin(theta)) d/dphi )
      !

      ! dr/dx = x/r, dr/dy = y/r, dr/dz = z/r
      tmp_grad(i,:) = ( dwdr * y_loc(i) + grid%w(i) * dy(i) ) * ylm(i)
      tmp_grad(i,1) = tmp_grad(i,1) * xx/norm
      tmp_grad(i,2) = tmp_grad(i,2) * yy/norm
      tmp_grad(i,3) = tmp_grad(i,3) * zz/norm

      ! dtheta/dx = x*z/(r^3*sqrt(1-z^2/r^2))
      ! dtheta/dy = y*z/(r^3*sqrt(1-z^2/r^2))
      ! dtheta/dz = -sqrt(rho^2/r^2)/r
      if (xx*zz .ne. 0._dp) then
      tmp_grad(i,1) = tmp_grad(i,1) + grid%w(i) * y_loc(i) * dylm(i, 1)&
         * xx*zz/( norm**3 * sqrt( 1._dp-(zz/norm)**2 ) )
      endif
      if (yy*zz .ne. 0._dp) then
      tmp_grad(i,2) = tmp_grad(i,2) + grid%w(i) * y_loc(i) * dylm(i, 1)&
         * yy*zz/( norm**3 * sqrt( 1._dp-(zz/norm)**2 ) )
      endif
      tmp_grad(i,3) = tmp_grad(i,3) + grid%w(i) * y_loc(i) * dylm(i, 1)&
         * sqrt((xx*xx+yy*yy)/norm**2)/norm

      ! dphi/dx = -y/(x^2+y^2)
      ! dphi/dy = x/(x^2+y^2)
      ! dphi/dz = 0
      if (yy .ne. 0._dp) then
         tmp_grad(i,1) = tmp_grad(i,1) - grid%w(i) * y_loc(i) * dylm(i, 2) * (yy/rho)
      endif
      if (xx .ne. 0._dp) then
         tmp_grad(i,2) = tmp_grad(i,2) + grid%w(i) * y_loc(i) * dylm(i, 2) * (xx/rho)
      endif
   enddo

   grad = 0._dp
   grad(1) = sum(tmp_grad(:,1))
   grad(2) = sum(tmp_grad(:,2))
   grad(3) = sum(tmp_grad(:,3))

   deallocate(y_loc)
   deallocate(dy)
   deallocate(ylm)
   deallocate(dylm)
   deallocate(tmp_grad)
end subroutine grad_onecenter_cart

! **********************************************
!> \brief Calculates the jacobian of (r, theta, phi)
! **********************************************
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

end module ni_gradients
