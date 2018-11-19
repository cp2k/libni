module grid
   USE lebedev, ONLY: dp, lebedev_grid
   implicit none
   type :: type_grid_point
      REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
      REAL(KIND=dp) :: weight = 0.0_dp
   end type type_grid_point
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   public :: type_grid_point, build_onecenter_grid, build_twocenter_grid, &
             build_threecenter_grid, radial_grid
    
contains

! Writes a radial grid of size N into r and the weights in to wr
subroutine radial_grid(r, wr, n, addr2, quadr)
   implicit none
   ! Input
   INTEGER, intent(in) :: n
   LOGICAL, OPTIONAL, intent(in) :: addr2
   !! 1: Gauss-Chebyshev
   !! 2: Gauss-Hermite
   INTEGER, intent(in) :: quadr
   ! Output
   REAL(KIND=dp), DIMENSION(:) :: r, wr
   ! Local variables
   INTEGER :: i
   REAL(KIND=dp) :: alpha, t, x
   REAL(KIND=dp), DIMENSION(n) :: her_r, her_wr

   alpha = pi/REAL(n+1, dp)
   if (quadr .eq. 1) then
      do i=1,n
         t = REAL(i, dp)*alpha
         x = COS(t)
         r(i) = (1.0_dp+x)/(1.0_dp-x)
         wr(i) = alpha*2.0_dp*SIN(t)/(1.0_dp-x)**2
      enddo
   else if (quadr .eq. 2) then
      call gauher(r=her_r, wr=her_wr, n=n)
      r = EXP(her_r)
      wr = her_wr * r
   else
      stop 'quadr'
   endif

   if (present(addr2) .and. (addr2 .eqv. .true.)) then
      ! dxdydz = dr r^2 dcos(theta) dphi
      wr = wr * r**2
   endif
   r = r(n:1:-1)
   wr = wr(n:1:-1)
end subroutine radial_grid

subroutine build_onecenter_grid(ileb, nshell, addr2, quadr, grid_r, grid_w)
   implicit none
   ! Input
   INTEGER, intent(in) :: ileb, nshell, quadr
   LOGICAL, OPTIONAL, intent(in) :: addr2
   ! Output
   REAL(KIND=dp), DIMENSION(:, :) :: grid_r
   REAL(KIND=dp), DIMENSION(:) :: grid_w
   ! Local variables
   INTEGER :: i, j, lower, upper
   REAL(KIND=dp), DIMENSION(nshell) :: radii, radii_w
   LOGICAL :: aa

   aa = .FALSE.
   if (present(addr2)) aa = addr2
   call radial_grid(r=radii, &
                    wr=radii_w, &
                    n=nshell, addr2=aa, quadr=quadr)

   do i=1, lebedev_grid(ileb)%n
      lower = 1+(i-1)*nshell
      upper = i*nshell

      do j=1,nshell
         grid_r(lower+j-1,:) = radii(j) * lebedev_grid(ileb)%r(:, i)
      enddo
      grid_w(lower:upper) = 4.0_dp*pi * radii_w * lebedev_grid(ileb)%w(i)
   enddo
end subroutine build_onecenter_grid


subroutine build_twocenter_grid(ileb, nshell, d12, addr2, grid_r, grid_w, grid_dw)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   LOGICAL, OPTIONAL, intent(in) :: addr2
   ! Output
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: grid_dw
   ! Local variables
   REAL(KIND=dp), DIMENSION(nshell(1)) :: radii1, radii_w1
   REAL(KIND=dp), DIMENSION(nshell(2)) :: radii2, radii_w2
   REAL(KIND=dp), DIMENSION(2) :: alpha
   INTEGER :: i, j, lower, upper, offset
   REAL(KIND=dp) :: R, r1, r2, mu, s1, s2, norm

   alpha(1) = pi/REAL(nshell(1)+1, dp)
   alpha(2) = pi/REAL(nshell(2)+1, dp)
   R = sqrt( sum( d12**2 ) )
   grid_w = 0._dp
   if (R .eq. 0.0_dp) then
      call build_onecenter_grid(ileb=ileb(1), nshell=nshell(1), addr2=addr2,&
                                grid_r=grid_r, grid_w=grid_w, quadr=1)
      return
   endif

   call radial_grid(r=radii1, &
                    wr=radii_w1, &
                    n=nshell(1), addr2=addr2, quadr=1)

   call radial_grid(r=radii2, &
                    wr=radii_w2, &
                    n=nshell(2), addr2=addr2, quadr=1)

   do i=1, lebedev_grid(ileb(1))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = 1+(i-1)*nshell(1)
      upper = i*nshell(1)

      do j=1,nshell(1)
         grid_r(lower+j-1, :) = radii1(j) * lebedev_grid(ileb(1))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w1 * lebedev_grid(ileb(1))%w(i)
   enddo

   offset = lebedev_grid(ileb(1))%n * nshell(1)
   do i=1, lebedev_grid(ileb(2))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(2)
      upper = offset + i*nshell(2)

      ! The second batch of grid points is displaced by `d12`
      do j=1,nshell(2)
         grid_r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i) + d12
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w2 * lebedev_grid(ileb(2))%w(i)
   enddo

   ! nuclear partition
   !! A = 0, B = d12
   !! r1 is the distance A -> grid_r = grid_r-0 = grid_r
   !! r2 is the distance B -> grid_r = grid_r-d12
   do i=1,offset
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu)
      s2 = s3(-mu)

      if (abs(s1+s2-1._dp) .ge. epsilon(1._dp)) stop 'Two-center nuclear partition'

      ! dw/dX = w_rad * ds1/dX
      if (present(grid_dw)) then
         ! dw/dX = w_rad * ds_3/dmu * dmu/dX
         grid_dw(i, 1) = (d12(1)*(r2-r1)/R**3 + (grid_r(i, 1)-d12(1))/(R*r2) )
         grid_dw(i, 2) = (d12(2)*(r2-r1)/R**3 + (grid_r(i, 2)-d12(2))/(R*r2) )
         grid_dw(i, 3) = (d12(3)*(r2-r1)/R**3 + (grid_r(i, 3)-d12(3))/(R*r2) )
         grid_dw(i, :) = grid_dw(i, :) * grid_w(i) * ds3dmu(mu)
      endif

      grid_w(i) = grid_w(i) * s1!/(s1+s2)
   enddo

   do i=1+offset,size(grid_w)
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu)
      s2 = s3(-mu)

      if (abs(s1+s2-1._dp) .ge. epsilon(1._dp)) stop 'Two-center nuclear partition'

      if (present(grid_dw)) then
         ! dw/dX = (1)w_rad * ds_3/dmu*dmu/dX + (2) dw_rad/dr * dr/dX * wpart
         ! (1): w_rad * ds_3/dmu*dmu/dX
         grid_dw(i, 1) = (d12(1)*(r2-r1)/R**3 + (grid_r(i, 1)-d12(1))/(R*r2) )
         grid_dw(i, 2) = (d12(2)*(r2-r1)/R**3 + (grid_r(i, 2)-d12(2))/(R*r2) )
         grid_dw(i, 3) = (d12(3)*(r2-r1)/R**3 + (grid_r(i, 3)-d12(3))/(R*r2) )
         grid_dw(i, :) = -grid_dw(i, :) * grid_w(i) * ds3dmu(-mu)

         ! TODO
         ! ! (2): dw_rad/dr * dr/dX * wpart
         ! grid_dw(i, 1) = grid_dw(i, 1)&
         !    + alpha(2)*(7._dp*r2**(2.5_dp) + 5._dp*r2**(1.5_dp)) &
         !    * (d12(1) - grid_r(i, 1))/r2 * s2
      endif
      
      grid_w(i) = grid_w(i) * s2!/(s1+s2)
   enddo
end subroutine build_twocenter_grid

subroutine build_threecenter_grid(ileb, nshell, d12, d13, addr2, grid_r, grid_w)
   implicit none
   ! Input
   INTEGER, DIMENSION(3), intent(in) :: ileb, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
   LOGICAL, OPTIONAL, intent(in) :: addr2
   ! Output
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   ! Local variables
   REAL(KIND=dp), DIMENSION(nshell(1)) :: radii1, radii_w1
   REAL(KIND=dp), DIMENSION(nshell(2)) :: radii2, radii_w2
   REAL(KIND=dp), DIMENSION(nshell(3)) :: radii3, radii_w3
   INTEGER :: i, j, lower, upper, offset, off1, off2
   REAL(KIND=dp) :: R12, R13, R23, r1, r2, r3,&
                    mu12, mu13, mu23, s12, s13, s23, s21, s31, s32

   REAL(KIND=dp) :: tP1, tP2, tP3, sP, p1, p2, p3

   call radial_grid(r=radii1, &
                    wr=radii_w1, &
                    n=nshell(1), addr2=.TRUE., quadr=1)

   call radial_grid(r=radii2, &
                    wr=radii_w2, &
                    n=nshell(2), addr2=.TRUE., quadr=1)

   call radial_grid(r=radii3, &
                    wr=radii_w3, &
                    n=nshell(3), addr2=.TRUE., quadr=1)

   R12 = sqrt(sum(d12**2))
   R13 = sqrt(sum(d13**2))
   R23 = sqrt(sum((d13-d12)**2))

   ! Center 1
   do i=1, lebedev_grid(ileb(1))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = 1+(i-1)*nshell(1)
      upper = i*nshell(1)

      do j=1,nshell(1)
          grid_r(lower+j-1,:) = radii1(j) * lebedev_grid(ileb(1))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w1 * lebedev_grid(ileb(1))%w(i)
   enddo

   ! Center 2
   offset = lebedev_grid(ileb(1))%n * nshell(1)
   do i=1, lebedev_grid(ileb(2))%n
       ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(2)
      upper = offset + i*nshell(2)

      do j=1,nshell(2)
          grid_r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i) + d12
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w2 * lebedev_grid(ileb(2))%w(i)
   enddo

   ! Center 3
   offset = offset + lebedev_grid(ileb(2))%n * nshell(2)
   do i=1, lebedev_grid(ileb(3))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(3)
      upper = offset + i*nshell(3)

      do j=1,nshell(3)
          grid_r(lower+j-1,:) = radii3(j) * lebedev_grid(ileb(3))%r(:, i) + d13
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w3 * lebedev_grid(ileb(3))%w(i)
   enddo

   ! nuclear partition
   !! A = 0, B = d12, C = d13
   !! r1 is the distance A -> grid_r = grid_r-0 = grid_r
   !! r2 is the distance B -> grid_r = grid_r-d12
   !! r3 is the distance C -> grid_r = grid_r-d13

   ! Lets us jump to the right index in grid_r, grid_w
   off1 = lebedev_grid(ileb(1))%n * nshell(1)
   off2 = off1 + lebedev_grid(ileb(2))%n * nshell(2)
   do i=1,off1
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid_r(i, :) - d13)**2 ))

      mu12 = (r1-r2)/R12; mu13 = (r1-r3)/R13; mu23 = (r2-r3)/R23

      s12 = s3(mu12)
      s21 = s3(-mu12)
      s13 = s3(mu13)
      s31 = s3(-mu13)
      s32 = s3(-mu23)
      s23 = s3(mu23)

      tP1 = s12*s13; tP2 = s21*s23; tP3 = s31*s32
      sP = tP1+tP2+tP3
      p1 = tP1/sP; p2 = tP2/sP; p3 = tP3/sP

      grid_w(i) = grid_w(i) * p1

      if (abs(mu12).gt.1.0_dp .or. abs(mu13).gt.1.0_dp .or. abs(mu23).gt.1.0_dp) then
         stop 'Nuclear partition - Three-center'
      endif
      ! if (abs(sP-1.0_dp) .gt. 0.1_dp) then
      !    print *, i
      !    print *, s3(-1.0_dp), s3(0.0_dp), s3(1.0_dp)
      !    print *, sP, grid_w(i), grid_r(i, :)
      ! endif
   enddo

   do i=off1+1,off2
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid_r(i, :) - d13)**2 ))

      mu12 = (r1-r2)/R12; mu13 = (r1-r3)/R13; mu23 = (r2-r3)/R23

      s12 = s3(mu12)
      s21 = s3(-mu12)
      s13 = s3(mu13)
      s31 = s3(-mu13)
      s32 = s3(-mu23)
      s23 = s3(mu23)

      tP1 = s12*s13; tP2 = s21*s23; tP3 = s31*s32
      sP = tP1+tP2+tP3
      p1 = tP1/sP; p2 = tP2/sP; p3 = tP3/sP
      
      grid_w(i) = grid_w(i) * p2
      if (abs(mu12).gt.1.0_dp .or. abs(mu13).gt.1.0_dp .or. abs(mu23).gt.1.0_dp) then
         stop 'Nuclear partition - Three-center'
      endif
   enddo

   do i=off2+1,size(grid_w)
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid_r(i, :) - d13)**2 ))

      mu12 = (r1-r2)/R12; mu13 = (r1-r3)/R13; mu23 = (r2-r3)/R23

      s12 = s3(mu12)
      s21 = s3(-mu12)
      s13 = s3(mu13)
      s31 = s3(-mu13)
      s32 = s3(-mu23)
      s23 = s3(mu23)

      tP1 = s12*s13; tP2 = s21*s23; tP3 = s31*s32
      sP = tP1+tP2+tP3
      p1 = tP1/sP; p2 = tP2/sP; p3 = tP3/sP
      
      grid_w(i) = grid_w(i) * p3
      if (abs(mu12).gt.1.0_dp .or. abs(mu13).gt.1.0_dp .or. abs(mu23).gt.1.0_dp) then
         stop 'Nuclear partition - Three-center'
      endif
   enddo
end subroutine build_threecenter_grid

   function h(mu)
      implicit none
      REAL(KIND=dp) :: mu, h
      h = 1.5_dp*mu-0.5_dp*mu**3
   end function h

   function z(mu)
      implicit none
      REAL(KIND=dp) :: mu, mua, mua2, mua4, mua6, z
      mua = (mu/0.64_dp)
      mua2 = mua*mua
      mua4 = mua2*mua2
      mua6 = mua4*mua2
      z = mua*(35.0_dp*(1-mua2)+21.0_dp*mua4-5.0_dp*mua6)/16.0_dp
   end function z

function s3(mu)
   implicit none
   REAL(KIND=dp) :: mu, a, a3, s3, mu3
   ! s3 = 0.5_dp*(1._dp - h(h(h(mu))) )
   mu3 = mu**3
   a = mu3 - 3*mu
   a3 = a**3

   s3 = (a3 - 12._dp*mu3 + 36._dp*mu)**3/16384._dp -&
         .046875_dp*a3 + 0.5625_dp*mu3 - 1.6875_dp*mu + 0.5_dp
end function s3

function ds3dmu(mu)
   implicit none
   REAL(KIND=dp) :: mu, mu2, mu3, a, a2, a3, ds3dmu
   mu2 = mu*mu; mu3 = mu2*mu
   a = mu3 - 3._dp*mu
   a2 = a*a; a3 = a2*a

   ds3dmu = 27._dp/16384._dp * ( a3 - 12._dp*mu3 + 36._dp*mu )**2&
             * (a2*(mu2-1._dp) - 4._dp*mu2 + 4._dp) &
             - 27._dp/64._dp * a2*(mu2-1._dp) + 27._dp/16._dp * (mu2 - 1)
end function ds3dmu

   subroutine grid_parameters(atom, nleb, nshell)
      implicit none
      INTEGER :: atom, nleb, nshell

      SELECT CASE (atom)
      CASE (0)
         nshell = 5
         nleb = 10
      CASE (1:2)
         nshell = 35
         nleb = 302
      CASE (3:10)
         nshell = 40
         nleb = 590
      CASE (11:18)
         nshell = 45
         nleb = 590
      CASE (19:)
         nshell = 50
         nleb = 590
      CASE DEFAULT
         nshell = 10
         nleb = 302
      END SELECT
   end subroutine grid_parameters

! Compute the `n`th Hermite polynomial y=H_n(x).
recursive subroutine hermite(n, x, y)
   implicit none
   INTEGER, intent(in) :: n
   REAL(KIND=dp), intent(in) :: x
   REAL(KIND=dp) :: y, a, b

   if (n .eq. 0) then
      y = 1.0_dp
      return
   else if (n .eq. 1) then
      y = 2.0_dp * x
      return
   endif

   call hermite(n=n-1, x=x, y=a)
   call hermite(n=n-2, x=x, y=b)
   y = 2.0_dp * (x*a - (n-1._dp)*b)
end subroutine hermite

! Return the nodes `r` and weights `wr` of Gauss-Hermite quadratures of order n
subroutine gauher(r, wr, n)
   implicit none
   ! Input
   INTEGER :: n
   ! Output
   REAL(KIND=dp), DIMENSION(n) :: r, wr
   ! Local variables
   INTEGER :: i, j, k, m, maxit
   REAL(KIND=dp) :: u, pim4, p1, p2, p3, p4, pp, z, z1

   ! Only find n/2 zeros
   m = (n+1)/2
   do i=1,m
      ! Initial guess for the root
      if (i .eq. 1) then
         z = sqrt( REAL(2*n+1, dp) ) - 1.85575_dp * (2.0_dp*n+1)**(-0.16667_dp)
      else if (i .eq. 2) then
         z = z-1.14_dp*n**0.426_dp / z
      else if (i .eq. 3) then
         z = 1.86_dp*z - 0.86_dp*r(1)
      else if (i .eq. 4) then
         z = 1.91_dp*z - 0.91_dp*r(2)
      else
         z = 2.0_dp*z - r(i-2)
      endif

      maxit = 1e6
      pim4 = pi**(-0.25_dp)
      do k=1,maxit
         p1 = pim4
         p2 = 0.0_dp
         do j=1,n
            p3 = p2
            p2 = p1
            p1 = z*sqrt(2.0_dp/REAL(j, dp))*p2 - sqrt(REAL(j-1, dp)/REAL(j, dp))*p3
         enddo
         pp = sqrt(2.0_dp * n)*p2
         z1 = z
         z = z1 - p1/pp
         if (abs(z-z1) .le. epsilon(z)) exit
      enddo
      r(i) = z
      r(n+1-i) = -z
      wr(i) = 2.0_dp/(pp**2) * exp(r(i)**2)
      wr(n+1-i) = wr(i)
   enddo
end subroutine gauher

end module grid