module ni_grid
use ni_types, only: dp, pi, type_grid, type_fun, ni_env
use lebedev, only: lebedev_grid
implicit none
public :: type_grid, build_onecenter_grid, build_twocenter_grid, &
          build_threecenter_grid, radial_grid
    
contains

! **********************************************
!> \brief Allocate the grid
! **********************************************
subroutine allocate_grid(grid, n)
   implicit none
   type(type_grid), pointer :: grid
   integer :: n
   if (associated(grid)) then
   allocate(grid%r(n, 3))
   allocate(grid%w(n))
   allocate(grid%dw(n, 3))
   else
   stop 'Grid is not associated!'
   endif
end subroutine allocate_grid

subroutine deallocate_grid(grid)
   implicit none
   type(type_grid), pointer :: grid

   deallocate(grid%r)
   deallocate(grid%w)
   deallocate(grid%dw)
end subroutine deallocate_grid

! **********************************************
!> \brief Write a radial grid with n points into r, wr
!>        quadr = 1: Chebyshev
!>        quadr = 2: Hermite
! **********************************************
subroutine radial_grid(r, wr, n, addr2, quadr)
   implicit none
   ! Input
   integer, intent(in) :: n
   logical, OPTIONAL, intent(in) :: addr2
   !! 1: Gauss-Chebyshev
   !! 2: Gauss-Hermite
   integer, intent(in) :: quadr
   ! Output
   real(kind=dp), dimension(:) :: r, wr
   ! Local variables
   integer :: i
   real(kind=dp) :: alpha, t, x
   real(kind=dp), dimension(n) :: her_r, her_wr

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

! **********************************************
!> \brief Generates the one-center grid
! **********************************************
subroutine build_onecenter_grid(ileb, nshell, addr2, quadr, grid)
   implicit none
   ! Input
   integer, intent(in) :: ileb, nshell, quadr
   logical, OPTIONAL, intent(in) :: addr2
   ! Output
   type(type_grid), pointer :: grid
   ! Local variables
   integer :: i, j, lower, upper, ngrid
   real(kind=dp), dimension(nshell) :: radii, radii_w
   logical :: aa

   aa = .false.
   if (present(addr2)) aa = addr2
   call radial_grid(r=radii, &
                    wr=radii_w, &
                    n=nshell, addr2=aa, quadr=quadr)
   ngrid = lebedev_grid(ileb)%n * nshell
   call allocate_grid(grid=grid, n=ngrid)

   do i=1, lebedev_grid(ileb)%n
      lower = 1+(i-1)*nshell
      upper = i*nshell

      do j=1,nshell
         grid%r(lower+j-1,:) = radii(j) * lebedev_grid(ileb)%r(:, i)
      enddo
      grid%w(lower:upper) = 4.0_dp*pi * radii_w * lebedev_grid(ileb)%w(i)
   enddo
   grid%dw = 0._dp
end subroutine build_onecenter_grid

! **********************************************
!> \brief Generates the two-center grid
! **********************************************
subroutine build_twocenter_grid(ileb, nshell, d12, addr2, grid)
   implicit none
   ! Input
   integer, dimension(2), intent(in) :: ileb, nshell
   real(kind=dp), dimension(3), intent(in) :: d12
   logical, OPTIONAL, intent(in) :: addr2
   ! Output
   type(type_grid), pointer :: grid
   ! Local variables
   real(kind=dp), dimension(nshell(1)) :: radii1, radii_w1
   real(kind=dp), dimension(nshell(2)) :: radii2, radii_w2
   integer :: i, j, c, lower, upper, offset, ngrid
   real(kind=dp) :: R, r1, r2, mu, s1, s2, alpha

   alpha = pi/(2._dp*REAL(nshell(2)+1, dp))
   R = sqrt( sum( d12**2 ) )
   if (R .eq. 0.0_dp) then
      call build_onecenter_grid(ileb=ileb(1), nshell=nshell(1), addr2=addr2,&
                                grid=grid, quadr=1)
      return
   endif

   ngrid = lebedev_grid(ileb(1))%n * nshell(1)&
            +lebedev_grid(ileb(2))%n * nshell(2)
   call allocate_grid(grid=grid, n=ngrid)

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
         grid%r(lower+j-1, :) = radii1(j) * lebedev_grid(ileb(1))%r(:, i)
      enddo
      grid%w(lower:upper) = radii_w1 * 4.0_dp*pi * lebedev_grid(ileb(1))%w(i)
   enddo

   offset = lebedev_grid(ileb(1))%n * nshell(1)
   do i=1, lebedev_grid(ileb(2))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(2)
      upper = offset + i*nshell(2)

      ! The second batch of grid points is displaced by `d12`
      do j=1,nshell(2)
         grid%r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i) + d12
      enddo

      grid%w(lower:upper) = radii_w2 * 4.0_dp*pi *  lebedev_grid(ileb(2))%w(i)
   enddo

   ! nuclear partition
   !! A = 0, B = d12
   !! r1 is the distance A -> grid_r = grid_r-0 = grid_r
   !! r2 is the distance B -> grid_r = grid_r-d12
   do i=1,offset
      r1 = sqrt(sum( grid%r(i, :)**2 ))
      r2 = sqrt(sum( (grid%r(i, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu)
      s2 = s3(-mu)

      if (abs(s1+s2-1._dp) .gt. epsilon(1._dp)) stop 'Two-center nuclear partition'

      grid%w(i) = grid%w(i) * s1!/(s1+s2)
   enddo

   do i=1+offset,size(grid%w)
      r1 = sqrt(sum( grid%r(i, :)**2 ))
      r2 = sqrt(sum( (grid%r(i, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu)
      s2 = s3(-mu)

      if (abs(s1+s2-1._dp) .gt. epsilon(1._dp)) stop 'Two-center nuclear partition'
      
      grid%w(i) = grid%w(i) * s2!/(s1+s2)
   enddo

   grid%dw = 0._dp
   c = 0
   do i=1, lebedev_grid(ileb(1))%n
   do j=1, nshell(1)
      c = c+1
      r1 = sqrt(sum( grid%r(c, :)**2 ))
      r2 = sqrt(sum( (grid%r(c, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu); s2 = s3(-mu)
      ! dw/dX = (1)w_rad * ds_3/dmu*dmu/dX + (2) dw_rad/dr * dr/dX * wpart
      ! (1):
      grid%dw(c, 1) = radii_w1(j) * ds3dmu(mu)&
         * (d12(1)/R**3*(r2-r1) + (grid%r(c, 1) - d12(1)))/(R*r2)

      grid%dw(c, 2) = radii_w1(j) * ds3dmu(mu)&
         * (d12(2)/R**3*(r2-r1) + (grid%r(c, 2) - d12(2)))/(R*r2)

      grid%dw(c, 3) = radii_w1(j) * ds3dmu(mu)&
         * (d12(3)/R**3*(r2-r1) + (grid%r(c, 3) - d12(3)))/(R*r2)

      grid%dw(c, :) = grid%dw(c, :) * 4._dp*pi*lebedev_grid(ileb(1))%w(i)
   enddo
   enddo

   do i=1, lebedev_grid(ileb(2))%n
   do j=1, nshell(2)
      c = c+1
      r1 = sqrt(sum( grid%r(c, :)**2 ))
      r2 = sqrt(sum( (grid%r(c, :) - d12)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu); s2 = s3(-mu)
      ! dw/dX = (1)w_rad * ds_3/dmu*dmu/dX + (2) dw_rad/dr * dr/dX * wpart
      ! (1):
      grid%dw(c, 1) = radii_w1(j) * ds3dmu(mu)&
         * (d12(1)/R**3*(r2-r1) + (grid%r(c, 1) - d12(1)))/(R*r2)& ! (2)
         - (grid%r(c, 1) - d12(1)) * alpha * (7._dp*r2**1.5_dp + 5._dp*sqrt(r2))*s2

      grid%dw(c, 2) = radii_w1(j) * ds3dmu(mu)&
         * (d12(2)/R**3*(r2-r1) + (grid%r(c, 2) - d12(2)))/(R*r2)&
         - (grid%r(c, 2) - d12(2)) * alpha * (7._dp*r2**1.5_dp + 5._dp*sqrt(r2))*s2

      grid%dw(c, 3) = radii_w1(j) * ds3dmu(mu)&
         * (d12(3)/R**3*(r2-r1) + (grid%r(c, 3) - d12(3)))/(R*r2)&
         - (grid%r(c, 3) - d12(3)) * alpha * (7._dp*r2**1.5_dp + 5._dp*sqrt(r2))*s2

      grid%dw(c, :) = grid%dw(c, :) * 4._dp*pi*lebedev_grid(ileb(2))%w(i)
   enddo
   enddo
end subroutine build_twocenter_grid

! **********************************************
!> \brief Generates the three-center grid
! **********************************************
subroutine build_threecenter_grid(ileb, nshell, d12, d13, addr2, grid)
   implicit none
   ! Input
   integer, dimension(3), intent(in) :: ileb, nshell
   real(kind=dp), dimension(3), intent(in) :: d12, d13
   logical, OPTIONAL, intent(in) :: addr2
   ! Output
   type(type_grid), pointer :: grid
   ! Local variables
   real(kind=dp), dimension(nshell(1)) :: radii1, radii_w1
   real(kind=dp), dimension(nshell(2)) :: radii2, radii_w2
   real(kind=dp), dimension(nshell(3)) :: radii3, radii_w3
   integer :: i, j, lower, upper, offset, off1, off2, ngrid
   real(kind=dp) :: R12, R13, R23, r1, r2, r3,&
                    mu12, mu13, mu23, s12, s13, s23, s21, s31, s32

   real(kind=dp) :: tP1, tP2, tP3, sP, p1, p2, p3
   logical :: myaddr2 = .true.

   ngrid = lebedev_grid(ileb(1))%n * nshell(1)&
            +lebedev_grid(ileb(2))%n * nshell(2)&
            +lebedev_grid(ileb(3))%n * nshell(3)
   call allocate_grid(grid=grid, n=ngrid)

   if(present(addr2)) myaddr2 = addr2

   call radial_grid(r=radii1, &
                    wr=radii_w1, &
                    n=nshell(1), addr2=myaddr2, quadr=1)

   call radial_grid(r=radii2, &
                    wr=radii_w2, &
                    n=nshell(2), addr2=myaddr2, quadr=1)

   call radial_grid(r=radii3, &
                    wr=radii_w3, &
                    n=nshell(3), addr2=myaddr2, quadr=1)

   R12 = sqrt(sum(d12**2))
   R13 = sqrt(sum(d13**2))
   R23 = sqrt(sum((d13-d12)**2))

   ! Center 1
   do i=1, lebedev_grid(ileb(1))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = 1+(i-1)*nshell(1)
      upper = i*nshell(1)

      do j=1,nshell(1)
          grid%r(lower+j-1,:) = radii1(j) * lebedev_grid(ileb(1))%r(:, i)
      enddo

      grid%w(lower:upper) = 4.0_dp*pi * radii_w1 * lebedev_grid(ileb(1))%w(i)
   enddo

   ! Center 2
   offset = lebedev_grid(ileb(1))%n * nshell(1)
   do i=1, lebedev_grid(ileb(2))%n
       ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(2)
      upper = offset + i*nshell(2)

      do j=1,nshell(2)
          grid%r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i) + d12
      enddo

      grid%w(lower:upper) = 4.0_dp*pi * radii_w2 * lebedev_grid(ileb(2))%w(i)
   enddo

   ! Center 3
   offset = offset + lebedev_grid(ileb(2))%n * nshell(2)
   do i=1, lebedev_grid(ileb(3))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(3)
      upper = offset + i*nshell(3)

      do j=1,nshell(3)
          grid%r(lower+j-1,:) = radii3(j) * lebedev_grid(ileb(3))%r(:, i) + d13
      enddo

      grid%w(lower:upper) = 4.0_dp*pi * radii_w3 * lebedev_grid(ileb(3))%w(i)
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
      r1 = sqrt(sum( grid%r(i, :)**2 ))
      r2 = sqrt(sum( (grid%r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid%r(i, :) - d13)**2 ))

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

      grid%w(i) = grid%w(i) * p1

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
      r1 = sqrt(sum( grid%r(i, :)**2 ))
      r2 = sqrt(sum( (grid%r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid%r(i, :) - d13)**2 ))

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
      
      grid%w(i) = grid%w(i) * p2
      if (abs(mu12).gt.1.0_dp .or. abs(mu13).gt.1.0_dp .or. abs(mu23).gt.1.0_dp) then
         stop 'Nuclear partition - Three-center'
      endif
   enddo

   do i=off2+1,size(grid%w)
      r1 = sqrt(sum( grid%r(i, :)**2 ))
      r2 = sqrt(sum( (grid%r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid%r(i, :) - d13)**2 ))

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
      
      grid%w(i) = grid%w(i) * p3
      if (abs(mu12).gt.1.0_dp .or. abs(mu13).gt.1.0_dp .or. abs(mu23).gt.1.0_dp) then
         stop 'Nuclear partition - Three-center'
      endif
   enddo
end subroutine build_threecenter_grid

! Not used anymore (?)
   function h(mu)
      implicit none
      real(kind=dp) :: mu, h
      h = 1.5_dp*mu-0.5_dp*mu**3
   end function h

   function z(mu)
      implicit none
      real(kind=dp) :: mu, mua, mua2, mua4, mua6, z
      mua = (mu/0.64_dp)
      mua2 = mua*mua
      mua4 = mua2*mua2
      mua6 = mua4*mua2
      z = mua*(35.0_dp*(1-mua2)+21.0_dp*mua4-5.0_dp*mua6)/16.0_dp
   end function z

! **********************************************
!> \brief The nuclear partition function (see Becke)
! **********************************************
function s3(mu)
   implicit none
   real(kind=dp) :: mu, a, a3, s3, mu3
   ! s3 = 0.5_dp*(1._dp - h(h(h(mu))) )
   mu3 = mu**3
   a = mu3 - 3*mu
   a3 = a**3

   s3 = (a3 - 12._dp*mu3 + 36._dp*mu)**3/16384._dp -&
         .046875_dp*a3 + 0.5625_dp*mu3 - 1.6875_dp*mu + 0.5_dp
end function s3

! **********************************************
!> \brief The derivative of the nuclear partition function wrt. mu
! **********************************************
function ds3dmu(mu)
   implicit none
   real(kind=dp) :: mu, mu2, mu3, a, a2, a3, ds3dmu
   mu2 = mu*mu; mu3 = mu2*mu
   a = mu3 - 3._dp*mu
   a2 = a*a; a3 = a2*a

   ds3dmu = 27._dp/16384._dp * ( a3 - 12._dp*mu3 + 36._dp*mu )**2&
             * (a2*(mu2-1._dp) - 4._dp*mu2 + 4._dp) &
             - 27._dp/64._dp * a2*(mu2-1._dp) + 27._dp/16._dp * (mu2 - 1)
end function ds3dmu


! **********************************************
!> \brief The size of the grid should somehow depend on the size of the atom.
! **********************************************
   subroutine grid_parameters(atom, nleb, nshell)
      implicit none
      integer :: atom, nleb, nshell

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

! **********************************************
!> \brief Compute the `n`th Hermite polynomial y=H_n(x).
! **********************************************
recursive subroutine hermite(n, x, y)
   implicit none
   integer, intent(in) :: n
   real(kind=dp), intent(in) :: x
   real(kind=dp) :: y, a, b

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
! This takes some time
! Numerical recipes
subroutine gauher(r, wr, n)
   implicit none
   ! Input
   integer :: n
   ! Output
   real(kind=dp), dimension(n) :: r, wr
   ! Local variables
   integer :: i, j, k, m, maxit
   real(kind=dp) :: pim4, p1, p2, p3, pp, z, z1

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

end module ni_grid
