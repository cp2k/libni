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
subroutine radial_grid(r, wr, n, addr2)
   implicit none
   INTEGER, intent(in) :: n
   LOGICAL, OPTIONAL, intent(in) :: addr2
   REAL(KIND=dp), DIMENSION(:) :: r, wr
   INTEGER :: i
   REAL(KIND=dp) :: alpha, t, x

   alpha = pi/REAL(n+1, dp)
   do i=1,n
      ! COORDINATE
      t = REAL(i, dp)*alpha
      x = COS(t)
      r(i) = (1.0_dp+x)/(1.0_dp-x)
      wr(i) = alpha*2.0_dp*SIN(t)/(1.0_dp-x)**2
   enddo
   if (present(addr2)) then
      ! dxdydz = dr r^2 dcos(theta) dphi
      wr = wr * r**2
   endif
end subroutine radial_grid

subroutine build_onecenter_grid(ileb, nshell, addr2, grid_r, grid_w)
   implicit none
   ! Input
   INTEGER, intent(in) :: ileb, nshell
   LOGICAL, OPTIONAL, intent(in) :: addr2
   ! Output
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   ! Local variables
   INTEGER :: i, j, lower, upper
   REAL(KIND=dp), DIMENSION(nshell) :: radii, radii_w

   call radial_grid(r=radii, &
                    wr=radii_w, &
                    n=nshell, addr2=.TRUE.)

   do i=1, lebedev_grid(ileb)%n
      lower = 1+(i-1)*nshell
      upper = i*nshell

      do j=1,nshell
         grid_r(lower+j-1,:) = radii(j) * lebedev_grid(ileb)%r(:, i)
      enddo
      grid_w(lower:upper) = 4.0_dp*pi * radii_w * lebedev_grid(ileb)%w(i)
   enddo
end subroutine build_onecenter_grid


subroutine build_twocenter_grid(ileb, nshell, displacement, addr2, grid_r, grid_w)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
   LOGICAL, OPTIONAL, intent(in) :: addr2
   ! Output
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   ! Local variables
   REAL(KIND=dp), DIMENSION(nshell(1)) :: radii1, radii_w1
   REAL(KIND=dp), DIMENSION(nshell(2)) :: radii2, radii_w2
   INTEGER :: i, j, lower, upper, offset
   REAL(KIND=dp) :: R, r1, r2, mu, s1, s2

   call radial_grid(r=radii1, &
                    wr=radii_w1, &
                    n=nshell(1), addr2=.TRUE.)

   call radial_grid(r=radii2, &
                    wr=radii_w2, &
                    n=nshell(2), addr2=.TRUE.)

   R = sqrt(sum(displacement**2))

   do i=1, lebedev_grid(ileb(1))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = 1+(i-1)*nshell(1)
      upper = i*nshell(1)

      ! the coordinates do not change, since we use `displacement` later on
      ! when evaluating the function(r)
      do j=1,nshell(1)
         grid_r(lower+j-1,:) = radii1(j) * lebedev_grid(ileb(1))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w1 * lebedev_grid(ileb(1))%w(i)
   enddo

   offset = lebedev_grid(ileb(1))%n * nshell(1)
   do i=1, lebedev_grid(ileb(2))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(2)
      upper = offset + i*nshell(2)

      ! the coordinates do not change, since we use `displacement` later on
      ! when evaluating the function(r)
      do j=1,nshell(2)
         grid_r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w2 * lebedev_grid(ileb(2))%w(i)
   enddo

   ! nuclear partition   
   do i=1,size(grid_w)
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - displacement)**2 ))
      mu = (r1-r2)/R
      s1 = s3(mu)
      s2 = s3(-mu)

      if (i .lt. (lebedev_grid(ileb(1))%n*nshell(1)+1)) then
         grid_w(i) = grid_w(i) * s1/(s1+s2)
      else
         grid_w(i) = grid_w(i) * s2/(s1+s2)
      endif
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
                    n=nshell(1), addr2=.TRUE.)

   call radial_grid(r=radii2, &
                    wr=radii_w2, &
                    n=nshell(2), addr2=.TRUE.)

   call radial_grid(r=radii3, &
                    wr=radii_w3, &
                    n=nshell(3), addr2=.TRUE.)

   R12 = sqrt(sum(d12**2))
   R13 = sqrt(sum(d13**2))
   R23 = sqrt(sum((d13-d12)**2))

   ! Let's us jump to the right index in grid_r, grid_w
   off1 = lebedev_grid(ileb(1))%n * nshell(1)+1
   off2 = off1 + lebedev_grid(ileb(2))%n * nshell(2)

   ! Center 1
   do i=1, lebedev_grid(ileb(1))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = 1+(i-1)*nshell(1)
      upper = i*nshell(1)

      ! the coordinates do not change, since we use `displacement` later on
      ! when evaluating the function(r)
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

      ! the coordinates do not change, since we use `displacement` later on
      ! when evaluating the function(r)
      do j=1,nshell(2)
          grid_r(lower+j-1,:) = radii2(j) * lebedev_grid(ileb(2))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w2 * lebedev_grid(ileb(2))%w(i)
   enddo

   ! Center 3
   offset = offset + lebedev_grid(ileb(2))%n * nshell(2)
   do i=1, lebedev_grid(ileb(3))%n
      ! lower:upper is the slice belonging to this angular point/ direction
      lower = offset + 1+(i-1)*nshell(3)
      upper = offset + i*nshell(3)

      ! the coordinates do not change, since we use `displacement` later on
      ! when evaluating the function(r)
      do j=1,nshell(3)
          grid_r(lower+j-1,:) = radii3(j) * lebedev_grid(ileb(3))%r(:, i)
      enddo

      grid_w(lower:upper) = 4.0_dp*pi * radii_w3 * lebedev_grid(ileb(3))%w(i)
   enddo

   ! nuclear partition
   do i=1,size(grid_w)
      r1 = sqrt(sum( grid_r(i, :)**2 ))
      r2 = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      r3 = sqrt(sum( (grid_r(i, :) - d13)**2 ))

      mu12 = (r1-r2)/R12
      mu13 = (r1-r3)/R13
      mu23 = (r2-r3)/R23

      s12 = 1.0_dp-z(mu12)
      s21 = 1.0_dp+z(mu12)
      s13 = 1.0_dp-z(mu13)
      s31 = 1.0_dp+z(mu13)
      s32 = 1.0_dp+z(mu23)
      s23 = 1.0_dp-z(mu23)

      tP1 = s12*s13
      tP2 = s21*s23
      tP3 = s31*s32

      sP = tP1+tP2+tP3

      p1 = tP1/sP
      p2 = tP2/sP
      p3 = tP3/sP

      if (i .lt. off1) then
         grid_w(i) = grid_w(i) * p1
      else if ((i .gt. (off1-1)) .and. (i .lt. off2) ) then
         grid_w(i) = grid_w(i) * p2
      else if (i .gt. (off2-1)) then
         grid_w(i) = grid_w(i) * p3
      else
         print *, i, '/', size(grid_w)
         stop 'Whoopsie'
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
      REAL(KIND=dp) :: mu, s3
      ! s3 = 0.5_dp*(1-h(h(h(mu))))
      s3 = 0.5_dp*(1-z(mu))
   end function s3

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

   subroutine print_grid(thegrid)
      implicit none
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      INTEGER :: i
      do i=1,size(thegrid)
         print *, thegrid(i)%r, thegrid(i)%weight
      enddo 
   end subroutine print_grid

end module grid