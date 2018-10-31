module nao_unit
USE eddi, ONLY: integration_onecenter, integration_twocenter, integration_threecenter, &
                kinetic_energy, coulomb_integral, spline, interpolation,&
                forward_derivative_weights, bisection, derivative_point
USE lebedev, ONLY: dp
USE grid, ONLY: radial_grid
implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi
contains

subroutine test_derivative_point()
   implicit none
   REAL(KIND=dp), DIMENSION(500) :: r, wr, y, y1_exact, y1_approx
   REAL(KIND=dp), DIMENSION(500-3) :: errors
   REAL(KIND=dp) :: alpha, r0
   INTEGER :: i
   INTEGER :: low, upper, atindex

   ! Set up the function to test
   call radial_grid(r=r, wr=wr, n=size(r), addr2=.FALSE., quadr=1)
   call RANDOM_NUMBER(alpha); alpha = alpha * 5
   y = exp(-alpha * r**2)
   y1_exact = -2.0_dp * alpha * r * y

   errors = 0._dp
   do i=1,size(errors)
      call bisection(r=r, r0=r(i), low=low, upper=upper)
      call derivative_point(r=r, y=y, r0=r(i), y1=y1_approx(i))
      if (y1_exact(i) .ne. 0._dp) then
         errors(i) = abs((y1_exact(i)-y1_approx(i))/y1_exact(i))
      endif
      ! print *, r(i), y1_exact(i), y1_approx(i), errors(i)
   enddo

   ! Oh my.
   if ( count(errors>0.2)/REAL(size(errors), dp) < 0.2_dp&
         .and. sum(errors)/size(errors) < 0.1_dp ) then
      ! The test passes, if less than 20% of all derivatives are ~bad~ awful
      print *, '👌 test_derivative_point – radial grid - passed'
   else
      print *, '💣 test_derivative_point – radial grid - failed'
      print *, count(errors>0.2)/REAL(size(errors), dp)
      print *, sum(errors)/size(errors)
   endif

   ! Should be much better on an equally spaced grid
   r = (/ ( 0.0001_dp*REAL(i, dp), i = 1,size(r) ) /)
   y = exp(-alpha * r**2)
   y1_exact = -2.0_dp * alpha * r * y
   errors = 0._dp
   do i=1,size(errors)
      call bisection(r=r, r0=r(i), low=low, upper=upper)
      call derivative_point(r=r, y=y, r0=r(i), y1=y1_approx(i))
      if (y1_exact(i) .ne. 0._dp) then
         errors(i) = abs((y1_exact(i)-y1_approx(i))/y1_exact(i))
      endif
   enddo
   if (sum(errors)/size(errors) .lt. 0.01) then
      print *, '👌 test_derivative_point – equally spaced grid - passed '
   else
      print *, '💣 test_derivative_point – equally spaced grid - failed '
      print *, count(errors>0.2)/REAL(size(errors), dp)
      print *, sum(errors)/size(errors)
   endif
end subroutine test_derivative_point

subroutine test_derivative(ntests, loud)
   implicit none
   ! Input
   INTEGER, intent(in) :: ntests
   LOGICAL, intent(in) :: loud
   ! Local variables
   REAL(KIND=dp), DIMENSION(100) :: r, wr, y, ys, y1_s, y2_s, y1_ex, y2_ex
   REAL(KIND=dp) :: alpha
   INTEGER :: i
   REAL(KIND=dp), DIMENSION(2,3,5) :: coeff

   ! Set up the gaussian function and its derivatives analytically
   call RANDOM_NUMBER(alpha)
   alpha = alpha * 5

   call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)
   ! r = (/ ( 0.1_dp*i**2, i = 1,size(r) ) /)

   y = exp(-alpha * r**2)
   y1_ex = -2.0_dp * alpha * r * y
   y2_ex = (4.0_dp * alpha**2 * r**2 - 2.0_dp * alpha)*y

   ! Get the derivatives via spline interpolation
   call spline(r, y, size(r), y2_s)
   do i=1,size(r)
      call interpolation(gr=r, gy=y, spline=y2_s, r=r(i), y=ys(i), yprime=y1_s(i))
      print *, 'r:' , r(i)
      print *, 'y_exact, y_spline:', y(i), ys(i), y(i)-ys(i)
      print *, 'y1_exact, y1_spline:', y1_ex(i), y1_s(i), abs(1.-y1_s(i)/y1_ex(i))/r(i)
      print *,
   enddo
   print *, alpha

   open(unit=100, file='poly')
   do i=1,size(r)
      write(100, *) r(i), y(i), y1_ex(i), y1_s(i), y2_ex(i), y2_s(i)
   enddo
   close(100)
end subroutine test_derivative

subroutine test_spline()
   implicit none
   ! Local variables
   REAL(KIND=dp), DIMENSION(30) :: r, wr, y, y1_exact, y2_exact, y2_spline
   REAL(KIND=dp) :: alpha
   INTEGER :: i

   call RANDOM_NUMBER(alpha)
   alpha = alpha * 5
   call radial_grid(r=r, wr=wr, n=size(r), addr2=.FALSE., quadr=1)
   y = exp(-alpha * r**2)
   y1_exact = -2.0_dp * alpha * r * y
   y2_exact = (4.0_dp * alpha**2 * r**2 - 2.0_dp * alpha)*y

   call spline(r=r, y=y, n=size(r), yspline=y2_spline)
   print *, 'alpha=', alpha
   do i=1,size(r)
      print *, r(i), y2_exact(i), y2_spline(i), abs(1-y2_exact(i)/y2_spline(i))
   enddo
end subroutine test_spline

subroutine test_forward_deriv_coeff()
   implicit none
   REAL(KIND=dp), DIMENSION(2,3,5) :: c
   REAL(KIND=dp), DIMENSION(7) :: tr
   INTEGER :: i
   LOGICAL :: failed = .FALSE.

   tr = (/(REAL(i, dp), i=0,6)/)
   ! call radial_grid(r=tr, wr=twr, n=10, addr2=.TRUE., quadr=1)

   call forward_derivative_weights(order=2, x0=0.0_dp, r=tr, coeff=c)
   failed = failed .or. all(c(1,1,1:4) .ne. (/-1._dp, 1._dp, 0._dp, 0._dp/))
   failed = failed .or. all(c(1,2,1:4) .ne. (/-1.5_dp, 2._dp, -0.5_dp, 0._dp/))
   failed = failed .or. all(c(2,1,1:4) .ne. (/1._dp, -2._dp, 1._dp, 0._dp/))
   failed = failed .or. all(c(2,2,1:4) .ne. (/2._dp, -5._dp, 4._dp, -1._dp/))
   if (failed .eqv. .TRUE.) then
      print *, 'test forward derivative coefficients - failed 💣'
      print *, c(1,1,:)
      print *, c(1,2,:)
      print *, c(2,1,:)
      print *, c(2,2,:)
   else
      print *, 'test forward derivative coefficients - passed 👌'
   endif

   print *, REPEAT('-', 20) // ' BACKWARDS ' // REPEAT('-', 20)
   call forward_derivative_weights(order=2, x0=6._dp, r=tr(7:1:-1), coeff=c)
   print *, c(1,1,:)
end subroutine test_forward_deriv_coeff

subroutine test_onecenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp) :: rand
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, y, spline1, wr
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing One-Center ' // REPEAT('-', 30)
   ngrid = 100

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y(ngrid))
   allocate(spline1(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)
   ! Perform the tests
   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand)
      rand = rand * 10.0_dp + 0.1_dp

      ! Prepare grids
      y = exp(-rand * r**2 )
      call spline(r, y, size(r), spline1)
      call integration_onecenter(nang=50, nshell=ngrid, r=r, y=y,&
                                 spline=spline1, quadr=1, integral=integral)

      ri = 4.0_dp*pi**1.5_dp/(4.0_dp*rand**(1.5_dp))

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.00001_dp)) then
         print *, 'Exponents: ', rand
         print *, 'Is: ', integral
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral-ri)
         print *, 'Relative Error: ', errors(j)
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err
   print *, REPEAT('-', 28) // ' End Testing One-Center ' // REPEAT('-', 28)
   print *, ''

   deallocate(r)
   deallocate(wr)
   deallocate(y)
   deallocate(spline1)
end subroutine test_onecenter

subroutine test_twocenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp), DIMENSION(2) :: rand2, nshell_rand
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, y1, y2, spline1, spline2, wr
   INTEGER, DIMENSION(2) :: nshell
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Two-Center ' // REPEAT('-', 30)
   ngrid = 5000

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(spline1(ngrid))
   allocate(spline2(ngrid))

   ! Discard wr
   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos)
      rand_pos = rand_pos * sqrt(5.0_dp)
      ! nshell
      CALL RANDOM_NUMBER(nshell_rand)
      nshell = (/ 75, 75 /) + INT(50 * nshell_rand)

      y1 = exp(-rand2(1) * r**2 )
      y2 = exp(-rand2(2) * r**2 )
      call spline(r, y1, size(r), spline1)
      call spline(r, y2, size(r), spline2)

      call integration_twocenter(nang=(/590, 590/), nshell=nshell, d12=rand_pos, &
                                 r1=r, y1=y1, r2=r, y2=y2,&
                                 spline1=spline1, spline2=spline2, integral=integral)

      ri = (pi/sum(rand2))**(1.5_dp)
      ri = ri * exp(-rand2(1)*rand2(2)*sum(rand_pos**2)/sum(rand2))

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.00001_dp)) then
         print *, 'Exponents: ', rand2
         print *, 'Distance: ', sqrt(sum(rand_pos**2))
         print *, 'Is: ', integral
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral-ri)
         print *, 'Relative Error: ', errors(j)
         print *, ''
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err
   print *, REPEAT('-', 28) // ' End Testing Two-Center ' // REPEAT('-', 28)
   print *, ''

   deallocate(r)
   deallocate(wr)
   deallocate(y1)
   deallocate(y2)
   deallocate(spline1)
   deallocate(spline2)
end subroutine test_twocenter

subroutine test_threecenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: abc, integral, ri, err
   REAL(KIND=dp), DIMENSION(3) :: rand3, rand_pos1, rand_pos2, exparg, nshell_rand
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr, y1, y2, y3, s1, s2, s3
   INTEGER, DIMENSION(3) :: nshell, nang
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Three-Center ' // REPEAT('-', 30)
   ngrid = 5000

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(y3(ngrid))
   allocate(s1(ngrid))
   allocate(s2(ngrid))
   allocate(s3(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand3)
      rand3 = rand3 * 5.0_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos1)
      rand_pos1 = rand_pos1 * sqrt(5.0_dp)
      CALL RANDOM_NUMBER(rand_pos2)
      rand_pos2 = rand_pos2 * sqrt(5.0_dp)
      ! nshell
      CALL RANDOM_NUMBER(nshell_rand)
      nshell = (/ 50, 50, 50 /) + INT(50 * nshell_rand)
      nang = (/590, 590, 590/)

      y1 = exp(-rand3(1) * r**2 )
      y2 = exp(-rand3(2) * r**2 )
      y3 = exp(-rand3(3) * r**2 )
      call spline(r, y1, size(r), s1)
      call spline(r, y2, size(r), s2)
      call spline(r, y3, size(r), s3)

      call integration_threecenter(nang=nang, nshell=nshell,&
         d12=rand_pos1, d13=rand_pos2, r1=r, y1=y1, r2=r, y2=y2, r3=r, y3=y3,&
         spline1=s1, spline2=s2, spline3=s3, integral=integral)

      abc = sum(rand3)

      exparg(1) = (rand3(2)*rand_pos1(1) + rand3(3)*rand_pos2(1))**2/abc
      exparg(1) = exparg(1) - (rand3(2)*rand_pos1(1)**2 + rand3(3)*rand_pos2(1)**2)

      exparg(2) = (rand3(2)*rand_pos1(2) + rand3(3)*rand_pos2(2))**2/abc
      exparg(2) = exparg(2) - (rand3(2)*rand_pos1(2)**2 + rand3(3)*rand_pos2(2)**2)

      exparg(3) = (rand3(2)*rand_pos1(3) + rand3(3)*rand_pos2(3))**2/abc
      exparg(3) = exparg(3) - (rand3(2)*rand_pos1(3)**2 + rand3(3)*rand_pos2(3)**2)

      ri = (pi/abc)**(1.5_dp) * exp(sum(exparg))

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.00001_dp)) then
         print *, 'Exponents: ', rand3
         print *, 'Distances: ', sqrt(sum(rand_pos1**2)), sqrt(sum(rand_pos2**2))
         print *, 'Is: ', integral
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral-ri)
         print *, 'Relative Error: ', errors(j)
         print *, ''
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err

   print *, REPEAT('-', 28) // ' End Testing Three-Center ' // REPEAT('-', 28)
   print *, ''
   ! Prepare the grid
   deallocate(r)
   deallocate(wr)
   deallocate(y1)
   deallocate(y2)
   deallocate(y3)
   deallocate(s1)
   deallocate(s2)
   deallocate(s3)
end subroutine test_threecenter

subroutine test_kinetic(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp), DIMENSION(2) :: rand2
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr, y1, y2, spline1, spline2, d2f2
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Kinetic energy ' // REPEAT('-', 30)
   ngrid = 25000

   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(spline1(ngrid))
   allocate(spline2(ngrid))
   allocate(d2f2(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp

      ! Prepare grids
      y1 = exp( -rand2(1) * r**2 )
      y2 = exp( -rand2(2) * r**2 )
      call spline(r, y1, size(r), spline1)
      call spline(r, y2, size(r), spline2)

      ! The result we get from subroutine kinetic_energy
      call kinetic_energy(nang=1, nshell=100,&
                          r1=r, y1=y1, r2=r, y2=y2,&
                          spline1=spline1, spline2=spline2, integral=integral)

      ! The result we want to have
      ri = 3.0_dp*rand2(2)*(pi/(sum(rand2)))**1.5_dp&
           - 3.0_dp*rand2(2)**2*sqrt(pi**3/(sum(rand2)**5))

      ! ! The result we want to have by one-center integration, analytically
      ! y1 = exp(-sum(rand2) * r**2)
      ! y1 = y1 * (3.0_dp*rand2(2) - 2.0_dp*rand2(2)**2 * r**2)
      ! call integration_onecenter(nang=590, nshell=100, r=r, y=y1,&
      !                            spline=spline1, integral=ri)
      ! ! 2.9644830114845719

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.000001_dp)) then
         print *, 'Exponents: ', rand2
         print *, 'Is: ', integral
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral-ri)
         print *, 'Relative Error: ', errors(j)
         print *, ''
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err

   print *, REPEAT('-', 28) // ' End Testing Kinetic energy ' // REPEAT('-', 28)

   deallocate(r)
   deallocate(wr)
   deallocate(y1)
   deallocate(y2)
   deallocate(spline1)
   deallocate(spline2)
   deallocate(d2f2)
end subroutine test_kinetic

subroutine test_coulomb(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err, alpha, norm
   REAL(KIND=dp), DIMENSION(2) :: rand2, nshell_rand
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, y1, y2, spline1, spline2, wr
   INTEGER, DIMENSION(2) :: nshell
   INTEGER :: j, ntests, ngrid, coul_n

   print *, REPEAT('-', 30) // ' Testing Coulomb ' // REPEAT('-', 30)
   ngrid = 5000

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(spline1(ngrid))
   allocate(spline2(ngrid))

   ! Discard wr
   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos)
      rand_pos = rand_pos * sqrt(10.0_dp) - 5_dp
      ! nshell
      CALL RANDOM_NUMBER(nshell_rand)
      nshell = (/ 75, 75 /) + INT(50 * nshell_rand) + 1000
      ! coul_n 
      CALL RANDOM_NUMBER(nshell_rand)
      coul_n = 1000 + INT(50*nshell_rand(1))

      y1 = exp(-rand2(1) * r**2 )
      y2 = exp(-rand2(2) * r**2 )
      call spline(r, y1, size(r), spline1)
      call spline(r, y2, size(r), spline2)

      call coulomb_integral(nang=(/590, 590/), nshell=nshell, coul_n=coul_n, d12=rand_pos,&
                            r1=r, y1=y1, r2=r, y2=y2, s1=spline1, s2=spline2,&
                            integral=integral)

      alpha = product(rand2)/sum(rand2)
      norm = sqrt(sum(rand_pos**2))
      ri = pi**3 / sqrt(product(rand2))**3
      if (norm .eq. 0.0_dp) then
      ri = ri * 2.0_dp*sqrt(alpha/pi) ! lim x->0 : erf(a*x)/x = 2a/sqrt(pi)
      else
      ri = ri * erf( sqrt(alpha)*norm )/norm 
      endif

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.00001_dp)) then
         print *, 'Exponents: ', rand2
         print *, 'Distance: ', sqrt(sum(rand_pos**2))
         print *, 'Is: ', integral
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral-ri)
         print *, 'Relative Error: ', errors(j)
         print *, ''
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err
   print *, REPEAT('-', 28) // ' End Testing Coulomb ' // REPEAT('-', 28)
   print *, ''

   deallocate(r)
   deallocate(wr)
   deallocate(y1)
   deallocate(y2)
   deallocate(spline1)
   deallocate(spline2)
end subroutine test_coulomb

subroutine test_radial_weight_pos(ntests)
   implicit none
   INTEGER :: ntests
   REAL(KIND=dp), DIMENSION(ntests) :: errors

   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr
   INTEGER :: ngrid, i
   LOGICAL :: failed

   failed = .FALSE.
   do i=1,ntests
      ngrid = i**2
      allocate(r(ngrid))
      allocate(wr(ngrid))

      call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)
      failed = failed .or. any(wr .lt. 0)
      call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=2)
      failed = failed .or. any(wr .lt. 0)

      deallocate(r)
      deallocate(wr)
   enddo
   if (.not. failed) then
      print *, 'test_radial - all weights positive - passed 👌'
   else
      print *, 'test_radial - all weights positive - failed 💣'
   endif
end subroutine test_radial_weight_pos

subroutine test_radial_chebyherm(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral_c, integral_h, ri, err
   REAL(KIND=dp) :: rand
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r_c, r_h, y, spline1, wr
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 26) // ' Testing Hermite-Chebyshev ' // REPEAT('-', 26)
   ngrid = 40

   ! Prepare the grid
   allocate(r_c(ngrid))
   allocate(r_h(ngrid))
   allocate(wr(ngrid))
   allocate(y(ngrid))
   allocate(spline1(ngrid))

   call radial_grid(r=r_c, wr=wr, n=ngrid, addr2=.TRUE., quadr=1)
   call radial_grid(r=r_h, wr=wr, n=ngrid, addr2=.TRUE., quadr=2)

   ! Perform the tests
   do j=1,ntests
      print *, j
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand)
      rand = rand * 5.0_dp + 0.1_dp

      ! Prepare grids
      y = exp(-rand * r_c**2 )
      call spline(r_c, y, size(r_c), spline1)
      call integration_onecenter(nang=50, nshell=ngrid, r=r_c, y=y,&
                                 spline=spline1, quadr=1, integral=integral_c)

      y = exp(-rand * r_h**2 )
      call spline(r_h, y, size(r_h), spline1)
      call integration_onecenter(nang=50, nshell=ngrid, r=r_h, y=y,&
                                 spline=spline1, quadr=2, integral=integral_h)

      ri = 4.0_dp*pi**1.5_dp/(4.0_dp*rand**(1.5_dp))

      errors(j) = abs(1.0_dp-integral_h/integral_c)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 0.01_dp)) then
         print *, 'Exponents: ', rand
         print *, 'Hermite: ', integral_h
         print *, 'Chebyshev: ', integral_c
         print *, 'Should:', ri
         print *, ''
         print *, 'Absolute Difference: ', abs(integral_h-integral_c)
         print *, 'Relative difference: ', errors(j)
         print *, REPEAT('-', 80)
      endif
   enddo

   err = sum(errors)/REAL(ntests, dp)
   print *, 'Mean rel. difference: ', err
   print *, REPEAT('-', 24) // ' End Testing Hermite-Chebyshev ' // REPEAT('-', 24)
   print *, ''

   deallocate(r_c)
   deallocate(r_h)
   deallocate(wr)
   deallocate(y)
   deallocate(spline1)
end subroutine test_radial_chebyherm

subroutine test_radial_weight_asc(ntests)
   implicit none
   INTEGER :: ntests
   REAL(KIND=dp), DIMENSION(ntests) :: errors

   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr
   INTEGER :: ngrid, i, j
   LOGICAL :: failed

   failed = .FALSE.
   do i=1,ntests
      ngrid = i**2
      allocate(r(ngrid))
      allocate(wr(ngrid))

      call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=1)
      do j=1,ngrid-1
         failed = r(j) .ge. r(j+1)
      enddo
      call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE., quadr=2)
      do j=1,ngrid-1
         failed = r(j) .ge. r(j+1)
      enddo
      failed = failed .or. any(wr .lt. 0)

      deallocate(r)
      deallocate(wr)
   enddo
   if (.not. failed) then
      print *, 'test_radial - radii ascending order - passed 👌'
   else
      print *, 'test_radial - radii ascending order - failed 💣'
   endif
end subroutine test_radial_weight_asc

end module nao_unit