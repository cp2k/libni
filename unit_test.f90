module nao_unit
USE eddi, ONLY: integration_onecenter, integration_twocenter, integration_threecenter, &
                kinetic_energy, coulomb_integral, spline, interpolation,&
                forward_derivative_weights, bisection, derivative_point,&
                derivatives, kah_sum, fun_grid
USE lebedev, ONLY: dp
USE grid, ONLY: radial_grid
implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi
contains

! --  Tests concerning the integrals -- !
! When we go about this stuff, the relative error does not matter
! all that much, since we want to compute values ~1. If we demand
! accuracy up to ~1meV we end up at something like 1e-4 = 1e-3/27 (adv. math)
! unless the total error should be 1meV and we have to make sure, that
! 1000 * our error < 1meV, assuming that the errors just add linearly. Then
! we can deal with an error of 1e-5. Ö
subroutine test_onecenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp) :: rand
   REAL(KIND=dp), DIMENSION(15000) :: r, y, spline1
   INTEGER :: j, ntests

   print *, REPEAT('-', 30) // ' Testing One-Center ' // REPEAT('-', 30)

   do j=1,ntests
      CALL RANDOM_NUMBER(rand)
      rand = rand * 5.0_dp + 0.1_dp
      call fun_grid(r=r, max=75._dp/rand)

      ! Prepare function
      y = exp(-rand * r**2 )
      call spline(r, y, size(r), spline1)
      call integration_onecenter(nang=1, nshell=1000, r=r, y=y,&
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
end subroutine test_onecenter

subroutine test_twocenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp), DIMENSION(2) :: rand2, nshell_rand
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(15000) :: r1, r2, y1, y2, spline1, spline2
   INTEGER, DIMENSION(2) :: nshell
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Two-Center ' // REPEAT('-', 30)
   ngrid = 100

   errors = 0._dp
   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 4.5_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos)
      rand_pos = rand_pos * sqrt(5.0_dp)
      ! nshell
      CALL RANDOM_NUMBER(nshell_rand)
      nshell = (/ 75, 75 /) + INT(50 * nshell_rand)

      call fun_grid(r=r1, max=75._dp/rand2(1))
      call fun_grid(r=r2, max=75._dp/rand2(2))
      y1 = exp(-rand2(1) * r1**2 )
      y2 = exp(-rand2(2) * r2**2 )
      call spline(r1, y1, size(r1), spline1)
      call spline(r2, y2, size(r2), spline2)

      ! Prepare the grid

      call integration_twocenter(l=(/0, 0/), m=(/0, 0/), nshell=nshell, d12=rand_pos, &
                                 r1=r1, y1=y1, r2=r2, y2=y2,&
                                 spline1=spline1, spline2=spline2, integral=integral)

      ri = (pi/sum(rand2))**(1.5_dp)
      ri = ri * exp(-rand2(1)*rand2(2)*sum(rand_pos**2)/sum(rand2))
      ri = ri / (4._dp * pi)

      if (abs(integral-ri) .gt. 1.e-13_dp) then
         errors(j) = abs(1.0_dp-integral/ri)
      endif
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
end subroutine test_twocenter

subroutine test_threecenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: abc, integral, ri, err
   REAL(KIND=dp), DIMENSION(3) :: rand3, rand_pos1, rand_pos2, exparg, nshell_rand
   REAL(KIND=dp), DIMENSION(10000) :: r, y1, y2, y3, s1, s2, s3
   INTEGER, DIMENSION(3) :: nshell, nang
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Three-Center ' // REPEAT('-', 30)
   ngrid = 5000

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

      call fun_grid(r=r, max=75._dp/maxval(rand3))

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
end subroutine test_threecenter

subroutine test_kinetic(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp), DIMENSION(2) :: rand2
   REAL(KIND=dp), DIMENSION(10000) :: r, wr, y1, y2, spline1, spline2, d2f2
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Kinetic energy ' // REPEAT('-', 30)
   ngrid = 2500


   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp
      call fun_grid(r=r, max=75._dp/maxval(rand2))

      ! Prepare grids
      y1 = exp( -rand2(1) * r**2 )
      y2 = exp( -rand2(2) * r**2 )
      call spline(r, y1, size(r), spline1)
      call spline(r, y2, size(r), spline2)

      ! The result we get from subroutine kinetic_energy
      call kinetic_energy(l=(/0, 0/), m=(/0, 0/), nshell=(/125,125/),&
                          r1=r, y1=y1, r2=r, y2=y2, d12=(/0._dp, 0._dp, 0._dp/),&
                          spline1=spline1, spline2=spline2, integral=integral)

      ! The result we want to have
      ri = 3.0_dp*rand2(2)*(pi/(sum(rand2)))**1.5_dp&
           - 3.0_dp*rand2(2)**2*sqrt(pi**3/(sum(rand2)**5))
      ri = ri/(4._dp*pi)
      ! ! The result we want to have by one-center integration, analytically
      ! y1 = exp(-sum(rand2) * r**2)
      ! y1 = y1 * (3.0_dp*rand2(2) - 2.0_dp*rand2(2)**2 * r**2)
      ! call integration_onecenter(nang=590, nshell=100, r=r, y=y1,&
      !                            spline=spline1, integral=ri)
      ! ! 2.9644830114845719

      errors(j) = abs(1.0_dp-integral/ri)
      if ((loud .eqv. .TRUE.) .or. (errors(j) .gt. 1.e-5_dp)) then
         print *, j
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

   err = kah_sum(errors)/REAL(ntests, dp)
   print *, 'Mean error: ', err

   print *, REPEAT('-', 28) // ' End Testing Kinetic energy ' // REPEAT('-', 28)

end subroutine test_kinetic

subroutine test_coulomb(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err, alpha, norm
   REAL(KIND=dp), DIMENSION(2) :: rand2, nshell_rand
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(10000) :: r, y1, y2, spline1, spline2
   INTEGER, DIMENSION(2) :: nshell
   INTEGER :: j, ntests, ngrid, coul_n

   print *, REPEAT('-', 30) // ' Testing Coulomb ' // REPEAT('-', 30)
   ngrid = 5000

   ! Discard wr

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

      call fun_grid(r=r, max=75._dp/maxval(rand2))

      y1 = exp(-rand2(1) * r**2 )
      y2 = exp(-rand2(2) * r**2 )
      call spline(r, y1, size(r), spline1)
      call spline(r, y2, size(r), spline2)

      call coulomb_integral(nshell=nshell, coul_n=coul_n, d12=rand_pos,&
                            r1=r, y1=y1, r2=r, y2=y2, s1=spline1, s2=spline2,&
                            l=(/0,0/), m=(/0,0/), integral=integral)

      alpha = product(rand2)/sum(rand2)
      norm = sqrt(sum(rand_pos**2))
      ri = pi**3 / sqrt(product(rand2))**3
      if (norm .eq. 0.0_dp) then
      ri = ri * 2.0_dp*sqrt(alpha/pi) ! lim x->0 : erf(a*x)/x = 2a/sqrt(pi)
      else
      ri = ri * erf( sqrt(alpha)*norm )/norm / (4._dp*pi)
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
end subroutine test_coulomb

! --  Tests concerning radial grids and functions on those grids -- !

! All weights have to be positive.
! We check this for radial grids of n=1,..,ntests**2
subroutine test_radial_weight_pos(ntests)
   implicit none
   INTEGER :: ntests

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
      print *, '👌 test_radial - all weights positive - passed'
   else
      print *, '💣 test_radial - all weights positive - failed'
   endif
end subroutine test_radial_weight_pos

! The code assumes that all radial points are given in ASCending order.
subroutine test_radial_weight_asc(ntests)
   implicit none
   INTEGER :: ntests

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
      print *, '👌 test_radial - radii ascending order - passed'
   else
      print *, '💣 test_radial - radii ascending order - failed'
   endif
end subroutine test_radial_weight_asc

! When integrating a function, the result should be the same for
! Gauss-Tschebyshev and Gauss-Hermite integration.
subroutine test_radial_chebyherm(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral_c, integral_h, ri, err
   REAL(KIND=dp) :: rand
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r_c, r_h, y, spline1, wr
   INTEGER :: j, ntests, ngrid

   ngrid = 100

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
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand)
      rand = rand * 5.0_dp + 0.1_dp

      ! Prepare grids
      y = exp(-rand * r_c**2 )
      call spline(r_c, y, size(r_c), spline1)
      call integration_onecenter(nang=1, nshell=ngrid, r=r_c, y=y,&
                                 spline=spline1, quadr=1, integral=integral_c)

      y = exp(-rand * r_h**2 )
      call spline(r_h, y, size(r_h), spline1)
      call integration_onecenter(nang=1, nshell=ngrid, r=r_h, y=y,&
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

   err = sum(errors)/ntests
   if (maxval(errors) .lt. 1e-5_dp .and. sum(errors)/ntests .lt. 1e-5_dp) then
      print *, '👌 test_radial - cheby.herm - passed'
   else
      print *, '💣 test_radial - cheby.herm - failed'
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/ntests
   endif
   deallocate(r_c)
   deallocate(r_h)
   deallocate(wr)
   deallocate(y)
   deallocate(spline1)
end subroutine test_radial_chebyherm

! The spline should in principle be exactly the second derivative of the function
subroutine test_spline(ntests)
   implicit none
   INTEGER, intent(in) :: ntests
   ! Local variables
   REAL(KIND=dp), DIMENSION(15000) :: r, y, y2_exact, y2_spline, errors
   REAL(KIND=dp), DIMENSION(ntests) :: tot_errors
   REAL(KIND=dp) :: alpha, error_cutoff, abs_error
   INTEGER :: i, t

   do t=1,ntests
      call RANDOM_NUMBER(alpha)
      alpha = alpha * 5._dp + 0.1_dp
      ! Get the exact Gaussian and its 2nd derivative wrt. r
      call fun_grid(r=r, max=sqrt(75._dp/alpha))
      y = exp(-alpha * r**2)
      y2_exact = (4.0_dp * alpha**2 * r**2 - 2.0_dp * alpha)*y

      ! There needs to be some cutoff as (1e-100-1e-120)/1e-120 > 1e100
      error_cutoff = maxval(abs(y2_exact))*1.e-16_dp

      ! So we compute the spline and compare it point for point with the exact result
      call spline(r=r, y=y, n=size(r), yspline=y2_spline)
      errors = 0._dp
      do i=1,size(r)
         abs_error = abs(y2_spline(i)-y2_exact(i))
         if (y2_exact(i) .ne. 0._dp .and. abs_error .gt. error_cutoff) then
            errors(i) = abs_error/y2_exact(i)
         endif
      enddo
      tot_errors(t) = kah_sum(errors)/size(errors)
      if (tot_errors(t) .gt. 1.e-5_dp) then
         print *, 'alpha', alpha
         print *, 'error', tot_errors(t)
      endif
   enddo

   if (maxval(tot_errors) .lt. 1._dp .and. kah_sum(tot_errors)/ntests .lt. 0.01_dp) then
      print *, '👌 test_spline - passed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', kah_sum(tot_errors)/ntests
   else
      print *, '💣 test_spline - failed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', kah_sum(tot_errors)/ntests
   endif
end subroutine test_spline

! 
subroutine test_interpolation(ntests)
   implicit none
   REAL(KIND=dp), DIMENSION(20000) :: r_orig, y_orig, s_orig
   REAL(KIND=dp), DIMENSION(500) :: r_interp, wr2, errors
   REAL(KIND=dp), DIMENSION(ntests) :: tot_errors
   REAL(KIND=dp) :: alpha, y_interp, y_exact, abs_error
   INTEGER :: i, ntests, t

   ! We want to now the function on r_interp...
   call radial_grid(r=r_interp, wr=wr2, n=size(r_interp), addr2=.FALSE., quadr=1)
   do t=1,ntests
      call RANDOM_NUMBER(alpha); alpha = alpha * 5._dp + 0.2_dp
      alpha = 1.1_dp
      ! ... while knowing it on r_orig.
      call fun_grid(r=r_orig, max=(75._dp/alpha))
      ! The exact function and its spline
      y_orig = exp(-alpha * r_orig**2)
      call spline(r=r_orig, y=y_orig, n=size(r_orig), yspline=s_orig)
      ! accurate: call spline(r=r_orig, y=y_orig, n=size(r_orig), yspline=s_orig)
      ! fast    : s_orig = (4.0_dp * alpha**2 * r_orig**2 - 2.0_dp * alpha)*y_orig

      errors = 0._dp
      ! For every point on r_interp...
      do i=1,size(r_interp)
         ! ... we interpolate ...
         call interpolation(gr=r_orig, gy=y_orig, spline=s_orig, r=r_interp(i), y=y_interp)
         ! ... and compute the exact value.
         y_exact = exp(-alpha * r_interp(i)**2)

         abs_error = abs(y_interp-y_exact)
         if (y_exact .ne. 0._dp .and. abs_error .gt. epsilon(1._dp)) then
            errors(i) = abs_error/y_exact
         endif
      enddo
      tot_errors(t) = sum(errors)/size(errors)
      if (tot_errors(t) .gt. 1._dp) print *, 'alpha ', alpha
   enddo

   if (maxval(tot_errors) .lt. 0.3_dp .and. sum(tot_errors)/ntests .lt. 0.3_dp) then
      print *, '👌 test_interpolation - passed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   else
      print *, '💣 test_interpolation - failed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   endif
end subroutine test_interpolation

! --  Tests concerning derivatives -- !

! When computing the forward differences coefficients
! on an unequally spaced grid the coefficients will be non-obvious:
subroutine test_forward_deriv_coeff()
   implicit none
   REAL(KIND=dp), DIMENSION(3,3,5) :: c
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
      print *, '💣 test forward derivative coefficients - failed'
      print *, c(1,1,:)
      print *, c(1,2,:)
      print *, c(2,1,:)
      print *, c(2,2,:)
   else
      print *, '👌 test forward derivative coefficients - passed'
   endif
end subroutine test_forward_deriv_coeff

! We need the derivatives of y for further calculations
subroutine test_derivatives()
   implicit none
   REAL(KIND=dp), DIMENSION(15000) :: r, wr, y, yy1, y1, yy2, y2, yy3, y3
   REAL(KIND=dp), DIMENSION(3,size(r)) :: errors
   REAL(KIND=dp) :: error_cutoff, len
   INTEGER :: i
   
   ! The exact function and derivatives
   call fun_grid(r=r, max=sqrt(75._dp/0.5_dp))

   ! Derivatives are more accurate on equally spaced grids
   wr = 0._dp
   ! call radial_grid(r=r, wr=wr, n=size(r), addr2=.FALSE., quadr=1)
   y = exp(-0.5_dp * r**2)
   yy1 = -r * y
   yy2 = y * ( r**2 - 1._dp )
   yy3 = y * ( 3._dp * r - r**3 )

   error_cutoff = maxval(abs(y)) * epsilon(1._dp)

   errors = 0._dp

   call derivatives(r=r, y=y, y1=y1, y2=y2, y3=y3)
   do i=1,size(r)-5
      errors(1, i) = abs(y1(i)-yy1(i))
      if (yy1(i) .ne. 0._dp .and. errors(1, i) .gt. error_cutoff) then
         errors(1, i) = errors(1, i)/abs(yy1(i))
      endif
      errors(2, i) = abs(y2(i)-yy2(i))
      if (yy2(i) .ne. 0._dp .and. errors(1, i) .gt. error_cutoff) then
         errors(2, i) = errors(2, i)/abs(yy2(i))
      endif
      errors(3, i) = abs(y3(i)-yy3(i))
      if (yy3(i) .ne. 0._dp .and. errors(1, i) .gt. error_cutoff) then
         errors(3, i) = errors(3, i)/abs(yy3(i))
      endif
   enddo

   len = size(r)-5._dp

   if (maxval( errors(1, :) ) .lt. 1._dp .and.&
       sum(errors(1, :))/len .lt. 1.e-3_dp) then
      print *, '👌 test derivative - first - passed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   else
      print *, '💣 test derivative - first - failed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   endif

   if (maxval( errors(2, :) ) .lt. 1._dp .and.&
       sum(errors(2, :))/len .lt. 1.e-3_dp) then
      print *, '👌 test derivative - second - passed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   else
      print *, '💣 test derivative - second - failed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   endif

   if (maxval( errors(3, :) ) .lt. 1._dp .and.&
       sum(errors(3, :))/len .lt. 1.e-3_dp) then
      print *, '👌 test derivative - third - passed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   else
      print *, '💣 test derivative - third - failed'
      print *, 'Average error', sum(errors(1, :))/len
      print *, 'Max error', maxval( errors(1, :) )
   endif
end subroutine test_derivatives

! Sometimes we want to know the derivative at some point which is on the
! radial grid, the function is tabulated at
subroutine test_derivative_point_on()
   implicit none
   REAL(KIND=dp), DIMENSION(15000) :: r, wr, y, y1_exact, y1_approx
   REAL(KIND=dp), DIMENSION(size(r)-3) :: errors
   REAL(KIND=dp) :: alpha, abs_error, error_cutoff, summer
   INTEGER :: i

   call RANDOM_NUMBER(alpha); alpha = alpha * 5

   ! Set up the function to test
   ! call radial_grid(r=r, wr=wr, n=size(r), addr2=.FALSE., quadr=1)
   wr = 0._dp
   call fun_grid(r=r, max=sqrt(75._dp/alpha))
   y = exp(-alpha * r**2)
   y1_exact = -2.0_dp * alpha * r * y

   error_cutoff = maxval(abs(y1_exact)) * epsilon(1._dp)
   errors = 0._dp

   do i=1,size(errors)
      call derivative_point(r=r, y=y, r0=r(i), y1=y1_approx(i))
      abs_error = abs(y1_exact(i)-y1_approx(i))
      if (y1_exact(i) .ne. 0._dp .and. abs_error .gt. error_cutoff) then
         errors(i) = abs(abs_error/y1_exact(i))
      endif
      ! print *, r(i), y1_exact(i), y1_approx(i), errors(i)
   enddo

   summer = abs( 1._dp - sum(abs(y1_approx(1:size(errors))))/sum(abs(y1_exact(1:size(errors)))) )

   if ( summer .lt. 1.e-4 ) then
      print *, '👌 test_derivative_point_on – radial grid - passed'
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   else
      print *, '💣 test_derivative_point_on – radial grid - failed'
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   endif

   ! Should be much better on an equally spaced grid
   r = (/ ( 0.0001_dp*REAL(i, dp), i = 1,size(r) ) /)
   y = exp(-alpha * r**2)
   y1_exact = -2.0_dp * alpha * r * y

   error_cutoff = maxval(abs(y1_exact))*1e-16_dp
   errors = 0._dp
   do i=1,size(errors)
      call derivative_point(r=r, y=y, r0=r(i), y1=y1_approx(i))
      abs_error = abs(y1_exact(i)-y1_approx(i))
      if (y1_exact(i) .ne. 0._dp .and. abs_error .gt. error_cutoff) then
         errors(i) = abs(abs_error/y1_exact(i))
      endif
      ! print *, r(i), y1_exact(i), y1_approx(i), errors(i)
   enddo

   summer = abs( 1._dp - sum(abs(y1_approx(1:size(errors))))/sum(abs(y1_exact(1:size(errors)))) )

   if (summer .lt. 1.e-4) then
      print *, '👌 test_derivative_point_on – equally spaced grid - passed '
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   else
      print *, '💣 test_derivative_point_on – equally spaced grid - failed '
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   endif
end subroutine test_derivative_point_on

! ... and sometimes this point is not on the grid.
subroutine test_derivative_point_off()
   implicit none
   REAL(KIND=dp), DIMENSION(15000) :: r_exact, wr, y_exact
   REAL(KIND=dp), DIMENSION(100) :: r_appr, y1_appr
   REAL(KIND=dp), DIMENSION(size(r_appr)-3) :: errors
   REAL(KIND=dp) :: alpha, y1_exact, error_cutoff, abs_error, y1_exact_sum, summer
   INTEGER :: i

   ! Set up the function to test
   call radial_grid(r=r_exact, wr=wr, n=size(r_exact), addr2=.FALSE., quadr=1)
   call radial_grid(r=r_appr, wr=wr, n=size(r_appr), addr2=.FALSE., quadr=1)
   call RANDOM_NUMBER(alpha); alpha = alpha * 5
   y_exact = exp(-alpha * r_exact**2)

   error_cutoff = maxval(abs(-2.0_dp * alpha * r_appr * exp(-alpha * r_appr**2)))*1.e-16_dp

   errors = 0._dp
   y1_exact_sum = 0._dp
   do i=1,size(errors)
      call derivative_point(r=r_exact, y=y_exact, r0=r_appr(i), y1=y1_appr(i))
      y1_exact = -2.0_dp * alpha * r_appr(i) * exp(-alpha * r_appr(i)**2)
      y1_exact_sum = y1_exact_sum + abs(y1_exact)

      abs_error = abs(y1_exact-y1_appr(i))
      if (y1_exact .ne. 0._dp .and. abs_error .gt. error_cutoff) then
         errors(i) = abs_error/abs(y1_exact)
      endif
      ! print *, r_appr(i), y1_exact, y1_appr(i), errors(i)
   enddo

   summer = abs( sum(abs(y1_appr(1:size(errors))) )-y1_exact_sum )/y1_exact_sum
   if ( summer .lt. 1.e-4_dp ) then
      print *, '👌 test_derivative_point_off – radial grid - passed'
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   else
      print *, '💣 test_derivative_point_off – radial grid - failed'
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   endif

   ! Should be much better on an equally spaced grid
   r_exact = (/ ( 0.012_dp*REAL(i, dp), i = 0,size(r_exact)-1 ) /)
   y_exact = exp(-alpha * r_exact**2)

   y1_exact_sum = 0._dp
   errors = 0._dp
   do i=1,size(errors)
      call derivative_point(r=r_exact, y=y_exact, r0=r_appr(i), y1=y1_appr(i))
      y1_exact = -2.0_dp * alpha * r_appr(i) * exp(-alpha * r_appr(i)**2)
      y1_exact_sum = y1_exact_sum + abs(y1_exact)

      abs_error = abs(y1_exact-y1_appr(i))
      if (y1_exact .ne. 0._dp .and. abs_error .gt. error_cutoff) then
         errors(i) = abs(abs_error/y1_exact)
      endif
      ! print *, r_appr(i), y1_exact, y1_appr(i), errors(i)
   enddo

   summer = abs( sum(abs(y1_appr(1:size(errors))))-y1_exact_sum )/y1_exact_sum
   if (summer .lt. 1.e-4_dp) then
      print *, '👌 test_derivative_point_off – equally spaced grid - passed '
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   else
      print *, '💣 test_derivative_point_off – equally spaced grid - failed '
      print *, 'max. error: ', maxval(errors)
      print *, 'mean error: ', sum(errors)/size(errors)
      print *, 'error of sum: ', summer
   endif
end subroutine test_derivative_point_off

! We can retrieve the first derivative of a function when doing the interpolation,
! ...when the point is on the tabulated grid.
subroutine test_derivative_on(ntests)
   implicit none
   ! Input
   INTEGER, intent(in) :: ntests
   ! Local variables
   REAL(KIND=dp), DIMENSION(15000) :: r, y, ys, y1_s, y2_s, y1_ex, errors
   REAL(KIND=dp), DIMENSION(ntests) :: tot_errors
   REAL(KIND=dp) :: alpha, abs_error, error_cutoff
   INTEGER :: i, t

   tot_errors = 0.0_dp
   do t=1,ntests
      ! Set up the gaussian function and its derivatives analytically
      call RANDOM_NUMBER(alpha); alpha = alpha * 5.0_dp + 0.1_dp

      call fun_grid(r=r, max=sqrt(75._dp/alpha))
      y = exp(-alpha * r**2)
      y1_ex = -2.0_dp * alpha * r * y

      error_cutoff = maxval(abs(y1_ex))*1.e-16_dp
      errors = 0._dp

      ! Get the derivatives via spline interpolation
      call spline(r, y, size(r), y2_s)
      ! y2_s = (4.0_dp * alpha**2 * r**2 - 2.0_dp * alpha)*y
      do i=1,size(r)
         call interpolation(gr=r, gy=y, spline=y2_s, r=(r(i)), y=ys(i), yprime=y1_s(i))
         abs_error = abs(y1_s(i)-y1_ex(i))
         if (y1_ex(i) .ne. 0.0_dp .and. abs_error .gt. error_cutoff) then
            errors(i) = abs_error/abs(y1_ex(i))
         endif     
      enddo

      tot_errors(t) = kah_sum(errors)/size(errors)
   enddo

   if (maxval(tot_errors) .lt. 1.e-1_dp .and. sum(tot_errors)/ntests .lt. 1.e-3_dp) then
      print *, '👌 test_derivative - on grid - passed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   else
      print *, '💣 test_derivative - on grid - failed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   endif
end subroutine test_derivative_on

! ...and when it's not.
subroutine test_derivative_off(ntests)
   implicit none
   ! Input
   INTEGER, intent(in) :: ntests
   ! Local variables
   REAL(KIND=dp), DIMENSION(15000) :: r, y, ys, y2_s
   REAL(KIND=dp), DIMENSION(1000) :: r_interp, wr_interp, y1_s, y1_ex, errors
   REAL(KIND=dp), DIMENSION(ntests) :: tot_errors
   REAL(KIND=dp) :: alpha, abs_error, error_cutoff
   INTEGER :: i, t

   call radial_grid(r=r_interp, wr=wr_interp, n=size(r_interp), addr2=.TRUE., quadr=1)

   tot_errors = 0.0_dp
   do t=1,ntests
      ! Set up the gaussian function and its derivatives analytically
      call RANDOM_NUMBER(alpha); alpha = alpha * 5.0_dp + 0.1_dp
      call fun_grid(r=r, max=75._dp/alpha)

      y = exp(-alpha * r**2)

      y1_ex = -2.0_dp * alpha * r_interp * exp(-alpha * r_interp**2)

      error_cutoff = maxval(abs(y1_ex))*1.e-14_dp
      errors = 0._dp

      ! Get the derivatives via spline interpolation
      call spline(r, y, size(r), y2_s)
      do i=1,size(r_interp)
         call interpolation(gr=r, gy=y, spline=y2_s, r=r_interp(i), y=ys(i), yprime=y1_s(i))
         abs_error = abs(y1_s(i)-y1_ex(i))
         if (y1_ex(i) .ne. 0.0_dp .and. abs_error .gt. error_cutoff) then
            errors(i) = abs_error/abs(y1_ex(i))
         endif     
      enddo

      tot_errors(t) = sum(errors)/size(errors)
   enddo

   if (maxval(tot_errors) .lt. 1e-1_dp .and. sum(tot_errors)/ntests .lt. 1e-4_dp) then
      print *, '👌  test_derivative - off grid - passed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   else
      print *, '💣 test_derivative - off grid - failed'
      print *, 'max. error: ', maxval(tot_errors)
      print *, 'mean error: ', sum(tot_errors)/ntests
   endif
end subroutine test_derivative_off
end module nao_unit