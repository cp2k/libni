module nao_unit
USE eddi, ONLY: integration_onecenter, integration_twocenter, integration_threecenter, kinetic_energy, spline
USE lebedev, ONLY: dp
USE grid, ONLY: radial_grid
implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi
contains

! Compute the integrals r^2 * exp(-a r^2)
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

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE.)
   r = r(ngrid:1:-1)

   ! Perform the tests
   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand)
      rand = rand * 10.0_dp + 0.1_dp

      ! Prepare grids
      y = exp(-rand * r**2 )
      call spline(r, y, size(r), 0.0_dp, 0.0_dp, spline1)

      call integration_onecenter(nang=50, nshell=ngrid, r=r, y=y,&
                                 spline=spline1, integral=integral)

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
end subroutine test_onecenter

subroutine test_twocenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: integral, ri, err
   REAL(KIND=dp), DIMENSION(2) :: rand2
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, y1, y2, spline1, spline2, wr
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Two-Center ' // REPEAT('-', 30)
   ngrid = 10000

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE.)
   r = r(ngrid:1:-1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos)
      rand_pos = rand_pos * sqrt(3.0_dp)

      y1 = exp(-rand2(1) * r**2 )
      y2 = exp(-rand2(2) * r**2 )
      call spline(r, y1, size(r), 0.0_dp, 0.0_dp, spline1)
      call spline(r, y2, size(r), 0.0_dp, 0.0_dp, spline2)

      call integration_twocenter(nang=(/590, 590/), nshell=(/75, 75/), d12=rand_pos, &
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
   print *, REPEAT('-', 28) // ' End Testing One-Center ' // REPEAT('-', 28)
   print *, ''

   deallocate(r)
   deallocate(wr)
   deallocate(y1)
   deallocate(y2)
end subroutine test_twocenter

subroutine test_threecenter(ntests, loud)
   implicit none
   LOGICAL :: loud
   REAL(KIND=dp), DIMENSION(ntests) :: errors
   REAL(KIND=dp) :: abc, exparg, integral, ri, err
   REAL(KIND=dp), DIMENSION(3) :: rand3, rand_pos1, rand_pos2
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr, y1, y2, y3, s1, s2, s3
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Three-Center ' // REPEAT('-', 30)
   ngrid = 5000

   ! Prepare the grid
   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(y3(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE.)
   r = r(ngrid:1:-1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand3)
      rand3 = rand3 * 5.0_dp + 0.5_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos1)
      rand_pos1 = rand_pos1 * sqrt(3.0_dp)
      CALL RANDOM_NUMBER(rand_pos2)
      rand_pos2 = rand_pos2 * sqrt(3.0_dp)

      y1 = exp(-rand3(1) * r**2 )
      y2 = exp(-rand3(2) * r**2 )
      y3 = exp(-rand3(3) * r**2 )
      call spline(r, y1, size(r), 0.0_dp, 0.0_dp, s1)
      call spline(r, y2, size(r), 0.0_dp, 0.0_dp, s2)
      call spline(r, y3, size(r), 0.0_dp, 0.0_dp, s3)

      call integration_threecenter(nang=(/590, 590, 590/), nshell=(/75, 75, 75/),&
         d12=rand_pos1, d13=rand_pos2, r1=r, y1=y1, r2=r, y2=y2, r3=r, y3=y3,&
         spline1=s1, spline2=s2, spline3=s3, integral=integral)

      abc = sum(rand3)

      ri = (pi/abc)**(1.5_dp)
      exparg = (rand3(2)*rand_pos1(1) + rand3(3)*rand_pos2(1))**2/abc
      exparg = exparg - (rand3(2)*rand_pos1(1)**2 + rand3(3)*rand_pos2(1)**2)
      ri = ri * exp(exparg)
      exparg = (rand3(2)*rand_pos1(2) + rand3(3)*rand_pos2(2))**2/abc
      exparg = exparg - (rand3(2)*rand_pos1(2)**2 + rand3(3)*rand_pos2(2)**2)
      ri = ri * exp(exparg)
      exparg = (rand3(2)*rand_pos1(3) + rand3(3)*rand_pos2(3))**2/abc
      exparg = exparg - (rand3(2)*rand_pos1(3)**2 + rand3(3)*rand_pos2(3)**2)
      ri = ri * exp(exparg)

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
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r, wr, y1, y2, spline1, spline2, d2f2
   INTEGER :: j, ntests, ngrid

   print *, REPEAT('-', 30) // ' Testing Kinetic energy ' // REPEAT('-', 30)
   ngrid = 10000

   allocate(r(ngrid))
   allocate(wr(ngrid))
   allocate(y1(ngrid))
   allocate(y2(ngrid))
   allocate(d2f2(ngrid))

   call radial_grid(r=r, wr=wr, n=ngrid, addr2=.FALSE.)
   r = r(ngrid:1:-1)

   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp + 0.5_dp

      ! Prepare grids
      y1 = exp( -rand2(1) * r**2 )
      y2 = exp( -rand2(2) * r**2 )
      call spline(r, y1, size(r), 0.0_dp, 0.0_dp, spline1)
      call spline(r, y2, size(r), 0.0_dp, 0.0_dp, spline2)

      ! The result we get from subroutine kinetic_energy
      call kinetic_energy(nang=10, nshell=75,&
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
   deallocate(d2f2)
end subroutine test_kinetic

end module nao_unit