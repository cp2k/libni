module nao_unit
USE eddi, ONLY: integration_twocenter, integration_threecenter, spline
USE lebedev, ONLY: dp
implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi
contains
subroutine test_twocenter(ntests)
   implicit none
   REAL(KIND=dp) :: dr, integral, ri, mean_rel_error, err
   REAL(KIND=dp), DIMENSION(2) :: rand2
   REAL(KIND=dp), DIMENSION(3) :: rand_pos
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gr1, gy1, gy2, spline1, spline2
   INTEGER :: i, j, ntests

   print *, REPEAT('-', 80)
   print *, REPEAT('-', 30) // ' Testing Two-Center ' // REPEAT('-', 30)
   print *, REPEAT('-', 80)
   print *, ''

   allocate(gr1(200))
   allocate(gy1(200))
   allocate(gy2(200))

   mean_rel_error = 0
   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2 * 5.0_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos)
      rand_pos = rand_pos * sqrt(3.0_dp)

      ! Prepare grids
      dr = 0.1_dp
      do i=1,200
         gr1(i) = i*dr
         gy1(i) = exp(-rand2(1) * gr1(i)**2 )
         gy2(i) = exp(-rand2(2) * gr1(i)**2 )
      enddo
      call spline(gr1, gy1, size(gr1), 0.0_dp, 0.0_dp, spline1)
      call spline(gr1, gy2, size(gr1), 0.0_dp, 0.0_dp, spline2)

      call integration_twocenter(nleb=(/590, 590/), nshell=(/100, 100/), d12=rand_pos, &
                                 gr1=gr1, gy1=gy1, gr2=gr1, gy2=gy2,&
                                 spline1=spline1, spline2=spline2, integral=integral)

      ri = (pi/sum(rand2))**(1.5_dp)
      ri = ri * exp(-rand2(1)*rand2(2)*sum(rand_pos**2)/(sum(rand2)))

      err = abs(1.0_dp-integral/ri)
      mean_rel_error = mean_rel_error+err
      ! print *, 'Is: ', integral
      ! print *, 'Should:', ri
      ! print *, ''
      ! print *, 'Absolute Difference: ', abs(integral-ri)
      print *, 'Exponents: ', rand2
      print *, 'Displacement: ', rand_pos
      print *, 'Relative Error: ', err
      print *, ''
   enddo

   mean_rel_error = mean_rel_error/REAL(ntests, dp)
   print *, 'Mean error in %: ', mean_rel_error*100.0_dp

   print *, REPEAT('-', 80)
   print *, REPEAT('-', 28) // ' End Testing Two-Center ' // REPEAT('-', 28)
   print *, REPEAT('-', 80)
   print *, ''
    
end subroutine test_twocenter

subroutine test_threecenter(ntests)
   implicit none
   REAL(KIND=dp) :: dr, abc, exparg, expargb, integral, ri, mean_rel_error, err
   REAL(KIND=dp), DIMENSION(3) :: rand3
   REAL(KIND=dp), DIMENSION(3) :: rand_pos1
   REAL(KIND=dp), DIMENSION(3) :: rand_pos2
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gr1, gy1, gy2, gy3, spline1, spline2, spline3
   INTEGER :: i, j, ntests

   print *, REPEAT('-', 80)
   print *, REPEAT('-', 30) // ' Testing Three-Center ' // REPEAT('-', 30)
   print *, REPEAT('-', 80)
   print *, ''

   allocate(gr1(200))
   allocate(gy1(200))
   allocate(gy2(200))
   allocate(gy3(200))

   mean_rel_error = 0
   do j=1,ntests
      ! Gaussian exponents
      CALL RANDOM_NUMBER(rand3)
      rand3 = rand3 * 5.0_dp
      ! Displacement
      CALL RANDOM_NUMBER(rand_pos1)
      rand_pos1 = rand_pos1 * sqrt(3.0_dp)
      CALL RANDOM_NUMBER(rand_pos2)
      rand_pos2 = rand_pos2 * sqrt(3.0_dp)

      ! Prepare grids
      dr = 0.1_dp
      do i=1,200
         gr1(i) = i*dr
         gy1(i) = exp(-rand3(1) * gr1(i)**2 )
         gy2(i) = exp(-rand3(2) * gr1(i)**2 )
         gy3(i) = exp(-rand3(3) * gr1(i)**2 )
      enddo
      call spline(gr1, gy1, size(gr1), 0.0_dp, 0.0_dp, spline1)
      call spline(gr1, gy2, size(gr1), 0.0_dp, 0.0_dp, spline2)
      call spline(gr1, gy3, size(gr1), 0.0_dp, 0.0_dp, spline3)

      call integration_threecenter(nleb=(/590, 590, 590/), nshell=(/100, 100, 100/),&
         d12=rand_pos1, d13=rand_pos2, gr1=gr1, gy1=gy1, gr2=gr1, gy2=gy2, gr3=gr1, gy3=gy3,&
         spline1=spline1, spline2=spline2, spline3=spline3, integral=integral)

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

      err = abs(1.0_dp-integral/ri)
      mean_rel_error = mean_rel_error+err
      ! print *, 'Is: ', integral
      ! print *, 'Should:', ri
      ! print *, ''
      ! print *, 'Absolute Difference: ', abs(integral-ri)
      print *, 'Exponents: ', rand3
      print *, 'Displacement: ', rand_pos1
      print *, 'Displacement: ', rand_pos2
      print *, 'Relative Error: ', err
      print *, ''
   enddo

   mean_rel_error = mean_rel_error/REAL(ntests, dp)
   print *, 'Mean error in %: ', mean_rel_error*100.0_dp

   print *, REPEAT('-', 80)
   print *, REPEAT('-', 28) // ' End Testing Two-Center ' // REPEAT('-', 28)
   print *, REPEAT('-', 80)
   print *, ''
    
end subroutine test_threecenter

end module nao_unit