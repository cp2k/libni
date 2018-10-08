program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: type_atom, integration_twocenter, &
                   read_nfun, pi, interpolation, spline, &
                   integration_threecenter, kinetic_energy, &
                   integration_onecenter, coulomb_integral, &
                   radial_integration
   USE grid, ONLY: grid_parameters
   USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic
   implicit none

   TYPE(type_atom), DIMENSION(3) :: atoms
   INTEGER, DIMENSION(:), ALLOCATABLE :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: d12, d13
   REAL(KIND=dp) :: integral, integral_low
   REAL(KIND=dp) :: start, finish, del, eps
   INTEGER :: i
   REAL(KIND=dp) :: y

   CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
   CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gr1, gr2, gy1, gy2, gr3, gy3,&
                                               spline1, spline2, spline3


! ! Get the gaussian on a useful grid
! allocate(gy1(5000))
! allocate(gr1(5000))

! del = 0.001_dp
! integral_low = 0.0_dp
! do i=1,size(gy1)
!    gr1(i) = REAL(i-1, dp)*del
!    gy1(i) = exp(-gr1(i)**2)
!    integral_low = integral_low + del*gy1(i)
! enddo
! call radial_integration(f=gy1, r=gr1, n=50, integral=integral)
! print *, integral, integral_low/2.0_dp
! print *, abs(1-integral/(sqrt(pi)/4))

! deallocate(gy1)
! deallocate(gr1)
! return

! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=50 , loud=.FALSE.)
call test_kinetic(ntests=100, loud=.FALSE.)
return
   ! ! Mean error in %:   0.41358588140703562
   ! return

   ! Build parameters
   atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
   atoms(1)%z = 1
   atoms(2)%r = (/ 100.0_dp, 0.0_dp, 0.0_dp /)
   atoms(2)%z = 1
   atoms(3)%r = (/ 0.0_dp, 2.0_dp, 0.0_dp /)
   atoms(3)%z = 1

   d12 = atoms(2)%r - atoms(1)%r
   d13 = atoms(3)%r - atoms(1)%r

   call read_nfun(fn1, gr1, gy1)
   call spline(gr1, gy1, size(gr1), 0.0_dp, 0.0_dp, spline1)
   call read_nfun(fn2, gr2, gy2)
   call spline(gr2, gy2, size(gr2), 0.0_dp, 0.0_dp, spline2)
   call read_nfun(fn2, gr3, gy3)
   call spline(gr3, gy3, size(gr3), 0.0_dp, 0.0_dp, spline3)

   allocate(nleb(2))
   allocate(nshell(2))
   call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   call grid_parameters(atoms(2)%z, nleb(2), nshell(2))

   ! Coulomb integral
   print *, '!' // REPEAT('-', 78) // '!'
   print *, '! Coulomb integral !'
   print *, '!' // REPEAT('-', 78) // '!'
   call coulomb_integral(nleb, nshell, gr1, gy1, gr2, gy2, spline1, spline2, d12,&
                         integral)
   print *, integral

   ! call integration_onecenter(590, 100, gr1, gy1, spline1, integral)
   ! print *, integral
   ! return

   ! ! Kinetic energy
   ! print *, '!' // REPEAT('-', 78) // '!'
   ! print *, '! Kinetic energy !'
   ! print *, '!' // REPEAT('-', 78) // '!'
   ! call cpu_time(start)
   ! call kinetic_energy(nleb(1), nshell(1), gr1, gy1, gr2, gy2,&
   !                     spline1, spline2, integral)
   ! print *, integral/(4.0_dp*pi)
   ! call cpu_time(finish)
   ! print *, 'took', finish-start
   ! ! End Kinteic energy
   ! ! return 

   ! ! 2 center
   ! print *, '!' // REPEAT('-', 78) // '!'
   ! print *, '! TWO CENTER !'
   ! print *, '!' // REPEAT('-', 78) // '!'
   ! call cpu_time(start)
   ! call integration_twocenter(nleb, nshell, d12, &
   !                            gr1, gy1, gr2, gy2, spline1, spline2, integral)
   ! print *, integral/(4.0_dp*pi)!, ',', pi**1.5_dp-integral
   ! deallocate(nleb)
   ! deallocate(nshell)
   ! call cpu_time(finish)
   ! print *, 'took', finish-start

   ! ! ! 3 center
   ! ! print *, '!' // REPEAT('-', 78) // '!'
   ! ! print *, '! THREE CENTER !'
   ! ! print *, '!' // REPEAT('-', 78) // '!'
   ! ! call cpu_time(start)
   ! ! allocate(nleb(3))
   ! ! allocate(nshell(3))
   ! ! call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   ! ! call grid_parameters(atoms(2)%z, nleb(2), nshell(2))
   ! ! call grid_parameters(atoms(3)%z, nleb(3), nshell(3))
   ! ! call integration_threecenter(nleb, nshell, d12, d13,&
   ! !                              gr1, gy1, gr2, gy2, gr3, gy3,&
   ! !                              spline1, spline2, spline3, integral)
   ! ! print *, integral
   ! ! call cpu_time(finish)
   ! ! print *, 'took', finish-start
   deallocate(nleb)
   deallocate(nshell)

end program hallo