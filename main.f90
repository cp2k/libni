program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: type_atom, integration_twocenter, &
                   read_nfun, pi, interpolation, spline, &
                   integration_threecenter
   USE grid, ONLY: grid_parameters
   implicit none

   TYPE(type_atom), DIMENSION(3) :: atoms
   INTEGER, DIMENSION(:), ALLOCATABLE :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: d12, d13
   REAL(KIND=dp) :: integral
   REAL(KIND=dp) :: start, finish, del
   INTEGER :: i
   REAL(KIND=dp) :: y

   CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
   CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian_alpha0_5.grid'
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gr1, gr2, gy1, gy2, gr3, gy3,&
                                               spline1, spline2, spline3

   ! Build parameters
   call cpu_time(start)
   atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
   atoms(1)%z = 10
   atoms(2)%r = (/ 1.0_dp, 1.0_dp, 1.0_dp /)
   atoms(2)%z = 10
   atoms(3)%r = (/ 0.0_dp, 2.0_dp, 0.0_dp /)
   atoms(3)%z = 10

   d12 = atoms(2)%r - atoms(1)%r
   d13 = atoms(3)%r - atoms(1)%r

   call read_nfun(fn1, gr1, gy1)
   call spline(gr1, gy1, size(gr1), 0.0_dp, 0.0_dp, spline1)
   call read_nfun(fn2, gr2, gy2)
   call spline(gr2, gy2, size(gr2), 0.0_dp, 0.0_dp, spline2)
   call read_nfun(fn2, gr3, gy3)
   call spline(gr3, gy3, size(gr3), 0.0_dp, 0.0_dp, spline3)

   call cpu_time(finish)
   print *, 'took', finish-start
   ! 2 center
   print *, '!' // REPEAT('-', 78) // '!'
   print *, '! TWO CENTER !'
   print *, '!' // REPEAT('-', 78) // '!'
   call cpu_time(start)
   allocate(nleb(2))
   allocate(nshell(2))
   call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   call grid_parameters(atoms(2)%z, nleb(2), nshell(2))
   call integration_twocenter(nleb, nshell, d12, &
                              gr1, gy1, gr2, gy2, spline1, spline2, integral)
   print *, integral!, ',', pi**1.5_dp-integral
   deallocate(nleb)
   deallocate(nshell)
   call cpu_time(finish)
   print *, 'took', finish-start

   ! 3 center
   print *, '!' // REPEAT('-', 78) // '!'
   print *, '! THREE CENTER !'
   print *, '!' // REPEAT('-', 78) // '!'
   call cpu_time(start)
   allocate(nleb(3))
   allocate(nshell(3))
   call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   call grid_parameters(atoms(2)%z, nleb(2), nshell(2))
   call grid_parameters(atoms(3)%z, nleb(3), nshell(3))
   call integration_threecenter(nleb, nshell, d12, d13,&
                                gr1, gy1, gr2, gy2, gr3, gy3,&
                                spline1, spline2, spline3, integral)
   print *, integral
   deallocate(nleb)
   deallocate(nshell)
   call cpu_time(finish)
   print *, 'took', finish-start

end program hallo