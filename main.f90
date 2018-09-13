program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: type_atom, integration_twocenter, grid_parameters, &
                   read_nfun, pi, interpolation
   implicit none

   TYPE(type_atom), DIMENSION(2) :: atoms
   INTEGER, DIMENSION(2) :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: displacement
   REAL(KIND=dp) :: integral
   REAL(KIND=dp) :: start, finish, del
   INTEGER :: i
   REAL(KIND=dp) :: y

   CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
   CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian_alpha0_5.grid'
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gridax1, gridax2, gridf1, gridf2

   call cpu_time(start)
   atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
   atoms(1)%z = 10
   atoms(2)%r = (/ 1.0_dp, 1.0_dp, 1.0_dp /)
   atoms(2)%z = 10

   displacement = atoms(2)%r - atoms(1)%r

   call read_nfun(fn1, gridax1, gridf1)
   call read_nfun(fn2, gridax2, gridf2)

   call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   call grid_parameters(atoms(2)%z, nleb(2), nshell(2))
   call integration_twocenter(nleb, nshell, displacement, &
                              gridax1, gridf1, gridax2, gridf2, integral)
   print *, integral!, ',', pi**1.5_dp-integral
   call cpu_time(finish)
   ! 1.1150509804006377 in 14ms
   ! Python: 1.1150476803117528 in 2min 37s

   print *, 'took', finish-start

end program hallo