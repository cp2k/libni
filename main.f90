program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: type_atom, integration_twocenter, grid_parameters, pi
   implicit none

   TYPE(type_atom), DIMENSION(2) :: atoms
   INTEGER, DIMENSION(2) :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: displacement
   REAL(KIND=dp) :: integral, time

   atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
   atoms(1)%z = 10
   atoms(2)%r = (/ 1.0_dp, 1.0_dp, 1.0_dp /)
   atoms(2)%z = 10

   displacement = atoms(2)%r - atoms(1)%r

   call grid_parameters(atoms(1)%z, nleb(1), nshell(1))
   call grid_parameters(atoms(2)%z, nleb(2), nshell(2))
   call integration_twocenter(nleb, nshell, displacement, integral)
   print *, integral!, ',', pi**1.5_dp-integral
   ! Should: 5.568327997
   ! Is:     5.5683279968305239
end program hallo