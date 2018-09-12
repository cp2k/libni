program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: integration_twoatom
   implicit none
   INTEGER, DIMENSION(2) :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: displacement
   REAL(KIND=dp) :: integral, time

   nleb = (/ 19, 19/)
   nshell = (/30, 30/)
   displacement = (/3.0_dp, 0.0_dp, 0.0_dp/)
   call integration_twoatom(nleb, nshell, displacement, integral)
   print *, 5.568327997
   print *, integral
   print *, 5.568327997-integral
   ! Should: 5.568327997
   ! Is:     5.5683279968305239
end program hallo