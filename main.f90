program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: integration_twoatom
   implicit none
   INTEGER, DIMENSION(2) :: nleb, nshell
   REAL(KIND=dp), DIMENSION(3) :: displacement
   REAL(KIND=dp) :: integral

   nleb = (/ 19, 19/)
   nshell = (/50, 50/)
   displacement = (/3.0_dp, 0.0_dp, 0.0_dp/)
   call integration_twoatom(nleb, nshell, displacement, integral)
   print *, integral
end program hallo