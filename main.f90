program hallo 
   USE lebedev, ONLY: dp
   USE eddi, ONLY: integration_oneatom
   implicit none

   call integration_oneatom(9, 5)
end program hallo