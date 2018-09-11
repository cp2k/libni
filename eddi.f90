module eddi
   USE lebedev, ONLY: lebedev_grid,&
                        get_number_of_lebedev_grid,&
                        dp
   USE grid_point, ONLY: type_grid_point
   implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   public :: integration_oneatom

   contains
      subroutine integration_oneatom(nleb, nshell)
         implicit none
         INTEGER, intent(in) :: nleb, nshell
         INTEGER :: ileb
         INTEGER :: ngrid = 0   ! number of points making up the grid
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: thefgrid

         ileb = get_number_of_lebedev_grid(nleb)
         ngrid = lebedev_grid(ileb)%n * nshell

         print *, '!' // REPEAT('-', 79)
         print *, '! Angular grid points: ', lebedev_grid(ileb)%n
         print *, '! Radial grid points: ', nshell
         print *, '! Total grid points:   ', ngrid
         print *, '!' // REPEAT('-', 79)
         allocate(thegrid(ngrid))
         allocate(thefgrid(ngrid))

         call build_grid(ileb, nshell, thegrid)
         call evaluate_atgrid(thegrid, thefgrid)
         ! call print_grid(thegrid, thefgrid)
         print *, sum(thefgrid)

         deallocate(thegrid)
         deallocate(thefgrid)

      end subroutine integration_oneatom

      subroutine evaluate_atgrid(thegrid, thefgrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: thefgrid
         REAL(KIND=dp) :: norm
         INTEGER :: i

         do i=1,size(thefgrid)
            norm = sqrt(sum((thegrid(i)%r)**2))
            thefgrid(i) = thegrid(i)%weight * gaussian(norm)
         enddo 
      end subroutine evaluate_atgrid

      subroutine build_grid(ileb, nshell, thegrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER, intent(in) :: ileb, nshell
         INTEGER :: cnt, iterrad, iterang
         REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ
         REAL(KIND=dp) :: alpha
         REAL(KIND=dp), DIMENSION(1:nshell) :: rrange

         if (.not. allocated(thegrid)) then
            print *, 'Allocating the grid in build_grid.' !CP
            allocate(thegrid(lebedev_grid(ileb)%n * nshell))
         end if

         alpha = pi/REAL(nshell, dp)
         rrange = range(nshell, -0.99_dp, +0.99_dp)
         cnt = 0
         do iterang=1, lebedev_grid(ileb)%n
            tangw = lebedev_grid(ileb)%w(iterang)
            do iterrad=1, nshell
               cnt = cnt+1
               targ = alpha*(2.0_dp*iterrad-1)/2.0_dp ! pi*(2i-1)/2n
               tsin = sin(targ)

               thegrid(cnt)%r = radial_mapping(rrange(iterrad))*lebedev_grid(ileb)%r(:, iterang)
               thegrid(cnt)%weight = lebedev_grid(ileb)%w(iterang) * tsin
            enddo !iterrad
         enddo !iterang
         
      end subroutine build_grid

      subroutine print_grid(thegrid, thefgrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: thefgrid
         INTEGER :: i
         do i=1,size(thegrid)
            print *, thegrid(i)%r, thegrid(i)%weight, thefgrid(i)
         enddo 
      end subroutine print_grid

      ! take values x in [-1, 1) and map them to (0, +infty)
      function radial_mapping(x) result(r)
         implicit none
         REAL(KIND=dp), INTENT(in) :: x
         REAL(KIND=dp) :: eps, alpha, r

         ! eps = 1.0_dp
         ! alpha = 0.6_dp
         ! r = eps/log(2.0_dp)*(1+x)**alpha*log(2/(1-x))
         r = (1+x)/(1-x)
          
      end function radial_mapping

      ! Evaluate a gaussian e^(-a.x^2)
      function gaussian(x, a) result(g)
         implicit none
         REAL(KIND=dp) :: x, g
         REAL(KIND=dp), OPTIONAL :: a
         if (PRESENT(a)) then
            g = exp(-a*x**2)
         else
            g = exp(-x**2)
         endif

      end function gaussian

      ! Take input N and return an array of N values evenly spaced
      ! between [-lower, upper]
      function range(N, lower, upper) result(r)
         implicit none
         INTEGER :: N, iter
         REAL(KIND=dp), DIMENSION(N) :: r
         REAL(KIND=dp)               :: spacing, lower, upper
         spacing = (upper-lower)/(N-1)
         do iter=0,N-1
            r(iter+1) = lower + spacing * iter
         enddo
      end function range

end module eddi