module eddi
   USE lebedev, ONLY: lebedev_grid,&
                        get_number_of_lebedev_grid,&
                        dp
   USE grid_point, ONLY: type_grid_point
   implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   public :: integration_twoatom

   contains
      subroutine integration_twoatom(nleb, nshell, displacement, integral)
         implicit none
         INTEGER, DIMENSION(2), intent(in) :: nleb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
         INTEGER, DIMENSION(2) :: ileb
         INTEGER :: ngrid
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp) :: integral

         ileb(1) = get_number_of_lebedev_grid(nleb(1))
         ileb(2) = get_number_of_lebedev_grid(nleb(2))
         ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
                 lebedev_grid(ileb(2))%n * nshell(2)

         print *, '!' // REPEAT('-', 78) // '!'
         print *, '! Angular grid1 points: ', lebedev_grid(ileb(1))%n
         print *, '! Radial grid1 points: ', nshell(1)
         print *, '! Angular grid2 points: ', lebedev_grid(ileb(2))%n
         print *, '! Radial grid2 points: ', nshell(2)
         print *, '! Total grid points:   ', ngrid
         print *, '!' // REPEAT('-', 78) // '!'

         allocate(thegrid(ngrid))

         call build_grid(ileb, nshell, displacement, thegrid)
         call evaluate_atgrid(thegrid, integral)
         integral = 4.0_dp*pi*integral

         deallocate(thegrid)

      end subroutine integration_twoatom

      subroutine build_grid(ileb, nshell, displacement, thegrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
         INTEGER :: cnt, iterrad, iterang
         REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
         REAL(KIND=dp) :: alpha

         if (.not. allocated(thegrid)) then
            print *, 'The grid has not been allocated!' !CP
         end if

         cnt = 0

         alpha = pi/REAL(nshell(1)+1, dp)
         do iterang=1, lebedev_grid(ileb(1))%n
            tangw = lebedev_grid(ileb(1))%w(iterang)
            do iterrad=1, nshell(1)
               cnt = cnt+1
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2

               thegrid(cnt)%r = tr*lebedev_grid(ileb(1))%r(:, iterang)
               thegrid(cnt)%weight = tangw * tradw
            enddo !iterrad
         enddo !iterang

         alpha = pi/REAL(nshell(2)+1, dp)
         do iterang=1, lebedev_grid(ileb(2))%n
            tangw = lebedev_grid(ileb(2))%w(iterang)
            do iterrad=1, nshell(2)
               cnt = cnt+1
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2

               thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang) +&
                                 displacement
               thegrid(cnt)%weight = tangw * tradw
            enddo !iterrad
         enddo !iterang
         
      end subroutine build_grid

      subroutine evaluate_atgrid(thegrid, integral)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp) :: norm, f1, f2, partition, integral
         INTEGER :: i

         partition = 0.5_dp
         integral = 0
         do i=1,size(thegrid)
            norm = sqrt(sum((thegrid(i)%r)**2))
            f1 = gaussian(norm)
            f2 = 1.0_dp
            integral = integral + thegrid(i)%weight * partition * f1 * f2
         enddo
      end subroutine evaluate_atgrid

      subroutine print_grid(thegrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER :: i
         do i=1,size(thegrid)
            print *, thegrid(i)%r, thegrid(i)%weight
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

      ! 

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

!! SUBROUTINE GRAVEYARD

      ! Take input N and return an array of N values evenly spaced
      ! between [-lower, upper]
      ! function range(N, lower, upper) result(r)
      !    implicit none
      !    INTEGER :: N, iter
      !    REAL(KIND=dp), DIMENSION(N) :: r
      !    REAL(KIND=dp)               :: spacing, lower, upper
      !    spacing = (upper-lower)/(N-1)
      !    do iter=0,N-1
      !       r(iter+1) = lower + spacing * iter
      !    enddo
      ! end function range

      ! subroutine integration_oneatom(nleb, nshell, integral)
      !    implicit none
      !    INTEGER, intent(in) :: nleb, nshell
      !    INTEGER :: ileb
      !    INTEGER :: ngrid = 0   ! number of points making up the grid
      !    TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      !    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: thefgrid
      !    REAL(KIND=dp) :: integral

      !    ileb = get_number_of_lebedev_grid(nleb)
      !    ngrid = lebedev_grid(ileb)%n * nshell

      !    print *, '!' // REPEAT('-', 78) // '!'
      !    print *, '! Angular grid points: ', lebedev_grid(ileb)%n
      !    print *, '! Radial grid points: ', nshell
      !    print *, '! Total grid points:   ', ngrid
      !    print *, '!' // REPEAT('-', 78) // '!'
      !    allocate(thegrid(ngrid))
      !    allocate(thefgrid(ngrid))

      !    call build_grid(ileb, nshell, thegrid)
      !    call evaluate_atgrid(thegrid, thefgrid)
      !    integral = 4.0_dp*pi*sum(thefgrid)

      !    deallocate(thegrid)
      !    deallocate(thefgrid)

      ! end subroutine integration_oneatom


end module eddi