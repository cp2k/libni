module eddi
   USE lebedev, ONLY: lebedev_grid,&
                        get_number_of_lebedev_grid,&
                        dp
   USE grid_point, ONLY: type_grid_point
   implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   type :: type_atom
      REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
      INTEGER :: z = 1
   end type type_atom

   public :: integration_twocenter

   contains
      subroutine integration_twocenter(nleb, nshell, displacement, &
                                       gr1, gy1, gr2, gy2, integral)
         implicit none
         INTEGER, DIMENSION(2), intent(in) :: nleb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                                 gr2, gy2
         INTEGER, DIMENSION(2) :: ileb
         INTEGER :: ngrid
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp) :: integral

         ileb(1) = get_number_of_lebedev_grid(n=nleb(1))
         ileb(2) = get_number_of_lebedev_grid(n=nleb(2))
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

         call integrate(ileb, nshell, displacement, thegrid,  &
                        gr1, gy1, gr2, gy2, integral)

         deallocate(thegrid)

      end subroutine integration_twocenter

      ! thegrid: the integration grid
      ! gr1, gy1: the grid `f1` is given on
      ! gr2, gy2: the grid `f2` is given on
      subroutine integrate(ileb, nshell, displacement, thegrid, &
                           gr1, gy1, gr2, gy2, integral)
         implicit none
         INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                                 gr2, gy2
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER :: cnt, iterrad, iterang, i
         REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
         REAL(KIND=dp) :: norm, f1, f2, f3, integral, distance, alpha

         cnt = 0
         distance = sqrt(sum(displacement**2))

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
               tradw = tr**2 * tradw * 2.0_dp/(1.0_dp-tcos)**2

               thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang) +&
                                 displacement
               thegrid(cnt)%weight = tangw * tradw
            enddo !iterrad
         enddo !iterang

         integral = 0
         do i=1,size(thegrid)
            norm = sqrt(sum((thegrid(i)%r)**2))
            !f1 = gaussian(norm)
            call interpolation(gr1, gy1, norm, f1)
            norm = sqrt(sum((thegrid(i)%r+displacement)**2))
            ! f2 = gaussian(norm, 0.5_dp)
            call interpolation(gr2, gy2, norm, f2)
            f3 = 1.0_dp
            integral = integral + thegrid(i)%weight * f1 * f3 * f2
         enddo
         integral = 4.0_dp*pi*integral * 0.5_dp

      end subroutine integrate

      subroutine grid_parameters(atom, nleb, nshell)
         implicit none
         INTEGER :: atom, nleb, nshell

         !nshell = ceiling(24.0_dp+2.0_dp*atoms/5)

         SELECT CASE (atom)
         CASE (0)
            nshell = 5
            nleb = 10
         CASE (1:2)
            nshell = 35
            nleb = 302
         CASE (3:10)
            nshell = 40
            nleb = 590
         CASE (11:18)
            nshell = 45
            nleb = 590
         CASE (19:)
            nshell = 50
            nleb = 590
         CASE DEFAULT
            nshell = 10
            nleb = 302
         END SELECT
      end subroutine grid_parameters

      subroutine print_grid(thegrid)
         implicit none
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER :: i
         do i=1,size(thegrid)
            print *, thegrid(i)%r, thegrid(i)%weight
         enddo 
      end subroutine print_grid

      subroutine read_nfun(fn, gridax, gridf)
         CHARACTER(len=*) :: fn
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gridax, gridf
         INTEGER :: dim, i

         open(unit=100, file=fn)
         read(100, *) dim
         allocate(gridax(dim))
         allocate(gridf(dim))
         do i=1, dim
            read(100, *) gridax(i), gridf(i)
         enddo
         close(100)
      end subroutine read_nfun

      ! Given a function `gridy` on a grid `gridr` and a requested
      ! function value y(r) interpolates the function value `y`
      subroutine interpolation(gridr, gridy, r, y)
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gridr, gridy
         REAL(KIND=dp), intent(in) :: r
         REAL(KIND=dp) :: y

         INTEGER :: low, high, mid
         REAL(KIND=dp) :: A, B
         ! find the closest grid point by bisection
         low = 0
         high = size(gridr)
         do while (high-low .gt. 1)
            mid = NINT((low+high)/2.0_dp)
            if (gridr(mid) .gt. r) then
               high = mid
            else
               low = mid
            endif
         enddo
         A = (gridr(high)-r)/(gridr(high)-gridr(low))
         y = A*gridy(low) + (1.0_dp-A)*gridy(high)

      end subroutine interpolation

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

      ! partition_0
      function p0(x) result(p)
         implicit none
         REAL(KIND=dp) :: x, p
         p = 2.0_dp*x**3 - 3.0_dp*x**2 + 1.0_dp
      end function p0

      function partition(x) result(p)
         implicit none
         REAL(KIND=dp) :: x, p
         if (x .gt. 1) then
            p = 0.0_dp
         else
            p = p0(p0(p0(x)))
         endif
      end function partition

end module eddi