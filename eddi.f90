module eddi
   USE lebedev, ONLY: lebedev_grid,&
                        get_number_of_lebedev_grid,&
                        dp
   USE grid, ONLY: type_grid_point, build_twocenter_grid, build_threecenter_grid
   implicit none
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   type :: type_atom
      REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
      INTEGER :: z = 1
   end type type_atom

   public :: integration_twocenter

   contains
      subroutine integration_twocenter(nleb, nshell, d12,&
                                       gr1, gy1, gr2, gy2,&
                                       spline1, spline2, integral)
         implicit none
         INTEGER, DIMENSION(2), intent(in) :: nleb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                                 gr2, gy2, &
                                                                 spline1, spline2
         INTEGER, DIMENSION(2) :: ileb
         INTEGER :: ngrid, i
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp) :: norm, f1, f2
         REAL(KIND=dp) :: integral

         ileb(1) = get_number_of_lebedev_grid(n=nleb(1))
         ileb(2) = get_number_of_lebedev_grid(n=nleb(2))
         ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
                 lebedev_grid(ileb(2))%n * nshell(2)

         ! print *, '!' // REPEAT('-', 78) // '!'
         ! print *, '! Angular grid1 points: ', lebedev_grid(ileb(1))%n
         ! print *, '! Radial grid1 points: ', nshell(1)
         ! print *, '! Angular grid2 points: ', lebedev_grid(ileb(2))%n
         ! print *, '! Radial grid2 points: ', nshell(2)
         ! print *, '!' // REPEAT('-', 78) // '!'
         ! print *, '! Total grid points:   ', ngrid
         ! print *, '!' // REPEAT('-', 78) // '!'

         allocate(thegrid(ngrid))
         call build_twocenter_grid(ileb, nshell, d12, thegrid, &
                                   gr1, gy1, gr2, gy2)

         integral = 0
         do i=1,size(thegrid)
            norm = sqrt(sum((thegrid(i)%r)**2))
            call interpolation(gr1, gy1, spline1, norm, f1)

            norm = sqrt(sum((thegrid(i)%r+d12)**2))
            call interpolation(gr2, gy2, spline2, norm, f2)

            integral = integral + thegrid(i)%weight * f1 * f2
         enddo
         ! lebedev integration: 4pi * sum_leb
         integral = 4.0_dp*pi*integral

         deallocate(thegrid)

      end subroutine integration_twocenter

      subroutine integration_threecenter(nleb, nshell, d12, d13, &
                                         gr1, gy1, gr2, gy2, gr3, gy3,&
                                         spline1, spline2, spline3, integral)
         implicit none
         INTEGER, DIMENSION(3), intent(in) :: nleb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                         gr2, gy2, gr3, gy3, spline1, spline2, spline3
         INTEGER, DIMENSION(3) :: ileb
         INTEGER :: ngrid, i
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp) :: norm, f1, f2, f3
         REAL(KIND=dp) :: integral

         ileb(1) = get_number_of_lebedev_grid(n=nleb(1))
         ileb(2) = get_number_of_lebedev_grid(n=nleb(2))
         ileb(3) = get_number_of_lebedev_grid(n=nleb(3))
         ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
                 lebedev_grid(ileb(2))%n * nshell(2) + &
                 lebedev_grid(ileb(3))%n * nshell(3)

         ! print *, '!' // REPEAT('-', 78) // '!'
         ! print *, '! Angular grid1 points: ', lebedev_grid(ileb(1))%n
         ! print *, '! Radial grid1 points: ', nshell(1)
         ! print *, '! Angular grid2 points: ', lebedev_grid(ileb(2))%n
         ! print *, '! Radial grid2 points: ', nshell(2)
         ! print *, '! Angular grid3 points: ', lebedev_grid(ileb(3))%n
         ! print *, '! Radial grid3 points: ', nshell(3)
         ! print *, '!' // REPEAT('-', 78) // '!'
         ! print *, '! Total grid points:   ', ngrid
         ! print *, '!' // REPEAT('-', 78) // '!'

         allocate(thegrid(ngrid))
         call build_threecenter_grid(ileb, nshell, d12, d13, thegrid, &
                                     gr1, gy1, gr2, gy2, gr3, gy3)

         integral = 0
         do i=1,size(thegrid)
            norm = sqrt(sum((thegrid(i)%r)**2))
            call interpolation(gr1, gy1, spline1, norm, f1)

            norm = sqrt(sum((thegrid(i)%r+d12)**2))
            call interpolation(gr2, gy2, spline2, norm, f2)

            norm = sqrt(sum((thegrid(i)%r+d13)**2))
            call interpolation(gr3, gy3, spline3, norm, f3)

            integral = integral + thegrid(i)%weight * f1 * f2 * f3
         enddo
         ! lebedev integration: 4pi * sum_leb
         integral = 4.0_dp*pi*integral

         deallocate(thegrid)

      end subroutine integration_threecenter

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

      subroutine spline(r, y, n, bound1, boundn, yspline)
         implicit none
         INTEGER, INTENT(in) :: n
         REAL(KIND=dp), DIMENSION(1:n), INTENT(in) :: r, y
         REAL(KIND=dp), INTENT(in) :: bound1, boundn
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: yspline, u

         INTEGER :: i
         REAL(KIND=dp) :: sig, p, un, qn
         
         if (.not. allocated(yspline)) then
            allocate(yspline(n))
         endif
         if (.not. allocated(u)) then
            allocate(u(n))
         endif

         ! ignore bound for now and just make it natural
         yspline(1) = 0
         u(1) = 0
         !

         do i=2,n-1
            sig = (r(i)-r(i-1))/(r(i+1)-r(i-1))
            p = sig*yspline(i-1)+2.0_dp
            yspline(i) = (sig-1.0_dp)/p

            u(i) = (6.0_dp * ( (y(i+1)-y(i))/(r(i+1)-r(i)) -&
                               (y(i)-y(i-1))/(r(i)-r(i-1)) )/&
                     (r(i+1)-r(i-1)) - sig*u(i-1)) / p
         enddo

         ! ignore bound for now and just make it natural
         qn = 0
         un = 0
         !

         yspline(n) = 0
         do i=n-1,1,-1
            yspline(i) = yspline(i)*yspline(i+1)+u(i)
         enddo
      end subroutine spline

      ! Given a function `gridy` on a grid `gridr` and a requested
      ! function value y(r) interpolates the function value `y`
      subroutine interpolation(gr, gy, spline, r, y)
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr, gy
         REAL(KIND=dp), DIMENSION(1:size(gr)) :: spline
         REAL(KIND=dp), intent(in) :: r
         REAL(KIND=dp) :: y

         INTEGER :: low, upper, mid
         REAL(KIND=dp) :: A, B, h
         ! find the closest grid point by bisection
         low = 0
         upper = size(gr)
         do while (upper-low .gt. 1)
            mid = NINT((low+upper)/2.0_dp)
            if (gr(mid) .gt. r) then
               upper = mid
            else
               low = mid
            endif
         enddo
         if (gr(upper) .eq. r) then
            y = gy(upper)
         else if (gr(low) .eq. r) then
            y = gy(low)
         else if ((gr(upper) .gt. r) .and. (gr(low) .lt. r)) then
            ! LINEAR INTERPOLATION
            ! A = (gr(upper)-r)/(gr(upper)-gr(low))
            ! y = A*gy(low) + (1.0_dp-A)*gy(upper)
            ! SPLINE INTERPOLATION
            h = gr(upper)-gr(low)
            A = (gr(upper)-r)/h
            B = (r-gr(low))/h
            y = A*gy(low) + B*gy(upper) + &
                  ((A**3-A)*spline(low)+(B**3-B)*spline(upper)) * (h**2)/6.0_dp
         else if (gr(upper) .lt. r) then
            y = gy(upper)
            ! print *, 'Extrapolation!'
            ! print *, upper, gr(upper), r, y
         else if (gr(low) .gt. r) then
            y = gy(low)
            ! print *, 'Extrapolation!'
         endif
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
end module eddi
