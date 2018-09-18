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

         print *, '!' // REPEAT('-', 78) // '!'
         print *, '! Angular grid1 points: ', lebedev_grid(ileb(1))%n
         print *, '! Radial grid1 points: ', nshell(1)
         print *, '! Angular grid2 points: ', lebedev_grid(ileb(2))%n
         print *, '! Radial grid2 points: ', nshell(2)
         print *, '!' // REPEAT('-', 78) // '!'
         print *, '! Total grid points:   ', ngrid
         print *, '!' // REPEAT('-', 78) // '!'

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

      ! thegrid: the integration grid
      ! gr1, gy1: the grid `f1` is given on
      ! gr2, gy2: the grid `f2` is given on
      subroutine build_twocenter_grid(ileb, nshell, displacement, thegrid, &
                                      gr1, gy1, gr2, gy2)
         implicit none
         INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                                 gr2, gy2
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER :: cnt, iterrad, iterang, i, iterileb
         REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
         REAL(KIND=dp) :: alpha, mu, s1, s2, p, ri, rj, R

         cnt = 0
         iterileb = 1
         R = sqrt(sum(displacement**2))

         alpha = pi/REAL(nshell(1)+1, dp)
         do iterang=1, lebedev_grid(ileb(1))%n
            tangw = lebedev_grid(ileb(1))%w(iterang)
            do iterrad=1, nshell(1)
               cnt = cnt+1
               ! COORDINATE
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)

               thegrid(cnt)%r = tr*lebedev_grid(ileb(1))%r(:, iterang)

               ! WEIGHTS
               ! nuclear partition
               ri = sqrt(sum((thegrid(cnt)%r**2)))
               rj = sqrt(sum(((thegrid(cnt)%r-displacement)**2)))

               mu = (ri-rj)/R
               s1 = s3(mu)
               s2 = s3(-mu)
               p = s1/(s1+s2)

               ! radial
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2
   
               thegrid(cnt)%weight = tangw * tradw * p
            enddo !iterrad
         enddo !iterang

         alpha = pi/REAL(nshell(2)+1, dp)
         do iterang=1, lebedev_grid(ileb(2))%n
            tangw = lebedev_grid(ileb(2))%w(iterang)
            do iterrad=1, nshell(2)
               cnt = cnt+1
               ! COORDINATE
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)

               thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang)

               ! WEIGHTS
               ! nuclear partition
               ri = sqrt(sum((thegrid(cnt)%r**2)))
               rj = sqrt(sum(((thegrid(cnt)%r-displacement)**2)))

               mu = (ri-rj)/R
               s1 = s3(mu)
               s2 = s3(-mu)
               p = s2/(s1+s2)

               ! radial
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2
   
               thegrid(cnt)%weight = tangw * tradw * p
            enddo !iterrad
         enddo !iterang
      end subroutine build_twocenter_grid

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

         print *, '!' // REPEAT('-', 78) // '!'
         print *, '! Angular grid1 points: ', lebedev_grid(ileb(1))%n
         print *, '! Radial grid1 points: ', nshell(1)
         print *, '! Angular grid2 points: ', lebedev_grid(ileb(2))%n
         print *, '! Radial grid2 points: ', nshell(2)
         print *, '! Angular grid3 points: ', lebedev_grid(ileb(3))%n
         print *, '! Radial grid3 points: ', nshell(3)
         print *, '!' // REPEAT('-', 78) // '!'
         print *, '! Total grid points:   ', ngrid
         print *, '!' // REPEAT('-', 78) // '!'

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

      ! thegrid: the integration grid
      subroutine build_threecenter_grid(ileb, nshell, d12, d13, thegrid, &
                                        gr1, gy1, gr2, gy2, gr3, gy3)
         implicit none
         INTEGER, DIMENSION(3), intent(in) :: ileb, nshell
         REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                                 gr2, gy2, &
                                                                 gr3, gy3
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         INTEGER :: cnt, iterrad, iterang, i, iterileb
         REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
         REAL(KIND=dp) :: alpha, tP1, tP2, tP3, sP, p1, p2, p3
         REAL(KIND=dp) :: mu12, mu13, mu23, s12, s13, s23, s21, s31, s32
         REAL(KIND=dp) :: r1, r2, r3, R12, R13, R23

         cnt = 0
         iterileb = 1
         R12 = sqrt(sum(d12**2))
         R13 = sqrt(sum(d13**2))
         R23 = sqrt(sum((d13-d12)**2))

         alpha = pi/REAL(nshell(1)+1, dp)
         do iterang=1, lebedev_grid(ileb(1))%n
            tangw = lebedev_grid(ileb(1))%w(iterang)
            do iterrad=1, nshell(1)
               cnt = cnt+1
               ! COORDINATE
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)

               thegrid(cnt)%r = tr*lebedev_grid(ileb(1))%r(:, iterang)

               ! WEIGHTS
               ! nuclear partition
               r1 = sqrt(sum((thegrid(cnt)%r**2)))
               r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
               r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

               mu12 = (r1-r2)/R12
               mu13 = (r1-r3)/R13
               mu23 = (r2-r3)/R23

               s12 = s3(mu12)
               s21 = s3(-mu12)
               s13 = s3(mu13)
               s31 = s3(-mu13)
               s32 = s3(-mu23)
               s23 = s3(mu23)

               tP1 = s12*s13
               tP2 = s21*s23
               tP3 = s31*s32

               sP = tP1+tP2+tP3

               p1 = tP1/sP
               p2 = tP2/sP
               p3 = tP3/sP

               ! radial
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2
   
               thegrid(cnt)%weight = tangw * tradw * p1
            enddo !iterrad
         enddo !iterang

         alpha = pi/REAL(nshell(2)+1, dp)
         do iterang=1, lebedev_grid(ileb(2))%n
            tangw = lebedev_grid(ileb(2))%w(iterang)
            do iterrad=1, nshell(2)
               cnt = cnt+1
               ! COORDINATE
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)

               thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang)

               ! WEIGHTS
               ! nuclear partition
               r1 = sqrt(sum((thegrid(cnt)%r**2)))
               r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
               r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

               mu12 = (r1-r2)/R12
               mu13 = (r1-r3)/R13
               mu23 = (r2-r3)/R23

               s12 = s3(mu12)
               s21 = s3(-mu12)
               s13 = s3(mu13)
               s31 = s3(-mu13)
               s32 = s3(-mu23)
               s23 = s3(mu23)

               tP1 = s12*s13
               tP2 = s21*s23
               tP3 = s31*s32

               sP = tP1+tP2+tP3

               p1 = tP1/sP
               p2 = tP2/sP
               p3 = tP3/sP

               ! radial
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2
   
               thegrid(cnt)%weight = tangw * tradw * p2
            enddo !iterrad
         enddo !iterang

         alpha = pi/REAL(nshell(3)+1, dp)
         do iterang=1, lebedev_grid(ileb(3))%n
            tangw = lebedev_grid(ileb(3))%w(iterang)
            do iterrad=1, nshell(3)
               cnt = cnt+1
               ! COORDINATE
               targ = REAL(iterrad, dp)*alpha
               tcos = cos(targ)
               tr = radial_mapping(tcos)

               thegrid(cnt)%r = tr*lebedev_grid(ileb(3))%r(:, iterang)

               ! WEIGHTS
               ! nuclear partition
               r1 = sqrt(sum((thegrid(cnt)%r**2)))
               r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
               r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

               mu12 = (r1-r2)/R12
               mu13 = (r1-r3)/R13
               mu23 = (r2-r3)/R23

               s12 = s3(mu12)
               s21 = s3(-mu12)
               s13 = s3(mu13)
               s31 = s3(-mu13)
               s32 = s3(-mu23)
               s23 = s3(mu23)

               tP1 = s12*s13
               tP2 = s21*s23
               tP3 = s31*s32

               sP = tP1+tP2+tP3

               p1 = tP1/sP
               p2 = tP2/sP
               p3 = tP3/sP

               ! radial
               tsin = alpha*sin(targ)**2
               tradw = tsin/sqrt(1.0_dp-tcos**2)
               tradw = 2.0_dp*tradw*tr**2/(1.0_dp-tcos)**2
   
               thegrid(cnt)%weight = tangw * tradw * p3
            enddo !iterrad
         enddo !iterang
      end subroutine build_threecenter_grid


      function h(mu)
         implicit none
         REAL(KIND=dp) :: mu, h
         h = 1.5_dp*mu-0.5_dp*mu**3
      end function h

      function s3(mu)
         implicit none
         REAL(KIND=dp) :: mu, s3
         s3 = 0.5_dp*(1-h(h(h(mu))))
      end function s3

      subroutine grid_parameters(atom, nleb, nshell)
         implicit none
         INTEGER :: atom, nleb, nshell

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

      subroutine spline(r, y, n, bound1, boundn, yspline)
         implicit none
         INTEGER, INTENT(in) :: n
         REAL(KIND=dp), DIMENSION(1:n), INTENT(in) :: r, y
         REAL(KIND=dp), INTENT(in) :: bound1, boundn
         REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: yspline, u

         INTEGER :: i
         REAL(KIND=dp) :: sig, p, un, qn
         
         allocate(yspline(n))
         allocate(u(n))

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