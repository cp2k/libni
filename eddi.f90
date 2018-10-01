module eddi
USE lebedev, ONLY: lebedev_grid,&
                     get_number_of_lebedev_grid,&
                     dp
USE grid, ONLY: type_grid_point, build_onecenter_grid, build_twocenter_grid, build_threecenter_grid
implicit none
REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

type :: type_atom
   REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
   INTEGER :: z = 1
end type type_atom

public :: integration_twocenter, integration_onecenter, integration_threecenter

contains
   subroutine integration_onecenter(nleb, nshell, gr, gy, spline, integral)
      implicit none
      INTEGER, intent(in) :: nleb, nshell
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr, gy, spline
      INTEGER :: ileb
      INTEGER :: ngrid, i
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      REAL(KIND=dp) :: norm, f
      REAL(KIND=dp) :: integral

      ileb = get_number_of_lebedev_grid(n=nleb)
      ngrid = lebedev_grid(ileb)%n * nshell

      ! print *, '!' // REPEAT('-', 78) // '!'
      ! print *, '! Angular grid points: ', lebedev_grid(ileb%n
      ! print *, '! Radial grid points: ', nshell
      ! print *, '!' // REPEAT('-', 78) // '!'
      ! print *, '! Total grid points:   ', ngrid
      ! print *, '!' // REPEAT('-', 78) // '!'

      allocate(thegrid(ngrid))
      call build_onecenter_grid(ileb, nshell, thegrid)

      integral = 0
      do i=1,size(thegrid)
         norm = sqrt(sum((thegrid(i)%r)**2))
         call interpolation(gr, gy, spline, norm, f)

         integral = integral + thegrid(i)%weight * f
      enddo
      ! lebedev integration: 4pi * sum_leb
      integral = 4.0_dp*pi*integral

      deallocate(thegrid)

   end subroutine integration_onecenter
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

   subroutine kinetic_energy(nleb, nshell,&
                             gr1, gy1, gr2, gy2,&
                             spline1, spline2, integral)
      implicit none
      INTEGER, intent(in) :: nleb, nshell
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                              gr2, gy2, &
                                                              spline1, spline2
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rf2, d2rf2, d2rf2_spline
      INTEGER :: ileb
      INTEGER :: ngrid, i
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      REAL(KIND=dp) :: norm, f1, f2
      REAL(KIND=dp) :: integral

      ileb = get_number_of_lebedev_grid(n=nleb)
      ngrid = lebedev_grid(ileb)%n * nshell

      allocate(thegrid(ngrid))
      call build_onecenter_grid(ileb, nshell, thegrid)
      ! Get the spline of the second derivative of r*f2 -> d2f2
      ! < f1 | -0.5*laplace | f2 >

      ! Compute (r*f2)
      rf2 = -0.5_dp*gr2*gy2
      ! Get the 2nd derivative d_r^2(r*f2)
      call spline(gr2, rf2, size(gr2), 0.0_dp, 0.0_dp, d2rf2)
      call spline(gr2, d2rf2, size(gr2), 0.0_dp, 0.0_dp, d2rf2_spline)

      ! Divide by r
      do i=1,size(gr2)
         if(gr2(i) .gt. 0.0_dp) then
            d2rf2(i) = d2rf2(i)/gr2(i)
         else
            d2rf2(i) = d2rf2(i)!0.0_dp
         endif
      enddo
      ! Laplace = 1/r * d_r^2(r*f2) ✓

      integral = 0
      do i=1,size(thegrid)
         norm = sqrt(sum((thegrid(i)%r)**2))
         call interpolation(gr1, gy1, spline1, norm, f1)
         call interpolation(gr2, d2rf2, d2rf2_spline, norm, f2)

         integral = integral + thegrid(i)%weight * f1 * f2
      enddo

      ! laplace = -0.5*🔺
      ! lebedev integration: 4pi * sum_leb
      integral = 4.0_dp*pi*integral

      deallocate(thegrid)
   end subroutine kinetic_energy

   subroutine coulomb_integral(nleb, nshell, r1, y1, r2, y2, s1, s2, d12,&
                               integral)
      implicit none
      INTEGER, DIMENSION(2), intent(in) :: nleb, nshell
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, &
                                                              r2, y2, s1, s2
      REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
      INTEGER, DIMENSION(2) :: ileb
      INTEGER :: ngrid, i, j
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pot
      REAL(KIND=dp) :: norm, f1, f2, integral, potential, r12, rnu1, rnu2

      ileb(1) = get_number_of_lebedev_grid(n=nleb(1))

      ngrid = lebedev_grid(ileb(1))%n * nshell(1)

      allocate(thegrid(ngrid))
      call build_onecenter_grid(ileb(1), nshell(1), thegrid)

      integral = 0
      do i=1,size(thegrid) ! r2
         potential = 0
         do j=1,size(thegrid) ! r1
            rnu1 = sqrt(sum(thegrid(i)%r**2))
            rnu2 = sqrt(sum(thegrid(j)%r**2))
            if (rnu1 .gt. rnu2) then
               r12 = rnu1
            else
               r12 = rnu2
            endif
            norm = sqrt(sum((thegrid(j)%r+d12)**2))
            call interpolation(r1, y1, s1, norm, f1)
            f1 = f1/r12

            potential = potential + thegrid(j)%weight * f1
         enddo

         norm = sqrt(sum((thegrid(j)%r)**2))
         call interpolation(r2, y2, s2, norm, f2)
         integral = integral + thegrid(i)%weight * potential * f2
         if (mod(i, 1000) == 0) then
            print *, integral
            print *, i, '/', size(thegrid)            
         endif
      enddo
      ! lebedev integration: 4pi * sum_leb
      integral = 4.0_dp*pi*integral

      deallocate(thegrid)

       
   end subroutine coulomb_integral

   subroutine read_nfun(fn, gridax, gridf)
      CHARACTER(len=*) :: fn
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: gridax, gridf
      INTEGER :: dim, i

      open(unit=100, file=fn)
      read(100, *) dim
      if (.not. allocated(gridax)) then
         allocate(gridax(dim))
      endif
      if (.not. allocated(gridf)) then
         allocate(gridf(dim))
      endif
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

      if (.TRUE. .or. bound1 .gt. 1e10) then
         yspline(1) = 0
         u(1) = 0
      else
         yspline(1) = -0.5_dp
         u(1) = (3.0_dp/(r(2)-r(1)))*((y(2)-y(1))/(r(2)-r(1))-bound1)
      endif

      do i=2,n-1
         sig = (r(i)-r(i-1))/(r(i+1)-r(i-1))
         p = sig*yspline(i-1)+2.0_dp
         yspline(i) = (sig-1.0_dp)/p

         u(i) = (6.0_dp * ( (y(i+1)-y(i))/(r(i+1)-r(i)) -&
                            (y(i)-y(i-1))/(r(i)-r(i-1)) )/&
                  (r(i+1)-r(i-1)) - sig*u(i-1)) / p
      enddo

      ! addendum: zero second derivative at r->infinity seems reasonable
      !! ignore bound for now and just make it natural
      qn = 0
      un = 0
      !

      yspline(n) = 0
      do i=n-1,1,-1
         yspline(i) = yspline(i)*yspline(i+1)+u(i)
      enddo
   end subroutine spline

   ! Given a function `gy` on a grid `gr` and a requested
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
end module eddi
