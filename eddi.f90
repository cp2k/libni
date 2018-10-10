module eddi
USE lebedev, ONLY: lebedev_grid,&
                   get_number_of_lebedev_grid,&
                   dp
USE grid, ONLY: build_onecenter_grid, build_twocenter_grid, build_threecenter_grid, &
                type_grid_point, radial_grid
implicit none
REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

type :: type_atom
   REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
   INTEGER :: z = 1
end type type_atom

public :: integration_twocenter, integration_onecenter, integration_threecenter,&
           radial_integration

contains
! **********************************************
!> \brief Computes the radial integral of f(r)
!> \param f(n): The tabulated function at n grid points
!> \param r(n): The tabulated grid points
!> \param n: The number of radial grid points
!> \param integral: The integral's value
!> \author 
! **********************************************
subroutine radial_integration(f, r, n, integral)
   implicit none
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: f, r
   INTEGER, intent(in)                                  :: n
   REAL(KIND=dp) :: integral
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rad, wr, d2f, fun
   INTEGER :: i

   allocate(rad(n))
   allocate(wr(n))
   allocate(d2f(n))
   allocate(fun(n))

   integral = 0.0_dp
   
   ! Put the radial grid points into `rad` and their weights into `wr`
   call radial_grid(r=rad, wr=wr, n=n, addr2=.TRUE.)

   ! Create the spline
   call spline(r=r, y=f, n=size(r), bound1=0.0_dp, boundn=0.0_dp, yspline=d2f)
   
   ! Sum over all radial grid points
   do i=1,n
      call interpolation(gr=r, gy=f, spline=d2f, r=rad(i), y=fun(i))
   enddo
   integral = sum( wr * fun )

   deallocate(rad)
   deallocate(wr)
   deallocate(d2f)
   deallocate(fun)
end subroutine radial_integration

subroutine integration_onecenter(nang, nshell, r, y, spline, integral)
   implicit none
   ! Input
   INTEGER, intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r, y, spline
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   INTEGER :: ileb, ngrid, i
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: int_i
   REAL(KIND=dp) :: norm
   ! End header

   ileb = get_number_of_lebedev_grid(n=nang)
   ngrid = lebedev_grid(ileb)%n * nshell
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(int_i(ngrid))
   int_i = 0.0_dp

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE.,&
                             grid_r=grid_r, grid_w=grid_w)

   do i=1,size(grid_w)
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r, y, spline, norm, int_i(i))
   enddo

   integral = sum(grid_w*int_i)

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(int_i)
end subroutine integration_onecenter

subroutine integration_twocenter(nang, nshell, d12, r1, y1, r2, y2, &
                                 spline1, spline2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, r2, y2, &
                                                           spline1, spline2
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f1, f2
   REAL(KIND=dp) :: norm
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid, i

   ileb(1) = get_number_of_lebedev_grid(n=nang(1))
   ileb(2) = get_number_of_lebedev_grid(n=nang(2))

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(f1(ngrid))
   allocate(f2(ngrid))

   call build_twocenter_grid(ileb=ileb, nshell=nshell, displacement=d12, &
                             addr2=.TRUE., grid_r=grid_r, grid_w=grid_w)

   do i=1,size(grid_w)
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))

      norm = sqrt(sum( (grid_r(i, :) + d12 )**2 ))
      call interpolation(r2, y2, spline2, norm, f2(i))
   enddo
   integral = sum(grid_w * f1*f2 )

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(f1)
   deallocate(f2)

end subroutine integration_twocenter

subroutine integration_threecenter(nang, nshell, d12, d13, &
                                   r1, y1, r2, y2, r3, y3, &
                                   spline1, spline2, spline3, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(3), intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, &
                   r2, y2, r3, y3, spline1, spline2, spline3
   ! Output
   REAL(KIND=dp) :: integral

   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f1, f2, f3
   REAL(KIND=dp) :: norm
   INTEGER, DIMENSION(3) :: ileb
   INTEGER :: ngrid, i

   ileb(1) = get_number_of_lebedev_grid(n=nang(1))
   ileb(2) = get_number_of_lebedev_grid(n=nang(2))
   ileb(3) = get_number_of_lebedev_grid(n=nang(3))
   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2) + &
           lebedev_grid(ileb(3))%n * nshell(3)

   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(f1(ngrid))
   allocate(f2(ngrid))
   allocate(f3(ngrid))

   call build_threecenter_grid(ileb=ileb, nshell=nshell, d12=d12, d13=d13, &
                               addr2=.TRUE., grid_r=grid_r, grid_w=grid_w)

   do i=1,size(grid_w)
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))

      norm = sqrt(sum( (grid_r(i, :) + d12 )**2 ))
      call interpolation(r2, y2, spline2, norm, f2(i))

      norm = sqrt(sum( (grid_r(i, :) + d13 )**2 ))
      call interpolation(r3, y3, spline3, norm, f3(i))
   enddo
   integral = sum(grid_w * f1*f2*f3 )

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(f1)
   deallocate(f2)
   deallocate(f3)
end subroutine integration_threecenter

subroutine kinetic_energy(nang, nshell, r1, y1, r2, y2,&
                          spline1, spline2, integral)
   implicit none
   ! Input
   INTEGER, intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, &
                      r2, y2, spline1, spline2
   ! Output
   REAL(KIND=dp) :: integral

   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rf2, d2rf2, d2rf2_spline, f1, f2
   REAL(KIND=dp) :: norm
   INTEGER :: ngrid, ileb, i

   ileb = get_number_of_lebedev_grid(n=nang)
   ngrid = lebedev_grid(ileb)%n * nshell

   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(f1(ngrid))
   allocate(f2(ngrid))

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE.,&
                             grid_r=grid_r, grid_w=grid_w)
   ! < f1 | -0.5*🔺 | f2 >
   rf2 = -0.5_dp*r2*y2
   ! Get the 2nd derivative d_r^2(r*f2) as well as its spline
   call spline(r2, rf2, size(r2), 0.0_dp, 0.0_dp, d2rf2)
   call spline(r2, d2rf2, size(r2), 0.0_dp, 0.0_dp, d2rf2_spline)

   ! Divide by r
   d2rf2 = d2rf2/r2
   ! Laplace = 1/r * d_r^2(r*f2) ✓

   do i=1,size(grid_w)
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))
      call interpolation(r2, d2rf2, d2rf2_spline, norm, f2(i))
   enddo

   integral = sum(grid_w * f1*f2)

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(f1)
   deallocate(f2)
end subroutine kinetic_energy

subroutine coulomb_integral(nang, nshell, d12, r1, y1, r2, y2, s1, s2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, r2, y2, s1, s2
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w, grid_r2
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid, i, l
   ! Local variables (potential)
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: a, b, G, H, unique_r
   INTEGER :: nunique

   l = 0 ! Quantum number

   ileb(1) = get_number_of_lebedev_grid(n=nang(1))
   ileb(2) = get_number_of_lebedev_grid(n=nang(2))

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))

   call build_twocenter_grid(ileb=ileb, nshell=nshell, displacement=d12, &
                             addr2=.TRUE., grid_r=grid_r, grid_w=grid_w)

   ! First sort the grid in such a way that the smallest distances come first
   ! R is the distance form the origin to a grid point
   allocate(grid_r2(ngrid)) ! we sort by r^2

   do i=1, lebedev_grid(ileb(1))%n * nshell(1)
      grid_r2(i) = sqrt(sum( (grid_r(i, :))**2 ))
   enddo
   do i=lebedev_grid(ileb(1))%n * nshell(1)+1,ngrid 
      grid_r2(i) = sqrt(sum( (grid_r(i, :)+d12)**2 ))
   enddo

   call qsort(grid_r2)

   ! The Coulomb potential is radially symmetric (for now).
   ! As such we need to evaluate it only once for every radius

   allocate(unique_r(ngrid))
   unique_r = 0.0_dp
   unique_r(1) = grid_r2(1)
   nunique = 1
   do i=2,size(grid_r2)
      if (any( unique_r == grid_r2(i) )) cycle
      nunique = nunique + 1
      unique_r(nunique) = grid_r2(i)
   enddo

   ! Evaluate the Coulomb potential at every unique r
   allocate(a(nunique))
   allocate(b(nunique))
   allocate(G(nunique))
   allocate(H(nunique))
   a = unique_r**(-l-1)
   b = unique_r**l
   do i=1,nunique
      
   enddo

   ! do i=1,size(grid_w)
   !    norm = sqrt(sum( grid_r(i, :)**2 ))
   !    call interpolation(r1, y1, spline1, norm, f1(i))

   !    norm = sqrt(sum( (grid_r(i, :) + d12 )**2 ))
   !    call interpolation(r2, y2, spline2, norm, f2(i))
   ! enddo
   ! integral = sum(grid_w * f1*f2 )

   deallocate(a)
   deallocate(b)
   deallocate(G)
   deallocate(H)

   deallocate(unique_r)

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(grid_r2)

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

      if (bound1 .gt. 1e10) then
         yspline(1) = 0.0_dp
         u(1) = 0.0_dp
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
      qn = 0.0_dp
      un = 0.0_dp
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
      REAL(KIND=dp), intent(out) :: y

      INTEGER :: low, upper, mid
      REAL(KIND=dp) :: A, B, h
      ! find the closest grid point by bisection
      low = 1
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
               ((A**3.0_dp-A)*spline(low)+(B**3.0_dp-B)*spline(upper)) * (h**2.0_dp)/6.0_dp
      else if (gr(upper) .lt. r) then
         y = gy(upper)
         ! print *, 'Extrapolation!'
         ! print *, upper, gr(upper), r, y
      else if (gr(low) .gt. r) then
         y = gy(low)
         ! print *, 'Extrapolation!'
      endif
   end subroutine interpolation

recursive subroutine qsort(arr)!, brr)
   implicit none
   ! Input
   REAL(KIND=dp), DIMENSION(:) :: arr
   ! REAL(KIND=dp), DIMENSION(:, :) :: brr
   ! Local variables
   INTEGER :: first, last, mid, i, j
   ! REAL(KIND=dp), DIMENSION(3) :: temp_arr
   REAL(KIND=dp) :: a, b, temp

   first = 1
   last = size(arr, 1)
   a = arr( (first+last)/2 )
   i = first
   j = last

   do
      do while (arr(i) .lt. a)
         i = i+1
      enddo
      do while (arr(j) .gt. a)
         j = j-1
      enddo
      if (j .le. i) exit
      temp = arr(i); arr(i) = arr(j); arr(j) = temp;
      ! temp_arr = brr(i, :); brr(i, :) = brr(j, :); brr(j, :) = temp_arr;
      i = i+1
      j = j-1
   enddo

   if (first .lt. (i-1)) call qsort( arr(first:i-1))!, brr(first:i-1, :) )
   if ((j+1) .lt. last) call qsort( arr(j+1:last))!, brr(j+1:last, :) )
end subroutine qsort
end module eddi
