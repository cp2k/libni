module eddi
USE lebedev, ONLY: lebedev_grid,&
                   get_number_of_lebedev_grid,&
                   dp
USE grid, ONLY: build_onecenter_grid, build_twocenter_grid, build_threecenter_grid, &
                type_grid_point, radial_grid
USE spherical_harmonics, ONLY: rry_lm
implicit none
REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

type :: type_atom
   REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
   INTEGER :: z = 1
end type type_atom

public :: integration_twocenter, integration_onecenter, integration_threecenter,&
           radial_integration, qsort

contains
! Get the derivatives of a function by finite differences
subroutine derivatives(r, y, y1, y2, y3)
   implicit none
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, y
   REAL(KIND=dp), DIMENSION(size(r)), intent(out) :: y1, y2, y3
   ! Local variables
   REAL(KIND=dp), DIMENSION(3,3,5) :: c
   INTEGER :: ir

   y1 = 0._dp; y2 = 0._dp; y3 = 0._dp
   do ir=1, size(r)-5
      ! [...] where coeff[derivative, accuracy, coefficients]
      call forward_derivative_weights(order=3, x0=r(ir), r=r(ir:ir+6), coeff=c)
      y1(ir) = sum( c(1,2,1:3) * y(ir:ir+2) )
      y2(ir) = sum( c(2,2,1:4) * y(ir:ir+3) )
      y3(ir) = sum( c(3,2,1:5) * y(ir:ir+4) )
   enddo
end subroutine derivatives

! **********************************************
!> \brief Computes the radial integral of f(r)
!> \param f(n): The tabulated function at n grid points
!> \param r(n): The tabulated grid points
!> \param n: The number of radial grid points
!> \param integral: The integral's value
!> \author 
! **********************************************
subroutine radial_integration(f, r, n, addr2, integral)
   implicit none
   REAL(KIND=dp), DIMENSION(:), intent(in) :: f, r
   INTEGER, intent(in) :: n
   LOGICAL, intent(in) :: addr2
   REAL(KIND=dp) :: integral
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rad, wr, d2f, fun
   INTEGER :: i

   allocate(rad(n))
   allocate(wr(n))
   allocate(d2f(n))
   allocate(fun(n))

   integral = 0.0_dp
   
   ! Put the radial grid points into `rad` and their weights into `wr`
   call radial_grid(r=rad, wr=wr, n=n, addr2=addr2, quadr=1)

   ! Create the spline
   call spline(r=r, y=f, n=size(r), yspline=d2f)
   
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

! Compute <Y_L | f>_w^ri for all r_i
subroutine pp_projector(l, m, r, f, s, d12, p)
   implicit none
   ! Inputs
   INTEGER, intent(in) :: l, m
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, f, s
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Outputs
   REAL(KIND=dp), DIMENSION(size(r)) :: p
   ! Local variables
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: ylm, funs
   REAL(KIND=dp), DIMENSION(3) :: gr
   REAL(KIND=dp) :: norm
   INTEGER :: ileb, iang, irad

   ileb = get_number_of_lebedev_grid(l=l+5)
   allocate(ylm(lebedev_grid(ileb)%n))
   do iang=1,lebedev_grid(ileb)%n
      call rry_lm(l=l, m=m, r=lebedev_grid(ileb)%r(:, iang), y=ylm(iang))
   enddo

   ! Compute the value of f at each point
   allocate(funs(lebedev_grid(ileb)%n))
   do irad=1,size(r)
      funs = 0.0_dp
      do iang=1,lebedev_grid(ileb)%n
         gr = r(irad) * lebedev_grid(ileb)%r(:, iang)
         norm = sqrt(sum( (gr-d12)**2 ))
         call interpolation(gr=r, gy=f, spline=s, r=norm, y=funs(iang))
      enddo
      p(irad) = sum(lebedev_grid(ileb)%w * ylm * funs)
   enddo
   deallocate(funs)
   deallocate(ylm)
   p = 4.0_dp * pi * p
end subroutine pp_projector

!  nloc                            L           L
! V     = Σ 1/Ω * Σ w (V_l(r_i) * P (alpha) * P (beta))
!         L       i                i           i
subroutine pp_nonloc(rv, v, rp1, p1, rp2, p2, d12, d13, lmax, nrad, integral)
   implicit none
   ! Input
   REAL(KIND=dp), DIMENSION(:), intent(in) :: rv, v, rp1, p1, rp2, p2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
   INTEGER, intent(in) :: lmax, nrad
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(nrad) :: r, wr, gv, gp1, gp2, gsv, gsp1, gsp2,&
                                     v_pp, proj1, proj2
   REAL(KIND=dp), DIMENSION(size(rv)) :: sv
   REAL(KIND=dp), DIMENSION(size(rp1)) :: sp1
   REAL(KIND=dp), DIMENSION(size(rp2)) :: sp2
   REAL(KIND=dp), DIMENSION((lmax+1)**2) :: integral_sub
   INTEGER :: i, il, im, h!elp
   ! End header

   ! First we transpose the three functions to a common grid
   call radial_grid(r=r, wr=wr, n=nrad, addr2=.TRUE., quadr=2) !2=hermite
   call spline(r=rv, y=v, n=size(rv), yspline=sv)
   call spline(r=rp1, y=p1, n=size(rp1), yspline=sp1)
   call spline(r=rp2, y=p2, n=size(rp2), yspline=sp2)

   do i=1,nrad
      call interpolation(rv , v , sv , r=r(i), y=gv(i) )
      call interpolation(rp1, p1, sp1, r=r(i), y=gp1(i))
      call interpolation(rp2, p2, sp2, r=r(i), y=gp2(i))
   enddo

   ! The interpolated functions need splines as well
   call spline(r=r, y=gv , n=size(r), yspline=gsv)
   call spline(r=r, y=gp1, n=size(r), yspline=gsp1)
   call spline(r=r, y=gp2, n=size(r), yspline=gsp2)

   ! Then we go over all L={l,m} where (l .le. lmax)
   h = 0
   do il=0,lmax
      do im=-il,+il
         h = h + 1
         ! gv – the interpolation of V_l – we can use as is
         call pp_projector(l=il, m=im, r=r, f=gp1, s=gsp1, d12=d12, p=proj1)
         call pp_projector(l=il, m=im, r=r, f=gp2, s=gsp2, d12=d13, p=proj2)
         integral_sub(h) = sum(wr * gv * proj1 * proj2)
      enddo
   enddo
   integral = sum(integral_sub)
end subroutine pp_nonloc

subroutine integration_onecenter(nang, nshell, r, y, spline, quadr, integral)
   implicit none
   ! Input
   INTEGER, intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r, y, spline
   INTEGER, intent(in) :: quadr
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

   call build_onecenter_grid(ileb=ileb, nshell=nshell, addr2=.TRUE., quadr=quadr,&
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

subroutine integration_twocenter(l, m, nshell, d12, r1, y1, r2, y2, &
                                 spline1, spline2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, r2, y2, &
                                              spline1, spline2
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: f1, f2
   REAL(KIND=dp) :: norm, ylm
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid, i

   ileb(1) = get_number_of_lebedev_grid(n=302)
   ileb(2) = get_number_of_lebedev_grid(n=302)

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(f1(ngrid))
   allocate(f2(ngrid))

   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.TRUE., grid_r=grid_r, grid_w=grid_w)

   do i=1,ngrid
      if (grid_w(i) .eq. 0.0_dp) cycle         
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))
      call rry_lm(l=l(1), m=m(1), r=grid_r(i, :)/norm, y=ylm)
      f1(i) = f1(i) * ylm

      norm = sqrt(sum( (grid_r(i, :) - d12 )**2 ))
      call interpolation(r2, y2, spline2, norm, f2(i))
      call rry_lm(l=l(2), m=m(2), r=(grid_r(i, :)-d12)/norm, y=ylm)
      f2(i) = f2(i) * ylm
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

   do i=1,ngrid
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))

      norm = sqrt(sum( (grid_r(i, :) - d12 )**2 ))
      call interpolation(r2, y2, spline2, norm, f2(i))

      norm = sqrt(sum( (grid_r(i, :) - d13 )**2 ))
      call interpolation(r3, y3, spline3, norm, f3(i))
   enddo
   integral = sum(grid_w * f1*f2*f3 )

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(f1)
   deallocate(f2)
   deallocate(f3)
end subroutine integration_threecenter

subroutine kinetic_energy(l, m, nshell, r1, y1, r2, y2, d12,&
                          spline1, spline2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: l, m, nshell
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r1, y1, &
                                              r2, y2, spline1, spline2
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: rf2, d2rf2, d2rf2_spline, f1, f2
   REAL(KIND=dp) :: norm, ylm, df2
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: ngrid, i

   ileb = get_number_of_lebedev_grid(n=590)  ! TODO
   if (all(d12 .eq. (/0._dp, 0._dp, 0._dp/))) then
      ngrid = lebedev_grid(ileb(1))%n * nshell(1)
   else
      ngrid = lebedev_grid(ileb(1))%n * nshell(1) &
               + lebedev_grid(ileb(2))%n * nshell(2)
   endif

   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))
   allocate(f1(ngrid))
   allocate(f2(ngrid))
   allocate(d2rf2(size(r2)))
   allocate(d2rf2_spline(size(r2)))

   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, addr2=.FALSE.,&
                             grid_r=grid_r, grid_w=grid_w)
   ! < f1 | -0.5*🔺 | f2 >
   ! Laplace_r = 1/r * d_r^2(r*f2)
   rf2 = r2*y2
   ! Get the 2nd derivative d_r^2(r*f2) as well as its spline
   call spline(r2, rf2, size(r2), d2rf2)
   call spline(r2, d2rf2, size(r2), d2rf2_spline)
   d2rf2 = d2rf2/r2

   do i=1,ngrid
      ! T = -0.5 * ∫f1*Ylm * (D_r f2(r) - l'(l'+1)*f2/r^2 ) * Yl'm'
      norm = sqrt(sum( grid_r(i, :)**2 ))
      call interpolation(r1, y1, spline1, norm, f1(i))
      call rry_lm(l=l(1), m=m(1), r=grid_r(i, :)/norm, y=ylm)
      f1(i) = f1(i) * ylm * norm**2

      norm = sqrt(sum( (grid_r(i, :) - d12)**2 ))
      call interpolation(r2, d2rf2, d2rf2_spline, norm, df2)  ! D_r f2
      call interpolation(r2, y2, spline2, norm, f2(i))        ! f2
      call rry_lm(l=l(2), m=m(2), r=(grid_r(i, :) - d12)/norm, y=ylm)

      ! (D_r f2(r) - l'(l'+1)*f2/r^2 ) * Yl'm'
      f2(i) = df2 - REAL(l(2)*(l(2)+1), dp)*f2(i)/norm**2
      f2(i) = f2(i) * ylm
   enddo

   integral = -0.5_dp*sum(grid_w * f1*f2)

   deallocate(grid_r)
   deallocate(grid_w)
   deallocate(f1)
   deallocate(f2)
   deallocate(d2rf2)
   deallocate(d2rf2_spline)
end subroutine kinetic_energy

subroutine coulomb_integral(nang, nshell, coul_n, d12, r1, y1, r2, y2, s1, s2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: nang, nshell
   INTEGER, intent(in) :: coul_n
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, r2, y2, s1, s2
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   INTEGER :: i, j, l
   ! Local variables (potential)
   REAL(KIND=dp), DIMENSION(coul_n) :: f, gi, hi, G, H, coul_w
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: coul_r, pot, pots

   l = 0 ! Quantum number TODO

   ! 1: Evaluate the Coulomb potential on a radial grid around A
   ! ! the integral is purely radial => addr2=False
   allocate(coul_r(coul_n))
   allocate(pot(coul_n))
   allocate(pots(coul_n))
   call radial_grid(r=coul_r, &
                    wr=coul_w, &
                    n=coul_n, addr2=.FALSE., quadr=1)

   do i=1,coul_n
      call interpolation(r1, y1, s1, coul_r(i), f(i))
   enddo
   gi = coul_w * coul_r**(l+2) * f
   hi = coul_w * coul_r**(1-l) * f

   G(coul_n) = sum(gi)
   H(coul_n) = 0.0_dp
   do j=coul_n,1,-1
      pot(j) = coul_r(j)**(-l-1) * G(j) + coul_r(j)**l * H(j)
      G(j-1) = G(j) - gi(j)
      H(j-1) = H(j) + hi(j)
   enddo

   pot = pot * 4.0_dp*pi/(2.0_dp*l+1.0_dp)
   call spline(coul_r, pot, coul_n, pots)

   ! 2: Calculate the overlap of y2(r-d12) and the coulomb potential
   call integration_twocenter(l=(/l,l/), m=(/0,0/), nshell=nshell, d12=d12, &
                              r1=coul_r, y1=pot, r2=r2, y2=y2,&
                              spline1=pots, spline2=s2, integral=integral)

   deallocate(coul_r)
   deallocate(pot)
   deallocate(pots)
end subroutine coulomb_integral

subroutine coulomb_integral_grid(nang, nshell, d12, r1, y1, r2, y2, s1, s2, integral)
   implicit none
   ! Input
   INTEGER, DIMENSION(2), intent(in) :: nang, nshell
   REAL(KIND=dp), DIMENSION(3), intent(in) :: d12
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: r1, y1, r2, y2, s1, s2
   ! Output
   REAL(KIND=dp) :: integral
   ! Local variables
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: grid_r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: grid_w, f, gi, hi, G, H
   REAL(KIND=dp) :: norm
   INTEGER, DIMENSION(2) :: ileb
   INTEGER :: i, j, idx, l, ngrid
   ! Local variables (potential)
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: coul_r, coul_w, grid_r2, pot, f1, f2
   REAL(KIND=dp) :: temp
   INTEGER :: coul_n 

   l = 0 ! Quantum number

   ! 1: Evaluate the Coulomb potential on the two-center grid
   ileb(1) = get_number_of_lebedev_grid(n=nang(1))
   ileb(2) = get_number_of_lebedev_grid(n=nang(2))

   ngrid = lebedev_grid(ileb(1))%n * nshell(1) + &
           lebedev_grid(ileb(2))%n * nshell(2)
   allocate(grid_r(ngrid, 3))
   allocate(grid_w(ngrid))

   ! TODO
   ! When integrating f1 - retrieving the coulomb potential - the radial weights
   ! need to be included in order to get the right potential. After that the
   ! regular weights are used for the overlap integral
   call build_twocenter_grid(ileb=ileb, nshell=nshell, d12=d12, &
                             addr2=.FALSE., grid_r=grid_r, grid_w=grid_w)

   ! First we need all unique distances
   !! together with the weights
   allocate(coul_r(ngrid))
   allocate(coul_w(ngrid))
   allocate(grid_r2(ngrid))
   do i=1,ngrid
      grid_r2(i) = sqrt(sum(grid_r(i, :)**2))
   enddo
   call qsort_sim2(grid_r2, grid_w)

   coul_n = 0
   do i=1,ngrid
      if (grid_r2(i) == temp) cycle
      coul_n = coul_n + 1
      temp = grid_r2(i)
      coul_r(coul_n) = temp
      coul_w(coul_n) = grid_w(i)
   enddo

   ! Next we evaluate the coulomb potential and spline at those distances
   allocate(pot(coul_n))
   allocate(f(coul_n))
   allocate(gi(coul_n))
   allocate(hi(coul_n))
   allocate(G(coul_n))
   allocate(H(coul_n))

   do i=1,coul_n
      call interpolation(r1, y1, s1, coul_r(i), f(i))
   enddo
   gi = coul_r**(l+2) * f * coul_w
   hi = coul_r**(1-l) * f * coul_w

   G(coul_n) = sum(gi)
   H(coul_n) = 0.0_dp
   do j=coul_n,1,-1
      pot(j) = coul_r(j)**(-l-1) * G(j) + coul_r(j)**l * H(j)
      G(j-1) = G(j) - gi(j)
      H(j-1) = H(j) + hi(j)
   enddo

   pot = pot * 4.0_dp*pi/(2.0_dp*l+1.0_dp)
   deallocate(f)
   deallocate(gi)
   deallocate(hi)
   deallocate(G)
   deallocate(H)

   ! 2: Finally calculate the overlap of y2(r-d12) and the coulomb potential
   allocate(f1(ngrid))
   allocate(f2(ngrid))
   do i=1,ngrid
      ! evaluate the potential
      norm = sqrt(sum( grid_r(i, :)**2 ))
      ! Look for 
      do j=1,coul_n
         f1(i) = j
         if ( norm .eq. coul_r(j) ) then
            f1(i) = pot(j)
            exit
         endif
      enddo
      if (f1(i) .eq. j) then
         print *, 'oh noes', i, j
      endif
      ! evaluate y2 at the same point
      norm = sqrt(sum( (grid_r(i, :) - d12 )**2 ))
      call interpolation(r2, y2, s2, norm, f2(i))
   enddo

   integral = sum(grid_w * f1*f2 )

   deallocate(f1)
   deallocate(f2)
   deallocate(pot)
   deallocate(grid_r2)
   deallocate(coul_r)
   deallocate(coul_w)
   deallocate(grid_r)
   deallocate(grid_w)
end subroutine coulomb_integral_grid

subroutine forward_derivative_weights(order, x0, r, coeff)
   implicit none
   ! Input
   INTEGER, intent(in) :: order
   REAL(KIND=dp), intent(in) :: x0
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r
   ! Output
   REAL(KIND=dp), DIMENSION(3,3,5) :: coeff
   ! Local variables
   INTEGER :: points, n, nu, m
   REAL(KIND=dp) :: c1, c2, c3
   REAL(KIND=dp), DIMENSION(0:order,0:size(r),0:size(r)) :: d

   points = size(r)-1

   d = 0.0_dp
   d(0,0,0) = 1.0_dp
   c1 = 1.0_dp

   do n=1,points
      c2 = 1.0_dp
      do nu=0,n-1
         c3 = r(n+1) - r(nu+1)
         c2 = c2 * c3
         do m=0,min(n,order)
            d(m,n,nu) = (r(n+1)-x0)*d(m,n-1,nu)
            if (m .ne. 0) then
               d(m,n,nu) = d(m,n,nu) - m*d(m-1,n-1,nu)
            endif
            d(m,n,nu) = d(m,n,nu)/c3
         enddo
      enddo
      do m=0,min(n,order)
         if (m .ne. 0) then
            d(m,n,n) = m*d(m-1,n-1,n-1)
         endif
         d(m,n,n) = c1/c2*(d(m,n,n) - (r(n)-x0)*d(m,n-1,n-1))
      enddo
      c1 = c2
   enddo

   ! d contains way more information than we need it to.
   ! instead we construct a smaller field `coeff`
   ! where coeff[derivative, accuracy, coefficients]
   ! First derivative
   coeff(1,1,:) = d(1,1,0:4)
   coeff(1,2,:) = d(1,2,0:4)
   coeff(1,3,:) = d(1,3,0:4)
   ! Second derivative
   coeff(2,1,:) = d(2,2,0:4)
   coeff(2,2,:) = d(2,3,0:4)
   coeff(2,3,:) = d(2,4,0:4)
   ! Third derivative
   coeff(3,1,:) = d(3,3,0:4)
   coeff(3,2,:) = d(3,4,0:4)
   ! coeff(3,3,:) = d(3,5,0:4) ! this one has 6 coefficients
end subroutine forward_derivative_weights

subroutine derivative_point(r, y, r0, y1)
   implicit none
   ! Input
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, y
   REAL(KIND=dp), intent(in) :: r0
   ! Output
   REAL(KIND=dp) :: y1
   ! Local variables
   REAL(KIND=dp), DIMENSION(3,3,5) :: coeff
   INTEGER :: low, upper

   call bisection(r=r, r0=r0, low=low, upper=upper)
   call forward_derivative_weights(order=2, x0=r0, r=r(low:size(r)), coeff=coeff)
   y1 = sum( coeff(1,2,1:3) * y(low:low+2) )  ! 2nd order accuracy
   ! this will crash/yield 0 if we want the derivative at the end points
end subroutine derivative_point

subroutine spline(r, y, n, yspline)
   implicit none
   ! Input
   INTEGER, INTENT(in) :: n
   REAL(KIND=dp), DIMENSION(:), INTENT(in) :: r, y
   ! Output
   REAL(KIND=dp), DIMENSION(n) :: yspline
   ! Local variables
   REAL(KIND=dp), DIMENSION(n) :: u
   REAL(KIND=dp), DIMENSION(3,3,5) :: coeff
   INTEGER :: i
   REAL(KIND=dp) :: sig, p, un, qn, der1, dern, h

   ! der1 is the first derivative at r(1)
   yspline(1) = -0.5_dp
   call forward_derivative_weights(order=2, x0=r(1), r=r, coeff=coeff)
   der1 = sum( coeff(1,2,1:4) * y(1:4) )  ! 2nd order accuracy
   u(1) = (3.0_dp/(r(2)-r(1))) * ((y(2)-y(1))/(r(2)-r(1))-der1)

   do i=2,n-1
      sig = (r(i)-r(i-1))/(r(i+1)-r(i-1))
      p = sig*yspline(i-1)+2.0_dp
      yspline(i) = (sig-1.0_dp)/p

      u(i) = (6.0_dp * ( (y(i+1)-y(i))/(r(i+1)-r(i)) -&
                         (y(i)-y(i-1))/(r(i)-r(i-1)) )/&
               (r(i+1)-r(i-1)) - sig*u(i-1)) / p
   enddo

   ! zero first derivative at r->infinity seems reasonable for our purposes
   qn = 0.5_dp
   dern = 0.0_dp
   un = (3.0_dp/(r(n)-r(n-1))) * (dern - (y(n)-y(n-1))/(r(n)-r(n-1)) )

   yspline(n) = (un-qn*u(n-1))/(qn*yspline(n-1)+1.0_dp);

   do i=n-1,1,-1
      yspline(i) = yspline(i)*yspline(i+1)+u(i)
   enddo
end subroutine spline

! Given a function `gy` on a grid `gr` and a requested
! function value y(r) interpolates the function value `y` using `spline`
subroutine interpolation(gr, gy, spline, r, y, yprime)
   ! Input
   REAL(KIND=dp), DIMENSION(:), intent(in) :: gr, gy
   REAL(KIND=dp), DIMENSION(size(gr)), intent(in) :: spline
   REAL(KIND=dp), intent(in) :: r
   ! Output
   REAL(KIND=dp) :: y
   REAL(KIND=dp), OPTIONAL :: yprime
   ! Local variables
   INTEGER :: low, upper
   REAL(KIND=dp) :: A, B, C, D, h

   ! find the closest grid point by bisection
   call bisection(r=gr, r0=r, low=low, upper=upper)

   if (gy(low).eq.0._dp .and. gy(upper).eq.0._dp) then
      y = 0._dp
   elseif (gr(upper) .eq. r) then
      y = gy(upper)
      if (present(yprime)) call derivative_point(r=gr, y=gy, r0=r, y1=yprime)
   else if (gr(low) .eq. r) then
      y = gy(low)
      if (present(yprime)) call derivative_point(r=gr, y=gy, r0=r, y1=yprime)
   else if ((gr(upper) .gt. r) .and. (gr(low) .lt. r)) then
      h = gr(upper)-gr(low)
      A = (gr(upper)-r)/h
      B = (r-gr(low))/h
      C = (A**3.0_dp-A) * (h**2.0_dp)/6.0_dp
      D = (B**3.0_dp-B) * (h**2.0_dp)/6.0_dp
      y = A*gy(low) + B*gy(upper) + C*spline(low) + D*spline(upper)
      if (present(yprime)) then
      yprime = ( gy(upper)-gy(low) )/h - (3._dp*A**2-1._dp)*h/6._dp*spline(low)&
                                       + (3._dp*B**2-1._dp)*h/6._dp*spline(upper)
      endif
   else if (gr(upper) .lt. r) then
      ! If the supplied r is higher than maxval(gr)
      y = gy(upper)
      if (present(yprime)) yprime = 0.0_dp
      ! print *, 'Extrapolation up!'
   else if (gr(low) .gt. r) then
      ! If the supplied r is lower than minval(gr)
      y = gy(low)
      if (present(yprime)) call derivative_point(r=gr, y=gy, r0=r, y1=yprime)
      ! print *, 'Extrapolation!'
   endif
end subroutine interpolation

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

subroutine bisection(r, r0, low, upper)
   implicit none
   ! Input
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r
   REAL(KIND=dp), intent(in) :: r0
   ! Output
   INTEGER :: low, upper
   ! Local variables
   INTEGER :: mid

   low = 1
   upper = size(r)
   do while (upper .gt. low+1)
      mid = NINT((low+upper)/2.0_dp)
      if (r(mid) .gt. r0) then
         upper = mid
      else
         low = mid
      endif
      if (r(low) .eq. r0) upper = low 
      if (r(upper) .eq. r0) low = upper 
   enddo
end subroutine bisection

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

recursive subroutine qsort_sim2(arr, brr)
   implicit none
   ! Input
   REAL(KIND=dp), DIMENSION(:) :: arr
   REAL(KIND=dp), DIMENSION(:) :: brr
   ! Local variables
   INTEGER :: first, last, mid, i, j
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
      temp = brr(i); brr(i) = brr(j); brr(j) = temp;
      i = i+1
      j = j-1
   enddo

   if (first .lt. (i-1)) call qsort_sim2( arr(first:i-1), brr(first:i-1) )
   if ((j+1) .lt. last) call qsort_sim2( arr(j+1:last), brr(j+1:last) )
end subroutine qsort_sim2

end module eddi
