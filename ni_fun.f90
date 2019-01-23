module ni_fun
USE lebedev, ONLY: dp
implicit none
REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp

type :: type_fun
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y, y1, y2, y3, y4, y5
end type type_fun

public :: type_fun, fun_grid, prepare_fun, prepare_gauss
contains
subroutine prepare_fun(r, f, fun)
   implicit none
   REAL(KIND=dp), DIMENSION(:) :: r, f
   TYPE(type_fun), POINTER :: fun

   fun%r = r
   fun%y = f
   ! We want to produce the first through fifth derivatives
   call derivatives(r=r, y=f, y1=fun%y1, y2=fun%y2, y3=fun%y3)
   call spline(r, fun%y2, size(r), fun%y4)
   call spline(r, fun%y3, size(r), fun%y5)

   ! Now in principle this is not a unique assignment.
   ! The spline of the function is also the second derivative,
   ! The spline of the first derivative is also the third derivative
   ! The fourth and fifth derivative could be retrieved from another call
   ! to `derivatives`.
end subroutine prepare_fun

subroutine prepare_gauss(r, alpha, fun)
   implicit none
   REAL(KIND=dp), DIMENSION(:) :: r
   REAL(KIND=dp) :: alpha
   TYPE(type_fun), POINTER :: fun

   fun%y = exp(-alpha * r**2)
   fun%y1 = -2._dp*alpha*r* fun%y
   fun%y2 = 2._dp*(2._dp*alpha*r**2 - 1._dp) * alpha * fun%y
   fun%y3 = -4._dp*(2._dp*alpha*r**2 - 3._dp) * alpha**2 * r * fun%y
   fun%y4 = (16._dp*alpha**2*r**4 - 48._dp*alpha*r**2 + 12._dp)*alpha**2*fun%y
   fun%y5 = -(32._dp*alpha**2*r**4 - 160._dp*alpha*r**2 + 120._dp)*alpha**3*r*fun%y
end subroutine prepare_gauss

! Get the derivatives of a function by finite differences
subroutine derivatives(r, y, y1, y2, y3)
   implicit none
   REAL(KIND=dp), DIMENSION(:), intent(in) :: r, y
   REAL(KIND=dp), DIMENSION(size(r)), OPTIONAL, intent(out) :: y1, y2, y3
   ! Local variables
   REAL(KIND=dp), DIMENSION(3,3,5) :: c
   INTEGER :: ir
   if(present(y1)) y1 = 0._dp
   if(present(y2)) y2 = 0._dp
   if(present(y3)) y3 = 0._dp

   do ir=1, size(r)-5
      ! [...] where coeff[derivative, accuracy, coefficients]
      call forward_derivative_weights(order=3, x0=r(ir), r=r(ir:ir+6), coeff=c)
      if(present(y1)) y1(ir) = sum( c(1,2,1:3) * y(ir:ir+2) )
      if(present(y2)) y2(ir) = sum( c(2,2,1:4) * y(ir:ir+3) )
      if(present(y3)) y3(ir) = sum( c(3,2,1:5) * y(ir:ir+4) )
   enddo
end subroutine derivatives

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
   REAL(KIND=dp) :: sig, p, un, qn, der1, dern

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

! Generate an equally spaced grid up to max
subroutine fun_grid(r, max)
   implicit none
   REAL(KIND=dp), DIMENSION(:), intent(out) :: r
   REAL(KIND=dp), intent(in) :: max
   INTEGER :: i
   do i=1,size(r)
      r(i) = REAL(i-1, dp)*max/size(r)
   enddo
end subroutine fun_grid

subroutine allocate_fun(fun, n)
   implicit none
   INTEGER :: n
   TYPE(type_fun), POINTER :: fun
   
   if (.not. allocated(fun%r)) allocate(fun%r(n))
   if (.not. allocated(fun%y)) allocate(fun%y(n))
   if (.not. allocated(fun%y1)) allocate(fun%y1(n))
   if (.not. allocated(fun%y2)) allocate(fun%y2(n))
   if (.not. allocated(fun%y3)) allocate(fun%y3(n))
   if (.not. allocated(fun%y4)) allocate(fun%y4(n))
   if (.not. allocated(fun%y5)) allocate(fun%y5(n))
end subroutine allocate_fun

subroutine deallocate_fun(fun)
   implicit none
   INTEGER :: n
   TYPE(type_fun), POINTER :: fun
   
   deallocate(fun%r)
   deallocate(fun%y)
   deallocate(fun%y1)
   deallocate(fun%y2)
   deallocate(fun%y3)
   deallocate(fun%y4)
   deallocate(fun%y5)
end subroutine deallocate_fun

end module
