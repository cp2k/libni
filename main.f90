program hallo 
USE lebedev, ONLY: dp, get_number_of_lebedev_grid
USE ni_fun, ONLY: allocate_fun, prepare_fun, deallocate_fun, type_fun, fun_grid
USE ni_grid, ONLY: deallocate_grid, build_onecenter_grid, type_grid
implicit none
REAL(KIND=dp), DIMENSION(50) :: r, f
TYPE(type_fun), POINTER :: pfun
TYPE(type_fun), TARGET :: fun

TYPE(type_grid), POINTER :: pgrid
TYPE(type_grid), target :: grid
INTEGER :: i

call fun_grid(r=r, max=10._dp)
f = exp(-r**2)

pfun => fun

call allocate_fun(fun=pfun, n=size(r))
call prepare_fun(r=r, f=f, fun=pfun)
do i=1,size(r)
  print *, fun%r(i), fun%y(i), fun%y1(i), fun%y2(i)
enddo
call deallocate_fun(fun=pfun)

pgrid => grid
call build_onecenter_grid(ileb=4, nshell=50, addr2=.TRUE.,&
                          quadr=1, grid=pgrid)
print *, 'grid'
print *, size(grid%r)
print *, size(grid%w)
print *, size(grid%dw)
print *, get_number_of_lebedev_grid(n=10)
call deallocate_grid(grid=pgrid)
end program hallo