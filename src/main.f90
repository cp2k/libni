program main_test 
use lebedev, only: dp, lebedev_grid, get_number_of_lebedev_grid
use ni_fun, only: allocate_fun, prepare_fun, deallocate_fun, type_fun, fun_grid
use ni_grid, only: deallocate_grid, build_onecenter_grid, type_grid
implicit none
real(kind=dp), dimension(50) :: r, f
type(type_fun), pointer :: pfun
type(type_fun), TARGET :: fun

type(type_grid), pointer :: pgrid
type(type_grid), target :: grid
integer :: i, j

i = get_number_of_lebedev_grid(l=3)

call fun_grid(r=r, max=10._dp)
f = exp(-r**2)

pfun => fun

call allocate_fun(fun=pfun, n=size(r))
call prepare_fun(r=r, f=f, fun=pfun)
do i=1,size(r)
  ! print *, fun%r(i), fun%y(i), fun%y1(i), fun%y2(i)
enddo
call deallocate_fun(fun=pfun)

pgrid => grid
call build_onecenter_grid(ileb=4, nshell=50, addr2=.true.,&
                          quadr=1, grid=pgrid)
print *, 'grid'
print *, size(grid%r)
print *, size(grid%w)
print *, size(grid%dw)
call deallocate_grid(grid=pgrid)
end program main_test
