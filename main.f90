program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, coulomb_integral_grid, &
                radial_integration
USE grid, ONLY: grid_parameters, radial_grid
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic, &
                    test_coulomb, test_radial_grid
implicit none

TYPE(type_atom), DIMENSION(3) :: atoms
INTEGER, DIMENSION(:), ALLOCATABLE :: nang, nshell
REAL(KIND=dp), DIMENSION(3) :: d12, d13
REAL(KIND=dp) :: integral

CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r1, r2, y1, y2, r3, y3,&
                                            spline1, spline2, spline3

! Temp variables (Coulomb)
REAL(KIND=dp), DIMENSION(500) :: d12_range
REAL(KIND=dp) :: ri, norm
INTEGER :: i, ngrid, coul_n


! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=100 , loud=.FALSE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! call test_coulomb(ntests=100, loud=.FALSE.)
call test_radial_grid(ntests=200)
return

! Build parameters
atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
atoms(1)%z = 1
atoms(2)%r = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
atoms(2)%z = 1
atoms(3)%r = (/ 0.0_dp, 2.0_dp, 0.0_dp /)
atoms(3)%z = 1

d12 = atoms(2)%r - atoms(1)%r
d13 = atoms(3)%r - atoms(1)%r

ngrid = 50000
allocate(r1(ngrid))
allocate(y1(ngrid))
allocate(spline1(ngrid))
allocate(r2(ngrid))
allocate(y2(ngrid))
allocate(spline2(ngrid))

call radial_grid(r=r1, wr=spline1, n=ngrid, addr2=.FALSE.)
r1 = r1(ngrid:1:-1)
r2 = r1

! call read_nfun(fn1, r1, y1)
y1 = exp(-r1**2.0_dp)
call spline(r1, y1, size(r1), 0.0_dp, 0.0_dp, spline1)

! call read_nfun(fn2, r2, y2)
y2 = exp(-r2**2.0_dp)
call spline(r2, y2, size(r2), 0.0_dp, 0.0_dp, spline2)

call read_nfun(fn2, r3, y3)
call spline(r3, y3, size(r3), 0.0_dp, 0.0_dp, spline3)

allocate(nang(2))
allocate(nshell(2))

call grid_parameters(atoms(1)%z, nang(1), nshell(1))
call grid_parameters(atoms(2)%z, nang(2), nshell(2))

! Coulomb integral
integral = 1e6
print *, REPEAT('-', 30) // '! Coulomb integral !' // REPEAT('-', 30)

nang = (/ 590, 590 /)
nshell = (/ 200, 200 /)
coul_n = 100

d12_range = [ (0.25_dp*REAL(i, dp) , i=0,499) ]
open(unit=100, file='coulomb_integral')
do i=1,50
   d12 = (/ d12_range(i), 0._dp, 0._dp /)
   norm = sqrt(sum(d12**2))
   call coulomb_integral(nang=nang, nshell=nshell, coul_n=coul_n, d12=d12, r1=r1, y1=y1,&
                         r2=r2, y2=y2, s1=spline1, s2=spline2, integral=integral)
   if (norm .eq. 0.0_dp) then
   ri = pi**3 / sqrt(1.0_dp * 1.0_dp)**3 * 2.0_dp*sqrt(1.0_dp*1.0_dp/(1.0_dp+1.0_dp))/sqrt(pi) ! lim x->0 : erf(a*x)/x 
   else
      ri = pi**3 / sqrt(1.0_dp * 1.0_dp)**3 * erf( sqrt(1.0_dp*1.0_dp/(1.0_dp+1.0_dp)) * norm )/norm 
   endif
   print *, d12_range(i), integral, ri
   ! call coulomb_integral_grid(nang=nang, nshell=nshell, d12=d12, r1=r1, y1=y1, r2=r2, &
   !                            y2=y2, s1=spline1, s2=spline2, integral=integral)
   ! print *, d12_range(i), integral
   ! return
   write(100, *) d12_range(i), integral
enddo
close(100)

deallocate(nang)
deallocate(nshell)

deallocate(r1)
deallocate(y1)
deallocate(spline1)
deallocate(r2)
deallocate(y2)
deallocate(spline2)

end program hallo