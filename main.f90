program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, &
                radial_integration
USE grid, ONLY: grid_parameters
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic
implicit none

TYPE(type_atom), DIMENSION(3) :: atoms
INTEGER, DIMENSION(:), ALLOCATABLE :: nang, nshell
REAL(KIND=dp), DIMENSION(3) :: d12, d13
REAL(KIND=dp) :: integral

CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r1, r2, y1, y2, r3, y3,&
                                            spline1, spline2, spline3

! Temp variables
REAL(KIND=dp), DIMENSION(500) :: d12_range
INTEGER :: i


! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=50 , loud=.FALSE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! return

! Build parameters
atoms(1)%r = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
atoms(1)%z = 1
atoms(2)%r = (/ 1.0_dp, 0.0_dp, 0.0_dp /)
atoms(2)%z = 1
atoms(3)%r = (/ 0.0_dp, 2.0_dp, 0.0_dp /)
atoms(3)%z = 1

d12 = atoms(2)%r - atoms(1)%r
d13 = atoms(3)%r - atoms(1)%r

call read_nfun(fn1, r1, y1)
call spline(r1, y1, size(r1), 0.0_dp, 0.0_dp, spline1)
call read_nfun(fn2, r2, y2)
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
nshell = (/ 1000, 1000 /)

d12_range = [ (0.1_dp*REAL(i, dp) , i=1,500) ]
open(unit=100, file='coulomb_integral')
do i=1,500
   d12 = (/ d12_range(i), 0._dp, 0._dp /)
   call coulomb_integral(nang=nang, nshell=nshell, d12=d12, r1=r1, y1=y1, r2=r2, &
                         y2=y2, s1=spline1, s2=spline2, integral=integral)
   write(100, *) d12_range(i), integral
enddo
close(100)

deallocate(nang)
deallocate(nshell)

end program hallo