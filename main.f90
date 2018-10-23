program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, coulomb_integral_grid, &
                radial_integration
USE grid, ONLY: grid_parameters, radial_grid, gauher
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic, &
                    test_coulomb, test_radial_grid, test_radial_quadrature
implicit none

INTEGER :: nang, nshell
REAL(KIND=dp) :: integral

CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r1, r2, y1, y2, r3, y3,&
                                            spline1, spline2, spline3

! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=1000 , loud=.TRUE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! call test_coulomb(ntests=100, loud=.FALSE.)
! call test_radial_grid(ntests=200)
! call test_radial_quadrature(ntests=100, loud=.FALSE.)
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
return

ngrid = 50000
allocate(r1(ngrid))
allocate(y1(ngrid))
allocate(spline1(ngrid))
allocate(r2(ngrid))
allocate(y2(ngrid))
allocate(spline2(ngrid))

call radial_grid(r=r1, wr=spline1, n=ngrid, addr2=.FALSE., quadr=1)
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


call grid_parameters(atoms(1)%z, nang(1), nshell(1))
call grid_parameters(atoms(2)%z, nang(2), nshell(2))



deallocate(r1)
deallocate(y1)
deallocate(spline1)
deallocate(r2)
deallocate(y2)
deallocate(spline2)

end program hallo