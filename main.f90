program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, coulomb_integral_grid, &
                radial_integration, pp_projector, pp_nonloc
USE grid, ONLY: grid_parameters, radial_grid, gauher
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic, &
                    test_coulomb, test_radial_grid, test_radial_quadrature
implicit none

CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'

! Local variables (Projector)
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r1, wr1, y1, s1
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: projector1, projector2
INTEGER :: i, nang, nshell, ngrid

! Local variables (TEST)

return

! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=1000 , loud=.TRUE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! call test_coulomb(ntests=100, loud=.FALSE.)
! call test_radial_grid(ntests=200)
! call test_radial_quadrature(ntests=100, loud=.FALSE.)
! return
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––

ngrid = 10

allocate(r1(ngrid))
allocate(wr1(ngrid))
allocate(y1(ngrid))
allocate(s1(ngrid))
allocate(projector1(ngrid)); allocate(projector2(ngrid))

call radial_grid(r=r1, wr=wr1, n=ngrid, addr2=.TRUE., quadr=1)
r1 = r1(ngrid:1:-1)
y1 = exp(-r1*r1)
call spline(r=r1, y=y1, n=size(r1), bound1=0.0_dp, boundn=0.0_dp, yspline=s1)
call pp_projector(l=1, m=0, r=r1, f=y1, s=s1, d12=(/1._dp,1._dp,1._dp/), p=projector1)
call pp_projector(l=1, m=0, r=r1, f=y1, s=s1, d12=(/2._dp,1._dp,1._dp/), p=projector2)
do i=1,ngrid
   print *, i, r1(i), y1(i), projector1(i), projector2(i), projector1(i)*projector2(i)
enddo


deallocate(r1)
deallocate(wr1)
deallocate(y1)
deallocate(s1)
deallocate(projector1); deallocate(projector2)

end program hallo