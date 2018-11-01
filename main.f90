program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, coulomb_integral_grid, &
                radial_integration, pp_projector, pp_nonloc,&
                forward_derivative_weights
USE grid, ONLY: grid_parameters, radial_grid, gauher
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic, &
                    test_coulomb,&
                    test_radial_weight_pos, test_radial_chebyherm, test_radial_weight_asc, &
                    test_forward_deriv_coeff, test_spline,&
                    test_derivative_point_on, test_derivative_point_off,&
                    test_interpolation,&
                    test_derivative_on, test_derivative_off
implicit none

CHARACTER(len=*), PARAMETER :: fn1 = 'gaussian.grid'
CHARACTER(len=*), PARAMETER :: fn2 = 'gaussian.grid'

! Local variables (Projector)
REAL(KIND=dp) :: integral, o
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r1, wr1, y1, s1
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: projector1, projector2
INTEGER :: i, nang, nshell, ngrid

! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
! call test_radial_weight_pos(ntests=9)
! call test_radial_weight_asc(ntests=9)
! call test_radial_chebyherm(ntests=25, loud=.FALSE.)
! call test_forward_deriv_coeff()
! call test_interpolation(ntests=100)
! call test_spline(ntests=100)
! call test_derivative_point_on()
! call test_derivative_point_off()

! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=100 , loud=.FALSE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! call test_coulomb(ntests=50, loud=.FALSE.)
call test_derivative_on(ntests=100)
call test_derivative_off(ntests=100)
return
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––

ngrid = 50

allocate(r1(ngrid))
allocate(wr1(ngrid))
allocate(y1(ngrid))
allocate(s1(ngrid))
allocate(projector1(ngrid)); allocate(projector2(ngrid))

call radial_grid(r=r1, wr=wr1, n=ngrid, addr2=.TRUE., quadr=2)
y1 = exp(-r1*r1)

! call spline(r=r1, y=y1, n=size(r1), bound1=0.0_dp, boundn=0.0_dp, yspline=s1)
! call pp_projector(l=1, m=0, r=r1, f=y1, s=s1, d12=(/1._dp,1._dp,1._dp/), p=projector1)
! call pp_projector(l=1, m=0, r=r1, f=y1, s=s1, d12=(/2._dp,1._dp,1._dp/), p=projector2)
! do i=1,ngrid
!    print *, i, r1(i), y1(i), projector1(i), projector2(i), projector1(i)*projector2(i)
! enddo
! o = 0
! do i=1,50
!    o = integral
!    call pp_nonloc(rv=r1, v=y1, rp1=r1, p1=y1, rp2=r1, p2=y1,&
!                   d12=(/1._dp, 1._dp, 1._dp/), d13=(/1._dp, 1._dp, 1._dp/),&
!                   lmax=0, nrad=i*10, integral=integral)
!    print *, i, i*10, integral, integral/o
   
! enddo

deallocate(r1)
deallocate(wr1)
deallocate(y1)
deallocate(s1)
deallocate(projector1); deallocate(projector2)

end program hallo