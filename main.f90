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
USE gradients, ONLY: jacobian, grad_onecenter, grad_twocenter
USE nao_grad_unit, ONLY: test_jacobian
implicit none

REAL(KIND=dp), DIMENSION(1000) :: r, y1, wr, y2
REAL(KIND=dp), DIMENSION(3) :: grad, d12
INTEGER :: l, m, i

call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)
y1 = exp(-r**2)
y2 = exp(-0.2_dp * r**2)

do l=0,3
   do m=-l,l
      d12 = (/ 1._dp, 1._dp, -0.25_dp /)
      call grad_twocenter(r1=r, y1=y1, r2=r, y2=y2, l=(/l,l/), m=(/m,m/),&
                          nshell=(/100, 100/), d12=d12, grad=grad)
      ! where(abs(grad) .lt. 1._dp*epsilon(1._dp)) grad = 0._dp
      ! call grad_onecenter(r=r, y=y, l=l, m=m, nshell=50, grad=grad)
      ! where(abs(grad) .lt. 10._dp*epsilon(1._dp)) grad = 0._dp
      print *, l, m, grad
   enddo
   print *,
enddo

! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
! call test_radial_weight_pos(ntests=9)
! call test_radial_weight_asc(ntests=9)
! call test_radial_chebyherm(ntests=25, loud=.FALSE.)
! call test_forward_deriv_coeff()
! call test_interpolation(ntests=100)
! call test_spline(ntests=100)
! call test_derivative_point_on()
! call test_derivative_point_off()

! call test_derivative_on(ntests=100)
! call test_derivative_off(ntests=100)
! call test_onecenter(ntests=100, loud=.FALSE.)
! call test_twocenter(ntests=100, loud=.FALSE.)
! call test_threecenter(ntests=100 , loud=.FALSE.)
! call test_kinetic(ntests=100, loud=.FALSE.)
! call test_coulomb(ntests=50, loud=.FALSE.)
return
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
end program hallo