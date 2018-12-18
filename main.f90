program hallo 
USE lebedev, ONLY: dp
USE eddi, ONLY: type_atom, integration_twocenter, &
                read_nfun, pi, interpolation, spline, &
                integration_threecenter, kinetic_energy, &
                integration_onecenter, coulomb_integral, coulomb_integral_grid, &
                radial_integration, pp_projector, pp_nonloc,&
                forward_derivative_weights
USE grid, ONLY: grid_parameters, radial_grid, gauher, s3, ds3dmu
USE nao_unit, ONLY: test_onecenter, test_twocenter, test_threecenter, test_kinetic, &
                    test_coulomb,&
                    test_radial_weight_pos, test_radial_chebyherm, test_radial_weight_asc, &
                    test_forward_deriv_coeff, test_spline,&
                    test_derivative_point_on, test_derivative_point_off,&
                    test_interpolation,&
                    test_derivative_on, test_derivative_off
USE gradients, ONLY: jacobian, grad_onecenter, grad_onecenter_cart, grad_twocenter_fd, &
                     grad_twocenter
USE nao_grad_unit, ONLY: test_jacobian, test_twocenter_grad
implicit none
REAL(KIND=dp), DIMENSION(1000) :: r, y1, wr, y2
REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
REAL(KIND=dp), DIMENSION(1000,3) :: error
INTEGER :: l1, l2, m1, m2, c, n

call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)

y1 = exp(-0.2_dp * r**2)
y2 = exp(-0.5_dp * r**2)

c = 0
error = 0._dp
d12 = (/ .1_dp, -.5_dp, .1_dp /)
n = 130
l1 = 0; m1 = 0; l2 = 1; m2 = 1;
   open(unit=222, file='b')
   write(222, *) r, y1, y2, (/l1,l2/), (/m1,m2/), (/n, n/), d12, grad1
   close(222)
call grad_twocenter(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                  nshell=(/n, n/), d12=d12, grad=grad1)
call grad_twocenter_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l1,l2/), m=(/m1,m2/),&
                     nshell=(/n, n/), d12=d12, grad=grad2)
error(c, :) = abs(grad1-grad2)
if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
! if ( any( error(c, :)  .gt. 0.1_dp  )) then
print *, '   ', l1, m1, l2, m2
print *, 'e  ', grad1
print *, 'fd ', grad2
print *, 'ra ', error(c, :), '/3 = ', sum(error(c, :))/3._dp

end program hallo