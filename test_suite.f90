program test_suite 
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
USE nao_grad_unit, ONLY: test_jacobian, test_twocenter_grad
implicit none

! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
call test_twocenter_grad()
call test_jacobian()
return

call test_radial_weight_pos(ntests=9)
call test_radial_weight_asc(ntests=9)
call test_radial_chebyherm(ntests=25, loud=.FALSE.)
call test_forward_deriv_coeff()
call test_interpolation(ntests=100)
call test_spline(ntests=100)
call test_derivative_point_on()
call test_derivative_point_off()


call test_derivative_on(ntests=100)
call test_derivative_off(ntests=100)
call test_onecenter(ntests=100, loud=.FALSE.)
call test_twocenter(ntests=100, loud=.FALSE.)
call test_threecenter(ntests=25 , loud=.FALSE.)
call test_kinetic(ntests=100, loud=.FALSE.)
call test_coulomb(ntests=50, loud=.FALSE.)
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
end program test_suite