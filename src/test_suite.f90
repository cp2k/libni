program test_suite 
use ni_types, only: dp, pi, type_grid, type_fun, ni_env
use ni_module, only: integration_twocenter,&
                read_nfun,&
                interpolation,&
                spline,&
                integration_threecenter,&
                kinetic_energy,&
                integration_onecenter,&
                coulomb_integral,&
                coulomb_integral_grid,&
                radial_integration,&
                pp_projector,&
                pp_nonloc,&
                forward_derivative_weights
use ni_grid, only: grid_parameters, radial_grid, gauher
use ni_fun, only: prepare_fun
use nao_unit, only: test_onecenter,&
                    test_twocenter,&
                    test_threecenter,&
                    test_kinetic,&
                    test_coulomb,&
                    test_radial_weight_pos,&
                    test_radial_chebyherm,&
                    test_radial_weight_asc,&
                    test_forward_deriv_coeff,&
                    test_spline,&
                    test_derivative_point_on,&
                    test_derivative_point_off,&
                    test_interpolation,&
                    test_derivative_on,&
                    test_derivative_off,&
                    test_derivatives,&
                    test_onecenter_acc,&
                    test_twocenter_acc,&
                    test_kinetic_acc
use nao_grad_unit, only: test_jacobian,&
                         test_twocenter_grad,&
                         test_kinetic_grad,&
                         test_kinetic_fd,&
                         test_twocenter_fd,&
                         test_coulomb_grad
implicit none
real(kind=dp) :: timer_start, timer_stop
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
! Benchmarks
call cpu_time(timer_start)

call test_onecenter_acc()
call test_twocenter_acc()
call test_kinetic_acc()

call test_twocenter_fd()

! --  Tests concerning radial grids and functions on those grids -- !
call test_radial_weight_pos(ntests=9)
call test_radial_weight_asc(ntests=9)
call test_radial_chebyherm(ntests=9, loud=.false.)
call test_spline(ntests=10)
call test_interpolation(ntests=1)

! ! --  Tests concerning derivatives -- !
call test_forward_deriv_coeff()
call test_derivatives()
call test_derivative_point_on()
call test_derivative_point_off()
call test_derivative_on(ntests=1)
call test_derivative_off(ntests=1)

! ! --  Tests concerning the integrals -- !
call test_onecenter(ntests=10, loud=.false.)
call test_twocenter(ntests=3, loud=.false.)
call test_threecenter(ntests=3 , loud=.false.)
call test_kinetic(ntests=3, loud=.false.)
call test_coulomb(ntests=3, loud=.false.)

! call test_twocenter_grad(loud=.false.)
! call test_kinetic_grad(loud=.true.)
! call test_coulomb_grad(loud=.false.)    
! call test_jacobian()

call cpu_time(timer_stop)
print *, 'Total time: ', timer_stop - timer_start
! ––––––––––––––––––––––––––––––––– Test suite –––––––––––––––––––––––––––––––––
end program test_suite

 ! 👌 test_radial - all weights positive - passed
 ! 👌 test_radial - radii ascending order - passed
 ! 👌 test_radial - cheby.herm - passed
 ! 👌 test_spline - passed
 ! max. error:    5.2546253072522031E-007
 ! mean error:    5.2535231079479282E-007
 ! 👌 test_interpolation - passed
 ! max. error:    3.8854068332578025E-012
 ! mean error:    3.8854068332578025E-012

 ! 👌 test forward derivative coefficients - passed
 ! 👌 test derivative - first - passed
 ! Average error   2.0734877921523409E-006
 ! Max error   1.1119069213850314E-005
 ! 👌 test derivative - second - passed
 ! Average error   2.0734877921523409E-006
 ! Max error   1.1119069213850314E-005
 ! 👌 test derivative - third - passed
 ! Average error   2.0734877921523409E-006
 ! Max error   1.1119069213850314E-005
 ! 👌 test_derivative_point_on – radial grid - passed
 ! max. error:    1.1351727411227652E-005
 ! mean error:    2.1391453623379814E-006
 ! error of sum:    2.2222213380551636E-007
 ! 👌 test_derivative_point_on – equally spaced grid - passed
 ! max. error:    8.9488852204707846E-008
 ! mean error:    3.3451227264683783E-008
 ! error of sum:    1.5149770549172104E-008
 ! 👌 test_derivative_point_off – radial grid - passed
 ! max. error:    4.8571690045776984E-004
 ! mean error:    9.2852862049533019E-006
 ! error of sum:    4.2502808304080525E-007
 ! 👌 test_derivative_point_off – equally spaced grid - passed
 ! max. error:    1.9034385756163846E-003
 ! mean error:    7.3650190031037660E-005
 ! error of sum:    1.1229037548912591E-005
 !  👌 test_derivative - on grid - passed
 ! max. error:    2.2464241946680770E-006
 ! mean error:    2.2464241946680770E-006
 !  👌  test_derivative - off grid - passed
 ! max. error:    1.7682772601954857E-008
 ! mean error:    1.7682772601954857E-008
