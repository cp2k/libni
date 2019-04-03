# libni

The module `ni_module` contains routines, which allow to evaluate common QM integrals, where the basis functions are tabulated on a grid. `src/test_suite.f90` serves as documentation:

```fortran
! test_suite.f90
USE ni_module, ONLY: integration_twocenter,&
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
```

## Integrals
The usual call signature:

```fortran
subroutine integration_twocenter(l, m, nshell, d12, r1, y1, r2, y2, &
                                 spline1, spline2, integral)
```

- `l, m`: angular and magnetic quantum number (used for multiplying the functions with spherical harmonics). Could also be used to decide on the size of the integration grid.
- `d12`: distance vector between the centers
- `r1, y1; r2, y2`: the function y1, y2 tabulated at radial points r1, r2
- `spline1, spline2`: the splines of the functions, used for interpolation
- `integral`: the output

## Gradients

```fortran
! test_gradients.f90
USE ni_gradients, ONLY: jacobian, grad_twocenter, grad_twocenter_fd,&
                     grad_kinetic, grad_kinetic_fd,&
                     grad_coulomb, grad_coulomb_fd
```

## Function structs

The tabulated function can be prepared as a `type_fun` type:

```fortran
! ni_fun.f90
subroutine prepare_fun(r, f, fun)
   implicit none
   REAL(KIND=dp), DIMENSION(:) :: r, f
   TYPE(type_fun), POINTER :: fun

   fun%r = r
   fun%y = f
   ! We need the first through fifth derivatives
   call derivatives(r=r, y=f, y1=fun%y1, y2=fun%y2, y3=fun%y3)
   call spline(r, fun%y2, size(r), fun%y4)
   call spline(r, fun%y3, size(r), fun%y5)

   ! Now in principle this is not a unique assignment.
   ! The spline of the function is also the second derivative,
   ! The spline of the first derivative is also the third derivative
   ! The fourth and fifth derivative could be retrieved from another call
   ! to `derivatives`.
end subroutine prepare_fun
```

## Tests

There are tests for all routines, comparing analytical Gaussian integrals with the numerical integration scheme. The latest output:

```
 👌 test_radial - all weights positive - passed
 👌 test_radial - radii ascending order - passed
 👌 test_radial - cheby.herm - passed
 👌 test_spline - passed
 max. error:    5.2551394519997629E-007
 mean error:    5.2532541405195017E-007
 👌 test_interpolation - passed
 max. error:    2.0330073160962031E-011
 mean error:    2.0330073160962031E-011
 👌 test forward derivative coefficients - passed
 👌 test derivative - first - passed
 Average error   2.0734877921523409E-006
 Max error   1.1119069213850314E-005
 👌 test derivative - second - passed
 Average error   2.0734877921523409E-006
 Max error   1.1119069213850314E-005
 👌 test derivative - third - passed
 Average error   2.0734877921523409E-006
 Max error   1.1119069213850314E-005
 👌 test_derivative_point_on – radial grid - passed
 max. error:    1.1351726254208810E-005
 mean error:    2.1391453761181597E-006
 error of sum:    2.2222216466971645E-007
 👌 test_derivative_point_on – equally spaced grid - passed
 max. error:    1.8984180834541634E-007
 mean error:    6.1935172909896709E-008
 error of sum:    1.9534793782582938E-008
 👌 test_derivative_point_off – radial grid - passed
 max. error:    1.1458953580938935E-004
 mean error:    3.2296242554542527E-006
 error of sum:    1.3098482794390939E-008
 👌 test_derivative_point_off – equally spaced grid - passed
 max. error:    4.8449707018411357E-002
 mean error:    1.7156257798068885E-003
 error of sum:    2.9853664109350093E-005
 💣 test_derivative - on grid - failed
 max. error:    8.6298361075689459E+017
 mean error:    8.6298361075689459E+017
 👌  test_derivative - off grid - passed
 max. error:    3.4527185031314547E-005
 mean error:    3.4527185031314547E-005
 
 ------------------------------ Testing One-Center ------------------------------
 Mean error:    3.8182790262908382E-013
 ---------------------------- End Testing One-Center ----------------------------

 ------------------------------ Testing Two-Center ------------------------------
 Mean error:    0.0000000000000000
 ---------------------------- End Testing Two-Center ----------------------------

 ------------------------------ Testing Three-Center ------------------------------
 Mean error:    2.1002187544638673E-008
 ---------------------------- End Testing Three-Center ----------------------------
 
 
 THE TEST SUITE WILL OUTPUT TEST CASES FOR WHICH THE ERROR IS HIGH. IN THIS CASE
 INTEGRATEING A COMBINATION OF HIGH α GAUSSIANS DOESN'T WORK PERFECTLY.
  ------------------------------ Testing Kinetic energy ------------------------------
           2
 Exponents:    3.2131044185975841        3.3375705296524085
 Is:   0.12980211413251669
 Should:  0.12980050036408144
 Absolute Difference:    1.6137684352435233E-006
 Relative Error:    1.2432682699348163E-005

 --------------------------------------------------------------------------------
           3
 Exponents:    2.1796953626752340        3.3250087368354655
 Is:   0.13551770967891422
 Should:  0.13551577206623969

 Absolute Difference:    1.9376126745296318E-006
 Relative Error:    1.4298060255102385E-005

 --------------------------------------------------------------------------------
 Mean error:    1.1632930605539235E-005
 ---------------------------- End Testing Kinetic energy ----------------------------
 ------------------------------ Testing Coulomb ------------------------------
 Mean error:    1.1597608758536400E-006
 ---------------------------- End Testing Coulomb ----------------------------
  👌 test_twocenter_grad - passed
```
