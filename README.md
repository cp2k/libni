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
