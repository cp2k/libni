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
USE nao_grad_unit, ONLY: test_jacobian
implicit none

REAL(KIND=dp), DIMENSION(1000) :: r, y1, wr, y2
REAL(KIND=dp), DIMENSION(3) :: grad1, grad2, d12
REAL(KIND=dp), DIMENSION(100,3) :: error
INTEGER :: l, m, i, c

call radial_grid(r=r, wr=wr, n=size(r), addr2=.TRUE., quadr=1)
open(unit=111, file='partition')
do i=1,size(r)
   write (111,*) i, r(i), s3( r(i) ), ds3dmu( r(i) )
enddo
close(111)
! y1 = exp(-r**2)
! y2 = exp(-0.2_dp * r**2)

! c = 0
! error = 0._dp
! do l=0,6
!    do m=-l,l
!       c = c+1
!       d12 = (/ .1_dp, .1_dp, .1_dp /)
!       call grad_twocenter(r1=r, y1=y1, r2=r, y2=y2, l=(/l,l/), m=(/m,m/),&
!                           nshell=(/100, 100/), d12=d12, grad=grad1)
!       print *, 'e  ', l, m, grad1

!       call grad_twocenter_fd(r1=r, y1=y1, r2=r, y2=y2, l=(/l,l/), m=(/m,m/),&
!                              nshell=(/100, 100/), d12=d12, grad=grad2)
!       error(c, :) = abs(grad1-grad2)
!       if(grad2(1) .ne. 0._dp) error(c, 1) = error(c, 1)/abs(grad2(1))
!       if(grad2(2) .ne. 0._dp) error(c, 2) = error(c, 2)/abs(grad2(2))
!       if(grad2(3) .ne. 0._dp) error(c, 3) = error(c, 3)/abs(grad2(3))
!       print *, 'fd ', l, m, grad2
!       print *, 'ra ', l, m, error(c, :)
!       print *,
!       ! where(abs(grad) .lt. 1._dp*epsilon(1._dp)) grad = 0._dp
!       ! call grad_onecenter_cart(r=r, y=y1, l=l, m=m, nshell=100, grad=grad)
!       ! where(abs(grad) .lt. 10._dp*epsilon(1._dp)) grad = 0._dp
!       ! print *, l, m, grad

!       ! call grad_onecenter(r=r, y=y1, l=l, m=m, nshell=100, grad=grad)
!       ! where(abs(grad) .lt. 10._dp*epsilon(1._dp)) grad = 0._dp
!       ! print *, l, m, grad
!    enddo
!    print *,
! enddo

! print *, 'error', sum(error, 1)/16._dp

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