module eddi
   USE lebedev, ONLY: lebedev_grid,&
                        get_number_of_lebedev_grid,&
                        dp
   USE grid_point, ONLY: type_grid_point
   implicit none

   public :: print_lebedev, print_spherical, print_total_grid

   contains
      subroutine print_lebedev(l)
         integer :: i, iter, l
         i = get_number_of_lebedev_grid(l=l)
         print *, '! Angular momentum qn = ', l
         print *, '! Lebedev grid          ', i

         ! print *, 'x', 'y', 'z', 'w'
         do iter=1,lebedev_grid(i)%n
            print *, lebedev_grid(i)%r(:, iter), lebedev_grid(i)%w(iter)
         enddo
      end subroutine print_lebedev

      subroutine print_spherical(nshells, router)
         implicit none
         REAL(KIND=dp) :: router
         INTEGER :: ishell, nshells
         do ishell=1,nshells
            ! print *, 'x', -1.0_dp+0.1_dp*ishell
            print *, radial_mapping(-1.0_dp+0.1_dp*ishell)
         enddo
      end subroutine print_spherical

      subroutine print_total_grid(nleb, nshell)
         implicit none
         INTEGER :: nleb, nshell, ileb
         INTEGER :: iterleb, itershell, cnt
         INTEGER :: ngrid = 0   ! number of points making up the grid
         TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
         REAL(KIND=dp), DIMENSION(nshell) :: rrange

         rrange = range(nshell)
         
         ileb = get_number_of_lebedev_grid(nleb)
         ngrid = lebedev_grid(ileb)%n * nshell
         print *, '! grid points: leb ', lebedev_grid(ileb)%n, '. shell', nshell
         print *, '!', ngrid
         print *, '!' // REPEAT('-', 79)
         ! print *, shape(lebedev_grid(ileb)%r(:, 1))
         allocate(thegrid(ngrid))
         cnt = 0
         do itershell=1,nshell
            do iterleb=1,lebedev_grid(ileb)%n
               cnt = cnt+1
               thegrid(cnt)%r = rrange(itershell)*lebedev_grid(ileb)%r(:, iterleb)
               thegrid(cnt)%weight = lebedev_grid(ileb)%w(iterleb)
            enddo
         enddo
         do iterleb=1,cnt
            print *, thegrid(iterleb)%r, thegrid(iterleb)%weight
         enddo
         deallocate(thegrid)

      end subroutine print_total_grid

      ! take values x in [-1, 1) and map them to (0, +infty)
      function radial_mapping(x) result(r)
         implicit none
         REAL(KIND=dp), INTENT(in) :: x
         REAL(KIND=dp) :: eps, alpha, r

         eps = 1.0_dp
         alpha = 0.6_dp
         ! r = eps/log(2.0_dp)*(1+x)**alpha*log(2/(1-x))
         r = eps*(1+x)/(1-x)
          
      end function radial_mapping

      ! Evaluate a gaussian e^(-a.x^2)
      function gaussian(a, x) result(g)
         implicit none
         REAL(KIND=dp) :: a, x, g
         g = exp(-a*x**2)

      end function gaussian

      ! Take input N and return an array of N values evenly spaced
      ! between (-1, 1)
      function range(N) result(r)
         implicit none
         INTEGER :: N, iter
         REAL(KIND=dp), DIMENSION(N) :: r
         REAL(KIND=dp)               :: spacing
         spacing = 2.0_dp/N
         do iter=1,N
            r(iter) = -1.0_dp + spacing * iter
         enddo
          
      end function range

end module eddi