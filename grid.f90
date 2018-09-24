module grid
   USE lebedev, ONLY: dp, lebedev_grid
   implicit none
   type :: type_grid_point
      REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
      REAL(KIND=dp) :: weight = 0.0_dp
   end type type_grid_point
   REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

   public :: type_grid_point, build_onecenter_grid, build_twocenter_grid, build_threecenter_grid
    
   contains

   subroutine build_onecenter_grid(ileb, nshell, thegrid)
      implicit none
      INTEGER, intent(in) :: ileb, nshell
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      INTEGER :: cnt, iterrad, iterang, i, iterileb
      REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
      REAL(KIND=dp) :: alpha

      cnt = 0

      alpha = pi/REAL(nshell+1, dp)
      do iterang=1, lebedev_grid(ileb)%n
         tangw = lebedev_grid(ileb)%w(iterang)
         do iterrad=1, nshell
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb)%r(:, iterang)

            ! WEIGHTS
            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw
         enddo !iterrad
      enddo !iterang
   end subroutine build_onecenter_grid

   ! thegrid: the integration grid
   ! gr1, gy1: the grid `f1` is given on
   ! gr2, gy2: the grid `f2` is given on
   subroutine build_twocenter_grid(ileb, nshell, displacement, thegrid, &
                                   gr1, gy1, gr2, gy2)
      implicit none
      INTEGER, DIMENSION(2), intent(in) :: ileb, nshell
      REAL(KIND=dp), DIMENSION(3), intent(in) :: displacement
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                              gr2, gy2
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      INTEGER :: cnt, iterrad, iterang, i, iterileb
      REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
      REAL(KIND=dp) :: alpha, mu, s1, s2, p, ri, rj, R

      cnt = 0
      iterileb = 1
      R = sqrt(sum(displacement**2))

      alpha = pi/REAL(nshell(1)+1, dp)
      do iterang=1, lebedev_grid(ileb(1))%n
         tangw = lebedev_grid(ileb(1))%w(iterang)
         do iterrad=1, nshell(1)
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb(1))%r(:, iterang)

            ! WEIGHTS
            ! nuclear partition
            if (R /= 0.0_dp) then         	
	            ri = sqrt(sum((thegrid(cnt)%r**2)))
	            rj = sqrt(sum(((thegrid(cnt)%r-displacement)**2)))

	            mu = (ri-rj)/R
	            s1 = s3(mu)
	            s2 = s3(-mu)
	            p = s1/(s1+s2)
            else
            	p = 0.5_dp
            endif
            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw * p
         	if (isnan(thegrid(cnt)%weight)) then
            	print *, p, ri, rj
         		stop '"x" is a NaN'
      		endif
         enddo !iterrad
      enddo !iterang

      alpha = pi/REAL(nshell(2)+1, dp)
      do iterang=1, lebedev_grid(ileb(2))%n
         tangw = lebedev_grid(ileb(2))%w(iterang)
         do iterrad=1, nshell(2)
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang)

            ! WEIGHTS
            ! nuclear partition
            if (R /= 0.0_dp) then   
	            ri = sqrt(sum((thegrid(cnt)%r**2)))
	            rj = sqrt(sum(((thegrid(cnt)%r-displacement)**2)))

	            mu = (ri-rj)/R
	            s1 = s3(mu)
	            s2 = s3(-mu)
	            p = s2/(s1+s2)
            else
            	p = 0.5_dp
	         endif

            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw * p
         enddo !iterrad
      enddo !iterang
   end subroutine build_twocenter_grid

   ! thegrid: the integration grid
   subroutine build_threecenter_grid(ileb, nshell, d12, d13, thegrid, &
                                     gr1, gy1, gr2, gy2, gr3, gy3)
      implicit none
      INTEGER, DIMENSION(3), intent(in) :: ileb, nshell
      REAL(KIND=dp), DIMENSION(3), intent(in) :: d12, d13
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, intent(in) :: gr1, gy1, &
                                                              gr2, gy2, &
                                                              gr3, gy3
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      INTEGER :: cnt, iterrad, iterang, i, iterileb
      REAL(KIND=dp) :: tradw, tangw, tsin, tcos, targ, tr
      REAL(KIND=dp) :: alpha, tP1, tP2, tP3, sP, p1, p2, p3
      REAL(KIND=dp) :: mu12, mu13, mu23, s12, s13, s23, s21, s31, s32
      REAL(KIND=dp) :: r1, r2, r3, R12, R13, R23

      cnt = 0
      iterileb = 1
      R12 = sqrt(sum(d12**2))
      R13 = sqrt(sum(d13**2))
      R23 = sqrt(sum((d13-d12)**2))

      alpha = pi/REAL(nshell(1)+1, dp)
      do iterang=1, lebedev_grid(ileb(1))%n
         tangw = lebedev_grid(ileb(1))%w(iterang)
         do iterrad=1, nshell(1)
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb(1))%r(:, iterang)

            ! WEIGHTS
            ! nuclear partition
            r1 = sqrt(sum((thegrid(cnt)%r**2)))
            r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
            r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

            mu12 = (r1-r2)/R12
            mu13 = (r1-r3)/R13
            mu23 = (r2-r3)/R23

            s12 = 1.0_dp-z(mu12)
            s21 = 1.0_dp+z(mu12)
            s13 = 1.0_dp-z(mu13)
            s31 = 1.0_dp+z(mu13)
            s32 = 1.0_dp+z(mu23)
            s23 = 1.0_dp-z(mu23)

            tP1 = s12*s13
            tP2 = s21*s23
            tP3 = s31*s32

            sP = tP1+tP2+tP3

            p1 = tP1/sP
            p2 = tP2/sP
            p3 = tP3/sP

            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw * p1
         enddo !iterrad
      enddo !iterang

      alpha = pi/REAL(nshell(2)+1, dp)
      do iterang=1, lebedev_grid(ileb(2))%n
         tangw = lebedev_grid(ileb(2))%w(iterang)
         do iterrad=1, nshell(2)
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb(2))%r(:, iterang)

            ! WEIGHTS
            ! nuclear partition
            r1 = sqrt(sum((thegrid(cnt)%r**2)))
            r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
            r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

            mu12 = (r1-r2)/R12
            mu13 = (r1-r3)/R13
            mu23 = (r2-r3)/R23

            s12 = 1.0_dp-z(mu12)
            s21 = 1.0_dp+z(mu12)
            s13 = 1.0_dp-z(mu13)
            s31 = 1.0_dp+z(mu13)
            s32 = 1.0_dp+z(mu23)
            s23 = 1.0_dp-z(mu23)

            tP1 = s12*s13
            tP2 = s21*s23
            tP3 = s31*s32

            sP = tP1+tP2+tP3

            p1 = tP1/sP
            p2 = tP2/sP
            p3 = tP3/sP

            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw * p2
         enddo !iterrad
      enddo !iterang

      alpha = pi/REAL(nshell(3)+1, dp)
      do iterang=1, lebedev_grid(ileb(3))%n
         tangw = lebedev_grid(ileb(3))%w(iterang)
         do iterrad=1, nshell(3)
            cnt = cnt+1
            ! COORDINATE
            targ = REAL(iterrad, dp)*alpha
            tcos = cos(targ)
            tr = (1+tcos)/(1-tcos)

            thegrid(cnt)%r = tr*lebedev_grid(ileb(3))%r(:, iterang)

            ! WEIGHTS
            ! nuclear partition
            r1 = sqrt(sum((thegrid(cnt)%r**2)))
            r2 = sqrt(sum(((thegrid(cnt)%r-d12)**2)))
            r3 = sqrt(sum(((thegrid(cnt)%r-d13)**2)))

            mu12 = (r1-r2)/R12
            mu13 = (r1-r3)/R13
            mu23 = (r2-r3)/R23

            s12 = 1.0_dp-z(mu12)
            s21 = 1.0_dp+z(mu12)
            s13 = 1.0_dp-z(mu13)
            s31 = 1.0_dp+z(mu13)
            s32 = 1.0_dp+z(mu23)
            s23 = 1.0_dp-z(mu23)

            tP1 = s12*s13
            tP2 = s21*s23
            tP3 = s31*s32

            sP = tP1+tP2+tP3

            p1 = tP1/sP
            p2 = tP2/sP
            p3 = tP3/sP

            ! radial
            tradw = 2.0_dp*alpha*sin(targ)*tr**2/(1.0_dp-tcos)**2

            thegrid(cnt)%weight = tangw * tradw * p3
         enddo !iterrad
      enddo !iterang
   end subroutine build_threecenter_grid

   function h(mu)
      implicit none
      REAL(KIND=dp) :: mu, h
      h = 1.5_dp*mu-0.5_dp*mu**3
   end function h

   function z(mu)
   	implicit none
   	REAL(KIND=dp) :: mu, mua, mua2, mua4, mua6, z
   	mua = (mu/0.64_dp)
   	mua2 = mua*mua
   	mua4 = mua2*mua2
   	mua6 = mua4*mua2
   	z = mua*(35.0_dp*(1-mua2)+21.0_dp*mua4-5.0_dp*mua6)/16.0_dp
	end function z

   function s3(mu)
      implicit none
      REAL(KIND=dp) :: mu, s3
      ! s3 = 0.5_dp*(1-h(h(h(mu))))
      s3 = 0.5_dp*(1-z(mu))
   end function s3

   subroutine grid_parameters(atom, nleb, nshell)
      implicit none
      INTEGER :: atom, nleb, nshell

      SELECT CASE (atom)
      CASE (0)
         nshell = 5
         nleb = 10
      CASE (1:2)
         nshell = 35
         nleb = 302
      CASE (3:10)
         nshell = 40
         nleb = 590
      CASE (11:18)
         nshell = 45
         nleb = 590
      CASE (19:)
         nshell = 50
         nleb = 590
      CASE DEFAULT
         nshell = 10
         nleb = 302
      END SELECT
   end subroutine grid_parameters

   subroutine print_grid(thegrid)
      implicit none
      TYPE(type_grid_point), DIMENSION(:), ALLOCATABLE :: thegrid
      INTEGER :: i
      do i=1,size(thegrid)
         print *, thegrid(i)%r, thegrid(i)%weight
      enddo 
   end subroutine print_grid

end module grid