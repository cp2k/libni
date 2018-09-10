module grid_point
    USE lebedev, ONLY: dp
    type :: type_grid_point
       REAL(KIND=dp), DIMENSION(3) :: r = 0.0_dp
       REAL(KIND=dp) :: weight = 0.0_dp
    end type type_grid_point

    public :: type_grid_point
    
end module grid_point