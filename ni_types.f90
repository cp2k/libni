module ni_types
implicit none
SAVE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp ! Pi

type :: type_fun
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: r
   REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: y, y1, y2, y3, y4, y5
end type type_fun

type :: type_grid
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: r
   REAL(KIND=dp), DIMENSION(:),    ALLOCATABLE :: w
   REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: dw
end type type_grid

type :: ni_env
    type(type_fun), POINTER :: fun
    type(type_grid), DIMENSION(:), POINTER :: grids
end type ni_env

contains
    
end module ni_types