module ni_types
implicit none
SAVE
integer, parameter :: dp = SELECTED_REAL_KIND(14, 200)
real(kind=dp), parameter :: pi = 3.14159265358979323846264338_dp ! Pi

type :: type_fun
   real(kind=dp), dimension(:), allocatable :: r
   real(kind=dp), dimension(:), allocatable :: y, y1, y2, y3, y4, y5
end type type_fun

type :: type_grid
   real(kind=dp), dimension(:, :), allocatable :: r
   real(kind=dp), dimension(:),    allocatable :: w
   real(kind=dp), dimension(:, :), allocatable :: dw
end type type_grid

type :: ni_env
    type(type_fun), pointer :: fun
    type(type_grid), dimension(:), pointer :: grids
end type ni_env

contains
    
end module ni_types
