gfortran -o out main.f90 eddi.f90 lebedev.f90 grid.f90;
if [ $? -eq 0 ]; then
    ./out
fi
