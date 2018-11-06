gfortran -C -o out main.f90 eddi.f90 gradients.f90 lebedev.f90 grid.f90 unit_test.f90 spherical_harmonics.f90;
if [ $? -eq 0 ]; then
    ./out
else
    echo 'Build failed ☹️'
fi
