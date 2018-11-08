gfortran -C -o out spherical_harmonics.f90 lebedev.f90 grid.f90 gradients.f90 eddi.f90 unit_test.f90 test_gradients.f90 main.f90;
if [ $? -eq 0 ]; then
    ./out
else
    echo 'Build failed ☹️'
fi
