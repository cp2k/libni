gfortran -C -o out_test spherical_harmonics.f90 lebedev.f90 grid.f90 gradients.f90 eddi.f90 unit_test.f90 test_gradients.f90 test_suite.f90;
if [ $? -eq 0 ]; then
    ./out_test
else
    echo 'Build failed ☹️'
fi
