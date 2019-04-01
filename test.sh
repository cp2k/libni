rm mods/*.o
gfortran -std=f2008 -J ./mods/ -I ./mods/ -c ni_types.f90 spherical_harmonics.f90 lebedev.f90 ni_fun.f90 ni_grid.f90 ni_module.f90 ni_gradients.f90 unit_test.f90 test_gradients.f90;
mv *.o mods/
gfortran -std=f2008 -J ./mods/ -I ./mods/ -C mods/*.o test_suite.f90 -o out_test;
if [ $? -eq 0 ]; then
    ./out_test | tee log.txt
else
    echo 'Build failed ☹️'
fi
