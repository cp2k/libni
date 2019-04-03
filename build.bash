gfortran -std=f2008 -J ./mods/ -I ./mods/ -c src/ni_types.f90 src/spherical_harmonics.f90 src/lebedev.f90 src/ni_fun.f90 src/ni_grid.f90 src/ni_gradients.f90 src/ni_module.f90 src/unit_test.f90 src/test_gradients.f90;
mv *.o mods/
gfortran -std=f2008 -J ./mods/ -I ./mods/ -C mods/*.o src/main.f90 -o out;
if [ $? -eq 0 ]; then
    ./out
else
    echo 'Build failed ☹️'
fi
