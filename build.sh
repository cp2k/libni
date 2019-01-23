gfortran -std=f2008 -J ./mods/ -I ./mods/ -c spherical_harmonics.f90 lebedev.f90 ni_fun.f90 grid.f90 gradients.f90 eddi.f90 unit_test.f90 test_gradients.f90;
mv *.o mods/
gfortran -std=f2008 -J ./mods/ -I ./mods/ -C mods/*.o main.f90 -o out;
if [ $? -eq 0 ]; then
    ./out
else
    echo 'Build failed ☹️'
fi
