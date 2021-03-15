clear
ls
rm *mod *exe *o
ls
gfortran -c MD.f95
ls
gfortran -c md_basic.f95
ls
gfortran MD.o md_basic.o -o MD_basic.exe
ls
./MD_basic.exe
