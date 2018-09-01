# AsFem
AsFem: **a** **s**imple **f**inite **e**lement **m**ethod program. A beginner
level's code for my personal interests.

# Installation
- Step-1: install MPI and PETSc, if you already have these two packages,
you can goto Step-2.
#### [mpich](https://www.mpich.org/) or [openmpi](https://www.open-mpi.org/)
```shell
make -j8
make install
```
#### [PETSc](https://www.mcs.anl.gov/petsc/)
you can use the default configuration for PETSc like:
```shell
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
```
then run
```shell
make all test
```

or compile PETSc with external packages like superlu_dist or mumps, you
can run:
```shell
./configure \
--prefix=$PETSC_DIR \
--download-hypre=1 \
--with-ssl=0 \
--with-debugging=no \
--with-pic=1 \
--with-shared-libraries=1 \
--with-cc=mpicc \
--with-cxx=mpicxx \
--with-fc=mpif90 \
--download-fblaslapack=1 \
--download-metis=1 \
--download-parmetis=1 \
--download-superlu_dist=1 \
--download-scalapack=1 \
--download-mumps=1 \
CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif77 F90=mpif90 \
CFLAGS='-fPIC -fopenmp -O3 -march=corei7-avx' \
CXXFLAGS='-fPIC -fopenmp -O3 -march=corei7-avx' \
FFLAGS='-fPIC -fopenmp -O3 -march=corei7-avx ' \
FCFLAGS='-fPIC -fopenmp -O3 -march=corei7-avx ' \
F90FLAGS='-fPIC -fopenmp -O3 -march=corei7-avx ' \
F77FLAGS='-fPIC -fopenmp -O3 -march=corei7-avx ' \
PETSC_DIR=`pwd`
```

- Step-2: set the correct path to PETSc and mpi

edit your CMakeLists.txt file, change the PETSc dir and MPI dir to
your own one, for example:
```cmake
set(PETSc "/your/path/to/petsc3.9.3")
```
then for mpi, you need to modify:
```cmake 
set(MPI "/your/path/to/mpi")
```
that's all you need to do.

- Step-3: compile AsFem

If all the path settings are correct, you can run:
```shell
CC=gcc CXX=g++ cmake CMakeLists.txt
```
after the configuration finished, it should generate a Makefile for you,
then run the following line:
```shell
make -j8
```
if everything is ok, you should find an executable file named 'ASFEM'
under the /bin folder.

- Step-4: run the simulation

Once you have the 'ASFEM' file, you can prepare an input file, for example
'test.i', then you can run it like:
```shell
./ASFEM test.i
```
or
```shell
mpirun -np 4 ./ASFEM test.i
```

# Contact
- Email: walkandthinker@gmail.com