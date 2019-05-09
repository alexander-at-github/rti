# Project rti

### Installing
````
cd <rti-base-directory>
mkdir build
cd build
cmake ..
cmake --build . --target rti_core -j <number-of-cores>
````

The build system will download and install OCE, Gmsh, TBB, and Embree locally and automatically. These external dependencies are saved in the directory `<rti-base-directory>/dependencies`.

### License
See the file [LICENSE.md].
