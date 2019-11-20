# Project rti

### Installing

Create a build directory:
````
cd <rti-base-directory>
mkdir build
cd build
````
The following commands will download and build Boost, OCE, Gmsh, and Embree in the folder `<rti-base-directory>/dependencies`. This may take some time. You need to build them only once:
````
cmake ..
cmake --build . --target boost-external
cmake --build . --target vtk-external
cmake --build . --target embree-external
````

Build the rti project:
````
cmake --build . --target rti
````
### License
See the [License file](./LICENSE).
