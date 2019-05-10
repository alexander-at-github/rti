# Project rti

### Installing

Create a build directory
````
cd <rti-base-directory>
mkdir build
cd build
````
The following commands download and build Boost, OCE, Gmsh, TBB, and Embree in the folder `<rti-base-directory>/dependencies`. You need to build them only once.
````
cmake ..
cmake --build . --target boost-external
cmake --build . --target gmsh-external
cmake --build . --target embree-external
````

Build the rti project
````
cmake --build . --target rti
````
### License
See the [License file](./LICENSE).
