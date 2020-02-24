# Project rti

### Installing

Create a build directory:
````
cd <rti-base-directory>
mkdir build
cd build
````
The following commands will download and build Boost, Embree, and VTK in the folder `<rti-base-directory>/dependencies`. This may take some time. You need to build them only once:
````
cmake ..
cmake --build . --target boost-external
cmake --build . --target embree-external
cmake --build . --target vtk-external
````

Build the rti library:
````
cmake --build . --target rti
````

The library and CMake files will be saved in directories under `build/lib`.
The API is declared in the header file `build/include/rti/device.hpp`.

In order to use the rti library in CMake use `find_package(librtidevice)` in your CMakeLists.txt to obtain the target `librtidevice::librtidevice`, for example:

````
find_package (
  librtidevice REQUIRED
  PATHS "/path/to/your/rti/build/directory"
  )
````

### License
See the [License file](./LICENSE).

