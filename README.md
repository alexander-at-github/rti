# rti &ndash; Flux Calculation Library

### Build Instructions

Create a build directory:
````
cd <rti-base-directory>
mkdir build
cd build
````
The following commands will download and build Embree and VTK in the folder `<rti-base-directory>/dependencies`. This may take some time. You need to build them only once:
````
cmake ..
cmake --build . --target embree-external
cmake --build . --target vtk-external
````

Build the rti library:
````
cmake --build . --target rti
````

If you want to use your own VTK installation, you can provide the path in your first invocation of cmake:

````
cmake -DVTK_DIR=/your/path/to/vtk ..
````

The library and CMake files will be saved in directories under `build/lib`.
The API is declared in the header file `build/include/rti/device.hpp`.

In order to use the rti library use `find_package(librtidevice)` in your CMakeLists.txt to obtain the target `librtidevice::librtidevice`, for example:

````
find_package (
  librtidevice REQUIRED
  PATHS "/path/to/your/rti/build/directory"
  )
````

### License
See the [License file](./LICENSE).

