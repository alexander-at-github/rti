## external project base directory for all projects in this hierarchy
set_property(DIRECTORY PROPERTY EP_BASE ${M_DEPENDENCIES_DIR})

add_subdirectory(boost)
# Alternatively one can use OCC instead of OCE
add_subdirectory(oce)
add_subdirectory(gmsh)
add_subdirectory(embree)
add_subdirectory(vtk)

add_subdirectory(gperftools)
