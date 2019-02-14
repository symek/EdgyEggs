# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( edgyeggs )

list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )
message (STATUS "Houdini HDK: ${Houdini_VERSION}")

set (edgyeggs_version_major       0)
set (edgyeggs_version_minor       2)
set (edgyeggs_version_patch       0)

add_subdirectory (src/SOP_IGLDeform)
add_subdirectory (src/SOP_IGLDiscreteGeo)
add_subdirectory (src/SOP_IGLUVproject)
add_subdirectory (src/SOP_ShapeOp)