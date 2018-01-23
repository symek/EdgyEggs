# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( edgyeggs )

list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )

# Includes headers based projects
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)
include_directories(${CMAKE_SOURCE_DIR}/external/libigl/include)
include_directories(${CMAKE_SOURCE_DIR}/external/ShapeOp.0.1.0)
# ?
include_directories(${CMAKE_SOURCE_DIR}/external/eltopo/eltopo3d)
include_directories(${CMAKE_SOURCE_DIR}/external/eltopo/common)

# Add a library with two source files.
set( library_name SOP_EdgyEggs )

add_library( ${library_name} SHARED
    src/converters.cpp
    src/SOP_IGLDeform.cpp
    src/SOP_IGLDiscreteGeo.cpp
    src/SOP_IGLUVproject.cpp
    src/SOP_ShapeOp.cpp
    src/SOP_IGLMain.cpp
)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )

# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )