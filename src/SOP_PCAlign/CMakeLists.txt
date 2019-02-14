# Includes headers based projects
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp/include)

# Add a library with two source files.
set( library_name SOP_PCAlign )

add_library( ${library_name} SHARED
    SOP_PCAlign.cpp
    SOP_PCAlign.hpp
)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )

# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )