# Add a library with two source files.
set( library_name SOP_Eltopo )

add_library( ${library_name} SHARED
    src/converters.cpp
    src/SOP_Eltopo.cpp
    src/SOP_Eltopo.hpp
)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )

# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )