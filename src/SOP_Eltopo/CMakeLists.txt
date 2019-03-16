# Add a library with two source files.
set( library_name SOP_Eltopo )


file(GLOB eltopo_src
        "../../external/eltopo/eltopo3d/*.h"
        "../../external/eltopo/eltopo3d/*.cpp"
        "../../external/eltopo/common/*.cpp"
        "../../external/eltopo/common/*.h"
        "../../external/eltopo/common/newsparse/*.h"
        "../../external/eltopo/common/tunicate/*.cpp"
        "../../external/eltopo/common/newsparse/*.h"
        "../../external/eltopo/common/tunicate/*.cpp"
        )

include_directories(
        "../../external/eltopo/eltopo3d/"
        "../../external/eltopo/common/"
        "../../external/eltopo/common/newsparse"
        "../../external/eltopo/common/tunicate"
        )

add_library( eltopo
        ${eltopo_src})

add_library( SOP_Eltopo SHARED
    SOP_Eltopo.cpp
    SOP_Eltopo.hpp
)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( eltopo "-framework Accelerate")
target_link_libraries( SOP_Eltopo Houdini )
target_link_libraries( SOP_Eltopo eltopo )

# Configure several common target properties, such as its output directory.
houdini_configure_target(SOP_Eltopo )