# Includes headers based projects
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp/include)
include_directories(${CMAKE_SOURCE_DIR}/external/cpd/include)

# Add a library with two source files.
set( library_name SOP_PCAlign )
# add_definitions(-DEIGEN_MALLOC_ALREADY_ALIGNED)
# Main sources
list(APPEND pcalign_sources
    SOP_PCAlign.cpp
    SOP_PCAlign.hpp
)
# Intel fgr related
if (BUILD_WITH_INTEL_FGR)
	include_directories(${PCL_INCLUDE_DIRS})
	link_directories(${PCL_LIBRARY_DIRS})
	add_definitions(${PCL_DEFINITIONS})
	list(APPEND pcalign_sources
	    intel_fgr.cpp
	    intel_fgr.hpp
	    pcl_fpfh.hpp
		)
endif()

add_library( ${library_name} SHARED
   ${pcalign_sources}
)

add_dependencies(${library_name} Cpd)

link_directories(${CMAKE_BINARY_DIR}/external/cpd)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )
#FIXME
target_link_libraries( ${library_name} ${CMAKE_BINARY_DIR}/external/cpd/libcpd.dylib )
# target_link_libraries( ${library_name} ${CMAKE_SOURCE_DIR}/external/cpd/build/libcpd.dylib )

if (BUILD_WITH_INTEL_FGR)
	target_link_libraries( ${library_name} ${PCL_COMMON_LIBRARIES} ${PCL_FEATURES_LIBRARIES})
endif()
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )