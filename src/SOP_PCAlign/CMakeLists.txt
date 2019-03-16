# Includes headers based projects
include_directories(${CMAKE_SOURCE_DIR}/external/eigen)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp)
include_directories(${CMAKE_SOURCE_DIR}/external/sparseicp/include)

# Add a library with two source files.
set( library_name SOP_PCAlign )

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
option(BUILD_WITH_CPD "Build with CPD" OFF)

IF ( BUILD_WITH_CPD )
	message(STATUS "Building with CPD and FGT.")
	# message(WARNING ${CMAKE_SOURCE_DIR}/install)
	# include_directories(${CMAKE_SOURCE_DIR}/external/cpd/include)
	# include_directories(${CMAKE_SOURCE_DIR}/external/fgt/include)
	# set (BUILD_WITH_CPD TRUE)

	# list( APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/install/cpd/lib/cmake/cpd")
	# list( APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/install/fgt/lib/cmake/fgt")

	# find_package(Fgt REQUIRED)
	# find_package(Cpd REQUIRED)

	# add_dependencies(${library_name} Fgt)
	# add_dependencies(${library_name} Cpd)
	# list(APPEND cpd_sources
	# 	${CMAKE_SOURCE_DIR}/external/fgt/src/cluster.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/direct.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/direct_tree.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/ifgt.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/openmp.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/transform.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/cluster-sequential.cpp
	#     ${CMAKE_SOURCE_DIR}/external/fgt/src/version.cpp.in
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/affine.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/gauss_transform.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/matrix.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/nonrigid.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/normalization.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/rigid.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/transform.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/utils.cpp
	#     ${CMAKE_SOURCE_DIR}/external/cpd/src/version.cpp.in
	# 	)
	# add_definitions(-DBUILD_WITH_CPD)
	# add_library(cpd_static STATIC ${cpd_sources} )
	# target_compile_options(cpd_static
	#         PUBLIC
	#         -std=c++11
	#         PRIVATE
	#         -Wall
	#         -pedantic
	#         # The no-nested-anon-types is for nanoflann
	#         -Wno-nested-anon-types
	#         )
endif()



if ( BUILD_WITH_CPD )
	# link_directories(${CMAKE_SOURCE_DIR}/install/cpd/lib)
	# link_directories(${CMAKE_SOURCE_DIR}/install/fgt/lib)
	# target_link_libraries(${library_name}  Cpd::Library-C++)
endif()

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name}  Houdini )

if (BUILD_WITH_INTEL_FGR)
	target_link_libraries( ${library_name} ${PCL_COMMON_LIBRARIES} ${PCL_FEATURES_LIBRARIES})
endif()
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )