# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( edgyeggs )

list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )
message (STATUS "Houdini HDK: ${Houdini_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Debug
         CACHE STRING "Choose the type of build, options are: ${CMAKE_CONFIGURATION_TYPES}" FORCE)
    set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${CMAKE_CONFIGURATION_TYPES})
endif ()

set (CMAKE_BUILD_TYPE DEBUG)
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    message("Buiding in debug mode")
endif (CMAKE_BUILD_TYPE MATCHES DEBUG)


# PCL time
find_package(PCL 1.3 REQUIRED COMPONENTS common features)
if (PCL_FOUND)
    message(STATUS "PCL found. Building PCAlign with Intel Fast Global Registration.")
	add_definitions(-DBUILD_WITH_INTEL_FGR)
	set(BUILD_WITH_INTEL_FGR TRUE)
endif()

# Cpd
include(ExternalProject)
ExternalProject_Add(Cpd
  SOURCE_DIR        "${CMAKE_SOURCE_DIR}/external/cpd"
  BINARY_DIR        "${CMAKE_SOURCE_DIR}/external/cpd/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)

# list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/cpd")
# find_package(Cpd REQUIRED)

set (edgyeggs_version_major       0)
set (edgyeggs_version_minor       2)
set (edgyeggs_version_patch       0)

# add_subdirectory (src/SOP_IGLDeform)
# add_subdirectory (src/SOP_IGLDiscreteGeo)
# add_subdirectory (src/SOP_IGLUVproject)
# add_subdirectory (src/SOP_ShapeOp)
add_subdirectory (${CMAKE_SOURCE_DIR}/external/cpd)
add_subdirectory (src/SOP_PCAlign)
