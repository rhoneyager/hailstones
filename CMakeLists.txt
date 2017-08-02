cmake_minimum_required (VERSION 3.1)
project(hailstones C CXX)
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/CMakeRules"
  "${CMAKE_CURRENT_SOURCE_DIR}"
)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (NOT WIN32)
	if (NOT CMAKE_BUILD_TYPE)
	  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
	endif (NOT CMAKE_BUILD_TYPE)
endif (NOT WIN32)

include(addBaseProject)
include(addBoostUniform)
include(addapp)
include(addlib)
include(generaterc)
include(GetGitRevisionDescription)

get_git_head_revision(GITREFSPECVAR GITHASHVAR)
set(MAJOR 0)
set(MINOR 0)
set(REVISION 1)
set(hailstones_VERSION ${MAJOR}.${MINOR}.${REVISION})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
addBaseProject()
configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake-settings.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake-settings.h"
    )

# Add in script-generated includes
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

#find_package(Eigen3 REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


set (boost_uniform_libs program_options filesystem system)
addBoostUniform(${boost_uniform_libs})
include_directories(BEFORE SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
if(VTK_LIBRARIES)
  set(app-libs ${VTK_LIBRARIES})
else()
  set(app-libs vtkHybrid vtkWidgets)
endif()
if(WIN32 AND NOT CYGWIN)
else()
	set (app-libs ${app-libs}
		${Boost_FILESYSTEM_LIBRARY}
		${Boost_SYSTEM_LIBRARY}
		${Boost_PROGRAM_OPTIONS_LIBRARY}
		)
endif()

generaterc("prog1" "exe"
	"snowflakes"
	"Program to convert Andy's surface meshes into solid volume representations"
	 "${CMAKE_CURRENT_BINARY_DIR}/prog1.rc"
	 )
set (prog1_srcs prog1.cpp)
IF (WIN32 AND NOT CYGWIN)
SET(prog1_srcs ${prog1_srcs} 
	"${CMAKE_CURRENT_BINARY_DIR}/prog1.rc")
ENDIF()

add_executable(prog1 MACOSX_BUNDLE ${prog1_srcs})
target_link_libraries(prog1 ${app-libs})
addapp(prog1 Main)
storebin(prog1)

generaterc("prog2a" "exe"
	"snowflakes"
	"Program to read solid volume meshes and decimate according to a table of known masses for the same particles"
	 "${CMAKE_CURRENT_BINARY_DIR}/prog2a.rc"
	 )
set (prog2a_srcs prog2a.cpp)
IF (WIN32 AND NOT CYGWIN)
SET(prog2a_srcs ${prog2a_srcs} 
	"${CMAKE_CURRENT_BINARY_DIR}/prog2a.rc")
ENDIF()

add_executable(prog2a MACOSX_BUNDLE ${prog2a_srcs})
target_link_libraries(prog2a ${app-libs})
addapp(prog2a Main)
storebin(prog2a)

generaterc("prog3a" "exe"
	"snowflakes"
	"Program to read volume meshes and write DDSCAT shape files"
	 "${CMAKE_CURRENT_BINARY_DIR}/prog3a.rc"
	 )
set (prog3a_srcs prog3a.cpp)
IF (WIN32 AND NOT CYGWIN)
SET(prog3a_srcs ${prog3a_srcs} 
	"${CMAKE_CURRENT_BINARY_DIR}/prog3a.rc")
ENDIF()

add_executable(prog3a MACOSX_BUNDLE ${prog3a_srcs})
target_link_libraries(prog3a ${app-libs})
addapp(prog3a Main)
storebin(prog3a)

generaterc("prog3b" "exe"
	"snowflakes"
	"Program to read volume meshes and write DDSCAT shape files"
	 "${CMAKE_CURRENT_BINARY_DIR}/prog3b.rc"
	 )
set (prog3b_srcs prog3b.cpp refract.cpp refract.h)
IF (WIN32 AND NOT CYGWIN)
SET(prog3b_srcs ${prog3b_srcs} 
	"${CMAKE_CURRENT_BINARY_DIR}/prog3b.rc")
ENDIF()

add_executable(prog3b MACOSX_BUNDLE ${prog3b_srcs})
target_link_libraries(prog3b ${app-libs})
addapp(prog3b Main)
storebin(prog3b)
