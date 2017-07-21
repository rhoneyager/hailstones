cmake_minimum_required (VERSION 3.1)
project(hailstones C CXX)
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/CMakeRules"
  "${CMAKE_CURRENT_SOURCE_DIR}"
)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(addBaseProject)
include(addBoostUniform)
include(addapp)
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

set (prog1_srcs prog1.cpp)
IF (WIN32 AND NOT CYGWIN)
SET(srcs ${srcs} 
	"${CMAKE_CURRENT_BINARY_DIR}/app.rc")
ENDIF()

add_executable(prog1 ${srcs})
#target_link_libraries(prog1 libicedb)
addapp(prog1 Main)
storebin(prog1)
