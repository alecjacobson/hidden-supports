cmake_minimum_required(VERSION 3.1)
set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
message(${ROOT})
project(shaderpipeline)
option(USE_SOLUTION "Use solution"  OFF)
# Add your project files
include_directories("include/")
# create executable
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# OpenGL library
if (NOT CMAKE_VERSION VERSION_LESS "3.11")
  cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)

if(USE_SOLUTION)
  add_definitions("-DUSE_SOLUTION")
endif()

if(NOT TARGET glfw)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory(${ROOT}/../glfw glfw)
endif()

add_subdirectory(${ROOT}/glad glad)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
#option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
#option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
#option(LIBIGL_WITH_EMBREE      		   "Use Embree"         ON)
option(GAUSS_USE_UI "Use Gauss' UI" OFF)

find_package(LIBIGL REQUIRED QUIET)
find_package(Gauss REQUIRED QUIET)

add_executable(${PROJECT_NAME} main-1.cpp)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${OPENGL_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/../eigen/)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/../glfw/include)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/glad/include)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${Gauss_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} glfw ${OPENGL_gl_LIBRARY} glad igl::core)
