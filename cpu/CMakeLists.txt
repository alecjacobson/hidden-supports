cmake_minimum_required(VERSION 3.1)
project(hidden-supports)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_EMBREE      		 "Use Embree"         ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)
list(APPEND SRCFILES main.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::opengl igl::embree)
