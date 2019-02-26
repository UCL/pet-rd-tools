cmake_minimum_required( VERSION 3.7.0 FATAL_ERROR)

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(BUILD_POSTFIX "_d")
else(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(BUILD_POSTFIX "")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

# Directory includes
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

# Set package information
set(PACKAGE_NAME "petmrrdtools")
set(PACKAGE_VERSION "1.1.1")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}-${PACKAGE_VERSION}")
set(PACKAGE_LIBNAME "${PACKAGE_NAME}-${PACKAGE_VERSION}")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(${PACKAGE_NAME} VERSION ${PACKAGE_VERSION})

make_directory(${PROJECT_BINARY_DIR}/config)

set(BUILDER_NAME ${CURRENT_USER})

configure_file(
  "${PROJECT_SOURCE_DIR}/EnvironmentInfo.h.in"
  "${PROJECT_BINARY_DIR}/config/EnvironmentInfo.h"
)

set(CMAKE_CXX_STANDARD 11)

#set(nmtools_DIR "${CMAKE_CURRENT_SOURCE_DIR}/inm")

include(nmtools-deps)

include_directories("${PROJECT_BINARY_DIR}/config")
include_directories(lib)
add_subdirectory(src)
