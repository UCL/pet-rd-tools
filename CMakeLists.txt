cmake_minimum_required( VERSION 3.21.0 FATAL_ERROR)

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set(BUILD_POSTFIX "_d")
else(CMAKE_BUILD_TYPE STREQUAL Debug)
    set(BUILD_POSTFIX "")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

# Directory includes
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

# Set package information
set(PACKAGE_NAME "petrdtools")
set(PACKAGE_VERSION "2.2.0")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME "${PACKAGE_NAME}-${PACKAGE_VERSION}")
set(PACKAGE_LIBNAME "${PACKAGE_NAME}-${PACKAGE_VERSION}")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(${PACKAGE_NAME} VERSION ${PACKAGE_VERSION})

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/config")

set(BUILDER_NAME ${CURRENT_USER})

configure_file(
  "${PROJECT_SOURCE_DIR}/EnvironmentInfo.h.in"
  "${PROJECT_BINARY_DIR}/config/EnvironmentInfo.h"
)

set(CMAKE_CXX_STANDARD 11)

#set(nmtools_DIR "${CMAKE_CURRENT_SOURCE_DIR}/inm")

include(nmtools-deps)

if(APPLE)
  # Homebrew ITK's usage requirements inject a raw -I<CommandLineTools SDK>/usr/include
  # that precedes the Xcode sysroot and breaks libc++/boost; drop SDK dirs.
  get_directory_property(_incdirs INCLUDE_DIRECTORIES)
  list(FILTER _incdirs EXCLUDE REGEX "/SDKs/")
  set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${_incdirs}")
endif()

include_directories("${PROJECT_BINARY_DIR}/config")
include_directories(lib)
add_subdirectory(src)
