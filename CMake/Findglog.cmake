# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2015 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: alexs.mac@gmail.com (Alex Stewart)
#

# FindGlog.cmake - Find Google glog logging library.
#
# This module defines the following variables:
#
# GLOG_FOUND: TRUE iff glog is found.
# GLOG_INCLUDE_DIRS: Include directories for glog.
# GLOG_LIBRARIES: Libraries required to link glog.
# FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION: True iff the version of glog found
#                                           was built & installed / exported
#                                           as a CMake package.
#
# The following variables control the behaviour of this module:
#
# GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION: TRUE/FALSE, iff TRUE then
#                           then prefer using an exported CMake configuration
#                           generated by glog > 0.3.4 over searching for the
#                           glog components manually.  Otherwise (FALSE)
#                           ignore any exported glog CMake configurations and
#                           always perform a manual search for the components.
#                           Default: TRUE iff user does not define this variable
#                           before we are called, and does NOT specify either
#                           GLOG_INCLUDE_DIR_HINTS or GLOG_LIBRARY_DIR_HINTS
#                           otherwise FALSE.
# GLOG_INCLUDE_DIR_HINTS: List of additional directories in which to
#                         search for glog includes, e.g: /timbuktu/include.
# GLOG_LIBRARY_DIR_HINTS: List of additional directories in which to
#                         search for glog libraries, e.g: /timbuktu/lib.
#
# The following variables are also defined by this module, but in line with
# CMake recommended FindPackage() module style should NOT be referenced directly
# by callers (use the plural variables detailed above instead).  These variables
# do however affect the behaviour of the module via FIND_[PATH/LIBRARY]() which
# are NOT re-called (i.e. search for library is not repeated) if these variables
# are set with valid values _in the CMake cache_. This means that if these
# variables are set directly in the cache, either by the user in the CMake GUI,
# or by the user passing -DVAR=VALUE directives to CMake when called (which
# explicitly defines a cache variable), then they will be used verbatim,
# bypassing the HINTS variables and other hard-coded search locations.
#
# GLOG_INCLUDE_DIR: Include directory for glog, not including the
#                   include directory of any dependencies.
# GLOG_LIBRARY: glog library, not including the libraries of any
#               dependencies.

# Reset CALLERS_CMAKE_FIND_LIBRARY_PREFIXES to its value when
# FindGlog was invoked.
macro(GLOG_RESET_FIND_LIBRARY_PREFIX)
  if (MSVC AND CALLERS_CMAKE_FIND_LIBRARY_PREFIXES)
    set(CMAKE_FIND_LIBRARY_PREFIXES "${CALLERS_CMAKE_FIND_LIBRARY_PREFIXES}")
  endif()
endmacro(GLOG_RESET_FIND_LIBRARY_PREFIX)

# Called if we failed to find glog or any of it's required dependencies,
# unsets all public (designed to be used externally) variables and reports
# error message at priority depending upon [REQUIRED/QUIET/<NONE>] argument.
macro(GLOG_REPORT_NOT_FOUND REASON_MSG)
  unset(GLOG_FOUND)
  unset(GLOG_INCLUDE_DIRS)
  unset(GLOG_LIBRARIES)
  # Make results of search visible in the CMake GUI if glog has not
  # been found so that user does not have to toggle to advanced view.
  mark_as_advanced(CLEAR GLOG_INCLUDE_DIR
                         GLOG_LIBRARY)

  glog_reset_find_library_prefix()

  # Note <package>_FIND_[REQUIRED/QUIETLY] variables defined by FindPackage()
  # use the camelcase library name, not uppercase.
  if (Glog_FIND_QUIETLY)
    message(STATUS "Failed to find glog - " ${REASON_MSG} ${ARGN})
  elseif (Glog_FIND_REQUIRED)
    message(FATAL_ERROR "Failed to find glog - " ${REASON_MSG} ${ARGN})
  else()
    # Neither QUIETLY nor REQUIRED, use no priority which emits a message
    # but continues configuration and allows generation.
    message("-- Failed to find glog - " ${REASON_MSG} ${ARGN})
  endif ()
  return()
endmacro(GLOG_REPORT_NOT_FOUND)

# Protect against any alternative find_package scripts for this library having
# been called previously (in a client project) which set GLOG_FOUND, but not
# the other variables we require / set here which could cause the search logic
# here to fail.
unset(GLOG_FOUND)

# -----------------------------------------------------------------
# By default, if the user has expressed no preference for using an exported
# glog CMake configuration over performing a search for the installed
# components, and has not specified any hints for the search locations, then
# prefer a glog exported configuration if available.
if (NOT DEFINED GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION
    AND NOT GLOG_INCLUDE_DIR_HINTS
    AND NOT GLOG_LIBRARY_DIR_HINTS)
  message(STATUS "No preference for use of exported glog CMake configuration "
    "set, and no hints for include/library directories provided. "
    "Defaulting to preferring an installed/exported glog CMake configuration "
    "if available.")
  set(GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION TRUE)
endif()

if (GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION)
  # Try to find an exported CMake configuration for glog, as generated by
  # glog versions > 0.3.4
  #
  # We search twice, s/t we can invert the ordering of precedence used by
  # find_package() for exported package build directories, and installed
  # packages (found via CMAKE_SYSTEM_PREFIX_PATH), listed as items 6) and 7)
  # respectively in [1].
  #
  # By default, exported build directories are (in theory) detected first, and
  # this is usually the case on Windows.  However, on OS X & Linux, the install
  # path (/usr/local) is typically present in the PATH environment variable
  # which is checked in item 4) in [1] (i.e. before both of the above, unless
  # NO_SYSTEM_ENVIRONMENT_PATH is passed).  As such on those OSs installed
  # packages are usually detected in preference to exported package build
  # directories.
  #
  # To ensure a more consistent response across all OSs, and as users usually
  # want to prefer an installed version of a package over a locally built one
  # where both exist (esp. as the exported build directory might be removed
  # after installation), we first search with NO_CMAKE_PACKAGE_REGISTRY which
  # means any build directories exported by the user are ignored, and thus
  # installed directories are preferred.  If this fails to find the package
  # we then research again, but without NO_CMAKE_PACKAGE_REGISTRY, so any
  # exported build directories will now be detected.
  #
  # To prevent confusion on Windows, we also pass NO_CMAKE_BUILDS_PATH (which
  # is item 5) in [1]), to not preferentially use projects that were built
  # recently with the CMake GUI to ensure that we always prefer an installed
  # version if available.
  #
  # NOTE: We use the NAMES option as glog erroneously uses 'google-glog' as its
  #       project name when built with CMake, but exports itself as just 'glog'.
  #       On Linux/OS X this does not break detection as the project name is
  #       not used as part of the install path for the CMake package files,
  #       e.g. /usr/local/lib/cmake/glog, where the <glog> suffix is hardcoded
  #       in glog's CMakeLists.  However, on Windows the project name *is*
  #       part of the install prefix: C:/Program Files/google-glog/[include,lib].
  #       However, by default CMake checks:
  #       C:/Program Files/<FIND_PACKAGE_ARGUMENT_NAME='glog'> which does not
  #       exist and thus detection fails.  Thus we use the NAMES to force the
  #       search to use both google-glog & glog.
  #
  # [1] http://www.cmake.org/cmake/help/v2.8.11/cmake.html#command:find_package
  find_package(glog QUIET
                    NAMES google-glog glog
                    NO_MODULE
                    NO_CMAKE_PACKAGE_REGISTRY
                    NO_CMAKE_BUILDS_PATH)
  if (glog_FOUND)
    message(STATUS "Found installed version of glog: ${glog_DIR}")
  else()
    # Failed to find an installed version of glog, repeat search allowing
    # exported build directories.
    message(STATUS "Failed to find installed glog CMake configuration, "
      "searching for glog build directories exported with CMake.")
    # Again pass NO_CMAKE_BUILDS_PATH, as we know that glog is exported and
    # do not want to treat projects built with the CMake GUI preferentially.
    find_package(glog QUIET
                      NAMES google-glog glog
                      NO_MODULE
                      NO_CMAKE_BUILDS_PATH)
    if (glog_FOUND)
      message(STATUS "Found exported glog build directory: ${glog_DIR}")
    endif(glog_FOUND)
  endif(glog_FOUND)

  set(FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION ${glog_FOUND})

  if (FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION)
    message(STATUS "Detected glog version: ${glog_VERSION}")
    set(GLOG_FOUND ${glog_FOUND})
    # glog wraps the include directories into the exported glog::glog target.
    set(GLOG_INCLUDE_DIR "")
    set(GLOG_LIBRARY glog::glog)
  else (FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION)
    message(STATUS "Failed to find an installed/exported CMake configuration "
      "for glog, will perform search for installed glog components.")
  endif (FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION)
endif(GLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION)

if (NOT GLOG_FOUND)
  # Either failed to find an exported glog CMake configuration, or user
  # told us not to use one.  Perform a manual search for all glog components.

  # Handle possible presence of lib prefix for libraries on MSVC, see
  # also GLOG_RESET_FIND_LIBRARY_PREFIX().
  if (MSVC)
    # Preserve the caller's original values for CMAKE_FIND_LIBRARY_PREFIXES
    # s/t we can set it back before returning.
    set(CALLERS_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
    # The empty string in this list is important, it represents the case when
    # the libraries have no prefix (shared libraries / DLLs).
    set(CMAKE_FIND_LIBRARY_PREFIXES "lib" "" "${CMAKE_FIND_LIBRARY_PREFIXES}")
  endif (MSVC)

  # Search user-installed locations first, so that we prefer user installs
  # to system installs where both exist.
  list(APPEND GLOG_CHECK_INCLUDE_DIRS
    /usr/local/include
    /usr/local/homebrew/include # Mac OS X
    /opt/local/var/macports/software # Mac OS X.
    /opt/local/include
    /usr/include)
  # Windows (for C:/Program Files prefix).
  list(APPEND GLOG_CHECK_PATH_SUFFIXES
    glog/include
    glog/Include
    Glog/include
    Glog/Include
    google-glog/include # CMake installs with project name prefix.
    google-glog/Include)

  list(APPEND GLOG_CHECK_LIBRARY_DIRS
    /usr/local/lib
    /usr/local/homebrew/lib # Mac OS X.
    /opt/local/lib
    /usr/lib)
  # Windows (for C:/Program Files prefix).
  list(APPEND GLOG_CHECK_LIBRARY_SUFFIXES
    glog/lib
    glog/Lib
    Glog/lib
    Glog/Lib
    google-glog/lib # CMake installs with project name prefix.
    google-glog/Lib)

  # Search supplied hint directories first if supplied.
  find_path(GLOG_INCLUDE_DIR
    NAMES glog/logging.h
    HINTS ${GLOG_INCLUDE_DIR_HINTS}
    PATHS ${GLOG_CHECK_INCLUDE_DIRS}
    PATH_SUFFIXES ${GLOG_CHECK_PATH_SUFFIXES})
  if (NOT GLOG_INCLUDE_DIR OR
      NOT EXISTS ${GLOG_INCLUDE_DIR})
    glog_report_not_found(
      "Could not find glog include directory, set GLOG_INCLUDE_DIR "
      "to directory containing glog/logging.h")
  endif (NOT GLOG_INCLUDE_DIR OR
    NOT EXISTS ${GLOG_INCLUDE_DIR})

  find_library(GLOG_LIBRARY NAMES glog
    HINTS ${GLOG_LIBRARY_DIR_HINTS}
    PATHS ${GLOG_CHECK_LIBRARY_DIRS}
    PATH_SUFFIXES ${GLOG_CHECK_LIBRARY_SUFFIXES})
  if (NOT GLOG_LIBRARY OR
      NOT EXISTS ${GLOG_LIBRARY})
    glog_report_not_found(
      "Could not find glog library, set GLOG_LIBRARY "
      "to full path to libglog.")
  endif (NOT GLOG_LIBRARY OR
    NOT EXISTS ${GLOG_LIBRARY})

  # Mark internally as found, then verify. GLOG_REPORT_NOT_FOUND() unsets
  # if called.
  set(GLOG_FOUND TRUE)

  # Glog does not seem to provide any record of the version in its
  # source tree, thus cannot extract version.

  # Catch case when caller has set GLOG_INCLUDE_DIR in the cache / GUI and
  # thus FIND_[PATH/LIBRARY] are not called, but specified locations are
  # invalid, otherwise we would report the library as found.
  if (GLOG_INCLUDE_DIR AND
      NOT EXISTS ${GLOG_INCLUDE_DIR}/glog/logging.h)
    glog_report_not_found(
      "Caller defined GLOG_INCLUDE_DIR:"
      " ${GLOG_INCLUDE_DIR} does not contain glog/logging.h header.")
  endif (GLOG_INCLUDE_DIR AND
    NOT EXISTS ${GLOG_INCLUDE_DIR}/glog/logging.h)
  # TODO: This regex for glog library is pretty primitive, we use lowercase
  #       for comparison to handle Windows using CamelCase library names, could
  #       this check be better?
  string(TOLOWER "${GLOG_LIBRARY}" LOWERCASE_GLOG_LIBRARY)
  if (GLOG_LIBRARY AND
      NOT "${LOWERCASE_GLOG_LIBRARY}" MATCHES ".*glog[^/]*")
    glog_report_not_found(
      "Caller defined GLOG_LIBRARY: "
      "${GLOG_LIBRARY} does not match glog.")
  endif (GLOG_LIBRARY AND
    NOT "${LOWERCASE_GLOG_LIBRARY}" MATCHES ".*glog[^/]*")

  glog_reset_find_library_prefix()

endif(NOT GLOG_FOUND)

# Set standard CMake FindPackage variables if found.
if (GLOG_FOUND)
  set(GLOG_INCLUDE_DIRS ${GLOG_INCLUDE_DIR})
  set(GLOG_LIBRARIES ${GLOG_LIBRARY})
endif (GLOG_FOUND)

# If we are using an exported CMake glog target, the include directories are
# wrapped into the target itself, and do not have to be (and are not)
# separately specified.  In which case, we should not add GLOG_INCLUDE_DIRS
# to the list of required variables in order that glog be reported as found.
if (FOUND_INSTALLED_GLOG_CMAKE_CONFIGURATION)
  set(GLOG_REQUIRED_VARIABLES GLOG_LIBRARIES)
else()
  set(GLOG_REQUIRED_VARIABLES GLOG_INCLUDE_DIRS GLOG_LIBRARIES)
endif()

# Handle REQUIRED / QUIET optional arguments.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(glog DEFAULT_MSG
  ${GLOG_REQUIRED_VARIABLES})

# Only mark internal variables as advanced if we found glog, otherwise
# leave them visible in the standard GUI for the user to set manually.
if (GLOG_FOUND)
  mark_as_advanced(FORCE GLOG_INCLUDE_DIR
                         GLOG_LIBRARY
                         glog_DIR) # Autogenerated by find_package(glog)
endif (GLOG_FOUND)
