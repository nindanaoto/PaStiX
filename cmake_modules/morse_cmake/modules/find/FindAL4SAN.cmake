###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2019 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# @copyright 2017-2018 King Abdullah University of Science and Technology (KAUST). All rights reserved.
###
#
# - Find AL4SAN include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(AL4SAN
#               [REQUIRED]             # Fail with error if al4san is not found
#              )
#
# Results are reported in variables:
#
#  AL4SAN_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = AL4SAN
#  <XPREFIX> = <PREFIX>        for common case
#  <XPREFIX> = <PREFIX>_STATIC for static linking
#  <XPREFIX>_FOUND          ... set to 1 if module(s) exist
#  <XPREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
#  <XPREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
#  <XPREFIX>_LDFLAGS        ... all required linker flags
#  <XPREFIX>_LDFLAGS_OTHER  ... all other linker flags
#  <XPREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
#  <XPREFIX>_CFLAGS         ... all required cflags
#  <XPREFIX>_CFLAGS_OTHER   ... the other compiler flags

#=============================================================================
# Copyright 2019 Inria
# Copyright 2018 Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)


# Use pkg-config to detect include/library dirs (if pkg-config is available)
# --------------------------------------------------------------------------
include(FindPkgConfig)
if (AL4SAN_FIND_REQUIRED)
  find_package(PkgConfig QUIET REQUIRED)
else()
  find_package(PkgConfig QUIET)
endif()
if(PKG_CONFIG_EXECUTABLE)
  pkg_search_module(AL4SAN al4san)
  if (NOT AL4SAN_FIND_QUIETLY)
    if (AL4SAN_FOUND AND AL4SAN_LIBRARIES)
      message(STATUS "Looking for AL4SAN - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for AL4SAN - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing al4san.pc"
            "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  # check version if a precise version is required
  if (AL4SAN_FIND_VERSION_EXACT)
    if( NOT (AL4SAN_FIND_VERSION_MAJOR STREQUAL AL4SAN_VERSION_MAJOR) OR
              NOT (AL4SAN_FIND_VERSION_MINOR STREQUAL AL4SAN_VERSION_MINOR) )
      if(NOT AL4SAN_FIND_QUIETLY)
              message(FATAL_ERROR
              "AL4SAN version found is ${AL4SAN_VERSION_STRING}"
              "when required is ${AL4SAN_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (AL4SAN_FIND_VERSION_MAJOR STRGREATER AL4SAN_VERSION_MAJOR) OR
              (AL4SAN_FIND_VERSION_MINOR STRGREATER AL4SAN_VERSION_MINOR) )
      if(NOT AL4SAN_FIND_QUIETLY)
              message(FATAL_ERROR
                "AL4SAN version found is ${AL4SAN_VERSION_STRING}"
                "when required is ${AL4SAN_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  if (AL4SAN_FOUND AND AL4SAN_LIBRARIES)
    set(AL4SAN_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(AL4SAN)
  else()
    set(AL4SAN_FOUND_WITH_PKGCONFIG "FALSE")
  endif()
  
endif(PKG_CONFIG_EXECUTABLE)

if (AL4SAN_LIBRARIES)
  if (AL4SAN_LIBRARY_DIRS)
    foreach(dir ${AL4SAN_LIBRARY_DIRS})
      if ("${dir}" MATCHES "al4san")
        set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET AL4SAN_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(AL4SAN_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of AL4SAN library" FORCE)
  else()
    set(AL4SAN_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of AL4SAN library" FORCE)
  endif()
endif()
mark_as_advanced(AL4SAN_DIR_FOUND)

# check that AL4SAN has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND AL4SAN_FOUND)
  find_package_handle_standard_args(AL4SAN DEFAULT_MSG
    AL4SAN_LIBRARIES)
else()
  find_package_handle_standard_args(AL4SAN DEFAULT_MSG
    AL4SAN_LIBRARIES
    AL4SAN_WORKS)
endif()
