cmake_minimum_required (VERSION 3.3)

include (CheckCCompilerFlag)
include (CheckCXXCompilerFlag)

set( SANITIZERS
  TSAN ASAN LSAN MSAN UBSAN )

# ThreadSanitizer
set( TSAN_OPTION "-fsanitize=thread" )
set( TSAN_FLAGS  "-g -O1" )

# AddressSanitizer
set( ASAN_OPTION "-fsanitize=address" )
set( ASAN_FLAGS  " -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer" )

# LeakSanitizer
set( LSAN_OPTION "-fsanitize=leak" )
set( LSAN_FLAGS  "-fno-omit-frame-pointer" )

# MemorySanitizer
set( MSAN_OPTION "-fsanitize=memory" )
set( MSAN_FLAGS  "-fno-optimize-sibling-calls -fsanitize-memory-track-origins=2 -fno-omit-frame-pointer" )

# UndefinedBehaviour
set( UBSAN_OPTION "-fsanitize=undefined" )
set( UBSAN_FLAGS  "" )

get_property( languages GLOBAL PROPERTY ENABLED_LANGUAGES )

foreach(sanitizer ${SANITIZERS} )
  set(CMAKE_REQUIRED_FLAGS "${${sanitizer}_OPTION}" ) # Also needs to be a link flag for test to pass
  if ( "C" IN_LIST languages )
    check_c_compiler_flag(   "${${sanitizer}_OPTION} ${${sanitizer}_FLAGS}" HAVE_C_${sanitizer} )
    if ( HAVE_C_${sanitizer} )
      set( CMAKE_C_FLAGS_${sanitizer} "${${sanitizer}_OPTION} ${${sanitizer}_FLAGS}" CACHE STRING
        "Flags used by the C compiler during ${sanitizer} builds."
        FORCE )
    endif()
  endif()
  if ( "CXX" IN_LIST languages )
    check_cxx_compiler_flag( "${${sanitizer}_OPTION} ${${sanitizer}_FLAGS}" HAVE_CXX_${sanitizer} )
    if ( HAVE_CXX_${sanitizer} )
      set( CMAKE_CXX_FLAGS_${sanitizer} "${${sanitizer}_OPTION} ${${sanitizer}_FLAGS}" CACHE STRING
        "Flags used by the C++ compiler during ${sanitizer} builds."
        FORCE )
    endif()
  endif()

  if ( HAVE_C_${sanitizer} OR HAVE_CXX_${sanitizer} )
    set( CMAKE_Fortran_FLAGS_${sanitizer} "" CACHE STRING
      "Flags used by the Fortran compiler during ${sanitizer} builds."
      FORCE )
    set( CMAKE_EXE_LINKER_FLAGS_${sanitizer}
      "${${sanitizer}_OPTION}" CACHE STRING
      "Flags used for linking binaries during ${sanitizer} builds."
      FORCE )
    set( CMAKE_MODULE_LINKER_FLAGS_${sanitizer}
      "" CACHE STRING
      "Flags used by the linker the creation of module during ${sanitizer} builds."
      FORCE )
    set( CMAKE_SHARED_LINKER_FLAGS_${sanitizer}
      "${${sanitizer}_OPTION}" CACHE STRING
      "Flags used by the shared libraries linker during ${sanitizer} builds."
      FORCE )
    set( CMAKE_STATIC_LINKER_FLAGS_${sanitizer}
      "" CACHE STRING
      "Flags used by the static libraries linker during ${sanitizer} builds."
      FORCE )

    mark_as_advanced(
      CMAKE_C_FLAGS_${sanitizer}
      CMAKE_CXX_FLAGS_${sanitizer}
      CMAKE_Fortran_FLAGS_${sanitizer}
      CMAKE_EXE_LINKER_FLAGS_${sanitizer}
      CMAKE_MODULE_LINKER_FLAGS_${sanitizer}
      CMAKE_SHARED_LINKER_FLAGS_${sanitizer}
      CMAKE_STATIC_LINKER_FLAGS_${sanitizer} )

    set( CMAKE_BUILD_TYPE_DROP_LIST ${CMAKE_BUILD_TYPE_DROP_LIST} "${sanitizer}")
  endif()
endforeach()

