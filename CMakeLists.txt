###
#
#  @copyright 2013-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.3
#  @author Mathieu Faverge
#  @date 2020-02-05
#
###
cmake_minimum_required (VERSION 2.8)
project (PASTIX C CXX Fortran)

include(CMakeDependentOption)
include(CheckFunctionExists)

# The current version number
set(PASTIX_VERSION_MAJOR 6)
set(PASTIX_VERSION_MINOR 1)
set(PASTIX_VERSION_MICRO 0)

set(PASTIX_VERSION "${PASTIX_VERSION_MAJOR}.${PASTIX_VERSION_MINOR}.${PASTIX_VERSION_MICRO}")

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
include(AddSourceFiles)

if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules)
  set( MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules )
  list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH})
  include(MorseInit)
else()
  message(FATAL_ERROR "Submodule cmake_morse not initialized - run `git submodule update --init`")
endif()

# Define precision supported by PaStiX
# ------------------------------------
set( RP_PASTIX_DICTIONNARY
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/local_subs.py )
set( RP_PASTIX_PRECISIONS  "p;s;d;c;z" )
include(RulesPrecisions)

### System parameter detection
include(CheckSystem)

# Set the RPATH config
# --------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

### Misc options
option(BUILD_SHARED_LIBS
  "Build shared libraries" OFF)
option(BUILD_64bits
  "Build 64 bits mode" ON)

### Distributed engine parameters
option(PASTIX_WITH_MPI
    "Build PaStiX for distributed memory with MPI" OFF)

### GPU engine parameters
option(PASTIX_WITH_CUDA
    "Enable GPU support using CUDA kernels" OFF)

cmake_dependent_option(PASTIX_CUDA_FERMI
    "Enable the Fermi kernel instead of the batched kernel" OFF "PASTIX_WITH_CUDA" OFF)

### Tracing parameters
option(PASTIX_WITH_EZTRACE
  "Enable Eztrace for tracing CPU and GPU kernels" OFF)

### Runtime engines
option(PASTIX_WITH_STARPU
    "Build PaStiX with StarPU runtime support" OFF)
option(PASTIX_WITH_PARSEC
  "Build PaStiX with PaRSEC runtime support" OFF)

# Internal options
option(PASTIX_GENERATE_MODEL
  "Enable performances profiling for model generation" OFF)

# Internal options
option(PASTIX_WITH_FORTRAN
  "Enable Fortran files/interface/examples to be compiled" ON)

cmake_dependent_option(PASTIX_DISTRIBUTED
 "Enable the distributed interface (requires PASTIX_WITH_MPI)" OFF "PASTIX_WITH_MPI" OFF)

cmake_dependent_option(PASTIX_THREAD_COMM
 "Enable the specific communication thread (requires PASTIX_WITH_MPI)" ON "PASTIX_WITH_MPI" OFF)

option(PASTIX_INT64
  "Choose between int32 and int64 for integer representation" ON)

# Precisions generated
if(NOT PASTIX_PRECISIONS)
  set(PASTIX_PRECISIONS "s;d;c;z" CACHE STRING "The precisions to compile in PaSTiX (accepts a colon separated list of s;d;c;z)" FORCE)
else()
  set(PASTIX_PRECISIONS "${PASTIX_PRECISIONS}" CACHE STRING "The precisions to compile in PaSTiX (accepts a colon separated list of s;d;c;z)" FORCE)
endif()

# Ordering step options
option(PASTIX_ORDERING_SCOTCH
  "Enable Scotch Ordering" ON)
option(PASTIX_ORDERING_METIS
  "Enable Metis ordering"  OFF)
cmake_dependent_option(PASTIX_ORDERING_PTSCOTCH
  "Activate the PT-scotch ordering (requires PASTIX_DISTRIBUTED and PASTIX_ORDERING_SCOTCH)" ON
  "PASTIX_DISTRIBUTED;PASTIX_ORDERING_SCOTCH" OFF)

# Symbolic factorization options
option(PASTIX_SYMBOL_DUMP_SYMBMTX
  "Dump the generated symbol matrix in a postscript file" OFF)
option(PASTIX_ORDER_DRAW_LASTSEP
  "Dump the last separator in a ivview file" OFF)

# Analyze step options
option(PASTIX_BLEND_DEEPEST_DISTRIB
  "Enable the candidate distribution to go to the deepest node matching the criterion" ON)
option(PASTIX_BLEND_PROPMAP_2STEPS
  "Enable the two stages proportionnal mapping" OFF)

# Numerical factorization options
option(PASTIX_NUMFACT_DUMP_SOLVER
  "Dump the generated factorized matrix in a postscript file with information on block's compressibility" OFF)
option(PASTIX_SUPERNODE_STATS
  "Enable statistics per original supernodes" OFF)

# Options to check
option(FORGET_PARTITION
  "Force to forget the partition generated by Scotch" OFF)
option(COMPACT_SMX
  "Optimization for solve computations (TODO: check if not obsolete because results don't converge)" OFF)

# Options to increase verbosity while debuging
option(PASTIX_DEBUG_GRAPH
  "Debug the graph data structure" OFF)
option(PASTIX_DEBUG_ORDERING
  "Debug the ordering step" OFF)
option(PASTIX_DEBUG_SYMBOL
  "Debug the symbol matrix factorization" OFF)
option(PASTIX_DEBUG_BLEND
  "Debug the analyze step" OFF)
option(PASTIX_DEBUG_DUMP_COEFTAB
  "Dump all cblk on disk after modifications" OFF)
cmake_dependent_option(PASTIX_DEBUG_FACTO
  "Debug the factorization step" OFF
  "PASTIX_WITH_MPI" OFF)
option(PASTIX_DEBUG_SOLVE
  "Debug the solve step" OFF)
cmake_dependent_option(PASTIX_DEBUG_PARSEC
  "Debug the parsec implementation of the factorization step" OFF
  "PASTIX_WITH_PARSEC" OFF)
option(PASTIX_DEBUG_LR
  "Debug the lowrank kernels (Valgrind)" OFF)
option(PASTIX_DEBUG_LR_NANCHECK
  "Debug the lowrank kernels (nancheck)" OFF)
option(PASTIX_DEBUG_GMRES
  "Debug the gmres refinment by enabling the computation of the true resiudal at each iteration" OFF)
option(PASTIX_DEBUG_EXIT_ON_SIGSEGV
  "Generate a segfault on exit calls to trakc in gdb" OFF)
cmake_dependent_option(PASTIX_DEBUG_MPI
  "Dump all the MPI messages" OFF
  "PASTIX_WITH_MPI" OFF)

###############################################################################
# Look for dependencies #
#########################

if (PASTIX_WITH_FORTRAN)
  include(FortranCInterface)
  FortranCInterface_HEADER(common/FCmangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE "FC_")
  link_directories( ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} )
endif()

# PaStiX depends on LAPACKE and CBLAS
#------------------------------------
find_package(CBLAS) # Should be REQUIRED for BLAS sequential only
if(CBLAS_FOUND)
    message(STATUS "cblas: ${CBLAS_INCLUDE_DIRS}")
    include_directories(${CBLAS_INCLUDE_DIRS})
endif()

find_package(LAPACKE) # Should be also REQUIRED
if(LAPACKE_FOUND)
  message(STATUS "lapacke: ${LAPACKE_INCLUDE_DIRS}")
  include_directories(${LAPACKE_INCLUDE_DIRS})
endif()

### Store dependencies not handled with pkg-config
set( DEPS_LIBRARIES
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP}
  )

list(APPEND CMAKE_INSTALL_RPATH
  ${LAPACKE_LIBRARY_DIRS_DEP}
  ${CBLAS_LIBRARY_DIRS_DEP}
  )

# PaStiX depends on HwLoc
#---------------------------
find_package(HWLOC)
set(HAVE_HWLOC ${HWLOC_FOUND})
if( HWLOC_FOUND )
  link_directories( ${HWLOC_LIBRARY_DIRS} )
  link_libraries( ${HWLOC_LIBRARIES} )
  include_directories( ${HWLOC_INCLUDE_DIRS} )
  list(APPEND CMAKE_INSTALL_RPATH
    ${HWLOC_LIBRARY_DIRS}
    )
endif (HWLOC_FOUND)

# PaStiX might depend on MPI
#------------------------------
if (PASTIX_WITH_MPI)
  # Force the detection of the C library
  find_package(MPI)
  string(REPLACE ";" " " MPI_C_COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
  if (MPI_C_FOUND)
    list(APPEND EXTRA_LIBS ${MPI_C_LIBRARIES} )
    include_directories( ${MPI_C_INCLUDE_PATH} )

    # Check to see if support for MPI 2.0 is available
    set(saved_include "${CMAKE_REQUIRED_INCLUDES}" )
    set(saved_libs    "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${MPI_C_INCLUDE_PATH}")
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${MPI_C_LIBRARIES}")
    check_function_exists(MPI_Type_create_resized HAVE_MPI_20)
    set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
    set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")
  else (MPI_C_FOUND)
    message(FATAL_ERROR "MPI is required but was not found. Please provide an MPI compiler in your environment or configure with -DPASTIX_WITH_MPI=OFF")
  endif (MPI_C_FOUND)
endif (PASTIX_WITH_MPI)

# PaStiX might depend on CUDA
#-------------------------------
if( PASTIX_WITH_CUDA )

    find_package(CUDA)

    set(HAVE_CUDA ${CUDA_FOUND})
    if (CUDA_FOUND)
        if(CUDA_VERSION VERSION_LESS "3.0")
            set(CUDA_HOST_COMPILATION_CPP OFF)
        endif(CUDA_VERSION VERSION_LESS "3.0")
        set(CUDA_BUILD_EMULATION OFF)
        include_directories(${CUDA_INCLUDE_DIRS})
        link_libraries(${CUDA_LIBRARIES}) # need CUDA libs to link.
        set(saved_include "${CMAKE_REQUIRED_INCLUDES}")
        set(saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
        set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${CUDA_INCLUDE_DIRS}")
        set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CUDA_LIBRARIES}")
        if(CUDA_VERSION VERSION_LESS "4.0")
            set(CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS 0)
        else()
            check_function_exists(cuDeviceCanAccessPeer CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS)
        endif()
        set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
        set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DCUDA_SM_VERSION=20")
    else(CUDA_FOUND)
        message(FATAL_ERROR "CUDA is required but was not found. Please provide a CUDA library in your environment or configure with -DPASTIX_WITH_CUDA=OFF")
    endif (CUDA_FOUND)

endif()

# PaStiX might depend on StarPU
#---------------------------------
if (PASTIX_WITH_STARPU)

  set(PASTIX_STARPU_VERSION "1.3" CACHE STRING "oldest STARPU version desired")

  # Create list of components in order to make a single call to find_package(starpu...)
  set(STARPU_COMPONENT_LIST "HWLOC")
  if(PASTIX_WITH_MPI)
    list(APPEND STARPU_COMPONENT_LIST "MPI")
  endif()
  if(PASTIX_WITH_CUDA)
    list(APPEND STARPU_COMPONENT_LIST "CUDA")
  endif()
  if(PASTIX_WITH_FXT)
    list(APPEND STARPU_COMPONENT_LIST "FXT")
  endif()

  find_package(STARPU ${PASTIX_STARPU_VERSION} REQUIRED
    COMPONENTS ${STARPU_COMPONENT_LIST})

  # Add definition and include_dir if found
  if(STARPU_FOUND)
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    link_directories(${STARPU_LIBRARY_DIRS_DEP})
    set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
    if(PASTIX_USE_FXT)
      # check if fxt profiling is accessible in starpu and activate it in magmamorse
      set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
      unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
      check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
      if ( STARPU_FXT_START_PROFILING_FOUND )
        message(STATUS "Set HAVE_STARPU_FXT_PROFILING - Activate FxT profiling through StarPU")
        add_definitions(-DHAVE_STARPU_FXT_PROFILING)
      endif()
    endif()
    if (PASTIX_WITH_MPI)
      # Check if a specific function exist
      set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
      unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
      check_function_exists(starpu_mpi_data_register STARPU_MPI_DATA_REGISTER_FOUND)
      if ( STARPU_MPI_DATA_REGISTER_FOUND )
        add_definitions(-DHAVE_STARPU_MPI_DATA_REGISTER)
      endif()
    endif()

    list(APPEND CMAKE_INSTALL_RPATH
      ${STARPU_LIBRARY_DIRS_DEP}
      )

  else(STARPU_FOUND)
    message(FATAL_ERROR "StarPU with ${STARPU_COMPONENT_LIST} components is required but was not found. Please provide the StarPU library in your environment or configure with -DPASTIX_WITH_STARPU=OFF")
  endif()

endif (PASTIX_WITH_STARPU)

# PaStiX might depend on Parsec
#------------------------------
if (PASTIX_WITH_PARSEC)

    # Create list of components in order to make a single call to find_package(starpu...)
    set(PARSEC_COMPONENT_LIST "HWLOC")
    if(PASTIX_WITH_MPI)
        list(APPEND PARSEC_COMPONENT_LIST "MPI")
    endif()
    if(PASTIX_WITH_CUDA)
        list(APPEND PARSEC_COMPONENT_LIST "CUDA")
    endif()

    find_package(PARSEC REQUIRED
      COMPONENTS ${PARSEC_COMPONENT_LIST})

    # Add definition and include_dir if found
    if(PARSEC_FOUND)
      include_directories(${PARSEC_INCLUDE_DIRS_DEP})
      link_directories(${PARSEC_LIBRARY_DIRS_DEP})
      list(APPEND CMAKE_INSTALL_RPATH ${PARSEC_LIBRARY_DIRS_DEP})
    else()
        message(FATAL_ERROR "Parsec with ${PARSEC_COMPONENT_LIST} components is required but was not found. Please provide the PaRSEC library in your environment or configure with -DPASTIX_WITH_PARSEC=OFF")
    endif()

    # Add the index-array dep management option to PTGCC compilation flags.
    include(RulesJDF)
    set(PARSEC_PTGPP_CFLAGS "--dep-management;index-array;${PARSEC_PTGPP_CFLAGS}")

endif (PASTIX_WITH_PARSEC)

### Add trace generation
if(PASTIX_WITH_EZTRACE)
  if (NOT EZTRACE_FOUND)
    find_package(EZTRACE REQUIRED)
  endif()

  include_directories(${EZTRACE_INCLUDE_DIRS})
  link_directories(${EZTRACE_LIBRARY_DIRS})
endif(PASTIX_WITH_EZTRACE)

# PaStiX might depend on Scotch/PT-Scotch
#----------------------------------------
if (PASTIX_ORDERING_PTSCOTCH)
  find_package(PTSCOTCH)
  if (PTSCOTCH_FOUND)
    include_directories( ${PTSCOTCH_INCLUDE_DIRS_DEP} )
    link_directories( ${PTSCOTCH_LIBRARY_DIRS_DEP} )
    link_libraries( ${PTSCOTCH_LIBRARIES_DEP} )

    list(APPEND CMAKE_INSTALL_RPATH
      ${PTSCOTCH_LIBRARY_DIRS_DEP}
      )
  else (PTSCOTCH_FOUND)
    message(FATAL_ERROR "PTSCOTCH is required if PASTIX_ORDERING_PTSCOCTH is on but was not found. Please provide an PtScotch library in your environment or configure with -DPASTIX_ORDERING_PTSCOTCH=OFF")
  endif (PTSCOTCH_FOUND)

  # Check coherency for integer size
  if(PASTIX_INT64 AND NOT PTSCOTCH_Num_8)
    message(FATAL_ERROR "PASTIX_INT64 is enabled and provided PT-Scotch is not compiled with int64 support.")
  endif()
  if(NOT PASTIX_INT64 AND NOT PTSCOTCH_Num_4)
    message(FATAL_ERROR "PASTIX_INT64 is disabled and provided PT-Scotch is not compiled with int32 support.")
  endif()
endif()

if (PASTIX_ORDERING_SCOTCH)
  find_package(SCOTCH)
  if (SCOTCH_FOUND)
    message(STATUS "Scotch inlude dirs: ${SCOTCH_INCLUDE_DIRS}" )
    include_directories( ${SCOTCH_INCLUDE_DIRS} )
    link_directories( ${SCOTCH_LIBRARY_DIRS} )
    link_libraries( ${SCOTCH_LIBRARIES} )

    list(APPEND CMAKE_INSTALL_RPATH
      ${SCOTCH_LIBRARY_DIRS}
      )
  else()
    message(FATAL_ERROR "Scotch is required but was not found. Please provide a Scotch library in your environment or configure with -DPASTIX_ORDERING_SCOTCH=OFF")
  endif()

  # Check coherency for integer size
  if(PASTIX_INT64 AND NOT SCOTCH_Num_8)
    message(FATAL_ERROR "PASTIX_INT64 is enabled and provided Scotch is not compiled with int64 support.")
  endif()
  if(NOT PASTIX_INT64 AND NOT SCOTCH_Num_4)
    message(FATAL_ERROR "PASTIX_INT64 is disabled and provided Scotch is not compiled with int32 support.")
  endif()
endif()

set( DEPS_LIBRARIES
  ${DEPS_LIBRARIES}
  ${PTSCOTCH_LIBRARIES_DEP}
  ${SCOTCH_LIBRARIES}
  )

# PaStiX might depend on Metis/ParMetis
#--------------------------------------
if (PASTIX_ORDERING_METIS)
  find_package(METIS)
  if (METIS_FOUND)
    include_directories( ${METIS_INCLUDE_DIRS} )
    link_directories( ${METIS_LIBRARY_DIRS} )
    link_libraries( ${METIS_LIBRARIES} )

    list(APPEND CMAKE_INSTALL_RPATH
      ${METIS_LIBRARY_DIRS}
      )
  else()
    message(FATAL_ERROR "Metis is required but was not found. Please provide a Metis library in your environment or configure with -DPASTIX_ORDERING_METIS=OFF")
  endif()

  # Check coherency for integer size
  if(PASTIX_INT64 AND NOT METIS_Idx_8)
    message(FATAL_ERROR "PASTIX_INT64 is enabled and provided Metis is not compiled with int64 support.")
  endif()
  if(NOT PASTIX_INT64 AND NOT METIS_Idx_4)
    message(FATAL_ERROR "PASTIX_INT64 is disabled and provided Metis is not compiled with int32 support.")
  endif()
endif()

set( DEPS_LIBRARIES
  ${DEPS_LIBRARIES}
  ${METIS_LIBRARIES}
  )

# PaStiX might depend on GTG
#---------------------------
find_package(GTG QUIET)

if (GTG_FOUND)
  # Symbolic factorization options
  option(PASTIX_BLEND_GENTRACE
    "Allow trace generation in Blend simulation" OFF)

  if (PASTIX_BLEND_GENTRACE)
    include_directories( ${GTG_INCLUDE_DIRS} )
    link_directories( ${GTG_LIBRARY_DIRS} )
  endif()
else()
  set( PASTIX_BLEND_GENTRACE OFF )
endif()

# valgrind detection
#-------------------
find_program( MEMORYCHECK_COMMAND valgrind )
set( MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )

#
##
###
# Finished detecting the system, lets do our own things now
###
##
#

# Compile the spm library
if ( PASTIX_WITH_FORTRAN )
  find_package( SPM QUIET COMPONENTS FORTRAN )
else()
  find_package( SPM QUIET )
endif()

if ( SPM_FOUND )
  if( (PASTIX_INT64 AND (NOT SPM_INTSIZE EQUAL 8)) OR (PASTIX_INT32 AND (NOT SPM_INTSIZE EQUAL 4)) )
    message( WARNING "PASTIX_INT64 is enabled and provided SPM does not have the same support. Internal SPM will be compiled.")
    unset( SPM_FOUND )
  endif()
endif()

if ( SPM_FOUND )
  message( STATUS "Use installed SPM" )
  include_directories(${SPM_INCLUDE_DIRS})
else()
  message( STATUS "Use internal SPM" )
  add_subdirectory(spm)
  set( SPM_LIBRARIES spm )
  set( SPM_LIBRARIES_FORTRAN spmf )
  include_directories("spm/include")
  include_directories("${CMAKE_CURRENT_BINARY_DIR}/spm/include")
endif()

# Disable restrict (temporary)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Drestrict=")

set(PROJECT_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories(BEFORE "${PROJECT_INCLUDE_DIR}")

STRING(COMPARE EQUAL ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} PROJECT_COMPILE_INPLACE)
if(NOT PROJECT_COMPILE_INPLACE)
  include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")
endif(NOT PROJECT_COMPILE_INPLACE)

#Configuration header
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/pastix/config.h.in"
  "${PROJECT_INCLUDE_DIR}/pastix/config.h")
install(FILES "${PROJECT_INCLUDE_DIR}/pastix/config.h" DESTINATION include/pastix)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/common")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/common")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/kernels")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/kernels")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/bcsc")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/bcsc")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/refinement")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/refinement")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/sopalin")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/blend")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/graph")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/order")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/symbol")

# Sub modules
set(modules kernels refinement) #csc kernels sopalin common)
foreach (module ${modules})
   add_subdirectory(${module})
endforeach()

if (PASTIX_WITH_PARSEC)
    include_directories("${CMAKE_CURRENT_BINARY_DIR}/sopalin")
    add_subdirectory(sopalin/parsec)
endif()

if (PASTIX_WITH_STARPU)
    include_directories("${CMAKE_CURRENT_BINARY_DIR}/sopalin")
    add_subdirectory(sopalin/starpu)
endif()

include(RulesPrecisions)

##
# Headers from bcsc directory
##
set(generated_bcsc_headers "")
set(HEADERS
  bcsc/bcsc_z.h
)

precisions_rules_py(generated_bcsc_headers
  "${HEADERS}"
  TARGETDIR  "bcsc"
  PRECISIONS "p;s;d;c;z")

add_custom_target(bcsc_headers_tgt
  DEPENDS ${generated_bcsc_headers} )

##
# Headers from sopalin directory
##
set(generated_sopalin_headers "")
set(HEADERS
  sopalin/coeftab_z.h
)

precisions_rules_py(generated_sopalin_headers
  "${HEADERS}"
  TARGETDIR  "sopalin"
  PRECISIONS "p;s;d;c;z")

add_custom_target( sopalin_headers
  DEPENDS ${generated_sopalin_headers} )

###############
# Source files

##
# Source files from common directory
##
set(generated_common_sources "")
set(SOURCES
  common/z_integer.c
  )
precisions_rules_py(generated_common_sources
  "${SOURCES}"
  TARGETDIR  "common"
  PRECISIONS "p;s;d;c;z"
  )

##
# Source files from bcsc directory
##
set(generated_bcsc_sources "")
set(SOURCES
  bcsc/bcsc_zinit.c
  bcsc/bcsc_znorm.c
  bcsc/bcsc_zspmv.c
  bcsc/bvec_zcompute.c
  )
precisions_rules_py(generated_bcsc_sources
  "${SOURCES}"
  TARGETDIR  "bcsc"
  PRECISIONS "p;s;d;c;z"
  )

##
# Source files from sopalin directory
##
set(generated_sopalin_sources "")
set(SOURCES
  sopalin/coeftab_z.c
  sopalin/coeftab_zinit.c
  sopalin/sequential_zdiag.c
  sopalin/sequential_zgetrf.c
  sopalin/sequential_zhetrf.c
  sopalin/sequential_zpotrf.c
  sopalin/sequential_zpxtrf.c
  sopalin/sequential_zsytrf.c
  sopalin/sequential_ztrsm.c
  )
precisions_rules_py(generated_sopalin_sources
  "${SOURCES}"
  TARGETDIR  "sopalin"
  PRECISIONS "p;s;d;c;z"
  )

##
# Source files from refinement directory
##
set(generated_refinement_sources "")
set(SOURCES
  refinement/z_refine_bicgstab.c
  refinement/z_refine_functions.c
  refinement/z_refine_gmres.c
  refinement/z_refine_grad.c
  refinement/z_refine_pivot.c
  )
precisions_rules_py(generated_refinement_sources
  "${SOURCES}"
  TARGETDIR  "refinement"
  PRECISIONS "p;s;d;c;z"
  )

#
# Build the pastix library
#
set(PASTIX_LIB_SRCS
  # Files that are fully converted
  #
  graph/graph.c
  graph/graph_compute_projection.c
  graph/graph_connected_components.c
  graph/graph_io.c
  graph/graph_isolate.c
  graph/graph_prepare.c
  graph/graph_symmetrize.c
  #
  order/order.c
  order/order_add_isolate.c
  order/order_amalgamate.c
  order/order_apply_level_order.c
  order/order_check.c
  order/order_find_supernodes.c
  order/order_grids.c
  order/order_io.c
  order/pastix_subtask_order.c
  #
  symbol/fax_csr.c
  symbol/fax_csr_amalgamate.c
  symbol/fax_csr_direct.c
  symbol/fax_csr_iluk.c
  symbol/symbol.c
  symbol/symbol_base.c
  symbol/symbol_check.c
  symbol/symbol_cost.c
  symbol/symbol_cost_flops.c
  symbol/symbol_cost_perfs.c
  symbol/symbol_draw.c
  symbol/symbol_expand.c
  symbol/symbol_fax_direct.c
  symbol/symbol_fax_iluk.c
  symbol/symbol_io.c
  symbol/symbol_reorder.c
  symbol/symbol_reordering.c
  symbol/pastix_subtask_reordering.c
  symbol/pastix_subtask_symbfact.c

  # Files that still require some changes
  blend/pastix_subtask_blend.c
  blend/pastix_task_analyze.c

  blend/blendctrl.c
  blend/cand.c
  blend/cand_gendot.c
  blend/cost.c
  blend/elimintree.c
  blend/extendVector.c
  blend/extracblk.c
  blend/propmap.c
  blend/simu.c
  blend/simu_run.c
  blend/simu_task.c
  blend/solver.c
  blend/solver_backup.c
  blend/solver_check.c
  blend/solver_copy.c
  blend/solver_draw.c
  blend/solver_io.c
  blend/solver_matrix_gen.c
  blend/solver_matrix_gen_utils.c
  blend/splitsymbol.c

  common/api.c
  common/integer.c
  common/isched.c
  common/models.c
  common/get_options.c
  common/getline.c
  #
  bcsc/bcsc.c
  bcsc/bvec.c
  #
  sopalin/coeftab.c
  sopalin/schur.c
  sopalin/diag.c
  sopalin/pastix_task_sopalin.c
  sopalin/pastix_task_solve.c
  sopalin/pastix.c
  #
  refinement/pastix_task_refine.c
  #
  ${generated_common_sources}
  ${generated_bcsc_sources}
  ${generated_sopalin_sources}
  ${generated_refinement_sources}
  )

set_source_files_properties(
  sopalin/coeftab.h
  PROPERTIES DEPENDS sopalin_headers
  )

if(PASTIX_ORDERING_SCOTCH)
  set(PASTIX_LIB_SRCS
    ${PASTIX_LIB_SRCS}
    graph/graph_compute_kway.c
    order/order_compute_scotch.c
    order/order_draw.c
    order/order_supernodes.c
    )
endif()
if(PASTIX_ORDERING_PTSCOTCH)
  set(PASTIX_LIB_SRCS
    ${PASTIX_LIB_SRCS}
    order/order_compute_ptscotch.c
    )
endif()
if(PASTIX_ORDERING_METIS)
  list(APPEND PASTIX_LIB_SRCS
    order/order_compute_metis.c
    )
endif()
if(HWLOC_FOUND)
  list(APPEND PASTIX_LIB_SRCS
    common/isched_hwloc.c
    )
else()
  list(APPEND PASTIX_LIB_SRCS
    common/isched_nohwloc.c
    )
endif()

add_library(pastix
  ${PASTIX_LIB_SRCS}
  )

set_target_properties(pastix PROPERTIES VERSION ${PASTIX_VERSION} SOVERSION ${PASTIX_VERSION_MAJOR}.${PASTIX_VERSION_MINOR})

add_dependencies(pastix
  kernels_headers_tgt
  bcsc_headers_tgt
  refinement_headers_tgt
  sopalin_headers
)

target_link_libraries(pastix
  ${SPM_LIBRARIES}
  pastix_kernels
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP}
)

if(PASTIX_ORDERING_PTSCOTCH)
  target_link_libraries(pastix ${PTSCOTCH_LIBRARIES_DEP} )
endif()
if(PASTIX_ORDERING_SCOTCH)
  target_link_libraries(pastix ${SCOTCH_LIBRARIES} )
endif()
if(PASTIX_ORDERING_METIS)
  target_link_libraries(pastix ${METIS_LIBRARIES} )
endif()

if(PASTIX_WITH_PARSEC)
  target_link_libraries(pastix
    pastix_parsec
  )

  add_dependencies(pastix
    parsec_headers_tgt
    )
endif()

if(PASTIX_WITH_STARPU)
  target_link_libraries(pastix
    pastix_starpu
    )

  add_dependencies(pastix
    starpu_headers_tgt
    )
endif()

if (PASTIX_BLEND_GENTRACE)
  target_link_libraries(pastix
    ${GTG_LIBRARIES}
    )
endif()

# foreach (module2 ${modules})
#   add_dependencies(pastix ${module2}_headers)
# endforeach()

install(TARGETS pastix
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

# Install header files
set(PASTIX_HDRS
  include/pastix/api.h
  include/pastix/old_api.h
  include/pastix/datatypes.h
  include/pastix/order.h
  graph/graph.h
  include/cblas.h
  include/lapacke.h
)
install(FILES
  include/pastix.h
  DESTINATION include )
install(FILES ${PASTIX_HDRS} DESTINATION include/pastix )

## Executable and tests
enable_testing()
include(CTest)
# Examples executables
add_subdirectory(example)
# Testing executables
add_subdirectory(test)

### Wrappers
add_subdirectory(wrappers)

### Build pkg-config and environment file
include(GenPkgConfig)

set(PASTIX_PKGCONFIG_LIBS pastix pastix_kernels)

# StarPU
if(PASTIX_WITH_STARPU)
  list(APPEND PASTIX_PKGCONFIG_LIBS pastix_starpu)
  if ( STARPU_FOUND_WITH_PKGCONFIG )
    if ( PASTIX_WITH_MPI )
      list(APPEND PASTIX_PKGCONFIG_REQUIRED starpumpi-${PASTIX_STARPU_VERSION})
    else()
      list(APPEND PASTIX_PKGCONFIG_REQUIRED starpu-${PASTIX_STARPU_VERSION})
    endif()
  else()
    list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${STARPU_LIBRARIES_DEP})
  endif()
endif()

# PaRSEC
if(PASTIX_WITH_PARSEC)
  list(APPEND PASTIX_PKGCONFIG_LIBS pastix_parsec)
  if ( PARSEC_FOUND_WITH_PKGCONFIG )
    list(APPEND PASTIX_PKGCONFIG_REQUIRED parsec)
  else()
    list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${PARSEC_LIBRARIES_DEP})
  endif()
endif()

if(PASTIX_WITH_CUDA)
  list(APPEND PASTIX_PKGCONFIG_LIBS pastix_kernels_cuda)
endif()

if(PASTIX_WITH_EZTRACE)
  list(APPEND PASTIX_PKGCONFIG_REQUIRED eztrace litl)
endif()

list(APPEND PASTIX_PKGCONFIG_INCS
  )
list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP}
  )
list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE
  ${EXTRA_LIBRARIES}
  )

# HwLoc
if ( HWLOC_FOUND_WITH_PKGCONFIG )
  list(APPEND PASTIX_PKGCONFIG_REQUIRED hwloc)
else()
  list(APPEND PASTIX_PKGCONFIG_LIBS_PRIVATE ${HWLOC_LIBRARIES})
endif()

generate_pkgconfig_files(
  "${CMAKE_CURRENT_SOURCE_DIR}/pastix.pc.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/pastixf.pc.in"
  PROJECTNAME PASTIX )

generate_env_file( PROJECTNAME PASTIX )

#############################################################
#
#      Build documentation
#
#############################################################
add_documented_files(
  README.md
  docs/doxygen/chapters/eztrace.doxy
  docs/doxygen/chapters/thread_bindings.doxy
  tools/binding_for_multimpi.c
  kernels/pastix_lowrank.h
  #
  include/pastix.h
  include/pastix/api.h
  common/api.c
  common/get_options.c
  common/models.c
  common/pastixdata.h
  #
  graph/graph.h
  graph/graph.c
  graph/graph_compute_kway.c
  graph/graph_compute_projection.c
  graph/graph_io.c
  graph/graph_isolate.c
  graph/graph_prepare.c
  graph/graph_symmetrize.c
  #
  include/pastix/order.h
  #order/order_scotch_strats.h TODO: include in developer version
  order/order.c
  order/order_add_isolate.c
  order/order_amalgamate.c
  order/order_apply_level_order.c
  order/order_check.c
  order/order_draw.c
  order/order_find_supernodes.c
  order/order_grids.c
  order/order_io.c
  order/order_compute_metis.c
  order/order_compute_ptscotch.c
  order/order_compute_scotch.c
  order/order_supernodes.c
  order/pastix_subtask_order.c
  #
  symbol/fax_csr.h
  symbol/symbol.h
  symbol/symbol_cost.h  #TODO: include in developer version
  symbol/symbol_fax.h   #TODO: include in developer version
  symbol/fax_csr.c
  symbol/fax_csr_amalgamate.c
  symbol/fax_csr_direct.c
  symbol/fax_csr_iluk.c
  symbol/symbol.c
  symbol/symbol_base.c
  symbol/symbol_check.c
  symbol/symbol_cost.c
  symbol/symbol_cost_flops.c
  symbol/symbol_cost_perfs.c
  symbol/symbol_draw.c
  symbol/symbol_expand.c
  #symbol/symbol_fax.c Not documented: template file
  symbol/symbol_fax_direct.c
  symbol/symbol_fax_iluk.c
  symbol/symbol_io.c
  symbol/symbol_reorder.c
  symbol/symbol_reordering.c
  symbol/pastix_subtask_reordering.c
  symbol/pastix_subtask_symbfact.c
  #
  blend/blend.h
  blend/blendctrl.h
  blend/blendctrl.c
  blend/cand.h
  blend/cand.c
  blend/cand_gendot.c
  blend/cost.h
  blend/cost.c
  blend/elimintree.h
  blend/elimintree.c
  blend/extendVector.h
  blend/extendVector.c
  blend/extracblk.h
  blend/extracblk.c
  blend/perf.h
  #
  blend/simu.h
  blend/simu_timer.h
  blend/simu.c
  blend/simu_run.c
  blend/simu_task.c
  #
  blend/solver.h
  blend/solver.c
  blend/solver_backup.c
  blend/solver_check.c
  blend/solver_copy.c
  blend/solver_draw.c
  blend/solver_io.c
  blend/solver_matrix_gen.c
  blend/solver_matrix_gen_utils.h
  blend/solver_matrix_gen_utils.c
  #
  blend/propmap.c
  blend/splitsymbol.c
  blend/pastix_subtask_blend.c
  blend/pastix_task_analyze.c
  #
  bcsc/bcsc.c
  bcsc/bcsc.h
  bcsc/bvec.c
  bcsc/bvec.h
  #
  sopalin/coeftab.h
  sopalin/coeftab.c
  sopalin/schur.c
  sopalin/pastix.c
  sopalin/pastix_task_sopalin.c
  sopalin/pastix_task_solve.c
  #
  refinement/pastix_task_refine.c
  #
  example/analyze.c
  example/compress.c
  example/personal.c
  example/reentrant.c
  example/refinement.c
  example/schur.c
  example/simple.c
  example/step-by-step.c
  )

add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_common_sources}
  ${generated_bcsc_headers}
  ${generated_bcsc_sources}
  ${generated_sopalin_headers}
  ${generated_sopalin_sources}
  ${generated_refinement_sources}
)

#
# Build documentation
#
add_subdirectory(docs)

#-- Add a custom target to run Doxygen when ever the project is built
add_custom_target (tags
  COMMAND etags ${PASTIX_LIB_SRCS}
  DEPENDS ${PASTIX_LIB_SRCS} )

#############################################################
# Print Options
#############################################################
include(PrintOpts)
