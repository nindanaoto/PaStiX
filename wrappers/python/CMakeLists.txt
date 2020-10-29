###
#
#  @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.3
#  @author Mathieu Faverge
#  @date 2019-12-04
#
###

# Configure enum.py
if (PASTIX_INT64)
  set(PASTIX_PYTHON_INTEGER c_int64)
else()
  set(PASTIX_PYTHON_INTEGER c_int)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/enum.py.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/enum.py" @ONLY)

# Install python package
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/__pastix__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/pastix.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/enum.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/pypastix/solver.py
  DESTINATION lib/python/pypastix )

# Install python examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/simple.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/simple_obj.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/schur.py
  ${CMAKE_CURRENT_SOURCE_DIR}/examples/schur_obj.py
  DESTINATION examples
  )

## CTest execution
find_package(PythonInterp QUIET)
if (PYTHONINTERP_FOUND)
  set( PYTHON_TESTS
    simple step-by-step schur simple_obj schur_obj )

  # Add path to SPM, and to the .so in the environment to run the tests
  if( NOT SPM_FOUND )
    set( TEST_ENV_LIST
      "PYTHONPATH=$ENV{PYTHONPATH}:${CMAKE_SOURCE_DIR}/spm/wrappers/python"
      "LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm"
      "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm"
      "DYLDLIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/spm"
      )
  endif()
  foreach(example ${PYTHON_TESTS} )
    add_test( shm_python_${example} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/examples/${example}.py)
    set_tests_properties( shm_python_${example} PROPERTIES
      ENVIRONMENT "${TEST_ENV_LIST}" )
  endforeach()
endif()

