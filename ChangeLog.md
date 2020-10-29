# pastix-6.1.0

- Add a new dynamic scheduler supported by the internal threads for numerical factorization
- Add MPI support for the numerical factorization and solve
        - Available for all schedulers: sequential, internal threads (static, dynamic), StarPU, and PaRSEC
        - WARNING: The RHS is not distributed yet, and must be replicated on all nodes
        - WARNING: The low-rank and Schur functionalities are not available in distributed yet
- Enable the use of an external SPM module
- Improve splitting strategy:
        - avoid unnecessary splits when using K-Way
        - reduce the range of possible split to limit the apparition of small blocks
- Change the preselected behavior to be:
        - never compressed in the JustInTime scenario
        - compress the preselected block just before applying the TRSM in the MinimalMemory scenario
        - the behavior can be change through IPARM_COMPRESS_PRESELECT
- Add a cmake summary
- Add coverity scan
- Update README.md and documentation

# pastix-6.0.3

- Update spm module to ada4963
- Update morse_cmake to ade4996
  - CMake: Update cmake_module to integrate the last version of the precision generator
- Change StarPU requirement to >= 1.3
- Refactor and extend the CI/CTests
- Update documentation
- Low-rank:
  - Add a new parameter IPARM_COMPRESS_RELTOL to switch between absolute and relative tolerance
  - Improved stability of low-rank kernels
  - Extend the number of tests
  - Add rotation QR kernels (unstable/work in progress)
  - Enable multiple low-rank factorization in a row
- Supports compilation with mpicc (no distributed solver yet)
- Add separated output directories for future distributed process, or for MPI multiple instances
- Octave: Fix issue with number of threads larger than the number of columns
- Octave: Fix compilation on Windows system
- Documentation: add documentation on process binding
- HwLoc: fix binding when already restricted through batch scheduler and/or MPI
- Fix issue solverstack/pastix#35, make pastix_task_analyze thread safe
- Add support for multi-dof in Fortran
- Fix issue in simulation, and a switch between cost and tree levels
- Refinement: Fix issue with gemv computation and PastixConjTrans
- CMake: Enable a round-robin selection of CMAKE_BUILD_TYPE  depending on the sanitizers provided by the compiler
- Homebrew: update formula

# pastix-6.0.2

- Integrate the clusting strategies developped for low-rank (See https://hal.inria.fr/hal-01961675)
- Restructure the ordering/symbolic factorization code to make sure with exit the ordering step with permutation, partition, and elimination tree.
- Relook the splitting/proportional mapping strategy
- Add new compression kernels: PQRCP, RQRCP, and TQRCP
- Fix inplace compilation (Issue #36)
- Fix issue when StarPU threads where fighting for ressources with PaStiX threads (Add pause/resume calls)
- Handle multi-dof static and/or variable in the analysis steps

# pastix-6.0.1

- Support for HWLOC 2.0.0
- Move the SPM library as a submodule
- Parallel (multithreaded) version of the reordering step
- Parallel (multithreaded) version of the refinement steps
- StarPU version of the solve step
- Fix Python/Fortran interface
- Fix Schur functions
- Fix multi-RHS solve
- Fix for METIS
- Update morse_cmake FindPACKAGE
- Update PaRSEC for release 6.0.1
- Improve PKG-CONFIG
- Improve documentation
- Better handle of disconnected graphs
- Add optimal reordering for grids
- Add more detailed statistics during analysis step
- Add detailed statistics about memory gain for the low-rank solver
- Add a function to compress the solver matrix outside the factorization step
- Add an example to dump the symbol matrix including the ranks of the block
- Add an refinement driver and testings
- Add a subtask_refine which does not perform the vector ordering
- Add a more complex testing based on example proposed by @andrea3.14
- Add an iparm IPARM_APPLYPERM_WS to enable/disable the use of an extra workspace to make the functions bvec_xlapmr thread-safe (by default, it is enabled, if disabled, the functions have no memory overhead but loose the thread-safe property)
- Remove the sparse-kit package to avoid conflict (the driver is replaced by HB)

# pastix-6.0.0

- low-rank compression (See https://hal.inria.fr/hal-01824275)
- static scheduler, PaRSEC and StarPU runtime support
- GPUs (Kepler) and KNL support through runtime systems
