!
! @file fmultidof.f90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.3
! @author Mathieu Faverge
! @date 2019-11-12
!
program fmultidof
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf
  ! use mpi_f08
  implicit none

  type(pastix_data_t),        pointer                      :: pastix_data
  type(pastix_order_t),       pointer                      :: order => null()
  type(spmatrix_t),           pointer                      :: spm
  type(spmatrix_t),           pointer                      :: spm2
  integer(kind=pastix_int_t), target                       :: iparm(iparm_size)
  real(kind=c_double),        target                       :: dparm(dparm_size)
  real(kind=c_double)                                      :: normA
  integer(c_int)                                           :: info
  integer(kind=pastix_int_t)                               :: nrhs
  real(kind=c_double), dimension(:,:), allocatable, target :: x0, x, b
  type(c_ptr)                                              :: x0_ptr, x_ptr, b_ptr

  !
  ! Initialize the problem
  !   1- The matrix
  ! See genLaplacian function for the parameters
  ! Here we generate a 10x101x0 laplacian with 5 degrees of freedom per unknown
  !
  allocate( spm )
  call spmReadDriver( SpmDriverLaplacian, "d:10:10:10:4.:.5:5", spm, info )

  !
  ! We do not check the matrix as it would expand the spm matrix, and
  ! the compressed version would not be used. If you reuse this
  ! testing, you have to make sure that the compressed form of the
  ! matrix is correct without calling the spmCheckAndCorrect. One way
  ! to do it, is to call the function on the matrix without the values
  ! to make sure its correct.
  !
  ! allocate( spm2 )
  ! call spmCheckAndCorrect( spm, spm2, info )
  ! if ( info .ne. 0 ) then
  !    call spmExit( spm )
  !    spm = spm2
  ! end if
  ! deallocate( spm2 )

  ! Print information about the matrix
  call spmPrintInfo( spm )

  ! Scale A for better stability with low-rank computations
  call spmNorm( SpmFrobeniusNorm, spm, normA )
  call spmScalMatrix( 1. / normA, spm )

  !
  ! Solve the problem
  !

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, 0, iparm, dparm )

  ! 2- Analyze the problem
  !   a- Split the analyze in substeps. The first three steps handle
  !      multi-dof spm, the last one does not yet
  call pastix_subtask_order( pastix_data, spm, order, info );
  call pastix_subtask_symbfact( pastix_data, info );
  call pastix_subtask_reordering( pastix_data, info );

  !   b- Expand the matrix and the associated substructure
  call pastixExpand( pastix_data, spm )

  !   c- Finish the analyze step
  call pastix_subtask_blend( pastix_data, info )

  ! 3- Factorize the matrix
  call pastix_task_numfact( pastix_data, spm, info )

  !
  ! We need to generate the right hand side once the spm has been
  ! expanded, spmGenRHS is not yet compatible with multi-dof
  !
  nrhs = 10
  allocate(x0(spm%n, nrhs))
  allocate(x( spm%n, nrhs))
  allocate(b( spm%n, nrhs))
  x0_ptr = c_loc(x0)
  x_ptr  = c_loc(x)
  b_ptr  = c_loc(b)

  call spmGenRHS( SpmRhsRndX, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, info )
  x = b

  ! 4- Solve the problem
  call pastix_task_solve( pastix_data, nrhs, x_ptr, spm%n, info )

  ! 5- Refine the solution
  call pastix_task_refine( pastix_data, spm%n, nrhs, b_ptr, spm%n, x_ptr, spm%n, info )

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  !
  ! Check the solution
  !
  call spmCheckAxb( dparm(DPARM_EPSILON_REFINEMENT), nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, x_ptr, spm%n, info )

  call spmExit( spm )
  deallocate(spm)
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program fmultidof
