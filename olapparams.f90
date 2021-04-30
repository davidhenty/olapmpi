module olapparams

  use mpi

  integer, parameter :: ndim = 3
  integer, parameter :: ndir = 2    ! down and up

  integer, parameter :: dndir = 1   ! ordering needs to be like this
  integer, parameter :: updir = 2   ! for neighourhood collectives

  integer :: ierr

end module olapparams
