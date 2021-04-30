program olapmpi

  use olapparams
  use benchclock
  use olapcomms
  use olapcalc

  implicit none

  integer, parameter :: nhalotest = 8
  integer :: ihalotest
  logical :: halocalled
  
  !
  ! Set main paramemers: size of each halo and number of repetitions
  !

  integer, parameter :: nbuf = 10000

  integer :: irep, commrep, calcrep
  
  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, recvbuf

  integer :: rank, size, comm, cartcomm, dblesize

  integer, dimension(ndir, ndim) :: neighbour
  integer, dimension(ndim) :: coords, dims

  integer, parameter :: mib = 1024*1024

  logical :: finished, reorder = .false.

  logical, dimension(:,:),   allocatable :: flag
  logical, dimension(:,:,:), allocatable :: allflag

  integer :: idir, idim

  ! Periodic boundaries

  logical, dimension(ndim) :: periods = [.true., .true., .true.]

  double precision :: t0, t1, time, iorate, mibdata, calcval
  double precision :: bwidth, latency, timetarget
  double precision :: msgtime, dummycalctime, commtime, calctime

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

  allocate(flag(ndir, ndim))
  allocate(allflag(ndir, ndim, size))

  ! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*nbuf*ndir*ndim)/float(mib)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple MPI overlap benchmark'
     write(*,*) '----------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', dims(1), ', ', dims(2), ', ', dims(3), ')'
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ' microseconds'
  end if

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

  ! Get all the data

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)
  
  call inithalodata(rank, sendbuf, recvbuf, nbuf)

  ! estimate bwidth (GiB/s) and latency (us)

  bwidth  = 10.0
  latency = 1.0

  ! estimate time for one call of calculation which is 100 sqrts (us)

  dummycalctime = 10.0

  ! Set target time

  timetarget = 1.0

  ! estimate commrep

  msgtime = ndir*ndim*(latency*1.0d-6 + 1.0d-9*dble(nbuf*dblesize)/bwidth)
  commrep = timetarget/msgtime
  

  if (rank == 0) write(*,*) "msgtime = ", msgtime, ", commrep = ", commrep

  finished = .false.

  do while (.not. finished .and. commrep > 1)

     finished = .true.

     call MPI_Barrier(comm, ierr)
     t0 = benchtime()

     do irep = 1, commrep

        call haloirecvisendwaitall(sendbuf, recvbuf, nbuf, neighbour, cartcomm)

     end do

     call MPI_Barrier(comm, ierr)
     t1 = benchtime()

     commtime = t1 - t0

     ! Make sure all ranks have the same time so bcast from 0

     call MPI_Bcast(commtime, 1, MPI_DOUBLE_PRECISION, 0, cartcomm, ierr)

     if (rank == 0) then
        write(*,*) "commrep = ", commrep, ", secs = ", commtime
     end if

     if (commtime > 2.0*timetarget) then
        commrep = commrep/2
        finished = .false.
     end if

     if (commtime < 0.5*timetarget) then
        commrep = commrep*2
        finished = .false.
     end if

  end do
  
  ! estimate calcrep

  calcrep = 1.0d6*commtime/dummycalctime

  calcval = 1.0

  finished = .false.

  do while (.not. finished .and. calcrep > 1)

     finished = .true.

     t0 = benchtime()

     call dummycalc(calcval, calcrep)

     t1 = benchtime()

     calctime = t1 - t0

     ! Make sure all ranks have the same time so bcast from 0

     call MPI_Bcast(calctime, 1, MPI_DOUBLE_PRECISION, 0, cartcomm, ierr)

     if (rank == 0) then
        write(*,*) "calcrep = ", calcrep, ", secs = ", calctime
     end if

     if (calctime > 1.0*commtime) then
        calcrep = calcrep/2
        finished = .false.
     end if

     if (calctime < 0.4*commtime) then
        calcrep = calcrep*2
        finished = .false.
     end if

  end do
  
  iorate = dble(commrep)*mibdata/time

  if (rank == 0) then
     write(*,*) "Final commrep = ", commrep, ", commtime = ", commtime
     write(*,*) "Final calcrep = ", calcrep, ", calctime = ", calctime
  end if

  ! Normalise to a single comms call

  calcrep = calcrep / commrep

  if (rank == 0) then
     write(*,*) "Normalised calcrep = ", calcrep
  end if

  call checkrecvdata(flag, recvbuf, nbuf, cartcomm)
        
  call MPI_Gather(flag,    ndir*ndim, MPI_LOGICAL, &
                  allflag, ndir*ndim, MPI_LOGICAL, &
                  0, comm, ierr)

  if (rank == 0) then
           
     if (any(allflag(:,:,:) .eqv. .false.)) then

        write(*,*)
        write(*,*) "ERROR: halo data did not verify"
        write(*,*)
              
        do idir = 1, ndir
           do idim = 1, ndim

              if (any(allflag(idir, idim,:) .eqv. .false.)) then
                 write(*,*) count(.not. allflag(idir, idim, :)), &
                      " processes failed for idir, idim = ", &
                      idir, ", ", idim
              end if

           end do
        end do
     end if

     write(*,*)
           
  end if
           
  if (rank == 0) then
     write(*,*)
     write(*,*) "Finished"
     write(*,*) "--------"
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program olapmpi
