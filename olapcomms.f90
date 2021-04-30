module olapcomms

  use olapparams
  use olapcalc
  
  implicit none

  integer, parameter :: nrequest = 2  ! receive and send

  integer, parameter :: rrequest = 1
  integer, parameter :: srequest = 2

  integer, parameter :: defsendtag = 1
  integer, parameter :: defrecvtag = 1

  integer, private, dimension(MPI_STATUS_SIZE) :: status

contains

subroutine haloirecvisendwaitall(sendbuf, recvbuf, n, neighbour, cartcomm)

  integer :: n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf
  integer, dimension(ndir, ndim) :: neighbour
  
  integer :: irep, idim, idir, irequest, reqid

  integer, dimension(MPI_STATUS_SIZE, ndir*ndim*nrequest) :: statuses
  integer, dimension(ndir*ndim*nrequest) :: requests

  ! Halo swap
  
  reqid = 0

  do irequest = 1, nrequest
     do idim = 1, ndim
        do idir = 1, ndir

           reqid = reqid + 1

           if (irequest == rrequest) then

              call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(ndir-idir+1, idim), defrecvtag, &
                             cartcomm, requests(reqid), ierr)

           else

              call MPI_Isend(sendbuf(1,idir,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(idir, idim), defsendtag, &
                             cartcomm, requests(reqid), ierr)

           end if

        end do
     end do
  end do

  call MPI_Waitall(ndir*ndim*nrequest, requests, &
                   statuses, ierr)


end subroutine haloirecvisendwaitall

subroutine olapirecvisendwaitall(x, calcrep, sendbuf, recvbuf, n, &
                                 neighbour, cartcomm)

  integer :: calcrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf
  integer, dimension(ndir, ndim) :: neighbour
  
  integer :: irep, idim, idir, irequest, reqid

  integer, dimension(MPI_STATUS_SIZE, ndir*ndim*nrequest) :: statuses
  integer, dimension(ndir*ndim*nrequest) :: requests

  double precision :: x

  ! Halo swap
  
  reqid = 0

  do irequest = 1, nrequest
     do idim = 1, ndim
        do idir = 1, ndir

           reqid = reqid + 1

           if (irequest == rrequest) then

              call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(ndir-idir+1, idim), defrecvtag, &
                             cartcomm, requests(reqid), ierr)

           else

              call MPI_Isend(sendbuf(1,idir,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(idir, idim), defsendtag, &
                             cartcomm, requests(reqid), ierr)

           end if

        end do
     end do
  end do

  call dummycalc(x, calcrep)

  call MPI_Waitall(ndir*ndim*nrequest, requests, &
                   statuses, ierr)


end subroutine olapirecvisendwaitall

subroutine nolapirecvisendwaitall(x, calcrep, sendbuf, recvbuf, n, &
                                  neighbour, cartcomm)

  integer :: calcrep, n, cartcomm
  double precision, dimension(n, ndir, ndim) :: sendbuf, recvbuf
  integer, dimension(ndir, ndim) :: neighbour
  
  integer :: irep, idim, idir, irequest, reqid

  integer, dimension(MPI_STATUS_SIZE, ndir*ndim*nrequest) :: statuses
  integer, dimension(ndir*ndim*nrequest) :: requests

  double precision :: x

  ! Halo swap
  
  reqid = 0

  do irequest = 1, nrequest
     do idim = 1, ndim
        do idir = 1, ndir

           reqid = reqid + 1

           if (irequest == rrequest) then

              call MPI_Irecv(recvbuf(1,ndir-idir+1,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(ndir-idir+1, idim), defrecvtag, &
                             cartcomm, requests(reqid), ierr)

           else

              call MPI_Isend(sendbuf(1,idir,idim), &
                             n, MPI_DOUBLE_PRECISION, &
                             neighbour(idir, idim), defsendtag, &
                             cartcomm, requests(reqid), ierr)

           end if

        end do
     end do
  end do

  call MPI_Waitall(ndir*ndim*nrequest, requests, &
                   statuses, ierr)

  call dummycalc(x, calcrep)

end subroutine nolapirecvisendwaitall

subroutine commdata(cartcomm, dims, periods, size, rank, coords, neighbour)

  integer :: cartcomm, size, rank

  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  integer :: dimsize, idim

  call MPI_Comm_size(cartcomm, size, ierr)
  call MPI_Comm_rank(cartcomm, rank, ierr)

  call MPI_Cartdim_get(cartcomm, dimsize, ierr)

  if (dimsize /= ndim) then

     if (rank == 0) write(*,*) "commdata: comm dimension ", dimsize, &
          " not equal to ndim = ", ndim
     call MPI_Finalize(ierr)
     stop

  end if

  call MPI_Cart_get(cartcomm, ndim, dims, periods, coords, ierr)

  !
  ! Find neighbours for this process
  !

  do idim = 1, ndim

     call MPI_Cart_shift(cartcomm, idim-1, 1, &
          neighbour(dndir, idim), &
          neighbour(updir, idim), ierr)

  end do

end subroutine commdata


subroutine checkrecvdata(flag, recvbuf, nbuf, cartcomm)

  logical, dimension(ndir, ndim) :: flag
  double precision, dimension(nbuf, ndir, ndim) :: recvbuf
  integer :: nbuf, cartcomm
  
  integer :: size, rank
  integer :: idim, idir
  
  integer, dimension(ndir, ndim) :: neighbour
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: coords, dims

  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, nullbuf

  ! Check that the correct data is sent and received

  call commdata(cartcomm, dims, periods, size, rank, coords, neighbour)
  
  flag(:,:) = .true.

  do idim = 1, ndim
     do idir = 1, ndir

        ! Compute the send halo data for this neighbour

        call inithalodata(neighbour(idir, idim), sendbuf, nullbuf, nbuf)

        if(any(recvbuf(:, idir, idim) /= sendbuf(:, ndir-idir+1, idim))) then
!           write(*,*) "checkrecvdata: error for rank, idim, idir = ", &
!                rank, idim, idir
!                "r/s = ", recvbuf(:,idir,idim), sendbuf(:,ndir-idir+1,idim) 

           flag(idir, idim) = .false.
           
        end if

     end do
  end do

end subroutine checkrecvdata


subroutine inithalodata(rank, sendbuf, recvbuf, nbuf)

  double precision, dimension(nbuf, ndir, ndim) :: sendbuf, recvbuf

  integer :: nbuf, rank

  integer :: idim, idir, ibuf

  ! set receive buffer to have illegal values

 recvbuf(:,:,:) = -1
  
! Set send buffer to have unique values

  do idim = 1, ndim
     do idir = 1, ndir
        do ibuf = 1, nbuf

           sendbuf(ibuf, idir, idim) =   rank*(nbuf*ndir*ndim) &
                                       + (idim-1)*nbuf*ndir  &
                                       + (idir-1)*nbuf + ibuf

!           write(*,*) 'rank = ', rank, ', x(', ibuf, ', ', idir, ', ', idim, &
!                ') = ', sendbuf(ibuf,idir,idim)

        end do
     end do
  end do

end subroutine inithalodata

subroutine printhaloerr(allflag)

  logical, dimension(:,:,:), allocatable :: allflag

  integer :: idir, idim
  
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

end subroutine printhaloerr
           
end module olapcomms
