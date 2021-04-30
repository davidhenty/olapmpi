module olapcalc

  implicit none

contains

subroutine dummycalc(x, nrep)

  integer :: nrep, irep, i
  double precision :: x

  do irep = 1, nrep
     do i = 1, 100
        x = sqrt(x)
     end do
  end do
  
end subroutine dummycalc

end module olapcalc
