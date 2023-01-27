!======================================================================
! Compaction algorithm
!======================================================================

subroutine compaction

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i
    double precision :: xraw(n), yraw(n), rraw(n), hraw(n)

    ! position of particles: need a one column prepared text file
    
    open(102, file = 'src/Xraw.dat')
    open(103, file = 'src/Yraw.dat')
    open(104, file = 'src/Rraw.dat')
    open(105, file = 'src/Hraw.dat')

    do i = 1, n
       
        read(102,*) xraw(i)
        read(103,*) yraw(i)
        read(104,*) rraw(i)                                            
        read(105,*) hraw(i)
        
    end do
    
    do i = 102, 105
        close(i)
    end do

end subroutine compaction

subroutine bubblesort (td, A, n)
  
    integer, intent(in) :: td
    integer, intent(in out), dimension(td) :: A
    integer, intent(in) :: n
    integer :: i, j, temp
    logical :: swapped

    do j = n-1, 1, -1

        swapped = .false.

        do i = 1, j
            if (A(i) > A(i+1)) then
                temp = A(i)
                A(i) = A(i+1)
                A(i+1) = temp
                swapped = .true.
            end if
        end do

        if (.not. swapped) exit
    end do
     
end subroutine bubblesort
