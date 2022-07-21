subroutine ini_get

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j

    ! position of particles: need a one column prepared text file
    
    open(102, file = 'src/X.dat')
    open(103, file = 'src/Y.dat')
    open(104, file = 'src/R.dat')
    open(105, file = 'src/M.dat')

    do i = 1, n

        read(102,*) xy(i,1)
        read(103,*) xy(i,2)
        read(104,*) xy(i,3)
        read(105,*) xy(i,4)
        
    end do

    do j = 102, 105
        close(j)
    end do

! initial velocity

    vu    =  0
    ome   =  0

! initial particle rotation
    teta  =  0 

!Allocate following array with zero value before main loop calculation
!(loop of contact force detection)
    fw      =  0d0
    f(:,2)  =  0d0
    f(:,1)  =  0d0
    s       =  0d0
    fsx     =  0d0
    fsy     =  0d0
    w       =  0d0
    fwx     =  0d0
    fwy     =  0d0
    d       =  0d0
    dx      =  0d0
    dy      =  0d0
    l       =  0d0
    lx      =  0d0
    ly      =  0d0

end subroutine ini_get