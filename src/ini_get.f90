subroutine ini_get

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i
    character Xfile*32, Yfile*32, Rfile*32, Hfile*32

    ! position of particles: need a one column prepared text file
    
    Xfile = 'src/X.dat'
    Yfile = 'src/Y.dat'
    Rfile = 'src/R.dat'
    Hfile = 'src/H.dat'

    open(102, file = Xfile, status='old')
    open(103, file = Yfile, status='old')
    open(104, file = Rfile, status='old')
    open(105, file = Hfile, status='old')

    do i = 1, n

        read(102,*) x(i)
        read(103,*) y(i)
        read(104,*) r(i)
        read(105,*) h(i)
        
    end do

    do i = 102, 105
        close(i)
    end do

    ! initial velocity

    u     =  0
    v     =  0
    omega =  0
    
    ! initial particle angle
    teta  =  0 

    ! initialize arrays
    tfx     =  0d0
    tfy     =  0d0
    fcn     =  0d0
    fct     =  0d0

end subroutine ini_get