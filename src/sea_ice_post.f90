subroutine sea_ice_post (tstep, expno)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i
    integer, intent(in) :: tstep, expno
    character(len=20) :: filename

    write (filename,'("output/x",i9.9,".", i2.2)') tstep, expno
    open (10, file = filename, status = 'unknown')

    write (filename,'("output/y",i9.9,".", i2.2)') tstep, expno
    open (11, file = filename, status = 'unknown')

    write (filename,'("output/r",i9.9,".", i2.2)') tstep, expno
    open (12, file = filename, status = 'unknown')

    write (filename,'("output/h",i9.9,".", i2.2)') tstep, expno
    open (13, file = filename, status = 'unknown')
    
    do i = 1, n
        write(10,*) x(i)
        write(11,*) y(i)
        write(12,*) r(i)
        write(13,*) h(i)
    end do

    do i = 10, 13
        close(i)
    end do

    
end subroutine sea_ice_post
