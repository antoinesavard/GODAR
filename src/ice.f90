program ice

    use datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    use datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    use datetime, ONLY: str2dt, datetime_str_6, now, delta_str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j, tstep
    type(datetime_type) :: tic, tac
    character(len=20) :: filename
    integer :: expno

    tic = now()

    print *, 'experiment #?'
    read  *, expno 

    call get_default
    call ini_get
    
    do tstep = 1, int(t / dt) + 1

        call stepper (tstep)
        write (filename,'("output/xy",i4.4,".", i2.2)') tstep, expno
        open (10, file = filename, status = 'unknown')
        do i = 1, n
            write(10,*) ( xy(i,j), j = 1, 4 )
        end do

    end do

    tac = now()

    print*, "Total execution time: ", delta_str(tac-tic)
    
end program ice