program ice

    use datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    use datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    use datetime, ONLY: str2dt, datetime_str_6, now, delta_str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: tstep, expno
    type(datetime_type) :: tic, tac

    tic = now()

    print *, 'experiment #?'
    read  *, expno 

    call ini_get
    call get_default
    
    do tstep = 1, int(t / dt) + 1

        call stepper (tstep)
        call sea_ice_post (tstep, expno)

    end do

    call execute_command_line("/aos/home/asavard/anaconda3/bin/python /storage/asavard/DEM/plots/video.py")

    tac = now()

    print*, "Total execution time: ", delta_str(tac-tic)
    
end program ice