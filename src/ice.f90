program ice

    use datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    use datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    use datetime, ONLY: str2dt, datetime_str_6, now, delta_str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: tstep
    integer :: expno, readnamelist
    type(datetime_type) :: tic, tac
    character(len=2) expno_str

    !-------------------------------------------------------------------
    !       Read run information
    !-------------------------------------------------------------------

    tic = now()

    print *, 'Read namelist?'
    read  *, readnamelist
    print *, readnamelist

    print *, 'experiment #?'
    read  *, expno 
    print *, expno 
	write(expno_str,'(i2.2)') expno

    call ini_get

    call get_default
    if (readnamelist .eq. 1) then
        call read_namelist      ! overwrite default based on namelist
    endif

	call clear_posts (expno_str)

    do tstep = 1, int(nt) + 1

        call stepper (tstep)

        if (MODULO(tstep, int(nt/comp)) .eq. 0) then

            print *, "Time step: ", tstep
            call sea_ice_post (expno_str)

        endif

    end do

    tac = now()

    print*, "Total simulation time: ", delta_str(tac - tic)

    call execute_command_line("/aos/home/asavard/anaconda3/bin/python /storage/asavard/DEM/plots/video.py "//expno_str)

    tic = now()

    print*, "Total video creation time: ", delta_str(tic - tac)
    
end program ice
