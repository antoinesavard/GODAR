program ice

    use datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    use datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    use datetime, ONLY: str2dt, datetime_str_6, now, delta_str

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: tstep
    integer :: expno, readnamelist, restart, expno_r
    type(datetime_type) :: tic, tac
    character(len=2) :: expno_str, expno_str_r

    !-------------------------------------------------------------------
    !       Read run information
    !-------------------------------------------------------------------

    tic = now()

    print *, 'Read namelist?'
    read  *, readnamelist
    print *, readnamelist

	print *, 'Restart experiment from previous one?'
	read  *, restart
	print *, restart

	if (restart .eq. 1) then
         read(*, '(i2)') expno_r
         write(*,*) "Restart experiment number is: ", expno_r
		 write(expno_str_r,'(i2.2)') expno_r
    endif

    print *, 'Experiment #?'
    read  *, expno 
    print *, expno 
    write(expno_str,'(i2.2)') expno

    call ini_get (restart, expno_str_r)

    call get_default
	! overwrite default based on namelist
    if (readnamelist .eq. 1) then
        call read_namelist
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
