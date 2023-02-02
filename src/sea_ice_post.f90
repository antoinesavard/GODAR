subroutine sea_ice_post (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, filer, fileh, filet, fileo

	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))

	open (10, file = filex, access = 'append', status = 'unknown')
	open (11, file = filey, access = 'append', status = 'unknown')
	open (12, file = filer, access = 'append', status = 'unknown')
	open (13, file = fileh, access = 'append', status = 'unknown')
	open (14, file = filet, access = 'append', status = 'unknown')
	open (15, file = fileo, access = 'append', status = 'unknown')

	write(10,*) ( x(i),                i=1, n )
	write(11,*) ( y(i),                i=1, n )
	write(12,*) ( r(i),                i=1, n )
	write(13,*) ( h(i),                i=1, n )
	write(14,*) ( theta(i),            i=1, n )
	write(15,*) ( sign(1d0, omega(i)), i=1, n )

    do i = 10, 15
        close(i)
    end do

end subroutine sea_ice_post


subroutine clear_posts (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, stat (6)
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, filer, fileh, filet, fileo

	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))

	open (10, file = filex, iostat = stat(1), status = 'old')
	open (11, file = filey, iostat = stat(2), status = 'old')
	open (12, file = filer, iostat = stat(3), status = 'old')
	open (13, file = fileh, iostat = stat(4), status = 'old')
	open (14, file = filet, iostat = stat(5), status = 'old')
	open (15, file = fileo, iostat = stat(6), status = 'old')

    do i = 10, 15
        if (stat(i-9) .eq. 0) then
			close(i, status = 'delete') 
		else 
			close(i)
		end if
    end do

end subroutine clear_posts


