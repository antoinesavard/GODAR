subroutine sea_ice_post (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

    integer :: i, j
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb

	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
    fileu = "output/u." // trim(adjustl(expno_str))
	filev = "output/v." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))
	fileb = "output/bond." // trim(adjustl(expno_str))

	open (10, file = filex, access = 'append', status = 'unknown')
	open (11, file = filey, access = 'append', status = 'unknown')
    open (12, file = fileu, access = 'append', status = 'unknown')
	open (13, file = filev, access = 'append', status = 'unknown')
	open (14, file = filer, access = 'append', status = 'unknown')
	open (15, file = fileh, access = 'append', status = 'unknown')
	open (16, file = filet, access = 'append', status = 'unknown')
	open (17, file = fileo, access = 'append', status = 'unknown')
	open (18, file = fileb, access = 'append', status = 'unknown')

	write(10,*) ( x(i),   	i=1, n )
	write(11,*) ( y(i),    	i=1, n )
    write(12,*) ( u(i),   	i=1, n )
	write(13,*) ( v(i),    	i=1, n )
	write(14,*) ( r(i),   	i=1, n )
	write(15,*) ( h(i),  	i=1, n )
	write(16,*) ( theta(i),	i=1, n )
	write(17,*) ( omega(i),	i=1, n )
	do i = 1, n
		write(18,*) ( bond(j,i), j=1, n )
	end do
	write(18,*)

    do i = 10, 18
        close(i)
    end do

end subroutine sea_ice_post


subroutine clear_posts (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, stat (9)
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb

	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
    fileu = "output/u." // trim(adjustl(expno_str))
	filev = "output/v." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))
	fileb = "output/bond." // trim(adjustl(expno_str))

	open (10, file = filex, iostat = stat(1), status = 'old')
	open (11, file = filey, iostat = stat(2), status = 'old')
    open (12, file = fileu, iostat = stat(3), status = 'old')
	open (13, file = filev, iostat = stat(4), status = 'old')
	open (14, file = filer, iostat = stat(5), status = 'old')
	open (15, file = fileh, iostat = stat(6), status = 'old')
	open (16, file = filet, iostat = stat(7), status = 'old')
	open (17, file = fileo, iostat = stat(8), status = 'old')
	open (18, file = fileb, iostat = stat(9), status = 'old')

    do i = 10, 18
        if (stat(i-9) .eq. 0) then
			close(i, status = 'delete') 
		else 
			close(i)
		end if
    end do

end subroutine clear_posts


