subroutine sea_ice_post (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

    integer :: i, j
    character(len=4), intent(in) :: expno_str
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

    integer :: i, stat (10)
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb, fileinfo

	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
    fileu = "output/u." // trim(adjustl(expno_str))
	filev = "output/v." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))
	fileb = "output/bond." // trim(adjustl(expno_str))
    fileinfo = "output/info." // trim(adjustl(expno_str))

	open (10, file = filex, iostat = stat(1), status = 'old')
	open (11, file = filey, iostat = stat(2), status = 'old')
    open (12, file = fileu, iostat = stat(3), status = 'old')
	open (13, file = filev, iostat = stat(4), status = 'old')
	open (14, file = filer, iostat = stat(5), status = 'old')
	open (15, file = fileh, iostat = stat(6), status = 'old')
	open (16, file = filet, iostat = stat(7), status = 'old')
	open (17, file = fileo, iostat = stat(8), status = 'old')
	open (18, file = fileb, iostat = stat(9), status = 'old')
    open (19, file = fileinfo, iostat = stat(10), status = 'old')

    do i = 10, 19
        if (stat(i-9) .eq. 0) then
			close(i, status = 'delete') 
		else 
			close(i)
		end if
    end do

end subroutine clear_posts


subroutine info (expno_str, restart)

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"

    character(len=4), intent(in) :: expno_str
    integer, intent(in) :: restart
    character(len=30) :: fileinfo, ffmt, lfmt, ifmt

    ffmt = "(es25.3e1)"
    lfmt = "(l25)"
    ifmt = "(i25)"

    fileinfo = 'output/info.' // trim(adjustl(expno_str))

    open (10, file = fileinfo, status = 'unknown')

    ! options of the model
    write(10,*) ('OPTIONS')
    write(10,*) ('')
    write(10,*) ('dynamics =')
    write(10,lfmt) (dynamics)
    write(10,*) ('thermodyn =')
    write(10,lfmt) (thermodyn)
    write(10,*) ('cohesion =') 
    write(10,lfmt) (cohesion)
    write(10,*) ('ridging =') 
    write(10,lfmt) (ridging)
    write(10,*) ('shelter =')
    write(10,lfmt) (shelter)

    ! numerical parameters
    write(10,*) ('')
    write(10,*) ('NUMERICAL PARAMETERS')
    write(10,*) ('')
    write(10,*) ('rtree =')
    write(10,ffmt) (rtree)
    write(10,*) ('dt =')
    write(10,ffmt) (dt)
    write(10,*) ('nt =') 
    write(10,ffmt) (nt)
    write(10,*) ('comp =') 
    write(10,ffmt) (comp)
    write(10,*) ('n =')
    write(10,ifmt) (n)
    write(10,*) ('nx =')
    write(10,ffmt) (nx)
    write(10,*) ('ny =')
    write(10,ffmt) (ny)

    ! physical parameters
    write(10,*) ('')
    write(10,*) ('PHYSICAL PARAMETERS')
    write(10,*) ('')
    write(10,*) ('Cdair =')
    write(10,ffmt) (Cdair)
    write(10,*) ('Csair =')
    write(10,ffmt) (Csair)
    write(10,*) ('Cdwater =')
    write(10,ffmt) (Cdwater)
    write(10,*) ('Cswater =')
    write(10,ffmt) (Cswater)
    write(10,*) ('z0w =')
    write(10,ffmt) (z0w)
    write(10,*) ('lat =')
    write(10,ffmt) (lat)
    write(10,*) ('rhoair =')
    write(10,ffmt) (rhoair)
    write(10,*) ('rhoice =')
    write(10,ffmt) (rhoice)
    write(10,*) ('rhowater =')
    write(10,ffmt) (rhowater)

    ! disk parameters
    write(10,*) ('')
    write(10,*) ('DISK PARAMETERS')
    write(10,*) ('')
    write(10,*) ('e_modul = ')
    write(10,ffmt) (e_modul)
    write(10,*) ('poiss_ratio =')
    write(10,ffmt) (poiss_ratio)
    write(10,*) ('friction_coeff =')
    write(10,ffmt) (friction_coeff)
    write(10,*) ('rest_coeff =')
    write(10,ffmt) (rest_coeff)
    write(10,*) ('sigmanc_crit =')
    write(10,ffmt) (sigmanc_crit)

    ! bond parameters
    write(10,*) ('')
    write(10,*) ('BOND PARAMETERS')
    write(10,*) ('')
    write(10,*) ('eb = ')
    write(10,ffmt) (eb)
    write(10,*) ('lambda_rb =')
    write(10,ffmt) (lambda_rb)
    write(10,*) ('lambda_lb =')
    write(10,ffmt) (lambda_lb)
    write(10,*) ('sigmatb_crit =')
    write(10,ffmt) (sigmatb_crit)
    write(10,*) ('sigmacb_crit =')
    write(10,ffmt) (sigmacb_crit)
    write(10,*) ('tau_crit =')
    write(10,ffmt) (tau_crit)
    write(10,*) ('gamma_d =')
    write(10,ffmt) (gamma_d)

    ! forcings
    write(10,*) ('')
    write(10,*) ('FORCINGS')
    write(10,*) ('')
    write(10,*) ('uw = ')
    write(10,ffmt) (uw)
    write(10,*) ('vw =')
    write(10,ffmt) (vw)
    write(10,*) ('ua =')
    write(10,ffmt) (ua)
    write(10,*) ('va =')
    write(10,ffmt) (va)

    ! plate forcigns
    write(10,*) ('')
    write(10,*) ('PLATES')
    write(10,*) ('')
    write(10,*) ('pfn =')
    write(10,ffmt) (pfn)
    write(10,*) ('pfs =')
    write(10,ffmt) (pfs)

    ! input files
    write(10,*) ('')
    write(10,*) ('INPUT FILES')
    write(10,*) ('')
    write(10,*) ('restart =')
    write(10,ifmt) (restart)

end subroutine info


