subroutine clear_posts (expno_str)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, stat(24)
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb 
    character(len=20) :: filetfx, filetfy, filefcx, filefcy, filefbx, &
                         filefby, filem, filemc, filemb
    character(len=20) :: filetsigxx, filetsigyy, filetsigxy, filetsigyx
    character(len=20) :: filetp, fileangle
    character(len=20) :: fileinfo

    ! position and state files
	filex = "output/x." // trim(adjustl(expno_str))
	filey = "output/y." // trim(adjustl(expno_str))
    fileu = "output/u." // trim(adjustl(expno_str))
	filev = "output/v." // trim(adjustl(expno_str))
	filer = "output/r." // trim(adjustl(expno_str))
	fileh = "output/h." // trim(adjustl(expno_str))
	filet = "output/theta." // trim(adjustl(expno_str))
	fileo = "output/omega." // trim(adjustl(expno_str))
	fileb = "output/bond." // trim(adjustl(expno_str))
    
    ! force files
    filetfx = "output/tfx." // trim(adjustl(expno_str))
	filetfy = "output/tfy." // trim(adjustl(expno_str))
    filefcx = "output/fcx." // trim(adjustl(expno_str))
	filefcy = "output/fcy." // trim(adjustl(expno_str))
	filefbx = "output/fbx." // trim(adjustl(expno_str))
	filefby = "output/fby." // trim(adjustl(expno_str))
	filem   = "output/mom." // trim(adjustl(expno_str))
	filemc  = "output/mc." // trim(adjustl(expno_str))
	filemb  = "output/mb." // trim(adjustl(expno_str))

    ! stress files
    filetsigxx = "output/tsigxx." // trim(adjustl(expno_str))
    filetsigyy = "output/tsigyy." // trim(adjustl(expno_str))
    filetsigxy = "output/tsigxy." // trim(adjustl(expno_str))
    filetsigyx = "output/tsigyx." // trim(adjustl(expno_str))

    ! pressure file
    filetp = "output/tp." // trim(adjustl(expno_str))
    fileangle = "output/angle." // trim(adjustl(expno_str))

    ! info file
    fileinfo = "output/info." // trim(adjustl(expno_str))

    ! open the files
    ! position and state
	open (10, file = filex, iostat = stat(1), status = 'old')
	open (11, file = filey, iostat = stat(2), status = 'old')
    open (12, file = fileu, iostat = stat(3), status = 'old')
	open (13, file = filev, iostat = stat(4), status = 'old')
	open (14, file = filer, iostat = stat(5), status = 'old')
	open (15, file = fileh, iostat = stat(6), status = 'old')
	open (16, file = filet, iostat = stat(7), status = 'old')
	open (17, file = fileo, iostat = stat(8), status = 'old')
	open (18, file = fileb, iostat = stat(9), status = 'old')
    ! forces
    open (19, file = filetfx, iostat = stat(10), status = 'old')
	open (20, file = filetfy, iostat = stat(11), status = 'old')
    open (21, file = filefcx, iostat = stat(12), status = 'old')
	open (22, file = filefcy, iostat = stat(13), status = 'old')
	open (23, file = filefbx, iostat = stat(14), status = 'old')
	open (24, file = filefby, iostat = stat(15), status = 'old')
    ! moments
	open (25, file = filem, iostat = stat(16), status = 'old')
	open (26, file = filemc, iostat = stat(17), status = 'old')
	open (27, file = filemb, iostat = stat(18), status = 'old')
    ! stresses
    open (28, file = filetsigxx, iostat = stat(19), status = 'old')
    open (29, file = filetsigyy, iostat = stat(20), status = 'old')
    open (30, file = filetsigxy, iostat = stat(21), status = 'old')
    open (31, file = filetsigyx, iostat = stat(22), status = 'old')
    ! pressure
    open (32, file = filetp, iostat = stat(23), status = 'old')
    ! angle
    open (33, file = fileangle, iostat = stat(24), status = 'old')
    ! info
    open (34, file = fileinfo, iostat = stat(25), status = 'old')

    ! delete them
    do i = 10, 34
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
    character(len=30) :: fileinfo, ffmt, lfmt, ifmt, cfmt

    ffmt = "(A, ES10.3)"
    lfmt = "(A, L10)"
    ifmt = "(A, I10)"
    cfmt = "(A, A)"

    fileinfo = 'output/info.' // trim(adjustl(expno_str))

    open (10, file = fileinfo, status = 'unknown')

    ! options of the model
    write(10,*) ('OPTIONS')
    write(10,*) ('')
    write(10,lfmt) 'dynamics  =', dynamics
    write(10,lfmt) 'slipping  =', slipping
    write(10,lfmt) 'thermodyn =', thermodyn
    write(10,lfmt) 'cohesion  =', cohesion
    write(10,lfmt) 'ridging   =', ridging
    write(10,lfmt) 'shelter   =', shelter
    ! numerical parameters
    write(10,*) ('')
    write(10,*) ('NUMERICAL PARAMETERS')
    write(10,*) ('')
    write(10,ffmt) 'rtree =', rtree
    write(10,ifmt) 'ntree =', int(ntree)
    write(10,ffmt) 'dt    =', dt
    write(10,ffmt) 'nt    =', nt
    write(10,ffmt) 'comp  =', comp
    write(10,ifmt) 'n     =', n
    write(10,ffmt) 'nx    =', nx
    write(10,ffmt) 'ny    =', ny

    ! physical parameters
    write(10,*) ('')
    write(10,*) ('PHYSICAL PARAMETERS')
    write(10,*) ('')
    write(10,ffmt) 'Cdair    =', Cdair
    write(10,ffmt) 'Csair    =', Csair
    write(10,ffmt) 'Cdwater  =', Cdwater
    write(10,ffmt) 'Cswater  =', Cswater
    write(10,ffmt) 'z0w      =', z0w
    write(10,ffmt) 'lat      =', lat
    write(10,ffmt) 'rhoair   =', rhoair
    write(10,ffmt) 'rhoice   =', rhoice
    write(10,ffmt) 'rhowater =', rhowater

    ! disk parameters
    write(10,*) ('')
    write(10,*) ('DISK PARAMETERS')
    write(10,*) ('')
    write(10,ffmt) 'e_modul        =', e_modul
    write(10,ffmt) 'poiss_ratio    =', poiss_ratio
    write(10,ffmt) 'friction_coeff =', friction_coeff
    write(10,ffmt) 'rest_coeff     =', rest_coeff
    write(10,ffmt) 'sigmanc_crit   =', sigmanc_crit
    ! bond parameters
    write(10,*) ('')
    write(10,*) ('BOND PARAMETERS')
    write(10,*) ('')
    write(10,ffmt) 'eb           =', eb
    write(10,ffmt) 'lambda_rb    =', lambda_rb
    write(10,ffmt) 'lambda_lb    =', lambda_lb
    write(10,ffmt) 'sigmatb_crit =', sigmatb_crit
    write(10,ffmt) 'sigmacb_crit =', sigmacb_crit
    write(10,ffmt) 'tau_crit     =', tau_crit
    write(10,ffmt) 'gamma_d      =', gamma_d
    write(10,ffmt) 'bond_lim     =', bond_lim

    ! forcings
    write(10,*) ('')
    write(10,*) ('FORCINGS')
    write(10,*) ('')
    write(10,ffmt) 'uw =', uw
    write(10,ffmt) 'vw =', vw
    write(10,ffmt) 'ua =', ua
    write(10,ffmt) 'va =', va

    ! plate forcigns
    write(10,*) ('')
    write(10,*) ('PLATES')
    write(10,*) ('')
    write(10,ffmt) 'pfn =', pfn
    write(10,ffmt) 'pfs =', pfs

    ! input files
    write(10,*) ('')
    write(10,*) ('INPUT FILES')
    write(10,*) ('')
    write(10,ifmt) 'restart =', restart
    write(10,cfmt) 'Xfile   =', Xfile
    write(10,cfmt) 'Yfile   =', Yfile
    write(10,cfmt) 'Rfile   =', Rfile
    write(10,cfmt) 'Hfile   =', Hfile
    write(10,cfmt) 'Tfile   =', Tfile
    write(10,cfmt) 'Ofile   =', Ofile

end subroutine info