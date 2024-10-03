subroutine sea_ice_post (tstep, expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

    integer, intent(in) :: tstep
    character(len=4), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb 
    character(len=20) :: filetfx, filetfy, filefcx, filefcy, filefbx, &
                         filefby, filem, filemc, filemb
    character(len=20) :: filetsigxx, filetsigyy, filetsigxy, filetsigyx

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
	filem   = "output/m." // trim(adjustl(expno_str))
	filemc  = "output/mc." // trim(adjustl(expno_str))
	filemb  = "output/mb." // trim(adjustl(expno_str))

    ! stress files
    filetsigxx = "output/tsigxx." // trim(adjustl(expno_str))
    filetsigyy = "output/tsigyy." // trim(adjustl(expno_str))
    filetsigxy = "output/tsigxy." // trim(adjustl(expno_str))
    filetsigyx = "output/tsigyx." // trim(adjustl(expno_str))

    ! writing in the files
    ! physical properties
    call write_file_1D(r, filer, "radius", "meters", expno_str, "r",   &
                    tstep)
    call write_file_1D(h, fileh, "thickness", "meters", expno_str, "h",&
                    tstep)
    call write_file_int_2D(bond, fileb, "bond pairs", "logical",       &
                    expno_str, "bond", tstep)
    
    ! position
    call write_file_1D(x, filex, "x position", "meters", expno_str,    &
                    "x", tstep)
    call write_file_1D(y, filey, "y position", "meters", expno_str,    &
                    "y", tstep)
    call write_file_1D(theta, filet, "angular position", "degrees",    &
                    expno_str, "theta", tstep)
    
    ! velocity
    call write_file_1D(u, fileu, "x velocity", "meters per second",    &
                    expno_str, "u", tstep)
    call write_file_1D(v, filev, "y velocity", "meters per second",    &
                    expno_str, "v", tstep)
    call write_file_1D(omega, fileo, "angular velocity",               &
                    "degrees per second", expno_str, "omega", tstep)
    
    ! forces
    call write_file_1D(tfx, filetfx, "total force in x", "Newtons",    &
                    expno_str, "tfx", tstep)
    call write_file_1D(tfy, filetfy, "total force in y", "Newtons",    &
                    expno_str, "tfy", tstep)
    call write_file_1D(fcx, filefcx, "contact force in x", "Newtons",  &
                    expno_str, "fcx", tstep)
    call write_file_1D(fcy, filefcy, "contact force in y", "Newtons",  &
                    expno_str, "fcy", tstep)
    call write_file_1D(fbx, filefbx, "bond force in x", "Newtons",     &
                    expno_str, "fbx", tstep)
    call write_file_1D(fby, filefby, "bond force in y", "Newtons",     &
                    expno_str, "fby", tstep)
    
    ! moments
    call write_file_1D(m, filem, "total moment", "Newton meters",      &
                    expno_str, "m", tstep)
    call write_file_1D(mc, filemc, "contact moment", "Newton meters",  &
                    expno_str, "mc", tstep)
    call write_file_1D(mb, filemb, "bond moment", "Newton meters",     &
                    expno_str, "mb", tstep)
    
    ! total stresses
    call write_file_1D(tsigxx, filetsigxx, "stress along x",           &
                    "Newtons per meter", expno_str, "tsigxx", tstep)
    call write_file_1D(tsigyy, filetsigyy, "stress along y",           &
                    "Newtons per meter", expno_str, "tsigyy", tstep)
    call write_file_1D(tsigxy, filetsigxy,                             &
                    "stress along x in + y-plane",                     &
                    "Newtons per meter", expno_str, "tsigxy", tstep)
    call write_file_1D(tsigyx, filetsigyx,                             &
                    "stress along y in + x-plane",                     &
                    "Newtons per meter", expno_str, "tsigyx", tstep)


end subroutine sea_ice_post


subroutine clear_posts (expno_str)

implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, stat (23)
	character(len=2), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb 
    character(len=20) :: filetfx, filetfy, filefcx, filefcy, filefbx, &
                         filefby, filem, filemc, filemb
    character(len=20) :: filetsigxx, filetsigyy, filetsigxy, filetsigyx
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
	filem   = "output/m." // trim(adjustl(expno_str))
	filemc  = "output/mc." // trim(adjustl(expno_str))
	filemb  = "output/mb." // trim(adjustl(expno_str))

    ! stress files
    filetsigxx = "output/tsigxx." // trim(adjustl(expno_str))
    filetsigyy = "output/tsigyy." // trim(adjustl(expno_str))
    filetsigxy = "output/tsigxy." // trim(adjustl(expno_str))
    filetsigyx = "output/tsigyx." // trim(adjustl(expno_str))

    ! info file
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
    open (19, file = filetfx, iostat = stat(10), status = 'old')
	open (20, file = filetfy, iostat = stat(11), status = 'old')
    open (21, file = filefcx, iostat = stat(12), status = 'old')
	open (22, file = filefcy, iostat = stat(13), status = 'old')
	open (23, file = filefbx, iostat = stat(14), status = 'old')
	open (24, file = filefby, iostat = stat(15), status = 'old')
	open (25, file = filem, iostat = stat(16), status = 'old')
	open (26, file = filemc, iostat = stat(17), status = 'old')
	open (27, file = filemb, iostat = stat(18), status = 'old')
    open (28, file = filetsigxx, iostat = stat(19), status = 'old')
    open (29, file = filetsigyy, iostat = stat(20), status = 'old')
    open (30, file = filetsigxy, iostat = stat(21), status = 'old')
    open (31, file = filetsigyx, iostat = stat(22), status = 'old')
    open (32, file = fileinfo, iostat = stat(23), status = 'old')

    do i = 10, 32
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


subroutine write_file_1D(array, path, description, units, expno_str, var, tstep)

    use netcdf

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    double precision, intent(in) :: array(n)
    character(len=20), intent(in) :: path
    character(len=*), intent(in) :: var, units, description
    character(len=4), intent(in) :: expno_str
    integer, intent(in) :: tstep


    integer :: file_id, xdim_id, time_id
    integer :: array_id
    integer, dimension(2) :: arrdims

    integer :: i
    integer :: ierr, time_index

    logical :: file_exists

    ! File status
    inquire(file=path, exist=file_exists)

    i = size(array,1)

    ! Create or open the NetCDF file
    if (file_exists) then
        ! Open the existing file in write mode
        ierr = nf90_open(path, NF90_WRITE, file_id)
    else
        ! Create a new file if it doesn't exist
        ierr = nf90_create(path=path, cmode=NF90_CLOBBER, ncid=file_id)
        
        ! Define the dimensions
        ierr = nf90_def_dim(file_id, 'i', i, xdim_id)
        ierr = nf90_def_dim(file_id, 'tstep', &
                            int(nt/comp), time_id)

        arrdims = (/time_id, xdim_id/)           

        ierr = nf90_def_var(file_id, var,  NF90_DOUBLE, arrdims, array_id)

        ! Assign attributes to the array
        ierr = nf90_put_att(file_id, array_id, 'units', units)
        ierr = nf90_put_att(file_id, array_id, 'dt', &
                            comp*dt)
        
        ! assign attributes to the file
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'expno', expno_str)
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'Description', &
                            description)


        ! End define mode
        ierr = nf90_enddef(file_id)
    endif

    ! which variable are you looking at
    ierr = nf90_inq_varid(file_id, var, array_id)

    ! Find the next available index for time dimension (appending)
    time_index = int(tstep/comp) 

    ! Write the array data to the appropriate time index
    ierr = nf90_put_var(file_id, array_id, array, &
                        start=(/time_index, 1/), count=(/1, i/))

    ! Close the NetCDF file
    ierr = nf90_close(file_id)

end subroutine write_file_1D


subroutine write_file_int_2D(array, path, description, units, expno_str, var, tstep)

    use netcdf

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: array(n,n)
    character(len=20), intent(in) :: path
    character(len=*), intent(in) :: var, units, description
    character(len=4), intent(in) :: expno_str
    integer, intent(in) :: tstep

    integer :: file_id, xdim_id, ydim_id, time_id
    integer :: array_id
    integer, dimension(3) :: arrdims

    integer :: i, j
    integer :: ierr, time_index

    logical :: file_exists

    ! File status
    inquire(file=path, exist=file_exists)

    i = size(array,1)
    j = size(array,2)

    ! Create or open the NetCDF file
    if (file_exists) then
        ! Open the existing file in write mode
        ierr = nf90_open(path, NF90_WRITE, file_id)
    else
        ! Create a new file if it doesn't exist
        ierr = nf90_create(path=path, cmode=NF90_CLOBBER, ncid=file_id)
        
        ! Define the dimensions
        ierr = nf90_def_dim(file_id, 'i', i, xdim_id)
        ierr = nf90_def_dim(file_id, 'j', j, ydim_id)
        ierr = nf90_def_dim(file_id, 'tstep', &
                            int(nt/comp), time_id)

        ! Define the variable based on the dimensions
        arrdims = (/time_id, xdim_id, ydim_id/)
        ierr = nf90_def_var(file_id, var,  NF90_INT, arrdims, array_id)

        ! Assign attributes to the array
        ierr = nf90_put_att(file_id, array_id, 'units', units)
        ierr = nf90_put_att(file_id, array_id, 'dt', &
                            comp*dt)

        ! assign attributes to the file
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'expno', expno_str)
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'Description', &
                            description)


        ! End define mode
        ierr = nf90_enddef(file_id)
    endif

    ! which variable are you looking at
    ierr = nf90_inq_varid(file_id, var, array_id)

    ! Find the next available index for time dimension (appending)
    time_index = int(tstep/comp)

    ! Write the array data to the appropriate time index
    ierr = nf90_put_var(file_id, array_id, array, &
                        start=(/time_index, 1, 1/), count=(/1, i, j/))

    ! Close the NetCDF file
    ierr = nf90_close(file_id)

end subroutine write_file_int_2D


