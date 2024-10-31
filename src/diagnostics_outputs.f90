subroutine sea_ice_post (tstep, expno_str)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"
    include "CB_diagnostics.h"

    integer :: i, j
    integer, intent(in) :: tstep
    character(len=4), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb 
    character(len=20) :: filetfx, filetfy, filefcx, filefcy, filefbx, &
                         filefby, filem, filemc, filemb
    character(len=20) :: filetsigxx, filetsigyy, filetsigxy, filetsigyx
    character(len=20) :: filetp

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

    ! pressure file
    filetp = "output/tp." // trim(adjustl(expno_str))

    ! ! writing in the files
    ! ! physical properties
    ! call write_file_1D(r, filer, "radius", "meters", expno_str, "r",   &
    !                 tstep)
    ! call write_file_1D(h, fileh, "thickness", "meters", expno_str, "h",&
    !                 tstep)
    ! call write_file_int_2D(bond, fileb, "bond pairs", "logical",       &
    !                 expno_str, "bond", tstep)
    
    ! ! position
    ! call write_file_1D(x, filex, "x position", "meters", expno_str,    &
    !                 "x", tstep)
    ! call write_file_1D(y, filey, "y position", "meters", expno_str,    &
    !                 "y", tstep)
    ! call write_file_1D(theta, filet, "angular position", "degrees",    &
    !                 expno_str, "theta", tstep)
    
    ! ! velocity
    ! call write_file_1D(u, fileu, "x velocity", "meters per second",    &
    !                 expno_str, "u", tstep)
    ! call write_file_1D(v, filev, "y velocity", "meters per second",    &
    !                 expno_str, "v", tstep)
    ! call write_file_1D(omega, fileo, "angular velocity",               &
    !                 "degrees per second", expno_str, "omega", tstep)
    
    ! ! forces
    ! call write_file_1D(tfx, filetfx, "total force in x", "Newtons",    &
    !                 expno_str, "tfx", tstep)
    ! call write_file_1D(tfy, filetfy, "total force in y", "Newtons",    &
    !                 expno_str, "tfy", tstep)
    ! call write_file_1D(fcx, filefcx, "contact force in x", "Newtons",  &
    !                 expno_str, "fcx", tstep)
    ! call write_file_1D(fcy, filefcy, "contact force in y", "Newtons",  &
    !                 expno_str, "fcy", tstep)
    ! call write_file_1D(fbx, filefbx, "bond force in x", "Newtons",     &
    !                 expno_str, "fbx", tstep)
    ! call write_file_1D(fby, filefby, "bond force in y", "Newtons",     &
    !                 expno_str, "fby", tstep)
    
    ! ! moments
    ! call write_file_1D(m, filem, "total moment", "Newton meters",      &
    !                 expno_str, "m", tstep)
    ! call write_file_1D(mc, filemc, "contact moment", "Newton meters",  &
    !                 expno_str, "mc", tstep)
    ! call write_file_1D(mb, filemb, "bond moment", "Newton meters",     &
    !                 expno_str, "mb", tstep)
    
    ! ! total stresses
    ! call write_file_1D(tsigxx, filetsigxx, "stress along x",           &
    !                 "Newtons per meter", expno_str, "tsigxx", tstep)
    ! call write_file_1D(tsigyy, filetsigyy, "stress along y",           &
    !                 "Newtons per meter", expno_str, "tsigyy", tstep)
    ! call write_file_1D(tsigxy, filetsigxy,                             &
    !                 "stress along x in + y-plane",                     &
    !                 "Newtons per meter", expno_str, "tsigxy", tstep)
    ! call write_file_1D(tsigyx, filetsigyx,                             &
    !                 "stress along y in + x-plane",                     &
    !                 "Newtons per meter", expno_str, "tsigyx", tstep)

    ! opening the files
    ! position and state
	open (10, file = filex, access = 'append', status = 'unknown')
	open (11, file = filey, access = 'append', status = 'unknown')
    open (12, file = fileu, access = 'append', status = 'unknown')
	open (13, file = filev, access = 'append', status = 'unknown')
	open (14, file = filer, access = 'append', status = 'unknown')
	open (15, file = fileh, access = 'append', status = 'unknown')
	open (16, file = filet, access = 'append', status = 'unknown')
	open (17, file = fileo, access = 'append', status = 'unknown')
	open (18, file = fileb, access = 'append', status = 'unknown')
    ! forces
    open (19, file = filetfx, access = 'append', status = 'unknown')
	open (20, file = filetfy, access = 'append', status = 'unknown')
    open (21, file = filefcx, access = 'append', status = 'unknown')
	open (22, file = filefcy, access = 'append', status = 'unknown')
	open (23, file = filefbx, access = 'append', status = 'unknown')
	open (24, file = filefby, access = 'append', status = 'unknown')
    ! moments
	open (25, file = filem, access = 'append', status = 'unknown')
	open (26, file = filemc, access = 'append', status = 'unknown')
	open (27, file = filemb, access = 'append', status = 'unknown')
    ! stresses
    open (28, file = filetsigxx, access = 'append', status = 'unknown')
    open (29, file = filetsigyy, access = 'append', status = 'unknown')
    open (30, file = filetsigxy, access = 'append', status = 'unknown')
    open (31, file = filetsigyx, access = 'append', status = 'unknown')
    ! pressure
    open (32, file = filetp, access = 'append', status = 'unknown')


    ! write on the files
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
    write(19,*) ( tfx(i),   i=1, n )
	write(20,*) ( tfy(i),   i=1, n )
    write(21,*) ( fcx(i),  	i=1, n )
	write(22,*) ( fcy(i),  	i=1, n )
	write(23,*) ( fbx(i),   i=1, n )
	write(24,*) ( fby(i),  	i=1, n )
	write(25,*) ( m(i),	    i=1, n )
	write(26,*) ( mc(i),	i=1, n )
    write(27,*) ( mb(i),	i=1, n )
    write(28,*) ( tsigxx(i),	i=1, n )
    write(29,*) ( tsigyy(i),	i=1, n )
    write(30,*) ( tsigxy(i),	i=1, n )
    write(31,*) ( tsigyx(i),	i=1, n )
    write(32,*) ( tp(i),	i=1, n )

    do i = 10, 32
        close(i)
    end do


end subroutine sea_ice_post

