subroutine sea_ice_post (tstep, expno_str)

    use parameters
    use variables
    use const
    use bonds
    use diagnostics

    implicit none


    integer :: i, j
    integer, intent(in) :: tstep
    character(len=4), intent(in) :: expno_str
    character(len=20) :: filex, filey, fileu, filev, filer, fileh, &
                         filet, fileo, fileb 
    character(len=20) :: filetfx, filetfy, filefcx, filefcy, filefbx, &
                         filefby, filem, filemc, filemb
    character(len=20) :: filetsigxx, filetsigyy, filetsigxy, filetsigyx
    character(len=20) :: filetp, fileangle, filedamage, filebeam

    ! bending angles of a bond and overlap of a pair
    double precision :: thetarelb_ji, thetarelb_ij, deltan_ji

    ! output step index written in the sparse files
    integer :: idx

    ! Explicit output formats. List-directed output (write(u,*)) is
    ! processor dependent. These formats give the
    ! same bytes with every compiler.
    !   ES24.16E3: 17 significant digits (exact round trip of a double),
    !              3-digit exponent so |exp| > 99 keeps its 'E'
    !   I0:        minimal-width integers
    character(len=*), parameter :: fmt_r    = 'ES24.16E3'
    character(len=*), parameter :: fmt_row  = '(*(1X,' // fmt_r // '))'
    character(len=*), parameter :: fmt_bond = '(I0,1X,I0,1X,I0)'
    character(len=*), parameter :: fmt_dmg  = '(I0,1X,I0,1X,I0,1X,' // fmt_r // ')'
    character(len=*), parameter :: fmt_beam = '(I0,1X,I0,1X,I0,5(1X,' // fmt_r // '))'

    ! diagnostics variables
    num_bonds = 0

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

    ! diagnostic files
    filetp = "output/tp." // trim(adjustl(expno_str))
    fileangle = "output/angle." // trim(adjustl(expno_str))
    filedamage = "output/damage." // trim(adjustl(expno_str))
    filebeam = "output/beam." // trim(adjustl(expno_str))

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
    call open_out(10, filex)
    call open_out(11, filey)
    call open_out(12, fileu)
    call open_out(13, filev)
    call open_out(14, filer)
    call open_out(15, fileh)
    call open_out(16, filet)
    call open_out(17, fileo)
    call open_out(18, fileb)
    ! forces
    call open_out(19, filetfx)
    call open_out(20, filetfy)
    call open_out(21, filefcx)
    call open_out(22, filefcy)
    call open_out(23, filefbx)
    call open_out(24, filefby)
    ! moments
    call open_out(25, filem)
    call open_out(26, filemc)
    call open_out(27, filemb)
    ! stresses
    call open_out(28, filetsigxx)
    call open_out(29, filetsigyy)
    call open_out(30, filetsigxy)
    call open_out(31, filetsigyx)
    ! pressure
    call open_out(32, filetp)
    ! angle
    call open_out(33, fileangle)
    ! damage
    call open_out(34, filedamage)
    ! beam state, for restarts
    call open_out(35, filebeam)

    ! write on the files: one record (line) per output step
    write(10, fmt_row) ( x(i),     i=1, n )
    write(11, fmt_row) ( y(i),     i=1, n )
    write(12, fmt_row) ( u(i),     i=1, n )
    write(13, fmt_row) ( v(i),     i=1, n )
    write(14, fmt_row) ( r(i),     i=1, n )
    write(15, fmt_row) ( h(i),     i=1, n )
    write(16, fmt_row) ( theta(i), i=1, n )
    write(17, fmt_row) ( omega(i), i=1, n )

    idx = int(tstep / comp)
    do i = 1, n
        do j = 1, n
            if (bond(j, i) == 1) then
                ! bending angles, as in bond_forces_timoshenko
                thetarelb_ji = -(theta(i) - alpha_total(j,i) - theta_offset(j,i))
                thetarelb_ij = -(theta(j) - alpha_total(j,i) - theta_offset(i,j))

                write(18, fmt_bond) idx, j, i
                write(34, fmt_dmg)  idx, j, i, damageb(j,i)
                write(35, fmt_beam) idx, j, i, lb(j,i), rb(j,i), &
                                    hb(j,i), thetarelb_ji, thetarelb_ij
                num_bonds = num_bonds + 1
            end if

            ! contacts between active particles (lower triangle),
            ! overlap recomputed from the positions
            if ( j .gt. i .and. active(i) .and. active(j) ) then
                deltan_ji = r(i) + r(j) - sqrt( ( x(i) - x(j) ) ** 2 + &
                                                ( y(i) - y(j) ) ** 2 )
                if ( deltan_ji .gt. 0 ) then
                    write(33, fmt_dmg) idx, j, i, &
                        atan2( y(j) - y(i), x(j) - x(i) ) * 180d0 / pi
                end if
            end if
        end do
    end do
    ! empty line separates output steps in the bond file
    write(18, '(A)') ''

    write(19, fmt_row) ( tfx(i),    i=1, n )
    write(20, fmt_row) ( tfy(i),    i=1, n )
    write(21, fmt_row) ( fcx(i),    i=1, n )
    write(22, fmt_row) ( fcy(i),    i=1, n )
    write(23, fmt_row) ( fbx(i),    i=1, n )
    write(24, fmt_row) ( fby(i),    i=1, n )
    write(25, fmt_row) ( m(i),      i=1, n )
    write(26, fmt_row) ( mc(i),     i=1, n )
    write(27, fmt_row) ( mb(i),     i=1, n )
    write(28, fmt_row) ( tsigxx(i), i=1, n )
    write(29, fmt_row) ( tsigyy(i), i=1, n )
    write(30, fmt_row) ( tsigxy(i), i=1, n )
    write(31, fmt_row) ( tsigyx(i), i=1, n )
    write(32, fmt_row) ( tp(i),     i=1, n )

    do i = 10, 35
        close(i)
    end do

    ! diagnostics stdout
    if (tstep .eq. 0) then
        num_bonds_nm1 = num_bonds
    end if
    write(*, '(A,ES3.5,A,F4.3,A)') " Mean velocity: ",          &
        sqrt((sum(u) / n) ** 2 + (sum(v) / n)** 2), " m/s at ", &
        atan2(sum(v) / n, sum(u) / n) * 180d0 / pi, " deg"
    write(*, '(A,ES3.5,A)') " RMS speed: ", sqrt(sum(u**2 + v**2) / n), " m/s"
    write(*, '(A,I0)') " Number of bonds: ", num_bonds
    write(*, '(A,I0)') " Number of bonds broken: ", num_bonds - num_bonds_nm1

    ! update number of bonds at tstep-1 for next tstep
    num_bonds_nm1 = num_bonds

contains

    subroutine open_out(unit, filename)
        ! Open an output file for appending

        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        integer :: ios

        open(unit, file = filename, status = 'unknown', &
             position = 'append', action = 'write',     &
             access = 'sequential', form = 'formatted', &
             recl = 25 * max(n, 1) + 256, iostat = ios)

        if (ios /= 0) then
            write(*, '(A,A,A,I0)') " Error opening ",   &
            trim(filename), ", iostat = ", ios
            stop 1
        end if

    end subroutine open_out

end subroutine sea_ice_post
