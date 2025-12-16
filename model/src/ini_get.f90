subroutine ini_get (restart, expno_str_r, nt_r)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_diagnostics.h"

	integer, intent(in) :: restart, nt_r
	character(2), intent(in) :: expno_str_r
    integer :: iostat

    integer :: i, j, k

	! load restart files
	if ( restart .eq. 1 ) then

        write(*,*) ('Reading last entry of restart version')

		Xfile = 'output/x.' // trim(adjustl(expno_str_r))
		Yfile = 'output/y.' // trim(adjustl(expno_str_r))
		Rfile = 'output/r.' // trim(adjustl(expno_str_r))
		Hfile = 'output/h.' // trim(adjustl(expno_str_r))
		Tfile = 'output/theta.' // trim(adjustl(expno_str_r))
		Ofile = 'output/omega.' // trim(adjustl(expno_str_r))
        Ufile = 'output/u.' // trim(adjustl(expno_str_r))
        Vfile = 'output/v.' // trim(adjustl(expno_str_r))
        Bfile = 'output/bond.' // trim(adjustl(expno_str_r))
		
        ! last timestep data index
		k = int(nt_r)

        ! read particle data at timestep k
        ! x, y, r, h, theta, omega, u, v, bond
        ! TODO: deltas to remember the stress state
        open(103, file = Xfile, status='old')
		do j = 1, k-1
			read (103, *, iostat=iostat)
            if (iostat /= 0) exit
        end do
        do j = k, k
			read (103, *, iostat=iostat) ( x(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(103)

        open(104, file = Yfile, status='old')
        do j = 1, k-1
			read (104, *, iostat=iostat)
            if (iostat /= 0) exit
        end do
        do j = k, k
			read (104, *, iostat=iostat) ( y(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(104)

        open(105, file = Rfile, status='old')
        do j = 1, k-1
			read (105, *, iostat=iostat)
            if (iostat /= 0) exit
        end do
        do j = k, k
			read (105, *, iostat=iostat) ( r(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(105)

        open(106, file = Hfile, status='old')
        do j = 1, k-1
			read (106, *, iostat=iostat)
            if (iostat /= 0) exit
        end do
        do j = k, k
			read (106, *, iostat=iostat) ( h(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(106)

        open(107, file = Tfile, status='old')
        do j = 1, k-1
			read (107, *, iostat=iostat)
            if (iostat /= 0) exit
        end do
        do j = k, k
			read (107, *, iostat=iostat) ( theta(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(107)

        open(108, file = Ofile, status='old')
        do j = 1, k-1
			read (108, *, iostat=iostat)
            if (iostat /= 0) exit
		end do
        do j = k, k
			read (108, *, iostat=iostat) ( omega(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(108)

        open(109, file = Ufile, status='old')
        do j = 1, k-1
			read (109, *, iostat=iostat)
            if (iostat /= 0) exit
		end do
        do j = k, k
			read (109, *, iostat=iostat) ( u(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(109)

        open(110, file = Vfile, status='old')
        do j = 1, k-1
			read (110, *, iostat=iostat)
            if (iostat /= 0) exit
		end do
        do j = k, k
			read (110, *, iostat=iostat) ( v(i), i = 1, n)
            if (iostat /= 0) exit
        end do
        close(110)

        ! read sparse bond data at timestep k
        call read_sparse_restart(Bfile, k)

        ! open(111, file = Bfile, status='old')
        ! do j = 1, (n + 1) * (k - 1)
        !     read (111, *, iostat=iostat)
        !     if (iostat /= 0) exit
        ! end do
        ! do j = (n + 1) * k - n, (n + 1) * k - 1
        !     read (111, *, iostat=iostat) ( bond(i, j - (n + 1) * (k - 1)),  i = 1, n )
        !     if (iostat /= 0) exit
        ! end do
        ! close(111)

	else

		! position of particles: need a one column prepared text file
		write(*,*) ('Reading prepared text files with given path')

		open(103, file = Xfile, status='old')
		open(104, file = Yfile, status='old')
		open(105, file = Rfile, status='old')
		open(106, file = Hfile, status='old')
        open(107, file = Tfile, status='old')
        open(108, file = Ofile, status='old')

		do i = 1, n

			read(103,*) x(i)
			read(104,*) y(i)
			read(105,*) r(i)
			read(106,*) h(i)
            read(107,*) theta(i)
            read(108,*) omega(i)
			
		end do

		do i = 103, 108
			close(i)
		end do

        do i = 1, n
            ! initial velocity
            u(i)      =  0d0
            v(i)      =  0d0
        end do

    end if

    do i = 1, n
        ! initial forces
        tfx(i)    =  0d0
        tfy(i)    =  0d0
        fcx(i)    =  0d0
        fcy(i)    =  0d0
        fbx(i)    =  0d0
        fby(i)    =  0d0
        
        ! initial moments
        m(i)      =  0d0
        mc(i)     =  0d0
        mb(i)     =  0d0

        ! initial stresses
        sigxx(i)  = 0d0
        sigyy(i)  = 0d0
        sigxy(i)  = 0d0
        sigyx(i)  = 0d0

        ! pressure
        tp(i)  = 0d0
        tac(i) = 0d0
        tab(i) = 0d0
        pc(i)  = 0d0
        pb(i)  = 0d0

        ! initial forces and moments in the contact plane
        do j = 1, n
            fcn(j,i)  =  0d0
            fct(j,i)  =  0d0
            fcr(j,i)  =  0d0
            mbb(j,i)  =  0d0
            mcc(j,i)  =  0d0
            fbn(j,i)  =  0d0
            fbt(j,i)  =  0d0
            ! relative angular position is 0 at first
            thetarelc(j,i) = 0d0
            thetarelb(j,i) = 0d0
            ! tangential compression at contact
            deltat(j,i) = 0d0
            ! bond elongation et deflection
            deltanb(j,i) = 0d0
            deltatb(j,i) = 0d0
            ! sheltering coefficient
            hsfa(j, i) = 1d0
            hsfw(j, i) = 1d0
            ! contact area
            ac(j, i) = 0d0
        end do

        ! boundary
        theta_bc1(i)  = 0d0
        theta_bc2(i)  = 0d0
        deltat_bc1(i) = 0d0
        deltat_bc2(i) = 0d0

        !---------------------------------------------------------------
        !           Disks physical properties
        !---------------------------------------------------------------

        call floe_properties(i)

    end do

end subroutine ini_get


subroutine read_sparse_restart(filename, last_tstep)
    
    implicit none

    include "parameter.h"
    include "CB_bond.h"

    character, intent(in) :: filename*32
    integer, intent(in) :: last_tstep

    integer :: i, j, tstep
    integer :: iostat
    integer :: current_tstep

    current_tstep = -1

    open(unit=111, file=filename, status='old', action='read')

    do
        read(111, *, iostat=iostat) tstep, j, i
        if (iostat /= 0) exit
        if (tstep == last_tstep + 1) exit

        ! New timestep detected
        if (tstep /= current_tstep) then
            bond = 0
            current_tstep = tstep
        end if

        bond(j, i) = 1

    end do

    close(111)
end subroutine read_sparse_restart