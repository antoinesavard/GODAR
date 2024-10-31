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
        Bfile = 'output/bond.' // trim(adjustl(expno_str_r))
		
		k = int(nt_r)

        open(102, file = Xfile, status='old')
		do j = 1, k-1
			read (102, *)
        end do
        do j = k, k
			read (102, *) ( x(i), 		i = 1, n)
        end do
        close(102)

        open(103, file = Yfile, status='old')
        do j = 1, k-1
			read (103, *)
        end do
        do j = k, k
			read (103, *) ( y(i), 		i = 1, n)
        end do
        close(103)

        open(104, file = Rfile, status='old')
        do j = 1, k-1
			read (104, *)
        end do
        do j = k, k
			read (104, *) ( r(i), 		i = 1, n)
        end do
        close(104)

        open(105, file = Hfile, status='old')
        do j = 1, k-1
			read (105, *)
        end do
        do j = k, k
			read (105, *) ( h(i), 		i = 1, n)
        end do
        close(105)

        open(106, file = Tfile, status='old')
        do j = 1, k-1
			read (106, *)
        end do
        do j = k, k
			read (106, *) ( theta(i),	i = 1, n)
        end do
        close(106)

        open(107, file = Ofile, status='old')
        do j = 1, k-1
			read (107, *)
		end do
        do j = k, k
			read (107, *) ( omega(i),	i = 1, n)
        end do
        close(107)

        open(108, file = Bfile, status='old')
        do j = 1, (n + 1) * (k - 1)
            read (108, *)
        end do
        do j = (n + 1) * k - n, (n + 1) * k - 1
            read (108, *) ( bond(i, j - (n + 1) * (k - 1)),  i = 1, n )
        end do
        close(108)

	else

		! position of particles: need a one column prepared text file
		write(*,*) ('Reading prepared text files with given path')

		open(102, file = Xfile, status='old')
		open(103, file = Yfile, status='old')
		open(104, file = Rfile, status='old')
		open(105, file = Hfile, status='old')
        open(106, file = Tfile, status='old')
        open(107, file = Ofile, status='old')

		do i = 1, n

			read(102,*) x(i)
			read(103,*) y(i)
			read(104,*) r(i)
			read(105,*) h(i)
            read(106,*) theta(i)
            read(107,*) omega(i)
			
		end do

		do i = 102, 107
			close(i)
		end do
        
    end if

    do i = 1, n
        ! initial velocity
        u(i)      =  0d0
        v(i)      =  0d0
        omega(i)  =  0d0

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