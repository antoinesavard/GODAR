subroutine ini_get (restart, expno_str_r, nt_r)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
	include "CB_bond.h"

	integer, intent(in) :: restart, nt_r
	character(2), intent(in) :: expno_str_r

    integer :: i, j, k
    character Tfile*32, Ofile*32, Bfile*32

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

        open(102, file = Xfile, status='old')
		open(103, file = Yfile, status='old')
		open(104, file = Rfile, status='old')
		open(105, file = Hfile, status='old')
		open(106, file = Tfile, status='old')
		open(107, file = Ofile, status='old')
        open(108, file = Bfile, status='old')
		
		k = int(nt_r)

		do j = 1, k-1
			read (102, *)
			read (103, *)
			read (104, *)
			read (105, *)
			read (106, *)
			read (107, *)
		end do

        do j = k, k
			read (102, *) ( x(i), 		i = 1, n)
			read (103, *) ( y(i), 		i = 1, n)
			read (104, *) ( r(i), 		i = 1, n)
			read (105, *) ( h(i), 		i = 1, n)
			read (106, *) ( theta(i), 	i = 1, n)
			read (107, *) ( omega(i), 	i = 1, n)
		end do

        do j = 1, (n + 1) * (k - 1)
            read (108, *)
        end do

        do i = 102, 107
			close(i)
		end do

		do j = k, k
			read (102, *) ( x(i), 		i = 1, n)
			read (103, *) ( y(i), 		i = 1, n)
			read (104, *) ( r(i), 		i = 1, n)
			read (105, *) ( h(i), 		i = 1, n)
			read (106, *) ( theta(i), 	i = 1, n)
			read (107, *) ( omega(i), 	i = 1, n)
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

		do i = 1, n

			read(102,*) x(i)
			read(103,*) y(i)
			read(104,*) r(i)
			read(105,*) h(i)
			
		end do

		do i = 102, 105
			close(i)
		end do
        
    end if

    do i = 1, n
        ! initial velocity
        u(i)      =  0d0
        v(i)      =  0d0
        omega(i)  =  0d0
        
        ! initial particle angle
        theta(i)  =  0d0

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
            ! bond elongation et deflection
            deltanb(j,i) = 0d0
            deltatb(j,i) = 0d0
        end do
    end do

    !-------------------------------------------------------------------
    !           Disks physical properties
    !-------------------------------------------------------------------

    ! mass of disk
    mass  =  rhoice * pi * h * r ** 2
    ! freeboard height
    hfa   =  h * (rhowater - rhoice) / rhowater
    ! drag from water height
    hfw   =  h * rhoice / rhowater

end subroutine ini_get