subroutine ini_get (restart, expno_r)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

	integer, intent(in) :: restart, expno_r

    integer :: i
    character Xfile*32, Yfile*32, Rfile*32, Hfile*32

	! load restart files
	if ( restart .eq. 1 ) then



	else

		! position of particles: need a one column prepared text file
		
		Xfile = 'src/X.dat'
		Yfile = 'src/Y.dat'
		Rfile = 'src/R.dat'
		Hfile = 'src/H.dat'

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

		do i = 1, n
			! initial velocity
			u(i)      =  0d0
			v(i)      =  0d0
			omega(i)  =  0d0
			
			! initial particle angle
			theta(i)  =  0d0 
		end do

	end if

end subroutine ini_get