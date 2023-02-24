subroutine ini_get (restart, expno_r)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

	integer, intent(in) :: restart
	character(2), intent(in) :: expno_r
	double precision :: xtemp (1000, n)

    integer :: i, j, k
    character Xfile*32, Yfile*32, Rfile*32, Hfile*32, Tfile*32, Ofile*32

	! load restart files
	if ( restart .eq. 1 ) then

		Xfile = 'output/x.' // trim(adjustl(expno_r))
		Yfile = 'output/y.' // trim(adjustl(expno_r))
		Rfile = 'output/r.' // trim(adjustl(expno_r))
		Hfile = 'output/h.' // trim(adjustl(expno_r))
		Tfile = 'output/theta.' // trim(adjustl(expno_r))
		Ofile = 'output/omega.' // trim(adjustl(expno_r))

        open(102, file = Xfile, status='old')
		open(103, file = Yfile, status='old')
		open(104, file = Rfile, status='old')
		open(105, file = Hfile, status='old')
		open(106, file = Tfile, status='old')
		open(107, file = Ofile, status='old')
		
		k = 1000

		do j = 1, k-1
			read (102, *)
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
		
		do i = 102, 107
			close(i)
		end do

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