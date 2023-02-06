subroutine stepper (tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j
    integer, intent(in) :: tstep

    ! reinitialize force arrays
    do i = 1, n
        do j = 1, n
            tfx = 0d0
            tfy = 0d0
            fcn = 0d0
            fct = 0d0
        end do
    end do
    
    ! put yourself in the referential of the ith particle
	! loop through all j particles and compute interactions
    do i = 1, n - 1 
        do j = i + 1, n
    
			! compute relative position and velocity
            call rel_pos_vel (i, j)

			! bond initialization
			if ( tstep .eq. 1 ) then
				if ( -deltan(i, j) .leq. 2 * r(i)) ! can be fancier
					bond (i, j) = 1
				end if
			end if

            ! verify if two particles are colliding
            if ( deltan(i,j) .ge. 0 ) then

                call contact_forces (i, j)

				! call bond_creation (i, j) to implement

				! update force on particle i by particle j
                tfx(i) = tfx(i) - fcn(i,j) * cosa(i,j)
                tfy(i) = tfy(i) - fcn(i,j) * sina(i,j)

            end if

			! compute forces from bonds between particle i and j
			if bond (i, j) .eq. 1 then

				call bond_forces (i, j)
				call bond_breaking (i, j)

			end if

			! compute forcing of air and water drag, and coriolis
            call forcing(i, j)
			call coriolis(i, j)

			! call moments at the end because you need all previous 
			! tangential forces to compute this (M = r x F_t)
            call moment(i, j)

			! sum all forces together
            tfx(i) = tfx(i) + fax(i) + fwx(i)
            tfy(i) = tfy(i) + fay(i) + fwy(i)

        end do
    end do

	! integration in time
    call velocity
    call euler

end subroutine stepper