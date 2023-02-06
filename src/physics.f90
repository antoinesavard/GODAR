subroutine moment (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    m(i) = r(i) * fct(i,j) + 	&	! tangential forces
			mw(i) + ma(i) !+ 	&	! forcing/drag
			!mbb(i) + (r(i) - deltan(i, j)) * fbt(i, j) ! bonds
    
end subroutine moment


subroutine coulomb (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

	! ensures slipping if force_t is too big
    if ( fct(i,j) > friction_coeff * fcn(i,j) ) then

        fct(i,j) = 0

    end if
    
end subroutine coulomb


subroutine rel_pos_vel (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    ! Components of unit vector ei=(cosa,sina) are:
	cosa(i,j) = ( x(j) - x(i) ) /                 	 &
				( sqrt(                              &
						( x(i) - x(j) ) ** 2 +       &
						( y(i) - y(j) ) ** 2         &
						)                            &
				)
	
	sina(i,j) = ( y(j) - y(i) ) /                 	 &
				( sqrt(                              &
						( x(i) - x(j) ) ** 2 +       &
						( y(i) - y(j) ) ** 2         &
						)                            &
				)

	! relative angular velocity
	omegarel(i,j) = ( omega(j) - omega(i) )

	! Normal components of the relative velocities:
	veln(i,j) = ( u(j) - u(i) ) * cosa(i,j) +     &
				( v(j) - v(i) ) * sina(i,j)

	! Tangential components of the relative velocities:
	velt(i,j) = ( u(j) - u(i) ) * sina(i,j) +   &
				( v(j) - v(i) ) * cosa(i,j) -   &
				( omega(i) * r(i) + omega(j) * r(j) )

	! normal overlap (displacement) deltan >=0
	deltan(i,j)  =  r(i) + r(j) -          & 
					sqrt(                  & 
					( x(i) - x(j) ) ** 2 + & 
					( y(i) - y(j) ) ** 2   & 
					)
    
end subroutine rel_pos_vel