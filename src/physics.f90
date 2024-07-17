subroutine coulomb (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

	! ensures slipping if force_t is too big
    if ( abs( fct(j,i) ) > friction_coeff * abs( fcn(j,i) ) ) then

        fct(j,i) = - friction_coeff * fcn(j,i) * sign(1d0, velt(j,i))

    end if
    
end subroutine coulomb


subroutine coulomb_bc (i, dir)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    integer, intent(in) :: dir

	! ensures slipping if force_t is too big
    if ( abs( ft_bc(i) ) > friction_coeff * abs( fn_bc(i) ) ) then

        ft_bc(i) = - friction_coeff * fn_bc(i) * &
                    sign(1d0, (1 - dir) * u(i) + dir * v(i) )

    end if
    
end subroutine coulomb_bc


subroutine rel_pos_vel (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    ! distance between particles are
    dist(j,i) = ( sqrt(                              &
						( x(i) - x(j) ) ** 2 +       &
						( y(i) - y(j) ) ** 2         &
						)                            &
				)

    ! Components of unit vector ei=(cosa,sina) are:
	cosa(j,i) = ( x(j) - x(i) ) / dist(j,i)
	
	sina(j,i) = ( y(j) - y(i) ) / dist(j,i)

	! relative angular velocity
	omegarel(j,i) = ( omega(j) - omega(i) )

	! Normal components of the relative velocities:
	veln(j,i) = ( u(j) - u(i) ) * cosa(j,i) +     &
				( v(j) - v(i) ) * sina(j,i)

	! Tangential components of the relative velocities:
	velt(j,i) = -( u(j) - u(i) ) * sina(j,i) +   &
				( v(j) - v(i) ) * cosa(j,i) -    &
				( omega(i) * r(i) + omega(j) * r(j) )

	! normal overlap (displacement) deltan >=0
	deltan(j,i)  =  r(i) + r(j) - dist(j,i)
    
    ! tangent overlap
    deltat(j,i)  =  0d0
    
end subroutine rel_pos_vel