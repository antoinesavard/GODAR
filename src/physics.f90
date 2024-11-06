subroutine coulomb (j, i, ktc, gamt)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_options.h"

    integer, intent(in) :: i, j
    double precision, intent(in) :: ktc, gamt

	! ensures slipping if fct is too big (kinetic friction)
    ! this is applied both on center of mass and moment
    if ( abs( fct(j,i) ) > friction_coeff * abs( fcn(j,i) ) ) then

        fct(j,i) = - friction_coeff * abs( fcn(j,i) ) * &
                    sign(1d0, velt(j,i))
        fcr(j,i) = fct(j,i)
        
        if ( slipping .eqv. .true. ) then
            deltat(j,i) = (fct(j,i) + gamt * velt(j,i)) / ktc
        end if

    ! static friction is not applied on center of mass, it only
    ! creates a moment
    else

        fcr(j,i) = fct(j,i)

    end if
    
end subroutine coulomb


subroutine coulomb_bc (i, velt_bc, ktc, gamt, deltat_bc)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_options.h" 

    integer, intent(in) :: i
    double precision, intent(in) :: velt_bc
    double precision, intent(in) :: ktc, gamt
    double precision, intent(out) :: deltat_bc

	! ensures slipping if force_t is too big (kinetic friction)
    ! this is applied both on center of mass and moment
    if ( abs( ft_bc(i) ) > friction_coeff * abs( fn_bc(i) ) ) then

        ft_bc(i) = - friction_coeff * abs( fn_bc(i) ) * &
                    sign(1d0, velt_bc )
        fr_bc(i) = ft_bc(i)
        
        if ( slipping .eqv. .true. ) then
            deltat_bc = (ft_bc(i) + gamt * velt_bc) / ktc
        end if

    ! static friction is not applied on center of mass, it only
    ! creates a moment
    else

        fr_bc(i) = ft_bc(i)

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
    
end subroutine rel_pos_vel


subroutine contact_local_to_global (j, i)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    ! update contact force on particle i by particle j
    fcx(i) = fcx(i) - fcn(j,i) * cosa(j,i)
    fcy(i) = fcy(i) - fcn(j,i) * sina(j,i)

    ! update moment on particule i by particule j due to tangent contact 
    mc(i) = mc(i) - r(i) * fct(j,i) - mcc(j,i)

    ! Newton's third law
    ! update contact force on particle j by particle i
    fcx(j) = fcx(j) + fcn(j,i) * cosa(j,i)
    fcy(j) = fcy(j) + fcn(j,i) * sina(j,i)

    ! update moment on particule j by particule i due to tangent contact 
    mc(j) = mc(j) - r(j) * fct(j,i) - mcc(j,i)

end subroutine contact_local_to_global


subroutine bond_local_to_global (j, i)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_bond.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    ! update force on particle i by j due to bond
    fbx(i) = fbx(i) - fbn(j,i) * cosa(j,i) +    &
                        fbt(j,i) * sina(j,i)
    fby(i) = fby(i) - fbn(j,i) * sina(j,i) -    &
                        fbt(j,i) * cosa(j,i)

    ! update moment on particule i by j to to bond
    mb(i) = mb(i) - r(i) * fbt(j,i) - mbb(j, i)

    ! Newton's third law
    ! update force on particle j by i due to bond
    fbx(j) = fbx(j) + fbn(j,i) * cosa(j,i) -    &
                        fbt(j,i) * sina(j,i)
    fby(j) = fby(j) + fbn(j,i) * sina(j,i) +    &
                        fbt(j,i) * cosa(j,i)


    ! update moment on particule j by i due to bond
    mb(j) = mb(j) - r(j) * fbt(j,i) - mbb(j, i)

end subroutine bond_local_to_global


subroutine floe_properties(i)

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_variables.h"

    integer, intent(in) :: i

    ! mass of disk
    mass(i)  =  rhoice * pi * h(i) * r(i) ** 2
    ! freeboard height
    hfa(i)   =  h(i) * (rhowater - rhoice) / rhowater
    ! drag from water height
    hfw(i)   =  h(i) * rhoice / rhowater

end subroutine