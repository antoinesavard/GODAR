subroutine coulomb (j, i, p, ktc, gamt)

    use parameters
    use variables
    use const
    use options
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: i, j
    type(pair_t), intent(inout) :: p
    double precision, intent(in) :: ktc, gamt

	! ensures slipping if fct is too big (kinetic friction)
    ! this is applied both on center of mass and moment
    if ( abs( p%fct ) > friction_coeff * abs( p%fcn ) ) then

        p%fct = - friction_coeff * abs( p%fcn ) * &
                    sign(1d0, p%velt)
        p%fcr = p%fct

        if ( slipping .eqv. .true. ) then
            deltat(j,i) = (p%fct + gamt * p%velt) / ktc
        end if

    ! static friction is not applied on center of mass, it only
    ! creates a moment
    else

        p%fcr = p%fct

    end if
    
end subroutine coulomb


subroutine coulomb_bc (i, velt_bc, ktc, gamt, deltat_bc)

    use parameters
    use variables
    use const
    use options

    implicit none


    integer, intent(in) :: i
    double precision, intent(in) :: velt_bc
    double precision, intent(in) :: ktc, gamt
    double precision, intent(inout) :: deltat_bc

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


subroutine rel_pos_vel (j, i, p)

    use parameters
    use variables
    use const
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: i, j
    type(pair_t), intent(inout) :: p

    ! distance between particles are
    p%dist = ( sqrt(                                 &
						( x(i) - x(j) ) ** 2 +       &
						( y(i) - y(j) ) ** 2         &
						)                            &
				)

    ! Components of unit vector ei=(cosa,sina) are:
	p%cosa = ( x(j) - x(i) ) / p%dist

	p%sina = ( y(j) - y(i) ) / p%dist

	! relative angular velocity
	p%omegarel = ( omega(j) - omega(i) )

	! Normal components of the relative velocities:
	p%veln = ( u(j) - u(i) ) * p%cosa +     &
				( v(j) - v(i) ) * p%sina

	! Center components of the relative velocities:
    p%veltb = -( u(j) - u(i) ) * p%sina +   &
                 ( v(j) - v(i) ) * p%cosa

    ! Tangential components of the relative velocities:
    p%velt = p%veltb - ( omega(i) * r(i) + omega(j) * r(j) )

	! normal overlap (displacement) deltan >=0
	p%deltan  =  r(i) + r(j) - p%dist
    
end subroutine rel_pos_vel


subroutine contact_local_to_global (j, i, p)

    use parameters
    use variables
    use const
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: i, j
    type(pair_t), intent(in) :: p

    ! update contact force on particle i by particle j
    fcx(i) = fcx(i) - p%fcn * p%cosa
    fcy(i) = fcy(i) - p%fcn * p%sina

    ! update moment on particule i by particule j due to tangent contact
    mc(i) = mc(i) - r(i) * p%fct - p%mcc

    ! Newton's third law
    ! update contact force on particle j by particle i
    fcx(j) = fcx(j) + p%fcn * p%cosa
    fcy(j) = fcy(j) + p%fcn * p%sina

    ! update moment on particule j by particule i due to tangent contact
    mc(j) = mc(j) - r(j) * p%fct - p%mcc

end subroutine contact_local_to_global


subroutine bond_local_to_global (j, i, p)

    use parameters
    use variables
    use bonds
    use const
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: i, j
    type(pair_t), intent(in) :: p

    ! update force on particle i by j due to bond
    fbx(i) = fbx(i) - p%fbn * p%cosa +    &
                        p%fbt * p%sina
    fby(i) = fby(i) - p%fbn * p%sina -    &
                        p%fbt * p%cosa

    ! update moment on particule i by j to to bond
    mb(i) = mb(i) - r(i) * p%fbt - p%mbb_ji

    ! Newton's third law
    ! update force on particle j by i due to bond
    fbx(j) = fbx(j) + p%fbn * p%cosa -    &
                        p%fbt * p%sina
    fby(j) = fby(j) + p%fbn * p%sina +    &
                        p%fbt * p%cosa


    ! update moment on particule j by i due to bond
    mb(j) = mb(j) - r(j) * p%fbt - p%mbb_ji

end subroutine bond_local_to_global


subroutine floe_properties(i)

    use parameters
    use const
    use variables

    implicit none


    integer, intent(in) :: i

    ! mass of disk
    mass(i)  =  rhoice * pi * h(i) * r(i) ** 2
    ! moment of inertia of disk
    inertia(i) = 0.5 * mass(i) * r(i) ** 2
    ! freeboard height
    hfa(i)   =  h(i) * (rhowater - rhoice) / rhowater
    ! drag from water height
    hfw(i)   =  h(i) * rhoice / rhowater

end subroutine floe_properties
