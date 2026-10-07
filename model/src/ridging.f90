subroutine plastic_contact (j, i, p, m_redu, hmin, ktc, krc, gamn, gamt, gamr)

    use parameters
    use variables
    use const
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: j, i
    type(pair_t), intent(inout) :: p
    double precision, intent(in) :: m_redu, hmin
    double precision, intent(out) :: ktc, krc, gamt, gamr
    double precision, intent(inout) :: gamn

    double precision :: knc

    ! compute the spring constants
    knc    = (sigmanc_crit * hmin ** 2 * p%delt_ridge &
                ) / p%deltan !+ gamn * veln(j,i)
    ktc    = 6d0 * gc / ec * knc
    krc    = knc * p%delt_ridge ** 2 / 12

    ! compute the dashpots constants
    gamn   = 2d0 * beta * sqrt( knc * m_redu )
    gamt   = 2d0 * beta * sqrt( 2d0/3d0 * ktc * m_redu )
    gamr   = gamn * p%delt_ridge ** 2 / 12

    ! compute the forces
    p%fcn = min(p%fcn, sigmanc_crit * hmin**2 * p%delt_ridge) !knc * deltan(j,i) - gamn * veln(j,i) !
    p%fct = ktc * deltat(j,i) - gamt * p%velt

    call update_shape (j, i, p)

end subroutine plastic_contact


subroutine update_shape (j, i, p)

    use parameters
    use variables
    use const
    use bonds
    use pairs, only: pair_t

    implicit none


    integer, intent(in) :: j, i
    type(pair_t), intent(in) :: p

    double precision :: hmin, Vol, Area
    double precision :: delta_ij, delta_ji, arg_ij, arg_ji
    integer :: k

    hmin = min(h(i), h(j))
    delta_ij = (p%dist ** 2 - r(j) ** 2 + r(i) ** 2) / &
                (2 * p%dist)
    delta_ji = (p%dist ** 2 - r(i) ** 2 + r(j) ** 2) / &
                (2 * p%dist)

    ! clamp arguments to [-1,1] to prevent NaN from -ffast-math
    arg_ij = max(min(delta_ij / r(i), 1d0), -1d0)
    arg_ji = max(min(delta_ji / r(j), 1d0), -1d0)

    Area = r(i) ** 2 * acos(arg_ij) - delta_ij * &
            p%delt_ridge / 2d0 + r(j) ** 2 * acos(arg_ji) - delta_ji * p%delt_ridge / 2d0

    ! take the max value to make sure we don't get negatives
    ! this is purely a compilation/numerical trick
    ! not doing this can result in particles exploding
    ! because when using aggressive optimization flags at
    ! compilation, and since asin/acos are extra sensitive on the
    ! distances, the area computation can become negative, even
    ! though it should not happen, in the substraction because 
    ! the distances are all super small. These negative volumes
    ! increase the radius -> which increase the contact forces ->
    ! leads to explosions. 
    Vol = max(Area * hmin / 2d0, 0d0)

    ! the thinner floe converts the overlap volume into thickness. The
    ! shape change is applied after the force loop (apply_ridging), so
    ! that every pair sees the same shapes during the step and no
    ! thread or rank modifies floes that belong to another one
    if ( hmin .eq. h(i) ) then
        k = i
    else
        k = j
    end if

    !$omp atomic update
    dvol(k) = dvol(k) + Vol

end subroutine update_shape


subroutine plastic_contact_bc (i, veln_bc, velt_bc, deltan_bc, deltat_bc, ktc, krc, gamn, gamt, gamr)

    use parameters
    use variables
    use const

    implicit none


    integer, intent(in) :: i
    double precision, intent(in) :: veln_bc, velt_bc
    double precision, intent(in) :: deltan_bc, deltat_bc
    double precision, intent(out) :: ktc, krc, gamt, gamr
    double precision, intent(inout) :: gamn

    double precision :: knc

    ! compute the spring constants
    knc    = (sigmanc_crit * h(i) ** 2 * delt_ridge_bc(i) &
                ) / deltan_bc ! + gamn * veln_bc
    ktc    = 6d0 * gc / ec * knc
    krc    = knc * delt_ridge_bc(i) ** 2 / 12

    ! compute the dashpots constants
    ! note the 1/2 factor in gamn and gamt: this is a choice.
    ! same reason as above for consistency
    gamn   = 2d0 * beta * sqrt( knc * mass(i) / 1d0 )
    gamt   = 2d0 * beta * sqrt( 2d0/3d0 * ktc * mass(i) / 1d0 )
    gamr   = gamn * delt_ridge_bc(i) ** 2 / 12

    ! compute the forces
    fn_bc(i) = min(fn_bc(i), sigmanc_crit * h(i)**2 * delt_ridge_bc(i)) !knc * deltan_bc - gamn * veln_bc !
    ft_bc(i) = ktc * deltat_bc - gamt * velt_bc

    call update_shape_bc (i, deltan_bc)

end subroutine plastic_contact_bc


subroutine update_shape_bc (i, deltan_bc)

    use parameters
    use variables
    use const

    implicit none 


    integer, intent(in) :: i
    double precision, intent(in) :: deltan_bc

    double precision :: Vol, Area
    double precision :: arg

    ! clamp argument to [-1,1] to prevent NaN from -ffast-math
    arg = max(min(delt_ridge_bc(i) / 2d0 / r(i), 1d0), -1d0)

    Area = r(i) ** 2 * asin(arg) - delt_ridge_bc(i) * ( r(i) - deltan_bc ) / 2

    ! take the max value to make sure we don't get negatives
    ! this is purely a compilation/numerical trick
    ! not doing this can result in particles exploding
    ! because when using aggressive optimization flags at
    ! compilation, and since asin/acos are extra sensitive on the
    ! distances, the area computation can become negative, even
    ! though it should not happen, in the substraction because 
    ! the distances are all super small. These negative volumes
    ! increase the radius -> which increase the contact forces ->
    ! leads to explosions. 
    Vol = max(Area * h(i), 0d0)

    ! the shape change is applied after the force loop (apply_ridging)
    !$omp atomic update
    dvol(i) = dvol(i) + Vol

end subroutine update_shape_bc


subroutine apply_ridging

    ! Converts the overlap volume accumulated during the step (dvol)
    ! into thickness. Mass is conserved: the thickness increase dh is
    ! compensated by a decrease in radius, so that r**2 * h is constant.
    ! Called after the force loop on every rank, with dvol summed over
    ! all ranks (broadcast_shape), so the shapes stay identical
    ! everywhere.

    use parameters
    use variables
    use const

    implicit none


    integer :: i
    double precision :: dh

    do i = 1, n
        if ( dvol(i) .gt. 0d0 ) then

            dh = dvol(i) / ( pi * r(i) ** 2d0 )

            r(i) = r(i) * sqrt(h(i) / (h(i) + dh) )

            h(i) = h(i) + dh

            ! recalculate floe freeboard and inertia
            call floe_properties(i)

        end if
    end do

end subroutine apply_ridging
