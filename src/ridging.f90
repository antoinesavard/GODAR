subroutine plastic_contact (j, i, m_redu, hmin, ktc, krc, gamn, gamt, gamr)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: j, i
    double precision, intent(in) :: m_redu, hmin
    double precision, intent(out) :: ktc, krc, gamt, gamr
    double precision, intent(inout) :: gamn

    double precision :: knc

    ! compute the spring constants
    knc    = (sigmanc_crit * hmin ** 2 * delt_ridge(j,i) &
                ) / deltan(j,i) !+ gamn * veln(j,i)
    ktc    = 6d0 * gc / ec * knc
    krc    = knc * delt_ridge(j,i) ** 2 / 12

    ! compute the dashpots constants
    gamn   = -beta * sqrt( 4d0 * knc * m_redu )
    gamt   = -2d0 * beta * sqrt( 2d0/3d0 * ktc * m_redu )
    gamr   = gamn * delt_ridge(j,i) ** 2 / 12

    ! compute the forces
    fcn(j,i) = min(fcn(j,i), sigmanc_crit * hmin**2 * delt_ridge(j,i)) !knc * deltan(j,i) - gamn * veln(j,i) !
    fct(j,i) = ktc * deltat(j,i) - gamt * velt(j,i)

    call update_shape (j, i)

end subroutine plastic_contact


subroutine update_shape (j, i)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: j, i

    double precision :: hmin, dh, Vol, Area
    double precision :: delta_ij, delta_ji

    hmin = min(h(i), h(j))
    delta_ij = (dist(j,i) ** 2 - r(j) ** 2 + r(i) ** 2) / &
                (2 * dist(j,i))
    delta_ji = (dist(j,i) ** 2 - r(i) ** 2 + r(j) ** 2) / &
                (2 * dist(j,i))

    Area = r(i) ** 2 * acos(delta_ij / r(i)) - delta_ij * &
            delt_ridge(j, i) / 2d0 + r(j) ** 2 * acos(delta_ji / r(j)) - delta_ji * delt_ridge(j,i) / 2d0

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

    if ( hmin .eq. h(i) ) then

        dh = Vol / ( pi * r(i) ** 2d0 )

        r(i) = r(i) * sqrt(h(i) / (h(i) + dh) )

        h(i) = h(i) + dh

        ! recalculate floe freeboard
        call floe_properties(i)

    else if ( hmin .eq. h(j) ) then

        dh = Vol / ( pi * r(j) ** 2d0 )

        r(j) = r(j) * sqrt(h(j) / (h(j) + dh) )

        h(j) = h(j) + dh

        ! recalculate floe freeboard
        call floe_properties(j)

    end if

    ! update the new tangent overlap for pressure
    delt_ridge(j,i) = 2 * sqrt( r(i) ** 2 - ( (dist(j,i) ** 2 - &
                    r(j) ** 2 + r(i) ** 2) / (2 * dist(j,i)) ) ** 2 )

end subroutine update_shape


subroutine plastic_contact_bc (i, veln_bc, velt_bc, deltan_bc, deltat_bc, ktc, krc, gamn, gamt, gamr)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

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
    gamn   = -beta * sqrt( 4d0 * knc * m(i) / 1d0 )
    gamt   = -2d0 * beta * sqrt( 2d0/3d0 * ktc * m(i) / 1d0 )
    gamr   = gamn * delt_ridge_bc(i) ** 2 / 12

    ! compute the forces
    fn_bc(i) = min(fn_bc(i), sigmanc_crit * h(i)**2 * delt_ridge_bc(i)) !knc * deltan_bc - gamn * veln_bc !
    ft_bc(i) = ktc * deltat_bc - gamt * velt_bc

    call update_shape_bc (i, deltan_bc)

end subroutine plastic_contact_bc


subroutine update_shape_bc (i, deltan_bc)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    double precision, intent(in) :: deltan_bc

    double precision :: dh, Vol, Area

    Area = r(i) ** 2 * asin(delt_ridge_bc(i) / 2 / r(i)) &
            - delt_ridge_bc(i) * ( r(i) - deltan_bc ) / 2

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

    dh = Vol / ( pi * r(i) ** 2d0 )

    r(i) = r(i) * sqrt(h(i) / (h(i) + dh) )

    h(i) = h(i) + dh

    ! recalculate the floe freeboard etc.
    call floe_properties(i)

    ! update the new tangent overlap for pressure
    delt_ridge_bc(i) = 2 * sqrt( r(i) ** 2 - ( r(i) - deltan_bc ) ** 2 )

end subroutine update_shape_bc