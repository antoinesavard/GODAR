subroutine plastic_contact (j, i, deltat_r, m_redu, hmin, krc)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: j, i
    double precision, intent(in) :: deltat_r, m_redu, hmin
    double precision, intent(out) :: krc

    double precision :: knc, ktc, gamn, gamt

    ! compute the spring constants
    knc    = sigmanc_crit * hmin ** 2 * deltat_r / deltan(j,i)
    ktc    = 6d0 * gc / ec * knc
    krc    = knc * deltat_r ** 2 / 12

    ! compute the dashpots constants
    gamn   = -beta * sqrt( 4d0 * knc * m_redu )
    gamt   = -beta * sqrt( 4d0 * ktc * m_redu )

    ! compute the forces
    fcn(j,i) = knc * deltan(j,i) - gamn * veln(j,i)
    fct(j,i) = ktc * deltat(j,i) - gamt * velt(j,i)

    call update_shape (j, i, deltat_r)

end subroutine plastic_contact


subroutine update_shape (j, i, deltat_r)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: j, i
    double precision, intent(in) :: deltat_r

    double precision :: hmin, dh, Vol, Area

    hmin = min(h(i), h(j))

    Area = r(i) ** 2 * acos((dist(j,i) ** 2 - r(j) ** 2 + r(i) ** 2) &
            / (2 * dist(j,i)) / r(i)) - (dist(j,i) ** 2 - r(j) ** 2  &
            + r(i) ** 2) / (2 * dist(j,i)) * deltat_r / 2d0 +        &
        r(j) ** 2 * acos((dist(j,i) ** 2 - r(i) ** 2 + r(j) ** 2)    &
            / (2 * dist(j,i)) / r(j)) - (dist(j,i) ** 2 - r(i) ** 2  &
            + r(j) ** 2) / (2 * dist(j,i)) * deltat_r / 2d0

    Vol = Area * hmin

    if ( hmin .eq. h(i) ) then

        dh = Vol / ( pi * r(i) ** 2d0 )

        r(i) = r(i) * sqrt(h(i) / (h(i) + dh) )

        h(i) = h(i) + dh

        call floe_properties(i)

    else if ( hmin .eq. h(j) ) then

        dh = Vol / ( pi * r(j) ** 2d0 )

        r(j) = r(j) * sqrt(h(j) / (h(j) + dh) )

        h(j) = h(j) + dh

        call floe_properties(j)

    end if

end subroutine update_shape


subroutine plastic_contact_bc (i, deltan_bc, deltat_bc, deltat_r_bc, krc, dir)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, dir
    double precision, intent(in) :: deltan_bc, deltat_bc, deltat_r_bc
    double precision, intent(out) :: krc

    double precision :: knc, ktc, gamn, gamt

    ! compute the spring constants
    knc    = sigmanc_crit * h(i) ** 2 * deltat_r_bc / deltan_bc
    ktc    = 6d0 * gc / ec * knc
    krc    = knc * deltat_r_bc ** 2 / 12

    ! compute the dashpots constants
    gamn   = -beta * sqrt( 4d0 * knc * m(i) )
    gamt   = -beta * sqrt( 4d0 * ktc * m(i) )

    ! compute the forces
    fn_bc(i) = knc * deltan_bc &
                - gamn * ( dir * u(i) + (1 - dir) * v(i) )

    ft_bc(i) = ktc * deltat_bc &
                - gamt * ( (1 - dir) * u(i) + dir * v(i) )

    call update_shape_bc (i, deltat_r_bc)

end subroutine plastic_contact_bc


subroutine update_shape_bc (i, deltat_r_bc)

    implicit none 

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i
    double precision, intent(in) :: deltat_r_bc

    double precision :: dh, Vol, Area

    Area = r(i) ** 2 * asin(deltat_r_bc / 2 / r(i))

    Vol = Area * h(i)

    dh = Vol / ( pi * r(i) ** 2d0 )

    r(i) = r(i) * sqrt(h(i) / (h(i) + dh) )

    h(i) = h(i) + dh

    call floe_properties(i)

end subroutine update_shape_bc