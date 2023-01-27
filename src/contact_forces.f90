subroutine contact_forces (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    double precision :: deltat(n,n), mc(n,n), rc(n,n), hmin(n,n)
    double precision :: fit
    double precision :: kn, kt, gamn, gamt

    deltat(i,j)  =  velt(i,j) * dt

    mc(i,j)      =  mass(i) * mass(j) / ( mass(i) + mass(j) )

    rc(i,j)      =  r(i) * r(j) / ( r(i) + r(j) )

    hmin(i,j)    =  min(h(i), h(j))

    kn     = pi * ec * hmin(i,j)  *             &
                fit( deltan(i,j) * rc(i,j) /    &
                ( 2 * hmin(i,j) ** 2 ) )

    kt     = 6 * gc / ec * kn

    gamn   = - beta * sqrt( 5 * kn * mc(i,j) )

    gamt   = - 2 * beta * sqrt( 5 * gc / ec * kn * mc(i,j) )

    ! compute the normal/tangent force

    fcn(i,j) = kn * deltan(i,j) - gamn * veln(i,j)

    fct(i,j) = kt * deltat(i,j) - gamt * velt(i,j)

    call coulomb (i,j)

    call moment (i,j)

end subroutine contact_forces

subroutine coulomb (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    if ( fct(i,j) > friction_coeff * fcn(i,j) ) then

        fct(i,j) = 0

    end if
    
end subroutine coulomb

subroutine moment (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    m(i) = r(i) * fct(i,j) + mw(i)
    
end subroutine moment

double precision function fit (xi)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    double precision, intent(in)  :: xi
    double precision :: p1, p2, p3, q1, q2

    p1 = 9117d-4
    p2 = 2722d-4
    p3 = 3324d-6
    q1 = 1524d-3
    q2 = 3159d-5

    fit = ( p1 * xi ** 2 + p2 * xi + p3 ) / ( xi ** 2 + q1 * xi + q2 )

end function fit
