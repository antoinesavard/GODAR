subroutine contact_forces

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i, j

    ! set values to zero at each timestep
    wafx  =  0d0
    wafy  =  0d0
    dafx  =  0d0
    dafy  =  0d0

    do i = 1, n
        do j = 1, n

            ! sum of global components of normal incremental force at contact point of each particle
            fsx(j,i)  =  fsx(j,i) + sx(j,i)
            fsy(j,i)  =  fsy(j,i) + sy(j,i)

            ! resultant shear spring force (damping included)
            fw(j,i)   =  fw(j,i) - w(j,i) - l(j,i)

            ! sum of global components of shear incremental force at contact point of each particle
            fwx(j,i)  =  fwx(j,i) + wx(j,i)
            fwy(j,i)  =  fwy(j,i) + wy(j,i)

            ! sum of global components of normal damping forces at contact point of each particle
            dafx(i)   =  dafx(i) + dx(j,i)
            dafy(i)   =  dafy(i) + dy(j,i)

            ! sum of global components of shear damping forces at contact point of each particle
            wafx(i)   =  wafx(i) + lx(j,i)
            wafy(i)   =  wafy(i) + ly(j,i)

        end do
    end do

    ! resultant force on each particle
    f = 0

    do i = 1, N
        do j = 1, N

            f(i,1) = f(i,1) + fsx(j,i) + fwx(j,i)
            f(i,2) = f(i,2) + fsy(j,i) + fwy(j,i)

        end do
    end do

    ! Resultant moment acting on each particle
    m = 0

    do i = 1, N
        do j = 1, N

            m(i) = m(i) + fw(j,i) * r(i)
            
        end do
    end do

end subroutine contact_forces


subroutine contact_forces2 (i, j, veln, velt)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    double precision, intent(in) :: veln, velt
    double precision :: deltan(n,n), deltat(n,n), mc(n,n), rc(n,n), hmin(n,n)
    double precision :: fit
    double precision :: kn, kt, gamn, gamt

    deltan(i,j)  =  r(i) + r(j) - sqrt(     & ! normal overlap
                    ( x(i) - x(j) ) ** 2 +  & ! (displacement)
                    ( y(i) - y(j) ) ** 2 )

    deltat(i,j)  =  velt * dt

    mc(i,j)      =  mass(i) * mass(j) / ( mass(i) + mass(j) )

    rc(i,j)      =  r(i) * r(j) / ( r(i) + r(j) )

    hmin(i,j)    =  min(h(i), h(j))

    kn     = pi * ec * hmin(i,j)  *             &
                fit( deltan(i,j) * rc(i,j) /     &
                ( 2 * hmin(i,j) ** 2 ) )

    kt     = 6 * gc / ec * kn

    gamn   = - beta * sqrt( 5 * kn * mc(i,j) )

    gamt   = - 2 * beta * sqrt( 5 * gc / ec * kn * mc(i,j) )

    ! compute the normal/tangent force

    fcn(i,j) = kn * deltan(i,j) - gamn * veln

    fct(i,j) = kt * deltat(i,j) - gamt * velt

    call coulomb (i,j)

    call moment (i,j)

end subroutine contact_forces2

subroutine coulomb (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    if ( fct(i,j) > friction_coeff * fcn(i,j) ) then

        fct(i,j) = friction_coeff * fcn(i,j)

    end if
    
end subroutine coulomb

subroutine moment (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    m(i) = r(i) * fct(i,j)
    
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
