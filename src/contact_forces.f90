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

subroutine viscous_terms

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer :: i

    ! Total springs-dashpots forces at each contact point
    do i = 1, N

        tfx(i)  =  wafx(i) + dafx(i) + f(i,1)
        tfy(i)  =  wafy(i) + dafy(i) + f(i,2)
        u(i)    =  u(i) + ( tfx(i) / mass(i) ) * dt
        v(i)    =  v(i) + ( tfy(i) / mass(i) ) * dt
        ome(i)  =  ome(i) + dt * m(i) / ( 4d-1 * mass(i) * r(i) ** 2 ) 

    end do

end subroutine

subroutine forcing (tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: tstep

    if (tstep == 1) then

        u(1) = 5d0
        u(2) = -5d0

    end if

end subroutine forcing