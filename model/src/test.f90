subroutine normal_forces (side, tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"

    character (*), intent(in) :: side
    integer, intent(in) :: tstep

    double precision :: tau
    integer :: i
    
    ! ~12 seconds is required to resolve the elastic wave travelling 30km
    tau = 50 / dt

    if ( side == "top" ) then
        do i = n, n-9, -1
            if ( tstep < tau ) then
                tfy(i) = tfy(i) + pfn/2 * sin( pi*tstep/tau - pi/2 ) + pfn/2
            else
                tfy(i) = tfy(i) + pfn
            end if
        end do
    end if 
    
    if ( side == "side" ) then
        do i = n, n, -1
            tfx(i) = tfx(i) + pfs * tanh( tstep / tau )
        end do
    end if

end subroutine normal_forces


subroutine gravity

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    double precision :: fg, angle
    integer :: i

    angle = 0d0 * 2d0 * pi / 360d0
    fg = mass(1) * 9.81

    do i = 1, n
        tfx(i) = tfx(i) + fg * sin(angle)
        tfy(i) = tfy(i) - fg * cos(angle)
    end do

end subroutine gravity