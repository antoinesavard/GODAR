subroutine normal_forces (side, tstep)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"

    character (*), intent(in) :: side
    integer, intent(in) :: tstep

    double precision :: ftmp, tau
    integer :: i
    
    ! ~12 seconds is required to resolve the elastic wave travelling 30km
    tau = 50 / dt

    ! forces on the right of the assembly
    ! forces on the right of the assembly
    if ( side == "right" ) then
        ftmp = 0d0!maxval(tfx(n-29:n))
        if (tstep*dt < 50d-2) then
!            if (tstep*dt > 25d-2) then
                !do i = n, n - 2, -1!n, n - 29 , -1
                    tfx(1) = tfx(1) + ftmp + pfn! * tanh( tstep / tau )
                !end do
!            end if       
        end if
     end if    

    ! no forces on the left of assembly (fixed)
    ! if ( side == "left" ) then
    !     do i = n - 30, n - 59 , -1
    !         tfy(i) = 0d0
    !         tfx(i) = 0d0
    !     end do
    ! end if


end subroutine normal_forces