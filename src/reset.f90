subroutine reset_contact (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: j, i

    thetarelc(j,i) = 0d0

end subroutine reset_contact


subroutine reset_bond (j, i)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: j, i

    thetarelb(j,i) = 0d0
    deltanb(j,i)   = 0d0
    deltatb(j,i)   = 0d0

end subroutine reset_bond


subroutine reset_forces

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"

    integer :: i

    ! reinitialize force arrays for contact and bonds
    do i = 1, n
        mc(i)    = 0d0
        mb(i)    = 0d0
        fcx(i)   = 0d0
        fcy(i)   = 0d0
        fbx(i)   = 0d0
        fby(i)   = 0d0
        ! mpi reduce
        mc_r(i)  = 0d0
        mb_r(i)  = 0d0
        fcx_r(i) = 0d0
        fcy_r(i) = 0d0
        fbx_r(i) = 0d0
        fby_r(i) = 0d0
        ! and total force arrays
        m(i)     = 0d0
		tfx(i)   = 0d0
        tfy(i)   = 0d0
        ! mpi reduce
        m_r(i)   = 0d0
		tfx_r(i) = 0d0
        tfy_r(i) = 0d0
    end do

    ! and shelter height
    hsfa = 0d0
    hsfw = 0d0

end subroutine reset_forces