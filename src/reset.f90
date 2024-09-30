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


subroutine reset_boundary (i, j)

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    integer, intent(in) :: i, j

    if (j .eq. 1) then

        theta_bc1(i) = 0d0
    
    else if (j .eq. 2) then

        theta_bc2(i) = 0d0

    else    
    
        theta_bc1(i) = 0d0
        theta_bc2(i) = 0d0

    end if
     

end subroutine reset_boundary



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
        sigxx_r(i) = 0d0
        sigyy_r(i) = 0d0
        sigxy_r(i) = 0d0
        sigyx_r(i) = 0d0
        ! and total force arrays
        m(i)     = 0d0
		tfx(i)   = 0d0
        tfy(i)   = 0d0
        ! mpi reduce
        m_r(i)   = 0d0
		tfx_r(i) = 0d0
        tfy_r(i) = 0d0
        tsigxx_r(i) = 0d0
        tsigyy_r(i) = 0d0
        tsigxy_r(i) = 0d0
        tsigyx_r(i) = 0d0
        ! boundary forces
        mc_bc(i) = 0d0
        m_bc(i)  = 0d0
        fn_bc(i) = 0d0
        ft_bc(i) = 0d0
        fx_bc(i) = 0d0
        fy_bc(i) = 0d0
        sigxx_bc(i) = 0d0
        sigyy_bc(i) = 0d0
        sigxy_bc(i) = 0d0
        sigyx_bc(i) = 0d0
        ! stresses
        sigxx(i) = 0d0
        sigyy(i) = 0d0
        sigxy(i) = 0d0
        sigyx(i) = 0d0
    end do

end subroutine reset_forces


subroutine reset_shelter

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    include "CB_mpi.h"

    ! reinitialize sheltering height
    hsfa = 1d0
    hsfw = 1d0

end subroutine reset_shelter