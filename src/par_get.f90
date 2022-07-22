subroutine get_default
    
    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

! set parameter for the run

    k       =  1.5d0        ! normal stiffness
    ks      =  1.5d0        ! shear stiffness
    dt      =  5d-2         ! length of time step
    t       =  150 * dt     ! length of time
    st      =  1            ! divide tot num of timestep to st
    eta     =  1.5d0        ! Normal damping coefficient
    eta2    =  1.5d0        ! shear damping coefficient

end subroutine get_default