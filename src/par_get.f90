subroutine get_default
    
    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

!-----------------------------------------------------------------------
! set parameter for the run
!-----------------------------------------------------------------------

    k        =  1.5d7              ! normal stiffness
    ks       =  1.5d7              ! shear stiffness
    dt       =  5d-3               ! length of time step
    t        =  1000 * dt          ! length of time
    st       =  1                  ! divide tot num of timestep to st
    eta      =  1.5d0              ! Normal damping coefficient
    eta2     =  1.5d0              ! shear damping coefficient

    pi       =  4d0 * datan(1d0)   ! pi

    Cpair    =  1d03               ! Specific heat of air.
    Cpwater  =  4d03               ! Specific heat of water

    rhoair   =  1.3d0              ! air density [kg/m3]
    rhoice   =  9d02               ! ice density [kg/m3]
    rhowater =  1026d0             ! water density [kg/m3]

!-----------------------------------------------------------------------
!           Disks physical parameters
!-----------------------------------------------------------------------

    e_modul         =  6d9            ! elastic moduli [N/m^2]
    poiss_ratio     =  33d-2          ! poisson ratio nu
    friction_coeff  =  7d-1           ! friction coefficient mu
    rest_coeff      =  88d-2          ! coefficient of restitution

    mass  =  rhoice * pi * h * r ** 2
    ec    =  e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    gc    =  e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    beta  =  log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )


end subroutine get_default