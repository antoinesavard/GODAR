subroutine get_default
    
    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"

    !-------------------------------------------------------------------
    ! set parameter for the run
    !-------------------------------------------------------------------

    dt       =  1d-3               ! length of time step [s]
    nt       =  1d3                ! number of tstep
    t        =  nt * dt            ! length of time [s]
    comp     =  1d0                ! output steps

    pi       =  4d0 * datan(1d0)   ! pi

    Cdair    =  3d-1 * pi / 4      ! body drag air
    Csair    =  5d-4               ! surface drag air
    Cdwater  =  3d-1 * pi / 4      ! body drag water
    Cswater  =  2d-3               ! surface drag water

	z0w      =  3.0d-4             ! viscosity limit of water  
    
    rhoair   =  1.3d0              ! air density [kg/m3]
    rhoice   =  9d02               ! ice density [kg/m3]
    rhowater =  1026d0             ! water density [kg/m3]

    !-------------------------------------------------------------------
    !           Disks physical parameters
    !-------------------------------------------------------------------

    e_modul         =  6d9            ! elastic moduli [N/m^2]
    poiss_ratio     =  33d-2          ! poisson ratio nu
    friction_coeff  =  7d-1           ! friction coefficient mu
    rest_coeff      =  88d-2          ! coefficient of restitution

    ! mass of disk
    mass  =  rhoice * pi * h * r ** 2
    ! effective contact modulus
    ec    =  e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    ! effective shear modulus
    gc    =  e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    ! damping ratio
    beta  =  log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )
	! freeboard height
	hfa   =  h * (rhowater - rhoice) / rhowater
	! drag from water height
	hfw   =  h * rhoice / rhowater


end subroutine get_default


subroutine read_namelist

    implicit none

    include 'parameter.h'
    include 'CB_const.h'
    
    integer :: nml_error, filenb
    character filename*32

    !---- namelist variables -------
    namelist /numerical_param_nml/ &
        nt, dt, comp

    filename ='namelistSIM'
    filenb = 10

    print *, 'Reading namelist values'
    
    open (filenb, file=filename, status='old',iostat=nml_error)
    if (nml_error /= 0) then
        nml_error = -1
    else
        nml_error =  1
    endif
        
    do while (nml_error > 0)
        print*,'Reading numerical parameters'
        read(filenb, nml=numerical_param_nml,iostat=nml_error)
        !if (nml_error /= 0) exit
        print *, nml_error
    enddo

    close(filenb)

    t = nt * dt

end subroutine read_namelist