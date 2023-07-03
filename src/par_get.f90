subroutine get_default
    
    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"

    !-------------------------------------------------------------------
    !           set parameter for the run
    !-------------------------------------------------------------------

    rtree    =  1d3                ! search radius in kd-tree

    dt       =  1d-3               ! length of time step [s]
    nt       =  1d3                ! number of tstep
    t        =  nt * dt            ! length of time [s]
    comp     =  1d0                ! output steps

    pi       =  4d0 * datan(1d0)   ! pi

    !-------------------------------------------------------------------
    !           physical parameters
    !-------------------------------------------------------------------

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

    e_modul         =  6d9         ! elastic moduli
    poiss_ratio     =  33d-2       ! poisson ratio nu
    friction_coeff  =  7d-1        ! friction coefficient mu
    rest_coeff      =  88d-2       ! coefficient of restitution

    ! effective contact modulus
    ec    =  e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    ! effective shear modulus
    gc    =  e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    ! damping ratio
    beta  =  log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )

	!-------------------------------------------------------------------
    !           Bonds physical parameters
    !-------------------------------------------------------------------

	eb			= 6d9
	lambda_rb	= 8d-1
	lambda_lb	= 1d0
	lambda_ns   = 2d0
	sigmatb_max	= 1d5
	sigmanb_max	= 1d6
	tau_max		= 1d6
	gamma_d		= 1d0

    !-------------------------------------------------------------------
    !           Winds and currents forcings
    !-------------------------------------------------------------------

    uw = 0d0
    vw = 0d0
    ua = 0d0
    va = 0d0

    !-------------------------------------------------------------------
    !           Default input files to use
    !-------------------------------------------------------------------

    Xfile = "input_files/X.dat"
    Yfile = "input_files/Y.dat"
    Rfile = "input_files/R.dat"
    Hfile = "input_files/H.dat"

end subroutine get_default


subroutine read_namelist

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    
    integer :: nml_error, filenb
    character filename*32

    !---- namelist variables -------------------------------------------
    namelist /numerical_param_nml/ &
        rtree, dt, nt, comp

    namelist /physical_param_nml/ &
        Cdair, Csair, Cdwater, Cswater, z0w, rhoair, rhoice, rhowater

    namelist /disk_param_nml/ &
        e_modul, poiss_ratio, friction_coeff, rest_coeff

    namelist /bond_param_nml/ &
        eb, lambda_rb, lambda_lb, lambda_ns, sigmatb_max, &
        sigmanb_max, tau_max, gamma_d

    namelist /input_files_nml/ &
        Xfile, Yfile, Rfile, Hfile
    !-------------------------------------------------------------------

    filename ='namelist.nml'
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
        if (nml_error /= 0) exit
        print*,'Reading input files names parameters'
        read(filenb, nml=physical_param_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading input files names parameters'
        read(filenb, nml=disk_param_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading input files names parameters'
        read(filenb, nml=bond_param_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading input files names parameters'
        read(filenb, nml=input_files_nml,iostat=nml_error)
        print *, nml_error
    enddo

    close(filenb)

    ! recompute compound variables with updated values of parameters
    t = nt * dt
    ec = e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    gc = e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    beta = log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )

end subroutine read_namelist
