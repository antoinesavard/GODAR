subroutine get_default
    
    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"


    !-------------------------------------------------------------------
    !           set run options
    !-------------------------------------------------------------------

    dynamics  = .true.             ! dynamical forcings
    slipping  = .true.             ! enables slipping
    thermodyn = .false.            ! thermo or not
    cohesion  = .true.             ! bonds/no bond
    ridging   = .true.             ! plastic behavior at contact
    shelter   = .true.             ! sheltering from forcings
    
    !-------------------------------------------------------------------
    !           set parameter for the run
    !-------------------------------------------------------------------

    rtree    =  1d3                ! search radius in kd-tree
    ntree    =  1d0                ! time steps between tree builds

    dt       =  1d-3               ! length of time step [s]
    nt       =  1d3                ! number of tstep
    t        =  nt * dt            ! length of time [s]
    comp     =  1d0                ! output steps

    pi       =  4d0 * datan(1d0)   ! pi

    !-------------------------------------------------------------------
    !           physical parameters
    !-------------------------------------------------------------------

    Cdair    =  3d-1               ! body drag air [Lupkes, 2012]
    Csair    =  5d-4               ! surface drag air
    Cdwater  =  3d-1               ! body drag water [Lupkes, 2012]
    Cswater  =  2d-3               ! surface drag water

    z0w      =  3.0d-4             ! viscosity limit of water
    lat      =  1.4                ! latitude ~80 degrees   
    
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
    sigmanc_crit    =  1d6         ! critical normal stress

    ! effective contact modulus
    ec    =  e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    ! effective shear modulus
    gc    =  e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    ! damping ratio
    beta =  1d1 * log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )

	!-------------------------------------------------------------------
    !           Bonds physical parameters
    !-------------------------------------------------------------------

	eb			 = 6d9
	lambda_rb	 = 8d-1
	lambda_lb	 = 1d0
	sigmatb_crit = 1d5
	sigmacb_crit = 1d6
	tau_crit	 = 1d6
	gamma_d		 = 1d0
    bond_lim     = 1d-2

    !-------------------------------------------------------------------
    !           Winds and currents forcings
    !-------------------------------------------------------------------

    uw = 0d0
    vw = 0d0
    ua = 0d0
    va = 0d0

    !-------------------------------------------------------------------
    !           Forcings on the plates
    !-------------------------------------------------------------------    
    pfn = 5d7
    pfs = 5d7

    !-------------------------------------------------------------------
    !           Default input files to use
    !-------------------------------------------------------------------

    Xfile = "files/x.dat"
    Yfile = "files/y.dat"
    Rfile = "files/r.dat"
    Hfile = "files/h.dat"
    Tfile = "files/theta.dat"
    Ofile = "files/omega.dat"

end subroutine get_default


subroutine read_namelist (namelist_name)

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_bond.h"
    include "CB_forcings.h"
    include "CB_options.h"
    
    integer :: nml_error, filenb
    logical :: exist
    character filename*32
    character(16), intent(in) :: namelist_name

    !---- namelist variables -------------------------------------------
    namelist /options_nml/ &
        dynamics, slipping, thermodyn, cohesion, ridging, shelter
    
    namelist /numerical_param_nml/ &
        rtree, ntree, dt, nt, comp

    namelist /physical_param_nml/ &
        Cdair, Csair, Cdwater, Cswater, z0w, lat, rhoair, &
        rhoice, rhowater

    namelist /disk_param_nml/ &
        e_modul, poiss_ratio, friction_coeff, rest_coeff, sigmanc_crit

    namelist /bond_param_nml/ &
        eb, lambda_rb, lambda_lb, sigmatb_crit, &
        sigmacb_crit, tau_crit, gamma_d, bond_lim

    namelist /forcings_nml/ &
        uw, vw, ua, va 

    namelist /plates_nml/ &
        pfn, pfs

    namelist /input_files_nml/ &
        Xfile, Yfile, Rfile, Hfile, Tfile, Ofile
    !-------------------------------------------------------------------

    filename = 'namelist/' // trim(adjustl(namelist_name))
    filenb = 10

    print '(a)', &
        '',&
        '|--------------------------------------------------------|',&
        '|                                                        |',&
        '|   Reading namelist values to update defaults.          |',&
        '|                                                        |',&
        '|--------------------------------------------------------|',&
        ''  
    

    inquire (file=filename, exist=exist)

    if (exist .eqv. .False.) then
        print *, "Error: input file does not exist: ", filename
        print *, "Default will be used"
        return
    end if

    ! open file
    open (filenb, file=filename, status='old', iostat=nml_error)
        
    ! read values inside file 
    print*,'Reading run options'
    read(filenb, nml=options_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading numerical parameters'
    read(filenb, nml=numerical_param_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading physical parameters'
    read(filenb, nml=physical_param_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading disk parameters'
    read(filenb, nml=disk_param_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading bond parameters'
    read(filenb, nml=bond_param_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading forcings'
    read(filenb, nml=forcings_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading plate forces'
    read(filenb, nml=plates_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    print*,'Reading input files names'
    read(filenb, nml=input_files_nml,iostat=nml_error)
    if (nml_error /= 0) then
        print *, '  error, default will be used', nml_error
    else
        print *, '  pass'
    end if

    close(filenb)

    ! warnings
    if (slipping .eqv. .false.) then
        print '(a)', &
    '',&
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',&
    '!                                                        !',&
    "!    Watch out! Slipping is required for                 !",&
    "!    numerical stability of the tangent forces!          !",&
    '!                                                        !',&
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',&
    ''
    end if

    ! recompute compound variables with updated values of parameters
    t = nt * dt
    ec = e_modul / ( 2 * ( 1 - poiss_ratio ** 2 ) )
    gc = e_modul / ( 4 * ( 1 - poiss_ratio ) * ( 2 + poiss_ratio ) )
    beta =  1d1 * log(rest_coeff) / sqrt( log(rest_coeff) ** 2 + pi ** 2 )

end subroutine read_namelist
