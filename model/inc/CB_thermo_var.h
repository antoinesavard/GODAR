!=======================================================================
!     Common block thermo_variables: Thermodynamic variables
!     (defined on the C-grid) 
!=======================================================================


!         double precision ::                  &
!             Ta       (0:nx+1,0:ny+1), &
!             Ti       (0:nx+1,0:ny+1), &
!             Tl       (0:nx+1,0:ny+1), &
!             To       (0:nx+1,0:ny+1), &
!             Qoa      (0:nx+1,0:ny+1), &
!             Qoa_f    (0:nx+1,0:ny+1), &
!             Qia      (0:nx+1,0:ny+1), &
!             Qsh_io   (0:nx+1,0:ny+1), &
!             Qadvdiff (0:nx+1,0:ny+1), &
!             Sh       (1:nx,1:ny)    , &
!             SA       (1:nx,1:ny)    , &
!             Pvap     (0:nx+1,0:ny+1)



! common/ThermoVariables/   &
!       Ta,             & ! air temperature
!       Ti,             & ! ice surface temperature
!       Tl,             & ! land surface temperature
!       To,             & ! ocean mixed layer temperature
!       Qoa,            & ! ocean-atmosphere heat flux
!       Qoa_f,          & ! ocean-atmosphere heat flux (ice growth)
!       Qia,            & ! conductive heat flux through ice
!       Qsh_io,         & ! sensible heat flux (ocn/ice)
!       Qadvdiff,       & ! Advection-Diffusion heat transfer
!       Sh,             & ! thermo source term (h, continuity equation)
!       SA,             & ! thermo source term (A, continuity equation) 
!       Pvap              ! atmospheric vapour pressure




    ! double precision ::	    & ! particle variables
    !     tice                  ! ice temperature
        

    ! double precision ::	    & ! thermodynamical forcings
    !     tocean         ,    & ! temperature of the ocean
    !     tair           ,    & ! temperature of the air
    !     tland

    ! double precision ::     & ! sensible heat and stuff
    !     source_h        ,   &
    !     source_r        ,   &
    !     diff_h          ,   &
    !     diff_r          ,   &
    !     Qsw             ,   &
    !     Qsens           ,   &
    !     Qlat            ,   &
    !     Qlw_down        ,   &
    !     humid_ice       ,   &
    !     humid_atm
    
    ! double precision ::     & ! constants
    !     Q0              ,   &
    !     Clat            ,   &
    !     Ls              ,   &
    !     Csens           ,   &
    !     Cpa             ,   &
    !     emi_atm         ,   &
    !     emi_ocn         ,   &
    !     stefanboltzmann ,   &
    !     ice_albedo      ,   &
    !     ocn_albedo      ,   &
    !     absorp_atm

    ! common/thermo_var/      & ! particle variables
    !     tice                  ! x positions                          [m]

    ! common/thermo_forcing/  & ! particle variables
    !     tocean         ,    & ! ocean temp                           [m]
    !     tair           ,    & ! ocean temp                           [m]
    !     tland