!=======================================================================
!   Module: thermodynamics variables (placeholder)
!
! Thermo variables and forcings are not yet implemented. This module
! exists so that thermo.f90 can `use` it now and pick up real fields
! later without further refactoring at the call sites.
!=======================================================================

module thermo_var

    use kind_parameter, only: dp
    use parameters, only: n

    implicit none

    ! TODO: declare per-particle thermodynamic state
    !       e.g. real(dp) :: tice(n), tocean(n), tair(n) ...

    !-------------------------------------------------------------------
    ! Sketches kept from the former include files (model/inc), none of
    ! them was in use. The gridded ones (0:nx+1, 0:ny+1) come from a
    ! C-grid model and would need a grid under the particles.
    !-------------------------------------------------------------------

    ! from CB_thermo_dym.h: dimensional thermodynamic coefficients
    ! (defined on the C-grid)
    !
    ! double precision :: Klat_ia, Klat_oa
    ! double precision :: Ksens_ai, Ksens_ao, Ksens_io
    ! double precision :: Kice, Kocn, Kadvo, Hocn
    ! double precision :: Kemis_i, Kemis_al, Kemis_o
    !
    ! common/thermo_dim/    &
    !         Klat_ia,      & ! LH transfer coefficient (ice/atm)
    !         Klat_oa,      & ! LH transfer coefficient (ocn/atm)
    !         Ksens_ai,     & ! SH transfer coefficient (ice/atm)
    !         Ksens_ao,     & ! SH transfer coefficient (ocn/atm)
    !         Ksens_io,     & ! SH transfer coefficient (ocn/ice)
    !         Kice,         & ! ice thermal conductivity
    !         Kocn,         & ! ocean diffusion coefficient [m2/s]
    !         Kadvo,        & ! adv heat transfer coefficient (ocn)
    !         Hocn,         & ! mixed layer depth
    !         Kemis_i,      & ! emissivity ice
    !         Kemis_al,     &
    !         Kemis_o

    ! from CB_thermo_forcing.h: thermodynamic forcing (air temp, ocn
    ! temp and shortwave radiation, defined on the C-grid)
    !
    ! double precision ::                   &
    !     To_clim     (0:nx+1,0:ny+1,0:13), &
    !     Qsw         (0:nx+1,0:ny+1)
    !
    ! common/ThermoForcing/   &
    !     To_clim,            & ! monthly clim ocean ML temperature
    !     Qsw                   ! shortave radiation from the sun

    ! from CB_thermo_var.h: thermodynamic variables (defined on the
    ! C-grid)
    !
    ! double precision ::                  &
    !     Ta       (0:nx+1,0:ny+1), &
    !     Ti       (0:nx+1,0:ny+1), &
    !     Tl       (0:nx+1,0:ny+1), &
    !     To       (0:nx+1,0:ny+1), &
    !     Qoa      (0:nx+1,0:ny+1), &
    !     Qoa_f    (0:nx+1,0:ny+1), &
    !     Qia      (0:nx+1,0:ny+1), &
    !     Qsh_io   (0:nx+1,0:ny+1), &
    !     Qadvdiff (0:nx+1,0:ny+1), &
    !     Sh       (1:nx,1:ny)    , &
    !     SA       (1:nx,1:ny)    , &
    !     Pvap     (0:nx+1,0:ny+1)
    !
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
    !
    ! double precision ::     & ! particle variables
    !     tice                  ! ice temperature
    !
    ! double precision ::     & ! thermodynamical forcings
    !     tocean         ,    & ! temperature of the ocean
    !     tair           ,    & ! temperature of the air
    !     tland
    !
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
    !
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
    !
    ! common/thermo_var/      & ! particle variables
    !     tice                  ! x positions                          [m]
    !
    ! common/thermo_forcing/  & ! particle variables
    !     tocean         ,    & ! ocean temp                           [m]
    !     tair           ,    & ! ocean temp                           [m]
    !     tland

end module thermo_var
