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

end module thermo_var
