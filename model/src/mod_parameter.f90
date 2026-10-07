!=======================================================================
!   Module: program-wide parameters
!
! `n` is set at runtime (read from the input file in godar.f90) and
! used to allocate every per-particle array via the init_<module>(n)
! procedures defined in each data module. `master` and the legacy
! rectangular-domain dimensions are compile-time constants.
!=======================================================================

module parameters

    use kind_parameter, only: dp

    implicit none

    ! number of particles (set at runtime)
    integer :: n = 0

    ! master mpi rank
    integer, parameter :: master = 0

    ! legacy rectangular-domain dimensions [m]; unused with mask
    real(dp), parameter :: nx = 10.0e3_dp
    real(dp), parameter :: ny = 10.0e3_dp

end module parameters
