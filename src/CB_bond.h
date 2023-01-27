!=======================================================================
!             COMMAND BLOCK: bond properties
!=======================================================================

    double precision :: E

    integer :: bond (n,n)

    common/bond_var/           &
                bond       ,   & ! bond between disks i and j (0 or 1)
                E