integer :: n, master

double precision :: nx, ny

parameter(                       &
            nx    = 70d3,        & ! x - size of the domain          [m]
            ny    = 30d3,        & ! y - size of the domain          [m]
            n     = 525,        & ! number of disks
            master= 0            & ! master process
        )
