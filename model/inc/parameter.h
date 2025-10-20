integer :: n, master

double precision :: nx, ny

parameter(                       &
            nx    = 5d3,        & ! x - size of the domain          [m]
            ny    = 5d3,        & ! y - size of the domain          [m]
            n     = 8,         & ! number of disks
            master= 0            & ! master process
        )
