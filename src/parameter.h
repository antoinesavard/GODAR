integer :: n, master

double precision :: nx, ny

parameter(                       &
            nx    = 1d5,         & ! x - size of the domain          [m]
            ny    = 1d5,         & ! y - size of the domain          [m]
            n     = 1000,        & ! number of disks
            master= 0            & ! master process
        )
