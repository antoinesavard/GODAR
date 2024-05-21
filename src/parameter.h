integer :: n, master

double precision :: nx, ny

parameter(                       &
            nx    = 35d3,        & ! x - size of the domain          [m]
            ny    = 35d3,        & ! y - size of the domain          [m]
            n     = 1012,        & ! number of disks
            master= 0            & ! master process
        )
