integer :: n

double precision :: nx, ny, rtree

parameter(                       &
            nx    = 1d5,         & ! x - size of the domain 
            ny    = 1d5,         & ! y - size of the domain 
            n     = 2500,        & ! number of disks
            rtree = 1d3          & ! search radius in kd-tree
        )
