subroutine write_file_1D(array, path, description, units, expno_str, var, tstep)

    use netcdf

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    double precision, intent(in) :: array(n)
    character(len=20), intent(in) :: path
    character(len=*), intent(in) :: var, units, description
    character(len=4), intent(in) :: expno_str
    integer, intent(in) :: tstep


    integer :: file_id, xdim_id, time_id
    integer :: array_id
    integer, dimension(2) :: arrdims

    integer :: i
    integer :: ierr, time_index

    logical :: file_exists

    ! File status
    inquire(file=path, exist=file_exists)

    i = size(array,1)

    ! Create or open the NetCDF file
    if (file_exists) then
        ! Open the existing file in write mode
        ierr = nf90_open(path, NF90_WRITE, file_id)
    else
        ! Create a new file if it doesn't exist
        ierr = nf90_create(path=path, cmode=NF90_CLOBBER, ncid=file_id)
        
        ! Define the dimensions
        ierr = nf90_def_dim(file_id, 'i', i, xdim_id)
        ierr = nf90_def_dim(file_id, 'tstep', &
                            int(nt/comp), time_id)

        arrdims = (/time_id, xdim_id/)           

        ierr = nf90_def_var(file_id, var,  NF90_DOUBLE, arrdims, array_id)

        ! Assign attributes to the array
        ierr = nf90_put_att(file_id, array_id, 'units', units)
        ierr = nf90_put_att(file_id, array_id, 'dt', &
                            comp*dt)
        
        ! assign attributes to the file
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'expno', expno_str)
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'Description', &
                            description)


        ! End define mode
        ierr = nf90_enddef(file_id)
    endif

    ! which variable are you looking at
    ierr = nf90_inq_varid(file_id, var, array_id)

    ! Find the next available index for time dimension (appending)
    time_index = int(tstep/comp) 

    ! Write the array data to the appropriate time index
    ierr = nf90_put_var(file_id, array_id, array, &
                        start=(/time_index, 1/), count=(/1, i/))

    ! Close the NetCDF file
    ierr = nf90_close(file_id)

end subroutine write_file_1D


subroutine write_file_int_2D(array, path, description, units, expno_str, var, tstep)

    use netcdf

    implicit none

    include "parameter.h"
    include "CB_variables.h"
    include "CB_const.h"
    include "CB_bond.h"

    integer, intent(in) :: array(n,n)
    character(len=20), intent(in) :: path
    character(len=*), intent(in) :: var, units, description
    character(len=4), intent(in) :: expno_str
    integer, intent(in) :: tstep

    integer :: file_id, xdim_id, ydim_id, time_id
    integer :: array_id
    integer, dimension(3) :: arrdims

    integer :: i, j
    integer :: ierr, time_index

    logical :: file_exists

    ! File status
    inquire(file=path, exist=file_exists)

    i = size(array,1)
    j = size(array,2)

    ! Create or open the NetCDF file
    if (file_exists) then
        ! Open the existing file in write mode
        ierr = nf90_open(path, NF90_WRITE, file_id)
    else
        ! Create a new file if it doesn't exist
        ierr = nf90_create(path=path, cmode=NF90_CLOBBER, ncid=file_id)
        
        ! Define the dimensions
        ierr = nf90_def_dim(file_id, 'i', i, xdim_id)
        ierr = nf90_def_dim(file_id, 'j', j, ydim_id)
        ierr = nf90_def_dim(file_id, 'tstep', &
                            int(nt/comp), time_id)

        ! Define the variable based on the dimensions
        arrdims = (/time_id, xdim_id, ydim_id/)
        ierr = nf90_def_var(file_id, var,  NF90_INT, arrdims, array_id)

        ! Assign attributes to the array
        ierr = nf90_put_att(file_id, array_id, 'units', units)
        ierr = nf90_put_att(file_id, array_id, 'dt', &
                            comp*dt)

        ! assign attributes to the file
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'expno', expno_str)
        ierr = nf90_put_att(file_id, NF90_GLOBAL, 'Description', &
                            description)


        ! End define mode
        ierr = nf90_enddef(file_id)
    endif

    ! which variable are you looking at
    ierr = nf90_inq_varid(file_id, var, array_id)

    ! Find the next available index for time dimension (appending)
    time_index = int(tstep/comp)

    ! Write the array data to the appropriate time index
    ierr = nf90_put_var(file_id, array_id, array, &
                        start=(/time_index, 1, 1/), count=(/1, i, j/))

    ! Close the NetCDF file
    ierr = nf90_close(file_id)

end subroutine write_file_int_2D