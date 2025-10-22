program godar

    use datetime, ONLY: datetime_init, datetime_str, datetime_type, delta_init
    use datetime, ONLY: datetime_delta_type, OPERATOR(+), OPERATOR(<), OPERATOR(==), OPERATOR(-)
    use datetime, ONLY: str2dt, datetime_str_6, now, delta_str
    use omp_lib
    use mpi_f08

    implicit none

    include "parameter.h"
    include "CB_const.h"
    include "CB_mpi.h"

    integer :: tstep
    integer :: expno, readnamelist, restart, expno_r, nt_r
    integer :: proc_num, thread_num, num_threads
    type(datetime_type) :: tic, tac
    character(len=4) :: expno_str, expno_str_r
    character(len=16) :: namelist_name
    character(10) :: n_str

    ! execute godar only on master rank
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, n_ranks, ierr)

    !---------------------------------------------------------------
    !       Read run information
    !---------------------------------------------------------------

    if ( rank .eq. master ) then

        print '(a)', &
        '',&
        '|========================================================|',&
        '|                                                        |',&
        '|   Welcome in GODAR (Granular flOes for Discrete        |',&
        '|   Arctic Representation). Please use the input file    |',&
        '|   as arguments to this program. If not, please provide |',&
        '|   the following:                                       |',&
        '|                                                        |',&
        '|========================================================|',&
        ''  

        ! read if inputing namelist or not
        print *, 'Read namelist? (0/1)'
        read  *, readnamelist
        print *, readnamelist

        if (readnamelist .eq. 1) then
            print *, 'Namelist name? (XXnamelist.nml)'
            read  *, namelist_name
            print *, namelist_name
        end if

        ! read if restarting or not
        print *, 'Restart experiment from previous one? (0/1)'
        read  *, restart
        print *, restart

        if (restart .eq. 1) then
            ! read restarting experience number
            print *, "Restart from experiment number: (XX)"
            read *, expno_r
            print *, expno_r
            write(expno_str_r,'(i2.2)') expno_r
            
            print *, "Last iteration of restart experiment is:"
            read *, nt_r
            print *, nt_r
        endif

        ! read experience number
        print *, 'This experiment number? (XX)'
        read  *, expno 
        print *, expno 
        write(expno_str,'(i2.2)') expno
        write(n_str,'(i0)') n

        ! number of particles
        print *, 'Number of particles is this simulation:'
        print *, n
        
    end if

    call get_default
    
    if ( rank .eq. master ) then
        ! overwrite default based on namelist
        if (readnamelist .eq. 1) then
            call read_namelist (namelist_name)
        endif
        call ini_get (restart, expno_str_r, nt_r)

        call clear_posts (expno_str)
        
        call info (expno_str, restart)

        !number of processor/threads
        print '(a)', &
        '',&
        '|--------------------------------------------------------|',&
        '|                                                        |',&
        '|   Checking the number of threads to use.               |',&
        '|                                                        |',&
        '|--------------------------------------------------------|',&
        '' 

        proc_num = omp_get_num_procs ( )
        thread_num = omp_get_max_threads ( )

        print *, 'Number of processes: ', n_ranks
        print *, 'Number of processors available: ', proc_num
        print *, 'Number of threads available:    ', thread_num
        print *, 'Number of threads to use?'
        
        ! number of threads to use per rank
        read  *, num_threads
        print *, num_threads

        if (num_threads .gt. thread_num) then
            print*, "Your requested thread number is going to be set to max value of available thread, which is: ", thread_num
            num_threads = thread_num
        end if



        print *, 'Total number of threads: ', n_ranks * num_threads

        if (n_ranks > 1) then

            print '(a)', &
        '',&
        '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',&
        '!                                                        !',&
        "!    Watch out! MPI is not working properly :'(          !",&
        '!                                                        !',&
        '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',&
        '' 
        end if

        print '(a)', &
        '',&
        '|--------------------------------------------------------|',&
        '|                                                        |',&
        '|   Entering main program now. Please wait...            |',&
        '|                                                        |',&
        '|--------------------------------------------------------|',&
        '' 
        
        print *, 'Broadcasting initialization to all processes'
        print *, 'Verify that you loop over all particles'

    end if

    ! broadcast the data to all the other ranks
    call broadcasting_ini (num_threads)
    
    ! set openmp in all ranks
    call mpi_barrier (mpi_comm_world, ierr)
    call omp_set_num_threads (num_threads)
    
    !-------------------------------------------------------------------
    ! Main (time) loop of the program
    !-------------------------------------------------------------------

    ! set mpi variables
    ! Compute the part of the array to loop over given rank
    iter_per_rank = n / n_ranks
    
    if ( mod(n, n_ranks) > 0 ) then
        iter_per_rank = iter_per_rank + 1
    end if

    first_iter = rank * iter_per_rank + 1
    
    last_iter  = first_iter + iter_per_rank - 1

    !first_iter = int(n * ( 1 - &
    !            sqrt(1. * (n_ranks - rank) / n_ranks) )) + 1

    !last_iter = int(n * ( 1 - &
    !            sqrt(1. * (n_ranks - rank - 1) / n_ranks) ))

    print*, "rank: ", rank, "particle id (start-stop): ", first_iter, last_iter
    
    if ( rank .eq. master ) then
        tic = now()
        print *, "Time step: ", 1, "/", int(nt)
    end if

    do tstep = 1, int(nt) + 1

        call stepper (tstep)

        if ( rank .eq. master ) then
            if (MODULO(tstep, int(comp)) .eq. 0) then

                call sea_ice_post (tstep, expno_str)
                if (MODULO(tstep, int(comp*10)) .eq. 0) then
                    print *, "Time step: ", tstep, "/", int(nt)
                end if

            endif
        end if

    end do

    if ( rank .eq. master ) then
        tac = now()

        print*, "Total simulation time: ", delta_str(tac - tic)

        print '(a)', &
        '',&
        '|--------------------------------------------------------|',&
        '|                                                        |',&
        '|   Exiting main program now.                            |',&
        '|                                                        |',&
        '|--------------------------------------------------------|',&
        '' 
    end if

    call mpi_finalize(ierr)
    
end program godar
