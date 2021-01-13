module pa_mod
    use mpi

    implicit none



    integer myid;
    character*50 file_name;             !output file name
    logical pa_operator_prepare_flag    !whether do time counting
    logical pa_space_operator_flag
    logical pa_reduce_flag
    logical pa_flag

    real*8    sub_time_start , sub_time_end , cal_time_start , cal_time_end , tran_time_start , tran_time_end ,step_time_start , step_time_end
    real*8    sub_time , cal_time , tran_time



contains

  subroutine pa_init(comm)
    
    integer, intent(in), optional :: comm

    integer ierr
    character*10 date
    character*10 time
    character*10 zone
    integer date_time(8)
    character*10 s_id 
    character*50 command_cd , command_mk
    logical istatus1,istatus2
    
    call mpi_comm_rank(comm,myid,ierr);
    write(s_id,"(i4.4)") myid 

    write(s_id,"(i4.4)") myid
    file_name(1:18)  = 'mpi_operator_time_'
    file_name(19:22) = s_id
    file_name(23:26) = '.txt'
    open(unit=(10000+myid),POSITION='APPEND',file=file_name)

    file_name(1:21)  = 'mpi_stateupdate_time_'
    file_name(22:25) = s_id
    file_name(26:29) = '.txt'
    open(unit=(20000+myid),POSITION='APPEND',file=file_name)
  
  end subroutine pa_init

  subroutine Get_Start_Time(get_time)
    real*8 , intent(inout) :: get_time

    get_time = mpi_wtime();
  end subroutine Get_Start_Time

  subroutine Get_End_Time(get_time)
    real*8 , intent(inout) :: get_time

    get_time = mpi_wtime();
  end subroutine Get_End_Time

  subroutine Get_Time_Init()
    cal_time = 0;
    tran_time = 0;
  end subroutine Get_Time_Init

  subroutine pa_final()
    close((10000+myid))
  end subroutine pa_final

end module pa_mod
