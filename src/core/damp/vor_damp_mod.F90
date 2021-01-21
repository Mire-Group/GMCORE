module vor_damp_mod
  use mpi
  use flogger
  use string
  use const_mod
  use process_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use tridiag_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: cv_full_lat(:,:)
  real(r8), allocatable :: cv_half_lat(:,:)

  real(r8) , allocatable :: rhs_all(:) , v_all(:)  ! send to solver , one member
  real(r8) , allocatable :: rhs_tran(:,:) , v_tran(:,:) ! mpi all member

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine vor_damp_init()

    integer j, j0, jr, k
    real(r8) a, b

    call vor_damp_final()

    ! Only do vorticity damping in reduced regions.
    ! First, find the interface when reduce starts.
    j0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat_deg(j) >= -vor_damp_lat0) then
        j0 = j
        exit
      end if
    end do
    if (is_root_proc()) then
      call log_notice('Vorticity damping control latitude index is ' // to_str(j0) // '.')
    end if

    allocate(cv_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cv_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    if (proc%NeedReduce) then
      allocate(rhs_all(global_mesh%num_full_lon))
      allocate(v_all(global_mesh%num_full_lon))
      allocate(rhs_tran(proc%member_num ,global_mesh%num_full_lon)) 
      allocate(v_tran(proc%member_num , global_mesh%num_full_lon))
    end if

    select case (vor_damp_order)
    case (2)
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          cv_full_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(vor_damp_decay) / j0**2) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          cv_half_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(vor_damp_decay) / j0**2) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_str(vor_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    allocate(use_implicit_solver(global_mesh%num_half_lat))
    use_implicit_solver = .false.
    allocate(zonal_solver(global_mesh%num_half_lat,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
        if (global_mesh%half_lat(j) <= 0) then
          jr = j - global_mesh%half_lat_ibeg_no_pole + 1
        else
          jr = global_mesh%half_lat_iend_no_pole - j + 1
        end if
        if (jr > size(reduce_factors)) cycle
        if (reduce_factors(jr) > 1) then
          use_implicit_solver(j) = .true.
          if (k > 1) then
            if (cv_half_lat(j,k) == cv_half_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -cv_half_lat(j,k) * (1 - beta) * dt_in_seconds / global_mesh%le_lat(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(global_mesh%num_full_lon, a, b)
        end if
      end do
    end do

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(cv_full_lat)) deallocate(cv_full_lat)
    if (allocated(cv_half_lat)) deallocate(cv_half_lat)

    if (allocated(rhs_all))  deallocate(rhs_all)
    if (allocated(rhs_tran)) deallocate(rhs_tran)
    if (allocated(v_all))    deallocate(v_all)
    if (allocated(v_tran))   deallocate(v_tran)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(member_num , block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer status(MPI_STATUS_SIZE), ierr
    integer i, j, k, im
    integer proc_length

    mesh => state%mesh
    proc_length = global_mesh%num_half_lon / proc%cart_dims(1)

    select case (vor_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            ! state%u(i,j,k) = state%u(i,j,k) - dt * cv_full_lat(j,k) * ( &
            !   state%vor(i,j+1,k) - state%vor(i,j,k)) / mesh%le_lon(j)
#else
            state%u(:,i,j,k) = state%u(:,i,j,k) - dt * cv_full_lat(j,k) * ( &
              state%vor(:,i,j,k) - state%vor(:,i,j-1,k)) / mesh%le_lon(j)
#endif
          end do
        end do
      end do
      call fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              rhs(:,i) = state%v(:,i,j,k) + &
                         cv_half_lat(j,k) * beta * dt / mesh%le_lat(j) * (       &
                         state%vor(:,i,j,k) - state%vor(:,i-1,j,k)                 &
                         ) +                                                     &
                       cv_half_lat(j,k) * (1 - beta) * dt / mesh%le_lat(j) * ( &
#ifdef V_POLE
                        !  state%u(i-1,j,k) - state%u(i-1,j-1,k) -               &
                        !  state%u(i  ,j,k) + state%u(i  ,j-1,k)                 &
#else
                         state%u(:,i-1,j+1,k) - state%u(:,i-1,j,k) -               &
                         state%u(:,i  ,j+1,k) + state%u(:,i  ,j,k)                 &
#endif
                       ) / mesh%de_lat(j)
            end do

            if (proc%cart_dims(1) == 1) then
              do im = 1 , member_num
                call zonal_solver(j,k)%solve(rhs(im,:), state%v(im , mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) )
              end do
            else
              if (proc%cart_coords(1) == 0) then 

                rhs_tran(:,1:proc_length) = rhs

                do i = 2 , proc%cart_dims(1) 
                  call MPI_RECV(rhs_tran(: , (i - 1) * proc_length + 1 : i  * proc_length ) , proc_length * member_num , MPI_DOUBLE , proc%id + proc%cart_dims(2) * (i-1)  , 100 , proc%comm , status,ierr)
                end do

                do im = 1 , member_num
                  rhs_all = rhs_tran(im , :)
                  call zonal_solver(j,k)%solve(rhs_all, v_all)
                  v_tran(im , :) = v_all
                end do

                do i = 2 , proc%cart_dims(1) 
                  call MPI_SEND(v_tran(: ,(i - 1) * proc_length + 1 : i  * proc_length) , proc_length * member_num , MPI_DOUBLE , proc%id + proc%cart_dims(2) * (i-1) , 100 , proc%comm , status,ierr)
                end do
                state%v(:,mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) = v_tran(: ,1:proc_length)
              else
                call MPI_SEND(rhs , size(rhs) , MPI_DOUBLE , proc%id - proc%cart_dims(2) * proc%cart_coords(1) , 100 , proc%comm , ierr)
                call MPI_RECV(state%v(:,mesh%full_lon_ibeg:mesh%full_lon_iend,j,k) , proc_length * member_num , MPI_DOUBLE , proc%id - proc%cart_dims(2) * proc%cart_coords(1) , 100 , proc%comm , status,ierr)
              end if
            end if
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              state%v(:,i,j,k) = state%v(:,i,j,k) + dt * cv_half_lat(j,k) * ( &
                state%vor(:,i,j,k) - state%vor(:,i-1,j,k)) / mesh%le_lat(j)
            end do
          end if
        end do
      end do
      call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    case (4)
    end select

  end subroutine vor_damp_run

end module vor_damp_mod
