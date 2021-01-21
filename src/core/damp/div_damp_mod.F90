module div_damp_mod
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

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: cd_full_lat(:,:)
  real(r8), allocatable :: cd_half_lat(:,:)
  real(r8), allocatable :: rhs_all(:)
  real(r8), allocatable :: u_all(:)
  real(r8), allocatable :: rhs_tran(:,:)
  real(r8), allocatable :: u_tran(:,:)

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine div_damp_init()

    integer j, k, r, jr, j0
    real(r8) a, b

    call div_damp_final()

    j0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        if (jr > size(reduce_factors)) exit
        if (reduce_factors(jr) > 1) j0 = jr
      end if
    end do
    j0 = max(div_damp_j0, j0)

    allocate(cd_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cd_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    if (proc%NeedReduce) then
      allocate(rhs_all(global_mesh%num_half_lon))
      allocate(u_all(global_mesh%num_half_lon))
      allocate(rhs_tran(proc%member_num ,global_mesh%num_half_lon))
      allocate(u_tran(proc%member_num , global_mesh%num_half_lon))
    end if

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          if (baroclinic) then
            cd_full_lat(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_upper * exp(k**2 * log(0.2_r8) / 3**2)) * &
              (global_mesh%full_cos_lat(j)**r + div_damp_polar * exp(jr**2 * log(div_damp_exp) / j0**2)) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            cd_full_lat(j,k) = div_damp_coef2 * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          if (baroclinic) then
            cd_half_lat(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_upper * exp(k**2 * log(0.2_r8) / 3**2)) * &
              (global_mesh%half_cos_lat(j)**r + div_damp_polar * exp(jr**2 * log(div_damp_exp) / j0**2)) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            cd_half_lat(j,k) = div_damp_coef2 * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_str(div_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    allocate(use_implicit_solver(global_mesh%num_full_lat))
    use_implicit_solver = .false.
    allocate(zonal_solver(global_mesh%num_full_lat,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
        if (global_mesh%full_lat(j) <= 0) then
          jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        else
          jr = global_mesh%full_lat_iend_no_pole - j + 1
        end if
        if (jr > size(reduce_factors)) cycle
        if (reduce_factors(jr) > 1) then
          use_implicit_solver(j) = .true.
          if (k > 1) then
            if (cd_full_lat(j,k) == cd_full_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -cd_full_lat(j,k) * (1 - beta) * dt_in_seconds / global_mesh%de_lon(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(global_mesh%num_half_lon, a, b)
        end if
      end do
    end do

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(cd_full_lat)) deallocate(cd_full_lat)
    if (allocated(cd_half_lat)) deallocate(cd_half_lat)

    if (allocated(rhs_all)) deallocate(rhs_all)
    if (allocated(u_all)) deallocate(u_all)
    if (allocated(rhs_tran)) deallocate(rhs_tran)
    if (allocated(u_tran)) deallocate(u_tran)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(member_num , block%mesh%half_lon_ibeg:block%mesh%half_lon_iend)
    integer status(MPI_STATUS_SIZE), ierr
    integer i, j, k, im
    integer proc_length

    mesh => state%mesh
    proc_length = global_mesh%num_half_lon / proc%cart_dims(1)

    select case (div_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            ! state%v(i,j,k) = state%v(i,j,k) + dt * cd_half_lat(j,k) * ( &
            !   state%div(i,j,k) - state%div(i,j-1,k)) / mesh%de_lat(j)
#else
            state%v(:,i,j,k) = state%v(:,i,j,k) + dt * cd_half_lat(j,k) * ( &
              state%div(:,i,j+1,k) - state%div(:,i,j,k)) / mesh%de_lat(j)
#endif
          end do
        end do
      end do
      call fill_halo_member(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              rhs(:,i) = state%u(:,i,j,k) + &
                         cd_full_lat(j,k) * beta * dt / mesh%de_lon(j) * (       &
                         state%div(:,i+1,j,k) - state%div(:,i,j,k)                 &
                       ) +                                                     &
                       cd_full_lat(j,k) * (1 - beta) * dt / mesh%de_lon(j) * ( &
#ifdef V_POLE
                        !  state%v(i+1,j+1,k) - state%v(i+1,j,k) -               &
                        !  state%v(i  ,j+1,k) + state%v(i  ,j,k)                 &
#else
                         state%v(:,i+1,j,k) - state%v(:,i+1,j-1,k) -               &
                         state%v(:,i  ,j,k) + state%v(:,i  ,j-1,k)                 &
#endif
                       ) / mesh%le_lon(j)
            end do
            if (proc%cart_dims(1) == 1) then
              do im = 1 , member_num
                call zonal_solver(j,k)%solve(rhs(im,:), state%u(im,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k))
              end do
            else
              if (proc%cart_coords(1) == 0) then 

                rhs_tran(:,1:proc_length) = rhs

                do i = 2 , proc%cart_dims(1) 
                  call MPI_RECV(rhs_tran(:, (i - 1) * proc_length + 1 : i  * proc_length ) , proc_length * member_num , MPI_DOUBLE , proc%id + proc%cart_dims(2) * (i-1)  , 100 , proc%comm , status,ierr)
                end do

                do im = 1 , member_num
                  rhs_all = rhs_tran(im , :)
                  call zonal_solver(j,k)%solve(rhs_all, u_all)
                  u_tran(im,:) = u_all
                end do

                do i = 2 , proc%cart_dims(1) 
                  call MPI_SEND(u_tran(:, (i - 1) * proc_length + 1 : i  * proc_length) , proc_length * member_num , MPI_DOUBLE , proc%id + proc%cart_dims(2) * (i-1) , 100 , proc%comm , status,ierr)
                end do

                state%u(:,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k) = u_tran(:,1:proc_length)
              else
                call MPI_SEND(rhs , size(rhs) , MPI_DOUBLE , proc%id - proc%cart_dims(2) * proc%cart_coords(1) , 100 , proc%comm , ierr)
                call MPI_RECV(state%u(:,mesh%half_lon_ibeg:mesh%half_lon_iend,j,k) , proc_length * member_num, MPI_DOUBLE , proc%id - proc%cart_dims(2) * proc%cart_coords(1) , 100 , proc%comm , status,ierr)
              end if
            end if 
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              state%u(:,i,j,k) = state%u(:,i,j,k) + dt * cd_full_lat(j,k) * ( &
                state%div(:,i+1,j,k) - state%div(:,i,j,k)) / mesh%de_lon(j)
            end do
          end if
        end do
      end do

      call fill_halo_member(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)
    case (4)
    end select

  end subroutine div_damp_run

end module div_damp_mod

