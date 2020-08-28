module process_mod

  use mpi
  use flogger
  use string
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_create_blocks
  public process_stop
  public process_final
  public proc
  public is_root_proc

  public decomp_1d
  public decomp_2d_simple
  public decomp_reduce_south_region
  public decomp_reduce_south_boundary
  public decomp_reduce_north_region
  public decomp_reduce_north_boundary
  public decomp_normal_region
  public decomp_normal_south_boundary
  public decomp_normal_north_boundary

  type process_neighbor_type
    integer :: id       = MPI_PROC_NULL ! ID in global comm
    integer :: local_id = MPI_PROC_NULL ! ID in local comm
    integer :: orient   = 0
    integer :: lon_ibeg = inf_i4
    integer :: lon_iend = inf_i4
    integer :: lat_ibeg = inf_i4
    integer :: lat_iend = inf_i4
  contains
    procedure :: init => process_neighbor_init
  end type process_neighbor_type

  integer, parameter :: decomp_1d = 1
  integer, parameter :: decomp_2d_simple = 2

  integer, parameter :: decomp_reduce_south_region   = 1
  integer, parameter :: decomp_reduce_south_boundary = 2
  integer, parameter :: decomp_reduce_north_region   = 3
  integer, parameter :: decomp_reduce_north_boundary = 4
  integer, parameter :: decomp_normal_region         = 5
  integer, parameter :: decomp_normal_south_boundary = 6
  integer, parameter :: decomp_normal_north_boundary = 7

  type process_type
    integer :: comm           = MPI_COMM_NULL
    integer :: local_comm     = MPI_COMM_NULL
    integer :: zonal_comm     = MPI_COMM_NULL
    integer :: group          = MPI_GROUP_NULL
    integer :: local_group    = MPI_GROUP_NULL
    integer :: zonal_group    = MPI_GROUP_NULL
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    integer :: id             = MPI_PROC_NULL          ! ID in global comm
    integer :: local_id       = MPI_PROC_NULL          ! ID in local comm
    integer idom                                       ! Nest domain index (root domain is 1)
    integer num_lon                                    ! Grid size along longitude
    integer num_lat                                    ! Grid size along latitude
    integer lon_ibeg                                   ! Index of beginning longitude grid
    integer lon_iend                                   ! Index of ending longitude grid
    integer lat_ibeg                                   ! Index of beginning latitude grid
    integer lat_iend                                   ! Index of ending latitude grid

    ! Working variables
    integer decomp_type
    integer decomp_loc
    integer np
    integer np_lon
    integer np_lat
    integer :: nd = 0
    integer :: nd2 = 0

    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes
    type(block_type), allocatable :: blocks(:)
  end type process_type

  type(process_type) proc

contains

  subroutine process_init()

    integer ierr

    call MPI_INIT(ierr)
    proc%comm = MPI_COMM_WORLD
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)

    call setup_mpi_2d()
    call decompose_domains()
    call setup_zonal_comm_for_reduce()
    call connect_parent() ! <-- FIXME: Needs implementation.

  end subroutine process_init

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer ierr

    if (allocated(proc%ngb   )) deallocate(proc%ngb   )
    if (allocated(proc%blocks)) deallocate(proc%blocks)
    if (proc%group       /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group      , ierr)
    if (proc%local_group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%local_group, ierr)
    if (proc%zonal_group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%zonal_group, ierr)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

  pure logical function is_root_proc()

    is_root_proc = proc%id == 0

  end function is_root_proc

  subroutine setup_mpi_1d()

    integer ierr, np2, tmp_comm, tmp_id(1), i
    logical periods(2)

    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)

    proc%decomp_type = decomp_1d
    proc%decomp_loc  = decomp_normal_region

    if (num_proc_lon(1) /= 0 .and. num_proc_lat(1) /= 0) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np2 = 0
      do i = 1, nest_max_dom
        np2 = np2 + num_proc_lon(i) * num_proc_lat(i)
      end do
      if (proc%np /= np2 .and. proc%id == 0) then
        call log_error('Namelist num_proc_lon and num_proc_lat are not compatible with MPI runtime!')
      end if
      ! Set the process topology into proc object.
      np2 = 0
      do i = 1, nest_max_dom
        np2 = np2 + num_proc_lon(i) * num_proc_lat(i)
        if (proc%id + 1 <= np2) then
          proc%cart_dims(1) = num_proc_lon(i)
          proc%cart_dims(2) = num_proc_lat(i)
          proc%idom = i
          exit
        end if
      end do
    else
      proc%cart_dims = 0
      call MPI_DIMS_CREATE(proc%np, 2, proc%cart_dims, ierr)
    end if
    periods = [.true.,.false.]
    if (proc%idom > 1 .and. (nest_lon_beg(proc%idom) /= 0 .or. nest_lon_end(proc%idom) /= 360)) then
      periods(1) = .false.
    end if
    ! Set MPI process topology.
    call MPI_COMM_SPLIT(proc%comm, proc%idom, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%local_comm, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_GROUP(proc%local_comm, proc%local_group, ierr)
    call MPI_COMM_RANK(proc%local_comm, proc%local_id, ierr)
    call MPI_CART_COORDS(proc%local_comm, proc%local_id, 2, proc%cart_coords, ierr)

  end subroutine setup_mpi_1d

  subroutine setup_mpi_2d()

    integer ierr, tag, tmp_comm, np_2d, i
    logical periods(2)

    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)

    proc%nd = count(reduce_factors > 1)
    ! If nd is odd, just use 1D decomposition.
    np_2d = proc%np - proc%nd
    if (np_2d <= 0 .or. mod(proc%nd, 2) /= 0) then
      proc%nd = 0
      call setup_mpi_1d()
      return
    end if
    proc%nd2 = proc%nd / 2
    ! Each process in reduce region gets two zonal circles.

    proc%decomp_type = decomp_2d_simple

    if (.false.) then ! For nesting
    else
      ! Separate processes for reduce regions and for normal regions.
      ! In reduce region, we use 1D merdional decomposition.
      ! In normal region, we use 2D Cartesian decomposition.
      if (proc%id < proc%nd2) then
        tag = 1
        proc%decomp_loc = decomp_reduce_south_region
        if (proc%id == proc%nd2 - 1) proc%decomp_loc = decomp_reduce_south_boundary
        proc%cart_dims = [1, proc%nd2]
      else if (proc%id >= proc%np - proc%nd2) then
        tag = 2
        proc%decomp_loc = decomp_reduce_north_region
        if (proc%id == proc%np - proc%nd2) proc%decomp_loc = decomp_reduce_north_boundary
        proc%cart_dims = [1, proc%nd2]
      else
        tag = 3
        proc%decomp_loc = decomp_normal_region
        call MPI_DIMS_CREATE(np_2d, 2, proc%cart_dims, ierr)
        proc%np_lon = proc%cart_dims(1)
        proc%np_lat = proc%cart_dims(2)
      end if
      ! Let processes in reduce region know the decomposition in normal region.
      call MPI_BCAST(proc%np_lon, 1, MPI_INT, proc%nd2, proc%comm, ierr)
      call MPI_BCAST(proc%np_lat, 1, MPI_INT, proc%nd2, proc%comm, ierr)
      periods = [.true.,.false.]
    end if

    call MPI_COMM_SPLIT(proc%comm, tag, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%local_comm, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_GROUP(proc%local_comm, proc%local_group, ierr)
    call MPI_COMM_RANK(proc%local_comm, proc%local_id, ierr)
    call MPI_CART_COORDS(proc%local_comm, proc%local_id, 2, proc%cart_coords, ierr)

    if (proc%decomp_loc == decomp_normal_region) then
      if (proc%cart_coords(2) == 0) then
        proc%decomp_loc = decomp_normal_south_boundary
      else if (proc%cart_coords(2) == proc%cart_dims(2) - 1) then
        proc%decomp_loc = decomp_normal_north_boundary
      end if
    end if

  end subroutine setup_mpi_2d

  subroutine decompose_domains()

    integer ierr, status(MPI_STATUS_SIZE), i, j, tmp_id(1)
    integer ngb_lon_ibeg, ngb_lon_iend
    integer num_lon, num_lat, res_num

    ! Set neighborhood.
    if (allocated(proc%ngb)) deallocate(proc%ngb)
    select case (proc%decomp_loc)
    case (decomp_normal_region, decomp_normal_south_boundary, decomp_normal_north_boundary)
      allocate(proc%ngb(4))
      call MPI_CART_SHIFT(proc%local_comm, 0, 1, proc%ngb(west )%local_id, proc%ngb(east )%local_id, ierr)
      call MPI_CART_SHIFT(proc%local_comm, 1, 1, proc%ngb(south)%local_id, proc%ngb(north)%local_id, ierr)
      if (proc%ngb(south)%local_id == MPI_PROC_NULL .and. proc%nd2 /= 0) then
        proc%ngb(south)%id = proc%nd2 - 1
      end if
      if (proc%ngb(north)%local_id == MPI_PROC_NULL .and. proc%nd2 /= 0) then
        proc%ngb(north)%id = proc%np - proc%nd2
      end if
    case (decomp_reduce_south_region, decomp_reduce_north_region)
      allocate(proc%ngb(4))
      call MPI_CART_SHIFT(proc%local_comm, 0, 1, proc%ngb(west )%local_id, proc%ngb(east )%local_id, ierr)
      call MPI_CART_SHIFT(proc%local_comm, 1, 1, proc%ngb(south)%local_id, proc%ngb(north)%local_id, ierr)
    case (decomp_reduce_south_boundary)
      allocate(proc%ngb(3 + proc%np_lon))
      call MPI_CART_SHIFT(proc%local_comm, 0, 1, proc%ngb(west )%local_id, proc%ngb(east )%local_id, ierr)
      call MPI_CART_SHIFT(proc%local_comm, 1, 1, proc%ngb(south)%local_id, proc%ngb(north)%local_id, ierr)
      ! Set neighbors from normal region mannually.
      do i = 1, proc%np_lon
        proc%ngb(3+i)%id = proc%nd2 + (i - 1) * proc%np_lat
      end do
    case (decomp_reduce_north_boundary)
      allocate(proc%ngb(3 + proc%np_lon))
      call MPI_CART_SHIFT(proc%local_comm, 0, 1, proc%ngb(west )%local_id, proc%ngb(east )%local_id, ierr)
      call MPI_CART_SHIFT(proc%local_comm, 1, 1, proc%ngb(south)%local_id, proc%ngb(north)%local_id, ierr)
      ! Set neighbors from normal region mannually.
      proc%ngb(south)%id = proc%np - proc%nd2 - proc%np_lon * proc%np_lat + 1
      do i = 2, proc%np_lon
        proc%ngb(3+i)%id = proc%ngb(south)%id + (i - 1) * proc%np_lat
      end do
    end select

    ! Translate local ID of neighbors to global ID.
    do i = 1, size(proc%ngb)
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%local_group, 1, [proc%ngb(i)%local_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    ! Set initial values for num_lon, num_lat, lon_ibeg, lat_ibeg.
    proc%num_lon = global_mesh%num_full_lon
    proc%lon_ibeg = 1
    proc%lat_ibeg = 1
    select case (proc%decomp_loc)
    case (decomp_reduce_south_region, decomp_reduce_south_boundary)
      proc%num_lat  = proc%nd
    case (decomp_reduce_north_region, decomp_reduce_north_boundary)
      proc%num_lat  = proc%nd
#ifdef V_POLE
      proc%lat_ibeg = global_mesh%num_half_lat - proc%nd + 1
#else
      proc%lat_ibeg = global_mesh%num_full_lat - proc%nd + 1
#endif
    case (decomp_normal_region, decomp_normal_south_boundary, decomp_normal_north_boundary)
#ifdef V_POLE
      proc%num_lat = global_mesh%num_half_lat - 2 * proc%nd
#else
      proc%num_lat = global_mesh%num_full_lat - 2 * proc%nd
#endif
    end select

    if (proc%idom > 1) then ! For nesting
      ! Get the start and end indices according to nest domain range.
      ! Zonal direction
      proc%lon_ibeg = 0; proc%lon_iend = 0
      do i = 1, global_mesh%num_full_lon
        if (global_mesh%full_lon_deg(i) >= nest_lon_beg(proc%idom-1)) then
          proc%lon_ibeg = i
          exit
        end if
      end do
      do i = 1, global_mesh%num_full_lon
        if (global_mesh%full_lon_deg(i) >= nest_lon_end(proc%idom-1)) then
          proc%lon_iend = i
          exit
        end if
      end do
      if (proc%lon_iend == 0) proc%lon_iend = global_mesh%num_full_lon
      proc%num_lon = proc%lon_iend - proc%lon_ibeg + 1
      ! Meridional direction
      proc%lat_ibeg = 0; proc%lat_iend = 0
#ifdef V_POLE
      do j = 1, global_mesh%num_half_lat
        if (global_mesh%half_lat_deg(j) >= nest_lat_beg(proc%idom-1)) then
          proc%lat_ibeg = j
          exit
        end if
      end do
      do j = 1, global_mesh%num_half_lat
        if (global_mesh%half_lat_deg(j) >= nest_lat_end(proc%idom-1)) then
          proc%lat_iend = j
          exit
        end if
      end do
      if (proc%lat_iend == 0) proc%lat_iend = global_mesh%num_full_lat
#else
      do j = 1, global_mesh%num_full_lat
        if (global_mesh%full_lat_deg(j) >= nest_lat_beg(proc%idom-1)) then
          proc%lat_ibeg = j
          exit
        end if
      end do
      do j = 1, global_mesh%num_full_lat
        if (global_mesh%full_lat_deg(j) >= nest_lat_end(proc%idom-1)) then
          proc%lat_iend = j
          exit
        end if
      end do
      if (proc%lat_iend == 0) proc%lat_iend = global_mesh%num_full_lat
#endif
      proc%num_lat = proc%lat_iend - proc%lat_ibeg + 1
    end if

    res_num = mod(proc%num_lon, proc%cart_dims(1))
    do i = 0, proc%cart_coords(1) - 1
      if (res_num /= 0 .and. i < res_num) then
        num_lon = proc%num_lon / proc%cart_dims(1) + 1
      else
        num_lon = proc%num_lon / proc%cart_dims(1)
      end if
      proc%lon_ibeg = proc%lon_ibeg + num_lon
    end do
    if (res_num /= 0 .and. proc%cart_coords(1) < res_num) then
      num_lon = proc%num_lon / proc%cart_dims(1) + 1
    else
      num_lon = proc%num_lon / proc%cart_dims(1)
    end if
    proc%num_lon = num_lon
    proc%lon_iend = proc%lon_ibeg + num_lon - 1

    res_num = mod(proc%num_lat, proc%cart_dims(2))
    do j = 0, proc%cart_coords(2) - 1
      if (res_num /= 0 .and. j < res_num) then
        num_lat = proc%num_lat / proc%cart_dims(2) + 1
      else
        num_lat = proc%num_lat / proc%cart_dims(2)
      end if
      proc%lat_ibeg = proc%lat_ibeg + num_lat
    end do
    if (res_num /= 0 .and. proc%cart_coords(2) < res_num) then
      num_lat = proc%num_lat / proc%cart_dims(2) + 1
    else
      num_lat = proc%num_lat / proc%cart_dims(2)
    end if
    proc%num_lat = num_lat
    proc%lat_iend = proc%lat_ibeg + num_lat - 1

    select case (proc%decomp_loc)
    case (decomp_normal_region)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
    case (decomp_normal_south_boundary)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call MPI_SEND(proc%lon_ibeg, 1, MPI_INT, proc%ngb(south)%id, 0, proc%comm, ierr)
      call MPI_SEND(proc%lon_iend, 1, MPI_INT, proc%ngb(south)%id, 1, proc%comm, ierr)
    case (decomp_normal_north_boundary)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call MPI_SEND(proc%lon_ibeg, 1, MPI_INT, proc%ngb(north)%id, 0, proc%comm, ierr)
      call MPI_SEND(proc%lon_iend, 1, MPI_INT, proc%ngb(north)%id, 1, proc%comm, ierr)
    case (decomp_reduce_south_region, decomp_reduce_north_region)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
    case (decomp_reduce_south_boundary)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      do i = 1, proc%np_lon
        call MPI_RECV(ngb_lon_ibeg, 1, MPI_INT, proc%ngb(south+i)%id, 0, proc%comm, status, ierr)
        call MPI_RECV(ngb_lon_iend, 1, MPI_INT, proc%ngb(south+i)%id, 1, proc%comm, status, ierr)
        call proc%ngb(south+i)%init(north, lon_ibeg=ngb_lon_ibeg, lon_iend=ngb_lon_iend)
      end do
    case (decomp_reduce_north_boundary)
      call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
      call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg, lon_iend=proc%lon_iend)
      call MPI_RECV(ngb_lon_ibeg, 1, MPI_INT, proc%ngb(south)%id, 0, proc%comm, status, ierr)
      call MPI_RECV(ngb_lon_iend, 1, MPI_INT, proc%ngb(south)%id, 1, proc%comm, status, ierr)
      call proc%ngb(south)%init(south, lon_ibeg=ngb_lon_ibeg, lon_iend=ngb_lon_iend)
      do i = 2, proc%np_lon
        call MPI_RECV(ngb_lon_ibeg, 1, MPI_INT, proc%ngb(south+i)%id, 0, proc%comm, status, ierr)
        call MPI_RECV(ngb_lon_iend, 1, MPI_INT, proc%ngb(south+i)%id, 1, proc%comm, status, ierr)
        call proc%ngb(south+i)%init(south, lon_ibeg=ngb_lon_ibeg, lon_iend=ngb_lon_iend)
      end do
    end select

  end subroutine decompose_domains

  subroutine setup_zonal_comm_for_reduce()

    integer ierr, i, j, jr
    integer, allocatable :: zonal_proc_id(:)

    ! Create zonal communicator for reduce algorithm.
    if (proc%idom == 1) then ! Only root domain has reduce region.
      jr = 0
      do j = 1, size(reduce_factors)
        if (reduce_factors(j) > 0) then
          jr = j
        else if (jr /= 0) then
          exit
        end if
      end do
      if (global_mesh%is_south_pole(proc%lat_ibeg) .or. global_mesh%is_north_pole(proc%lat_iend) .or. &
#ifdef V_POLE
          proc%lat_ibeg <= jr .or. proc%lat_iend > global_mesh%num_half_lat - jr) then
#else
          proc%lat_ibeg <= jr .or. proc%lat_iend > global_mesh%num_full_lat - jr) then
#endif
        call log_notice('Create zonal communicator on process ' // to_string(proc%id) // '.')
        allocate(zonal_proc_id(proc%cart_dims(1)))
        do i = 1, proc%cart_dims(1)
          call MPI_CART_RANK(proc%local_comm, [i-1,proc%cart_coords(2)], zonal_proc_id(i), ierr)
        end do
        call MPI_GROUP_INCL(proc%local_group, size(zonal_proc_id), zonal_proc_id, proc%zonal_group, ierr)
        call MPI_COMM_CREATE_GROUP(proc%local_comm, proc%zonal_group, sum(zonal_proc_id), proc%zonal_comm, ierr)
        deallocate(zonal_proc_id)
      end if
    end if

  end subroutine setup_zonal_comm_for_reduce

  subroutine connect_parent()

  end subroutine connect_parent

  subroutine process_create_blocks()

    integer i, dtype

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init(proc%id, global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                             proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    select case (r8)
    case (4)
      dtype = MPI_REAL
    case (8)
      dtype = MPI_DOUBLE
    case (16)
      dtype = MPI_REAL16
    case default
      call log_error('Unsupported parameter r8!')
    end select

    ! Setup halos (only normal halos for the time being).
    allocate(proc%blocks(1)%halo(size(proc%ngb)))
    do i = 1, size(proc%ngb)
      select case (proc%ngb(i)%orient)
      case (west, east)
        call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype, ngb_proc_id=proc%ngb(i)%id, &
                                         global_lat_ibeg=proc%ngb(i)%lat_ibeg, global_lat_iend=proc%ngb(i)%lat_iend, proc_id=proc%id)
      case (south, north)
        call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype, ngb_proc_id=proc%ngb(i)%id, &
                                         global_lon_ibeg=proc%ngb(i)%lon_ibeg, global_lon_iend=proc%ngb(i)%lon_iend, proc_id=proc%id)
      end select
    end do

  end subroutine process_create_blocks

  subroutine process_neighbor_init(this, orient, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(process_neighbor_type), intent(inout) :: this
    integer, intent(in) :: orient
    integer, intent(in), optional :: lon_ibeg
    integer, intent(in), optional :: lon_iend
    integer, intent(in), optional :: lat_ibeg
    integer, intent(in), optional :: lat_iend

    this%orient = orient

    select case (orient)
    case (west, east)
      this%lat_ibeg = lat_ibeg
      this%lat_iend = lat_iend
    case (south, north)
      this%lon_ibeg = lon_ibeg
      this%lon_iend = lon_iend
    end select

  end subroutine process_neighbor_init

end module process_mod
