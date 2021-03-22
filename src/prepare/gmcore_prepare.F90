program gmcore_prepare

  use fiona
  use string
  use topo_mod
  use bkg_mod
  use mesh_mod
  use process_mod
  use block_mod
  use vert_coord_mod
  use initial_mod
  use namelist_mod
  use file_date_mod

  implicit none

  character(256) namelist_file
  character(256) topo_file
  character(256) bkg_file
  character(30) :: bkg_type = 'era5'
  character(30) :: initial_time = '1970-01-01'

  real(r8) :: zero_min_lon(100) = -1.0e33
  real(r8) :: zero_max_lon(100) = -1.0e33
  real(r8) :: zero_min_lat(100) = -1.0e33
  real(r8) :: zero_max_lat(100) = -1.0e33
  real(r8) :: smth_min_lon(100) = -1.0e33
  real(r8) :: smth_max_lon(100) = -1.0e33
  real(r8) :: smth_min_lat(100) = -1.0e33
  real(r8) :: smth_max_lat(100) = -1.0e33
  integer :: smth_steps(100) = 1

  integer iblk, i , iter_file
  integer ini_minute , ini_hour , ini_day , ini_month , ini_year
  integer mpas_len , gmcore_len
  character(30) :: s_id

  namelist /gmcore_prepare_params/ &
    initial_time                 , &
    topo_file                    , &
    bkg_type                     , &
    bkg_file                     , &
    initial_file                 , &
    num_lon                      , &
    num_lat                      , &
    num_lev                      , &
    coarse_pole_mul              , &
    coarse_pole_decay            , &
    vert_coord_scheme            , &
    vert_coord_template          , &
    zero_min_lon                 , &
    zero_max_lon                 , &
    zero_min_lat                 , &
    zero_max_lat                 , &
    smth_min_lon                 , &
    smth_max_lon                 , &
    smth_min_lat                 , &
    smth_max_lat                 , &
    smth_steps                   , &
    file_num                     , &
    initial_interval 

  call get_command_argument(1, namelist_file)

  open(10, file=namelist_file)
  read(10, nml=gmcore_prepare_params)
  close(10)

  write(*, *) '=================== GMCORE Parameters ==================='
  write(*, *) 'num_lon              = ', to_str(num_lon)
  write(*, *) 'num_lat              = ', to_str(num_lat)
  write(*, *) 'num_lev              = ', to_str(num_lev)
  if (coarse_pole_mul /= 0) then
  write(*, *) 'coarse_pole_mul      = ', to_str(coarse_pole_mul, 2)
  write(*, *) 'coarse_pole_decay    = ', to_str(coarse_pole_decay, 2)
  end if
  write(*, *) 'vert_coord_scheme    = ', trim(vert_coord_scheme)
  write(*, *) 'vert_coord_template  = ', trim(vert_coord_template)
  write(*, *) 'initial_time         = ', trim(initial_time)
  write(*, *) 'namelist_file        = ', trim(namelist_file)
  write(*, *) 'topo_file            = ', trim(topo_file)
  write(*, *) 'bkg_file             = ', trim(bkg_file)
  write(*, *) 'initial_file         = ', trim(initial_file)
  write(*, *) 'bkg_type             = ', trim(bkg_type)
  write(*, *) 'file_num            = ', to_str(file_num)
  if (file_num > 1) then
    write(*, *) 'initial_interval    = ', to_str(initial_interval)
  end if
  write(*, *) '========================================================='

  time_scheme = 'N/A'

  call fiona_init(start_time=initial_time, time_units='hours')

  call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=2, lat_halo_width=2)
  call process_init()
  call vert_coord_init(num_lev, scheme=vert_coord_scheme, template=vert_coord_template)
  call process_create_blocks()

  call topo_read(topo_file)

  do iblk = 1, size(proc%blocks)
    call topo_regrid(proc%blocks(iblk))
    do i = 1, size(zero_min_lon)
      if (zero_min_lon(i) /= -1.0e33 .and. zero_max_lon(i) /= -1.0e33 .and. &
          zero_min_lat(i) /= -1.0e33 .and. zero_max_lat(i) /= -1.0e33) then
        call topo_zero(proc%blocks(iblk), zero_min_lon(i), zero_max_lon(i), zero_min_lat(i), zero_max_lat(i))
      end if
    end do
    do i = 1, size(smth_min_lon)
      if (smth_min_lon(i) /= -1.0e33 .and. smth_max_lon(i) /= -1.0e33 .and. &
          smth_min_lat(i) /= -1.0e33 .and. smth_max_lat(i) /= -1.0e33) then
        call topo_smth(proc%blocks(iblk), smth_min_lon(i), smth_max_lon(i), smth_min_lat(i), smth_max_lat(i), smth_steps(i))
      end if
    end do
  end do

  mpas_len = len_trim(bkg_file)
  gmcore_len = len_trim(initial_file)
  read(bkg_file(mpas_len - 21 : mpas_len - 18) , *) ini_year
  read(bkg_file(mpas_len - 16 : mpas_len - 15) , *) ini_month
  read(bkg_file(mpas_len - 13 : mpas_len - 12) , *) ini_day
  read(bkg_file(mpas_len - 10 : mpas_len -  9) , *) ini_hour
  read(bkg_file(mpas_len -  7 : mpas_len -  6) , *) ini_minute

  do iter_file = 1 , file_num

    initial_file(gmcore_len - 21 : gmcore_len - 18) = bkg_file(mpas_len - 21 : mpas_len - 18)
    initial_file(gmcore_len - 16 : gmcore_len - 15) = bkg_file(mpas_len - 16 : mpas_len - 15)
    initial_file(gmcore_len - 13 : gmcore_len - 12) = bkg_file(mpas_len - 13 : mpas_len - 12)
    initial_file(gmcore_len - 10 : gmcore_len - 9 ) = bkg_file(mpas_len - 10 : mpas_len - 9 )
    initial_file(gmcore_len - 7  : gmcore_len - 6 ) = bkg_file(mpas_len - 7  : mpas_len - 6 )

    call bkg_read(bkg_type, bkg_file)

    call bkg_regrid_phs()
    call bkg_calc_ph()
    call bkg_regrid_pt()
    call bkg_regrid_u()
    call bkg_regrid_v()

    call initial_write(initial_file)

    call date_advance(ini_year , ini_month , ini_day , ini_hour , ini_minute)

    write(s_id,"(i4.4)") ini_year
    bkg_file(mpas_len - 21 : mpas_len - 18) = s_id
    s_id = 'N/A'

    write(s_id,"(i2.2)") ini_month
    bkg_file(mpas_len - 16 : mpas_len - 15) = s_id
    s_id = 'N/A'
    
    write(s_id,"(i2.2)") ini_day
    bkg_file(mpas_len - 13 : mpas_len - 12) = s_id
    s_id = 'N/A'

    write(s_id,"(i2.2)") ini_hour
    bkg_file(mpas_len - 10 : mpas_len - 9)  = s_id
    s_id = 'N/A'

    write(s_id,"(i2.2)") ini_minute
    bkg_file(mpas_len - 7  : mpas_len - 6)  = s_id
    s_id = 'N/A'



  end do

  call topo_final()
  call bkg_final()
  call process_final()

end program gmcore_prepare
