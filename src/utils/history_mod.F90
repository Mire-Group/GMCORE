module history_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write_state
  public history_write_debug

contains

  subroutine history_init()

    character(10) time_value, time_units
    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4), lon_dims_2d(3)
    character(4) lat_dims(4), lat_dims_2d(3)
    character(4) vtx_dims(4), vtx_dims_2d(3)
    character(4) cell_lev_dims(4), lon_lev_dims(4), lat_lev_dims(4), lev_dims(4)
    real(8) seconds

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select

        cell_dims(1) =  'lon';     cell_dims(2) =  'lat';     cell_dims(3) =  'lev';     cell_dims(4) = 'time'
    cell_lev_dims(1) =  'lon'; cell_lev_dims(2) =  'lat'; cell_lev_dims(3) = 'ilev'; cell_lev_dims(4) = 'time'
         lon_dims(1) = 'ilon';      lon_dims(2) =  'lat';      lon_dims(3) =  'lev';      lon_dims(4) = 'time'
     lon_lev_dims(1) = 'ilon';  lon_lev_dims(2) =  'lat';  lon_lev_dims(3) = 'ilev';  lon_lev_dims(4) = 'time'
         lat_dims(1) =  'lon';      lat_dims(2) = 'ilat';      lat_dims(3) =  'lev';      lat_dims(4) = 'time'
     lat_lev_dims(1) =  'lon';  lat_lev_dims(2) = 'ilat';  lat_lev_dims(3) = 'ilev';  lat_lev_dims(4) = 'time'
         vtx_dims(1) = 'ilon';      vtx_dims(2) = 'ilat';      vtx_dims(3) =  'lev';      vtx_dims(4) = 'time'
         lev_dims(1) =  'lon';      lev_dims(2) =  'lat';      lev_dims(3) = 'ilev';      lev_dims(4) = 'time'
     cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'
      lon_dims_2d(1) = 'ilon';   lon_dims_2d(2) =  'lat';   lon_dims_2d(3) = 'time'
      lat_dims_2d(1) =  'lon';   lat_dims_2d(2) = 'ilat';   lat_dims_2d(3) = 'time'
      vtx_dims_2d(1) = 'ilon';   vtx_dims_2d(2) = 'ilat';   vtx_dims_2d(3) = 'time'

    call fiona_init(time_units, start_time_str)

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm)
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('h0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
      call fiona_add_dim('h0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
      call fiona_add_var('h0', 't'    , long_name='temperature'                 , units='K'      , dim_names=cell_dims)
      call fiona_add_var('h0', 't850' , long_name='temperature on 850hPa'       , units='K'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 't700' , long_name='temperature on 700hPa'       , units='K'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'pt'   , long_name='potential temperature'       , units='K'      , dim_names=cell_dims)
      call fiona_add_var('h0', 'phs'  , long_name='surface hydrostatic pressure', units='Pa'     , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'ph'   , long_name='hydrostatic pressure'        , units='Pa'     , dim_names=cell_dims)
      call fiona_add_var('h0', 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims)
      call fiona_add_var('h0', 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims)
      call fiona_add_var('h0', 'z'    , long_name='height'                      , units='m'      , dim_names=cell_dims)
      call fiona_add_var('h0', 'zs'   , long_name='surface height'              , units='m'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'         , units='m-1 s-1', dim_names=vtx_dims)
      call fiona_add_var('h0', 'vor'  , long_name='relative vorticity'          , units='s-1'    , dim_names=vtx_dims)
      call fiona_add_var('h0', 'tm'   , long_name='total mass'                  , units='m'      , dim_names=['time'])
      call fiona_add_var('h0', 'te'   , long_name='total energy'                , units='m4 s-4' , dim_names=['time'], data_type='real(8)')
      call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
      call fiona_add_var('h0', 'tpv'  , long_name='total potential vorticity'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
    else
      call fiona_add_var('h0', 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims_2d)
      call fiona_add_var('h0', 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims_2d)
      call fiona_add_var('h0', 'z'    , long_name='height'                      , units='m'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'zs'   , long_name='surface height'              , units='m'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'         , units='m-1 s-1', dim_names=vtx_dims_2d)
      call fiona_add_var('h0', 'vor'  , long_name='relative vorticity'          , units='s-1'    , dim_names=vtx_dims_2d)
      call fiona_add_var('h0', 'tm'   , long_name='total mass'                  , units='m'      , dim_names=['time'])
      call fiona_add_var('h0', 'te'   , long_name='total energy'                , units='m4 s-4' , dim_names=['time'], data_type='real(8)')
      call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
      call fiona_add_var('h0', 'tpv'  , long_name='total potential vorticity'   , units='m2 s-5' , dim_names=['time'], data_type='real(8)')
    end if

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm)
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('h1', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
      call fiona_add_dim('h1', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
      call fiona_add_var('h1', 'dphs'     , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'pt_lev'   , long_name='potential temperature on half levels'          , units='', dim_names=cell_lev_dims)
      call fiona_add_var('h1', 'ph_lev'   , long_name='hydrostatic pressure on half levels'           , units='', dim_names=cell_lev_dims)
      call fiona_add_var('h1', 'dptfdlon' , long_name='zonal potential temperature flux gradient'     , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dptfdlat' , long_name='meridional potential temperature flux gradient', units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dptfdlev' , long_name='vertical potential temperature flux gradient'  , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'wedphdlev', long_name='vertical coordinate velocity'                  , units='', dim_names=lev_dims)
      call fiona_add_var('h1', 'wedphdlev_lon', long_name='vertical coordinate velocity'              , units='', dim_names=lon_lev_dims)
      call fiona_add_var('h1', 'wedphdlev_lat', long_name='vertical coordinate velocity'              , units='', dim_names=lat_lev_dims)
      call fiona_add_var('h1', 'dpdlon'   , long_name='zonal pressure gradient force'                 , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'wedudlev' , long_name='vertical advection of u'                       , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'dpdlat'   , long_name='meridional pressure gradient force'            , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'wedvdlev' , long_name='vertical advection of v'                       , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'qhv'      , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'qhu'      , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dpedlon'  , long_name='zonal geopotential energy gradient force'      , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'dkedlon'  , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'dpedlat'  , long_name='meridional geopotential energy gradient force' , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dkedlat'  , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dmfdlon'  , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dmfdlat'  , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'mf_lon_n' , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'mf_lon_t' , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'mf_lat_n' , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'mf_lat_t' , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'm'        , long_name='dph on full levels'                            , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'm_lon'    , long_name='dph on lon edges on full levels'               , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'm_lat'    , long_name='dph on lat edges on full levels'               , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'm_vtx'    , long_name='dph on vertices on full levels'                , units='', dim_names=vtx_dims)
      call fiona_add_var('h1', 'pv_lon'   , long_name='pv on U grid'                                  , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'pv_lat'   , long_name='pv on V grid'                                  , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'ke'       , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'wp'       , long_name='omega'                                         , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'div'      , long_name='divergence'                                    , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'div2'     , long_name='laplacian of divergence'                       , units='', dim_names=cell_dims)
    else
      call fiona_add_var('h1', 'qhv'      , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'qhu'      , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dpedlon'  , long_name='zonal geopotential energy gradient force'      , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'dkedlon'  , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'dpedlat'  , long_name='meridional geopotential energy gradient force' , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dkedlat'  , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dmfdlon'  , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'dmfdlat'  , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'mf_lon_n' , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'mf_lon_t' , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'mf_lat_n' , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'mf_lat_t' , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'm_lon'    , long_name='mass on lon edges'                             , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'm_lat'    , long_name='mass on lat edges'                             , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'm_vtx'    , long_name='mass on vertices'                              , units='', dim_names=vtx_dims_2d)
      call fiona_add_var('h1', 'pv_lon'   , long_name='pv on U grid'                                  , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'pv_lat'   , long_name='pv on V grid'                                  , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'ke'       , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d)
    end if

    call time_add_alert('history_write', seconds=seconds)

  end subroutine history_init

  subroutine history_final()

  end subroutine history_final

  subroutine history_write_state(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    call fiona_start_output('h0', dble(elapsed_seconds), new_file=time_step==0)
    call fiona_output('h0', 'lon' , dble(global_mesh%full_lon_deg(1:global_mesh%num_full_lon)))
    call fiona_output('h0', 'lat' , dble(global_mesh%full_lat_deg(1:global_mesh%num_full_lat)))
    call fiona_output('h0', 'ilon', dble(global_mesh%half_lon_deg(1:global_mesh%num_half_lon)))
    call fiona_output('h0', 'ilat', dble(global_mesh%half_lat_deg(1:global_mesh%num_half_lat)))
    if (baroclinic) then
      call fiona_output('h0', 'lev' , dble(global_mesh%full_lev))
      call fiona_output('h0', 'ilev', dble(global_mesh%half_lev))
    end if

    do iblk = 1, size(blocks)
      mesh => blocks(iblk)%mesh
      state => blocks(iblk)%state(itime)
      static => blocks(iblk)%static

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('h0', 'zs' , dble(static%gzs(is:ie,js:je      ) / g), start=start, count=count)
      call fiona_output('h0', 'z'  , dble(state %gz (is:ie,js:je,ks:ke) / g), start=start, count=count)

      if (baroclinic) then
        call fiona_output('h0', 't'     , dble(state%t     (is:ie,js:je,ks:ke)), start=start, count=count)
        call fiona_output('h0', 't850'  , dble(state%t850  (is:ie,js:je      )), start=start, count=count)
        call fiona_output('h0', 't700'  , dble(state%t700  (is:ie,js:je      )), start=start, count=count)
        call fiona_output('h0', 'pt'    , dble(state%pt    (is:ie,js:je,ks:ke)), start=start, count=count)
        call fiona_output('h0', 'phs'   , dble(state%phs   (is:ie,js:je      )), start=start, count=count)
        call fiona_output('h0', 'ph'    , dble(state%ph    (is:ie,js:je,ks:ke)), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('h0', 'u'  , dble(state %u  (is:ie,js:je,ks:ke)), start=start, count=count)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output('h0', 'v'  , dble(state %v  (is:ie,js:je,ks:ke)), start=start, count=count)

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output('h0', 'pv' , dble(state %pv (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h0', 'vor', dble(state %vor(is:ie,js:je,ks:ke)), start=start, count=count)

      call fiona_output('h0', 'tm' , dble(state %tm ))
      call fiona_output('h0', 'te' , dble(state %te ))
      call fiona_output('h0', 'tpe', dble(state %tpe))
      call fiona_output('h0', 'tpv', dble(state %tav))
    end do
    call fiona_end_output('h0')

  end subroutine history_write_state

  subroutine history_write_debug(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(tend_type), pointer :: tend

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    mesh => blocks(1)%mesh
    state => blocks(1)%state(itime)
    tend => blocks(1)%tend(itime)

    call fiona_start_output('h1', dble(elapsed_seconds), new_file=time_step==0)
    call fiona_output('h1', 'lon'   , dble(global_mesh%full_lon_deg(1:global_mesh%num_full_lon)))
    call fiona_output('h1', 'lat'   , dble(global_mesh%full_lat_deg(1:global_mesh%num_full_lat)))
    call fiona_output('h1', 'ilon'  , dble(global_mesh%half_lon_deg(1:global_mesh%num_half_lon)))
    call fiona_output('h1', 'ilat'  , dble(global_mesh%half_lat_deg(1:global_mesh%num_half_lat)))
    if (baroclinic) then
      call fiona_output('h1', 'lev' , dble(global_mesh%full_lev))
      call fiona_output('h1', 'ilev', dble(global_mesh%half_lev))
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

    call fiona_output('h1', 'dmfdlon' , dble(tend%dmfdlon  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'dmfdlat' , dble(tend%dmfdlat  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'ke'      , dble(state%ke      (is:ie,js:je,ks:ke)), start=start, count=count)
    if (baroclinic) then
      call fiona_output('h1', 'm'       , dble(state%m      (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'dphs'    , dble(tend%dphs    (is:ie,js:je      )), start=start, count=count)
      call fiona_output('h1', 'dptfdlon', dble(tend%dptfdlon(is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'dptfdlat', dble(tend%dptfdlat(is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'dptfdlev', dble(tend%dptfdlev(is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'wp'      , dble(state%wp     (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'div'     , dble(state%div    (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'div2'    , dble(state%div2   (is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

    call fiona_output('h1', 'qhv'     , dble(tend%qhv      (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'dpedlon' , dble(tend%dpedlon  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'dkedlon' , dble(tend%dkedlon  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'm_lon'   , dble(state%m_lon   (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'mf_lon_n', dble(state%mf_lon_n(is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'mf_lon_t', dble(state%mf_lon_t(is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'pv_lon'  , dble(state%pv_lon  (is:ie,js:je,ks:ke)), start=start, count=count)

    if (baroclinic) then
      call fiona_output('h1', 'dpdlon'  , dble(tend%dpdlon  (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'wedudlev', dble(tend%wedudlev(is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

    call fiona_output('h1', 'qhu'     , dble(tend%qhu      (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'dpedlat' , dble(tend%dpedlat  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'dkedlat' , dble(tend%dkedlat  (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'm_lat'   , dble(state%m_lat   (is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'mf_lat_n', dble(state%mf_lat_n(is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'mf_lat_t', dble(state%mf_lat_t(is:ie,js:je,ks:ke)), start=start, count=count)
    call fiona_output('h1', 'pv_lat'  , dble(state%pv_lat  (is:ie,js:je,ks:ke)), start=start, count=count)

    if (baroclinic) then
      call fiona_output('h1', 'dpdlat'  , dble(tend%dpdlat  (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'wedvdlev', dble(tend%wedvdlev(is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]

    call fiona_output('h1', 'm_vtx', dble(state%m_vtx(is:ie,js:je,ks:ke)), start=start, count=count)

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]

    if (baroclinic) then
      call fiona_output('h1', 'wedphdlev', dble(state%wedphdlev(is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'pt_lev'   , dble(state%pt_lev   (is:ie,js:je,ks:ke)), start=start, count=count)
      call fiona_output('h1', 'ph_lev'   , dble(state%ph_lev   (is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_half_lev]

    if (baroclinic) then
      call fiona_output('h1', 'wedphdlev_lon', dble(state%wedphdlev_lon(is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_half_lev]

    if (baroclinic) then
      call fiona_output('h1', 'wedphdlev_lat', dble(state%wedphdlev_lat(is:ie,js:je,ks:ke)), start=start, count=count)
    end if

    call fiona_end_output('h1')

  end subroutine history_write_debug

end module history_mod
