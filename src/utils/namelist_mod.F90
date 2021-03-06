module namelist_mod

  use string
  use const_mod
  use time_mod, start_time => start_time_array, end_time => end_time_array

  implicit none

  character(256)  :: case_desc            = 'N/A'
  character(256)  :: case_name            = 'N/A'
  character(30 )  :: test_case            = 'N/A'
  character(30 )  :: history_interval(1)  = 'N/A'
  character(30 )  :: restart_interval     = 'N/A'
  character(30 )  :: print_interval       = '1 hours'
  character(256)  :: initial_file         = 'N/A'
  character(256)  :: restart_file         = 'N/A'

  logical         :: restart              = .false.

  integer         :: num_lon
  integer         :: num_lat
  integer         :: num_lev              = 1

  logical         :: baroclinic           = .false.
  logical         :: hydrostatic          = .true.
  logical         :: nonhydrostatic       = .false.

  integer         :: num_proc_lon(20)     = 0
  integer         :: num_proc_lat(20)     = 0

  character(30)   :: tangent_wgt_scheme   = 'classic'

  character(30)   :: vert_coord_scheme    = 'hybrid'
  character(30)   :: vert_coord_template  = 'N/A'
  character(30)   :: refer_state_scheme   = 'wrf'
  real(r8)        :: ptop                 = 2.194e2_r8

  integer         :: ke_scheme            = 1
  real(r8)        :: ke_cell_wgt          = 3.0_r8 / 8.0_r8

  integer         :: pv_scheme            = 1
  logical         :: pv_pole_stokes       = .true.
  character(8)    :: pgf_scheme           = 'lin97'
  integer         :: coriolis_scheme      = 1
  integer         :: upwind_order         = -1 ! -1, 1, 3
  real(r8)        :: upwind_wgt           = 1.0_r8
  real(r8)        :: upwind_wgt_pt        = 0.25_r8

  integer         :: fast_cycles          = 1
  character(30)   :: split_scheme         = ''
  character(30)   :: time_scheme          = 'pc2'

  real(r8)        :: coarse_pole_mul      = 0
  real(r8)        :: coarse_pole_decay    = 100.0

  ! Reduce settings
  integer         :: reduce_factors(100)  = 0
  logical         :: reduce_pv_directly   = .true.
  logical         :: do_reduce_ke         = .true.

  ! Damping settings
  logical         :: use_polar_damp       = .false.
  integer         :: polar_damp_order     = 4
  real(r8)        :: polar_damp_lat0      = 80.0_r8
  logical         :: use_div_damp         = .false.
  integer         :: div_damp_order       = 2
  integer         :: div_damp_j0          = 0
  real(r8)        :: div_damp_imp_lat0    = 80
  real(r8)        :: div_damp_upper       = 2.0_r8
  real(r8)        :: div_damp_polar       = 0.5_r8
  real(r8)        :: div_damp_exp         = 0.01_r8
  real(r8)        :: div_damp_coef2       = 1.0_r8 / 128.0_r8
  real(r8)        :: div_damp_coef4       = 0.01_r8
  logical         :: use_vor_damp         = .false.
  integer         :: vor_damp_order       = 2
  real(r8)        :: vor_damp_imp_lat0    = 80
  real(r8)        :: vor_damp_lat0        = 70.0_r8
  real(r8)        :: vor_damp_decay       = 0.2_r8
  real(r8)        :: vor_damp_coef2       = 0.001_r8
  real(r8)        :: vor_damp_coef4       = 0.01_r8
  real(r8)        :: vor_damp_bkg_coef    = 0.0001_r8
  logical         :: use_rayleigh_damp    = .false.
  logical         :: use_smag_damp        = .false.
  real(r8)        :: smag_damp_coef       = 0.2

  ! Output settings
  character(8)    :: output_h0_vars(100)  = ''
  logical         :: output_debug         = .false.

  ! Nest settings
  character(30)   :: nest_time_scheme     = 'pc2'
  integer         :: nest_max_dom         = 1
  integer         :: nest_parent_id(20)   = 1
  real(r8)        :: nest_lon_beg(20)     = inf
  real(r8)        :: nest_lon_end(20)     = inf
  real(r8)        :: nest_lat_beg(20)     = inf
  real(r8)        :: nest_lat_end(20)     = inf

  !muti prepare
  integer         :: file_num = 1
  integer         :: initial_interval = 0

  namelist /gmcore_control/     &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    num_lon                   , &
    num_lat                   , &
    num_lev                   , &
    nonhydrostatic            , &
    num_proc_lon              , &
    num_proc_lat              , &
    start_time                , &
    end_time                  , &
    dt_in_seconds             , &
    run_hours                 , &
    run_days                  , &
    history_interval          , &
    restart_interval          , &
    print_interval            , &
    initial_file              , &
    restart_file              , &
    restart                   , &
    tangent_wgt_scheme        , &
    vert_coord_scheme         , &
    vert_coord_template       , &
    refer_state_scheme        , &
    ptop                      , &
    ke_scheme                 , &
    ke_cell_wgt               , &
    pv_scheme                 , &
    pv_pole_stokes            , &
    pgf_scheme                , &
    coriolis_scheme           , &
    upwind_order              , &
    upwind_wgt                , &
    upwind_wgt_pt             , &
    split_scheme              , &
    time_scheme               , &
    fast_cycles               , &
    coarse_pole_mul           , &
    coarse_pole_decay         , &
    reduce_factors            , &
    reduce_pv_directly        , &
    do_reduce_ke              , &
    use_polar_damp            , &
    polar_damp_order          , &
    polar_damp_lat0           , &
    use_div_damp              , &
    div_damp_order            , &
    div_damp_imp_lat0         , &
    div_damp_polar            , &
    div_damp_upper            , &
    div_damp_j0               , &
    div_damp_exp              , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    use_vor_damp              , &
    vor_damp_order            , &
    vor_damp_imp_lat0         , &
    vor_damp_lat0             , &
    vor_damp_decay            , &
    vor_damp_coef2            , &
    vor_damp_coef4            , &
    vor_damp_bkg_coef         , &
    use_rayleigh_damp         , &
    use_smag_damp             , &
    smag_damp_coef            , &
    output_h0_vars            , &
    output_debug              , &
    nest_time_scheme          , &
    nest_max_dom              , &
    nest_parent_id            , &
    nest_lon_beg              , &
    nest_lon_end              , &
    nest_lat_beg              , &
    nest_lat_end              , &
    file_num                  , &
    initial_interval       

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_control)
    close(10)

    hydrostatic = .not. nonhydrostatic

  end subroutine parse_namelist

  subroutine print_namelist()

      write(*, *) '=================== GMCORE Parameters ==================='
      write(*, *) 'num_lon             = ', to_str(num_lon)
      write(*, *) 'num_lat             = ', to_str(num_lat)
      write(*, *) 'num_lev             = ', to_str(num_lev)
      write(*, *) 'num_proc_lon        = ', to_str(num_proc_lon)
      write(*, *) 'num_proc_lat        = ', to_str(num_proc_lat)
    if (coarse_pole_mul /= 0) then
      write(*, *) 'coarse_pole_mul     = ', to_str(coarse_pole_mul, 3)
      write(*, *) 'coarse_pole_decay   = ', to_str(coarse_pole_decay, 3)
    end if
      write(*, *) 'hydrostatic         = ', to_str(hydrostatic)
      write(*, *) 'nonhydrostatic      = ', to_str(nonhydrostatic)
      write(*, *) 'vert_coord_scheme   = ', trim(vert_coord_scheme)
      write(*, *) 'vert_coord_template = ', trim(vert_coord_template)
      write(*, *) 'ptop                = ', to_str(ptop, 4)
      write(*, *) 'dt_in_seconds       = ', to_str(dt_in_seconds, 2)
      write(*, *) 'pgf_scheme          = ', trim(pgf_scheme)
      write(*, *) 'ke_scheme           = ', to_str(ke_scheme)
    if (ke_scheme == 2) then
      write(*, *) 'ke_cell_wgt         = ', to_str(ke_cell_wgt, 2)
    end if
      write(*, *) 'pv_scheme           = ', to_str(pv_scheme)
      write(*, *) 'pv_pole_stokes      = ', to_str(pv_pole_stokes)
      write(*, *) 'time_scheme         = ', trim(time_scheme)
      write(*, *) 'upwind_order        = ', to_str(upwind_order)
    if (upwind_order > 0) then
      write(*, *) 'upwind_wgt          = ', to_str(upwind_wgt, 2)
      write(*, *) 'upwind_wgt_pt       = ', to_str(upwind_wgt_pt, 2)
    end if
      write(*, *) 'reduce_pv_directly  = ', to_str(reduce_pv_directly)
      write(*, *) 'do_reduce_ke        = ', to_str(do_reduce_ke)
      write(*, *) 'use_div_damp        = ', to_str(use_div_damp)
    if (use_div_damp) then
      write(*, *) 'div_damp_coef2      = ', to_str(div_damp_coef2, 3)
    end if
      write(*, *) 'use_vor_damp        = ', to_str(use_vor_damp)
    if (use_vor_damp) then
      write(*, *) 'vor_damp_lat0       = ', to_str(vor_damp_lat0, 1)
      write(*, *) 'vor_damp_decay      = ', to_str(vor_damp_decay, 1)
      write(*, *) 'vor_damp_coef2      = ', to_str(vor_damp_coef2, 3)
    end if
      write(*, *) 'use_polar_damp      = ', to_str(use_polar_damp)
      write(*, *) 'use_rayleigh_damp   = ', to_str(use_rayleigh_damp)
      write(*, *) 'use_smag_damp       = ', to_str(use_smag_damp)
    if (use_smag_damp) then
      write(*, *) 'smag_damp_coef      = ', to_str(smag_damp_coef, 1)
    end if
      write(*, *) 'history_interval    = ', history_interval(1)
      write(*, *) 'print_interval      = ', print_interval
      write(*, *) '========================================================='

  end subroutine print_namelist

end module namelist_mod
