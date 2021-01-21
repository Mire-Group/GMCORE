module refer_state_mod

  use flogger
  use namelist_mod
  use formula_mod
  use process_mod
  use static_mod
  use refer_state_wrf_mod
  use refer_state_types_mod

  implicit none

  private

  type(refer_state_type) refer_state

contains

  subroutine refer_state_init(static)

    type(static_type), intent(in) :: static

    integer i, j, k, im

    select case (refer_state_scheme)
    case ('wrf')
      call refer_state_wrf_init(static, refer_state)
    case default
      if (is_root_proc()) call log_error('Unknown refer_state_scheme ' // trim(refer_state_scheme) // '!')
    end select

    do k = refer_state%mesh%full_lev_ibeg, refer_state%mesh%full_lev_iend
      do j = refer_state%mesh%full_lat_ibeg, refer_state%mesh%full_lat_iend
        do i = refer_state%mesh%full_lon_ibeg, refer_state%mesh%full_lon_iend
          do im = 1 , member_num
            refer_state%pt(im,i,j,k) = potential_temperature(refer_state%t(im,i,j,k), refer_state%ph(im,i,j,k))
            refer_state%rhod(im,i,j,k) = dry_air_density(refer_state%pt(im,i,j,k), refer_state%ph(im,i,j,k))
          end do
        end do
      end do
    end do

  end subroutine refer_state_init

end module refer_state_mod
